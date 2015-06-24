/*
 * InputStreamMgr.cpp
 *
 *  Created on: Mar 21, 2013
 *      Author: nek3d
 */

#include "InputStreamMgr.h"
#include <cstring> //for memset
#include "gzstream.h"
#include "CompressionTools.h"

const char *InputStreamMgr::FIFO_STRING_LITERAL = "/dev/fd";

InputStreamMgr::InputStreamMgr(const QuickString &filename, bool buildScanBuffer)
:
 _filename(filename),
 _pushBackStreamBuf(NULL),
 _inputFileStream(NULL),
 _infStreamBuf(NULL),
 _oldInputStream(NULL),
 _isStdin(false),
 _isGzipped(false),
 _isBam(false),
 _isBgzipped(false),
 _tmpZipBuf(NULL),
 _bamRuledOut(false),
 _streamFinished(false),
 _numBytesInBuffer(0),
 _bamReader(NULL),
 _bgStream(NULL),
 _eofHit(false)
{
	_possibleBamCode.resize(4, 0);
}


InputStreamMgr::~InputStreamMgr() {
	delete _pushBackStreamBuf;
	_pushBackStreamBuf = NULL;

	delete _inputFileStream;
	_inputFileStream = NULL;

	delete _oldInputStream;
	_oldInputStream = NULL;

	delete _infStreamBuf;
	_infStreamBuf = NULL;

	delete _bamReader;
	_bgStream = NULL;

	delete _bgStream;
	_bgStream = NULL;

	delete _finalInputStream;
	_finalInputStream = NULL;

	delete [] _tmpZipBuf;
	_tmpZipBuf = NULL;
}

bool InputStreamMgr::init()
{
	if (_filename == "-" || _filename == "stdin") { //stdin
		_isStdin = true;
		//peek at the first char of stdin to see if this is gzipped.
		if ((unsigned char)cin.peek()  == 0x1f) {
			_isGzipped = true;
		}
		_pushBackStreamBuf = new PushBackStreamBuf(cin.rdbuf());
	} else {
		if (strncmp(_filename.c_str(), FIFO_STRING_LITERAL, strlen(FIFO_STRING_LITERAL)) == 0) {
			_isStdin = true;
		}
		_inputFileStream = new ifstream(_filename.c_str());
		if (_inputFileStream->fail()) {
			cerr << "Error: Unable to open file " << _filename << ". Exiting." << endl;
			delete _inputFileStream;
			_inputFileStream = NULL;
			exit(1);
		}
		//peek at the first char of stdin to see if this is gzipped.
		if ((unsigned char)_inputFileStream->peek()  == 0x1f) {
			_isGzipped = true;
		}
		_pushBackStreamBuf = new PushBackStreamBuf(_inputFileStream->rdbuf());
	}
	//now we have a PushBackStreamBuf. Make a new stream.
	_finalInputStream = new istream(_pushBackStreamBuf);
	populateScanBuffer();
	return true;
}

int InputStreamMgr::read(char *data, size_t dataSize)
{
	size_t origRead = 0;
	if (!_saveDataStr.empty()) {
		//must first copy contents of savedData into requested data read buffer.
		if (dataSize >= _saveDataStr.size()) {
			//They asked for the same amount of data or more than we saved. Give them all the saved data,
			//then decrement the requested data size accordingly.
			origRead = _saveDataStr.size();
			memcpy(data, _saveDataStr.c_str(), origRead);
			data += origRead;
			dataSize -= origRead;
			_saveDataStr.clear();
		} else {
			//This part is tricky. They want less data than we saved. Give them what they
			//requested, then delete from the front of the saveDataStr by using it's substr method.
			memcpy(data, _saveDataStr.c_str(), dataSize);
			QuickString newDataStr;
			_saveDataStr.substr(newDataStr, dataSize, _saveDataStr.size() - dataSize);
			_saveDataStr = newDataStr;
			return dataSize;
		}
	}
	if (_streamFinished) {
		return origRead;
	}
	if (_isBgzipped) {
		return (int)(origRead + _bgStream->Read(data, dataSize));
	}
	_finalInputStream->read(data, dataSize);
	return origRead + _finalInputStream->gcount();
}

bool InputStreamMgr::populateScanBuffer()
{
	_scanBuffer.clear();
	_saveDataStr.clear();
	int numChars=0;
	int currChar = 0;
	while (1) {
		if (_isGzipped && _bamRuledOut) {
			return readZipChunk();
		}
		currChar = _pushBackStreamBuf->sbumpc();
		//Stop when EOF hit.
		if (currChar == EOF) {
			_eofHit = true;
			break;
		}
		numChars++;
		_scanBuffer.push_back(currChar);
		if (_isGzipped) {
			if (!_bamRuledOut && detectBamOrBgzip(numChars, currChar)) {
				//we now know the file is in bgzipped format
				return true;
			}
			if (numChars == 0) {
				continue; //this will only happen when we've just discovered that this
				//is definitely not BAM, and want to start over.
			}
		}

		//Stop if we have the minimum number of bytes and newline is hit.
		//For gzip, stop at SCAN_BUFFER_SIZE.
		if (currChar == '\n' && numChars >= MIN_SCAN_BUFFER_SIZE ){
			break;
		}
	}
	_numBytesInBuffer = _scanBuffer.size();

	//append it to the savedDataStr.
	_scanBuffer.toStr(_saveDataStr, true);
	if (_numBytesInBuffer == 0) return false;
	return true;

}

bool InputStreamMgr::detectBamOrBgzip(int &numChars, int currChar)
{
	//Look for the BAM magic string "BAM\1" in the first fouur characters of the input stream.
	//In compressed form, the first char is the gzip signifier, which was already found.
	//The next three are the integers 139, 8, and 4.
	if (numChars < 5) {
		_possibleBamCode[numChars -1] = currChar;
		//special: test for BAM
		if (numChars == 4 && _possibleBamCode[1] == 139 && _possibleBamCode[2] == 8 && _possibleBamCode[3] == 4) {
			//BAM magic string detected.This is either a BAM or bgzip file. To find out which, we have to try and
			//open the file as BAM, with a BAM reader, and see if the header and references are both non-empty.
			//However, if they are empty, we will have had to save all bytes consumed in the attempt, meaning still
			//fill the scanBuffer and push it back onto the pushBackStream as normal.
			for (; numChars < BAM_SCAN_BUFFER_SIZE; numChars++) {
				currChar = _pushBackStreamBuf->sbumpc();
				//Stop when EOF hit.
				if (currChar == EOF) {
					_eofHit = true;
					break;
				}
				_scanBuffer.push_back(currChar);

			}
			_pushBackStreamBuf->pushBack(_scanBuffer);

			//ok, now all the data read so far is saved in the scan buffer, and pushbackstream is reset.
			//now we make a BamReader and try to open the file.


			_bamReader = new BamTools::BamReader();
			if (!_bamReader->OpenStream(_finalInputStream)) {
				//This is NOT a bam file, but it is bgzipped.
				_pushBackStreamBuf->clear();
				//Put all bytes read so far back onto the scan buffer, then reset
				//everything so that we're effectively starting over.
				_pushBackStreamBuf->pushBack(_scanBuffer);
				_scanBuffer.clear();
				numChars = 0;
				_isBam = false;
				_isBgzipped = true;
				_bamRuledOut = true;
				_numBytesInBuffer = 0;
				delete _bamReader;
				_bamReader = NULL;

				//Alter the finalInputSream to become a bgzfReader.
				_bgStream = new BamTools::Internal::BgzfStream();
				_bgStream->OpenStream(_finalInputStream, BamTools::IBamIODevice::ReadOnly);

				return false;
			}
			//This is a BAM file.
			_isBam = true;
			_numBytesInBuffer = _scanBuffer.size();
			return true;
		} else if (numChars == 4) {
			//This is a gzipped file, and it is not bgzipped or BAM.
			_pushBackStreamBuf->clear();
			_pushBackStreamBuf->pushBack(_scanBuffer);
			_scanBuffer.clear();
			numChars = 0;
			_isBam = false;
			_isBgzipped = false;
			_bamRuledOut = true;
			_numBytesInBuffer = 0;
			_infStreamBuf = new InflateStreamBuf(_finalInputStream);
			delete _oldInputStream;
			_oldInputStream = _finalInputStream;
			_finalInputStream = new istream(_infStreamBuf);
			return false;
		}
	}
	return false;
}

bool InputStreamMgr::readZipChunk()
{
	if (_tmpZipBuf == NULL) {
		_tmpZipBuf = new char[SCAN_BUFFER_SIZE +1];
	}
	memset(_tmpZipBuf, 0, SCAN_BUFFER_SIZE +1);
	size_t numCharsRead = read(_tmpZipBuf, (size_t)SCAN_BUFFER_SIZE);
	_saveDataStr.append(_tmpZipBuf);
	_numBytesInBuffer = _saveDataStr.size();
	if ((int)numCharsRead < SCAN_BUFFER_SIZE) {
		_streamFinished = true;
	}
	if (numCharsRead == 0) return false;
	return true;
}

bool InputStreamMgr::resetStream()
{
	_saveDataStr.clear();
	if (!_isBam && !_isStdin && !_isGzipped) {
		//For non-compressed, non-stdin file input, just re-open the file.
		delete _finalInputStream;
		_finalInputStream = new ifstream(_filename.c_str());
		return true;
	}
	return false;
}


