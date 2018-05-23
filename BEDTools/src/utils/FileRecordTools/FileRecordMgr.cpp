
#include "FileRecordMgr.h"
#include "FreeList.h"
#include "Record.h"
#include "NewGenomeFile.h"

FileRecordMgr::FileRecordMgr(const QuickString &filename)
: _fileIdx(-1),
  _filename(filename),
  _bufStreamMgr(NULL),
  _fileReader(NULL),
  _fileType(FileRecordTypeChecker::UNKNOWN_FILE_TYPE),
  _recordType(FileRecordTypeChecker::UNKNOWN_RECORD_TYPE),
  _recordMgr(NULL),
  _isSortedInput(false),
  _freeListBlockSize(512),
  _useFullBamTags(false),
  _prevStart(INT_MAX),
  _prevChromId(-1),
  _mustBeForward(false),
  _mustBeReverse(false),
  _totalRecordLength(0),
  _totalMergedRecordLength(0),
  _blockMgr(NULL),
  _bamReader(NULL),
  _hasGenomeFile(false),
  _genomeFile(NULL),
  _ioBufSize(0),
  _noEnforceCoordSort(false)
 {
}

FileRecordMgr::~FileRecordMgr(){

	delete _bufStreamMgr;
	_bufStreamMgr = NULL;

	close(); //just make sure file was closed.
	delete _fileReader;
	_fileReader = NULL;

	delete _recordMgr;
	_recordMgr = NULL;
}

bool FileRecordMgr::open(bool inheader){
	_bufStreamMgr = new BufferedStreamMgr(_filename);
	_bufStreamMgr->getTypeChecker().setInHeader(inheader);

	if (_ioBufSize > 0) _bufStreamMgr->setIoBufSize(_ioBufSize);
	if (!_bufStreamMgr->init()) {
		cerr << "Error: unable to open file or unable to determine types for file " << _filename << endl;
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
		exit(1);
	}

	_fileType = _bufStreamMgr->getTypeChecker().getFileType();
	_recordType = _bufStreamMgr->getTypeChecker().getRecordType();
	if (_fileType == FileRecordTypeChecker::UNKNOWN_FILE_TYPE || _recordType == FileRecordTypeChecker::UNKNOWN_RECORD_TYPE) {
		cerr << "Error: Unable to determine type for file " << _filename << endl;
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
		exit(1);
	}
	allocateFileReader(inheader);
	_recordMgr = new RecordMgr(_recordType, _freeListBlockSize);

	_fileReader->setFileName(_filename.c_str());
	_fileReader->setInputStream(_bufStreamMgr);
	if (!_fileReader->open()) {
		cerr << "Error: Types determined but can't open file " << _filename << endl;
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
		exit(1);
	}

	return true;
}

void FileRecordMgr::close(){
	delete _bufStreamMgr;
	_bufStreamMgr = NULL;

	if (_fileReader != NULL) {
		_fileReader->close();
		delete _fileReader;
		_fileReader = NULL;
	}
}

bool FileRecordMgr::eof(){
	return _fileReader->eof();
}

Record *FileRecordMgr::getNextRecord(RecordKeyVector *keyList)
{
	if (!_fileReader->isOpen()) {
		return NULL;
	}
	if (!_fileReader->readEntry()) {
		return NULL;
	}
	Record *record = NULL;
	record = _recordMgr->allocateRecord();
	if (!record->initFromFile(_fileReader)) {
		_recordMgr->deleteRecord(record);
		return NULL;
	}

	// If the record is unmapped, don't test for valid coords or sort order,
	// but still return it so the -v (noHit) option and the like will still
	// see it.

	if (!record->isUnmapped() ) {
		if (!record->coordsValid() && (record->getType() != FileRecordTypeChecker::NO_POS_PLUS_RECORD_TYPE)) {
			cerr << "Error: Invalid record in file " << _filename << ". Record is " << endl << *record << endl;
			exit(1);
		}

		//test for sorted order, if necessary.
		if (_isSortedInput && !_noEnforceCoordSort) {
			testInputSortOrder(record);
		}
	}
	assignChromId(record);
	_totalRecordLength += (unsigned long)(record->getEndPos() - record->getStartPos());
	if (keyList != NULL) {
		keyList->setKey(record);
	}
	return record;
}

void FileRecordMgr::assignChromId(Record *record) {
	const QuickString &currChrom = record->getChrName();
	if (currChrom != _prevChrom  && _hasGenomeFile) {
		_prevChromId = _genomeFile->getChromId(currChrom);
		record->setChromId(_prevChromId);
	} else {
		record->setChromId(_prevChromId);
	}
}

void FileRecordMgr::testInputSortOrder(Record *record)
{

	// Special: For BAM records that aren't mapped, we actually don't want
	// to test the sort order. Another ugly hack sponsored by the letters B, A, and M.
	if (record->isUnmapped()) {
		return;
	}


	const QuickString &currChrom = record->getChrName();
	int currStart = record->getStartPos();
	if (record->isZeroLength()) {
		currStart++;
	}
	if (currChrom != _prevChrom) {
		if ( _foundChroms.find(currChrom) != _foundChroms.end()) {
			//this is a different chrom than the last record had, but we've already seen this chrom.
			sortError(record, false);
		} else {
			//new chrom has not been seen before.
			//TBD: test genome file for ChromId.
			if (_hasGenomeFile) {
				int currChromId = _genomeFile->getChromId(currChrom);
				if (currChromId < _prevChromId) {
					sortError(record, true);
				} else {
					_prevChromId = currChromId;
				}
			}
			_foundChroms.insert(currChrom);
			_prevChrom = currChrom;
			_prevStart = INT_MAX;
			record->setChromId(_prevChromId);
		}
	} else if (currStart < _prevStart) { //same chrom as last record, but with lower startPos, so still out of order.
		sortError(record, false);
	}
	_prevStart = currStart;

}

void FileRecordMgr::sortError(const Record *record, bool genomeFileError)
{
	if (genomeFileError) {
		cerr << "Error: Sorted input specified, but the file " << _filename << " has the following record with a different sort order than the genomeFile " <<
				_genomeFile->getGenomeFileName() << endl;
	} else {
		cerr << "Error: Sorted input specified, but the file " << _filename << " has the following out of order record" << endl;
	}
	cerr << *record << endl;
	exit(1);
}


void FileRecordMgr::deleteRecord(const Record *record) {
	_recordMgr->deleteRecord(record);
}

void FileRecordMgr::deleteRecord(RecordKeyVector *keyList) {
	_recordMgr->deleteRecord(keyList->getKey());
}

void FileRecordMgr::allocateFileReader(bool inheader)
{
	switch (_fileType) {
	case FileRecordTypeChecker::EMPTY_FILE_TYPE:
	case FileRecordTypeChecker::SINGLE_LINE_DELIM_TEXT_FILE_TYPE:
	case FileRecordTypeChecker::VCF_FILE_TYPE:
		_fileReader = new SingleLineDelimTextFileReader(_bufStreamMgr->getTypeChecker().getNumFields(), _bufStreamMgr->getTypeChecker().getDelimChar());
		static_cast<SingleLineDelimTextFileReader *>(_fileReader)->setInHeader(inheader);
		break;

	case FileRecordTypeChecker::BAM_FILE_TYPE:
		_fileReader = new BamFileReader();
		(static_cast<BamFileReader *>(_fileReader))->setUseTags(_useFullBamTags);
		(static_cast<BamFileReader *>(_fileReader))->setBamReader(_bufStreamMgr->getBamReader());
		break;
	default:
		break;
	}
	_fileReader->setFileIdx(_fileIdx);
}

const BamTools::RefVector & FileRecordMgr::getBamReferences() {
	// exta safety check to insure user checked the file type first.
	if (_fileType != FileRecordTypeChecker::BAM_FILE_TYPE) {
		cerr << "Error: Attempted to get BAM references from file " << _filename << ", which is NOT a BAM file." << endl;
		exit(1);
	}
	return static_cast<BamFileReader *>(_fileReader)->getReferences();
}
