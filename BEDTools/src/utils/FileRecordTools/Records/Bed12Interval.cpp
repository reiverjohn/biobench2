#include "Bed12Interval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstdlib>

Bed12Interval::Bed12Interval()
:_thickStart(-1),
 _thickEnd(-1),
 _blockCount(-1)

{
}

Bed12Interval::~Bed12Interval()
{
}
const Bed12Interval &Bed12Interval::operator=(const Bed12Interval &other) {
	Bed6Interval::operator=(other);

	_thickStart = other._thickStart;
	_thickEnd = other._thickEnd;
	_itemRGB = other._itemRGB;
	_blockCount = other._blockCount;
	_blockSizes = other._blockSizes;
	_blockStarts = other._blockStarts;

	return *this;
}

bool Bed12Interval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	bool baseRetFlag = Bed6Interval::initFromFile(fileReader);

	 fileReader->getField(6, _thickStartStr);
	 fileReader->getField(7, _thickEndStr);
	 fileReader->getField(8, _itemRGB);
	 fileReader->getField(9, _blockCountStr);
	 fileReader->getField(10, _blockSizes);
	 fileReader->getField(11, _blockStarts);

	 _thickStart = str2chrPos(_thickStartStr);
	 _thickEnd = str2chrPos(_thickEndStr);
	 _blockCount = str2chrPos(_blockCountStr);
	return baseRetFlag;
}

void Bed12Interval::clear() {
	Bed6Interval::clear();

	 _thickStart = -1;
	 _thickEnd = -1;
	 _itemRGB.clear();
	 _blockCount = -1;
	 _blockSizes.clear();
	 _blockStarts.clear();
	 _thickStartStr.clear();
	 _thickEndStr.clear();
	 _blockCountStr.clear();

}

void Bed12Interval::print(QuickString &outBuf) const
{
	Bed6Interval::print(outBuf);

	outBuf.append('\t');
	outBuf.append(_thickStartStr);
	outBuf.append('\t');
	outBuf.append(_thickEndStr);
	outBuf.append('\t');
	outBuf.append(_itemRGB);
	outBuf.append('\t');
	outBuf.append(_blockCountStr);
	outBuf.append('\t');
	outBuf.append(_blockSizes);
	outBuf.append('\t');
	outBuf.append(_blockStarts);
}

void Bed12Interval::print(QuickString &outBuf, int start, int end) const
{
	Bed6Interval::print(outBuf, start, end);

	outBuf.append('\t');
	outBuf.append(_thickStartStr);
	outBuf.append('\t');
	outBuf.append(_thickEndStr);
	outBuf.append('\t');
	outBuf.append(_itemRGB);
	outBuf.append('\t');
	outBuf.append(_blockCountStr);
	outBuf.append('\t');
	outBuf.append(_blockSizes);
	outBuf.append('\t');
	outBuf.append(_blockStarts);
}

void Bed12Interval::print(QuickString &outBuf, const QuickString & start, const QuickString & end) const
{
	Bed6Interval::print(outBuf, start, end);
	outBuf.append('\t');
	outBuf.append(_thickStartStr);
	outBuf.append('\t');
	outBuf.append(_thickEndStr);
	outBuf.append('\t');
	outBuf.append(_itemRGB);
	outBuf.append('\t');
	outBuf.append(_blockCountStr);
	outBuf.append('\t');
	outBuf.append(_blockSizes);
	outBuf.append('\t');
	outBuf.append(_blockStarts);
}


void Bed12Interval::printNull(QuickString &outBuf) const
{
	Bed6Interval::printNull(outBuf);

	outBuf.append("\t.\t.\t.\t.\t.\t.", 12);
}

const QuickString &Bed12Interval::getField(int fieldNum) const
{
	switch (fieldNum) {
	case 7:
		return _thickStartStr;
		break;
	case 8:
		return _thickEndStr;
		break;
	case 9:
		return _itemRGB;
		break;
	case 10:
		return _blockCountStr;
		break;
	case 11:
		return _blockSizes;
		break;
	case 12:
		return _blockStarts;
		break;
	default:
		return Bed6Interval::getField(fieldNum);
		break;
	}
}

bool Bed12Interval::isNumericField(int fieldNum) {
	switch (fieldNum) {
	case 7:
		return true;
		break;
	case 8:
		return true;
		break;
	case 9:
		return false;
		break;
	case 10:
		return true;
		break;
	case 11:
		return false;
		break;
	case 12:
		return false;
		break;
	default:
		return Bed6Interval::isNumericField(fieldNum);
		break;
	}
}

int Bed12Interval::getLength(bool obeySplits) const {
	//only bed12 and BAM need to check splits
	if (!obeySplits || _blockCount <=0) {
		return _endPos - _startPos;
	} else {
		int length = 0;
    	//parse the blockSizes string.
		char numBuf[16];
		const char *startPtr = _blockSizes.c_str();
		const char *endPtr = startPtr;
	    for (int i=0; i < _blockCount; i++) {
	    	memset(numBuf, 0, 16);
	    	endPtr = strchr(endPtr, ',');
	    	memcpy(numBuf, startPtr, endPtr - startPtr);
	    	length += str2chrPos(numBuf);
	    	startPtr = ++endPtr;
	    }
	    return length;
	}
}
