
#include "Record.h"
#include <cstdio>

Record::Record()
: _fileIdx(-1),
  _chrId(-1),
  _startPos(-1),
  _endPos(-1),
  _strandVal(UNKNOWN),
  _zeroLength(false),
  _isUnmapped(false),
  _isMateUnmapped(false)
{
}

Record::~Record() {
}

const Record &Record::operator=(const Record &other)
{
	_fileIdx = other._fileIdx;
	_chrName = other._chrName;
	_chrId = other._chrId;
	_startPos = other._startPos;
	_endPos = other._endPos;
	_strand = other._strand;
	_strandVal = other._strandVal;
	_name = other._name;
	return *this;
}

void Record::clear() {
	_fileIdx = -1;
	_chrName.clear();
	_chrId = -1;
	_startPos = -1;
	_endPos = -1;
	_name.clear();
	_score.clear();
	_strand.clear();
	_strandVal = UNKNOWN;
	_startPosStr.clear();
	_endPosStr.clear();
	_zeroLength = false;
	_isUnmapped = false;
	_isMateUnmapped = false;
}

bool Record::operator < (const Record &other) const
{

	if (!sameChrom(&other)) {
		return chromBefore(&other);
	}
	if (_startPos != other._startPos) {
		return _startPos < other._startPos;
	}
	if (_endPos != other._endPos) {
		return _endPos < other._endPos;
	}
	return false;
}

bool Record::operator > (const Record &other) const {
	if (!sameChrom(&other)) {
		return chromAfter(&other);
	}
	if (_startPos != other._startPos) {
		return _startPos > other._startPos;
	}
	if (_endPos != other._endPos) {
		return _endPos > other._endPos;
	}
	return false;
}

bool Record::lessThan(const Record *other) const
{
	if (!sameChrom(other)) {
		return chromBefore(other);
	}
	if (_startPos != other->_startPos) {
		return _startPos < other->_startPos;
	}
	if (_endPos != other->_endPos) {
		return _endPos < other->_endPos;
	}
	return false;
}

bool Record::greaterThan(const Record *other) const
{
	if (!sameChrom(other)) {
		return chromAfter(other);
	}
	if (_startPos != other->_startPos) {
		return _startPos > other->_startPos;
	}
	if (_endPos != other->_endPos) {
		return _endPos > other->_endPos;
	}
	return false;

}

bool Record::sameChrom(const Record *other) const {
	return (_chrId == -1 || other->_chrId == -1) ? ( _chrName == other->_chrName) : (_chrId == other->_chrId);
}

bool Record::chromBefore(const Record *other) const
{
	return (_chrId == -1 || other->_chrId == -1) ? ( _chrName < other->_chrName) : (_chrId < other->_chrId);
}

bool Record::chromAfter(const Record *other) const
{
	return (_chrId == -1 || other->_chrId == -1) ? ( _chrName > other->_chrName) : (_chrId > other->_chrId);
}



bool Record::after(const Record *other) const
{
	return (sameChrom(other) && _startPos >= other->_endPos);
}

bool Record::intersects(const Record *record,
                        bool sameStrand,
                        bool diffStrand,
                        float overlapFractionA,
                        float overlapFractionB,
                        bool reciprocalFraction,
                        bool eitherFraction) const
{
	//must be on same chromosome
	if (!sameChrom(record)) {
		return false;
	}
	return sameChromIntersects(record,
                               sameStrand, diffStrand,
                               overlapFractionA, overlapFractionB,
                               reciprocalFraction, eitherFraction);
}

bool Record::sameChromIntersects(const Record *record,
                                 bool sameStrand,
                                 bool diffStrand,
                                 float overlapFractionA,
                                 float overlapFractionB,
                                 bool reciprocalFraction,
                                 bool eitherFraction) const
{
	// Special: For records that are unmapped, intersect should automatically return false
	if (_isUnmapped || record->isUnmapped()) {
		return false;
	}

	//If user requested hits only on same strand, or only on different strands,
	//rule out different strandedness first.
	//If the strand is unknown in either case, then queries regarding strandedness
	//can not be answered, so we return false;
	bool isSameStrand = (_strandVal == record->_strandVal && _strandVal != UNKNOWN);
	bool isDiffStrand = ( _strandVal != UNKNOWN && record->_strandVal != UNKNOWN && _strandVal != record->_strandVal);

	if (sameStrand && !isSameStrand) {
		return false; //want same, but they're not same.
	}
	if (diffStrand && !isDiffStrand) {
		return false; //want different, but they're not different.
	}

	int otherStart = record->getStartPos();
	int otherEnd = record->getEndPos();

	bool otherZeroLen = (otherStart - otherEnd == 0);
	int maxStart = max(_startPos, otherStart);
	int minEnd = min(_endPos, otherEnd);

	bool localZeroLen = (_endPos - _startPos == 0);
	//rule out all cases of no intersection at all
	if (minEnd < maxStart) {
		return false;
	}


	if ((overlapFractionA == 0.0) && (overlapFractionB == 0.0))
    {
		//don't care about amount of overlap.
		//however, if minEnd and maxStart are equal, and
		//neither record is zeroLen, return false.
		if (minEnd == maxStart && !otherZeroLen && !localZeroLen) {
			return false;
		}
		return true;
	}

	int overlapBases = minEnd - maxStart;
	int aLen = _endPos - _startPos;
	int bLen = otherEnd - otherStart;

    float overlapA = (float)overlapBases / (float)aLen;
    float overlapB = (float)overlapBases / (float)bLen;

    bool sufficentFractionA = (overlapA >= overlapFractionA);
    bool sufficentFractionB = (overlapB >= overlapFractionB);

    if (!eitherFraction)
    {
        if (sufficentFractionA && sufficentFractionB) { return true; }
        return false;
    }
    else {
        if (sufficentFractionA || sufficentFractionB) { return true; }
        return false;
    }

	return false;
}

bool Record::coordsValid() {
	if (_startPos < 0 || _endPos < 0 || _endPos < _startPos) {
		return false;
	}
	adjustZeroLength();
	return true;
}

void Record::adjustZeroLength()
{
	if (_startPos == _endPos) {
		_zeroLength = true;
		_startPos--;
		_endPos++;
	}
}

void Record::undoZeroLength()
{
	if (_zeroLength) {
		_startPos++;
		_endPos--;
		_zeroLength = false;
	}
}

ostream &operator << (ostream &out, const Record &record)
{
	QuickString outBuf;
	record.print(outBuf);
	out << outBuf;
	return out;
}

const QuickString &Record::getField(int fieldNum) const
{
    cerr << endl << "*****" << endl
         << "*****ERROR: requested column " << fieldNum <<
         " , but record only has fields 1 - " << getNumFields() << ". Exiting." << endl
          << endl << "*****" << endl;
    exit(1);
}


bool Record::hasChrInChromName() const {
	const char *str = _chrName.c_str();
	//if the chrom name has at least 3 characters,
	//and the first 3 are c, h, r, case-insensitive,
	//return true. Otherwise, return false.
	return ((_chrName.size() >= 3) &&
			(str[0] == 'c' || str[0] == 'C') &&
			(str[1] == 'h' || str[1] == 'H') &&
			(str[2] == 'r' || str[2] == 'R'));
}

bool Record::hasLeadingZeroInChromName(bool chrKnown) const {
	const char *str = _chrName.c_str();
	// only check to see if the digit zero follows the occurance of
	// "chr", case insensitive.
	return (_chrName.size() >= 4 && str[3] == '0' && (chrKnown || hasChrInChromName()));
}

void Record::print(FILE *fp, bool newline) const {
	QuickString buf;
	print(buf);
	fprintf(fp, "%s", buf.c_str());
	if(newline) fprintf(fp, "\n");
}

int Record::getLength(bool obeySplits) const {
	//only bed12 and BAM need to check splits
	return _endPos - _startPos;
}
