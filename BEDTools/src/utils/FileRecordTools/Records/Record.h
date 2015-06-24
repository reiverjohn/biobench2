/*
 * Record.h
 *
 *  Created on: Nov 8, 2012
 *      Author: nek3d
 */

#ifndef RECORD_H_
#define RECORD_H_

#include <string>
#include "FreeList.h"
#include "QuickString.h"
#include "FileRecordTypeChecker.h"

using namespace std;

class FileRecordMgr;
class FileReader;
class ChromIdLookup;

class Record {
public:
	friend class RecordMgr;
	friend class RecordOutputMgr;

	friend class FreeList<Record>;

	typedef enum { FORWARD, REVERSE, UNKNOWN } strandType;
	Record();
	virtual bool initFromFile(FileReader *) =0;
	virtual void clear();
	virtual void print(QuickString &outBuf) const {}
	virtual void print(QuickString &outBuf, int start, int end) const {}
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const {}
	virtual void print(FILE *fp, bool newline = false) const;
	virtual void printNull(QuickString &outBuf) const {}
	friend ostream &operator << (ostream &out, const Record &record);

	virtual const Record & operator=(const Record &);

	virtual const QuickString &getChrName() const { return _chrName; }
	virtual void setChrName(const QuickString &chr) { _chrName = chr; }
	virtual void setChrName(const string &chr) { _chrName = chr; }
	virtual void setChrName(const char *chr) { _chrName = chr; }

	virtual int getFileIdx() const { return _fileIdx; }
	virtual void setFileIdx(int fileIdx) { _fileIdx = fileIdx; }

	virtual int getChromId() const { return _chrId; }
	virtual void setChromId(int id) { _chrId = id; }

	virtual int getStartPos() const { return _startPos; }
	virtual void setStartPos(int startPos) { _startPos = startPos; }
	virtual const QuickString &getStartPosStr() const { return _startPosStr; }
	virtual void setStartPosStr(const QuickString &str) { _startPosStr = str; }

	virtual int getEndPos() const { return _endPos; }
	virtual void setEndPos(int endPos) { _endPos = endPos; }
	virtual const QuickString &getEndPosStr() const { return _endPosStr; }
	virtual void setEndPosStr(const QuickString &str) { _endPosStr = str; }

	virtual bool getZeroLength() const { return _zeroLength; }
	virtual void setZeroLength(bool val) { _zeroLength = val; }

	virtual const QuickString &getStrand() const { return _strand; }
	virtual void setStrand(const QuickString &val) { _strand = val;
		_strandVal = (val == "+" ? FORWARD : (val == "-" ? REVERSE : UNKNOWN));
	}
	virtual void setStrand(char val) { _strand = val;
		_strandVal = (val == '+' ? FORWARD : (val == '-' ? REVERSE : UNKNOWN));
	}
	virtual void adjustStrandVal() {
		_strandVal = (_strand == "+" ? FORWARD : (_strand == "-" ? REVERSE : UNKNOWN));
	}

	virtual strandType getStrandVal() const {return _strandVal; }

	virtual const QuickString &getName() const { return _name; }
	virtual void setName(const QuickString &chr) { _name = chr; }
	virtual void setName(const string &chr) { _name = chr; }
	virtual void setName(const char *chr) { _name = chr; }

	virtual const QuickString &getScore() const { return _score; }
	virtual void setScore(const QuickString &chr) { _score = chr; }
	virtual void setScore(const string &chr) { _score = chr; }
	virtual void setScore(const char *chr) { _score = chr; }

	virtual const QuickString &getField(int fieldNum) const;
	virtual int getNumFields() const  = 0;

	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::UNKNOWN_RECORD_TYPE; }

	virtual bool coordsValid(); //test that no coords negative, end not less than start, check zeroLength (see below).

	//Some files can have insertions of the form 2,2. If found this should translate to cover the base before and after,
	//thus meaning the startPos is decremented and the endPos is incremented. This method will find and handle that case.
	//Don't adjust the startPosStr and endPosStr strings because they aren't used in
	//calculation. They're only used in output, and it would be slower to change them
	//and then change them back.
	virtual void adjustZeroLength();
	virtual void undoZeroLength(); //change it back just before output;
	virtual bool isZeroLength() const { return _zeroLength; }

	// "Unmapped" only applies to BamRecord, but for design reasons, it has to be here,
	// because we want to short circuit the intersects method if either record is an unmapped
	// Bam record.
	bool isUnmapped() const { return _isUnmapped; }
	bool isMateUnmapped() const { return _isMateUnmapped; }
	virtual void printUnmapped(QuickString &outBuf) const {}



	virtual bool operator < (const Record &other) const;
	virtual bool operator > (const Record &other) const;
	virtual bool lessThan(const Record *other) const;
	virtual bool greaterThan(const Record *other) const;

	//is this on the same chromosome as another record?
	bool sameChrom(const Record *other) const;
	bool chromBefore(const Record *other) const;
	bool chromAfter(const Record *other) const;

	//is this record after the other one?
	virtual bool after(const Record *other) const;

	//does this record intersect with another record?
	virtual bool intersects(const Record *otherRecord,
			bool sameStrand, bool diffStrand, float overlapFraction, bool reciprocal) const;

	// *** WARNING !!! ** sameChromIntersects is a faster version of the intersects method,
	// BUT the caller MUST ensure that the records are on the same
	//chromosome. If you're not absolutely sure, use the regular intersects method.
	virtual bool sameChromIntersects(const Record *otherRecord,
			bool sameStrand, bool diffStrand, float overlapFraction, bool reciprocal) const;

//	virtual static bool isNumericField(int fieldNum) const = 0;

	bool hasChrInChromName() const;
	bool hasLeadingZeroInChromName(bool chrKnown = false) const;
	virtual int getLength(bool obeySplits) const;


protected:
	virtual ~Record(); //by making the destructor protected, only the friend class(es) can actually delete Record objects, or objects derived from Record.

	int _fileIdx; //associated file the record came from
	QuickString _chrName;
	int _chrId;
	int _startPos;
	int _endPos;
	//It is actually faster to also store the start and end positions as their original strings than to
	//have to convert their integer representations back to strings when printing them.
	QuickString _startPosStr;
	QuickString _endPosStr;
	QuickString _name;
	QuickString _score;
	QuickString _strand;
	strandType _strandVal;
	bool _zeroLength;
	bool _isUnmapped;
	bool _isMateUnmapped;
};

class RecordPtrSortAscFunctor {
public:
	bool operator()(const Record *rec1, const Record *rec2) const { return *rec1 < *rec2; }
};

class RecordPtrSortDescFunctor {
public:
	bool operator()(const Record *rec1, const Record *rec2) const { return *rec1 > *rec2; }
};
#endif /* RECORD_H_ */
