/*
 * BedPlusInterval.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef BEDPLUSINTERVAL_H_
#define BEDPLUSINTERVAL_H_

#include "Bed3Interval.h"
#include "PlusFields.h"

class SingleLineDelimTextFileReader;

class BedPlusInterval : public Bed3Interval {
public:
	friend class FreeList<BedPlusInterval>;

	BedPlusInterval();
	virtual ~BedPlusInterval() {}
	void setNumFixedFields(int numFields);
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void clear();
	virtual void print(QuickString &outBuf) const;
	virtual void print(QuickString &outBuf, int start, int end) const;
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED_PLUS_RECORD_TYPE; }

	virtual const QuickString &getField(int fieldNum) const;
	virtual int getNumFields() const  { return _numFixedFields + _plusFields.size(); }

	virtual void setNumPrintFields(int num) { _numPrintFields = num; }
	virtual int getNumPrintFields() const { return _numPrintFields; }
	static bool isNumericField(int fieldNum);


protected:
	int _numFixedFields; //first fields have names, and are not stored in otherIdxs.
	static const int defaultNumFixedFields = 3;
	PlusFields _plusFields;
	int _numPrintFields;

	void printBed6PlusFields(QuickString &outBuf) const;
	void printBed6PlusNullFields(QuickString &outBuf) const;

};



#endif /* BEDPLUSINTERVAL_H_ */
