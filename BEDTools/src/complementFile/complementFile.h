/*
 * complementFile.h
 *
 *  Created on: Jun 16, 2015
 *      Author: nek3d
 */

#ifndef COMPLEMENTFILE_H_
#define COMPLEMENTFILE_H_

#include "ToolBase.h"
#include "ContextComplement.h"


class BlockMgr;
class BinTree;
class NewGenomeFile;
class RecordOutputMgr;

class ComplementFile : public ToolBase {

public:
    ComplementFile(ContextComplement *context);
    virtual ~ComplementFile();
	virtual bool init();
	virtual bool findNext(RecordKeyVector &hits);
	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits);
	virtual void cleanupHits(RecordKeyVector &hits);
	virtual bool finalizeCalculations();
	virtual void  giveFinalReport(RecordOutputMgr *outputMgr) {}


protected:
	FileRecordMergeMgr *_frm;
	Bed3Interval _outRecord;
	QuickString _currChrom;
	const NewGenomeFile *_genomeFile;
	int _currStartPos;
	RecordOutputMgr *_outputMgr;
	const vector<QuickString> &_chromList;
	int _currPosInGenomeList;

	virtual ContextComplement *upCast(ContextBase *context) { return static_cast<ContextComplement *>(context); }

	void outPutLastRecordInPrevChrom();
	bool fastForward(const QuickString &newChrom);
	void printRecord(int endPos);

};




#endif /* COMPLEMENTFILE_H_ */
