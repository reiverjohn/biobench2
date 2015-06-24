/*
 * coverageFile.h
 *
 *  Created on: May 8, 2015
 *      Author: nek3d
 */

#ifndef COVERAGEFILE_H_
#define COVERAGEFILE_H_

#include "intersectFile.h"
#include "ContextCoverage.h"

class CoverageFile : public IntersectFile {
public:
	CoverageFile(ContextCoverage *);
	~CoverageFile();
	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits);
	virtual void cleanupHits(RecordKeyVector &hits);
	virtual void  giveFinalReport(RecordOutputMgr *outputMgr);


protected:
	QuickString _finalOutput;

	size_t *_depthArray;
	size_t _depthArrayCapacity;
	size_t _queryLen;
	size_t _totalQueryLen;
	int _queryOffset;
	static const int DEFAULT_DEPTH_CAPACITY = 1024;
	char *_floatValBuf;
	static const int floatValBufLen = 16;

	typedef map<size_t, size_t> depthMapType;
	depthMapType _currDepthMap;
	depthMapType _finalDepthMap;

	virtual ContextCoverage *upCast(ContextBase *context) { return static_cast<ContextCoverage*>(context); }
	void makeDepthCount(RecordKeyVector &hits);

	size_t countBasesAtDepth(size_t depth);

	void doCounts(RecordOutputMgr *outputMgr, RecordKeyVector &hits);
	void doPerBase(RecordOutputMgr *outputMgr, RecordKeyVector &hits);
	void doHist(RecordOutputMgr *outputMgr, RecordKeyVector &hits);
	void doDefault(RecordOutputMgr *outputMgr, RecordKeyVector &hits);

	void format(float val);

};


#endif /* COVERAGEFILE_H_ */
