/*
 * jaccard.cpp
 *
 *  Created on: Apr 24, 2015
 *      Author: nek3d
 */

#include "jaccard.h"

Jaccard::Jaccard(ContextJaccard *context)
: IntersectFile(context),
 _queryUnion(0),
 _dbUnion(0),
 _intersectionVal(0),
 _unionVal(0),
 _numIntersections(0)
{

}

bool Jaccard::findNext(RecordKeyVector &hits) {
	if (nextSortedFind(hits)) {
		checkSplits(hits);
		_intersectionVal += getTotalIntersection(hits);
		return true;
	}
	return false;
}

void Jaccard::cleanupHits(RecordKeyVector &hits)
{
//	upCastFRM(_queryFRM)->deleteMergedRecord(hits);
}


bool Jaccard::finalizeCalculations() {
	_sweep->closeOut();
	_queryUnion = _sweep->getQueryTotalRecordLength();
	_dbUnion = _sweep->getDatabaseTotalRecordLength();

	_unionVal = _queryUnion + _dbUnion;
	return true;
}

void  Jaccard::giveFinalReport(RecordOutputMgr *outputMgr) {
	// header
	outputMgr->checkForHeader();

	cout << "intersection\tunion-intersection\tjaccard\tn_intersections" << endl;

	unsigned long adjustedUnion = _unionVal - _intersectionVal;

	cout << _intersectionVal << "\t" << adjustedUnion << "\t" <<
			(float) _intersectionVal / (float)adjustedUnion << "\t" << _numIntersections << endl;
}

unsigned long Jaccard::getTotalIntersection(RecordKeyVector &hits)
{
	unsigned long intersection = 0;
	const Record *key = hits.getKey();
	int keyStart = key->getStartPos();
	int keyEnd = key->getEndPos();

	int hitIdx = 0;
	for (RecordKeyVector::const_iterator_type iter = hits.begin(); iter != hits.end(); iter = hits.next()) {
		const Record *currRec = *iter;
		int maxStart = max(currRec->getStartPos(), keyStart);
		int minEnd = min(currRec->getEndPos(), keyEnd);
		if (_context->getObeySplits()) {
			intersection += upCast(_context)->getSplitBlockInfo()->getOverlapBases(hitIdx);
			hitIdx++;
		} else {
			intersection += (unsigned long)(minEnd - maxStart);
		}
	}
	_numIntersections += (int)hits.size();
	return intersection;
}

