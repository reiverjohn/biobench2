/*
 * CloseSweep.h
 *
 *  Created on: Sep 25, 2014
 *      Author: nek3d
 */

#ifndef CLOSESWEEP_H_
#define CLOSESWEEP_H_

#include "NewChromsweep.h"
#include <list>
#include <set>

class ContextClosest;

class distanceTuple {
public:
	distanceTuple() : _dist(0), _rec(NULL), _isNeg(false) {}
	distanceTuple(int dist, const Record *rec, bool isNeg = false) : _dist(dist), _rec(rec), _isNeg(isNeg) {}
//	bool operator < (const distanceTuple & other) const { return (_dist < other._dist); }
	int _dist;
	const Record *_rec;
	bool _isNeg;
};

class DistanceTupleSortAscFunctor {
public:
	bool operator()(const distanceTuple & d1, const distanceTuple & d2) const {
//		return ((d1._dist < d2._dist) ? true : (d1._dist == d2._dist ? (d1._rec->lessThan(d2._rec)) :  false)); }

		return (d1._dist < d2._dist ? true : (d1._dist == d2._dist ? d1._rec->lessThan(d2._rec) : false));
//		if (d1._dist < d2._dist) {
//			return true;
//		} else if (d1._dist == d2._dist) {
//			if () {
//				return true;
//			}
//		}
//		return false;
	}
};


class RecDistList {
public:
    typedef enum { LEFT, OVERLAP, RIGHT } chromDirType;
	RecDistList(int maxSize);
	~RecDistList();
	bool empty() const { return _empty; }
	void clear();
	int uniqueSize() const { return _currNumIdxs; }
	size_t totalSize() const { return _totalRecs; }
	bool addRec(int dist, const Record *, chromDirType chromDir);
	bool exists(int dist) const {
		int dummyVal = 0;
		return find(dist, dummyVal);
	}
	typedef pair<chromDirType, const Record *> elemPairType;
	typedef vector<elemPairType>elemsType;
	typedef pair<int, int> indexType;

	int getMaxDist() const { return _empty ? 0 : _distIndex[_currNumIdxs-1].first; }
	typedef int constIterType; //used to be a map iter, trying not to change interface too much.
	constIterType begin() const { return 0; }
	constIterType end() const { return _currNumIdxs; }
	int currDist(constIterType iter) const { return _distIndex[iter].first; }
	size_t currNumElems(constIterType iter) const { return allElems(iter)->size(); }
	const elemsType *allElems(constIterType iter) const { return _allRecs[_distIndex[iter].second]; }
	int getMaxLeftEndPos() const;

private:

	void insert(int dist, const Record *, chromDirType chromDir);


	//if true, pos will be the idx the distance is at.
	//if false, pos will be the idx to insert at.
	bool find(int dist, int &pos) const;


	int _kVal; //max unique allowed
	bool _empty;
	int _currNumIdxs;
	int _totalRecs;

	vector<elemsType *> _allRecs;
	indexType * _distIndex;
};

class CloseSweep : public NewChromSweep {
public:
	CloseSweep(ContextClosest *context);
	~CloseSweep(void);
	bool init();
	const vector<int> &getDistances() { return _finalDistances; }

private:
   ContextClosest *_context;
   int _kClosest; // how many closest hits we want to each query.
	vector<RecDistList *> _minUpstreamRecs;
	vector<RecDistList *> _minDownstreamRecs;
	vector<RecDistList *> _overlapRecs;
	vector<int> _maxPrevLeftClosestEndPos;
	vector<int> _maxPrevLeftClosestEndPosReverse;

	vector<int> _finalDistances;


	//structs to help with finding closest among all of multiple dbs.
	RecordKeyVector _copyRetList;
	vector<int> _copyDists;

	//override these methods from chromsweep
	void masterScan(RecordKeyVector &retList);
    void scanCache(int dbIdx, RecordKeyVector &retList);
    bool chromChange(int dbIdx, RecordKeyVector &retList, bool wantScan);

 	typedef enum { IGNORE, DELETE } rateOvlpType;
    rateOvlpType considerRecord(const Record *cacheRec, int dbIdx, bool &stopScanning);
    void finalizeSelections(int dbIdx, RecordKeyVector &retList);
    void checkMultiDbs(RecordKeyVector &retList);

    typedef enum { LEFT, OVERLAP, RIGHT } chromDirType;
    typedef enum { UPSTREAM, INTERSECT, DOWNSTREAM } streamDirType;

    void setLeftClosestEndPos(int dbIdx);
    bool beforeLeftClosestEndPos(int dbIdx, const Record *rec);
    void clearClosestEndPos(int dbIdx);
    bool canStopScan(const Record *cacheRec, bool ignored, streamDirType streamDir);
    int addRecsToRetList(const RecDistList::elemsType *recs, int currDist, RecordKeyVector &retList);
    void addSingleRec(const Record *rec, int currDist, int &hitsUsed, RecordKeyVector &retList);
    rateOvlpType tryToAddRecord(const Record *cacheRec, int dist, int dbIdx, bool &stopScanning, chromDirType chromDir, streamDirType streamDir);
    bool purgePointException(int dbIdx);

};


#endif /* CLOSESWEEP_H_ */
