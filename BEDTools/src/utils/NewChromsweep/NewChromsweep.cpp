/*****************************************************************************
  NewChromsweep.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include "NewChromsweep.h"
#include "ContextIntersect.h"
#include "FileRecordMgr.h"

NewChromSweep::NewChromSweep(ContextIntersect *context)
:	_context(context),
 	_queryFRM(NULL),
 	_numDBs(_context->getNumDatabaseFiles()),
 	_numFiles(_context->getNumInputFiles()),
 	_queryRecordsTotalLength(0),
 	_databaseRecordsTotalLength(0),
    _queryTotalRecords(0),
    _databaseTotalRecords(0),
 	_wasInitialized(false),
 	_currQueryRec(NULL),
 	_runToQueryEnd(false),
 	_lexicoDisproven(false),
 	_lexicoAssumed(false),
 	_lexicoAssumedFileIdx(-1),
 	_testLastQueryRec(false)
{
}


bool NewChromSweep::init() {
    
	//Create new FileRecordMgrs for the input files.
	//Open them, and get the first record from each.
	//otherwise, return true.
    _queryFRM = _context->getFile(_context->getQueryFileIdx());
    
    _dbFRMs.resize(_numDBs, NULL);
    for (int i=0; i < _numDBs; i++) {
    	_dbFRMs[i] = _context->getDatabaseFile(i);
    }

    _currDbRecs.resize(_numDBs, NULL);
    if (!_context->hasGenomeFile()) {
    	_fileTracks.resize(_numFiles, NULL);
    	for (int i=0; i < _numFiles; i++) {
    		_fileTracks[i] = new _orderTrackType;
    	}
    }


    for (int i=0; i < _numDBs; i++) {
    	nextRecord(false, i);
    	testChromOrder(_currDbRecs[i]);
    }

    _caches.resize(_numDBs);

    //determine whether to stop when the database end is hit, or keep going until the
    //end of the query file is hit as well.

    if (_context->getNoHit() || _context->getWriteCount() || _context->getWriteOverlap() || _context->getWriteAllOverlap() || _context->getLeftJoin()) {
    	_runToQueryEnd = true;
    }
    _wasInitialized = true;
    return true;
 }

void NewChromSweep::closeOut(bool testChromOrderVal) {
	if (_testLastQueryRec) {
		testChromOrder(_currQueryRec);
	}
	while (!_queryFRM->eof()) {
		nextRecord(true);
		testChromOrder(_currQueryRec);
	}
	if (testChromOrderVal) testChromOrder(_currQueryRec);

    for (int i=0; i < _numDBs; i++) {
    	while (!_dbFRMs[i]->eof()) {
    		if (testChromOrderVal) testChromOrder(_currDbRecs[i]);
    		nextRecord(false, i);
    	}
   		if (testChromOrderVal) testChromOrder(_currDbRecs[i]);
    }
}

NewChromSweep::~NewChromSweep(void) {
	if (!_wasInitialized) {
		return;
	}
	testThatAllDbChromsExistInQuery();

	_queryFRM->deleteRecord(_currQueryRec);
	_currQueryRec = NULL;


    for (int i=0; i < _numDBs; i++) {
    	_dbFRMs[i]->deleteRecord(_currDbRecs[i]);
    	_currDbRecs[i] = NULL;
    }

	_queryFRM->close();

   for (int i=0; i < _numDBs; i++) {
	   _dbFRMs[i]->close();
   }

	if (!_context->hasGenomeFile()) {
		for (int i=0; i < _numFiles; i++) {
			delete _fileTracks[i];
		}
	}
}


void NewChromSweep::scanCache(int dbIdx, RecordKeyVector &retList) {
	recListIterType cacheIter = _caches[dbIdx].begin();
    while (cacheIter != _caches[dbIdx].end())
    {
    	const Record *cacheRec = cacheIter->value();
        if (_currQueryRec->sameChrom(cacheRec) && !_currQueryRec->after(cacheRec)) {
            if (intersects(_currQueryRec, cacheRec)) {
                retList.push_back(cacheRec);
            } else if (cacheRec->after(_currQueryRec)) break; // cacheRec is after the query rec, stop scanning.
            cacheIter = _caches[dbIdx].next();
        }
        else {
            cacheIter = _caches[dbIdx].deleteCurrent();
    		_dbFRMs[dbIdx]->deleteRecord(cacheRec);
        }
    }
}

void NewChromSweep::clearCache(int dbIdx)
{
	//delete all objects pointed to by cache
	recListType &cache = _caches[dbIdx];
	for (recListIterType iter = cache.begin(); iter != cache.end(); iter = cache.next()) {
		_dbFRMs[dbIdx]->deleteRecord(iter->value());
	}
	cache.clear();
}

void NewChromSweep::masterScan(RecordKeyVector &retList) {

	for (int i=0; i < _numDBs; i++) {
		if (dbFinished(i) || chromChange(i, retList, true)) {
			continue;
		} else {

			// scan the database cache for hits
			scanCache(i, retList);
			//skip if we hit the end of the DB
			// advance the db until we are ahead of the query. update hits and cache as necessary
			while (_currDbRecs[i] != NULL &&
					_currQueryRec->sameChrom(_currDbRecs[i]) &&
					!(_currDbRecs[i]->after(_currQueryRec))) {
				if (intersects(_currQueryRec, _currDbRecs[i])) {
					retList.push_back(_currDbRecs[i]);
				}
				if (_currQueryRec->after(_currDbRecs[i])) {
					_dbFRMs[i]->deleteRecord(_currDbRecs[i]);
					_currDbRecs[i] = NULL;
				} else {
					_caches[i].push_back(_currDbRecs[i]);
					_currDbRecs[i] = NULL;
				}
				nextRecord(false, i);
			}
		}
	}
}

bool NewChromSweep::chromChange(int dbIdx, RecordKeyVector &retList, bool wantScan)
{
	const Record *dbRec = _currDbRecs[dbIdx];

	if (_currQueryRec != NULL && _currQueryChromName != _prevQueryChromName) {
		_context->testNameConventions(_currQueryRec);
		testChromOrder(_currQueryRec);
	}

	if (dbRec != NULL) {
		_context->testNameConventions(dbRec);
		testChromOrder(dbRec);
	}

	// If the query rec and db rec are on the same chrom, stop.
	if (dbRec != NULL && _currQueryRec != NULL && _currQueryRec->sameChrom(dbRec)) return false;


	if (dbRec == NULL || _currQueryRec == NULL) return false;

	if (queryChromAfterDbRec(dbRec)) {
		// the query is ahead of the database. fast-forward the database to catch-up.
		QuickString oldDbChrom(dbRec->getChrName());
		while (dbRec != NULL &&
				queryChromAfterDbRec(dbRec)) {
				_dbFRMs[dbIdx]->deleteRecord(dbRec);
			if (!nextRecord(false, dbIdx)) break;
			dbRec =  _currDbRecs[dbIdx];
			const QuickString &newDbChrom = dbRec->getChrName();
			if (newDbChrom != oldDbChrom) {
				testChromOrder(dbRec);
				oldDbChrom = newDbChrom;
			}
		}
		clearCache(dbIdx);
        return false;
    } else {
        // the database is ahead of the query.
        // scan the cache for remaining hits on the query's current chrom.
    	if (wantScan) scanCache(dbIdx, retList);
        return true;
    }
}


bool NewChromSweep::next(RecordKeyVector &retList) {
	retList.clearVector();


	//make sure the first read of the query file is tested for chrom sort order.
	bool needTestSortOrder = false;
	if (_currQueryRec != NULL) {
		_queryFRM->deleteRecord(_currQueryRec);
	} else {
		needTestSortOrder = true;
	}

	if (!nextRecord(true)) return false; // query EOF hit
	retList.setKey(_currQueryRec);

	if (needTestSortOrder) testChromOrder(_currQueryRec);

	if (allCurrDBrecsNull() && allCachesEmpty() && !_runToQueryEnd) {
		_testLastQueryRec = true;
		return false;
	}
	_currQueryChromName = _currQueryRec->getChrName();

	masterScan(retList);

	if (_context->getSortOutput()) {
		retList.sortVector();
	}


	_prevQueryChromName = _currQueryChromName;
	return true;
}

bool NewChromSweep::nextRecord(bool query, int dbIdx) {
	if (query) {
		_currQueryRec = _queryFRM->getNextRecord();
		if (_currQueryRec != NULL) {
			_queryRecordsTotalLength += (unsigned long)(_currQueryRec->getLength(_context->getObeySplits()));
            _queryTotalRecords++;
			return true;
		}
		return false;
	} else { //database
		Record *rec = _dbFRMs[dbIdx]->getNextRecord();
		_currDbRecs[dbIdx] = rec;
		if (rec != NULL) {
			_databaseRecordsTotalLength += (unsigned long)(rec->getLength(_context->getObeySplits()));
			_databaseTotalRecords++;
			return true;
		}
		return false;
	}
}

bool NewChromSweep::intersects(const Record *rec1, const Record *rec2) const
{
	return rec1->sameChromIntersects(rec2, _context->getSameStrand(), _context->getDiffStrand(),
			_context->getOverlapFraction(), _context->getReciprocal());
}


bool NewChromSweep::allCachesEmpty() {
	for (int i=0; i < _numDBs; i++) {
		if (!_caches[i].empty()) {
			return false;
		}
	}
	return true;
}

bool NewChromSweep::allCurrDBrecsNull() {
	for (int i=0; i < _numDBs; i++) {
		if (_currDbRecs[i] != NULL) {
			return false;
		}
	}
	return true;
}

bool NewChromSweep::dbFinished(int dbIdx) {
	if (_currDbRecs[dbIdx] == NULL && _caches[dbIdx].empty()) {
		return true;
	}
	return false;
}

void NewChromSweep::testChromOrder(const Record *rec)
{
	// Only use this method if we don't have a genome file
	// and the record is valid
	if (_context->hasGenomeFile() || rec == NULL) return;

	int fileIdx = rec->getFileIdx();

	const QuickString &chrom = rec->getChrName();

	findChromOrder(rec);

	//determine what the previous chrom was for this file.
	map<int, QuickString>::iterator prevIter = _filePrevChrom.find(fileIdx);
	if (prevIter == _filePrevChrom.end()) {
		_filePrevChrom[fileIdx] = chrom;
		return; //no previously stored chrom for this file.
	}
	const QuickString prevChrom(prevIter->second);
	_filePrevChrom[fileIdx] = chrom;

	if (chrom == prevChrom) return;

	if (verifyChromOrderMismatch(chrom, prevChrom, fileIdx)) {
		fprintf(stderr, "ERROR: chromomsome sort ordering for file %s is inconsistent with other files. Record was:\n", _context->getInputFileName(fileIdx).c_str());
		rec->print(stderr, true);
		exit(1);
	}

	if (!_lexicoDisproven && chrom < prevChrom) {
		if (_lexicoAssumed) {
			// ERROR.
			fprintf(stderr, "ERROR: Sort order was unspecified, and file %s is not sorted lexicographically.\n",
					_context->getInputFileName(fileIdx).c_str());
			fprintf(stderr, "       Please re-reun with the -g option for a genome file.\n       See documentation for details.\n");
			exit(1);
		}
		_lexicoDisproven = true;
	}
}

bool NewChromSweep::queryChromAfterDbRec(const Record *dbRec)
{
	//If using a genome file, compare chrom ids.
	//Otherwise, compare global order, inserting as needed.
	if (_context->hasGenomeFile()) {
		return (_currQueryRec->getChromId() > dbRec->getChromId()) ;
	}
	//see if query has both
	const QuickString &qChrom = _currQueryRec->getChrName();
	const QuickString &dbChrom = dbRec->getChrName();
	const _orderTrackType *track = _fileTracks[_currQueryRec->getFileIdx()];
	_orderTrackType::const_iterator iter = track->find(qChrom);


	int qOrder = iter->second;
	iter = track->find(dbChrom);
	if (iter == track->end()) {
		//query file does not contain the dbChrom.
		//try a lexicographical comparison, if possible.
		return testLexicoQueryAfterDb(_currQueryRec, dbRec);
	}
	int dbOrder = iter->second;

	return (qOrder > dbOrder);
}


int NewChromSweep::findChromOrder(const Record *rec) {
	const QuickString &chrom = rec->getChrName();
	int fileIdx = rec->getFileIdx();
	_orderTrackType *track = _fileTracks[fileIdx];

	_orderTrackType::const_iterator iter = track->find(chrom);
	if (iter == track->end()) {
		//chrom never seen before. Enter into map.
		int val = (int)track->size();
		track->insert(pair<QuickString, int>(chrom, val));
		return val;
	} else {
		return iter->second;
	}
}

bool NewChromSweep::verifyChromOrderMismatch(const QuickString & chrom, const QuickString &prevChrom, int skipFile) {
	//for every file except the one being checked,
	//find the current and previous chrom. If a given file
	//is missing either, skip it and go on. If it has both,
	//and the curr has a lower order num than the prev, return true.

	//if that never happens, return false.
	for (int i=0; i < _numFiles; i++) {
		if (i == skipFile) continue;
		const _orderTrackType *track = _fileTracks[i];
		_orderTrackType::const_iterator iter = track->find(chrom);
		if (iter == track->end()) continue; //this file does not contain the curr chrom
		int currOrder = iter->second;
		iter = track->find(prevChrom);
		if (iter == track->end()) continue; //this file does not contain the prevChrom.
		int prevOrder = iter->second;

		if (currOrder < prevOrder) return true;
	}
	return false;
}

void NewChromSweep::testThatAllDbChromsExistInQuery()
{
	if (_context->hasGenomeFile() || !_lexicoDisproven) return;

	int queryIdx = _context->getQueryFileIdx();
	//get the query file track. Then check that every chrom in every db exists in it.
	const _orderTrackType *qTrack = _fileTracks[queryIdx];

	for (int i=0; i < _numFiles; i++) {
		if (i == queryIdx) continue;
		const _orderTrackType *dbTrack = _fileTracks[i];
		for (_orderTrackType::const_iterator iter = dbTrack->begin(); iter != dbTrack->end(); iter++) {
			const QuickString &chrom = iter->first;
			if (qTrack->find(chrom) == qTrack->end()) {
				fprintf(stderr, "ERROR: Database file %s contains chromosome %s, but the query file does not.\n",
						_context->getInputFileName(i).c_str(), chrom.c_str());
				fprintf(stderr, "       Please re-reun with the -g option for a genome file.\n       See documentation for details.\n");
				exit(1);
			}
		}
	}
}


bool NewChromSweep::testLexicoQueryAfterDb(const Record *queryRec, const Record *dbRec)
{
	if (_lexicoDisproven) return false;

	bool queryGreater = queryRec->getChrName() > dbRec->getChrName();
	if (!_lexicoAssumed && queryGreater) {
		_lexicoAssumed = true;
		_lexicoAssumedFileIdx = dbRec->getFileIdx();
		_lexicoAssumedChromName = dbRec->getChrName();
	}
	return queryGreater;
}
