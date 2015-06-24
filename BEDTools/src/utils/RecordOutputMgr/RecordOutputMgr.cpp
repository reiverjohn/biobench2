/*
 * RecordOutputMgr.cpp
 *
 *  Created on: May 28, 2013
 *      Author: nek3d
 */

#include "RecordOutputMgr.h"
#include "ContextBase.h"
#include "ContextIntersect.h"
#include "ContextClosest.h"
#include "BlockMgr.h"
#include "Bed3Interval.h"
#include "Bed4Interval.h"
#include "BedGraphInterval.h"
#include "Bed5Interval.h"
#include "Bed6Interval.h"
#include "BedPlusInterval.h"
#include "Bed12Interval.h"
#include "BamRecord.h"
#include "VcfRecord.h"
#include "GffRecord.h"



#include <cstdio>


RecordOutputMgr::RecordOutputMgr()
: _context(NULL),
  _printable(true),
  _bamWriter(NULL),
  _currBamBlockList(NULL),
  _bamBlockMgr(NULL)
{
	_bamBlockMgr = new BlockMgr();
}

RecordOutputMgr::~RecordOutputMgr()
{
	// In the rare case when a file had a header but was otherwise empty,
	// we'll need to make a last check to see if the header still needs to be printed.
	checkForHeader();

	if (_outBuf.size() > 0) {
		flush();
	}
	if (_bamWriter != NULL) {
		_bamWriter->Close();
		delete _bamWriter;
		_bamWriter = NULL;
	}
	delete _bamBlockMgr;
	_bamBlockMgr = NULL;

}

void RecordOutputMgr::init(ContextBase *context) {
	_context = context;
	if (_context->getOutputFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
		//set-up BAM writer.
		_bamWriter = new BamTools::BamWriter();
		_bamWriter->SetCompressionMode(_context->getUncompressedBam() ?  BamTools::BamWriter::Uncompressed : BamTools::BamWriter::Compressed);

		int bamFileIdx = _context->getBamHeaderAndRefIdx();
		_bamWriter->Open("stdout", _context->getFile(bamFileIdx)->getHeader().c_str(), _context->getFile(bamFileIdx)->getBamReferences());
	} else {
		//for everything but BAM, we'll copy output to an output buffer before printing.
		_outBuf.reserve(MAX_OUTBUF_SIZE);
	}
	if (_context->getProgram() == ContextBase::INTERSECT) {
		if ((static_cast<ContextIntersect *>(_context))->getAnyHit() || (static_cast<ContextIntersect *>(_context))->getNoHit() ||
				(static_cast<ContextIntersect *>(_context))->getWriteCount()) {
			_printable = false;
		}
	}
}

bool RecordOutputMgr::printKeyAndTerminate(RecordKeyVector &keyList) {
	if (_context->getProgram() == ContextBase::MERGE) {
		//when printing merged records, we want to force the printing into
		//bed3 format, which is surprisingly difficult to do. Had to use the following:
		const Bed3Interval *bed3 = static_cast<const Bed3Interval *>(keyList.getKey());
		bed3->Bed3Interval::print(_outBuf);
		return false;
	}
	printBamType bamCode = printBamRecord(keyList);
	if (bamCode == BAM_AS_BAM) {
		return true;
	} else if (bamCode == NOT_BAM) {
		keyList.getKey()->print(_outBuf);
		return false;
	}
	//otherwise, it was BAM_AS_BED, and the key was printed.
	return false;

}

RecordOutputMgr::printBamType RecordOutputMgr::printBamRecord(RecordKeyVector &keyList, bool bamOutputOnly)
{
	const Record *record = keyList.getKey();
	if (record->getType() == FileRecordTypeChecker::BAM_RECORD_TYPE) {
		if (_context->getOutputFileType() == FileRecordTypeChecker::BAM_FILE_TYPE) {
			_bamWriter->SaveAlignment(static_cast<const BamRecord *>(record)->getAlignment());
			return BAM_AS_BAM;
		} else {
			if (!bamOutputOnly) {
				if (record->isUnmapped()) {
					record->printUnmapped(_outBuf);
				} else {
					static_cast<const BamRecord *>(record)->print(_outBuf, _currBamBlockList);
				}
			}
			return BAM_AS_BED;
		}
	}
	return NOT_BAM;
}

void RecordOutputMgr::printRecord(const Record *record)
{
	RecordKeyVector keyList(record);
	printRecord(keyList);
}

void RecordOutputMgr::printRecord(const Record *record, const QuickString & value)
{	
	_afterVal = value;
	bool recordPrinted = false;
	if (record != NULL) {
		printRecord(record);
		recordPrinted = true;
	}
	if (!value.empty()) {
		if (recordPrinted) tab();
		_outBuf.append(value);
	}
	newline();

	if (needsFlush()) flush();
}

void RecordOutputMgr::printClosest(RecordKeyVector &keyList, const vector<int> *dists) {

	//The first time we print a record is when we print any header, because the header
	//hasn't been read from the query file until after the first record has also been read.
	checkForHeader();

	const ContextClosest *context = static_cast<const ContextClosest *>(_context);
	bool deleteBlocks = false;
	const Record *keyRec = keyList.getKey();
	RecordKeyVector blockList(keyRec);
	if (keyRec->getType() == FileRecordTypeChecker::BAM_RECORD_TYPE) {
		_bamBlockMgr->getBlocks(blockList, deleteBlocks);
		_currBamBlockList = &blockList;
	}
	if (!keyList.empty()) {
		int distCount = 0;
		for (RecordKeyVector::const_iterator_type iter = keyList.begin(); iter != keyList.end(); iter = keyList.next()) {
			const Record *hitRec = *iter;
			printKey(keyRec, keyRec->getStartPosStr(), keyRec->getEndPosStr());
			tab();
			addDbFileId(hitRec->getFileIdx());
			printKey(hitRec, hitRec->getStartPosStr(), hitRec->getEndPosStr());
			if (dists != NULL) {
				tab();
				int dist = (*dists)[distCount];
				//if not using sign distance, use absolute value instead.
				dist = context->signDistance() ? dist : abs(dist);
				_outBuf.append(dist);
				distCount++;
			}
			newline();
			if (needsFlush()) flush();
		}
	} else {
		printKey(keyRec, keyRec->getStartPosStr(), keyRec->getEndPosStr());
		tab();
		null(false, true);
		if (context->reportDistance()) {
			tab();
			_outBuf.append(-1);
		}
		newline();
	}
	if (deleteBlocks) {
		_bamBlockMgr->deleteBlocks(blockList);
		_currBamBlockList = NULL;
	}
	return;
}


void RecordOutputMgr::printRecord(RecordKeyVector &keyList) {
	if (keyList.getKey()->getType() == FileRecordTypeChecker::BAM_RECORD_TYPE) {
		RecordKeyVector blockList(keyList.getKey());
		bool deleteBlocks = false;
		_bamBlockMgr->getBlocks(blockList, deleteBlocks);
		printRecord(keyList, &blockList);
		if (deleteBlocks) {
			_bamBlockMgr->deleteBlocks(blockList);
		}
		return;
	}
    printRecord(keyList, NULL);

}

void RecordOutputMgr::printRecord(RecordKeyVector &keyList, RecordKeyVector *blockList)
{
	if (needsFlush()) {
		flush();
	}

	//The first time we print a record is when we print any header, because the header
	//hasn't been read from the query file until after the first record has also been read.
	checkForHeader();

	const_cast<Record *>(keyList.getKey())->undoZeroLength();
	_currBamBlockList = blockList;

	if (_context->getProgram() == ContextBase::INTERSECT || _context->getProgram() == ContextBase::SUBTRACT) {
		if (_printable) {
			if (keyList.empty()) {
				if ((static_cast<ContextIntersect *>(_context))->getWriteAllOverlap()) {
					// -wao the user wants to force the reporting of 0 overlap
					if (printKeyAndTerminate(keyList)) {

						_currBamBlockList = NULL;
						const_cast<Record *>(keyList.getKey())->adjustZeroLength();

						return;
					}
					tab();
					null(false, true);
					tab();
					_outBuf.append('0');
					newline();
					if (needsFlush()) flush();
				}
				else if ((static_cast<ContextIntersect *>(_context))->getLeftJoin()) {
					if (printKeyAndTerminate(keyList)) {
						_currBamBlockList = NULL;

						const_cast<Record *>(keyList.getKey())->adjustZeroLength();
						return;
					}
					tab();
					null(false, true);
					newline();
					if (needsFlush()) flush();
					_currBamBlockList = NULL;

					return;
				}
			} else {
				if (printBamRecord(keyList, true) == BAM_AS_BAM) {
					_currBamBlockList = NULL;

					const_cast<Record *>(keyList.getKey())->adjustZeroLength();
					return;
				}
				int hitIdx = 0;
				for (RecordKeyVector::const_iterator_type iter = keyList.begin(); iter != keyList.end(); iter = keyList.next()) {
					reportOverlapDetail(keyList.getKey(), *iter, hitIdx);
					hitIdx++;
				}
			}
		} else { // not printable
			reportOverlapSummary(keyList);
		}
		_currBamBlockList = NULL;
	} else if (_context->getProgram() == ContextBase::SAMPLE) {
		if (!printKeyAndTerminate(keyList)) {
			newline();
		}
	} else { // if (_context->getProgram() == ContextBase::MAP || _context->getProgram() == ContextBase::MERGE) {
		printKeyAndTerminate(keyList);
	}
	_currBamBlockList = NULL;
	const_cast<Record *>(keyList.getKey())->adjustZeroLength();

}

void RecordOutputMgr::checkForHeader() {
	// Do we need to print a header?
	if (!_context->getPrintHeader()) return;

	//if the program is based on intersection, we want the header from the query file.
	if (_context->hasIntersectMethods()) {
		int queryIdx = (static_cast<ContextIntersect *>(_context))->getQueryFileIdx();
		const QuickString &header  = _context->getFile(queryIdx)->getHeader();
		_outBuf.append(header);
	} else {
		_outBuf.append(_context->getFile(0)->getHeader());
	}
	_context->setPrintHeader(false);
	flush();
}

void RecordOutputMgr::reportOverlapDetail(const Record *keyRecord, const Record *hitRecord, int hitIdx)
{
	//get the max start and min end as strings.
	const_cast<Record *>(hitRecord)->undoZeroLength();


	const QuickString *startStr = NULL;
	const QuickString *endStr = NULL;
	int maxStart = 0;
	int minEnd = 0;

	int keyStart = keyRecord->getStartPos();
	int keyEnd = keyRecord->getEndPos();
	int hitStart = hitRecord->getStartPos();
	int hitEnd = hitRecord->getEndPos();

	if (  keyStart>= hitStart) {
		//the key start is after the hit start, but we need to check and make sure the hit end is at least after the keyStart.
		//The reason for this is that, in some rare cases, such as both the key and hit having been zero length intervals,
		//the normal process for intersection that allows us to simply report the maxStart and minEnd do not necessarily apply.
		if (hitEnd >= keyStart) {
			//this is ok. We have a normal intersection where the key comes after the hit.

			maxStart = keyStart;
			startStr = &(keyRecord->getStartPosStr());

			minEnd = min(keyEnd, hitEnd);
			endStr = keyRecord->getEndPos() < hitRecord->getEndPos() ? &(keyRecord->getEndPosStr()) : &(hitRecord->getEndPosStr());

		} else {
			//this is the weird case of not a "real" intersection. The keyStart is greater than the hitEnd. So just report the key as is.
			maxStart = keyStart;
			minEnd = keyEnd;
			startStr = &(keyRecord->getStartPosStr());
			endStr = &(keyRecord->getEndPosStr());
		}

	} else {
		//all of the above, but backwards. keyStart is before hitStart.
		if (keyEnd >= hitStart) {
			//normal intersection, key first
			maxStart = hitStart;
			startStr = &(hitRecord->getStartPosStr());
			minEnd = min(keyEnd, hitEnd);
			endStr = keyRecord->getEndPos() < hitRecord->getEndPos() ? &(keyRecord->getEndPosStr()) : &(hitRecord->getEndPosStr());
		} else {
			//this is the weird case of not a "real" intersection. The hitStart is greater than the keyEnd. So just report the hit as is.
			maxStart = hitStart;
			minEnd = hitEnd;
			startStr = &(hitRecord->getStartPosStr());
			endStr = &(hitRecord->getEndPosStr());

		}
	}


	if (!(static_cast<ContextIntersect *>(_context))->getWriteA() && !(static_cast<ContextIntersect *>(_context))->getWriteB()
			&& !(static_cast<ContextIntersect *>(_context))->getWriteOverlap() && !(static_cast<ContextIntersect *>(_context))->getLeftJoin()) {
		printKey(keyRecord, *startStr, *endStr);
		newline();
		if (needsFlush()) flush();
	}
	else if (((static_cast<ContextIntersect *>(_context))->getWriteA() &&
			(static_cast<ContextIntersect *>(_context))->getWriteB()) || (static_cast<ContextIntersect *>(_context))->getLeftJoin()) {
		printKey(keyRecord);
		tab();
		addDbFileId(hitRecord->getFileIdx());
		hitRecord->print(_outBuf);
		newline();
		if (needsFlush()) flush();
	}
	else if ((static_cast<ContextIntersect *>(_context))->getWriteA()) {
		printKey(keyRecord);
		newline();
		if (needsFlush()) flush();
	}
	else if ((static_cast<ContextIntersect *>(_context))->getWriteB()) {
		printKey(keyRecord, *startStr, *endStr);
		tab();
		addDbFileId(hitRecord->getFileIdx());
		hitRecord->print(_outBuf);
		newline();
		if (needsFlush()) flush();
	}
	else if ((static_cast<ContextIntersect *>(_context))->getWriteOverlap()) {
		int printOverlapBases = 0;
		if (_context->getObeySplits()) {
			printOverlapBases = _context->getSplitBlockInfo()->getOverlapBases(hitIdx);
		} else {
			printOverlapBases = minEnd - maxStart;
		}
		printKey(keyRecord);
		tab();
		addDbFileId(hitRecord->getFileIdx());
		hitRecord->print(_outBuf);
		tab();
		int2str(printOverlapBases, _outBuf, true);
		newline();
		if (needsFlush()) flush();
	}
	const_cast<Record *>(hitRecord)->adjustZeroLength();
}

void RecordOutputMgr::reportOverlapSummary(RecordKeyVector &keyList)
{
	int numOverlapsFound = (int)keyList.size();
	if ((static_cast<ContextIntersect *>(_context))->getAnyHit() && numOverlapsFound > 0) {
		if (printKeyAndTerminate(keyList)) {
			return;
		}
		newline();
		if (needsFlush()) flush();
	} else if ((static_cast<ContextIntersect *>(_context))->getWriteCount()) {
		if (printKeyAndTerminate(keyList)) {
			return;
		}
		tab();
		int2str(numOverlapsFound, _outBuf, true);
		newline();
		if (needsFlush()) flush();
	} else if ((static_cast<ContextIntersect *>(_context))->getNoHit() && numOverlapsFound == 0) {
		if (printKeyAndTerminate(keyList)) {
			return;
		}
		newline();
		if (needsFlush()) flush();
	}
}

void RecordOutputMgr::addDbFileId(int fileId) {
	if ((static_cast<ContextIntersect *>(_context))->getNumDatabaseFiles()  == 1) return;
	if (!_context->getUseDBnameTags() && (!_context->getUseDBfileNames())) {
		_outBuf.append(fileId);
	} else if (_context->getUseDBnameTags()){
		_outBuf.append((static_cast<ContextIntersect *>(_context))->getDatabaseNameTag((static_cast<ContextIntersect *>(_context))->getDbIdx(fileId)));
	} else {
		_outBuf.append(_context->getInputFileName(fileId));
	}
	tab();
}

void RecordOutputMgr::null(bool queryType, bool dbType)
{
	FileRecordTypeChecker::RECORD_TYPE recordType = FileRecordTypeChecker::UNKNOWN_RECORD_TYPE;
	if (_context->hasIntersectMethods()) {
		if (queryType) {
			recordType = (static_cast<ContextIntersect *>(_context))->getQueryRecordType();
		} else if (dbType) {
			recordType = (static_cast<ContextIntersect *>(_context))->getDatabaseRecordType(0);
		}
	} else  {
		recordType = _context->getFile(0)->getRecordType();
	}
	//This is kind of a hack. Need an instance of the correct class of record in order to call it's printNull method.
	Record *dummyRecord = NULL;

	switch (recordType) {
	case FileRecordTypeChecker::BED3_RECORD_TYPE:
		dummyRecord = new Bed3Interval();
		break;
	case FileRecordTypeChecker::BED4_RECORD_TYPE:
		dummyRecord = new Bed4Interval();
		break;
	case FileRecordTypeChecker::BEDGRAPH_RECORD_TYPE:
		dummyRecord = new BedGraphInterval();
		break;
	case FileRecordTypeChecker::BED5_RECORD_TYPE:
		dummyRecord = new Bed5Interval();
		break;
	case FileRecordTypeChecker::BED6_RECORD_TYPE:
		dummyRecord = new Bed6Interval();
		break;
	case FileRecordTypeChecker::BED12_RECORD_TYPE:
		dummyRecord = new Bed12Interval();
		break;
	case FileRecordTypeChecker::BED_PLUS_RECORD_TYPE:
		dummyRecord = new BedPlusInterval();
		(static_cast<BedPlusInterval *>(dummyRecord))->setNumPrintFields((static_cast<ContextIntersect *>(_context))->getMaxNumDatabaseFields());
		break;
	case FileRecordTypeChecker::VCF_RECORD_TYPE:
		dummyRecord = new VcfRecord();
		(static_cast<VcfRecord *>(dummyRecord))->setNumPrintFields((static_cast<ContextIntersect *>(_context))->getMaxNumDatabaseFields());
		break;
	case FileRecordTypeChecker::BAM_RECORD_TYPE:
		dummyRecord = new BamRecord();
		break;
	case FileRecordTypeChecker::GFF_RECORD_TYPE:
		dummyRecord = new GffRecord();
		(static_cast<GffRecord *>(dummyRecord))->setNumFields((static_cast<ContextIntersect *>(_context))->getMaxNumDatabaseFields());
		break;
	default:
		break;
	}

	dummyRecord->printNull(_outBuf);
	delete dummyRecord;

}

void RecordOutputMgr::printKey(const Record *key, const QuickString & start, const QuickString & end)
{
	if (key->getType() != FileRecordTypeChecker::BAM_RECORD_TYPE) {
		key->print(_outBuf, start, end);
	} else {
		static_cast<const BamRecord *>(key)->print(_outBuf, start, end, _currBamBlockList);
	}
}

void RecordOutputMgr::printKey(const Record *key)
{
	if (key->getType() != FileRecordTypeChecker::BAM_RECORD_TYPE) {
		key->print(_outBuf);
	} else {
		static_cast<const BamRecord *>(key)->print(_outBuf, _currBamBlockList);
	}
}

void RecordOutputMgr::flush() {
	fwrite(_outBuf.c_str(), 1, _outBuf.size(), stdout);
	_outBuf.clear();
}
