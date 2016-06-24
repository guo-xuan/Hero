/*
 * SubjectRead.h
 *
 *  Created on: Oct 16, 2015
 *      Author: xgo
 */

#ifndef SUBJECTREAD_H_
#define SUBJECTREAD_H_

#include "Config.h"
#include "AlignmentRecord.h"
#include "QueryDataset.h"
using namespace std;

class PairKeyPosition;

class SubjectRead {
public:
	string sReadSequence;
	string sReadSequenceReverseComplement;
	INT64 iReadUniqueId;
	string sReadName;
	bool bIsGoodRead;
	vector<AlignmentRecord * > vAlignment;
	vector<PairKeyPosition * > vKeyPosition;
	string sAlignment;

	SubjectRead();
	~SubjectRead();

	void addAlignment(AlignmentRecord * _record);
	void clearKeyPosition();
	void clearAlignment();
	bool getAlignment(stringstream & _sstream);
	bool getAlignmentLongFormat(stringstream & _sstream);
	INT64 getReadUniqueId();
	INT64 getReadLength();
	void getReverseComplement(string & _sReverseSequence);
	void getSequence(string & _sReadSequence);
	void setReadUniqueId(UINT64 _iReadId);
	void setSequence(string & _sReadSequence);
	void sortKeyPosition();
	void writeAlignment(stringstream & _sstream);
	void writeAlignmentLongFormat(stringstream & _sstream);

	void updateQueryCoverage(QueryDataset * _query);
	void updateQueryCoverageConsiderMerging(QueryDataset * _query);

	void chooseQueryWithBestAlignment();

private:
	UINT64 iLargestQueryId;
	bool bBeginReverse;
};

class PairKeyPosition{
public:
	INT64 iQueryReadInnerId;
	INT64 iPositionOnSubject;
	INT64 iPositionOnQuery;

	PairKeyPosition();
	~PairKeyPosition();

};

#endif /* SUBJECTREAD_H_ */
