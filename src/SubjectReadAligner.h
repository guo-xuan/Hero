/*
 * SubjectReadAligner.h
 *
 *  Created on: Oct 12, 2015
 *      Author: xgo
 */

#ifndef SUBJECTREADALIGNER_H_
#define SUBJECTREADALIGNER_H_

#include "Config.h"
#include "HashTableMethod.h"
#include "SubjectDataset.h"
#include "BandAligner.h"
#include "HashTableLongKey.h"

class SubjectReadAligner {
public:

	double dMaxMismatchRate;
	double dMaxIndelRate;
	UINT64 iMinimumOverlapLength;
	HashTableMethod * hashTable;

	double dTimeHashing;
	double dTimeAligning;
	double dTimeLoading;
	double dTimeWriting;
	double dTimeAlignmentConvertToString;

	UINT64 iNumberReadsProcessed;
	uint iIndexShard;

	SubjectReadAligner(HashTableMethod * _hashTable);
	virtual ~SubjectReadAligner();
	void getAlignment(SubjectRead & _subjectRead, string & _subjectSequence, uint & _iOrientation,
			vector<PairKeyPosition *> * _keyArray, BandAligner * _aligner);
	void getAlignment2(SubjectRead & _subjectRead, string & _subjectSequence, uint & _iOrientation,
			vector<PairKeyPosition *> * _keyArray, BandAligner * _aligner);
	bool isFromTheSameAlignment(PairKeyPosition * _keyToBeAdded, PairKeyPosition * _keyBegin,
			PairKeyPosition * & _keyEnd);
	bool start();
	bool searchHashTable(SubjectRead & _read, uint _iOrientation);
	bool processOneSubjectRead(SubjectRead & _subjectRead, uint _iOrientation,
			vector<vector<PairKeyPosition *> *> * _vvKeyArray, BandAligner * _aligner);
	bool searchHashTable2(SubjectRead & _read, uint _iOrientation);
	bool processOneSubjectRead2(SubjectRead & _subjectRead, uint _iOrientation,
			vector<vector<PairKeyPosition *> *> * _vvKeyArray, BandAligner * _aligner);
	bool writeToFile(SubjectDataset * _pSubjectData, string _fileName);
	bool writeToFile2(SubjectDataset * _pSubjectData, string _fileName);
};

#endif /* SUBJECTREADALIGNER_H_ */
