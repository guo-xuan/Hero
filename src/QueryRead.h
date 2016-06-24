/*
 * QueryRead.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYREAD_H_
#define QUERYREAD_H_

#include "Config.h"

class QueryRead {

public:
	string sReadSequence;
	string sReadSequenceReverseComplement;
	INT64 iReadUniqueId;
	UINT32 * aCoverageMatrix;
	UINT32 iRawReadLength;
	UINT32 iMin;
	UINT32 iMed;
	UINT32 iMax;
	double dGcContent;
	double dActualMismatchRate;
	bool bMerged;

	INT64 iReadInnerID; // Unique Identification of the read. start from one. zero means a new read.This is used in hashtable.
	string sReadName;
	bool bIsGoodRead;
	QueryRead();
	~QueryRead();
	bool correctErrors();
	void setReadUniqueId(INT64 name);
	void setSequence(string & sequence);
	void setReadInnerId(INT64 id);

	string * getSequence();
	string * getReverseComplement();
	INT64 getReadInnerId();
	INT64 getReadUniqueId();
	INT64 getRawReadLength();
	INT64 getCurrentReadLength();
	void setGcContent();

	static string reverseComplement(const string & read);
	string reverseComplement();

	//call consensus
	bool initilizeCoverageMatrix();
	void deleteCoverageMatrix();
	void updateCoveragematrix(UINT64 _pos, char _c);
	void updateCoveragematrix(UINT64 _pos, int _len);
	void callConsensus();
	void showMatrix();
	void calculateCoverageStatistic(vector<UINT32> * _vector);
	bool isGoodRead();

	//merging
	static int iLength4ExtraMatrix;
	UINT32 * aExtraCoverageMatrix;
	void mergingPairedReads(QueryRead * _rightRead, string & _sequence1, string & _sequence2);
	void mergingPairedReadsCoverage(QueryRead * _rightRead, vector<double> & _data1, vector<double> & _data2);
	void updateExtraCoverageMatrix(UINT64 _posOnSequence, string & _sequence);
	char getBase(int _pos);
	char getComplementBase(int _pos);
	void getExtendSequence(string & _sequence);
	void getExtendComplementSequence(string & _sequence);
	UINT32 getCoverageAt(int _pos);
};

#endif /* QUERYREAD_H_ */
