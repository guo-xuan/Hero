/*
 * AlignmentRecord.h
 *
 *  Created on: Oct 16, 2015
 *      Author: xgo
 */

#ifndef ALIGNMENTRECORD_H_
#define ALIGNMENTRECORD_H_

#include "Config.h"

class AlignmentRecord {
public:
	AlignmentRecord();
	~AlignmentRecord();

	uint iOrientation; //the orientation of reference (subject read)
	static uint FORWARD;
	static uint REVERSECOMPLEMENT;

	vector<int> vCigarCount;
	vector<char> vCigarType;
	INT64 iStartPositionOnRef;
	INT64 iEndPositionOnRef;
	INT64 iPositionOfQueryOnSubject;
	INT64 iLengthOnRef;

	double dMismatchRate;
	double dIndelRate;
	double dNumMismatch;
	double dNumIndel;
	double dAlignmentLength;
	UINT64 iQueryReadUniqueId;

	string * querySequence;
	string * pQueryReadName;

	void addAlignment(const vector<char> * _cigarType, const vector<int> * _cigarCount);
	void addCigar(INT64 _count, char _type);
	void calculateScore();
	bool getCigar(stringstream & _sstream);
	bool getCigarOnRef(stringstream & _sstream);
	INT64 getCountDeletion();
	INT64 getCountInsert();
	INT64 getCountMismatch();
	INT64 getClipLeft();
	INT64 getClipRight();
	INT64 getEditDistance();
	char getRelativePosition(INT64 _subjectLength);
	INT64 getEndPositionOnRef(INT64 _subjectLength);
	INT64 getStartPositionOnRef(INT64 _subjectLength);
	char getStrand();
	bool isGoodAlignment();
	bool isGoodAlignment(const INT64 & _numMisMatch, const INT64 & _numIndel);
	void refineAlignment();
	void setPositionOfQueryOnSubject(INT64 _pos);

	bool isSubjectCoveredByQuery(INT64 _iSubjectLength);
	double getAccuracy();
};

#endif /* ALIGNMENTRECORD_H_ */
