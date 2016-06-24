/*
 * BandAligner.h
 *
 *  Created on: Oct 16, 2015
 *      Author: xgo
 */

#ifndef BANDALIGNER_H_
#define BANDALIGNER_H_

#include "Config.h"
#include "AlignmentRecord.h"

using namespace std;

#define MAXLENGTHSUBJECT 700
#define MAXLENGTHQUERY 700
#define MAXDUL 100000

class BandAligner {
public:
	vector<int> vCigarCount;
	vector<char> vCigarType;
	static char CWILD;
	static UINT16 IAFFINE;
	static UINT16 IGLOBAL;

	//char cSubjectSeq[MAXLENGTHSUBJECT];
	//char cQuerySeq[MAXLENGTHQUERY];
	char * cSubjectSeq;
	char * cQuerySeq;
	int lengthSubjectSeq;
	int lengthQuerySeq;

	int iPositionOnSubjectSequence;
	int iPositionOnQuerySequence;

	INT64 m_band;
	double dx;
	double dy;
	double sx;
	double ex; // affine gap penalty
	double miss;
	double hit;
	double unobserved_match;

	BandAligner();
	~BandAligner();

	void bandedAlign();
	void bandedAlignAffineGap();
	void calculateAlignment();
	void calculateCigarAndPosition();
	void calculateCigarAndPositionReversed();
	void createMallocArray();
	void deleteMallocArray();
	void getAlignment();
	vector<int> * getCigarCount();
	vector<char> * getCigarType();
	void setAlignLeft();
	void setAlignRight();
	void setBand(INT64 _band);
	void setCost();
	void setQuerySequence(const string & _querySequence, const INT64 & _queryStart,
			const INT64 & _queryEnd);
	void setSequences(const string & _subjectSequence, const INT64 & _subjectStart,
			const INT64 & _subjectEnd, const string & _querySequence, const INT64 & _queryStart,
			const INT64 & _queryEnd); //end is exclusive
	void setSubjectSequence(const string & _subjectSequence, const INT64 & _subjectStart,
			const INT64 & _subjectEnd);
	void printAlignment();
	void printMatrix();

private:
	INT64 left_margin;
	INT64 right_margin;
	double **FM1;
	double **IX;
	double **IY;
	double **BM1;
	double vals[3];
	double temp[3];
	INT64 ei;
	int * edit;
	int si[3];
	int sj[3];
	INT64 icurrentyCigarCount;
	char cCurrentCigarType;
	bool isSubjectQueryReversed;
	bool isRightAlign;
	bool isLeftAlign;
};

#endif /* BANDALIGNER_H_ */
