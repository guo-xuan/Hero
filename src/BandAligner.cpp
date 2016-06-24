/*
 * BandAligner.cpp
 *
 *  Created on: Oct 16, 2015
 *      Author: xgo
 */

#include "BandAligner.h"

char BandAligner::CWILD = 'N';

UINT16 BandAligner::IAFFINE = 0;
UINT16 BandAligner::IGLOBAL = 1;

BandAligner::BandAligner() {
	lengthSubjectSeq = 0;
	lengthQuerySeq = 0;
	this->cQuerySeq = NULL;
	this->cSubjectSeq = NULL;
	this->FM1 = NULL;
	this->IX = NULL;
	this->IY = NULL;
	this->BM1 = NULL;
	this->edit = NULL;
	iPositionOnSubjectSequence = 0;
	iPositionOnQuerySequence = 0;

	m_band = 0;
	dx = 0.0;
	dy = 0.0;
	sx = 0.0;
	ex = 0.0;
	miss = 0.0;
	hit = 0.0;
	unobserved_match = 0.0;

	left_margin = 0;
	right_margin = 0;

	ei = 0;

	icurrentyCigarCount = 0;
	cCurrentCigarType = 0;

	this->setCost();
	this->createMallocArray();
	this->isSubjectQueryReversed = false;
	this->isLeftAlign = false;
	this->isRightAlign = false;
}

BandAligner::~BandAligner() {
	this->deleteMallocArray();
}

void BandAligner::bandedAlign() {
	FM1[0][0] = 0;
	int i = 0, j = 0, id = 0;
	id = this->lengthSubjectSeq - this->lengthQuerySeq;
	if (id > m_band) {
		m_band = id;
	}
	for (i = 0; i <= this->lengthSubjectSeq; i++) {
		left_margin = 0 > i - m_band ? 0 : i - m_band;
		right_margin = this->lengthQuerySeq < i + m_band ? this->lengthQuerySeq : i + m_band;
		for (j = left_margin; j <= right_margin; j++) {
			if (i == 0 && j == 0)
				continue;
			vals[0] = -MAXDUL;
			vals[1] = -MAXDUL;
			vals[2] = -MAXDUL;
			if (i > 0 && j < i + m_band) {
				vals[0] = FM1[i - 1][j];
				if (j > 0 && j < this->lengthQuerySeq) {
					vals[0] = vals[0] + dx; //insertion in the Sb
				} else {
					if (j == 0 && this->isRightAlign) {
						vals[0] = vals[0] + sx;
					} else if (j >= this->lengthQuerySeq && this->isLeftAlign) {
						vals[0] = vals[0] + sx;
					} else {
						vals[0] = vals[0] + dx; //insertion outside the Sb
					}
				}
			}
			if (j > left_margin) {
				vals[1] = dy + FM1[i][j - 1];
			}
			if (i > 0 && j > 0) {
				vals[2] = FM1[i - 1][j - 1]
						+ (cSubjectSeq[i - 1] != cQuerySeq[j - 1] ?
								(cSubjectSeq[i - 1] == BandAligner::CWILD
										|| cQuerySeq[j - 1] == BandAligner::CWILD) ?
										unobserved_match : miss
								:
								hit);
			}
			//max3
			if (vals[2] >= vals[1] && vals[2] >= vals[0]) {
				FM1[i][j] = vals[2];
				BM1[i][j] = 2;
			} else {
				if (vals[1] >= vals[0]) {
					FM1[i][j] = vals[1];
					BM1[i][j] = 1;
				} else {
					FM1[i][j] = vals[0];
					BM1[i][j] = 0;
				}
			}
		}
	}

	ei = this->lengthQuerySeq + this->lengthSubjectSeq;
	i = this->lengthSubjectSeq;
	j = this->lengthQuerySeq;
	id = BM1[i][j];
	while (true) {
		edit[--ei] = (sj[id] - si[id]);
		i = i + si[id];
		j = j + sj[id];
		id = BM1[i][j];
		if (i == 0 && j == 0)
			break;
	}
}

double max3(double _x, double _y, double _z) {
	if (_x >= _y && _x >= _z) {
		return _x;
	} else if (_y >= _x && _y >= _z) {
		return _y;
	} else {
		return _z;
	}

}

void BandAligner::bandedAlignAffineGap() {
	FM1[0][0] = 0;
	IX[0][0] = dx - ex;
	IY[0][0] = dx - ex;
	int i = 0, j = 0, id = 0;
	id = this->lengthSubjectSeq - this->lengthQuerySeq;
	//m_band *= 2;
	if (id > m_band) {
		m_band = id;
	}
	//initialize
	right_margin = lengthSubjectSeq < i + m_band ? lengthSubjectSeq : i + m_band;
	j = 0;
	for (i = 1; i <= right_margin; i++) {
		if (isRightAlign) {
			IX[i][j] = FM1[i - 1][j] + sx;
			IY[i][j] = dy;
			FM1[i][j] = IX[i][j];
		} else {
			IX[i][j] = IX[i - 1][j] + ex;
			IY[i][j] = IX[i][j];
			FM1[i][j] = IX[i][j];
		}
		BM1[i][j] = 0;
	}
	i = 0;
	right_margin = lengthQuerySeq < i + m_band ? lengthQuerySeq : i + m_band;
	for (j = 1; j <= right_margin; j++) {
		IY[i][j] = IY[i][j - 1] + ex;
		IX[i][j] = IY[i][j];
		FM1[i][j] = IY[i][j];
		BM1[i][j] = 1;
	}
	//banded align
	for (i = 1; i <= this->lengthSubjectSeq; i++) {
		left_margin = 0 > i - m_band ? 0 : i - m_band;
		right_margin = this->lengthQuerySeq < i + m_band ? this->lengthQuerySeq : i + m_band;
		for (j = left_margin; j <= right_margin; j++) {
			if (j == 0) {
				continue;
			}
			if (j == left_margin) { //&& j != 0) {
				IY[i][left_margin] = -MAXDUL;
			} else {
				//if (j > 0) {
				temp[0] = FM1[i][j - 1] + dy;
				temp[1] = IY[i][j - 1] + ex;
				IY[i][j] = temp[0] >= temp[1] ? temp[0] : temp[1];
				//}
			}
			if (j == right_margin && i + m_band <= right_margin) {
				IX[i][j] = -MAXDUL;
			} else {
				if (isLeftAlign && j >= lengthQuerySeq) {
					temp[0] = FM1[i - 1][j] + sx;
					temp[1] = IX[i - 1][j] + sx;
					IX[i][j] = temp[1] > temp[0] ? temp[1] : temp[0];
				} else { //if (j > 0) {
					temp[0] = FM1[i - 1][j] + dx;
					temp[1] = IX[i - 1][j] + ex;
					IX[i][j] = temp[1] > temp[0] ? temp[1] : temp[0];
				}
			}
			//if (j > 0) {
			temp[2] = (cSubjectSeq[i - 1] != cQuerySeq[j - 1] ? miss : hit);
			temp[0] = IX[i - 1][j - 1] + temp[2];
			temp[1] = IY[i - 1][j - 1] + temp[2];
			temp[2] += FM1[i - 1][j - 1];
			FM1[i][j] = max3(temp[2], temp[0], temp[1]);
			//}

			//max3
			/*if(j==0){
			 continue;
			 }*/
			if (FM1[i][j] > IX[i][j] && FM1[i][j] >= IY[i][j]) {
				BM1[i][j] = 2;
			} else {
				if (IY[i][j] > IX[i][j]) {
					BM1[i][j] = 1;
				} else {
					BM1[i][j] = 0;
				}
			}

		}
	}

	ei = this->lengthQuerySeq + this->lengthSubjectSeq;
	i = this->lengthSubjectSeq;
	j = this->lengthQuerySeq;
	id = BM1[i][j];
	while (true) {
		edit[--ei] = (sj[id] - si[id]);
		i = i + si[id];
		j = j + sj[id];
		id = BM1[i][j];
		if (i == 0 && j == 0)
			break;
	}
}

void BandAligner::calculateAlignment() {
	if (Config::iAlignerType == IAFFINE) {
		this->bandedAlignAffineGap();
	} else {
		this->bandedAlign();
	}
	/*if(Config::bDebug){
	 printAlignment();
	 }*/
//printAlignment();
	if (!this->isSubjectQueryReversed) {
		this->calculateCigarAndPosition();
	} else {
		this->calculateCigarAndPositionReversed();
	}
}

void BandAligner::calculateCigarAndPosition() {
	this->vCigarCount.clear();
	this->vCigarType.clear();
	int editLen = this->lengthSubjectSeq + this->lengthQuerySeq;
	int subjectIndex = 0;
	int queryIndex = 0;
	int i = ei;
	switch (edit[i]) {
	case -1:
		this->cCurrentCigarType = 'I';
		this->icurrentyCigarCount = 1;
		queryIndex++;
		break;
	case 0:
		if (this->cQuerySeq[queryIndex] == this->cSubjectSeq[subjectIndex]) {
			this->cCurrentCigarType = 'M';
			this->icurrentyCigarCount = 1;
		} else {
			this->icurrentyCigarCount = 1;
			this->cCurrentCigarType = this->cQuerySeq[queryIndex];
		}
		queryIndex++;
		subjectIndex++;
		break;
	case 1:
		this->cCurrentCigarType = 'D';
		this->icurrentyCigarCount = 1;
		subjectIndex++;
		break;
	}
	i++;
//find the start position
	for (; i < editLen; i++) {
		switch (edit[i]) {
		case -1:
			if (this->cCurrentCigarType == 'I') {
				this->icurrentyCigarCount++;
			} else {
				this->vCigarCount.push_back(icurrentyCigarCount);
				this->vCigarType.push_back(cCurrentCigarType);
				this->cCurrentCigarType = 'I';
				this->icurrentyCigarCount = 1;
			}
			queryIndex++;
			break;
		case 0:
			if (this->cQuerySeq[queryIndex] == this->cSubjectSeq[subjectIndex]) {
				if (this->cCurrentCigarType == 'M') {
					this->icurrentyCigarCount++;
				} else {
					this->vCigarCount.push_back(icurrentyCigarCount);
					this->vCigarType.push_back(cCurrentCigarType);
					this->cCurrentCigarType = 'M';
					this->icurrentyCigarCount = 1;
				}
			} else {
				this->vCigarCount.push_back(icurrentyCigarCount);
				this->vCigarType.push_back(cCurrentCigarType);
				this->icurrentyCigarCount = 1;
				this->cCurrentCigarType = this->cQuerySeq[queryIndex];
			}
			queryIndex++;
			subjectIndex++;
			break;
		case 1:
			if (this->cCurrentCigarType == 'D') {
				this->icurrentyCigarCount++;
			} else {
				this->vCigarCount.push_back(icurrentyCigarCount);
				this->vCigarType.push_back(cCurrentCigarType);
				this->cCurrentCigarType = 'D';
				this->icurrentyCigarCount = 1;
			}
			subjectIndex++;
			break;
		}
	}
	this->vCigarCount.push_back(icurrentyCigarCount);
	this->vCigarType.push_back(cCurrentCigarType);
}

void BandAligner::calculateCigarAndPositionReversed() {
	this->vCigarCount.clear();
	this->vCigarType.clear();
	int editLen = this->lengthSubjectSeq + this->lengthQuerySeq;
	int subjectIndex = 0;
	int queryIndex = 0;
	int i = ei;
	switch (edit[i]) {
	case -1:
		this->cCurrentCigarType = 'D';
		this->icurrentyCigarCount = 1;
		queryIndex++;
		break;
	case 0:
		if (this->cQuerySeq[queryIndex] == this->cSubjectSeq[subjectIndex]) {

			this->cCurrentCigarType = 'M';
			this->icurrentyCigarCount = 1;
		} else {
			this->icurrentyCigarCount = 1;
			this->cCurrentCigarType = this->cSubjectSeq[subjectIndex];
		}
		queryIndex++;
		subjectIndex++;
		break;
	case 1:
		this->cCurrentCigarType = 'I';
		this->icurrentyCigarCount = 1;
		subjectIndex++;
		break;
	}
	i++;
//record the cigar
	for (; i < editLen; i++) {
		switch (edit[i]) {
		case -1:
			if (this->cCurrentCigarType == 'D') {
				this->icurrentyCigarCount++;
			} else {
				this->vCigarCount.push_back(icurrentyCigarCount);
				this->vCigarType.push_back(cCurrentCigarType);
				this->cCurrentCigarType = 'D';
				this->icurrentyCigarCount = 1;
			}
			queryIndex++;
			break;
		case 0:
			if (this->cQuerySeq[queryIndex] == this->cSubjectSeq[subjectIndex]) {
				if (this->cCurrentCigarType == 'M') {
					this->icurrentyCigarCount++;
				} else {
					this->vCigarCount.push_back(icurrentyCigarCount);
					this->vCigarType.push_back(cCurrentCigarType);
					this->cCurrentCigarType = 'M';
					this->icurrentyCigarCount = 1;
				}
			} else {
				this->vCigarCount.push_back(icurrentyCigarCount);
				this->vCigarType.push_back(cCurrentCigarType);
				this->icurrentyCigarCount = 1;
				this->cCurrentCigarType = this->cSubjectSeq[subjectIndex];
			}
			queryIndex++;
			subjectIndex++;
			break;
		case 1:
			if (this->cCurrentCigarType == 'I') {
				this->icurrentyCigarCount++;
			} else {
				this->vCigarCount.push_back(icurrentyCigarCount);
				this->vCigarType.push_back(cCurrentCigarType);
				this->cCurrentCigarType = 'I';
				this->icurrentyCigarCount = 1;
			}
			subjectIndex++;
			break;
		}
	}
	this->vCigarCount.push_back(icurrentyCigarCount);
	this->vCigarType.push_back(cCurrentCigarType);
}

void BandAligner::createMallocArray() {
	this->cSubjectSeq = (char *) malloc(sizeof(char) * Config::iLongestRead);
	this->cQuerySeq = (char *) malloc(sizeof(char) * Config::iLongestRead);
	this->FM1 = (double **) malloc(sizeof(double*) * Config::iLongestRead);
	this->BM1 = (double **) malloc(sizeof(double*) * Config::iLongestRead);
	this->edit = (int *) malloc(sizeof(int) * Config::iLongestRead * 2);
	INT64 i = 0;
	for (i = 0; i < Config::iLongestRead; i++) {
		FM1[i] = (double *) malloc(sizeof(double) * Config::iLongestRead);
		BM1[i] = (double *) malloc(sizeof(double) * Config::iLongestRead);
	}
	if (Config::iAlignerType==BandAligner::IAFFINE) {
		this->IX = (double **) malloc(sizeof(double*) * Config::iLongestRead);
		this->IY = (double **) malloc(sizeof(double*) * Config::iLongestRead);
		for (i = 0; i < Config::iLongestRead; i++) {
			IX[i] = (double *) malloc(sizeof(double) * Config::iLongestRead);
			IY[i] = (double *) malloc(sizeof(double) * Config::iLongestRead);
		}
	}
}

void BandAligner::deleteMallocArray() {
	free(this->cQuerySeq);
	free(this->cSubjectSeq);
	INT64 i = 0;
	for (i = 0; i < Config::iLongestRead; i++) {
		free(FM1[i]);
		free(BM1[i]);
	}
	free(FM1);
	free(BM1);
	free(edit);
	if (Config::iAlignerType==BandAligner::IAFFINE) {
		for (i = 0; i < Config::iLongestRead; i++) {
			free(IX[i]);
			free(IY[i]);
		}
		free(IX);
		free(IY);
	}
}

void BandAligner::getAlignment() {
	this->bandedAlignAffineGap();
	this->calculateCigarAndPosition();
	this->printAlignment();
}

vector<int> * BandAligner::getCigarCount() {
	return &this->vCigarCount;
}

vector<char> * BandAligner::getCigarType() {
	return &this->vCigarType;
}

void BandAligner::printAlignment() {
	int editLen = this->lengthSubjectSeq + this->lengthQuerySeq;
	int _index = 0;
	for (int i = ei; i < editLen; i++) {
		switch (edit[i]) {
		case -1:
			cout << '-';
			break;
		case 0:
			cout << this->cSubjectSeq[_index++];
			break;
		case 1:
			cout << this->cSubjectSeq[_index++];
			break;
		}
	}
	cout << endl;
	_index = 0;
	for (int l = ei; l < editLen; l++) {
		switch (edit[l]) {
		case -1:
			cout << this->cQuerySeq[_index++];
			break;
		case 0:
			cout << this->cQuerySeq[_index++];
			break;
		case 1:
			cout << '-';
			break;
		}
	}
	cout << endl;
	for (uint i = 0; i < vCigarCount.size(); i++) {
		cout << vCigarCount.at(i);
		cout << vCigarType.at(i);
	}
	cout << endl;
}

void BandAligner::setAlignLeft() {
	this->isLeftAlign = true;
}

void BandAligner::setAlignRight() {
	this->isRightAlign = true;
}

void BandAligner::setBand(INT64 _band) {
	this->m_band = _band;
}

void BandAligner::setCost() {
	this->dx = -4;
	this->dy = -4;
	this->miss = -4;
	this->hit = 1;
	this->unobserved_match = -1;
	this->sx = 0;
	this->ex = -1;
	si[0] = -1;
	si[1] = 0;
	si[2] = -1;
	sj[0] = 0;
	sj[1] = -1;
	sj[2] = -1;
}

void BandAligner::setQuerySequence(const string & _querySequence, const INT64 & _queryStart,
		const INT64 & _queryEnd) {
	if (_queryStart < 0) {
		cout << "error in set subject query sequence" << endl;
		exit(1);
	} else {
		this->iPositionOnQuerySequence = _queryStart;
	}
	this->lengthQuerySeq = 0;
	for (int i = iPositionOnQuerySequence; i < _queryEnd; i++) {
		this->cQuerySeq[lengthQuerySeq] = _querySequence.at(i);
		this->lengthQuerySeq++;
	}
}

void BandAligner::setSequences(const string & _subjectSequence, const INT64 & _subjectStart,
		const INT64 & _subjectEnd, const string & _querySequence, const INT64 & _queryStart,
		const INT64 & _queryEnd) {
	if (_queryEnd - _queryStart <= _subjectEnd - _subjectStart) {
		this->setQuerySequence(_querySequence, _queryStart, _queryEnd);
		this->setSubjectSequence(_subjectSequence, _subjectStart, _subjectEnd);
		this->isSubjectQueryReversed = false;
	} else {
		this->setQuerySequence(_subjectSequence, _subjectStart, _subjectEnd);
		this->setSubjectSequence(_querySequence, _queryStart, _queryEnd);
		this->isSubjectQueryReversed = true;
	}
	this->isLeftAlign = false;
	this->isRightAlign = false;
}

void BandAligner::setSubjectSequence(const string & _subjectSequence, const INT64 & _subjectStart,
		const INT64 & _subjectEnd) {
	if (_subjectStart < 0) {
		cout << "error in set subject sequence" << endl;
		exit(1);
	} else {
		this->iPositionOnSubjectSequence = _subjectStart;
	}
	this->lengthSubjectSeq = 0;
	for (int i = _subjectStart; i < _subjectEnd; i++) {
		this->cSubjectSeq[this->lengthSubjectSeq] = _subjectSequence.at(i);
		this->lengthSubjectSeq++;
	}
}

void BandAligner::printMatrix() {
	for (int i = 0; i <= this->lengthSubjectSeq; i++) {
		for (int j = 0; j <= this->lengthQuerySeq; j++) {
			cout << FM1[i][j] << "\t";
		}
		cout << endl;
	}

	for (int i = 0; i <= this->lengthSubjectSeq; i++) {
		for (int j = 0; j <= this->lengthQuerySeq; j++) {
			cout << IX[i][j] << "\t";
		}
		cout << endl;
	}
	for (int i = 0; i <= this->lengthSubjectSeq; i++) {
		for (int j = 0; j <= this->lengthQuerySeq; j++) {
			cout << IY[i][j] << "\t";
		}
		cout << endl;
	}
	for (int i = 0; i <= this->lengthSubjectSeq; i++) {
		for (int j = 0; j <= this->lengthQuerySeq; j++) {
			cout << BM1[i][j] << "\t";
		}
		cout << endl;
	}
}

