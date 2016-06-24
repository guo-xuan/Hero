/*
 * QueryRead.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "QueryRead.h"
#include "Config.h"

char aDNA[] = { 'A', 'C', 'T', 'G' };

QueryRead::QueryRead() {
	sReadSequence = "";
	iReadUniqueId = -1;
	iReadInnerID = 0;
	bIsGoodRead = true;
	this->aCoverageMatrix = NULL;
	this->aExtraCoverageMatrix = NULL;
	iRawReadLength = 0;
	iMin = 0;
	iMax = 0;
	iMed = 0;
	dActualMismatchRate = 0;
	bMerged = false;
	dGcContent = 0;
}

QueryRead::~QueryRead() {
	deleteCoverageMatrix();
}

void QueryRead::setReadUniqueId(INT64 name) {
	iReadUniqueId = name;
}
void QueryRead::setSequence(string & sequence) {
	sReadSequence = sequence;
	iRawReadLength = sequence.length();
	if (!Config::bCallConsensus) {
		this->sReadSequenceReverseComplement = QueryRead::reverseComplement(this->sReadSequence);
	}
	this->setGcContent();
}

void QueryRead::setReadInnerId(INT64 id) {
	iReadInnerID = id;
}

string * QueryRead::getSequence() {
	return &sReadSequence;
}

string * QueryRead::getReverseComplement() {
	return &sReadSequenceReverseComplement;
}

INT64 QueryRead::getReadInnerId() {
	return iReadInnerID;
}

INT64 QueryRead::getReadUniqueId() {
	return this->iReadUniqueId;
}

INT64 QueryRead::getRawReadLength() {
	return iRawReadLength;
}

INT64 QueryRead::getCurrentReadLength() {
	return this->sReadSequence.length();
}

void QueryRead::setGcContent() {
	this->dGcContent = 0;
	for (size_t i = 0; i < this->sReadSequence.length(); i++) {
		if (this->sReadSequence.at(i) == 'C' || this->sReadSequence.at(i) == 'G') {
			dGcContent++;
		}
	}
	dGcContent = dGcContent / (double) sReadSequence.length();
}

string QueryRead::reverseComplement(const string & read) {
	UINT64 stringLength = read.length();
	string reverse(stringLength, '0');
	for (UINT64 i = 0; i < stringLength; i++)	// Then complement the string. Change A to T, C to G, G to C and T to A.
			{
		if (read[i] & 0X02) // C <==> G
			reverse.at(stringLength - i - 1) = read[i] ^ 0X04;
		else
			// A <==> T
			reverse.at(stringLength - i - 1) = read[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
}
string QueryRead::reverseComplement() {
	return this->QueryRead::reverseComplement(this->sReadSequence);
}

bool QueryRead::initilizeCoverageMatrix() {
	UINT32 len = this->getRawReadLength();
	this->aCoverageMatrix = new UINT32[len * 4];
	len *= 4;
	for (UINT32 i = 0; i < len; i++) {
		aCoverageMatrix[i] = 0;
	}
	if (Config::bMerging) {
		len = QueryRead::iLength4ExtraMatrix * 4;
		this->aExtraCoverageMatrix = new UINT32[len];
		for (UINT32 i = 0; i < len; i++) {
			aExtraCoverageMatrix[i] = 0;
		}
	}
	return true;
}

void QueryRead::deleteCoverageMatrix() {
	delete[] aCoverageMatrix;
	if (Config::bMerging) {
		delete[] aExtraCoverageMatrix;
	}
}

void QueryRead::updateCoveragematrix(UINT64 _pos, char _c) {
	int rowId = (_c >> 1) & 0X03;
	UINT32 * p = NULL;
	p = &aCoverageMatrix[_pos * 4 + rowId];
#pragma omp atomic
	(*p)++;
	//this->aCoverage[_pos][rowId]++;
}

void QueryRead::updateCoveragematrix(UINT64 _pos, int _len) {
	char _c;
	int rowId = 0;
	UINT32 * p = NULL;
	for (INT64 i = 0; i < _len; i++) {
		_c = this->sReadSequence.at(_pos + i);
		rowId = (_c >> 1) & 0X03;
		p = &aCoverageMatrix[(_pos + i) * 4 + rowId];
#pragma omp atomic
		(*p)++;
		//this->aCoverage[_pos+i][rowId]++;
	}

}

void QueryRead::callConsensus() {
	this->updateCoveragematrix(0, (int) this->sReadSequence.length());
	//this->showMatrix();
	int iLargestId = -1, j = 0;
	double fLargest = 0, fLarger = 0, sum = 0, fTemp = 0;
	this->dActualMismatchRate = 0;
	for (size_t i = 0; i < this->sReadSequence.length(); i++) {
		iLargestId = -1;
		fLargest = 0;
		fLarger = 0;
		sum = 0;
		for (j = 0; j < 4; j++) {
			sum += this->aCoverageMatrix[i * 4 + j];
		}
		for (j = 0; j < 4; j++) {
			fTemp = ((double) this->aCoverageMatrix[i * 4 + j]) / sum;
			if (fLargest < fTemp) {
				fLarger = fLargest;
				fLargest = fTemp;
				iLargestId = j;
			}
		}
		if (fLargest - fLarger >= Config::dMinFrequencyCallConsensus) {
			//if (fLargest >= 0.6) {
			//if (fLargest >= Config::dMinFrequencyCallConsensus) {
			if (sReadSequence.at(i) != aDNA[iLargestId]) {
				dActualMismatchRate++;
			}
			this->sReadSequence.at(i) = aDNA[iLargestId];
		}
	}
	dActualMismatchRate = dActualMismatchRate / ((double) sReadSequence.length());
}

void QueryRead::showMatrix() {
	for (int j = 0; j < 4; j++) {
		for (size_t i = 0; i < this->sReadSequence.length(); i++) {
			cout << this->aCoverageMatrix[i * 4 + j] << "\t";
		}
		cout << endl;
	}

}

void QueryRead::calculateCoverageStatistic(vector<UINT32> * _vector) {
	size_t len = this->getRawReadLength();
	_vector->clear();
	UINT32 sum = 0;
	int j = 0;
	for (size_t i = 0; i < len; i++) {
		sum = 0;
		for (j = 0; j < 4; j++) {
			sum += this->aCoverageMatrix[i * 4 + j];
		}
		_vector->push_back(sum);
	}
	sort((*_vector).begin(), (*_vector).end());
	iMin = _vector->at(0);
	iMax = _vector->back();
	if (len % 2 == 0) {
		iMed = (_vector->at(len / 2) + _vector->at(len / 2 - 1)) / 2;
	} else {
		iMed = _vector->at((len - 1) / 2);
	}
}

bool QueryRead::isGoodRead() {
	if (!bIsGoodRead) {
		return false;
	}
	if (iMin < Config::iCoverageDepthFilter) {
		return false;
	}
	return true;
}

int QueryRead::iLength4ExtraMatrix = 20;

double searchsorted(vector<double> & _data, double _value) {
	int iUpper = _data.size();
	int iLower = 0;
	int iMiddle = (iUpper + iLower) / 2;
	while (iLower < iUpper) {
		if (_data.at(iMiddle) > _value) {
			iUpper = iMiddle - 1;
		} else if (_data.at(iMiddle) == _value) {
			iUpper = iMiddle;
		} else {
			iLower = iMiddle + 1;
		}
		iMiddle = (iUpper + iLower) / 2;
	}
	return iMiddle;
}

double KsTest(vector<double> & _data1, vector<double> & _data2) {
	std::sort(_data1.begin(), _data1.end());
	std::sort(_data2.begin(), _data2.end());
	int iTotal = _data1.size() + _data2.size();
	vector<double> cdf1(iTotal);
	vector<double> cdf2(iTotal);
	double dMax = 0, temp = 0;
	for (size_t i = 0; i < _data1.size(); i++) {
		cdf1.at(i) = i;
		cdf2.at(i) = searchsorted(_data2, _data1.at(i));
	}
	for (size_t i = 0; i < _data2.size(); i++) {
		cdf2.at(i + _data1.size()) = i;
		cdf1.at(i + _data1.size()) = searchsorted(_data1, _data2.at(i));
	}
	for (int i = 0; i < iTotal; i++) {
		temp = fabs(cdf1.at(i) - cdf2.at(i));
		if (temp > dMax) {
			dMax = temp;
		}
	}
	return dMax;
}

void QueryRead::mergingPairedReads(QueryRead * _rightRead, string & _sequenceLeft, string & _sequenceRight) {
	this->getExtendSequence(_sequenceLeft);
	_rightRead->getExtendComplementSequence(_sequenceRight);
	//find the start position for merging
	int lenLeftSequence = this->getRawReadLength();
	int lenRightSequence = _rightRead->getRawReadLength();
	int iPositionLeft = 0, iPositionRight = 0, iOld;
	bool bMatch = true;
	for (; iPositionLeft < lenLeftSequence; iPositionLeft++) {
		iOld = iPositionLeft;
		iPositionRight = lenRightSequence - 1;
		bMatch = true;
		for (; iPositionRight >= 0 && iPositionLeft < lenLeftSequence; iPositionRight--, iPositionLeft++) {
			if (this->getBase(iPositionLeft) != _rightRead->getComplementBase(iPositionRight)) {
				//cout << this->getBase(iPositionLeft) << endl;
				//cout << _rightRead->getComplementBase(iPositionRight) << endl;
				bMatch = false;
				break;
			}
		}
		iPositionLeft = iOld;
		if (bMatch) {
			break;
		}
	}
	int iNumMismatch = 0;
	if (bMatch && iPositionLeft < lenLeftSequence) {
		if ((_sequenceLeft.length() + _sequenceRight.length() + (lenLeftSequence - iPositionLeft)) < 40) {
			return;
		}
		//now compare the extend sequences
		for (int i = (lenRightSequence - lenLeftSequence + iPositionLeft - 1), j = 0;
				i >= 0 && j < (int) _sequenceLeft.length(); j++, i--) {
			if (_rightRead->getComplementBase(i) != _sequenceLeft.at(j)) {
				iNumMismatch++;
			}
		}
		for (int i = iPositionLeft - 1, j = 0; i >= 0 && j < (int) _sequenceRight.length(); j++, i--) {
			if (this->getBase(i) != _sequenceRight.at(j)) {
				iNumMismatch++;
			}
		}
		if (iNumMismatch <= 0) {			//merge two reads
			for (int i = (lenRightSequence - lenLeftSequence + iPositionLeft - 1); i >= 0; i--) {
				this->sReadSequence.push_back(_rightRead->getComplementBase(i));
			}
			this->bMerged = true;
			_rightRead->bMerged = true;
			_rightRead->sReadSequence = "";
		}
	}

}

void QueryRead::mergingPairedReadsCoverage(QueryRead * _rightRead, vector<double> & _data1, vector<double> & _data2) {
	//find the start position for merging
	int lenLeftSequence = this->getRawReadLength();
	int lenRightSequence = _rightRead->getRawReadLength();
	int iPositionLeft = 0, iPositionRight = 0, iOld;
	bool bMatch = true;
	for (; iPositionLeft < lenLeftSequence; iPositionLeft++) {
		iOld = iPositionLeft;
		iPositionRight = lenRightSequence;
		bMatch = true;
		for (; iPositionRight >= 0 && iPositionLeft < lenLeftSequence; iPositionRight--, iPositionLeft++) {
			if (this->getBase(iPositionLeft) != this->getComplementBase(iPositionRight)) {
				bMatch = false;
				break;
			}
		}
		iPositionLeft = iOld;
		if (bMatch) {
			break;
		}
	}
	//check if these paired end reads can be merged
	if (lenLeftSequence - iPositionLeft > Config::iMinMergeOverlap) {
		_data1.clear();
		_data2.clear();

	}
	//merge the reads or keep them unchanged

}

void QueryRead::updateExtraCoverageMatrix(UINT64 _posOnSequence, string & _sequence) { //_posOnSequence included
	int count = QueryRead::iLength4ExtraMatrix - 1;
	int index = 0;
	char _c;
	int rowId = 0;
	UINT32 * p = NULL;
	UINT64 len = _sequence.length();
	while (count >= 0 && _posOnSequence < len) {
		_c = _sequence.at(_posOnSequence);
		rowId = (_c >> 1) & 0X03;
		p = &aExtraCoverageMatrix[index * 4 + rowId];
#pragma omp atomic
		(*p)++;
		count--;
		index++;
		_posOnSequence++;
	}
}

char QueryRead::getBase(int _pos) {
	return this->sReadSequence.at(_pos);
}

char QueryRead::getComplementBase(int _pos) {
	//_pos = this->sReadSequence.length() - _pos - 1;
	if (sReadSequence[_pos] & 0X02) // C <==> G
		return sReadSequence[_pos] ^ 0X04;
	else
		// A <==> T
		return sReadSequence[_pos] ^ 0X15;
}

void QueryRead::getExtendSequence(string & _sequence) {
	_sequence.clear();
	int iLargestId = -1, j = 0;
	double fLargest = 0, fTemp = 0;
	for (int i = 0; i < iLength4ExtraMatrix; i++) {
		iLargestId = -1;
		fLargest = 0;
		for (j = 0; j < 4; j++) {
			fTemp = ((double) this->aExtraCoverageMatrix[i * 4 + j]);
			if (fLargest < fTemp) {
				fLargest = fTemp;
				iLargestId = j;
			}
		}
		if (fLargest < 1) {
			break;
		}
		_sequence.push_back(aDNA[iLargestId]);
	}
}

void QueryRead::getExtendComplementSequence(string & _sequence) {
	_sequence.clear();
	int iLargestId = -1, j = 0;
	double fLargest = 0, fTemp = 0;
	for (int i = 0; i < iLength4ExtraMatrix; i++) {
		iLargestId = -1;
		fLargest = 0;
		for (j = 0; j < 4; j++) {
			fTemp = ((double) this->aExtraCoverageMatrix[i * 4 + j]);
			if (fLargest < fTemp) {
				fLargest = fTemp;
				iLargestId = j;
			}
		}
		if (fLargest < 1) {
			break;
		}
		_sequence.push_back(aDNA[iLargestId]);
	}
	//reverse(_sequence.begin(), _sequence.end());
	for (UINT64 i = 0; i < _sequence.length(); i++) { // Then complement the string. Change A to T, C to G, G to C and T to A.
		if (_sequence[i] & 0X02) // C <==> G
			_sequence[i] = _sequence[i] ^ 0X04;
		else
			// A <==> T
			_sequence[i] = _sequence[i] ^ 0X15;
	}
}

UINT32 QueryRead::getCoverageAt(int _pos) {
	UINT32 sum = 0;
	int j = 0;
	for (j = 0; j < 4; j++) {
		sum += this->aCoverageMatrix[_pos * 4 + j];
	}
	return sum;
}
