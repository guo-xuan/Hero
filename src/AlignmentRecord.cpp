/*
 * AlignmentRecord.cpp
 *
 *  Created on: Oct 16, 2015
 *      Author: xgo
 */

#include "AlignmentRecord.h"

uint AlignmentRecord::FORWARD = 0;
uint AlignmentRecord::REVERSECOMPLEMENT = 1;

AlignmentRecord::AlignmentRecord() {

	iQueryReadUniqueId = 0;
	iPositionOfQueryOnSubject = -1000000000;
	this->iStartPositionOnRef = -1000000000;
	this->iEndPositionOnRef = -1000000000;
	iLengthOnRef = -10000000;
	dMismatchRate = 100;
	dIndelRate = 100;
	dNumMismatch = 0;
	dNumIndel = 0;
	dAlignmentLength = 0;
	iOrientation = 10000;
	this->pQueryReadName = NULL;
	querySequence = NULL;
}

AlignmentRecord::~AlignmentRecord() {
	this->vCigarCount.clear();
	this->vCigarType.clear();
}

void AlignmentRecord::addCigar(INT64 count, char _type) {
	int _count = (int) count;
	if (this->vCigarCount.size() > 0) {
		if (this->vCigarType.back() == _type) {
			this->vCigarCount.back() += _count;
		} else {
			this->vCigarCount.push_back(_count);
			this->vCigarType.push_back(_type);
		}
	} else {
		this->vCigarCount.push_back(_count);
		this->vCigarType.push_back(_type);
	}
}

void AlignmentRecord::addAlignment(const vector<char> * _cigarType, const vector<int> * _cigarCount) {
	for (size_t i = 0; i < _cigarType->size(); i++) {
		this->addCigar(_cigarCount->at(i), _cigarType->at(i));
	}
}

void AlignmentRecord::calculateScore() {
	dNumMismatch = 0;
	dNumIndel = 0;
	dAlignmentLength = 0;
	if (this->vCigarType.size() > 0) {
		if (this->vCigarType.at(0) != 'I' && this->vCigarType.at(0) != 'D') {
			if (this->vCigarType.at(0) == 'M') {
				dAlignmentLength += this->vCigarCount.at(0);
			} else {
				dNumMismatch += this->vCigarCount.at(0);
				dAlignmentLength += this->vCigarCount.at(0);
			}
		}
		for (size_t i = 1; i < this->vCigarType.size() - 1; i++) {
			if (this->vCigarType.at(i) == 'M') {
				dAlignmentLength += this->vCigarCount.at(i);
			} else if (this->vCigarType.at(i) == 'D' || this->vCigarType.at(i) == 'I') {
				dNumIndel += this->vCigarCount.at(i);
				dAlignmentLength += this->vCigarCount.at(i);
			} else {
				dNumMismatch += this->vCigarCount.at(i);
				dAlignmentLength += this->vCigarCount.at(i);
			}
		}
	}

	if (this->vCigarType.size() > 1) {
		if ((this->vCigarType.back()) != 'I' && (this->vCigarType.back()) != 'D') {
			if ((this->vCigarType.back()) == 'M') {
				dAlignmentLength += (this->vCigarCount.back());
			} else {
				dNumMismatch += (this->vCigarCount.back());
				dAlignmentLength += (this->vCigarCount.back());
			}
		}
	}
	if (dAlignmentLength == 0) {
		this->dIndelRate = 100;
		this->dMismatchRate = 100;
	} else {
		this->dIndelRate = dNumIndel / dAlignmentLength;
		this->dMismatchRate = dNumMismatch / dAlignmentLength;
	}

}

bool AlignmentRecord::getCigar(stringstream & _sstream) {
	_sstream << this->iQueryReadUniqueId << ",";
	_sstream << this->iOrientation << ",";
	_sstream << this->iPositionOfQueryOnSubject << ",";
	for (size_t i = 0; i < this->vCigarType.size(); i++) {
		_sstream << this->vCigarCount.at(i);
		_sstream << this->vCigarType.at(i);
	}
	_sstream << ",";
	if (Config::storeOverhang) {
		size_t index = 0;
		for (size_t i = 0; i < this->vCigarType.size(); i++) {
			if (this->vCigarType.at(i) == 'I') {
				for (int j = 0; j < this->vCigarCount.at(i); j++) {
					_sstream << this->querySequence->at(index);
					index++;
				}
				continue;
			}
			if (this->vCigarType.at(i) == 'D') {
				continue;
			} else {
				index += this->vCigarCount.at(i);
			}
		}
		_sstream << ",";
	} else {
		size_t index = 0;
		for (size_t i = 0; i < this->vCigarType.size(); i++) {
			if (this->vCigarType.at(i) == 'I') {
				for (int j = 0; j < this->vCigarCount.at(i); j++) {
					if (i != 0 && i != this->vCigarCount.size() - 1) {
						_sstream << this->querySequence->at(index);
					}
					index++;
				}
				continue;
			}
			if (this->vCigarType.at(i) == 'D') {
				continue;
			} else {
				index += this->vCigarCount.at(i);
			}
		}
		_sstream << ",";
	}
	return true;
}

bool AlignmentRecord::getCigarOnRef(stringstream & _sstream) {
	if (this->vCigarType.at(0) != 'I') {
		_sstream << this->vCigarCount.at(0);
		_sstream << this->vCigarType.at(0);
	}
	int end = ((int) this->vCigarType.size()) - 1;
	for (int i = 1; i < end; i++) {
		_sstream << this->vCigarCount.at(i);
		if (this->vCigarType.at(i) == 'D') {
			_sstream << 'I';
		} else if (this->vCigarType.at(i) == 'I') {
			_sstream << 'D';
		} else {
			_sstream << this->vCigarType.at(i);
		}
	}
	if (this->vCigarType.size() > 1) {
		if (this->vCigarType.back() != 'I') {
			_sstream << this->vCigarCount.back();
			if (this->vCigarType.back() == 'D') {
				_sstream << 'I';
			} else {
				_sstream << this->vCigarType.back();
			}

		}
	}
	_sstream << ",";
	size_t index = 0;
	for (size_t i = 0; i < this->vCigarType.size(); i++) {
		if (this->vCigarType.at(i) == 'I') {
			for (int j = 0; j < this->vCigarCount.at(i); j++) {
				if (i != 0 && i != this->vCigarCount.size() - 1) {
					_sstream << this->querySequence->at(index);
				}
				index++;
			}
			continue;
		}
		if (this->vCigarType.at(i) == 'D') {
			continue;
		} else {
			index += this->vCigarCount.at(i);
		}
	}
	return true;
}

INT64 AlignmentRecord::getCountDeletion() {
	INT64 count = 0;
	for (size_t i = 1; i < this->vCigarType.size(); i++) {
		if (this->vCigarType.at(i) == 'I') {
			count += this->vCigarCount.at(i);
		}
	}
	if (this->vCigarType.back() == 'I') {
		count -= this->vCigarCount.back();
	}
	return count;
}

INT64 AlignmentRecord::getCountInsert() {
	INT64 count = 0;
	for (size_t i = 0; i < this->vCigarType.size(); i++) {
		if (this->vCigarType.at(i) == 'D') {
			count += this->vCigarCount.at(i);
		}
	}

	return count;
}

INT64 AlignmentRecord::getCountMismatch() {
	INT64 count = 0;
	for (size_t i = 0; i < this->vCigarType.size(); i++) {
		if (this->vCigarType.at(i) != 'D' && this->vCigarType.at(i) != 'I' && this->vCigarType.at(i) != 'M') {
			count += this->vCigarCount.at(i);
		}
	}
	return count;
}

INT64 AlignmentRecord::getClipLeft() {
	INT64 clip = 0;
	if (this->vCigarType.at(0) == 'I') {
		clip = this->vCigarCount.at(0);
	} else {
		clip = 0;
	}

	return clip;
}

INT64 AlignmentRecord::getClipRight() {
	INT64 clip = 0;
	if (this->vCigarType.back() == 'I') {
		clip = this->vCigarCount.back();
	} else {
		clip = 0;
	}

	return clip;
}

INT64 AlignmentRecord::getEditDistance() {
	INT64 count = 0;
	count = this->getCountDeletion();
	count += this->getCountInsert();
	count += this->getCountMismatch();
	return count;
}

char AlignmentRecord::getRelativePosition(INT64 _subjectLength) {
	this->getStartPositionOnRef(_subjectLength);
	this->getEndPositionOnRef(_subjectLength);
	if (this->iStartPositionOnRef < 0 && this->iEndPositionOnRef < _subjectLength) {
		return 'L';
	} else if (this->iStartPositionOnRef < 0 && this->iEndPositionOnRef >= _subjectLength) {
		return 'O';
	} else if (this->iStartPositionOnRef >= 0 && this->iEndPositionOnRef < _subjectLength) {
		return 'I';
	}
	return 'R';
}

INT64 AlignmentRecord::getEndPositionOnRef(INT64 _subjectLength) {
	this->iEndPositionOnRef = this->iPositionOfQueryOnSubject + this->iLengthOnRef;
	if (this->vCigarType.back() == 'I') {
		this->iEndPositionOnRef += this->vCigarCount.back();
	}
	this->iEndPositionOnRef--;
	return this->iEndPositionOnRef;
}

INT64 AlignmentRecord::getStartPositionOnRef(INT64 _subjectLength) {
	if (this->vCigarType.at(0) == 'I') {
		this->iStartPositionOnRef = this->iPositionOfQueryOnSubject - this->vCigarCount.at(0);
	} else {
		this->iStartPositionOnRef = this->iPositionOfQueryOnSubject;
	}
	return this->iStartPositionOnRef;
}

char AlignmentRecord::getStrand() {
	if (this->iOrientation == FORWARD) {
		return 'F';
	} else {
		return 'R';
	}
}

bool AlignmentRecord::isGoodAlignment() {
	this->calculateScore();
	if (this->dMismatchRate <= Config::dMaxMismatchRate && this->dIndelRate <= Config::dMaxIndelRate) {
		return true;
	}
	return false;
}

bool AlignmentRecord::isGoodAlignment(const INT64 & _numMisMatch, const INT64 & _numIndel) {
	this->calculateScore();
	if (this->dNumMismatch > _numMisMatch || this->dNumIndel > _numIndel) {
		return false;
	}
	return true;
}

void AlignmentRecord::refineAlignment() {
	if (this->vCigarType.at(0) == 'D') {
		this->iPositionOfQueryOnSubject += this->vCigarCount.at(0);
		this->vCigarType.erase(this->vCigarType.begin());
		this->vCigarCount.erase(this->vCigarCount.begin());
	}
	if (this->vCigarType.back() == 'D') {
		this->vCigarType.pop_back();
		this->vCigarCount.pop_back();
	}
	this->iLengthOnRef = 0;
	if (this->vCigarType.at(0) != 'I') {
		iLengthOnRef += this->vCigarCount.at(0);
	}
	for (size_t i = 1; i < this->vCigarType.size(); i++) {
		if (this->vCigarType.at(i) != 'I') {
			iLengthOnRef += this->vCigarCount.at(i);
		}
	}
}

void AlignmentRecord::setPositionOfQueryOnSubject(INT64 _pos) {
	this->iPositionOfQueryOnSubject = _pos;
}

bool AlignmentRecord::isSubjectCoveredByQuery(INT64 _len) {
	if (this->iPositionOfQueryOnSubject != 0) {
		return false;
	}
	INT64 len = 0;
	if (this->vCigarType.at(0) == 'D') {
		return false;
	} else if (this->vCigarType.at(0) != 'I') {
		len += this->vCigarCount.at(0);
	}
	if (this->vCigarType.size() > 1) {
		for (size_t i = 1; i < this->vCigarType.size() - 1; i++) {
			if (this->vCigarType.at(i) != 'I') {
				len += this->vCigarCount.at(i);
			}
		}
		if (this->vCigarType.back() == 'D') {
			return false;
		} else if (this->vCigarType.back() != 'I') {
			len += this->vCigarCount.back();
		}
	}
	if (len < _len) {
		return false;
	}

	return true;
}

double AlignmentRecord::getAccuracy() {
	double dSubsitution = 0;
	for (size_t i = 0; i < this->vCigarType.size(); i++) {
		if (this->vCigarType.at(i) != 'I' && this->vCigarType.at(i) != 'D' && this->vCigarType.at(i) != 'M') {
			dSubsitution += this->vCigarCount.at(i);
		}
	}

	return dSubsitution;
}

