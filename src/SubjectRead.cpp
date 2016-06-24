/*
 * SubjectRead.cpp
 *
 *  Created on: Oct 16, 2015
 *      Author: xgo
 */

#include "SubjectRead.h"

SubjectRead::SubjectRead() {
	this->iReadUniqueId = -1;
	this->bIsGoodRead = true;
	this->iLargestQueryId = 0;
	this->bBeginReverse = false;
}

SubjectRead::~SubjectRead() {
	this->clearKeyPosition();
	for (size_t i = 0; i < this->vAlignment.size(); i++) {
		delete this->vAlignment.at(i);
	}
	this->vAlignment.clear();
}

void SubjectRead::addAlignment(AlignmentRecord * _record) {
	if (Config::bSpeicialVersion) {
		this->vAlignment.push_back(_record);
		return;
	}
	if (this->vAlignment.empty()) {
		this->iLargestQueryId = _record->iQueryReadUniqueId;
		this->vAlignment.push_back(_record);
		return;
	}
	if ((this->vAlignment.back())->iQueryReadUniqueId == _record->iQueryReadUniqueId) {
		if ((this->vAlignment.back())->dIndelRate + (this->vAlignment.back())->dMismatchRate
				< _record->dIndelRate + _record->dMismatchRate) {
			delete _record;
		} else if ((this->vAlignment.back())->dIndelRate + (this->vAlignment.back())->dMismatchRate
				== _record->dIndelRate + _record->dMismatchRate) {
			if ((this->vAlignment.back())->dAlignmentLength >= _record->dAlignmentLength) {
				delete _record;
			} else {
				delete (this->vAlignment.back());
				(vAlignment.back()) = _record;
			}
		} else {
			delete (this->vAlignment.back());
			(vAlignment.back()) = _record;
		}
	} else if (iLargestQueryId >= _record->iQueryReadUniqueId) {
		for (size_t i = 0; i < this->vAlignment.size(); i++) {
			if (this->vAlignment.at(i)->iQueryReadUniqueId > _record->iQueryReadUniqueId) {
				break;
			} else if (this->vAlignment.at(i)->iQueryReadUniqueId == _record->iQueryReadUniqueId) {
				if (this->vAlignment.at(i)->dIndelRate + this->vAlignment.at(i)->dMismatchRate
						< _record->dIndelRate + _record->dMismatchRate) {
					delete _record;
				} else if (this->vAlignment.at(i)->dIndelRate + this->vAlignment.at(i)->dMismatchRate
						== _record->dIndelRate + _record->dMismatchRate) {
					if (this->vAlignment.at(i)->dAlignmentLength >= _record->dAlignmentLength) {
						delete _record;
					} else {
						delete (this->vAlignment.at(i));
						this->vAlignment.at(i) = _record;
					}
				} else {
					delete (this->vAlignment.at(i));
					(this->vAlignment.at(i)) = _record;
				}
				return;
			}
		}
		this->vAlignment.push_back(_record);
	} else {
		this->iLargestQueryId = _record->iQueryReadUniqueId;
		this->vAlignment.push_back(_record);
	}
}

void SubjectRead::clearKeyPosition() {
	for (size_t i = 0; i < this->vKeyPosition.size(); i++) {
		delete this->vKeyPosition.at(i);
	}
	this->vKeyPosition.clear();
}

void SubjectRead::clearAlignment() {
	for (size_t i = 0; i < this->vAlignment.size(); i++) {
		delete this->vAlignment.at(i);
	}
	this->vAlignment.clear();
}

bool SubjectRead::getAlignment(stringstream & _sstream) {
	for (size_t i = 0; i < this->vAlignment.size(); i++) {
		this->vAlignment.at(i)->getCigar(_sstream);
	}
	return true;
}

bool SubjectRead::getAlignmentLongFormat(stringstream & _sstream) {
	for (size_t i = 0; i < this->vAlignment.size(); i++) {
		this->vAlignment.at(i)->refineAlignment();
		_sstream << *(this->vAlignment.at(i)->pQueryReadName) << "\t";
		_sstream << this->vAlignment.at(i)->getRelativePosition(this->getReadLength()) << "\t";
		_sstream << this->vAlignment.at(i)->getEditDistance() << "\t";
		_sstream << this->vAlignment.at(i)->getStrand() << "\t";
		_sstream << this->vAlignment.at(i)->getClipLeft() << "\t";
		_sstream << this->vAlignment.at(i)->getClipRight() << "\t";
		_sstream << this->vAlignment.at(i)->querySequence->length() << "\t";
		_sstream << this->vAlignment.at(i)->getStartPositionOnRef(this->getReadLength()) << "\t";
		_sstream << this->vAlignment.at(i)->getEndPositionOnRef(this->getReadLength()) << "\t";
		_sstream << this->vAlignment.at(i)->getCountMismatch() << "\t";
		_sstream << this->vAlignment.at(i)->getCountInsert() << "\t";
		_sstream << this->vAlignment.at(i)->getCountDeletion() << "\t";
		this->vAlignment.at(i)->getCigarOnRef(_sstream);
		_sstream << "\n";
	}
	return true;
}

INT64 SubjectRead::getReadUniqueId() {
	return iReadUniqueId;
}

INT64 SubjectRead::getReadLength() {
	return this->sReadSequence.length();
}

void SubjectRead::getReverseComplement(string & _sReverseSequence) {
	UINT64 stringLength = sReadSequence.length();
	_sReverseSequence.resize(stringLength);
	for (UINT64 i = 0; i < stringLength; i++) { // Then complement the string. Change A to T, C to G, G to C and T to A.
		if (sReadSequence[i] & 0X02) // C <==> G
			_sReverseSequence.at(stringLength - i - 1) = sReadSequence[i] ^ 0X04;
		else
			// A <==> T
			_sReverseSequence.at(stringLength - i - 1) = sReadSequence[i] ^ 0X15;
	}
}

void SubjectRead::getSequence(string & _sReadSequence) {
	_sReadSequence = this->sReadSequence;
}

void SubjectRead::setReadUniqueId(UINT64 _iReadId) {
	this->iReadUniqueId = _iReadId;
}

void SubjectRead::setSequence(string & _sReadSequence) {
	this->sReadSequence = _sReadSequence;
	if (Config::bCallConsensus)
		getReverseComplement(sReadSequenceReverseComplement);
}

bool comparePairKeyPosition(const PairKeyPosition * _a, const PairKeyPosition * _b) {
	if (_a->iQueryReadInnerId < _b->iQueryReadInnerId) {
		return true;
	} else if (_a->iQueryReadInnerId == _b->iQueryReadInnerId) {
		if (_a->iPositionOnQuery < _b->iPositionOnQuery) {
			return true;
		} else if (_a->iPositionOnQuery == _b->iPositionOnQuery) {
			if (_a->iPositionOnSubject < _b->iPositionOnSubject) {
				return true;
			} else {
				return false;
			}
		}
		return false;
	}
	return false;
}

void SubjectRead::sortKeyPosition() {
	sort(vKeyPosition.begin(), vKeyPosition.end(), comparePairKeyPosition);
}

void SubjectRead::writeAlignment(stringstream & _sstream) {
	this->sAlignment = "";
	_sstream.clear();
	_sstream.str(std::string());
	_sstream << ">";
	_sstream << iReadUniqueId;
	_sstream << "\n";
	getAlignment(_sstream);
	_sstream << "\n";
	sAlignment = _sstream.str();
	this->clearAlignment();
	/*if(this->iReadUniqueId==36082){
	 cout << sAlignment << endl;
	 }*/
}

void SubjectRead::writeAlignmentLongFormat(stringstream & _sstream) {
	this->sAlignment = "";
	_sstream.clear();
	_sstream.str(std::string());
	_sstream << ">";
	_sstream << sReadName << "\t" << sReadSequence.length();
	_sstream << "\n";
	getAlignmentLongFormat(_sstream);
	_sstream << "\n";
	sAlignment = _sstream.str();
	this->clearAlignment();
}

void SubjectRead::updateQueryCoverage(QueryDataset * _query) {
	UINT64 id = 0;
	UINT64 iPosOnSubjectRead = 0;
	UINT64 iPosOnQueryRead = 0;
	QueryRead * pQueryRead;
	AlignmentRecord * pAlignmentRecord;
	for (size_t i = 0; i < this->vAlignment.size(); i++) {
		pAlignmentRecord = this->vAlignment.at(i);
		//get Query Id
		id = pAlignmentRecord->iQueryReadUniqueId;
		pQueryRead = _query->vQueryReads->at(id);
		//get variants
		iPosOnSubjectRead = pAlignmentRecord->iPositionOfQueryOnSubject;
		iPosOnQueryRead = 0;
		//process the beginning
		if (pAlignmentRecord->vCigarType.at(0) == 'I') {
			iPosOnQueryRead += pAlignmentRecord->vCigarCount.at(0);
		} else if (pAlignmentRecord->vCigarType.at(0) == 'D') {
			iPosOnSubjectRead += pAlignmentRecord->vCigarCount.at(0);
		} else if (pAlignmentRecord->vCigarType.at(0) == 'M') {
			pQueryRead->updateCoveragematrix(iPosOnQueryRead, pAlignmentRecord->vCigarCount.at(0));
			iPosOnQueryRead += pAlignmentRecord->vCigarCount.at(0);
			iPosOnSubjectRead += pAlignmentRecord->vCigarCount.at(0);
		} else {
			if (pAlignmentRecord->vCigarCount.at(0) > 1) {
				for (int k = 0; k < pAlignmentRecord->vCigarCount.at(0); k++) {
					if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
						pQueryRead->updateCoveragematrix(iPosOnQueryRead, this->sReadSequence.at(iPosOnSubjectRead));
					} else {
						pQueryRead->updateCoveragematrix(iPosOnQueryRead,
								this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
					}
					iPosOnQueryRead++;
					iPosOnSubjectRead++;
				}
			} else {
				if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
					pQueryRead->updateCoveragematrix(iPosOnQueryRead, this->sReadSequence.at(iPosOnSubjectRead));
				} else {
					pQueryRead->updateCoveragematrix(iPosOnQueryRead,
							this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
				}
				iPosOnQueryRead++;
				iPosOnSubjectRead++;
			}
		}
		//processing the middle
		if (pAlignmentRecord->vCigarType.size() > 2) {
			for (size_t j = 1; j < pAlignmentRecord->vCigarType.size() - 1; j++) {
				if (pAlignmentRecord->vCigarType.at(j) == 'I') {
					iPosOnQueryRead += pAlignmentRecord->vCigarCount.at(j);
				} else if (pAlignmentRecord->vCigarType.at(j) == 'D') {
					iPosOnSubjectRead += pAlignmentRecord->vCigarCount.at(j);
				} else if (pAlignmentRecord->vCigarType.at(j) == 'M') {
					pQueryRead->updateCoveragematrix(iPosOnQueryRead, pAlignmentRecord->vCigarCount.at(j));
					iPosOnQueryRead += pAlignmentRecord->vCigarCount.at(j);
					iPosOnSubjectRead += pAlignmentRecord->vCigarCount.at(j);
				} else {
					if (pAlignmentRecord->vCigarCount.at(j) > 1) {
						for (int k = 0; k < pAlignmentRecord->vCigarCount.at(j); k++) {
							if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
								pQueryRead->updateCoveragematrix(iPosOnQueryRead,
										this->sReadSequence.at(iPosOnSubjectRead));
							} else {
								pQueryRead->updateCoveragematrix(iPosOnQueryRead,
										this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
							}
							iPosOnQueryRead++;
							iPosOnSubjectRead++;
						}
					} else {
						if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
							pQueryRead->updateCoveragematrix(iPosOnQueryRead,
									this->sReadSequence.at(iPosOnSubjectRead));
						} else {
							pQueryRead->updateCoveragematrix(iPosOnQueryRead,
									this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
						}
						iPosOnQueryRead++;
						iPosOnSubjectRead++;
					}
				}
			}
		}
		if (pAlignmentRecord->vCigarType.size() == 1) {
			continue;
		}
		//processing the end
		if (pAlignmentRecord->vCigarType.back() == 'I') {
			iPosOnQueryRead += pAlignmentRecord->vCigarCount.back();
		} else if (pAlignmentRecord->vCigarType.back() == 'D') {
			iPosOnSubjectRead += pAlignmentRecord->vCigarCount.back();
		} else if (pAlignmentRecord->vCigarType.back() == 'M') {
			pQueryRead->updateCoveragematrix(iPosOnQueryRead, pAlignmentRecord->vCigarCount.back());
			iPosOnQueryRead += pAlignmentRecord->vCigarCount.back();
			iPosOnSubjectRead += pAlignmentRecord->vCigarCount.back();
		} else {
			if (pAlignmentRecord->vCigarCount.back() > 1) {
				for (int k = 0; k < pAlignmentRecord->vCigarCount.back(); k++) {
					if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
						pQueryRead->updateCoveragematrix(iPosOnQueryRead, this->sReadSequence.at(iPosOnSubjectRead));
					} else {
						pQueryRead->updateCoveragematrix(iPosOnQueryRead,
								this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
					}
					iPosOnQueryRead++;
					iPosOnSubjectRead++;
				}
			} else {
				if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
					pQueryRead->updateCoveragematrix(iPosOnQueryRead, this->sReadSequence.at(iPosOnSubjectRead));
				} else {
					pQueryRead->updateCoveragematrix(iPosOnQueryRead,
							this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
				}
				iPosOnQueryRead++;
				iPosOnSubjectRead++;
			}
		}
		//pQueryRead->showMatrix();
	}
	this->clearAlignment();
}

void SubjectRead::updateQueryCoverageConsiderMerging(QueryDataset * _query) {
	UINT64 id = 0;
	UINT64 iPosOnSubjectRead = 0;
	UINT64 iPosOnQueryRead = 0;
	QueryRead * pQueryRead;
	AlignmentRecord * pAlignmentRecord;
	for (size_t i = 0; i < this->vAlignment.size(); i++) {
		pAlignmentRecord = this->vAlignment.at(i);
		//get Query Id
		id = pAlignmentRecord->iQueryReadUniqueId;
		pQueryRead = _query->vQueryReads->at(id);
		//get variants
		iPosOnSubjectRead = pAlignmentRecord->iPositionOfQueryOnSubject;
		iPosOnQueryRead = 0;
		//process the beginning
		if (pAlignmentRecord->vCigarType.at(0) == 'I') {
			iPosOnQueryRead += pAlignmentRecord->vCigarCount.at(0);
		} else if (pAlignmentRecord->vCigarType.at(0) == 'D') {
			iPosOnSubjectRead += pAlignmentRecord->vCigarCount.at(0);
		} else if (pAlignmentRecord->vCigarType.at(0) == 'M') {
			pQueryRead->updateCoveragematrix(iPosOnQueryRead, pAlignmentRecord->vCigarCount.at(0));
			iPosOnQueryRead += pAlignmentRecord->vCigarCount.at(0);
			iPosOnSubjectRead += pAlignmentRecord->vCigarCount.at(0);
		} else {
			if (pAlignmentRecord->vCigarCount.at(0) > 1) {
				for (int k = 0; k < pAlignmentRecord->vCigarCount.at(0); k++) {
					if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
						pQueryRead->updateCoveragematrix(iPosOnQueryRead, this->sReadSequence.at(iPosOnSubjectRead));
					} else {
						pQueryRead->updateCoveragematrix(iPosOnQueryRead,
								this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
					}
					iPosOnQueryRead++;
					iPosOnSubjectRead++;
				}
			} else {
				if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
					pQueryRead->updateCoveragematrix(iPosOnQueryRead, this->sReadSequence.at(iPosOnSubjectRead));
				} else {
					pQueryRead->updateCoveragematrix(iPosOnQueryRead,
							this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
				}
				iPosOnQueryRead++;
				iPosOnSubjectRead++;
			}
		}
		//processing the middle
		if (pAlignmentRecord->vCigarType.size() > 2) {
			for (size_t j = 1; j < pAlignmentRecord->vCigarType.size() - 1; j++) {
				if (pAlignmentRecord->vCigarType.at(j) == 'I') {
					iPosOnQueryRead += pAlignmentRecord->vCigarCount.at(j);
				} else if (pAlignmentRecord->vCigarType.at(j) == 'D') {
					iPosOnSubjectRead += pAlignmentRecord->vCigarCount.at(j);
				} else if (pAlignmentRecord->vCigarType.at(j) == 'M') {
					pQueryRead->updateCoveragematrix(iPosOnQueryRead, pAlignmentRecord->vCigarCount.at(j));
					iPosOnQueryRead += pAlignmentRecord->vCigarCount.at(j);
					iPosOnSubjectRead += pAlignmentRecord->vCigarCount.at(j);
				} else {
					if (pAlignmentRecord->vCigarCount.at(j) > 1) {
						for (int k = 0; k < pAlignmentRecord->vCigarCount.at(j); k++) {
							if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
								pQueryRead->updateCoveragematrix(iPosOnQueryRead,
										this->sReadSequence.at(iPosOnSubjectRead));
							} else {
								pQueryRead->updateCoveragematrix(iPosOnQueryRead,
										this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
							}
							iPosOnQueryRead++;
							iPosOnSubjectRead++;
						}
					} else {
						if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
							pQueryRead->updateCoveragematrix(iPosOnQueryRead,
									this->sReadSequence.at(iPosOnSubjectRead));
						} else {
							pQueryRead->updateCoveragematrix(iPosOnQueryRead,
									this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
						}
						iPosOnQueryRead++;
						iPosOnSubjectRead++;
					}
				}
			}
		}
		if (pAlignmentRecord->vCigarType.size() == 1) {
			continue;
		}
		//processing the end
		if (pAlignmentRecord->vCigarType.back() == 'I') {
			iPosOnQueryRead += pAlignmentRecord->vCigarCount.back();
		} else if (pAlignmentRecord->vCigarType.back() == 'D') {
			iPosOnSubjectRead += pAlignmentRecord->vCigarCount.back();
		} else if (pAlignmentRecord->vCigarType.back() == 'M') {
			pQueryRead->updateCoveragematrix(iPosOnQueryRead, pAlignmentRecord->vCigarCount.back());
			iPosOnQueryRead += pAlignmentRecord->vCigarCount.back();
			iPosOnSubjectRead += pAlignmentRecord->vCigarCount.back();
		} else {
			if (pAlignmentRecord->vCigarCount.back() > 1) {
				for (int k = 0; k < pAlignmentRecord->vCigarCount.back(); k++) {
					if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
						pQueryRead->updateCoveragematrix(iPosOnQueryRead, this->sReadSequence.at(iPosOnSubjectRead));
					} else {
						pQueryRead->updateCoveragematrix(iPosOnQueryRead,
								this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
					}
					iPosOnQueryRead++;
					iPosOnSubjectRead++;
				}
			} else {
				if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
					pQueryRead->updateCoveragematrix(iPosOnQueryRead, this->sReadSequence.at(iPosOnSubjectRead));
				} else {
					pQueryRead->updateCoveragematrix(iPosOnQueryRead,
							this->sReadSequenceReverseComplement.at(iPosOnSubjectRead));
				}
				iPosOnQueryRead++;
				iPosOnSubjectRead++;
			}
		}
		if (iPosOnSubjectRead < this->sReadSequence.length()) {
			if (pAlignmentRecord->iOrientation == AlignmentRecord::FORWARD) {
				pQueryRead->updateExtraCoverageMatrix(iPosOnSubjectRead, this->sReadSequence);
			} else {
				pQueryRead->updateExtraCoverageMatrix(iPosOnSubjectRead, this->sReadSequenceReverseComplement);
			}
		}
		//pQueryRead->showMatrix();
	}
	this->clearAlignment();
}

void SubjectRead::chooseQueryWithBestAlignment() {
	if (this->vAlignment.size() <= 0) {
		return;
	}
	vector<int> viQuarlifiedAlignment;
	for (size_t i = 0; i < vAlignment.size(); i++) {
		if (vAlignment.at(i)->isSubjectCoveredByQuery(this->getReadLength())) {
			viQuarlifiedAlignment.push_back(i);
		}
	}
	if ((int)viQuarlifiedAlignment.size() <= 0) {
		this->clearAlignment();
		return;
	} else if (viQuarlifiedAlignment.size() > 1) {
		int iBestAlignmentIndex = -1;
		double dMinimumSubsitution = this->getReadLength();
		double dTemp = 0;
		for (size_t i = 0; i < viQuarlifiedAlignment.size(); i++) {
			dTemp = this->vAlignment.at(viQuarlifiedAlignment.at(i))->getAccuracy();
			if (dTemp < dMinimumSubsitution) {
				dMinimumSubsitution = dTemp;
				iBestAlignmentIndex = viQuarlifiedAlignment.at(i);
			}
		}
		if (iBestAlignmentIndex == -1) {
			cout << "error in SubjectRead (chooseQueryWithBestAlignment)." << endl;
		}

		for (size_t i = 0; i < this->vAlignment.size(); i++) {
			if ((int) i != iBestAlignmentIndex) {
				delete this->vAlignment.at(i);
			}
		}
		AlignmentRecord * pAlignmentRecord = vAlignment.at(iBestAlignmentIndex);
		this->vAlignment.clear();
		this->vAlignment.push_back(pAlignmentRecord);
	}
}

PairKeyPosition::PairKeyPosition() {
	this->iQueryReadInnerId = -1;
	this->iPositionOnQuery = -1;
	this->iPositionOnSubject = -1;
}

PairKeyPosition::~PairKeyPosition() {

}

