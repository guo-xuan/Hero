/*
 * SubjectReadAligner.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: xgo
 */

#include "SubjectReadAligner.h"

SubjectReadAligner::SubjectReadAligner(HashTableMethod * _hashTable) {
	this->dMaxMismatchRate = Config::dMaxMismatchRate;
	this->dMaxIndelRate = Config::dMaxIndelRate;
	this->iMinimumOverlapLength = Config::minimumOverlapLength;
	this->hashTable = _hashTable;
	this->dTimeHashing = 0;
	this->dTimeAligning = 0;
	dTimeLoading = 0;
	dTimeWriting = 0;
	dTimeAlignmentConvertToString = 0;
	iNumberReadsProcessed = 0;
	iIndexShard = 1;
}

SubjectReadAligner::~SubjectReadAligner() {
	// TODO Auto-generated destructor stub
}

void SubjectReadAligner::getAlignment(SubjectRead & _subjectRead, string & _subjectSequence, uint & _iOrientation,
		vector<PairKeyPosition *> * _keyArray, BandAligner * _aligner) {
	string * _querySequence =
			this->hashTable->dataSet->getReadFromID(_keyArray->at(0)->iQueryReadInnerId)->getSequence();
	INT64 positionOnQuery = _keyArray->at(0)->iPositionOnQuery;
	INT64 positionOnSubject = _keyArray->at(0)->iPositionOnSubject;
	INT64 positionOnQueryCurrent = 0;
	INT64 positionOnSubjectCurrent = 0;
	INT64 lenOverlapEstimated = 0;
	lenOverlapEstimated += positionOnQuery < positionOnSubject ? positionOnQuery : positionOnSubject;
	lenOverlapEstimated +=
			_querySequence->size() - positionOnQuery < _subjectSequence.size() - positionOnSubject ?
					_querySequence->size() - positionOnQuery : _subjectSequence.size() - positionOnSubject;
	INT64 extraResidues = round(lenOverlapEstimated * Config::dMaxIndelRate);
	INT64 maxMisMatch = round(lenOverlapEstimated * Config::dMaxMismatchRate);
	INT64 maxIndel = extraResidues;
	AlignmentRecord * alignmentRecord = new AlignmentRecord();
	alignmentRecord->iOrientation = _iOrientation;
	alignmentRecord->iQueryReadUniqueId =
			this->hashTable->dataSet->getReadFromID(_keyArray->at(0)->iQueryReadInnerId)->getReadUniqueId();
	alignmentRecord->querySequence = _querySequence;
	_aligner->setBand(extraResidues * 2);
	//cout << "begin" << endl;
	//the start part needs to be aligned
	if ((positionOnQuery != 0 && positionOnSubject != 0)) {
		if (positionOnQuery < positionOnSubject) {
			_aligner->setSequences(_subjectSequence,
					(positionOnSubject - positionOnQuery - extraResidues < 0 ?
							0 : positionOnSubject - positionOnQuery - extraResidues), positionOnSubject,
					(*_querySequence), 0, positionOnQuery);
			if (positionOnSubject - positionOnQuery - extraResidues < 0) {
				alignmentRecord->setPositionOfQueryOnSubject(0);
			} else {
				alignmentRecord->setPositionOfQueryOnSubject(positionOnSubject - positionOnQuery - extraResidues);
			}
		} else if (positionOnQuery == positionOnSubject) {
			_aligner->setSequences(_subjectSequence, 0, positionOnSubject, (*_querySequence), 0, positionOnQuery);
			alignmentRecord->setPositionOfQueryOnSubject(0);
		} else {
			_aligner->setSequences(_subjectSequence, 0, positionOnSubject, (*_querySequence),
					(positionOnQuery - positionOnSubject - extraResidues < 0 ?
							0 : positionOnQuery - positionOnSubject - extraResidues), positionOnQuery);
			alignmentRecord->setPositionOfQueryOnSubject(0);
			if ((positionOnQuery - positionOnSubject - extraResidues) > 0) {
				alignmentRecord->addCigar(positionOnQuery - positionOnSubject - extraResidues, 'I');
			}
		}
		_aligner->setAlignRight();
		_aligner->calculateAlignment();
		alignmentRecord->addAlignment(_aligner->getCigarType(), _aligner->getCigarCount());
		alignmentRecord->addCigar(Config::iKmer, 'M');
	} else {
		if (positionOnSubject == 0 && positionOnQuery == 0) {
			alignmentRecord->setPositionOfQueryOnSubject(0);
		} else if (positionOnSubject == 0) {
			alignmentRecord->setPositionOfQueryOnSubject(0);
			alignmentRecord->addCigar(positionOnQuery, 'I');
		} else {
			alignmentRecord->setPositionOfQueryOnSubject(positionOnSubject);
		}
		alignmentRecord->addCigar(Config::iKmer, 'M');
	}
	if (!(alignmentRecord->isGoodAlignment(maxMisMatch, maxIndel))) {
		delete alignmentRecord;
		return;
	}
	//cout << "middle" << endl;
	//the middle part needs to be aligned
	for (size_t i = 1; i < _keyArray->size(); i++) {
		positionOnQueryCurrent = _keyArray->at(i)->iPositionOnQuery;
		positionOnSubjectCurrent = _keyArray->at(i)->iPositionOnSubject;
		if (positionOnQueryCurrent - positionOnQuery == positionOnSubjectCurrent - positionOnSubject
				&& positionOnSubject + Config::iKmer >= positionOnSubjectCurrent) {
			if (positionOnSubject + Config::iKmer > positionOnSubjectCurrent) {
				alignmentRecord->addCigar(positionOnSubjectCurrent - positionOnSubject, 'M');
			} else {
				alignmentRecord->addCigar(Config::iKmer, 'M');
			}
		} else {
			if (positionOnSubject + Config::iKmer == positionOnSubjectCurrent) {
				if (positionOnQueryCurrent - positionOnQuery - Config::iKmer > 0) {
					alignmentRecord->addCigar(positionOnQueryCurrent - positionOnQuery - Config::iKmer, 'I');
					alignmentRecord->addCigar(Config::iKmer, 'M');
				} else {
					alignmentRecord->addCigar(-positionOnQueryCurrent + positionOnQuery + Config::iKmer, 'D');
					alignmentRecord->addCigar(positionOnQueryCurrent - positionOnQuery, 'M');
				}
			} else if (positionOnSubject + Config::iKmer > positionOnSubjectCurrent) {
				if (positionOnQueryCurrent - positionOnQuery < positionOnSubjectCurrent - positionOnSubject) {
					alignmentRecord->addCigar(
							positionOnSubjectCurrent - positionOnSubject - positionOnQueryCurrent + positionOnQuery,
							'D');
					alignmentRecord->addCigar(positionOnQueryCurrent - positionOnQuery, 'M');
				} else {
					alignmentRecord->addCigar(
							positionOnQueryCurrent - positionOnQuery + positionOnSubject - positionOnSubjectCurrent,
							'I');
					alignmentRecord->addCigar(positionOnSubjectCurrent - positionOnSubject, 'M');
				}
			} else {
				if (positionOnQuery + Config::iKmer == positionOnQueryCurrent) {
					alignmentRecord->addCigar(positionOnSubjectCurrent - positionOnSubject - Config::iKmer,
							'D');
					alignmentRecord->addCigar(Config::iKmer, 'M');
				} else if (positionOnQuery + Config::iKmer > positionOnQueryCurrent) {
					alignmentRecord->addCigar(
							positionOnSubjectCurrent - positionOnSubject + positionOnQuery - positionOnQueryCurrent,
							'D');
					alignmentRecord->addCigar(Config::iKmer + positionOnQuery - positionOnQueryCurrent, 'M');
				} else {
					_aligner->setSequences(_subjectSequence, positionOnSubject + Config::iKmer,
							positionOnSubjectCurrent, (*_querySequence), positionOnQuery + Config::iKmer,
							positionOnQueryCurrent);
					_aligner->calculateAlignment();
					alignmentRecord->addAlignment(_aligner->getCigarType(), _aligner->getCigarCount());
					alignmentRecord->addCigar(Config::iKmer, 'M');
				}
			}
		}
		positionOnQuery = positionOnQueryCurrent;
		positionOnSubject = positionOnSubjectCurrent;
		if (!(alignmentRecord->isGoodAlignment(maxMisMatch, maxIndel))) {
			delete alignmentRecord;
			return;
		}
	}
	//cout << "end" << endl;
	//the end part needs to be aligned
	positionOnQuery = (_keyArray->back())->iPositionOnQuery + Config::iKmer;
	positionOnSubject = (_keyArray->back())->iPositionOnSubject + Config::iKmer;
	INT64 tmpLenQuery = _querySequence->length();
	INT64 tempLenSub = _subjectSequence.length();
	if (positionOnQuery < tmpLenQuery) {
		if (tempLenSub - positionOnSubject < tmpLenQuery - positionOnQuery) {
			_aligner->setSequences(_subjectSequence, positionOnSubject, tempLenSub, (*_querySequence), positionOnQuery,
					(tempLenSub - positionOnSubject + positionOnQuery + extraResidues < tmpLenQuery ?
							tempLenSub - positionOnSubject + positionOnQuery + extraResidues : tmpLenQuery));
			_aligner->setAlignLeft();
			_aligner->calculateAlignment();
			alignmentRecord->addAlignment(_aligner->getCigarType(), _aligner->getCigarCount());
			if (tempLenSub - positionOnSubject + positionOnQuery + extraResidues < tmpLenQuery) {
				alignmentRecord->addCigar(
						tmpLenQuery - (tempLenSub - positionOnSubject + positionOnQuery + extraResidues), 'I');
			}
		} else if (tempLenSub - positionOnSubject == tmpLenQuery - positionOnQuery) {
			_aligner->setSequences(_subjectSequence, positionOnSubject, tempLenSub, (*_querySequence), positionOnQuery,
					tmpLenQuery);
			_aligner->setAlignLeft();
			_aligner->calculateAlignment();
			alignmentRecord->addAlignment(_aligner->getCigarType(), _aligner->getCigarCount());
		} else {
			_aligner->setSequences(_subjectSequence, positionOnSubject,
					(positionOnSubject + (tmpLenQuery - positionOnQuery) + extraResidues < tempLenSub ?
							positionOnSubject + (tmpLenQuery - positionOnQuery) + extraResidues : tempLenSub),
					(*_querySequence), positionOnQuery, tmpLenQuery);
			_aligner->setAlignLeft();
			_aligner->calculateAlignment();
			alignmentRecord->addAlignment(_aligner->getCigarType(), _aligner->getCigarCount());
		}
	}
	//cout << "last" << endl;
	if (!(alignmentRecord->isGoodAlignment())) {
		delete alignmentRecord;
	} else {
		_subjectRead.addAlignment(alignmentRecord);
	}
	//cout << "final" << endl;
}

void SubjectReadAligner::getAlignment2(SubjectRead & _subjectRead, string & _subjectSequence, uint & _iOrientation,
		vector<PairKeyPosition *> * _keyArray, BandAligner * _aligner) {
	string * _querySequence;
	if (Config::bCallConsensus) {
		_querySequence = this->hashTable->dataSet->getReadFromID(_keyArray->at(0)->iQueryReadInnerId)->getSequence();
	} else {
		if (_iOrientation == AlignmentRecord::FORWARD) {
			_querySequence =
					this->hashTable->dataSet->getReadFromID(_keyArray->at(0)->iQueryReadInnerId)->getSequence();
		} else {
			_querySequence =
					this->hashTable->dataSet->getReadFromID(_keyArray->at(0)->iQueryReadInnerId)->getReverseComplement();
		}
	}

	INT64 positionOnQuery = _keyArray->at(0)->iPositionOnQuery;
	INT64 positionOnSubject = _keyArray->at(0)->iPositionOnSubject;
	INT64 positionOnQueryCurrent = 0;
	INT64 positionOnSubjectCurrent = 0;
	INT64 lenOverlapEstimated = 0;
	lenOverlapEstimated += positionOnQuery < positionOnSubject ? positionOnQuery : positionOnSubject;
	lenOverlapEstimated +=
			_querySequence->size() - positionOnQuery < _subjectSequence.size() - positionOnSubject ?
					_querySequence->size() - positionOnQuery : _subjectSequence.size() - positionOnSubject;
	if (lenOverlapEstimated < Config::minimumOverlapLength) {
	//if (lenOverlapEstimated < (_querySequence->size()*0.3)) {
		return;
	}
	INT64 extraResidues = round(lenOverlapEstimated * Config::dMaxIndelRate);
	INT64 maxMisMatch = round(lenOverlapEstimated * Config::dMaxMismatchRate);
	INT64 maxIndel = extraResidues;
	AlignmentRecord * alignmentRecord = new AlignmentRecord();
	alignmentRecord->iOrientation = _iOrientation;
	if (Config::bCallConsensus) {
		alignmentRecord->iQueryReadUniqueId = _keyArray->at(0)->iQueryReadInnerId;
	} else {
		alignmentRecord->iQueryReadUniqueId = this->hashTable->dataSet->getReadFromID(
				_keyArray->at(0)->iQueryReadInnerId)->getReadUniqueId();
	}
	alignmentRecord->querySequence = _querySequence;
	if (Config::bSpeicialVersion) {
		alignmentRecord->pQueryReadName =
				&(this->hashTable->dataSet->getReadFromID(_keyArray->at(0)->iQueryReadInnerId)->sReadName);
	}
	_aligner->setBand(extraResidues * 2);
	//cout << "begin" << endl;
	//the start part needs to be aligned
	if ((positionOnQuery != 0 && positionOnSubject != 0)) {
		if (positionOnQuery < positionOnSubject) {
			_aligner->setSequences(_subjectSequence,
					(positionOnSubject - positionOnQuery - extraResidues < 0 ?
							0 : positionOnSubject - positionOnQuery - extraResidues), positionOnSubject,
					(*_querySequence), 0, positionOnQuery);
			if (positionOnSubject - positionOnQuery - extraResidues < 0) {
				alignmentRecord->setPositionOfQueryOnSubject(0);
			} else {
				alignmentRecord->setPositionOfQueryOnSubject(positionOnSubject - positionOnQuery - extraResidues);
			}
		} else if (positionOnQuery == positionOnSubject) {
			_aligner->setSequences(_subjectSequence, 0, positionOnSubject, (*_querySequence), 0, positionOnQuery);
			alignmentRecord->setPositionOfQueryOnSubject(0);
		} else {
			_aligner->setSequences(_subjectSequence, 0, positionOnSubject, (*_querySequence),
					(positionOnQuery - positionOnSubject - extraResidues < 0 ?
							0 : positionOnQuery - positionOnSubject - extraResidues), positionOnQuery);
			alignmentRecord->setPositionOfQueryOnSubject(0);
			if ((positionOnQuery - positionOnSubject - extraResidues) > 0) {
				alignmentRecord->addCigar(positionOnQuery - positionOnSubject - extraResidues, 'I');
			}
		}
		_aligner->setAlignRight();
		_aligner->calculateAlignment();
		alignmentRecord->addAlignment(_aligner->getCigarType(), _aligner->getCigarCount());
		alignmentRecord->addCigar(Config::iKmer, 'M');
	} else {
		if (positionOnSubject == 0 && positionOnQuery == 0) {
			alignmentRecord->setPositionOfQueryOnSubject(0);
		} else if (positionOnSubject == 0) {
			alignmentRecord->setPositionOfQueryOnSubject(0);
			alignmentRecord->addCigar(positionOnQuery, 'I');
		} else {
			alignmentRecord->setPositionOfQueryOnSubject(positionOnSubject);
		}
		alignmentRecord->addCigar(Config::iKmer, 'M');
	}
	if (!(alignmentRecord->isGoodAlignment(maxMisMatch, maxIndel))) {
		delete alignmentRecord;
		return;
	}
	//cout << "middle" << endl;
	//the middle part needs to be aligned
	for (size_t i = 1; i < _keyArray->size(); i++) {
		positionOnQueryCurrent = _keyArray->at(i)->iPositionOnQuery;
		positionOnSubjectCurrent = _keyArray->at(i)->iPositionOnSubject;
		if (positionOnQueryCurrent - positionOnQuery == positionOnSubjectCurrent - positionOnSubject
				&& positionOnSubject + Config::iKmer >= positionOnSubjectCurrent) {
			if (positionOnSubject + Config::iKmer > positionOnSubjectCurrent) {
				alignmentRecord->addCigar(positionOnSubjectCurrent - positionOnSubject, 'M');
			} else {
				alignmentRecord->addCigar(Config::iKmer, 'M');
			}
		} else {
			if (positionOnSubject + Config::iKmer == positionOnSubjectCurrent) {
				if (positionOnQueryCurrent - positionOnQuery - Config::iKmer > 0) {
					alignmentRecord->addCigar(positionOnQueryCurrent - positionOnQuery - Config::iKmer, 'I');
					alignmentRecord->addCigar(Config::iKmer, 'M');
				} else {
					alignmentRecord->addCigar(-positionOnQueryCurrent + positionOnQuery + Config::iKmer, 'D');
					alignmentRecord->addCigar(positionOnQueryCurrent - positionOnQuery, 'M');
				}
			} else if (positionOnSubject + Config::iKmer > positionOnSubjectCurrent) {
				if (positionOnQueryCurrent - positionOnQuery < positionOnSubjectCurrent - positionOnSubject) {
					alignmentRecord->addCigar(
							positionOnSubjectCurrent - positionOnSubject - positionOnQueryCurrent + positionOnQuery,
							'D');
					alignmentRecord->addCigar(positionOnQueryCurrent - positionOnQuery, 'M');
				} else {
					alignmentRecord->addCigar(
							positionOnQueryCurrent - positionOnQuery + positionOnSubject - positionOnSubjectCurrent,
							'I');
					alignmentRecord->addCigar(positionOnSubjectCurrent - positionOnSubject, 'M');
				}
			} else {
				if (positionOnQuery + Config::iKmer == positionOnQueryCurrent) {
					alignmentRecord->addCigar(positionOnSubjectCurrent - positionOnSubject - Config::iKmer,
							'D');
					alignmentRecord->addCigar(Config::iKmer, 'M');
				} else if (positionOnQuery + Config::iKmer > positionOnQueryCurrent) {
					alignmentRecord->addCigar(
							positionOnSubjectCurrent - positionOnSubject + positionOnQuery - positionOnQueryCurrent,
							'D');
					alignmentRecord->addCigar(Config::iKmer + positionOnQuery - positionOnQueryCurrent, 'M');
				} else {
					_aligner->setSequences(_subjectSequence, positionOnSubject + Config::iKmer,
							positionOnSubjectCurrent, (*_querySequence), positionOnQuery + Config::iKmer,
							positionOnQueryCurrent);
					_aligner->calculateAlignment();
					alignmentRecord->addAlignment(_aligner->getCigarType(), _aligner->getCigarCount());
					alignmentRecord->addCigar(Config::iKmer, 'M');
				}
			}
		}
		positionOnQuery = positionOnQueryCurrent;
		positionOnSubject = positionOnSubjectCurrent;
		if (!(alignmentRecord->isGoodAlignment(maxMisMatch, maxIndel))) {
			delete alignmentRecord;
			return;
		}
	}
	//cout << "end" << endl;
	//the end part needs to be aligned
	positionOnQuery = (_keyArray->back())->iPositionOnQuery + Config::iKmer;
	positionOnSubject = (_keyArray->back())->iPositionOnSubject + Config::iKmer;
	INT64 tmpLenQuery = _querySequence->length();
	INT64 tempLenSub = _subjectSequence.length();
	if (positionOnQuery < tmpLenQuery) {
		if (positionOnSubject == tempLenSub) {
			alignmentRecord->addCigar(tmpLenQuery - positionOnQuery, 'I');
		} else if (tempLenSub - positionOnSubject < tmpLenQuery - positionOnQuery) {
			_aligner->setSequences(_subjectSequence, positionOnSubject, tempLenSub, (*_querySequence), positionOnQuery,
					(tempLenSub - positionOnSubject + positionOnQuery + extraResidues < tmpLenQuery ?
							tempLenSub - positionOnSubject + positionOnQuery + extraResidues : tmpLenQuery));
			_aligner->setAlignLeft();
			_aligner->calculateAlignment();
			alignmentRecord->addAlignment(_aligner->getCigarType(), _aligner->getCigarCount());
			if (tempLenSub - positionOnSubject + positionOnQuery + extraResidues < tmpLenQuery) {
				alignmentRecord->addCigar(
						tmpLenQuery - (tempLenSub - positionOnSubject + positionOnQuery + extraResidues), 'I');
			}
		} else if (tempLenSub - positionOnSubject == tmpLenQuery - positionOnQuery) {
			_aligner->setSequences(_subjectSequence, positionOnSubject, tempLenSub, (*_querySequence), positionOnQuery,
					tmpLenQuery);
			_aligner->setAlignLeft();
			_aligner->calculateAlignment();
			alignmentRecord->addAlignment(_aligner->getCigarType(), _aligner->getCigarCount());
		} else {
			_aligner->setSequences(_subjectSequence, positionOnSubject,
					(positionOnSubject + (tmpLenQuery - positionOnQuery) + extraResidues < tempLenSub ?
							positionOnSubject + (tmpLenQuery - positionOnQuery) + extraResidues : tempLenSub),
					(*_querySequence), positionOnQuery, tmpLenQuery);
			_aligner->setAlignLeft();
			_aligner->calculateAlignment();
			alignmentRecord->addAlignment(_aligner->getCigarType(), _aligner->getCigarCount());
		}
	}
	//cout << "last" << endl;
	if (!(alignmentRecord->isGoodAlignment())) {
		delete alignmentRecord;
	} else {
		_subjectRead.addAlignment(alignmentRecord);
	}
	//cout << "final" << endl;
}

bool SubjectReadAligner::isFromTheSameAlignment(PairKeyPosition * _keyToBeAdded, PairKeyPosition * _keyBegin,
		PairKeyPosition * & _keyEnd) {
	if (_keyToBeAdded->iPositionOnSubject <= _keyEnd->iPositionOnSubject) {
		return false;
	}
	if (_keyToBeAdded->iPositionOnQuery <= _keyEnd->iPositionOnQuery) {
		return false;
	}
	double lengthQuery = fabs(_keyToBeAdded->iPositionOnQuery - _keyBegin->iPositionOnQuery);
	double lengthSubject = fabs(_keyToBeAdded->iPositionOnSubject - _keyBegin->iPositionOnSubject);
	double totalLength = (lengthQuery >= lengthSubject ? lengthQuery : lengthSubject);
	totalLength += Config::iKmer;
	double difference = fabs(lengthQuery - lengthSubject);
	if ((difference / totalLength) > Config::dMaxIndelRate) {
		return false;
	}
	return true;
}

bool SubjectReadAligner::searchHashTable(SubjectRead & _read, uint _iOrientation) {
	string sequence;
	if (_iOrientation == AlignmentRecord::FORWARD) {
		_read.getSequence(sequence);
	} else {
		_read.getReverseComplement(sequence);
	}
	_read.clearKeyPosition();
	INT64 end = ((INT64) sequence.size()) - Config::iKmer;
	for (INT64 i = 0; i < end; i++) {
		vector<INT64> * listOfReads = this->hashTable->getListOfReads(sequence, i, Config::iKmer);
		if (listOfReads != NULL && !listOfReads->empty()) { // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
			for (UINT64 k = 0; k < listOfReads->size(); k++) { // For each such reads.
				INT64 data = listOfReads->at(k); // We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				INT64 queryReadInnerId = data & (hashTable->iMask);
				INT64 queryReadUniqueId = this->hashTable->dataSet->getReadFromID(queryReadInnerId)->getReadUniqueId();
				if (queryReadUniqueId == _read.getReadUniqueId())
					continue;
				INT64 queryPosition = data >> (this->hashTable->iBitOfReadId);
				PairKeyPosition * pkp = new PairKeyPosition();
				pkp->iQueryReadInnerId = queryReadInnerId;
				pkp->iPositionOnSubject = i;
				pkp->iPositionOnQuery = queryPosition;
				_read.vKeyPosition.push_back(pkp);
			}
		}
	}
	if (!_read.vKeyPosition.empty()) {
		_read.sortKeyPosition();
		return true;
	}
	return false;
}

bool SubjectReadAligner::searchHashTable2(SubjectRead & _read, uint _iOrientation) {
	string sequence;
	if (_iOrientation == AlignmentRecord::FORWARD) {
		_read.getSequence(sequence);
	} else {
		_read.getReverseComplement(sequence);
	}
	_read.clearKeyPosition();
	INT64 end = ((INT64) sequence.size()) - Config::iKmer;
	INT64 len = (INT64) sequence.size();
	/*if (_read.iReadUniqueId == 2574) {
	 cout << "check" << endl;
	 }*/
	for (INT64 i = 0; i <= end; i++) {
		vector<INT64> * listOfReads = this->hashTable->getListOfReads(sequence, i, Config::iKmer);
		if (listOfReads != NULL && !listOfReads->empty()) { // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
			for (UINT64 k = 0; k < listOfReads->size(); k++) { // For each such reads.
				INT64 data = listOfReads->at(k); // We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				INT64 queryReadInnerId = data & (hashTable->iMask);
				INT64 queryReadUniqueId = this->hashTable->dataSet->getReadFromID(queryReadInnerId)->getReadUniqueId();
				if ((!Config::bSpeicialVersion) && queryReadUniqueId == _read.getReadUniqueId())
					continue;
				INT64 queryPosition = data >> (this->hashTable->iBitOfReadId);
				PairKeyPosition * pkp = new PairKeyPosition();
				pkp->iQueryReadInnerId = queryReadInnerId;
				pkp->iPositionOnSubject = i;
				pkp->iPositionOnQuery = queryPosition;
				/*if (queryReadUniqueId == 8544) {
				 cout << "check 2" <<endl;
				 }*/
				if (Config::bCallConsensus) { //it is for getting the consensus on the fly, so use the read in Query as the target, always in forward direction
					_read.vKeyPosition.push_back(pkp);
					continue;
				} else if (_iOrientation != AlignmentRecord::FORWARD) {
					pkp->iPositionOnSubject = len - (pkp->iPositionOnSubject + Config::iKmer);
					pkp->iPositionOnQuery = (this->hashTable->dataSet->getReadFromID(queryReadInnerId)->getRawReadLength())
							- (pkp->iPositionOnQuery + Config::iKmer);
				}
				_read.vKeyPosition.push_back(pkp);
			}
		}
	}
	if (!_read.vKeyPosition.empty()) {
		_read.sortKeyPosition();
		return true;
	}
	return false;
}

bool SubjectReadAligner::start() {

	int count = 1;
	SubjectDataset * pSubjectData = new SubjectDataset();
	CLOCKSTART
	MEMORYSTART
	double beginTime = omp_get_wtime();
	while (pSubjectData->loadNextChunkParallel(Config::numberOfThreads)) {
		double endTime = omp_get_wtime();
		dTimeLoading += endTime - beginTime;
		if (pSubjectData->iLongestReadLength < this->hashTable->dataSet->longestReadLength) {
			Config::iLongestRead = pSubjectData->iLongestReadLength;
		} else {
			Config::iLongestRead = this->hashTable->dataSet->longestReadLength;
		}
		//cout << "Loading Success: " << count << endl;
		UINT64 numOfReads = pSubjectData->iCountSubjectRead;
		if (Config::isNumberOfThreadsSet) {
			omp_set_dynamic(0);
			omp_set_num_threads(Config::numberOfThreads);
		}
		UINT64 numMaxThreads = omp_get_max_threads();
		vector<double> vTimeHashing(numMaxThreads, 0);
		vector<double> vTimeAligning(numMaxThreads, 0);
		vector<double> vTimeAlignmentConvertToString(numMaxThreads, 0);
#pragma omp parallel
		{
			vector<vector<PairKeyPosition *> *> * vvKeyArray = new vector<vector<PairKeyPosition *> *>();
			BandAligner * aligner = new BandAligner();
			stringstream * sstream = new stringstream();
#pragma omp for schedule(dynamic)
			for (UINT64 i = 0; i < numOfReads; i++) {
				SubjectRead * pSubjectRead = pSubjectData->vpSubjectRead->at(i);
				if (!(*pSubjectRead).bIsGoodRead) {
					continue;
				}
				double beginTime = omp_get_wtime();
				double endTime = 0;
				if (this->searchHashTable2((*pSubjectRead), AlignmentRecord::FORWARD)) {
					endTime = omp_get_wtime();
					vTimeHashing.at(omp_get_thread_num()) += (endTime - beginTime);
					beginTime = omp_get_wtime();
					this->processOneSubjectRead2((*pSubjectRead), AlignmentRecord::FORWARD, vvKeyArray, aligner);
					endTime = omp_get_wtime();
					vTimeAligning.at(omp_get_thread_num()) += (endTime - beginTime);
				}
				beginTime = omp_get_wtime();
				if (this->searchHashTable2((*pSubjectRead), AlignmentRecord::REVERSECOMPLEMENT)) {
					endTime = omp_get_wtime();
					vTimeHashing.at(omp_get_thread_num()) += (endTime - beginTime);
					beginTime = omp_get_wtime();
					this->processOneSubjectRead2((*pSubjectRead), AlignmentRecord::REVERSECOMPLEMENT, vvKeyArray,
							aligner);
					endTime = omp_get_wtime();
					vTimeAligning.at(omp_get_thread_num()) += (endTime - beginTime);
				}
				beginTime = omp_get_wtime();
				(*pSubjectRead).clearKeyPosition();
				endTime = omp_get_wtime();
				vTimeHashing.at(omp_get_thread_num()) += (endTime - beginTime);
				beginTime = omp_get_wtime();
				if (Config::bCallConsensus) {
					if (Config::bMerging) {
						(*pSubjectRead).updateQueryCoverageConsiderMerging(this->hashTable->dataSet);
					} else {
						if(Config::bMappingVersion){
							(*pSubjectRead).chooseQueryWithBestAlignment();
						}
						(*pSubjectRead).updateQueryCoverage(this->hashTable->dataSet);
					}
				} else {
					if (Config::bSpeicialVersion) {
						(*pSubjectRead).writeAlignmentLongFormat(*sstream);
					} else {
						(*pSubjectRead).writeAlignment(*sstream);
					}
				}
				endTime = omp_get_wtime();
				vTimeAlignmentConvertToString.at(omp_get_thread_num()) += (endTime - beginTime);
			}
			delete aligner;
			delete sstream;
			for (size_t iv = 0; iv < vvKeyArray->size(); iv++) {
				vvKeyArray->at(iv)->clear();
				delete vvKeyArray->at(iv);
			}
			vvKeyArray->clear();
			delete vvKeyArray;
		}
		//cout << "Aligning Success: " << count << endl;
		//dump alignment to hard drive
		beginTime = omp_get_wtime();
		if (!Config::bCallConsensus) {
			if (Config::bSpeicialVersion) {
				this->writeToFile2(pSubjectData, Config::outputfilename);
			} else {
				this->writeToFile(pSubjectData, Config::outputfilename);
			}
		}
		endTime = omp_get_wtime();
		dTimeWriting += endTime - beginTime;
		//cout << "Writing Success: " << count << endl;
		count++;
		double max = 0;
		for (size_t i = 0; i < vTimeHashing.size(); i++) {
			if (max < vTimeHashing.at(i)) {
				max = vTimeHashing.at(i);
			}
		}
		this->dTimeHashing += max;
		max = 0;
		for (size_t i = 0; i < vTimeAligning.size(); i++) {
			if (max < vTimeAligning.at(i)) {
				max = vTimeAligning.at(i);
			}
		}
		this->dTimeAligning += max;
		max = 0;
		for (size_t i = 0; i < vTimeAlignmentConvertToString.size(); i++) {
			if (max < vTimeAlignmentConvertToString.at(i)) {
				max = vTimeAlignmentConvertToString.at(i);
			}
		}
		this->dTimeAlignmentConvertToString += max;
		beginTime = omp_get_wtime(); // to record the loading time
	}
	if (Config::bWriteFinishFlag) {
		ofstream filePointer;
		char str[4];
		snprintf(str, 4, "%03u", this->iIndexShard);
		string sTemp = "";
		sTemp.append(Config::sFlagFolder);
		sTemp.append(str);
		sTemp.append(1, '_');
		sTemp.append(Config::sSuffix);
		filePointer.open(sTemp.c_str(), ios_base::app);
		filePointer.close();
	}

	MEMORYSTOP
	CLOCKSTOP

	cout << this->dTimeLoading << " seconds (Time for loading)" << endl;
	cout << this->dTimeHashing << " seconds (Time for hashing)" << endl;
	cout << this->dTimeAligning << " seconds (Time for aligning)" << endl;
	cout << this->dTimeAlignmentConvertToString << " seconds (Time for converting alignment to string or updating matrix)" << endl;
	if (!Config::bCallConsensus)
		cout << this->dTimeWriting << " seconds (Time for writing results)" << endl;

	delete pSubjectData;

	return true;
}

/*
 bool SubjectReadAligner::start() {

 CLOCKSTART
 MEMORYSTART
 int count = 1;
 SubjectDataset * pSubjectData = new SubjectDataset();
 while (pSubjectData->loadNextChunkParallel(Config::numberOfThreads)) {
 cout << "Loading Success: " << count << endl;
 UINT64 numOfReads = pSubjectData->iCountSubjectRead;
 if (Config::isNumberOfThreadsSet) {
 omp_set_dynamic(0);
 omp_set_num_threads(Config::numberOfThreads);
 }
 #pragma omp parallel
 {
 UINT64 id, i, Nthrds, istart, iend;
 id = omp_get_thread_num();
 Nthrds = omp_get_num_threads();
 istart = id * numOfReads / Nthrds;
 iend = (id + 1L) * numOfReads / Nthrds;
 if (id == (Nthrds - 1L))
 iend = numOfReads;
 vector<vector<PairKeyPosition *> *> * vvKeyArray = new vector<
 vector<PairKeyPosition *> *>();
 BandAligner * aligner = new BandAligner();
 for (i = istart; i < iend; i++) {
 if (!(*(pSubjectData->vpSubjectRead->at(i))).isGoodRead) {
 continue;
 }
 if (this->searchHashTable((*(pSubjectData->vpSubjectRead->at(i))),
 AlignmentRecord::FORWARD)) {
 this->processOneSubjectRead((*(pSubjectData->vpSubjectRead->at(i))),
 AlignmentRecord::FORWARD, vvKeyArray, aligner);
 }
 if (this->searchHashTable((*(pSubjectData->vpSubjectRead->at(i))),
 AlignmentRecord::REVERSECOMPLEMENT)) {
 this->processOneSubjectRead((*(pSubjectData->vpSubjectRead->at(i))),
 AlignmentRecord::REVERSECOMPLEMENT, vvKeyArray, aligner);
 }
 }
 delete aligner;
 for (size_t iv = 0; iv < vvKeyArray->size(); iv++) {
 vvKeyArray->at(iv)->clear();
 delete vvKeyArray->at(iv);
 }
 vvKeyArray->clear();
 delete vvKeyArray;
 }
 cout << "Aligning Success: " << count << endl;
 //dump alignment to hard drive
 this->writeToFile(pSubjectData, Config::outputfilename);
 cout << "Writinging Success: " << count << endl;
 count++;
 }

 delete pSubjectData;

 MEMORYSTOP
 CLOCKSTOP

 return true;
 }*/

INT64 calKeyLength(vector<PairKeyPosition *> * _vKeyArray) {
	INT64 length = 0;
	if (_vKeyArray->size() == 1) {
		length = Config::iKmer;
		return length;
	}
	for (size_t i = 0; i < _vKeyArray->size() - 1; i++) {
		if (_vKeyArray->at(i)->iPositionOnSubject + Config::iKmer
				<= _vKeyArray->at(i + 1)->iPositionOnSubject) {
			length += Config::iKmer;
		} else {
			length += (_vKeyArray->at(i + 1)->iPositionOnSubject - _vKeyArray->at(i)->iPositionOnSubject);
		}
	}
	length += Config::iKmer;
	return length;
}

void markTheBestPossibleAlignment(vector<vector<PairKeyPosition *> *> * _vvKeyArray) {
	INT64 longestLength = 0, temp = 0;
	for (size_t i = 0; i < _vvKeyArray->size(); i++) {
		if (_vvKeyArray->at(i)->size() == 0) {
			continue;
		} else {
			temp = calKeyLength(_vvKeyArray->at(i));
			if (temp > longestLength) {
				longestLength = temp;
				for (size_t j = 0; j < i; j++) {
					_vvKeyArray->at(j)->clear();
				}
			} else if (temp < longestLength) {
				_vvKeyArray->at(i)->clear();
			}
		}
	}
	/*if (longestLength < Config::minimumOverlapLength) {
	 for (size_t i = 0; i < _vvKeyArray->size(); i++) {
	 _vvKeyArray->at(i)->clear();
	 }
	 }*/
}

bool SubjectReadAligner::processOneSubjectRead(SubjectRead & _subjectRead, uint _iOrientation,
		vector<vector<PairKeyPosition *> *> * _vvKeyArray, BandAligner * _aligner) {
	string subjectSequence;
	if (_iOrientation == AlignmentRecord::FORWARD) {
		_subjectRead.getSequence(subjectSequence);
	} else {
		_subjectRead.getReverseComplement(subjectSequence);
	}
	for (size_t i = 0; i < _vvKeyArray->size(); i++) {
		_vvKeyArray->at(i)->clear();
	}
	bool bAdded = false;
	for (size_t i = 0; i < _subjectRead.vKeyPosition.size(); i++) {
		bAdded = false;
		size_t j = 0;
		for (j = 0; j < _vvKeyArray->size(); j++) {
			if (_vvKeyArray->at(j)->size() == 0) {
				break;
			}
			if (this->isFromTheSameAlignment(_subjectRead.vKeyPosition.at(i), _vvKeyArray->at(j)->at(0),
					_vvKeyArray->at(j)->back())) {
				_vvKeyArray->at(j)->push_back(_subjectRead.vKeyPosition.at(i));
				bAdded = true;
				break;
			}
		}
		if (!bAdded) {
			if (j >= _vvKeyArray->size()) {
				vector<PairKeyPosition *> * temp = new vector<PairKeyPosition *>();
				_vvKeyArray->push_back(temp);
			}
			_vvKeyArray->at(j)->push_back(_subjectRead.vKeyPosition.at(i));
		}

		if (i + 1 == _subjectRead.vKeyPosition.size()
				|| _subjectRead.vKeyPosition.at(i)->iQueryReadInnerId
						!= _subjectRead.vKeyPosition.at(i + 1)->iQueryReadInnerId) { // get detailed alignment
			markTheBestPossibleAlignment(_vvKeyArray);
			for (j = 0; j < _vvKeyArray->size(); j++) {
				if (_vvKeyArray->at(j)->size() != 0) {
					/*if(_subjectRead.getIdentifier()==1866&&_vvKeyArray->at(j)->at(0)->iQueryReadId==143931){
					 Config::bDebug = true;
					 }*/
					this->getAlignment(_subjectRead, subjectSequence, _iOrientation, _vvKeyArray->at(j), _aligner);
					_vvKeyArray->at(j)->clear();
				}
			}
		}
	}
	return true;
}

bool SubjectReadAligner::processOneSubjectRead2(SubjectRead & _subjectRead, uint _iOrientation,
		vector<vector<PairKeyPosition *> *> * _vvKeyArray, BandAligner * _aligner) {
	string subjectSequence;
	if (Config::bCallConsensus && _iOrientation != AlignmentRecord::FORWARD) {
		subjectSequence = _subjectRead.sReadSequenceReverseComplement;
	} else {
		_subjectRead.getSequence(subjectSequence);
	}
	for (size_t i = 0; i < _vvKeyArray->size(); i++) {
		_vvKeyArray->at(i)->clear();
	}
	bool bAdded = false;
	for (size_t i = 0; i < _subjectRead.vKeyPosition.size(); i++) {
		bAdded = false;
		size_t j = 0;
		for (j = 0; j < _vvKeyArray->size(); j++) {
			if (_vvKeyArray->at(j)->size() == 0) {
				break;
			}
			if (this->isFromTheSameAlignment(_subjectRead.vKeyPosition.at(i), _vvKeyArray->at(j)->at(0),
					_vvKeyArray->at(j)->back())) {
				_vvKeyArray->at(j)->push_back(_subjectRead.vKeyPosition.at(i));
				bAdded = true;
				break;
			}
		}
		if (!bAdded) {
			if (j >= _vvKeyArray->size()) {
				vector<PairKeyPosition *> * temp = new vector<PairKeyPosition *>();
				_vvKeyArray->push_back(temp);
			}
			_vvKeyArray->at(j)->push_back(_subjectRead.vKeyPosition.at(i));
		}

		if (i + 1 == _subjectRead.vKeyPosition.size()
				|| _subjectRead.vKeyPosition.at(i)->iQueryReadInnerId
						!= _subjectRead.vKeyPosition.at(i + 1)->iQueryReadInnerId) { // get detailed alignment
			markTheBestPossibleAlignment(_vvKeyArray);
			for (j = 0; j < _vvKeyArray->size(); j++) {
				if (_vvKeyArray->at(j)->size() != 0) {
					this->getAlignment2(_subjectRead, subjectSequence, _iOrientation, _vvKeyArray->at(j), _aligner);
					_vvKeyArray->at(j)->clear();
				}
			}
		}
	}
	return true;
}

bool SubjectReadAligner::writeToFile(SubjectDataset * _pSubjectData, string _fileName) {
	ofstream filePointer;
	string _fileNameBackup = _fileName;
	if (Config::iSizeLimit > 0) {
		_fileName.append(1, '/');
		char str[4];
		snprintf(str, 4, "%03u", this->iIndexShard);
		_fileName.append(str);
		_fileName.append(1, '_');
		_fileName.append(Config::sSuffix);
	} else if (!Config::sSuffix.empty()) {
		_fileName.append(1, '/');
		_fileName.append(Config::sSuffix);
	}
	_fileName.append(".aln");
	filePointer.open(_fileName.c_str(), ios_base::app);
	if (filePointer == NULL) {
		cout << "Unable to open file: " << endl;
		return false;
	}
	UINT64 i = 0;
	for (i = 0; i < _pSubjectData->iCountSubjectRead; i++) {
		if (Config::iSizeLimit > 0 && this->iNumberReadsProcessed >= Config::iSizeLimit) {
			filePointer.close();
			if (Config::bWriteFinishFlag) {
				char str[4];
				snprintf(str, 4, "%03u", this->iIndexShard);
				string sTemp = "";
				sTemp.append(Config::sFlagFolder);
				sTemp.append(str);
				sTemp.append(1, '_');
				sTemp.append(Config::sSuffix);
				filePointer.open(sTemp.c_str(), ios_base::app);
				filePointer.close();
			}
			this->iIndexShard++;
			_fileName = _fileNameBackup;
			_fileName.append(1, '/');
			char str[4];
			snprintf(str, 4, "%03u", this->iIndexShard);
			_fileName.append(str);
			_fileName.append(1, '_');
			_fileName.append(Config::sSuffix);
			_fileName.append(".aln");
			filePointer.open(_fileName.c_str(), ios_base::app);
			this->iNumberReadsProcessed = 0;
		}
		if (_pSubjectData->vpSubjectRead->at(i)->sAlignment.empty()) {
			filePointer << ">" << _pSubjectData->vpSubjectRead->at(i)->iReadUniqueId << endl << endl;
		} else {
			filePointer << _pSubjectData->vpSubjectRead->at(i)->sAlignment;
		}
		iNumberReadsProcessed++;
	}
	filePointer.close();
	return true;
}

bool SubjectReadAligner::writeToFile2(SubjectDataset * _pSubjectData, string _fileName) {
	ofstream filePointer;
	_fileName.append(".aln");
	filePointer.open(_fileName.c_str(), ios_base::app);
	if (filePointer == NULL) {
		cout << "Unable to open file: " << endl;
		return false;
	}
	UINT64 i = 0;
	for (i = 0; i < _pSubjectData->iCountSubjectRead; i++) {
		/*filePointer << ">";
		 filePointer << _pSubjectData->vpSubjectRead->at(i)->sReadName << "\t"
		 << _pSubjectData->vpSubjectRead->at(i)->getReadLength() << endl;
		 sstream.clear();
		 sstream.str(std::string());
		 _pSubjectData->vpSubjectRead->at(i)->getAlignmentLongFormat(sstream);
		 filePointer << sstream.str() << endl;*/
		if (_pSubjectData->vpSubjectRead->at(i)->sAlignment.empty()) {
			filePointer << ">" << _pSubjectData->vpSubjectRead->at(i)->iReadUniqueId << endl << endl;
		} else {
			filePointer << _pSubjectData->vpSubjectRead->at(i)->sAlignment;
		}
	}
	filePointer.close();
	return true;
}

