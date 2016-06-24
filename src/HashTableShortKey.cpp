/*
 * HashTableShortKey.cpp
 *
 *  Created on: Oct 28, 2015
 *      Author: xgo
 */

#include "HashTableShortKey.h"

HashTableShortKey::HashTableShortKey(QueryDataset * qDataset) {
	this->hashKeyLength = Config::iKmer;
	this->dataSet = qDataset;
	this->vviHashTable = NULL;
	this->hashTableSize = 0;
	this->iBitOfReadId = 52;
	this->iMask = (1L << this->iBitOfReadId) - 1;
	this->maxSingleHashCollision = 0;
	this->numberOfHashCollision = 0;

}

HashTableShortKey::~HashTableShortKey() {
	if (this->vviHashTable != NULL) {
		for (size_t i = 0; i < this->vviHashTable->size(); i++) {
			vector<INT64> * temp = this->vviHashTable->at(i);
			if (temp != NULL) {
				temp->clear();
				delete temp;
			}
		}

		this->vviHashTable->clear();
		delete this->vviHashTable;
		this->vviHashTable = NULL;
	}
}

bool HashTableShortKey::createHashTables() {

	if (this->dataSet == NULL) {
		cout << "no data set" << endl;
		return false;
	} else if (this->vviHashTable != NULL) {
		cout << "Hash Table already exists." << endl;
		return false;
	} else {
		CLOCKSTART
		MEMORYSTART
		hashTableSize = pow(4ull, hashKeyLength);
		setHashTableSizeAndInitialize(hashTableSize);
		MEMORYSTOP
		CLOCKSTOP
		return true;
	}
}

vector<INT64> * HashTableShortKey::getListOfReads(const string & _sequence, INT64 _start,
		INT64 _length) {
	INT64 index = -1;
	index = this->hashFunction(_sequence, _start, _length);
	return this->vviHashTable->at(index);
}

UINT64 HashTableShortKey::getPrimeLargerThanNumber(UINT64 number) {
	return 0;
}

UINT64 HashTableShortKey::hashFunction(const string & _subString) {
	UINT64 hashValue = 0;
	UINT64 length = _subString.length();
	for (UINT64 i = 0; i < length; i++) { // We take the bit representation of the string. A = 00, C = 01, G = 11 and T = 10
		// sum1 is for the first 32 bp. bit shifted to left 2 bits.
		// Change the character to integer. A=Ox41=01000001
		//                                  C=0x43=01000011
		//                                  G=0x47=01000111
		//                                  T=0x54=01010100
		// Then, shift to right way 1 bit.
		// Then, bit and operation with 00000011
		// Then, it just have A=00, C=01, G=11,T=10
		hashValue = (hashValue << 2ull);
		hashValue += (((UINT64) (_subString.at(i)) >> 1) & 0X03ull);
	}
	return hashValue;
}

UINT64 HashTableShortKey::hashFunction(const string & _sequence, INT64 _start, INT64 _length) {
	UINT64 hashValue = 0;
	UINT64 length = _length; // initial value of 1 to avoid zero when doing the multiplication.
	for (UINT64 i = 0; i < length; i++) { // We take the bit representation of the string. A = 00, C = 01, G = 11 and T = 10
		// sum1 is for the first 32 bp. bit shifted to left 2 bits.
		// Change the character to integer. A=Ox41=01000001
		//                                  C=0x43=01000011
		//                                  G=0x47=01000111
		//                                  T=0x54=01010100
		// Then, shift to right way 1 bit.
		// Then, bit and operation with 00000011
		// Then, it just have A=00, C=01, G=11,T=10
		hashValue = (hashValue << 2ull);
		hashValue += (((UINT64) (_sequence.at(i + _start)) >> 1) & 0X03ull);
	}
	return hashValue;
}

bool HashTableShortKey::insertQueryDataset(QueryDataset* _qDataset) {

	CLOCKSTART
	MEMORYSTART

	if (this->vviHashTable == NULL) {
		cout << "Hash Table hasn't been created yet." << endl;
		return false;
	} else {
		UINT64 datasetsize = this->dataSet->getNumberOfReads();
		UINT64 currentID = 0;
		while (currentID < datasetsize) {
			if ((currentID + 1L) % 1000000L == 0)
				cout << "Number of reads inserted in the hash table: " << (currentID + 1L) << "\r";
			QueryRead * read = this->dataSet->getReadFromID(currentID);
			if (read == NULL) {
				cout << endl << "Empty read: " << currentID << " ";
				currentID++;
				continue;
			}
			if (!read->bIsGoodRead) {
				currentID++;
				continue;
			}
			string * forwardRead = read->getSequence();
			uint i = 0;
			for (i = 0; i <= read->getRawReadLength() - this->hashKeyLength; i +=
					this->hashKeyLength) {
				string subString = forwardRead->substr(i, this->hashKeyLength);
				this->insertQueryRead(read, subString, i);
			}
			if (i < read->getRawReadLength()) {
				i = read->getRawReadLength() - this->hashKeyLength;
				string subString = forwardRead->substr(i, this->hashKeyLength);
				this->insertQueryRead(read, subString, i);
			}
			currentID++;
		}
		MEMORYSTOP
		CLOCKSTOP

		return true;
	}

	return false;
}

bool HashTableShortKey::insertQueryRead(QueryRead *read, string subString, INT64 postion) {

	UINT64 index = this->hashFunction(subString);
	if (this->vviHashTable->at(index) == NULL) {
		vector<INT64> * newList = new vector<INT64>;
		newList->resize(newList->size());
		this->vviHashTable->at(index) = newList;
	}
	INT64 positionShifted = (postion << this->iBitOfReadId);
	INT64 iReadInnerId = read->getReadInnerId();
	INT64 key = positionShifted | iReadInnerId;
	this->vviHashTable->at(index)->push_back(key); // Add the string in the list.
	return true;
}

bool HashTableShortKey::isEmptyAt(INT64 _hashTableIndex) {
	if (this->vviHashTable->at(_hashTableIndex) == NULL) {
		return true;
	}
	return this->vviHashTable->at(_hashTableIndex)->empty();
}

void HashTableShortKey::setHashTableSizeAndInitialize(INT64 size) {
	cout << "Hash Table size set to: " << size << endl;
	this->vviHashTable = new vector<vector<INT64> *>(size, NULL);
}
