/*
 * HashTable.h
 *
 *  Created on: Oct 6, 2015
 *      Author: xgo
 */

#ifndef HASHTABLELONGKEY_H_
#define HASHTABLELONGKEY_H_

#include "QueryDataset.h"
#include "HashTableMethod.h"

class HashTableLongKey : public HashTableMethod{

public:
	/*QueryDataset *dataSet;

	UINT64 hashTableSize; 		// Size of hash table, primer number usually.
	vector<vector<UINT64> *> *vviHashTable; // Main structure of hashtable, storing read identifiers.
	UINT8 iBitOfReadId;
	UINT64 iMask;

	UINT16 hashKeyLength;

	UINT64 numberOfHashCollision;// Counted total number of hash collisions. For debugging only.
	UINT64 maxSingleHashCollision;// Counted maximal number of hash collision for a single case.*/

	HashTableLongKey();
	HashTableLongKey(QueryDataset * qDataset);
	~HashTableLongKey();
	bool createHashTables();
	vector<INT64> * getListOfReads(const string & _sequence, INT64 _start, INT64 _length);
	UINT64 getPrimeLargerThanNumber(UINT64 number);
	UINT64 hashFunction(const string & _subString);
	UINT64 hashFunction(const string & _sequence, INT64 _start, INT64 _length);
	bool insertQueryDataset(QueryDataset* _qDataset);
	bool insertQueryRead(QueryRead *read, string subString, INT64 position);
	bool isEmptyAt(INT64 _hashTableIndex);
	void setHashTableSizeAndInitialize(INT64 size);
};

#endif /* HASHTABLELONGKEY_H_ */
