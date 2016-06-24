/*
 * HashTableMethod.h
 *
 *  Created on: Oct 28, 2015
 *      Author: xgo
 */

#ifndef HASHTABLEMETHOD_H_
#define HASHTABLEMETHOD_H_

#include "QueryDataset.h"

class HashTableMethod {
public:
	QueryDataset *dataSet;

	UINT64 hashTableSize; 		// Size of hash table, primer number usually.
	vector<vector<INT64> *> *vviHashTable; // Main structure of hashtable, storing read identifiers.
	INT64 iBitOfReadId;
	INT64 iMask;

	INT64 hashKeyLength;

	INT64 numberOfHashCollision; // Counted total number of hash collisions. For debugging only.
	UINT64 maxSingleHashCollision; // Counted maximal number of hash collision for a single case.

	HashTableMethod();
	virtual ~HashTableMethod()=0;

	virtual bool createHashTables()=0;
	virtual vector<INT64> * getListOfReads(const string & _sequence, INT64 _start, INT64 _length)=0;
	virtual UINT64 getPrimeLargerThanNumber(UINT64 number)=0;
	virtual UINT64 hashFunction(const string & _subString)=0;
	virtual UINT64 hashFunction(const string & _sequence, INT64 _start, INT64 _length)=0;
	virtual bool insertQueryDataset(QueryDataset* _qDataset)=0;
	virtual bool insertQueryRead(QueryRead *read, string subString, INT64 position)=0;
	virtual bool isEmptyAt(INT64 _hashTableIndex)=0;
	virtual void setHashTableSizeAndInitialize(INT64 size)=0;

};

#endif /* HASHTABLEMETHOD_H_ */
