/*
 * HashTableShortKey.h
 *
 *  Created on: Oct 28, 2015
 *      Author: xgo
 */

#ifndef HASHTABLESHORTKEY_H_
#define HASHTABLESHORTKEY_H_

#include "QueryDataset.h"
#include "HashTableMethod.h"

class HashTableShortKey: public HashTableMethod {
public:
	HashTableShortKey();
	HashTableShortKey(QueryDataset * qDataset);
	virtual ~HashTableShortKey();

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

#endif /* HASHTABLESHORTKEY_H_ */
