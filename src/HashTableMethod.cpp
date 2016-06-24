/*
 * HashTableMethod.cpp
 *
 *  Created on: Oct 28, 2015
 *      Author: xgo
 */

#include "HashTableMethod.h"

HashTableMethod::HashTableMethod() {
	this->hashKeyLength = Config::iKmer;
	this->vviHashTable = NULL;
	this->hashTableSize = 0;
	this->iBitOfReadId = 52ull;
	this->iMask = (1ull << this->iBitOfReadId) - 1;
	this->maxSingleHashCollision = 0;
	this->numberOfHashCollision = 0;
	this->dataSet = NULL;
}

HashTableMethod::~HashTableMethod() {
	// TODO Auto-generated destructor stub
}

