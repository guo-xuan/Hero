/*
 * QueryDataset.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYDATASET_H_
#define QUERYDATASET_H_

#include "Config.h"
#include "QueryRead.h"

class QueryDataset {

	UINT64 numberOfReads;								// Number of total reads present in the dataset.
	UINT64 numberOfUniqueReads; 						// number of unique reads in the dataset.

	UINT64 shortestReadLength;


	UINT16 dataset_minimumoverlaplength;

	string dataset_QueryFilename;
	vector<string> vDataset_QueryFilename;




	void sortReads();//in order to assign IDs and facilitate binary sequence search


public:
	UINT64 longestReadLength;
	vector<QueryRead *>* vQueryReads;
	QueryDataset(const string & QueryFilename);
	QueryDataset(const vector<string> & QueryFilename);
	~QueryDataset();
	static bool qualityFilter(string & sequence);
	bool loadDatasetParallel(const string & _QueryFilename, UINT16 _numberOfThreads);
	bool loadDatasetParallel(const vector<string> & _QueryFilename, UINT16 _numberOfThreads);
	bool buildDataset();
	UINT64 getNumberOfReads(); 						// Get the number of total reads in the database.
	QueryRead * getReadFromID(INT64 ID); 					// Find a read in the database given the ID in constant time.

	void callConsensus();
	void writeQueryOut(string & _sFilename);

	void mergePairedReads();

	void getCoverageAndWriteResults(string & _sFilename);
};

#endif /* QUERYDATASET_H_ */
