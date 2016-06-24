/*
 * main.cpp
 *
 *  Created on: Oct 6, 2015
 *      Author: xgo
 */
#include "Config.h"
#include "HashTableLongKey.h"
#include "HashTableShortKey.h"
#include "QueryDataset.h"
#include "SubjectReadAligner.h"

#define KEYLENTHROSHOLD 0

void callConcensusAndMergeAndWriteOuput(QueryDataset* queryDataset);

bool bDebug = false;
void debug();
void debugAffineGapAlignRight(BandAligner * _aligner, string & _subject, string & _query);
void debugAffineGapAlignLeft(BandAligner * _aligner, string & _subject, string & _query);
void debugAffineGap(BandAligner * _aligner, string & _subject, string & _query);
void debugGlobalAlignRight(BandAligner * _aligner, string & _subject, string & _query);
void debugGlobalAlignLeft(BandAligner * _aligner, string & _subject, string & _query);
void debugGlobal(BandAligner * _aligner, string & _subject, string & _query);

int main(int argc, char **argv) {

	double dTimePointBegin = omp_get_wtime();
	if (bDebug) {
		debug();
		return 0;
	}

	// parse command line options:
	if (!Config::setConfig(argc, argv)) {
		cout << "Please follow the above help information." << endl;
		return false;
	}

	QueryDataset* queryDataset;
	if (Config::sQueryFileName.empty()) {
		queryDataset = new QueryDataset(Config::vQueryFiles);
	} else {
		queryDataset = new QueryDataset(Config::getQueryDatasetFilename());
	}
	if (!queryDataset->buildDataset()) {
		cout << "Error: cannot build query dataset" << endl;
		return false;
	} else if (queryDataset->getNumberOfReads() == 0) {
		cout << "Unfortunately all the query data are filtered out. " << endl << "Further processing is aborted. "
				<< endl << "Please double check with the data or command settings." << endl;
		return false;
	}

	bool bDebugTemp = false;
	if (bDebugTemp) {
		return true;
	}

	HashTableMethod* hashTable;
	if (Config::iKmer <= KEYLENTHROSHOLD) {
		hashTable = new HashTableShortKey(queryDataset);
	} else {
		hashTable = new HashTableLongKey(queryDataset);
	}
	hashTable->createHashTables();
	hashTable->insertQueryDataset(queryDataset);

	SubjectReadAligner * aligner = new SubjectReadAligner(hashTable);
	aligner->start();

	callConcensusAndMergeAndWriteOuput(queryDataset);

	double dTimePointEnd = omp_get_wtime();
	cout << (dTimePointEnd - dTimePointBegin) << " seconds (Total time used)" << endl;
	cout << "DONE";
	delete queryDataset;
	delete hashTable;
	delete aligner;
}

void callConcensusAndMergeAndWriteOuput(QueryDataset* queryDataset) {
	if(Config::bMappingVersion){
		string fileName = Config::outputfilename;
		queryDataset->getCoverageAndWriteResults(fileName);
		return;
	}
	if (Config::bCallConsensus) {
		double beginTime = omp_get_wtime();
		queryDataset->callConsensus();
		double endTime = omp_get_wtime();
		cout << (endTime - beginTime) << " seconds (Time for calling consensus)" << endl;

		if(Config::bMerging){
			queryDataset->mergePairedReads();
			beginTime = omp_get_wtime();
			cout << (beginTime - endTime) << " seconds (Time for merging paired reads.)" << endl;
		}

		beginTime = omp_get_wtime();
		string fileName = Config::outputfilename + Config::sSuffix;
		queryDataset->writeQueryOut(fileName);
		endTime = omp_get_wtime();
		cout << (endTime - beginTime) << " seconds (Time for writing results)" << endl;
	}
}

void debug() {
	/*UINT64 iMask = (1L << 52) - 1;
	 cout << iMask << endl;
	 iMask = 0xFFFFFFFFFFFFF;
	 cout << iMask << endl;*/
	//BandAligner aligner;
	string sRead =
			"GAGTTGTTTGCTTTCATGCATGACAGCCTTTGTAAAGCAGCGGCAAACTTGATTTTTTCGAGGTATTCCAGGAAGATTTTGATGCCGCCGAAATTGGTCGCGTTCTGTAGGGAAAATTCGGTCTTGATTTTGCTGATGGGTGTGGTAGACTTCACCTAAAAGGTGCTCCTCTCTTTGGGTGTAGTGTTCTCGATAAACACCATTCTACCAAAGATTTGGGCACTTTTTTATGTTTTCATAGGATCGCACTTTCGAAATTCAGGATAAATATCAACATACAAAGGGCTCAACATCAAAGTTAGTGCT";
	string sRev = QueryRead::reverseComplement(sRead);
	cout << sRev << endl;
	BandAligner * aligner = new BandAligner();
	string subject = "AAGCCAGCTTGACTGCGAGCGTGACGGCGCGAGCAGGTGCGAAAGCAGGTCATAGTGATCCGGTGGTTCTGAATGGAAGGGCCATCGC";
	string query = "TTGACTGCGAGATCGACGGATCAAGCAGGTACGAAAGAAGGACTTAGTGATCCGGTGGTTCTGTATGGAAGGGCCATCGC";
	INT64 band = 9;
	aligner->setBand(band);
	debugGlobalAlignRight(aligner, subject, query);
	debugGlobal(aligner, subject, query);
	debugAffineGapAlignRight(aligner, subject, query);
	debugAffineGap(aligner, subject, query);
	/*subject = "ATTCTGATCTAGACGGCGATCCGGTGCAGGTCTTCG";
	 query = "TGGCGATCCGGTGCAGGTCTTCG";
	 debugGlobalAlignRight(aligner, subject, query);
	 debugGlobal(aligner, subject, query);
	 debugAffineGapAlignRight(aligner, subject, query);
	 debugAffineGap(aligner, subject, query);

	 subject = "ATTCTGATCTAGACGGCGATCCGGTGCAGGTCTTCG";
	 query = "TGGGGATCCGGTGCAGGTCTTCG";
	 debugGlobalAlignRight(aligner, subject, query);
	 debugGlobal(aligner, subject, query);
	 debugAffineGapAlignRight(aligner, subject, query);
	 debugAffineGap(aligner, subject, query);

	 //aligner->printMatrix();
	 subject = "ATTCTGATCTAGACGGCGATCCGGTGCAGGTCTTCG";
	 query = "ATTCTGATCTAGACGGTGCAGGTCTTCG";
	 debugGlobalAlignRight(aligner, subject, query);
	 debugGlobal(aligner, subject, query);
	 debugAffineGapAlignRight(aligner, subject, query);
	 debugAffineGap(aligner, subject, query);

	 subject = "CGGCATTCTGATCTAGACGGCGATCCGGTGCAGGTCTTCG";
	 query = "ATTCTTGACGGCGATCCGGCTTCG";
	 debugGlobalAlignRight(aligner, subject, query);
	 debugGlobal(aligner, subject, query);
	 debugAffineGapAlignRight(aligner, subject, query);
	 debugAffineGap(aligner, subject, query);*/

	//subject = "CGGCATTCTGATCTAGACGGCGATCCGGTGCA";
	//query = "CGGGATCTAGTCGGCGATC";
	debugGlobalAlignLeft(aligner, subject, query);
	debugGlobal(aligner, subject, query);
	debugAffineGapAlignLeft(aligner, subject, query);
	debugAffineGap(aligner, subject, query);

	//subject = "CGGCATTCTGATCTAGACGGCGATCCGGTGCA";
	//query = "CGGGATCTAGTCGGCTATG";
	debugGlobalAlignLeft(aligner, subject, query);
	debugGlobal(aligner, subject, query);
	debugAffineGapAlignLeft(aligner, subject, query);
	debugAffineGap(aligner, subject, query);

	//aligner->printMatrix();
	delete aligner;
}

void debugAffineGapAlignRight(BandAligner * _aligner, string & _subject, string & _query) {
	cout << "AffineGap Align Right:" << endl;
	_aligner->setSequences(_subject, 0, _subject.length(), _query, 0, _query.length());
	_aligner->setAlignRight();
	_aligner->bandedAlignAffineGap();
	_aligner->calculateCigarAndPosition();
	_aligner->printAlignment();
}

void debugAffineGapAlignLeft(BandAligner * _aligner, string & _subject, string & _query) {
	cout << "AffineGap Align Left:" << endl;
	_aligner->setSequences(_subject, 0, _subject.length(), _query, 0, _query.length());
	_aligner->setAlignLeft();
	_aligner->bandedAlignAffineGap();
	_aligner->calculateCigarAndPosition();
	_aligner->printAlignment();
}

void debugAffineGap(BandAligner * _aligner, string & _subject, string & _query) {
	cout << "AffineGap:" << endl;
	_aligner->setSequences(_subject, 0, _subject.length(), _query, 0, _query.length());
	_aligner->bandedAlignAffineGap();
	_aligner->calculateCigarAndPosition();
	_aligner->printAlignment();
}

void debugGlobalAlignRight(BandAligner * _aligner, string & _subject, string & _query) {
	cout << "Global Align Right:" << endl;
	_aligner->setSequences(_subject, 0, _subject.length(), _query, 0, _query.length());
	_aligner->setAlignRight();
	_aligner->bandedAlign();
	_aligner->calculateCigarAndPosition();
	_aligner->printAlignment();
}

void debugGlobalAlignLeft(BandAligner * _aligner, string & _subject, string & _query) {
	cout << "Global Align Left:" << endl;
	_aligner->setSequences(_subject, 0, _subject.length(), _query, 0, _query.length());
	_aligner->setAlignLeft();
	_aligner->bandedAlign();
	_aligner->calculateCigarAndPosition();
	_aligner->printAlignment();
}

void debugGlobal(BandAligner * _aligner, string & _subject, string & _query) {
	cout << "Global:" << endl;
	_aligner->setSequences(_subject, 0, _subject.length(), _query, 0, _query.length());
	_aligner->bandedAlign();
	_aligner->calculateCigarAndPosition();
	_aligner->printAlignment();
}

