/*
 * QueryDataset.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "QueryDataset.h"

QueryDataset::QueryDataset(const string & QueryFilename) {
	numberOfReads = 0;			// Number of total reads present in the dataset.
	numberOfUniqueReads = 0; 		// number of unique reads in the dataset.

	shortestReadLength = 0xFFFFFFFFFFFFFFFF;
	longestReadLength = 0;
	vQueryReads = NULL;

	dataset_QueryFilename = QueryFilename;
	dataset_minimumoverlaplength = Config::getminimumOverlapLength();

}

QueryDataset::QueryDataset(const vector<string> & QueryFilename) {
	numberOfReads = 0;			// Number of total reads present in the dataset.
	numberOfUniqueReads = 0; 		// number of unique reads in the dataset.

	shortestReadLength = 0xFFFFFFFFFFFFFFFF;
	longestReadLength = 0;
	vQueryReads = NULL;
	dataset_QueryFilename = "";
	vDataset_QueryFilename = QueryFilename;
	dataset_minimumoverlaplength = Config::getminimumOverlapLength();
}

QueryDataset::~QueryDataset() {
	if (vQueryReads != NULL) {
		for (UINT64 i = 0; i < numberOfReads; i++) {
			delete vQueryReads->at(i);
		}
		vQueryReads->clear();
		delete vQueryReads;
	}
}

UINT64 QueryDataset::getNumberOfReads() {
	return this->numberOfReads;
}

QueryRead * QueryDataset::getReadFromID(INT64 ID) {
	if (vQueryReads != NULL) {
		if (ID < 0 || ID >= (INT64) numberOfReads) { // ID outside the range.
			cout << "ID " << ID << " out of bound." << endl;
			return NULL;
		} else
			return vQueryReads->at(ID);
	} else
		return NULL;
}

bool QueryDataset::qualityFilter(string & sequence) {

	//Returns true if the read contains only {A,C,G,T} and does not contain more than 80% of the same nucleotide
	UINT64 cnt[4] = { 0, 0, 0, 0 };
	UINT64 readLength = sequence.length();
	for (UINT64 i = 0; i < readLength; i++) { // Count the number of A's, C's , G's and T's in the string.
		if (sequence[i] != 'A' && sequence[i] != 'C' && sequence[i] != 'G' && sequence[i] != 'T')
			return false;
		cnt[(sequence[i] >> 1) & 0X03]++; // Least significant 2nd and 3rd bits of ASCII value used here
	}
	UINT64 threshold = sequence.length() * .8;	// 80% of the length.
	if (cnt[0] >= threshold || cnt[1] >= threshold || cnt[2] >= threshold || cnt[3] >= threshold)
		return false;	// If 80% bases are the same base.
	return true;

}

bool compareReads(QueryRead *read1, QueryRead *read2) {
	return read1->getSequence() < read2->getSequence();
}

void QueryDataset::sortReads() {
	if (vQueryReads != NULL)
		std::sort(vQueryReads->begin(), vQueryReads->end(), compareReads);	// Sort the reads lexicographically.
	else
		cout << "no reads in the dataset";

}

bool QueryDataset::buildDataset() {
	if (this->dataset_QueryFilename.empty()) {
		cout << "Loading multiple query files" << endl;
		return this->loadDatasetParallel(this->vDataset_QueryFilename, Config::numberOfThreads);
	}
	cout << "Loading one query file" << endl;
	return this->loadDatasetParallel(this->dataset_QueryFilename, Config::numberOfThreads);
	/*
	 if (Config::numberOfThreads == 1)
	 return this->buildDataset(this->dataset_QueryFilename);
	 else
	 return this->buildDatasetParallel(this->dataset_QueryFilename, Config::numberOfThreads);*/
}

bool QueryDataset::loadDatasetParallel(const string & QueryFilename, UINT16 numberOfThreads) {
	CLOCKSTART
	MEMORYSTART

	cout << "Reading data set: " << " from file: " << QueryFilename << endl;
	ifstream myFile;
	myFile.open(QueryFilename.c_str());
	if (!myFile.is_open()) {
		cout << "Unable to open file: " << QueryFilename << endl;
		return false;
	}
	enum FileType {
		FASTA, FASTQ, UNDEFINED
	};
	FileType fileType = UNDEFINED;
	int numberOfReads = 0;
	if (fileType == UNDEFINED) {
		string text;
		getline(myFile, text);
		if (text[0] == '>')
			fileType = FASTA;
		else if (text[0] == '@')
			fileType = FASTQ;
		else
			cout << "Unknown input file format." << endl;
		myFile.seekg(0, ios::beg);
		while (!myFile.eof()) {
			if (fileType == FASTA) {
				getline(myFile, text);
				getline(myFile, text, '>');
			} else if (fileType == FASTQ) {
				for (UINT64 i = 0; i < 4; i++) {
					getline(myFile, text);
				}
			}
			numberOfReads++;
		}
		if (text == "") {
			numberOfReads--;
		}
		myFile.close();
		myFile.open(QueryFilename.c_str());
		if (myFile == NULL) {
			cout << "Unable to open file: " << QueryFilename << endl;
			return false;
		}
	}
	if (vQueryReads == NULL)
		vQueryReads = new vector<QueryRead *>(numberOfReads);
	vector<QueryRead *> * thisqueryReadList = this->vQueryReads;
	this->numberOfReads = numberOfReads;
	cout << "Number of reads in query:  " << numberOfReads << endl;
	UINT64 lineOffset;
	if (fileType == FASTA)
		lineOffset = 2;
	else if (fileType == FASTQ)
		lineOffset = 4;

	UINT64 indexAnchor = 0;
	while (!myFile.eof()) {
		UINT64 chunkCounter;
		vector<string>* lineData = new vector<string>();
		string text = "";
		for (chunkCounter = 0; !myFile.eof() && chunkCounter < Config::queryChunkSize; chunkCounter++) {
			if (fileType == FASTA) {			// Fasta file
				getline(myFile, text);
				lineData->push_back(text);
				getline(myFile, text, '>');
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				lineData->push_back(text);
			} else if (fileType == FASTQ) {			// Fastq file.
				for (UINT64 i = 0; i < 4; i++) {// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
					getline(myFile, text);
					lineData->push_back(text);
				}
			}
		}
		if (text == "") {
			if (fileType == FASTA) {
				lineData->pop_back();
				lineData->pop_back();
			} else if (fileType == FASTQ) {
				for (UINT64 i = 0; i < 4; i++)
					lineData->pop_back();
			}
			chunkCounter--;
		}
		if (Config::isNumberOfThreadsSet) {
			omp_set_dynamic(0);
			omp_set_num_threads(numberOfThreads);
		}
		//cout << indexAnchor << ", " << chunkCounter << endl;
#pragma omp parallel shared(thisqueryReadList)
		{
#pragma omp for schedule(dynamic)
			for (UINT64 i = 0; i < chunkCounter; i++) {
				string line0, line1;
				line0 = lineData->at(i * lineOffset);
				line1 = lineData->at(i * lineOffset + 1);
				string readname;
				if (line0[0] == '>' || line0[0] == '@')
					readname = line0.substr(1);
				else {
					readname = line0;
				}
				for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
					*p = toupper(*p);
				QueryRead *r1 = new QueryRead();
				if (Config::bSpeicialVersion) {
					r1->sReadName = readname;
					r1->setReadUniqueId(indexAnchor + i);
				} else {
					r1->setReadUniqueId(atoi(readname.c_str()));
				}
				r1->setSequence(line1);
				if (line1.length() > this->longestReadLength) {
					this->longestReadLength = line1.length();
				}
				r1->setReadInnerId(indexAnchor + i);
				thisqueryReadList->at(indexAnchor + i) = r1;
				if (Config::bCallConsensus) {
					r1->initilizeCoverageMatrix();
				}
				if (!(line1.length() > dataset_minimumoverlaplength && qualityFilter(line1))) { // Test the read is of good quality.
					r1->bIsGoodRead = false;
				}
			}
		}
		indexAnchor += chunkCounter;
		cout << "Number of reads has been read in: " << indexAnchor << "\r";
		lineData->clear();
		delete lineData;
	}
	myFile.close();
	cout << endl;
//sort the reads by its alphabetical order; and remove the duplicate reads so then the readID can be assigned uniquely
//sortReads();
//duplicateFilter();
	MEMORYSTOP
	CLOCKSTOP
	return true;
}

bool QueryDataset::loadDatasetParallel(const vector<string> & QueryFilename, UINT16 numberOfThreads) {
	CLOCKSTART
	MEMORYSTART
	ifstream myFile;
	enum FileType {
		FASTA, FASTQ, UNDEFINED
	};
	int numberOfReads = 0;
	FileType fileType = UNDEFINED;
	for (size_t i = 0; i < QueryFilename.size(); i++) {
		cout << "Reading data set: " << " from file: " << QueryFilename.at(i) << endl;
		myFile.open(QueryFilename.at(i).c_str());
		if (!myFile.is_open()) {
			cout << "Unable to open file: " << QueryFilename.at(i) << endl;
			return false;
		}
		string text;
		if (fileType == UNDEFINED) {
			getline(myFile, text);
			if (text[0] == '>')
				fileType = FASTA;
			else if (text[0] == '@')
				fileType = FASTQ;
			else
				cout << "Unknown input file format." << endl;
		}
		myFile.seekg(0, ios::beg);
		while (!myFile.eof()) {
			if (fileType == FASTA) {
				getline(myFile, text);
				getline(myFile, text, '>');
			} else if (fileType == FASTQ) {
				for (UINT64 i = 0; i < 4; i++) {
					getline(myFile, text);
				}
			}
			numberOfReads++;
		}
		if (text == "") {
			numberOfReads--;
		}
		myFile.close();
	}
	if (vQueryReads == NULL)
		vQueryReads = new vector<QueryRead *>(numberOfReads);
	vector<QueryRead *> * thisqueryReadList = this->vQueryReads;
	this->numberOfReads = numberOfReads;
	cout << "Number of reads in query:  " << numberOfReads << endl;
	UINT64 lineOffset;
	UINT64 indexAnchor = 0;
	for (size_t i = 0; i < QueryFilename.size(); i++) {
		myFile.open(QueryFilename.at(i).c_str());
		if (myFile == NULL) {
			cout << "Unable to open file: " << QueryFilename.at(i) << endl;
			return false;
		}
		if (fileType == FASTA)
			lineOffset = 2;
		else if (fileType == FASTQ)
			lineOffset = 4;

		while (!myFile.eof()) {
			UINT64 chunkCounter;
			vector<string>* lineData = new vector<string>();
			string text = "";
			for (chunkCounter = 0; !myFile.eof() && chunkCounter < Config::queryChunkSize; chunkCounter++) {
				if (fileType == FASTA) {			// Fasta file
					getline(myFile, text);
					lineData->push_back(text);
					getline(myFile, text, '>');
					text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
					lineData->push_back(text);
				} else if (fileType == FASTQ) {			// Fastq file.
					for (UINT64 i = 0; i < 4; i++) {// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
						getline(myFile, text);
						lineData->push_back(text);
					}
				}
			}
			if (text == "") {
				if (fileType == FASTA) {
					lineData->pop_back();
					lineData->pop_back();
				} else if (fileType == FASTQ) {
					for (UINT64 i = 0; i < 4; i++)
						lineData->pop_back();
				}
				chunkCounter--;
			}
			if (Config::isNumberOfThreadsSet) {
				omp_set_dynamic(0);
				omp_set_num_threads(numberOfThreads);
			}
			//cout << indexAnchor << ", " << chunkCounter << endl;
#pragma omp parallel shared(thisqueryReadList)
			{
#pragma omp for schedule(dynamic)
				for (UINT64 i = 0; i < chunkCounter; i++) {
					string line0, line1;
					line0 = lineData->at(i * lineOffset);
					line1 = lineData->at(i * lineOffset + 1);
					string readname;
					if (line0[0] == '>' || line0[0] == '@')
						readname = line0.substr(1);
					else {
						readname = line0;
					}
					for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
						*p = toupper(*p);
					QueryRead *r1 = new QueryRead();
					if (Config::bSpeicialVersion || Config::bMappingVersion) {
						r1->sReadName = readname;
						r1->setReadUniqueId(indexAnchor + i);
					} else {
						r1->setReadUniqueId(atoi(readname.c_str()));
					}
					r1->setSequence(line1);
					if (line1.length() > this->longestReadLength) {
						this->longestReadLength = line1.length();
					}
					r1->setReadInnerId(indexAnchor + i);
					thisqueryReadList->at(indexAnchor + i) = r1;
					if (Config::bCallConsensus || Config::bMappingVersion) {
						r1->initilizeCoverageMatrix();
					}
					if (!(line1.length() > dataset_minimumoverlaplength && qualityFilter(line1))) { // Test the read is of good quality.
						r1->bIsGoodRead = false;
					}
				}
			}
			indexAnchor += chunkCounter;
			cout << "Number of reads has been read in: " << indexAnchor << "\r";
			lineData->clear();
			delete lineData;
		}
		myFile.close();
	}

	cout << endl;
//sort the reads by its alphabetical order; and remove the duplicate reads so then the readID can be assigned uniquely
//sortReads();
//duplicateFilter();
	MEMORYSTOP
	CLOCKSTOP
	return true;
}

void QueryDataset::callConsensus() {
	if (Config::isNumberOfThreadsSet) {
		omp_set_dynamic(0);
		omp_set_num_threads(Config::numberOfThreads);
	}
	UINT64 iNumIteration = this->vQueryReads->size();
#pragma omp parallel
	{
		vector<UINT32> * vs = new vector<UINT32>();
#pragma omp for schedule(dynamic)
		for (UINT64 i = 0; i < iNumIteration; i++) {
			if (!vQueryReads->at(i)->bIsGoodRead) {
				continue;
			}
			vQueryReads->at(i)->callConsensus();
			vQueryReads->at(i)->calculateCoverageStatistic(vs);
		}
		delete vs;
	}
}

void QueryDataset::writeQueryOut(string & _sFilename) {
	string sFile;
	ofstream filePointer;
	if (Config::bMerging) {
		INT64 iNumReadsMerged = 0;
		INT64 iNumReadsUnmerged = 0;
		sFile = _sFilename + "_merged.fasta";
		filePointer.open(sFile.c_str(), ios_base::app);
		ofstream filePointerUnmerged;
		sFile = _sFilename + "_unmerged.fasta";
		filePointerUnmerged.open(sFile.c_str(), ios_base::app);
		for (size_t i = 0; i < this->vQueryReads->size(); i += 2) {
			if (!vQueryReads->at(i)->isGoodRead()) {
				if (vQueryReads->at(i + 1)->isGoodRead()) {
					filePointer << ">";
					filePointer << this->vQueryReads->at(i + 1)->iReadUniqueId << endl;
					filePointer << this->vQueryReads->at(i + 1)->sReadSequence << endl;
					iNumReadsUnmerged++;
				}
			} else {
				if (!vQueryReads->at(i + 1)->isGoodRead()) {
					filePointer << ">";
					filePointer << this->vQueryReads->at(i)->iReadUniqueId << endl;
					filePointer << this->vQueryReads->at(i)->sReadSequence << endl;
					iNumReadsUnmerged++;
				} else {
					if (vQueryReads->at(i)->bMerged && vQueryReads->at(i + 1)->bMerged) {
						filePointer << ">";
						filePointer << this->vQueryReads->at(i)->iReadUniqueId << endl;
						filePointer << this->vQueryReads->at(i)->sReadSequence << endl;
						iNumReadsMerged += 2;
					} else {
						filePointerUnmerged << ">";
						filePointerUnmerged << this->vQueryReads->at(i)->iReadUniqueId << endl;
						filePointerUnmerged << this->vQueryReads->at(i)->sReadSequence << endl;
						filePointerUnmerged << ">";
						filePointerUnmerged << this->vQueryReads->at(i + 1)->iReadUniqueId << endl;
						filePointerUnmerged << this->vQueryReads->at(i + 1)->sReadSequence << endl;
						iNumReadsUnmerged += 2;
					}
				}
			}
		}
		cout << "Merged Reads: " << iNumReadsMerged << endl;
		cout << "Not merged Reads: " << iNumReadsUnmerged << endl;
		filePointer.close();
		filePointerUnmerged.close();
	} else {
		sFile = _sFilename + ".fasta";
		filePointer.open(sFile.c_str(), ios_base::app);
		if (Config::bMergedReads) {
			for (size_t i = 0; i < this->vQueryReads->size(); i++) {
				if (!vQueryReads->at(i)->isGoodRead()) {
					continue;
				}
				filePointer << ">";
				filePointer << this->vQueryReads->at(i)->iReadUniqueId << endl;
				filePointer << this->vQueryReads->at(i)->sReadSequence << endl;
			}
		} else {
			for (size_t i = 0; i < this->vQueryReads->size(); i += 2) {
				if (vQueryReads->at(i)->isGoodRead() && vQueryReads->at(i + 1)->isGoodRead()) {
					filePointer << ">";
					filePointer << this->vQueryReads->at(i)->iReadUniqueId << endl;
					filePointer << this->vQueryReads->at(i)->sReadSequence << endl;
					filePointer << ">";
					filePointer << this->vQueryReads->at(i + 1)->iReadUniqueId << endl;
					filePointer << this->vQueryReads->at(i + 1)->sReadSequence << endl;
				}
			}
		}

		filePointer.close();
	}

	if (Config::bCoverageStatistic) {
		sFile = _sFilename + ".info";
		filePointer.open(sFile.c_str(), ios_base::app);
		if (Config::bMergedReads) {
			for (size_t i = 0; i < this->vQueryReads->size(); i++) {
				if (!vQueryReads->at(i)->isGoodRead()) {
					continue;
				}
				filePointer << this->vQueryReads->at(i)->iReadUniqueId << "\t";
				filePointer << this->vQueryReads->at(i)->iMin << "\t";
				filePointer << this->vQueryReads->at(i)->iMed << "\t";
				filePointer << this->vQueryReads->at(i)->iMax << "\t";
				filePointer << this->vQueryReads->at(i)->getRawReadLength();
				if (Config::bRealMismatchRate) {
					filePointer << "\t" << this->vQueryReads->at(i)->dActualMismatchRate;
				}
				filePointer << "\t" << this->vQueryReads->at(i)->dGcContent;
				filePointer << endl;
			}
			filePointer.close();
		} else {
			for (size_t i = 0; i < this->vQueryReads->size(); i += 2) {
				if (vQueryReads->at(i)->isGoodRead() && vQueryReads->at(i + 1)->isGoodRead()) {
					filePointer << this->vQueryReads->at(i)->iReadUniqueId << "\t";
					filePointer << this->vQueryReads->at(i)->iMin << "\t";
					filePointer << this->vQueryReads->at(i)->iMed << "\t";
					filePointer << this->vQueryReads->at(i)->iMax << "\t";
					filePointer << this->vQueryReads->at(i)->getRawReadLength();
					if (Config::bRealMismatchRate) {
						filePointer << "\t" << this->vQueryReads->at(i)->dActualMismatchRate;
					}
					filePointer << "\t" << this->vQueryReads->at(i)->dGcContent;
					filePointer << endl;
					filePointer << this->vQueryReads->at(i + 1)->iReadUniqueId << "\t";
					filePointer << this->vQueryReads->at(i + 1)->iMin << "\t";
					filePointer << this->vQueryReads->at(i + 1)->iMed << "\t";
					filePointer << this->vQueryReads->at(i + 1)->iMax << "\t";
					filePointer << this->vQueryReads->at(i + 1)->getRawReadLength();
					if (Config::bRealMismatchRate) {
						filePointer << "\t" << this->vQueryReads->at(i + 1)->dActualMismatchRate;
					}
					filePointer << "\t" << this->vQueryReads->at(i + 1)->dGcContent;
					filePointer << endl;
				}
			}
			filePointer.close();
		}
	}
}

void QueryDataset::mergePairedReads() {
	if (Config::isNumberOfThreadsSet) {
		omp_set_dynamic(0);
		omp_set_num_threads(Config::numberOfThreads);
		//omp_set_num_threads(1);
	}
	UINT64 iNumIteration = this->vQueryReads->size();
#pragma omp parallel
	{
		string sequenceLeft = "";
		string sequenceRight = "";
#pragma omp for schedule(dynamic)
		for (UINT64 i = 0; i < iNumIteration; i += 2) {
			if ((!vQueryReads->at(i)->isGoodRead()) || (!vQueryReads->at(i + 1)->isGoodRead())) {
				continue;
			}
			vQueryReads->at(i)->mergingPairedReads(vQueryReads->at(i + 1), sequenceLeft, sequenceRight);
		}
	}
}

void QueryDataset::getCoverageAndWriteResults(string & _sFilename) {
	UINT32 dMaxCoverage = 0;
	UINT32 dTemp = 0;
	INT64 j, iReadLength;
	for (size_t i = 0; i < this->vQueryReads->size(); i++) {
		iReadLength = this->vQueryReads->at(i)->getRawReadLength();
		for (j = 0; j < iReadLength; j++) {
			dTemp = this->vQueryReads->at(i)->getCoverageAt(j);
			if (dTemp > dMaxCoverage) {
				dMaxCoverage = dTemp;
			}
		}
	}
	dMaxCoverage++;
	vector<UINT32> vCoverageBp(dMaxCoverage);
	for (UINT32 i = 0; i < dMaxCoverage; i++) {
		vCoverageBp.at(i) = 0;
	}
	//UINT32 iCheck = 0;
	for (size_t i = 0; i < this->vQueryReads->size(); i++) {
		iReadLength = this->vQueryReads->at(i)->getRawReadLength();
		//cout << iReadLength << endl;
		for (j = 0; j < iReadLength; j++) {
			dTemp = this->vQueryReads->at(i)->getCoverageAt(j);
			vCoverageBp.at(dTemp) += 1;
			//iCheck++;
		}
	}
	//cout << iCheck << endl;
	//write out the results
	ofstream filePointer;
	filePointer.open(_sFilename.c_str(), ios_base::app);
	filePointer << "COV\tFRE" << endl;
	for (UINT32 i = 0; i < dMaxCoverage; i++) {
		if (vCoverageBp.at(i) > 0) {
			filePointer << i << "\t" << vCoverageBp.at(i) << endl;
		}
	}
	filePointer.flush();
	filePointer.close();
}
