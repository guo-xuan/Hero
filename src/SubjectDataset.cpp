/*
 * SubjectDataset.cpp
 *
 *  Created on: Oct 16, 2015
 *      Author: xgo
 */

#include "SubjectDataset.h"

SubjectDataset::SubjectDataset() {
	this->vpSubjectRead = new vector<SubjectRead *>(Config::streamChunkSize, NULL);
	this->iInputFileIndex = 0;
	this->iCountSubjectRead = 0;
	this->eFileType = UNDEFINED;
	this->iReadIdAnchor = 0;
	sInputFileStreamer.open(Config::vSubjectFiles.at(iInputFileIndex).c_str());
	if (eFileType == UNDEFINED) {
		string text;
		getline(sInputFileStreamer, text);
		if (text[0] == '>')
			eFileType = FASTA;
		else if (text[0] == '@')
			eFileType = FASTQ;
		else
			cout << "Unknown input file format." << endl;
		sInputFileStreamer.seekg(0, ios::beg);
	}
	this->iLongestReadLength = 0;
}

SubjectDataset::~SubjectDataset() {
	this->clearSubjectRead();
	vpSubjectRead->clear();
	delete vpSubjectRead;
}

bool SubjectDataset::clearSubjectRead() {
	for (size_t i = 0; i < this->iCountSubjectRead; i++) {
		delete vpSubjectRead->at(i);
	}
	this->iCountSubjectRead = 0;
	return true;
}

bool SubjectDataset::loadNextChunkParallel(UINT16 _numberOfThreads) {
	this->iLongestReadLength = 0;
	this->clearSubjectRead();
	if (this->iInputFileIndex >= Config::vSubjectFiles.size() && this->sInputFileStreamer.eof()) {
		if (this->sInputFileStreamer.is_open())
			this->sInputFileStreamer.close();
		return false;
	}
	UINT64 lineOffset;
	if (!this->sInputFileStreamer.is_open()) {
		sInputFileStreamer.open(Config::vSubjectFiles.at(iInputFileIndex).c_str());
		if (!this->sInputFileStreamer.is_open()) {
			cout << "Unable to open file: 1 " << sInputFileStreamer << endl;
			exit(1);
		}
		if (eFileType == UNDEFINED) {
			string text;
			getline(sInputFileStreamer, text);
			if (text[0] == '>')
				eFileType = FASTA;
			else if (text[0] == '@')
				eFileType = FASTQ;
			else
				cout << "Unknown input file format." << endl;
			sInputFileStreamer.seekg(0, ios::beg);
		}
	}
	if (eFileType == FASTA)
		lineOffset = 2;
	else if (eFileType == FASTQ)
		lineOffset = 4;
	vector<string>* lineData = new vector<string>();
	UINT64 chunkCounter;
	string text = "";
	for (chunkCounter = 0; chunkCounter < Config::streamChunkSize; chunkCounter++) {
		if (sInputFileStreamer.eof()) {
			sInputFileStreamer.close();
			iInputFileIndex++;
			if (iInputFileIndex >= Config::vSubjectFiles.size()) {
				break;
			} else {
				sInputFileStreamer.open(Config::vSubjectFiles.at(iInputFileIndex).c_str());
				if (!this->sInputFileStreamer.is_open()) {
					cout << "Unable to open file: 2 " << sInputFileStreamer << endl;
					exit(1);
				}
				if (eFileType == UNDEFINED) {
					string text;
					getline(sInputFileStreamer, text);
					if (text[0] == '>')
						eFileType = FASTA;
					else if (text[0] == '@')
						eFileType = FASTQ;
					else
						cout << "Unknown input file format." << endl;
					sInputFileStreamer.seekg(0, ios::beg);
				}
				if (eFileType == FASTA)
					lineOffset = 2;
				else if (eFileType == FASTQ)
					lineOffset = 4;
			}
		}

		if (eFileType == FASTA) { // Fasta file
			getline(sInputFileStreamer, text);
			lineData->push_back(text);
			getline(sInputFileStreamer, text, '>');
			text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
			lineData->push_back(text);
		} else if (eFileType == FASTQ) { // Fastq file.
			for (UINT64 i = 0; i < 4; i++) { // Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
				getline(sInputFileStreamer, text);
				lineData->push_back(text);
			}
		}
	}
	if (text == "") {
		chunkCounter--;
	}
	//cout << "Number of reads in line data: " << chunkCounter << endl;
	if (Config::isNumberOfThreadsSet) {
		omp_set_dynamic(0);
		omp_set_num_threads(_numberOfThreads);
	}
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for (UINT64 i = 0; i < chunkCounter; i++) {
			string line0, line1;
			line0 = lineData->at(i * lineOffset);
			line1 = lineData->at(i * lineOffset + 1);
			//remove the > if it appears in the title or name
			string readname = "";
			if (line0[0] == '>' || line0[0] == '@')
				readname = line0.substr(1);
			else
				readname = line0;

			for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
				*p = toupper(*p);
			SubjectRead * subjectRead = new SubjectRead();
			subjectRead->setSequence(line1);
			if (line1.length() > this->iLongestReadLength) {
				this->iLongestReadLength = line1.length();
			}
			if (Config::bSpeicialVersion) {
				subjectRead->sReadName = readname;
			} else {
				subjectRead->setReadUniqueId(atoi(readname.c_str()));
			}
			this->vpSubjectRead->at(i) = subjectRead;
			if (!(line1.length() < Config::getminimumOverlapLength()
					|| QueryDataset::qualityFilter(line1))) { // Test the read is of good quality.
				subjectRead->bIsGoodRead = false;
			} else {
				subjectRead->bIsGoodRead = true;
			}
		}
	}

	this->iCountSubjectRead = chunkCounter;
	int oldvalue = floor(iReadIdAnchor / 5000000);
	iReadIdAnchor += chunkCounter;
	int newvalue = floor(iReadIdAnchor / 5000000);
	if (newvalue - oldvalue >= 1) {
		CURMEM
		cout << iReadIdAnchor << " good reads have been streamed from the data file(s)." << endl;
		CURTIME
	}
	lineData->clear();
	delete lineData;

	return true;
}
