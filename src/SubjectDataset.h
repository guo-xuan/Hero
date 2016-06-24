/*
 * SubjectDataset.h
 *
 *  Created on: Oct 16, 2015
 *      Author: xgo
 */

#ifndef SUBJECTDATASET_H_
#define SUBJECTDATASET_H_

#include "Config.h"
#include "SubjectRead.h"
#include "QueryDataset.h"

using namespace std;

enum FileType {FASTA, FASTQ, UNDEFINED};

class SubjectDataset {
public:

	vector<SubjectRead *> * vpSubjectRead;
	UINT64 iCountSubjectRead; // the true number of subject read

	ifstream sInputFileStreamer;
	UINT64 iInputFileIndex;
	FileType eFileType;

	UINT64 iReadIdAnchor;

	UINT64 iLongestReadLength;

	SubjectDataset();
	virtual ~SubjectDataset();
	bool clearSubjectRead();
	bool loadNextChunkParallel(UINT16 _numberOfThreads);

	bool writeToFile(string & _outputFileName);
};

#endif /* SUBJECTDATASET_H_ */
