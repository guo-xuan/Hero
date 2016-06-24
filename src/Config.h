/*
 * Config.h
 *
 *  Created on: Oct 6, 2015
 *      Author: xgo
 */

#ifndef CONFIG_H_
#define CONFIG_H_

// Common headers:

//multi-thread library OPENMP
#include <omp.h>

// C headers:
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

// C++ headers:
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <utility>
#include <limits>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <map>
#include <exception>

using namespace std;

//define CLOCKSTART clock_t begin = clock(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
//define CLOCKSTOP clock_t end = clock(); cout << "Function " << __FUNCTION__ << "() finished in " << double(end - begin) / CLOCKS_PER_SEC<< " Seconds." << endl << endl;
#define CLOCKSTART double begin = omp_get_wtime(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
#define CLOCKSTOP double end = omp_get_wtime(); cout << "Function " << __FUNCTION__ << "() finished in " << double(end - begin) << " Seconds." << endl << endl;

#define CURTIME time_t curr=time(0);cout << "current time is: " << ctime(&curr) <<endl;
#define CURMEM INT64 mem = checkMemoryUsage();cout<<"Memory used: " << mem <<  " MB."<< endl;

// To keep time information of functions.
#define MEMORYSTART INT64 mem_start = checkMemoryUsage(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
#define MEMORYSTOP INT64 mem_end = checkMemoryUsage(); cout << "Function " << __FUNCTION__ << "() finished. " << "Memory used: " << mem_end << " - " <<  mem_start << " = "<< mem_end - mem_start << " MB."<< endl;
// Get the memory usage with a Linux kernel.
inline unsigned int checkMemoryUsage() {
// get KB memory into count
	unsigned int count = 0;

#if defined(__linux__)
	ifstream f("/proc/self/status"); // read the linux file
	while (!f.eof()) {
		string key;
		f >> key;
		if (key == "VmData:") {     // size of data
			f >> count;
			break;
		}

	}
	f.close();
#endif

// return MBs memory (size of data)
	return (count / 1024);
}
;

#define MAX_VECTOR 500000
typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef unsigned int UINT32;
typedef int INT32;
typedef unsigned long long UINT64;
typedef long long INT64;

#include "BandAligner.h"

class Config {

public:
	static string sQueryFileName;
	static vector<string> vSubjectFiles;
	static vector<string> vQueryFiles;
	static string outputfilename;
	static string sWorkingDirectory;

	static INT64 iKmer;
	static UINT64 streamChunkSize;
	static UINT64 queryChunkSize; //used for parallel I/O
	static UINT16 maxMismatch;
	static UINT16 maxInsertion;
	static UINT16 maxDeletion;
	static UINT16 maxIndel;
	static bool useErrorRate;
	static double dMaxMismatchRate;
	static double maxInsertionRate;
	static double maxDeletionRate;
	static double dMaxIndelRate;
	static bool isNumberOfThreadsSet;
	static UINT16 numberOfThreads;
	static bool storeOverhang;
	static INT64 minimumOverlapLength; //it's left side minimumoverlaplength when used in the paired end merging process
	static bool bDebug;
	static INT64 iLongestRead; // the shorter read among the longest reads from query or subject
	static UINT16 iAlignerType; // 0 affine gap; 1 global;
	static bool bSpeicialVersion; // output alignment in a human readable way
	static UINT64 iHashTableFactor; // hash table size
	static string sSuffix; // the suffix for the output alignment file
	static UINT64 iSizeLimit; // the maximum number reads per alignment file
	static bool bWriteFinishFlag; // once finish writing an alignment file, write out a empty file with that name
	static string sFlagFolder;

	static bool bCallConsensus;
	static bool bMergedReads;
	static double dMinFrequencyCallConsensus;
	static bool bCoverageStatistic;
	static UINT32 iCoverageDepthFilter;

	static bool bMerging;
	static int iMinMergeOverlap;

	static bool bRealMismatchRate;

	static bool bMappingVersion;

	Config();
	virtual ~Config();

	static inline unsigned int checkMemoryUsage();
	static UINT16 getminimumOverlapLength();
	static string getQueryDatasetFilename();	// query is database A loaded to memory
	static vector<string> getSubjectDatasetFilenames(); // subject is dataset B streamed from hard drive
	static void printHelp();
	static bool setConfig(int argc, char **argv);
};

#endif /* CONFIG_H_ */
