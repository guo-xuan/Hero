/*
 * Config.cpp
 *
 *  Created on: Oct 6, 2015
 *      Author: xgo
 */

#include "Config.h"

//Here is to create and initialize the static members
string Config::sQueryFileName = "";
vector<string> Config::vSubjectFiles;
vector<string> Config::vQueryFiles;
string Config::outputfilename = "";
string Config::sWorkingDirectory = "./";

INT64 Config::iKmer = 39ll; //key needs to be smaller than minimumoverlap, if double key, 2*key<minoverlap
UINT64 Config::streamChunkSize = 400;
UINT64 Config::queryChunkSize = 800;
bool Config::useErrorRate = false;
UINT16 Config::maxMismatch = 1; //only valid when perfectMatch is set to false
UINT16 Config::maxInsertion = 0; //only valid when perfectMatch is set to false
UINT16 Config::maxDeletion = 0;
UINT16 Config::numberOfThreads = 1;
double Config::dMaxMismatchRate = 1; //only valid when useErrorRate is set to true, in percentage
double Config::maxInsertionRate = 1; //only valid when useErrorRate is set to true, in percentage
double Config::maxDeletionRate = 1; //only valid when useErrorRate is set to true, in percentage
UINT16 Config::maxIndel = 1;
double Config::dMaxIndelRate = 1;
bool Config::isNumberOfThreadsSet = false;
bool Config::storeOverhang = false;
INT64 Config::minimumOverlapLength = 40ll;
bool Config::bDebug = false;
INT64 Config::iLongestRead = 700;
UINT16 Config::iAlignerType = BandAligner::IGLOBAL;
bool Config::bSpeicialVersion = false;
UINT64 Config::iHashTableFactor = 32;

string Config::sSuffix = "";
UINT64 Config::iSizeLimit = 0;

bool Config::bWriteFinishFlag = false;
string Config::sFlagFolder = "";

bool Config::bCallConsensus = false;
bool Config::bMergedReads = false;
double Config::dMinFrequencyCallConsensus = 0.6;
bool Config::bCoverageStatistic = false;
UINT32 Config::iCoverageDepthFilter = 1;

bool Config::bMerging = false;
int Config::iMinMergeOverlap = 10;

bool Config::bRealMismatchRate = true;

bool Config::bMappingVersion = false;

Config::Config() {
	// TODO Auto-generated constructor stub

}

Config::~Config() {
	// TODO Auto-generated destructor stub
}

UINT16 Config::getminimumOverlapLength() {
	return Config::minimumOverlapLength;
}

string Config::getQueryDatasetFilename() {
	return Config::sQueryFileName;
}

void Config::printHelp() {

	cout << endl << "Usage:" << endl << "PairwiseAliger [OPTION]... <PARAM>..." << endl
			<< "[OPTION]                                                    " << endl
			<< "  -h/--help\t only print out the help contents" << endl
			<< "  -s/--subject\t list of subject file name(s) (comma separated)" << endl
			<< "  -q/--query\t query file name (usually it's a small piece of large subject file(s))" << endl
			<< "  -o/--out\t output file name (default:query.aln)" << endl
			<< "  -c\t call the consensus when frequency larger than the given value (default: 0.6)" << endl
			<< "  --covfilter\t the minimum coverage for a read (default: 1)" << endl
			<< "  --covinfo\t storing the coverage info (default: false)" << endl
			<< "  --merge\t enable merging feature (default: false)           " << endl
			<< "  --singleReads\t indicate the input are unpaired reads       " << endl << endl

			<< "<PARAM: aligning setting>                                      " << endl
			<< "  -k\t k-mer length (default:39)                             " << endl
			<< "  -m\t maximum allowed mismatch rate in percentage [1 means 1% of the overlap length] " << endl
			<< "  -i\t maximum allowed INDEL rate in percentage [1 means 1% of the overlap length] " << endl
			<< "  -l\t minimum overlap length (default:40)                     " << endl << endl

			<< "<PARAM: multithreading setting>                                " << endl
			<< "  -t\t number of threads (default:1 [single thread])" << endl
			<< "  -y\t stream chunk size of query read file (default:400)" << endl
			<< "  -z\t stream chunk size of subject read file (default:400) " << endl << endl

			<< "<PARAM: setting for long reads ( > 1000 bp )>                   " << endl
			<< "  -a\t choose aligner type [a: affine gap alignment, g: global alignment(default)]" << endl
			<< "  -v\t special version (show the alignment in a human readable way)" << endl
			<< "  -f\t hash table factor (32 default, for high coverage dataset, set it to a large integer)" << endl
			<< "  -mapping\t enable mapping (only support substitution) to generate base pair coverage info." << endl
			<< endl

			//<< "<PARAM: shard and split setting>                                " << endl
			//<< "  --suffix\t suffix for the output file                         " << endl
			//<< "  --size\t the number of reads per shard                        " << endl
			//<< "  --notification\t the folder store the notification when the alignment completed" << endl << endl

			<< "Example: " << endl
			<< "  ./PairwiseAliger -q qreads.fasta -s sreads1.fasta,sreads2.fasta -out outreads -m 5 -i 5 -k 40 -t 4 -y 4000 -z 4000 -l 30 -v -a g -f 40"
			<< endl << endl;

}

bool Config::setConfig(int argc, char **argv) {

	Config::vSubjectFiles.clear();

	vector<string> argumentsList;
	cout << "PRINTING ARGUMENTS" << endl;
	for (int i = 0; i < argc; i++) {
		cout << argv[i] << ' ';
	}
	cout << endl;
	while (argc--)
		argumentsList.push_back(*argv++);

	if (argumentsList.size() == 1) {
		Config::printHelp();
		return false;
	}

	for (UINT64 i = 1; i <= argumentsList.size() - 1; i++) {
		if (argumentsList[i] == "-h" || argumentsList[i] == "--help") {
			Config::printHelp();
			return false;
		} else if (argumentsList[i] == "-s" || argumentsList[i] == "--subject") {
			string inputFilenames = argumentsList[++i];
			stringstream ss(inputFilenames);
			string item;
			while (getline(ss, item, ',')) {
				Config::vSubjectFiles.push_back(item);
			}
		} else if (argumentsList[i] == "-q" || argumentsList[i] == "--query") {
			string inputFilename = argumentsList[++i];
			Config::sQueryFileName = inputFilename;
			stringstream ss(inputFilename);
			string item;
			while (getline(ss, item, ',')) {
				Config::vQueryFiles.push_back(item);
			}
			if (Config::vQueryFiles.size() > 1) {
				Config::sQueryFileName = "";
			} else {
				Config::vQueryFiles.clear();
			}
		} else if (argumentsList[i] == "-o" || argumentsList[i] == "--out") {
			string outputFilename = argumentsList[++i];
			Config::outputfilename = outputFilename;
		} else if (argumentsList[i] == "-k")
			Config::iKmer = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-m") {
			Config::useErrorRate = true;
			Config::dMaxMismatchRate = atof(argumentsList[++i].c_str()) / 100.0;
		} else if (argumentsList[i] == "-i") {
			Config::useErrorRate = true;
			Config::dMaxIndelRate = atof(argumentsList[++i].c_str()) / 100.0;
		} else if (argumentsList[i] == "-t") {
			Config::numberOfThreads = atoi(argumentsList[++i].c_str());
			Config::isNumberOfThreadsSet = true;
		} else if (argumentsList[i] == "-z") {
			Config::streamChunkSize = atoi(argumentsList[++i].c_str());
		} else if (argumentsList[i] == "-y") {
			Config::queryChunkSize = atoi(argumentsList[++i].c_str());
		} else if (argumentsList[i] == "-l") {
			Config::minimumOverlapLength = atoi(argumentsList[++i].c_str());
		} else if (argumentsList[i] == "-a") {
			string sTemp = argumentsList[++i];
			if (sTemp == "a") {
				Config::iAlignerType = BandAligner::IAFFINE;
			} else if (sTemp == "g") {
				Config::iAlignerType = BandAligner::IGLOBAL;
			}
		} else if (argumentsList[i] == "-v") {
			Config::bSpeicialVersion = true;
		} else if (argumentsList[i] == "--mapping") {
			Config::bMappingVersion = true;
		} else if (argumentsList[i] == "-f") {
			Config::iHashTableFactor = atoi(argumentsList[++i].c_str());
		} else if (argumentsList[i] == "--suffix") {
			string sTemp = argumentsList[++i];
			Config::sSuffix = sTemp;
		} else if (argumentsList[i] == "--size") {
			Config::iSizeLimit = atoi(argumentsList[++i].c_str());
		} else if (argumentsList[i] == "-c") {
			Config::bCallConsensus = true;
			if (i + 1 < argumentsList.size() && argumentsList[i + 1].at(0) >= '0'
					&& argumentsList[i + 1].at(0) <= '9') {
				Config::dMinFrequencyCallConsensus = atof(argumentsList[++i].c_str());
			}
		} else if (argumentsList[i] == "--covfilter") {
			Config::iCoverageDepthFilter = atoi(argumentsList[++i].c_str());
		} else if (argumentsList[i] == "--singleReads") {
			Config::bMergedReads = true;
		} else if (argumentsList[i] == "--notification") {
			Config::bWriteFinishFlag = true;
			Config::sFlagFolder = argumentsList[++i];
			if (Config::sFlagFolder.at(Config::sFlagFolder.length() - 1) != '/') {
				Config::sFlagFolder.append("/");
			}
		} else if (argumentsList[i] == "--covinfo") {
			Config::bCoverageStatistic = true;
		} else if (argumentsList[i] == "--merge") {
			Config::bMerging = true;
		} else {
			Config::printHelp();
			return false;
		}
	}

	if ((Config::outputfilename == "")) {
		Config::outputfilename = "query";
	} else {
		int len = Config::outputfilename.length();
		if (len >= 4) {
			int pos = Config::outputfilename.find_last_of(".", 0);
			if (Config::outputfilename.substr(pos + 1) == "aln") {
				Config::outputfilename = Config::outputfilename.substr(0, pos);
			}
		}
	}

	if ((Config::vSubjectFiles.size() == 0) || (Config::sQueryFileName == "" && Config::vQueryFiles.size() == 0)) {
		cout << "missed -subject or -query input files!" << endl;
		printHelp();
		return false;
	}

	return true;
}
