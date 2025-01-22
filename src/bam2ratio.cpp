#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <thread>
#include "zlib.h"
#include <regex>
#include <getopt.h>
#include "../include/bam2ratio.hpp"

using namespace std;
void help_bam2ratio(char** argv);

int main_bam2ratio(int argc, char** argv) {
    // Input File
	string bamFileName;

    // Read filtering threshold
	float matchTd = 0.2;

    // Number of read splits
    int readSplitNum = 1000;

    // Whether to calculate the ratio relative to the reference genome
	string calRefCovBool = "false";

    // Output File
    string outputFileName;

	// Number of threads
	int threadsNum = 10;

	// Input Parameters
    int c;
    
    while (true) {
        static const struct option long_options[] = 
		{
			{"bam", required_argument, 0, 'b'},

            {"matchTd", required_argument, 0, 'm'},
            {"number", required_argument, 0, 'n'},
            {"cal", required_argument, 0, 'c'},
            {"output", required_argument, 0, 'o'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "b:m:n:c:o:t:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c) {
        case 'b':
            bamFileName = optarg;
            break;
        case 'm':
            matchTd = stof(optarg);
            break;
        case 'n':
            readSplitNum = stoi(optarg);
            break;
		case 'c':
            calRefCovBool = optarg;
            break;
        case 'o':
            outputFileName = optarg;
            break;
        case 't':
            threadsNum = stoi(optarg);
            break;
        case 'h':
        case '?':
            help_bam2ratio(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_bam2ratio(argv);
        return 1;
    }

	// Determine whether the parameters are correct
	if (bamFileName.empty())
	{
		cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -b\n\n";
		help_bam2ratio(argv);
        return 1;
	}

	cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";

    bam2ratio::bam2ratio(
        (char*)bamFileName.c_str(),
        matchTd, 
        readSplitNum, 
        outputFileName, 
        threadsNum, 
        calRefCovBool
    );

	cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// Help Documentation
void help_bam2ratio(char** argv) {
  cerr << "usage: " << argv[0] << " " << argv[1] << " -b [options]" << endl
       << "Build minimizer bam2ratio for file" << endl
       << endl
       << "required arguments:" << endl
	   << "    -b, --bam          FILE     input BAM" << endl
	   << endl
	   << "optional arguments:" << endl
       << "    -m, --matchTd      FLOAT    read filtering threshold (0-1) [0.2]" << endl
       << "    -n, --number       INT      number of reads used by each thread [1000]" << endl
       << "    -c, --cal          BOOL     calculate coverage relative to refgenome (true/false) [false]" << endl
       << "    -o, --output       FILE     output filename [stdout]" << endl
	   << "    -t, --threads      INT      number of compute threads to use [10]" << endl
       << endl
       << "    -h, --help                  print this help document" << endl;
}