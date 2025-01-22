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
#include "../include/genotype.hpp"

using namespace std;
void help_index(char** argv);

int main_index(int argc, char** argv) {
    // Input File
	string inputBaseFile;

    // Output File
    string outputFileName;
	
	// The length of kmer and minimizer
	int kmerLen = 21;

	// Number of threads
	int threadsNum = 10;

	// Input Parameters
    int c;
    
    while (true)
    {
        static const struct option long_options[] = 
		{
			{"fasta", required_argument, 0, 'i'},

            {"kmer", required_argument, 0, 'k'},
            {"output", required_argument, 0, 'o'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:k:o:t:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 'i':
            inputBaseFile = optarg;
            break;
        case 'k':
            kmerLen = stoi(optarg);
            break;
        case 'o':
            outputFileName = optarg;
            break;
        case 't':
            threadsNum = stoi(optarg);
            break;
        case 'h':
        case '?':
            help_index(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_index(argv);
        return 1;
    }

	// Determine whether the parameters are correct
	if (inputBaseFile.empty()) {
		cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -i\n\n";
		help_index(argv);
        return 1;
	}

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
	cerr << "[" << __func__ << "::" << getTime() << "] " << "kmer size: " << kmerLen << endl;

	// Constructing a haplotype index
    unordered_map<uint64_t, vector<uint64_t>> kmerBaseMap;
    unordered_map<uint64_t, tuple<string, uint32_t>> chromosomeIdLenMap;
	tie(kmerBaseMap, chromosomeIdLenMap) = genotype::build_base_index(inputBaseFile, kmerLen);
	
    // Output
    if (outputFileName.empty()) {  // If no output file name is specified, it is printed to standard output.
        for (auto it1 : kmerBaseMap) {
            for (auto it2 : it1.second) {
                cout << it1.first << "\t" << it2 << endl;
            }
        }
    } else {  // Save to file
        string OutTxt;

        for (auto it1 : kmerBaseMap) {
            for (auto it2 : it1.second) {
                OutTxt += to_string(it1.first) + "\t" + to_string(it2) + "\n";
            }
        }

        // Save to file
		gzFile gzfp = gzopen(outputFileName.c_str(), "wb");
		gzwrite(gzfp, OutTxt.c_str(), OutTxt.length());

		// Freeing up memory
		gzclose(gzfp);
    }

	cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// Help Documentation
void help_index(char** argv) {
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i [options]" << endl
       << "Build minimizer index for file" << endl
       << endl
       << "required arguments:" << endl
	   << "    -i, --fasta        FILE     input file for indexing" << endl
	   << endl
	   << "optional arguments:" << endl
       << "    -k, --kmer         INT      k-mer size [21]" << endl
       << "    -o, --output       FILE     output filename [stdout]" << endl
	   << "    -t, --threads      INT      number of compute threads to use [10]" << endl
       << endl
       << "    -h, --help                  print this help document" << endl;
}