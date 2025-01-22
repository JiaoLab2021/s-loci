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
#include "../include/filter_bloomfilter.hpp"

using namespace std;
void help_augment(char** argv);

int main_augment(int argc, char** argv) {
    // Bitmap File
	string bitmapFile;

	// Input File
	string inputBaseFile;

	// Output file name prefix
	string prefix = "augment";
	
	// kmer length
	int kmerLen = 21;

	// Number of threads
	int threadsNum = 10;

	// Input Parameters
    int c;
    
    while (true) {
        static const struct option long_options[] = 
		{
			{"fasta", required_argument, 0, 'i'},
			{"bitmapFile", required_argument, 0, 'b'},

            {"kmer", required_argument, 0, 'k'},
			{"prefix", required_argument, 0, 'p'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:b:k:p:t:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c) {
        case 'i':
            inputBaseFile = optarg;
            break;
        case 'b':
            bitmapFile = optarg;
            break;
		case 'p':
            prefix = optarg;
            break;
        case 'k':
            kmerLen = stoi(optarg);
            break;
        case 't':
            threadsNum = stoi(optarg);
            break;
        case 'h':
        case '?':
            help_augment(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_augment(argv);
        return 1;
    }

	// Determine whether the parameters are correct
	if (
		bitmapFile.empty() || 
	    inputBaseFile.empty()
	) {
		cerr << "[" << __func__ << "::" << getTime() << "] " 
			 << "Parameter error: -i -b\n\n";
		help_augment(argv);
        return 1;
	}

	cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
	cerr << "[" << __func__ << "::" << getTime() << "] " << "kmer size: " << kmerLen << endl;

	// Loading a bitmap
	cerr << "[" << __func__ << "::" << getTime() << "] " << "Load the bitmap file.\n";

	// Loading Bloom filter
	BloomFilter<string, HashFun1, HashFun2, HashFun3, 
				HashFun4, HashFun5, HashFun6, 
				HashFun7, HashFun8, HashFun9, 
				HashFun10, HashFun11, HashFun12> bf(bitmapFile);

	// Building the base index
	vector<string> allKmerVectorBase = filter_bloomfilter::build_base_index(
		inputBaseFile, 
		kmerLen,
		threadsNum
	);

	// augment
	// Adding a minimizer to a bloom filter
	for (auto it1 : allKmerVectorBase) {
		bf.set(it1);
	}

	// Clear the vector and release the memory
	allKmerVectorBase.clear();
	vector<string>().swap(allKmerVectorBase);

	// Save the bitmap to memory
	string bitMapFilename = prefix + bitmapFile;
	bf.save(bitMapFilename);

	// Print bitmap usage
	size_t bitMapNum = bf.get_num();
	size_t bitMapCap = bf.get_cap();
	cerr << "[" << __func__ << "::" << getTime() << "] " 
		<< "Bitmap used/capacity: " 
		<< setprecision(2) << bitMapNum/(float)bitMapCap
		<< "(" << to_string(bitMapNum) << "/" << to_string(bitMapCap) << ")" 
		<< endl;

	cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// Help Documentation
void help_augment(char** argv) {
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i -b [options]" << endl
       << "Add haplotype informations to bitmap" << endl
       << endl
       << "required arguments:" << endl
	   << "    -i, --fasta        FILE     input haplotype informations for building bitmap" << endl
	   << "    -b, --bitmapFile   FILE     bitmap file" << endl
	   << endl
	   << "optional arguments:" << endl
       << "    -k, --kmer         INT      k-mer size [21]" << endl
	   << "    -p, --prefix       STRING   output prefix [augment]" << endl
	   << "    -t, --threads      INT      number of compute threads to use [10]" << endl
       << endl
       << "    -h, --help                  print this help document" << endl;
}