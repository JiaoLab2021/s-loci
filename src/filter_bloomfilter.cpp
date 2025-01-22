#include <fstream>
#include <string>
#include <iostream>
#include <getopt.h>
#include "../include/filter_bloomfilter.hpp"

using namespace std;
void help_filter_bloomfilter(char** argv);

int main_filter_bloomfilter(int argc, char** argv) {
    // Bitmap size or file
	long long int bitmapSize = 2857660380;
	string bitmapFile;

	// Input File
	string inputBaseFile;
	string inputFastqFile1;
	string inputFastqFile2;
	string inputFastaFile1;
	string inputFastaFile2;

	// Output file name prefix
	string prefix = "sample";
	
	// kmer length
	int kmerLen = 21;

	// Read filtering threshold
	float matchTd = 0.85;

	// Number of threads
	int threadsNum = 10;

	// Input Parameters
    int c;
    
    while (true) {
        static const struct option long_options[] = 
		{
			{"fasta", required_argument, 0, 'i'},
			{"bitmapSize", required_argument, 0, 's'},
            {"bitmapFile", required_argument, 0, 'b'},
			{"FASTQ", required_argument, 0, 'f'},
			{"FASTA", required_argument, 0, 'F'},

            {"matchTd", required_argument, 0, 'm'},
            {"kmer", required_argument, 0, 'k'},
			{"prefix", required_argument, 0, 'p'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:s:b:f:F:m:k:p:t:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 'i':
            inputBaseFile = optarg;
            break;
        case 's':
            bitmapSize = stol(optarg);
            break;
        case 'b':
            bitmapFile = optarg;
            break;
        case 'f':
			if (inputFastqFile1.empty()) {
				inputFastqFile1 = optarg;
			} else {
				inputFastqFile2 = optarg;
			}
            break;
        case 'F':
			if (inputFastaFile1.empty()) {
				inputFastaFile1 = optarg;
			} else {
				inputFastaFile2 = optarg;
			}
            break;
        case 'm':
            matchTd = stof(optarg);
            break;
        case 'k':
            kmerLen = stoi(optarg);
            break;
        case 'p':
            prefix = optarg;
            break;
        case 't':
            threadsNum = stoi(optarg);
            break;
        case 'h':
        case '?':
            help_filter_bloomfilter(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_filter_bloomfilter(argv);
        return 1;
    }

	// Determine whether the parameters are correct
	if ((bitmapFile.empty() && (inputBaseFile.empty() || bitmapSize == 0)) || 
		(inputFastqFile1.empty() && inputFastqFile2.empty() && inputFastaFile1.empty() && inputFastaFile2.empty()) || 
		matchTd <= 0 || 
		matchTd > 1) {
		cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: (-f / -F) (-i -s / -b)\n\n";
		help_filter_bloomfilter(argv);
        return 1;
	}

	cerr << "[" << __func__ << "::" << getTime() << "] "<< "Running.\n";
	cerr << "[" << __func__ << "::" << getTime() << "] "<< "kmer size: " << kmerLen << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] "<< "Filter threshold: " << matchTd << endl;

	// Determine whether to reconstruct the Bloom filter or load it from memory
	if (inputBaseFile.length() > 0 && bitmapFile.length() == 0) {
		BloomFilter<string, HashFun1, HashFun2, HashFun3, 
					HashFun4, HashFun5, HashFun6, 
					HashFun7, HashFun8, HashFun9, 
					HashFun10, HashFun11, HashFun12> bf(bitmapSize);
	
		// Build base index and Bloom filter
		vector<string> allKmerVectorBase = filter_bloomfilter::build_base_index(
			inputBaseFile, 
			kmerLen,
			threadsNum
		);
		cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Kmer number: " 
             << allKmerVectorBase.size() 
             << endl;

        // Adding a minimizer to a bloom filter
		for (auto it1 : allKmerVectorBase) {
			bf.set(it1);
		}

		// Clear the vector and release the memory
		allKmerVectorBase.clear();
		vector<string>().swap(allKmerVectorBase);

        // Save the bitmap to memory
		string bitMapFilename = "BitMap._s" + to_string(bitmapSize) + "_k" + to_string(kmerLen) + ".dat";
		bf.save(bitMapFilename);

		// Print bitmap usage
		size_t bitMapNum = bf.get_num();
		size_t bitMapCap = bf.get_cap();
		cerr << "[" << __func__ << "::" << getTime() << "] "
			<< "Bitmap used/capacity: " 
			<< setprecision(2) << bitMapNum/(float)bitMapCap
			<< "(" << to_string(bitMapNum) << "/" << to_string(bitMapCap) << ")" 
			<< endl;

		// Filter fastq files
		cerr << "[" << __func__ << "::" << getTime() << "] "<< "Filtering.\n";
		filter_bloomfilter::build_sequence_index(
			inputFastqFile1, 
			inputFastqFile2, 
			inputFastaFile1, 
			inputFastaFile2, 
			kmerLen, 
			matchTd, 
			bf,
			threadsNum,
			prefix
		);
	} else if (inputBaseFile.length() == 0 && bitmapFile.length() > 0) {  // Load from memory
		cerr << "[" << __func__ << "::" << getTime() << "] " << "Load the bitmap file.\n";
		// Loading Bloom filter
		BloomFilter<string, HashFun1, HashFun2, HashFun3, 
					HashFun4, HashFun5, HashFun6, 
					HashFun7, HashFun8, HashFun9, 
					HashFun10, HashFun11, HashFun12> bf(bitmapFile);

        // Print bitmap usage
		size_t bitMapNum = bf.get_num();
		size_t bitMapCap = bf.get_cap();
		cerr << "[" << __func__ << "::" << getTime() << "] "
			<< "Bitmap used/capacity: " 
			<< setprecision(2) << bitMapNum/(float)bitMapCap
			<< "(" << to_string(bitMapNum) << "/" << to_string(bitMapCap) << ")" 
			<< endl;

		// Filter fastq files
		cerr << "[" << __func__ << "::" << getTime() << "] "<< "Filtering.\n";
		filter_bloomfilter::build_sequence_index(
			inputFastqFile1, 
			inputFastqFile2, 
			inputFastaFile1, 
			inputFastaFile2, 
			kmerLen, 
			matchTd, 
			bf,
			threadsNum,
			prefix
		);
	} else {  // Unclear whether to load the bitmap or reconstruct it
		cerr << "[" << __func__ << "::" << getTime() << "] "
			 << "Not clear whether to load or reconstruct the bitmap \n";
		exit(1);
	}

	cerr << "[" << __func__ << "::" << getTime() << "] "<< "Done.\n";

    return 0;
}

// Help Documentation
void help_filter_bloomfilter(char** argv) {
  cerr << "usage: " << argv[0] << " " << argv[1] << " (-f / -F) (-i -s / -b) [options]" << endl
       << "Filter fastq/a files using haplotype informations" << endl
       << endl
       << "required arguments:" << endl
	   << "    -i, --fasta        FILE     input haplotype informations for building bitmap" << endl
	   << "    -s, --bitmapSize   INT      size of the bitmap [2857660380]" << endl
	   << "    -b, --bitmapFile   FILE     loads bitmap file from memory (-f -s or -b)" << endl
	   << "    -f, --FASTQ        FILE     input fastq, possibly compressed, two are allowed, one for each mate" << endl
	   << "    -F, --FASTA        FILE     input fasta, possibly compressed, two are allowed, one for each mate" << endl
	   << endl
	   << "optional arguments:" << endl
	   << "    -m, --matchTd      FLOAT    read filtering threshold (0-1) [0.85]" << endl
       << "    -k, --kmer         INT      k-mer size [21]" << endl
	   << "    -p, --prefix       STRING   output prefix [sample]" << endl
	   << "    -t, --threads      INT      number of compute threads to use [10]" << endl
       << endl
       << "    -h, --help                  print this help document" << endl;
}