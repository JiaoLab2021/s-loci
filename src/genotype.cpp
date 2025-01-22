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
#include "../include/hmm.hpp"

using namespace std;
void help_genotype(char** argv);

int main_genotype(int argc, char** argv) {
    // Input File
	string inputBaseFile;
    vector<string> inputReadVec;

    // Output File
    string outputFilePrefix = "out";
	
	// The length of kmer and minimizer
	int kmerLen = 27;

	// Number of threads
	int threadsNum = 10;
    
    // Non-recombination rate
    long double noRecombProb = 0.99;

    // Debugging Code
    bool debug = false;

	// Input Parameters
    int c;
    
    while (true) {
        static const struct option long_options[] = 
		{
			{"fasta", required_argument, 0, 'f'},
            {"read", required_argument, 0, 'r'},

            {"kmer", required_argument, 0, 'k'},
            {"threads", required_argument, 0, 't'},
            {"prefix", required_argument, 0, 'p'},
            {"debug", no_argument, 0, '1'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "f:r:k:t:p:1h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c) {
        case 'f':
            inputBaseFile = optarg;
            break;
        case 'r':
            inputReadVec.push_back(optarg);
            break;
        case 'k':
            kmerLen = stoi(optarg);
            break;
        case 't':
            threadsNum = stoi(optarg);
            break;
        case 'p':
            outputFilePrefix = optarg;
            break;
        case '1':
            debug = true;
            break;
        case 'h':
        case '?':
            help_genotype(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_genotype(argv);
        return 1;
    }

	// Determine whether the parameters are correct
	if (inputBaseFile.empty() || inputReadVec.size() == 0) {
		cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -f -r \n\n";
		help_genotype(argv);
        return 1;
	}

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
	cerr << "[" << __func__ << "::" << getTime() << "] " << "kmer size: " << kmerLen << endl;
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Haplotype information file: " << inputBaseFile << endl;


	// Constructing a haplotype index
    unordered_map<uint64_t, vector<uint64_t> > kmerBaseMap;  // unordered_map<hash, vector<chr/end/strand>>
    unordered_map<uint64_t, tuple<string, uint32_t> > chromosomeIdLenMap;  // map<ID, tuple<chr, length> >
	tie(kmerBaseMap, chromosomeIdLenMap) = genotype::build_base_index(inputBaseFile, kmerLen);

    // Input sequencing file for typing
    // Record the chromosome number of the haplotype corresponding to the kmer matched by the sequencing file, and the number of times END corresponds
    // unordered_map<chrId, map<end, frequency> >
    unordered_map<uint64_t, map<uint32_t, uint32_t> > HapkmerMatMap; // unordered_map<chrId, map<end, frequency> >
    unordered_map<uint64_t, uint64_t> readAllKmerMap;  // map<kmerHash, frequence>
    tie(HapkmerMatMap, readAllKmerMap) = genotype::genotype(
        inputReadVec,  
        kmerLen, 
        kmerBaseMap, 
        chromosomeIdLenMap, 
        threadsNum
    );

    // Forward algorithm to calculate probability
    HMM HMMClass(readAllKmerMap, kmerBaseMap, chromosomeIdLenMap, noRecombProb, debug);
    HMMClass.forward();
    HMMClass.get_path();
    vector<string> outHapVec = HMMClass.get_result(2);

    // Calculate jaccard
    string jaccardTxt = genotype::Jaccard_score(kmerLen, chromosomeIdLenMap, HapkmerMatMap);
    // Calculate the most likely combination and save the file
    genotype::cal_pro(chromosomeIdLenMap, jaccardTxt, outHapVec, outputFilePrefix);

    // Calculate site coverage
    genotype::site_cov(kmerLen, chromosomeIdLenMap, HapkmerMatMap, outputFilePrefix);

	cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// Help Documentation
void help_genotype(char** argv) {
  cerr << "usage: " << argv[0] << " " << argv[1] << " -f FILE -r FILE [options]" << endl
       << "Genotype using haplotypes information" << endl
       << endl
       << "required arguments:" << endl
	   << "    -f, --fasta        FILE     input haplotype informations for genotyping" << endl
       << "    -r, --read         FILE     input fastq/a, possibly compressed, two are allowed, one for each mate" << endl
	   << endl
	   << "optional arguments:" << endl
       << "    -k, --kmer         INT      k-mer size [27]" << endl
	   << "    -t, --threads      INT      number of compute threads to use [10]" << endl
       << "    --prefix           STRING   output prefix [out]" << endl
       << "    --debug                     debug code (false by default)" << endl
       << endl
       << "    -h, --help                  print this help document" << endl;
}