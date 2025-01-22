#include <fstream>
#include <string>
#include <iostream>
#include <getopt.h>
#include "../include/filter_kmer.hpp"

using namespace std;
void help_filter_kmer(char** argv);

int main_filter_kmer(int argc, char** argv) {
    // Input File
	vector<string> inputBaseVec;

	// Files to filter
	vector<string> inputReadVec;

	// Output file name prefix
	string prefix = "sample";
	
	// The length of kmer and minimizer
	unsigned int kmerLen = 27;

	// Whether to calculate the ratio relative to the reference genome
	bool calRefCovBool = false;

	// Read filtering threshold
	float matchTd = 0.2;

	// Number of threads
	int threadsNum = 10;

	// Input Parameters
    int c;
    
    while (true) {
        static const struct option long_options[] = 
		{
			{"input", required_argument, 0, 'i'},
			{"read", required_argument, 0, 'r'},

            {"matchTd", required_argument, 0, 'm'},
            {"minimizerK", required_argument, 0, 'k'},
			{"cal", required_argument, 0, '1'},
			{"prefix", required_argument, 0, '2'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:r:m:k:1:2:t:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c)
        {
        case 'i':
			inputBaseVec.push_back(optarg);
            break;
        case 'r':
			inputReadVec.push_back(optarg);
            break;
        case 'm':
            matchTd = stof(optarg);
            break;
        case 'k':
            kmerLen = stoi(optarg);
            break;
		case '1':
            calRefCovBool = true;
            break;
        case '2':
            prefix = optarg;
            break;
        case 't':
            threadsNum = stoi(optarg);
            break;
        case 'h':
        case '?':
            help_filter_kmer(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_filter_kmer(argv);
        return 1;
    }

	// Determine whether the parameters are correct
	if (
		inputBaseVec.size() == 0 || 
		inputReadVec.size() == 0 || 
		matchTd <= 0 || 
		matchTd > 1
	) {
		cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -i (-f / -F)\n\n";
		help_filter_kmer(argv);
        return 1;
	}

	cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";
	cerr << "[" << __func__ << "::" << getTime() << "] " << "kmer size: " << kmerLen << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] " << "Filter threshold: " << matchTd << endl;
	cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculate coverage relative to refgenome: " << calRefCovBool << endl;
	
	// Build base file index
	// Record read kmer map<kmerScore, map<readId, kmerNum>>
	unordered_map<uint64_t, unordered_map<uint32_t, uint32_t>> kmerIdMap;
	// Record the number of kmers of each ref to determine coverage <readId, kmerNum (number)>
	unordered_map<uint32_t, uint32_t> refkmerNumMap;
	
	tie(kmerIdMap, refkmerNumMap) = filter_kmer::build_base_index(
		inputBaseVec, 
		kmerLen, 
		threadsNum
	);

	malloc_trim(0);	// 0 is for heap memory

	// Filtering Sequencing Files
	filter_kmer::build_read_index(
		inputReadVec, 
		kmerLen,
		matchTd, 
		kmerIdMap, 
		refkmerNumMap,
		threadsNum,
		prefix,
		calRefCovBool
	);

	cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// Help Documentation
void help_filter_kmer(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i FILE -r FILE [options]" << endl
       << "Filter fastq/a files using haplotype informations" << endl
       << endl
       << "required arguments:" << endl
	   << "    -i, --input       FILE     fasta/q file for index building, two are allowed, one for each mate" << endl
	   << "    -r, --read        FILE     input fastq/a, possibly compressed, two are allowed, one for each mate" << endl
	   << endl
	   << "optional arguments:" << endl
	   << "    -m, --matchTd      FLOAT    read filtering threshold (0-1] [0.2]" << endl
       << "    -k, --kmer         INT      k-mer size (no larger than 28) [27]" << endl
	   << "    --cal                       calculate coverage relative to refgenome (default: false)" << endl
	   << "    --prefix           STRING   output prefix [sample]" << endl
	   << "    -t, --threads      INT      number of compute threads to use [10]" << endl
       << endl
       << "    -h, --help                  print this help document" << endl;
}