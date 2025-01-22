#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <regex>
#include <getopt.h>
#include "../include/assembly.hpp"

using namespace std;
void help_assembly(char** argv);

int main_assembly(int argc, char** argv) {
    // Input File
	string inputFastqFile1;
	string inputFastqFile2;
	string inputBaseFile;

    // Sequencing Type
    int sequenceType = 0;
    uint32_t insLen;

	// Output file name prefix
	string prefix = "out";

	// kmer length
	int kmerLen = 31;

	// Number of threads
	int threadsNum = 10;

    // Software Path
    string SOAPdenovoPath = "";
    string minimap2Path = "";

	// Input Parameters
    int c;

    while (true) {
        static const struct option long_options[] = 
		{
			{"fastq", required_argument, 0, 'i'},
			{"avg_ins", required_argument, 0, '1'},

            {"type", required_argument, 0, '2'},
            {"kmer", required_argument, 0, 'k'},
			{"prefix", required_argument, 0, 'p'},
			{"fasta", required_argument, 0, 'f'},
            {"threads", required_argument, 0, 't'},

            {"SOAPdenovo-63mer", required_argument, 0, '3'},
            {"minimap2", required_argument, 0, '4'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "i:1:2:k:p:f:t:3:4:h", long_options, &option_index);

        if (c == -1)
            break;
        
        switch (c) {
        case 'i':
			if (inputFastqFile1.empty()) {
				inputFastqFile1 = optarg;
			} else {
				inputFastqFile2 = optarg;
			}
            break;
        case '1':
            insLen = stoi(optarg);
            break;
        case '2':
            sequenceType = stoi(optarg);
            break;
        case 'k':
            kmerLen = stoi(optarg);
            break;
        case 'p':
            prefix = optarg;
            break;
        case 'f':
            inputBaseFile = optarg;
            break;
        case 't':
            threadsNum = stoi(optarg);
            break;

        case '3':
            SOAPdenovoPath = optarg;
            break;
        case '4':
            minimap2Path = optarg;
            break;
        case 'h':
        case '?':
            help_assembly(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_assembly(argv);
        return 1;
    }

	// Determine whether the parameters are correct
	if (inputFastqFile1.empty() && inputFastqFile2.empty()) {
		cerr << "[" << __func__ << "::" << getTime() << "] " 
			 << "Parameter error: -i\n\n";
		help_assembly(argv);
        return 1;
	}

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running.\n";

    // Record the shortest read, longest read, average read length, and total read length
    uint32_t minLen;
    uint32_t maxLen;
    float aveLen;
    uint64_t readBase;
    uint64_t readNum;
    
    // Calculate read length of sequencing data
    tie(minLen, maxLen, aveLen, readBase, readNum) = assembly::fastq_count(inputFastqFile1, inputFastqFile2);

    // Generate configuration files
    assembly::make_config(aveLen, insLen, inputFastqFile1, inputFastqFile2, sequenceType);

    // assembly
    assembly::SOAPdenovo(SOAPdenovoPath, kmerLen, prefix, threadsNum);

    // alignment
    if (inputBaseFile.size() > 0) {
        assembly::minimap2(minimap2Path, inputBaseFile, prefix, threadsNum);
    }
    
	cerr << "[" << __func__ << "::" << getTime() << "] " << "Done.\n";

    return 0;
}

// Help Documentation
void help_assembly(char** argv) {
  cerr << "usage: " << argv[0] << " " << argv[1] << " -i --avg_ins [options]" << endl
       << "Assemble haplotype sequences based on filtered reads" << endl
       << endl
       << "required arguments:" << endl
	   << "    -i, --fastq         FILE     input fastq, possibly compressed, two are allowed, one for each mate" << endl
	   << "    --avg_ins           INT      average insert length of library" << endl
	   << endl
	   << "optional arguments:" << endl
       << "    --type              INT      0-pair-end, 1-mate-pair [0]" << endl
       << "    -k, --kmer          INT      k-mer size for assembly [31]" << endl
	   << "    -p, --prefix        STRING   output prefix [out]" << endl
	   << "    -f, --fasta         FILE     gene sequences for alignment (minimap2) [None]" << endl
	   << "    -t, --threads       INT      number of compute threads to use [10]" << endl
       << endl
       << "software path:" << endl
       << "    --SOAPdenovo-63mer  PATH     SOAPdenovo-63mer path [None]" << endl
       << "    --minimap2          PATH     minimap2 path [None]" << endl
       << endl
       << "    -h, --help                  print this help document" << endl;
}