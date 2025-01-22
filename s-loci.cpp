// g++ s-loci.cpp -o s-loci -lpthread -lz -lhts -O2 -std=c++17
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <mutex>
#include <iomanip>
#include "src/filter_bloomfilter.cpp"
#include "src/filter_kmer.cpp"
#include "src/bam2ratio.cpp"
#include "src/genotype.cpp"
#include "src/assembly.cpp"
#include "src/index.cpp"
#include "src/augment.cpp"

using namespace std;


// define data
#define PROGRAM_DATA "2025/01/22"
// define version
#define PROGRAM_VERSION "1.0.1"
// define author
#define PROGRAM_AUTHOR "Zezhen Du"
// define E-mail
#define PROGRAM_E_MAIL "dzz0539@gmail.com or dzz0539@163.com"


void help(char** argv);

int main(int argc, char** argv) {
	// Print Help Document
	if (argc == 1) {
        help(argv);
        return 1;
    }

	// Selected sub-function
	string subcommand = argv[1];
	
	if (subcommand == "-h" || subcommand == "--help") {
		help(argv);
        return 1;
	} else if (subcommand == "bloomfilter") {
		main_filter_bloomfilter(argc, argv);	
	} else if (subcommand == "kmerfilter") {
		main_filter_kmer(argc, argv);	
	} else if (subcommand == "bam2ratio") {
		main_bam2ratio(argc, argv);	
	} else if (subcommand == "genotype") {
		main_genotype(argc, argv);	
	} else if (subcommand == "assembly") {
		main_assembly(argc, argv);
	} else if (subcommand == "augment") {
		main_augment(argc, argv);
	} else if (subcommand == "index") {
		main_index(argc, argv);	
	} else {
		cerr << "Error: ["<< argv[0] << "] command " << subcommand << " not found" << endl;
		help(argv);
        return 1;
	}

    return 0;
}

// Help Documentation
void help(char** argv)
{
  cerr << "usage: " << argv[0] << " <command> [options]" << endl
	   << endl
	   << "Data: " << PROGRAM_DATA << endl
	   << "Version: " << PROGRAM_VERSION << endl
	   << "Author: " << PROGRAM_AUTHOR << endl
	   << endl
       << "subcommands:" << endl
       << "   bloomfilter      filter fastq/a files based on haplotypes information (fasta)" << endl
	   << "   kmerfilter       filter fastq/a files based on haplotypes information (fasta)" << endl
	   << "   bam2ratio        filter fastq/a files based on alignment score (bam)" << endl
	   << "   genotype         genotyping based on haplotypes information (fasta)" << endl
	   << "   assembly         assemble haplotype sequences based on filtered reads" << endl
	   << "   augment          add haplotype information to bitmap" << endl
	   << "   index            build minimizer index (fasta)" << endl
       << endl
       << "   -h, --help       print this help document" << endl;
}