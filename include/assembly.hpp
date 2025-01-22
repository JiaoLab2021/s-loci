#ifndef assembly_hpp
#define assembly_hpp
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include <regex>
#include <tuple>
#include "kmer.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "shell.hpp"

using namespace std;

namespace assembly{
    // Open fastq.gz file
    // tuple<minLen, maxLen, aveLen, readBase, readNum>
    tuple<uint32_t, uint32_t, float, uint64_t, uint64_t> fastq_count(string inputFileName1, string inputFileName2) {
        // log
        cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Calculate the base number and length of fastq." 
                    << endl;

        // Record the shortest read, longest read, average read length, and total read length
        uint32_t minLen = INT32_MAX;
        uint32_t maxLen = 0;
        uint64_t readBase = 0;
        uint64_t readNum = 0;
        float aveLen = 0.0;

        if (inputFileName1.length() > 0 && inputFileName2.length() > 0) {  // Paired-end sequencing
            // Input file stream
            gzFile gzfp1 = gzopen(inputFileName1.c_str(), "rb");
            gzFile gzfp2 = gzopen(inputFileName2.c_str(), "rb");

            // Opening a file
            if(!gzfp1 || !gzfp2) {
                cerr << "[" << __func__ << "::" << getTime() << "] '"
                    << inputFileName1 << "' or '" << inputFileName2
                    << "': No such file or directory." 
                    << endl;
                exit(1);
            } else {
                kseq_t *ks1;
                kseq_t *ks2;
                ks1 = kseq_init(gzfp1);
                ks2 = kseq_init(gzfp2);
            
                while(kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0)
                {
                    uint32_t length1 = ks1->seq.l;
                    uint32_t length2 = ks2->seq.l;

                    minLen = min(minLen, min(length1, length2));
                    maxLen = max(minLen, max(length1, length2));

                    readBase += length1 + length2;
                    readNum += 2;
                }

                // Release memory and close the file
                kseq_destroy(ks1);
                kseq_destroy(ks2);
                gzclose(gzfp1);
                gzclose(gzfp2);
            }
        } else if (inputFileName1.length() > 0 && inputFileName2.length() == 0) {  // Single-end sequencing
            // Input file stream
            gzFile gzfp = gzopen(inputFileName1.c_str(), "rb");

            // Opening a file
            if(!gzfp) {
                cerr << "[" << __func__ << "::" << getTime() << "] '"
                    << inputFileName1 
                    << "': No such file or directory." 
                    << endl;
                exit(1);
            } else {
                kseq_t *ks;
                ks = kseq_init(gzfp);
            
                while(kseq_read(ks) >= 0) {
                    uint32_t length = ks->seq.l;

                    minLen = min(minLen, length);
                    maxLen = max(minLen, length);

                    readBase += length;
                    readNum++;
                }

                // Release memory and close the file
                kseq_destroy(ks);
                gzclose(gzfp);
            }
        } else {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Input file error." 
                    << endl;
            exit(1);
        }
        aveLen = (float)readBase / readNum;

        return make_tuple(minLen, maxLen, aveLen, readBase, readNum);
    }

    // Generate configuration files
    void make_config(const uint32_t & aveLen, const uint32_t & insLen, const string & inputFileName1, const string & inputFileName2, const int & sequenceType) {
        /*
        max_rd_len=150
        [LIB]
        avg_ins=200
        reverse_seq=0
        asm_flags=3
        rd_len_cutoff=100
        rank=1
        q1=/home/hujianbing/dzz/test/PN02/0.85-new-new/sample.PN02_R1_clean.fq.gz
        q2=/home/hujianbing/dzz/test/PN02/0.85-new-new/sample.PN02_R2_clean.fq.gz
        */

       // log
        cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Build configuration file." 
                    << endl;

       string outTxt;
       outTxt += "max_rd_len=" + to_string(aveLen) + 
                 "\n[LIB]\navg_ins=" + to_string(insLen) + 
                 "\nreverse_seq=" + to_string(sequenceType) + 
                 "\nasm_flags=3\nrd_len_cutoff=100\nrank=1\nq1=" + inputFileName1 + 
                 "\nq2=" + inputFileName2 + "\n";

        // Output file stream
        ofstream outputFile;
        outputFile.open("config_A", ios::out);
        outputFile << outTxt;

        // Close File
        outputFile.close();
    }

    // SOAPdenovo
    void SOAPdenovo(const string & SOAPdenovoPath, const int & kmerLen, const string & prefix, const int & threadsNum) {
        // log
        cerr << "[" << __func__ << "::" << getTime() << "] "
             << "Running SOAPdenovo-63mer." 
             << endl;

        shell::IN string cmd;
        if (SOAPdenovoPath.empty()) {
            cmd = "SOAPdenovo-63mer all -s ./config_A -K " + to_string(kmerLen) + " -R -o " + prefix + " -p " + to_string(threadsNum);
        } else {
            cmd = SOAPdenovoPath + " all -s ./config_A -K " + to_string(kmerLen) + " -R -o " + prefix + " -p " + to_string(threadsNum);
        }
        
        // Execute shell script
        int cmdReturnValue = 0;
        bool exitBool; // Program exit status
        string exitTxt; // Program exit string
        tie(exitBool, exitTxt) = shell::exeShellCmd(cmd, & cmdReturnValue);

        // Check whether the program runs successfully
        if (!exitBool) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "SOAPdenovo-63mer error: " << exitBool << " " << exitTxt 
                 << endl;
            exit(1); 
        }
    }

    // minimap2
    void minimap2(const string & minimap2Path, const string & inputBaseFile, const string & prefix, const int & threadsNum) {
        // log
        cerr << "[" << __func__ << "::" << getTime() << "] "
             << "Running minimap2." 
             << endl;

        shell::IN string cmd;
        if (minimap2Path.empty())
        {
            cmd = "minimap2 -ax asm5 -t " + to_string(threadsNum) + " " + inputBaseFile + " " + prefix + ".scafSeq > " + prefix + ".sam";
        }
        else
        {
            cmd = minimap2Path + " -ax asm5 -t " + to_string(threadsNum) + " " + inputBaseFile + " " + prefix + ".scafSeq > " + prefix + ".sam";
        }
        
        // Execute shell script
        int cmdReturnValue = 0;
        bool exitBool; // Program exit status
        string exitTxt; // Program exit string
        tie(exitBool, exitTxt) = shell::exeShellCmd(cmd, & cmdReturnValue);

        // Check whether the program runs successfully
        if (!exitBool) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "minimap2 error: " << exitBool << " " << exitTxt 
                 << endl;
            exit(1); 
        }

        // awk -F '\t' 'BEGIN{FS=OFS="\t"} {if($3!~/*/ && $0!~/@/) {print ">"$1" Haplotype:"$3";MAPQ:"$5";CIGAR:"$6"\n"$10} else {pass}}' out.sam > out.fa
        // Convert sam file to fasta
        cmd = "awk -F '\\t' 'BEGIN{FS=OFS=\"\\t\"} {if($3!~/*/ && $0!~/@/) {print \">\"$1\" Haplotype:\"$3\";MAPQ:\"$5\";CIGAR:\"$6\"\\n\"$10} else {pass}}' " + prefix + ".sam" + " > " + prefix + ".fa";
        
        // Execute shell script
        tie(exitBool, exitTxt) = shell::exeShellCmd(cmd, & cmdReturnValue);

        // Check whether the program runs successfully
        if (!exitBool) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "awk error: " << exitBool << " " << exitTxt 
                 << endl;
            exit(1); 
        }
    }
}

#endif