#ifndef bam2ratio_hpp
#define bam2ratio_hpp
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include <regex>
#include <getopt.h>
#include "htslib/sam.h"
#include <stdio.h>
#include <stdlib.h>
#include "kmer.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "ThreadPool.hpp"

using namespace std;

namespace bam2ratio {
    struct alignmentInfo {
        vector<uint32_t> refLenVec;
        vector<string> cigarVec;
    };
    
    /*
    @discussion In the CIGAR array, each element is a 32-bit integer. The
    lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
    length of a CIGAR.
    */
    string getCigar(const bam1_t *b)
    {
        uint32_t *data = (uint32_t *)bam_get_cigar(b);
        int cigarNum = b->core.n_cigar;
        stringstream ss;
        for(int i=0; i<cigarNum; i++) {
            uint32_t val = data[i];
            char op = bam_cigar_opchr(val);
            uint32_t len = bam_cigar_oplen(val);
            ss << len << op;
        }
        return ss.str();
    }

    //Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.
    char fourbits2base(uint8_t val) {
        switch(val) {
            case 1:
                return 'A';
            case 2:
                return 'C';
            case 4:
                return 'G';
            case 8:
                return 'T';
            case 15:
                return 'N';
            default:
                cerr << "ERROR: Wrong base with value "<< (int)val << endl ;
                return 'N';
        }
    }

    // get seq The sequence information is recorded in an 8-bit data structure. The first 4 bits are the previous base and the last 4 bits are the next base.
    string getSeq(const bam1_t *b) {
        uint8_t *data = bam_get_seq(b);
        int len = b->core.l_qseq;
        string s(len, '\0');
        for(int i=0; i<len; i++) {
            char base;
            if(i%2 == 1)
                base = fourbits2base(data[i/2] & 0xF); 
            else
                base = fourbits2base((data[i/2]>>4) & 0xF);
            s[i] = base;
        }
        return s;
    }

    // get seq quality
    string getQual(const bam1_t *b) {
        uint8_t *data = bam_get_qual(b);
        int len = b->core.l_qseq;
        string s(len, '\0');
        for(int i=0; i<len; i++) {
            s[i] = (char)(data[i] + 33); // Convert to printable ascii
        }
        return s;
    }

    string getAux(const bam1_t *b,const char tag[2]) {
        kstring_t res = KS_INITIALIZE;    // Need to initialize
        if(bam_aux_get_str(b,tag,&res) == 1) {  // The string buffer of kstring is not terminated by \0
            int len = ks_len(&res);
            char *ks_s = ks_str(&res);
            string s(len, '\0');
            for (int i = 0;i<len;i++ ){
                s[i] = ks_s[i];
            }
            ks_free(&res); // Release resources
            return s;
        } else  {
            cerr << "no tag :" << tag << '\n';
            ks_free(&res);
            return "";
        }
    }


    // Split bam file
    int bam_split(
        samFile *bam_in, 
        bam_hdr_t *bam_header, 
        bam1_t *aln, 
        int splitNum, 
        map<string, alignmentInfo> & readCigraMap
    ) {
        while (sam_read1(bam_in,bam_header,aln) >= 0) 
        {
            uint32_t pos = aln->core.pos ;
            string chr = "*";
            uint32_t chrLen = 0;
            if (aln->core.tid != -1) {  // There are reads that cannot be mapped to the genome
                chr = bam_header->target_name[aln->core.tid]; // config name(chromosome)
                chrLen = bam_header->target_len[aln->core.tid];
            }
                
            string queryname = bam_get_qname(aln);
            string cigar = getCigar(aln);
            string seq = getSeq(aln);
            string qual = getQual(aln);

            // Reads that do not match are skipped
            if (cigar != "*" &&  chr != "*") {
                readCigraMap[queryname].refLenVec.push_back(chrLen);
                readCigraMap[queryname].cigarVec.push_back(cigar);
            }
            
            // Return once for every splitNum reads
            if (readCigraMap.size() >= splitNum) {
                return 0;
            }
        }

        return -1;
    }


    // Regular Expression Functions
    vector<string> Search(string & str, regex & pattern) {
        vector<string> outVec;

        sregex_token_iterator p(str.begin(), str.end(), pattern, 0);
        sregex_token_iterator end;
        vector<string> vec;
        while (p != end)
            outVec.push_back(*p++);

        return outVec;
    }


    // cigar转ratio
    int cigar2ratio(
        map<string, alignmentInfo> readCigraMap, 
        string & outReads, 
        const float & ratio, 
        const string & calRefCovBool
    ) {
        // Save the read information that meets the conditions
        string outRead;

        // Regular Expressions
        regex pattern1("\\d+");
        regex pattern2("[A-Z]+");

        for (auto forIter1 : readCigraMap) {
            for (size_t i = 0; i < forIter1.second.cigarVec.size(); i++) {
                // Record comparison ratio
                uint32_t alignmentLen = 0;

                // Calculate the length of read
                uint32_t readLen = 0;

                // Compare cigar
                string cigra = forIter1.second.cigarVec[i];
                vector<string> numberVec = Search(cigra, pattern1);
                vector<string> cigarVec = Search(cigra, pattern2);

                for (size_t j = 0; j < cigarVec.size(); j++) {
                    if (cigarVec[j] != "H" && cigarVec[j] != "S") {
                        alignmentLen += stol(numberVec[j]);
                    }
                    readLen += stol(numberVec[j]);
                }
                
                // Calculate the ratio
                float refRatio = alignmentLen/(float)forIter1.second.refLenVec[i];
                float readRatio = alignmentLen/(float)readLen;

                // If the reference genome ratio is not calculated, the value is 0.
                if (calRefCovBool == "false") {
                    refRatio = 0;
                }

                if (refRatio >= ratio || readRatio >= ratio) {
                    outRead += forIter1.first + 
                               "\trefgenome:" + to_string(refRatio) + 
                               "\tread:" + to_string(readRatio) + "\n";
                    break;
                }
            }
        }
        
        // Add to total reads
        // Multithreaded data lock
        std::lock_guard<std::mutex> mtx_locker(mtx);
        outReads += outRead;
        return 0;
    }


    // Output
    int out_save(const string & outputFileName, const string & outReads) {
        if (outputFileName.empty()) {
            cout << outReads << endl;
        } else {
            // Output
            ofstream outputFile;
            outputFile.open(outputFileName, ios::out);
            outputFile << outReads;
            // Close File
            outputFile.close();
        }
        
        return 0;
    }


    // bam to ratio
    int bam2ratio(
        char * bamFileName, 
        const float & matchTd, 
        const int & readSplitNum, 
        const string & outputFileName, 
        const int & threadsNum, 
        const string & calRefCovBool
    ) {
        // Process Pool
        ThreadPool pool(threadsNum);

        // Initialize the thread pool
        pool.init();

        // Save the split read information
        map<string, alignmentInfo> readCigraMap;

        // Save reads that meet the threshold
        string outReads;

        // Open bam file
        samFile *bam_in = sam_open(bamFileName, "r"); // open bam file
        hts_idx_t *bam_index = sam_index_load(bam_in, bamFileName); //load index 
        bam_hdr_t *bam_header = sam_hdr_read(bam_in); // read header
        bam1_t *aln = NULL;
        aln = bam_init1(); //initialize an alignment

        // Split read to speed up thread running
        while (
            bam_split(bam_in, bam_header, aln, readSplitNum, readCigraMap) == 0 || readCigraMap.size() > 0
        ) {
            pool.submit(cigar2ratio, readCigraMap, ref(outReads), ref(matchTd), ref(calRefCovBool));
            
            // Freeing up memory
            readCigraMap.clear();
            map<string, alignmentInfo>().swap(readCigraMap);
            
            malloc_trim(0);	// 0 is for heap memory

            // Check if the task queue exceeds the threshold, and wait if it exceeds the threshold to prevent loading the data into memory all at once
            int maxRetries = 120; // Set the maximum number of retries, for example 120 times (60 seconds)
            int retryCount = 0;
            while (pool.get_queue() >= threadsNum*100) {
                if (retryCount >= maxRetries) {
                    cerr << "[" << __func__ << "::" << getTime() << "] Task queue exceeded threshold for too long, continuing with execution." << endl;
                    break; // Jump out of the loop and continue to execute the following code
                }
                // Check every 0.5 seconds
                sleep(0.5);
                retryCount++;
            }
        }

        // Check whether the task queue has been executed. If it is completed, close the thread pool. Otherwise, check it every 0.5s.
        int maxRetries = 120; // Set the maximum number of retries, for example 120 times (60 seconds)
        int retryCount = 0;
        while (pool.get_queue() > 0) {
            if (retryCount >= maxRetries) {
                cerr << "[" << __func__ << "::" << getTime() << "] Task queue exceeded threshold for too long, continuing with execution." << endl;
                break; // Jump out of the loop and continue to execute the following code
            }
            // Check every 0.5 seconds
            sleep(0.5);
            retryCount++;
        }

        // Shutdown thread pool
        pool.shutdown();
        
        bam_destroy1(aln); // Recycling
        bam_hdr_destroy(bam_header);
        sam_close(bam_in);

        // 输出结果
        out_save(outputFileName, outReads);

        return 0;
    }
}

#endif