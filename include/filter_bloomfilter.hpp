#ifndef filter_bloomfilter_hpp
#define filter_bloomfilter_hpp
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <thread>
#include "zlib.h"
#include <malloc.h>
#include <regex>
#include "kmer.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "ThreadPool.hpp"
#include "bitmap.hpp"
#include "bloomfilter.hpp"
#include <mutex>

std::mutex mtx;

using namespace std;

// kseq.h Opening a file
KSEQ_INIT(gzFile, gzread)

namespace filter_bloomfilter{
    /*#############################################################################################*/
    vector<string> fasta_open(const string & inputFileName, const int & kmerLen, const int & threadsNum);

    // Open fastq.gz file Paired-end sequencing
    int fastq_open_p(
        const string & fastqFileName1, 
        const string & fastqFileName2, 
        const int & kmerLen, 
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        const int & threadsNum,
        const string & prefix
    );

    // Fastq.gz multithreaded function double-end sequencing
    int fastq_filter_p(
        const int & kmerLen,
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        string readName1,
        string readName2,
        string readSeq1,
        string readSeq2,
        string readQual1, 
        string readQual2, 
        string & outRead1, 
        string & outRead2, 
        gzFile & gzfp1, 
        gzFile & gzfp2
    );

    // Open fastq.gz file Single-end sequencing
    int fastq_open_s(
        const string & fastqFileName, 
        const int & kmerLen, 
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        const int & threadsNum,
        const string & prefix
    );

    // Fastq.gz multithreaded function single end sequencing
    int fastq_filter_s(
        const int & kmerLen,
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        string readName,
        string readSeq,
        string readQual, 
        string & outRead, 
        gzFile & gzfp
    );

    // Open the fasta.gz file and perform double-end sequencing
    int fasta_open_p(
        const string & fastaFileName1, 
        const string & fastaFileName2, 
        const int & kmerLen, 
        const float & matchTd, 
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9,
        HashFun10, HashFun11, HashFun12> & bf,
        const int & threadsNum, 
        const string & prefix
    );

    // fasta.gz multithreaded function, double-end sequencing
    int fasta_filter_p(
        const int & kmerLen,
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        string readName1,
        string readName2,
        string readSeq1,
        string readSeq2,
        string & outRead1, 
        string & outRead2, 
        gzFile & gzfpO1,
        gzFile & gzfpO2
    );

    //Open the fasta.gz file and perform single-end sequencing
    int fasta_open_s(
        const string & fastaFileName, 
        const int & kmerLen, 
        const float & matchTd, 
        BloomFilter<string, HashFun1, HashFun2, HashFun3,
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf, 
        const int & threadsNum, 
        const string & prefix
    );

    // fasta.gz multithreaded function, single-end sequencing
    int fasta_filter_s(
        const int & kmerLen,
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3,
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        string readName,
        string readSeq,
        string & outRead, 
        gzFile & gzfpO
    );

    /*#############################################################################################*/


    vector<string> build_base_index(
        const string & inputFileName, 
        const int & kmerLen, 
        const int & threadsNum
    ) {
        // Check if it is a fasta file
        if (
            inputFileName.find(".fa") == string::npos &&
            inputFileName.find(".FA") == string::npos &&
            inputFileName.find(".fasta") == string::npos &&
            inputFileName.find(".FASTA") == string::npos
        ) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "The file must be in fasta format. "
                 << "(.fa / .fasta / .FA / .FASTA)\n";
            exit(1);
        }

        // Save kmers of all sequences
        vector<string> allKmerVecBase;
            
        allKmerVecBase = fasta_open(inputFileName, kmerLen, threadsNum);

        return allKmerVecBase;
    }

    int build_sequence_index(
        const string & fastqFileName1, 
        const string & fastqFileName2, 
        const string & inputFastaFile1, 
        const string & inputFastaFile2,
        const int & kmerLen,
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        const int & threadsNum,
        const string & prefix
    ) {
        // Paired-end sequencing
        if (
            fastqFileName1.length() > 0 && fastqFileName2.length() > 0 
            && inputFastaFile1.empty() && inputFastaFile2.empty()
        ) {
            fastq_open_p(
                fastqFileName1, 
                fastqFileName2, 
                kmerLen, 
                matchTd,
                bf,
                threadsNum,
                prefix
            );
        } else if (  // Single-end sequencing
            fastqFileName1.length() > 0 && fastqFileName2.empty() && inputFastaFile1.empty() && inputFastaFile2.empty()
        ) {
            fastq_open_s(fastqFileName1, 
                         kmerLen, 
                         matchTd,
                         bf,
                         threadsNum,
                         prefix);
        } else if (
            fastqFileName1.empty() && fastqFileName2.empty() && inputFastaFile1.size() > 0 && inputFastaFile2.empty()
        ) {
            fasta_open_s(inputFastaFile1, 
                         kmerLen, 
                         matchTd,
                         bf,
                         threadsNum,
                         prefix);
        } else if (
            fastqFileName1.empty() && fastqFileName2.empty() && inputFastaFile1.size() > 0 && inputFastaFile2.size() > 0
        ) {
            fasta_open_p(inputFastaFile1, 
                         inputFastaFile2, 
                         kmerLen, 
                         matchTd,
                         bf,
                         threadsNum,
                         prefix);
        } else {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                 << "Please submit sequencing files."
                 << endl;
            exit(1);
        }
        

        return 0;
    }
    
    vector<string> fasta_open(
        const string & inputFileName, 
        const int & kmerLen,  
        const int & threadsNum
    ) {
        // print log
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Collect kmer.\n";

        // Record all kmers
        vector<string> allKmerVecBase;

        // open base fastafile
        gzFile gzfp = gzopen(inputFileName.c_str(), "rb");

        // Opening a file
        if(!gzfp) {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                    << inputFileName 
                    << "': No such file or directory." 
                    << endl;
            exit(1);
        } else {
            kseq_t *ks;
            ks = kseq_init(gzfp);
        
            while( kseq_read(ks) >= 0) {
                // ks->name.s The name is recorded
                // ks->seq.s The sequence is recorded
                // Capitalize
                string sequence = ks->seq.s;
                transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

                // Building kmer index
                vector<string> kmerVecBase = build_kmer_index(sequence, kmerLen);

                // Merge the kmer of this sequence into the total vector
                allKmerVecBase.insert(allKmerVecBase.end(),kmerVecBase.begin(),kmerVecBase.end());
            }

            // Release memory and close the file
            kseq_destroy(ks);
            gzclose(gzfp);
        }

        // Vector deduplication
        sort(allKmerVecBase.begin(), allKmerVecBase.end());
        allKmerVecBase.erase(unique(allKmerVecBase.begin(), allKmerVecBase.end()), allKmerVecBase.end());

        return allKmerVecBase;
    }


    // Open fastq.gz file Paired-end sequencing
    int fastq_open_p(
        const string & fastqFileName1,
        const string & fastqFileName2, 
        const int & kmerLen, 
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6,
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        const int & threadsNum,
        const string & prefix
    ) {
        // Output File
        string outputFile1 = prefix + "." +  split(fastqFileName1, "/")[split(fastqFileName1, "/").size()-1];
        if (outputFile1.find(".gz") == string::npos && outputFile1.find(".GZ") == string::npos) {
            outputFile1 += ".gz";
        }
        string outputFile2 = prefix + "." +  split(fastqFileName2, "/")[split(fastqFileName2, "/").size()-1];
        if (outputFile2.find(".gz") == string::npos && outputFile2.find(".GZ") == string::npos) {
            outputFile2 += ".gz";
        }
        // Output file stream
        gzFile gzfpO1 = gzopen(outputFile1.c_str(), "wb");
        gzFile gzfpO2 = gzopen(outputFile2.c_str(), "wb");

        if(!gzfpO1 || !gzfpO2) {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << outputFile1 << "' or '" << outputFile2
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        }

        // Save reads whose matching degree is greater than the threshold
        string outRead1;
        string outRead2;

        // Input file stream
        gzFile gzfpI1 = gzopen(fastqFileName1.c_str(), "rb");
        gzFile gzfpI2 = gzopen(fastqFileName2.c_str(), "rb");

        // Process Pool
        ThreadPool pool(threadsNum);

        // Initialize the thread pool
        pool.init();

        // Opening a file 
        if(!gzfpI1 || !gzfpI2) {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << fastqFileName1 << "' or '" << fastqFileName2
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        } else {
            kseq_t *ks1;
            ks1 = kseq_init(gzfpI1);
            kseq_t *ks2;
            ks2 = kseq_init(gzfpI2);
        
            while(kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
                string readName1 = ks1->name.s;
                readName1 = "@" + readName1;

                string readName2 = ks2->name.s;
                readName2 = "@" + readName2;

                // Determine whether fastq has adapter sequence information
                string readNameLast1;
                if (ks1->comment.s) {
                    readNameLast1 = ' ' + string(ks1->comment.s);
                    readName1 += readNameLast1; // read name plus connector information
                } else {
                    readNameLast1 = "";
                    readName1 += " 1"; // read name plus label
                }

                string readNameLast2;
                if (ks2->comment.s) {
                    readNameLast2 = ' ' + string(ks2->comment.s);
                    readName2 += readNameLast2; // read name plus connector information
                } else {
                    readNameLast2 = "";
                    readName2 += " 2"; // read name plus label
                }

                // The value must be reassigned here, otherwise multithreading will report an error.
                string readSeq1 = ks1->seq.s;
                string readQual1 = ks1->qual.s;
                string readSeq2 = ks2->seq.s;
                string readQual2 = ks2->qual.s;

                // Capitalize
                transform(readSeq1.begin(),readSeq1.end(),readSeq1.begin(),::toupper);
                transform(readSeq2.begin(),readSeq2.end(),readSeq2.begin(),::toupper);
                
                pool.submit(
                    fastq_filter_p, 
                    kmerLen, 
                    matchTd,
                    ref(bf),
                    readName1,
                    readName2,
                    readSeq1,
                    readSeq2,
                    readQual1,
                    readQual2,
                    ref(outRead1), 
                    ref(outRead2), 
                    ref(gzfpO1), 
                    ref(gzfpO2)
                );

                // Clear the string and release the memory            
                readName1.clear();
                readNameLast1.clear();
                readSeq1.clear();
                readQual1.clear();
                string().swap(readName1);
                string().swap(readNameLast1);
                string().swap(readSeq1);
                string().swap(readQual1);
                readName2.clear();
                readNameLast2.clear();
                readSeq2.clear();
                readQual2.clear();
                string().swap(readName2);
                string().swap(readNameLast2);
                string().swap(readSeq2);
                string().swap(readQual2);

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

            // Release memory and close the file
            kseq_destroy(ks1);
            gzclose(gzfpI1);
            kseq_destroy(ks2);
            gzclose(gzfpI2);
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

        // Write the file one last time
        if (outRead1.length() > 0 || outRead2.length() > 0) {
           gzwrite(gzfpO1, outRead1.c_str(), outRead1.length());
           gzwrite(gzfpO2, outRead2.c_str(), outRead2.length());
        }
        gzclose(gzfpO1);
        gzclose(gzfpO2);

        // Clear a string
        outRead1.clear();
        string().swap(outRead1);
        outRead2.clear();
        string().swap(outRead2);

        // Freeing up memory
        malloc_trim(0);	// 0 is for heap memory
        
        return 0;
    }

    // Fastq.gz multithreaded function double-end sequencing
    int fastq_filter_p(
        const int & kmerLen,
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6,
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        string readName1,
        string readName2,
        string readSeq1,
        string readSeq2,
        string readQual1, 
        string readQual2, 
        string & outRead1, 
        string & outRead2, 
        gzFile & gzfpO1, 
        gzFile & gzfpO2
    ) {
        // Fastq.gz multithreaded function double-end sequencing
        vector<string> kmerVector1 = build_kmer_index(readSeq1, kmerLen);
        vector<string> kmerVector2 = build_kmer_index(readSeq2, kmerLen);

        // Vector deduplication
        sort(kmerVector1.begin(), kmerVector1.end());
        kmerVector1.erase(unique(kmerVector1.begin(), kmerVector1.end()), kmerVector1.end());
        sort(kmerVector2.begin(), kmerVector2.end());
        kmerVector2.erase(unique(kmerVector2.begin(), kmerVector2.end()), kmerVector2.end());

        // See if read1 is in the Bloom filter
        int kmerVectorSize1 = kmerVector1.size();
        float minMatchKmerNum1 = kmerVectorSize1*matchTd;
        float maxMisMatchKmerNum1 = kmerVectorSize1*(1-matchTd);
        int matchKmerNum1 = 0;
        int misMatchKmerNum1 = 0;

        for (auto it : kmerVector1) {
            // First determine whether the number of misMatchKmerNum exceeds the maximum value
            // If it exceeds the limit, exit the function directly and judge the situation of read2
            if (misMatchKmerNum1 > maxMisMatchKmerNum1) {
                // Clear Memory
                vector<string>().swap(kmerVector1);
                break;
            }

            if (bf.find(it)) {  // If it exists, the function returns true
                matchKmerNum1++; // The number of kmers matched is increased by 1
            } else {
                misMatchKmerNum1++; // The number of unmatched kmers is increased by 1
            }
            
            // Determine whether the number of MatchKmerNum1 exceeds the minimum value
            // If it exceeds the limit, it will be added directly to outRead1 and the function will be exited.
            if (matchKmerNum1 >= minMatchKmerNum1) {
                // Multithreaded data lock
                std::lock_guard<std::mutex> mtx_locker(mtx);

                // Print progress and save results to map
                cerr << "[" << __func__ << "::" << getTime() << "] pass: " << readName1 << endl;
                outRead1 += readName1 + "\n" + readSeq1 + "\n+\n" + readQual1 + "\n";
                outRead2 += readName2 + "\n" + readSeq2 + "\n+\n" + readQual2 + "\n";

                // Clear the string and release the memory            
                readName1.clear();
                readSeq1.clear();
                readQual1.clear();
                string().swap(readName1);
                string().swap(readSeq1);
                string().swap(readQual1);
                readName2.clear();
                readSeq2.clear();
                readQual2.clear();
                string().swap(readName2);
                string().swap(readSeq2);
                string().swap(readQual2);

                // Check if outRead1 is greater than 10Mb, if so, output and clear outRead1 and outRead2
                if (outRead1.size() > 10 * 1024 * 1024) {
                    gzwrite(gzfpO1, outRead1.c_str(), outRead1.length());
                    gzwrite(gzfpO2, outRead2.c_str(), outRead2.length());

                    // Clear a string
                    outRead1.clear();
                    string().swap(outRead1);
                    outRead2.clear();
                    string().swap(outRead2);
                }

                return 0;
            }
        }

        // If read1 is not in the bloom filter, check read2
        int kmerVectorSize2 = kmerVector2.size();
        float minMatchKmerNum2 = kmerVectorSize2*matchTd;
        float maxMisMatchKmerNum2 = kmerVectorSize2*(1-matchTd);
        int matchKmerNum2 = 0;
        int misMatchKmerNum2 = 0;

        for (auto it : kmerVector2) {
            // First determine whether the number of misMatchKmerNum exceeds the maximum value
            // If it exceeds the limit, exit the function directly without going to the next step
            if (misMatchKmerNum2 > maxMisMatchKmerNum2) {
                // Clear Memory
                vector<string>().swap(kmerVector2);
                return 0;
            }

            if (bf.find(it)) {  // If it exists, the function returns true
                matchKmerNum2++; // The number of kmers matched is increased by 1
            } else {
                misMatchKmerNum2++; // The number of unmatched kmers is increased by 1
            }
            
            // Determine whether the number of MatchKmerNum2 exceeds the minimum value
            // If it exceeds the limit, it will be added directly to outRead2 and the function will be exited.
            if (matchKmerNum2 >= minMatchKmerNum2) {
                // Multithreaded data lock
                std::lock_guard<std::mutex> mtx_locker(mtx);

                // Print progress and save results to map
                cerr << "[" << __func__ << "::" << getTime() << "] pass: " << readName2 << endl;
                outRead1 += readName1 + "\n" + readSeq1 + "\n+\n" + readQual1 + "\n";
                outRead2 += readName2 + "\n" + readSeq2 + "\n+\n" + readQual2 + "\n";

                // Clear the string and release the memory            
                readName1.clear();
                readSeq1.clear();
                readQual1.clear();
                string().swap(readName1);
                string().swap(readSeq1);
                string().swap(readQual1);
                readName2.clear();
                readSeq2.clear();
                readQual2.clear();
                string().swap(readName2);
                string().swap(readSeq2);
                string().swap(readQual2);

                // Check if outRead1 is greater than 10Mb, if so, output and clear outRead1 and outRead2
                if (outRead1.size() > 10 * 1024 * 1024) {
                    gzwrite(gzfpO1, outRead1.c_str(), outRead1.length());
                    gzwrite(gzfpO2, outRead2.c_str(), outRead2.length());

                    // Clear a string
                    outRead1.clear();
                    string().swap(outRead1);
                    outRead2.clear();
                    string().swap(outRead2);
                }

                return 0;
            }
        }

        // Clear Memory
        vector<string>().swap(kmerVector1);
        vector<string>().swap(kmerVector2);

        // Freeing up memory
        malloc_trim(0);	// 0 is for heap memory

        return 0;
    }


    // Open fastq.gz file Single-end sequencing
    int fastq_open_s(
        const string & fastqFileName, 
        const int & kmerLen, 
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        const int & threadsNum,
        const string & prefix
    ) {
        // Output File
        string outputFile = prefix + "." +  split(fastqFileName, "/")[split(fastqFileName, "/").size()-1];
        if (outputFile.find(".gz") == string::npos && outputFile.find(".GZ") == string::npos) {
            outputFile += ".gz";
        }
        // Output file stream
        gzFile gzfpO = gzopen(outputFile.c_str(), "wb");

        if(!gzfpO) {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << outputFile 
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        }

        // Save reads whose matching degree is greater than the threshold
        string outRead;

        // Input file stream
        gzFile gzfpI = gzopen(fastqFileName.c_str(), "rb");

        // Process Pool
        ThreadPool pool(threadsNum);

        // Initialize the thread pool
        pool.init();

        // Opening a file 
        if(!gzfpI) {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << fastqFileName
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        } else {
            kseq_t *ks;
            ks = kseq_init(gzfpI);
        
            while(kseq_read(ks) >= 0) {
                string readName = ks->name.s;
                readName = "@" + readName;

                // Determine whether fastq has adapter sequence information
                string readNameLast;
                if (ks->comment.s) {
                    readNameLast = ' ' + string(ks->comment.s);
                } else {
                    readNameLast = "";
                }

                // Read name plus label and connector information
                readName += readNameLast;

                // The value must be reassigned here, otherwise multithreading will report an error.
                string readSeq = ks->seq.s;
                string readQual = ks->qual.s;

                // Capitalize
                transform(readSeq.begin(),readSeq.end(),readSeq.begin(),::toupper);
                
                pool.submit(
                    fastq_filter_s, 
                    kmerLen, 
                    matchTd,
                    ref(bf),
                    readName,
                    readSeq,
                    readQual,
                    ref(outRead), 
                    ref(gzfpO)
                );

                // Clear the string and release the memory            
                readName.clear();
                readNameLast.clear();
                readSeq.clear();
                readQual.clear();
                string().swap(readName);
                string().swap(readNameLast);
                string().swap(readSeq);
                string().swap(readQual);

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

            // Release memory and close the file
            kseq_destroy(ks);
            gzclose(gzfpI);
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

        // Write the file one last time
        if (outRead.length() > 0) {
           gzwrite(gzfpO, outRead.c_str(), outRead.length());
        }
        gzclose(gzfpO);

        // Clear a string
        outRead.clear();
        string().swap(outRead);

        // Freeing up memory
        malloc_trim(0);	// 0 is for heap memory
        
        return 0;
    }

    // Fastq.gz multithreaded function single end sequencing
    int fastq_filter_s(
        const int & kmerLen,
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        string readName,
        string readSeq,
        string readQual, 
        string & outRead, 
        gzFile & gzfpO
    ) {
        // Building kmer index
        vector<string> kmerVector = build_kmer_index(readSeq, kmerLen);

        // Vector deduplication
        sort(kmerVector.begin(), kmerVector.end());
        kmerVector.erase(unique(kmerVector.begin(), kmerVector.end()), kmerVector.end());

        // See if read is in the Bloom filter
        int kmerVectorSize = kmerVector.size();
        float minMatchKmerNum = kmerVectorSize*matchTd;
        float maxMisMatchKmerNum = kmerVectorSize*(1-matchTd);
        int matchKmerNum = 0;
        int misMatchKmerNum = 0;

        for (auto it : kmerVector) {
            // First determine whether the number of misMatchKmerNum exceeds the maximum value
            // If it exceeds the limit, exit the function directly without going to the next step
            if (misMatchKmerNum > maxMisMatchKmerNum) {
                // Clear Memory
                vector<string>().swap(kmerVector);
                return 0;
            }

            if (bf.find(it)) {  // If it exists, the function returns true
                matchKmerNum++; // The number of kmers matched is increased by 1
            } else {
                misMatchKmerNum++; // The number of unmatched kmers is increased by 1
            }
            
            // Determine whether the number of MatchKmerNum exceeds the minimum value
            // If it exceeds the limit, it will be added directly to outRead and the function will be exited.
            if (matchKmerNum >= minMatchKmerNum) {
                // Multithreaded data lock
                std::lock_guard<std::mutex> mtx_locker(mtx);

                // Print progress and save results to map
                cerr << "[" << __func__ << "::" << getTime() << "] pass: " << readName << endl;
                outRead += readName + "\n" + readSeq + "\n+\n" + readQual + "\n";

                // Clear the string and release the memory            
                readName.clear();
                readSeq.clear();
                readQual.clear();
                string().swap(readName);
                string().swap(readSeq);
                string().swap(readQual);

                // Check if outRead1 is greater than 10Mb, if so, output and clear outRead
                if (outRead.size() > 10 * 1024 * 1024) {
                    gzwrite(gzfpO, outRead.c_str(), outRead.length());

                    // Clear a string
                    outRead.clear();
                    string().swap(outRead);
                }

                return 0;
            }
        }

        // Clear Memory
        vector<string>().swap(kmerVector);

        // Freeing up memory
        malloc_trim(0);	// 0 is for heap memory

        return 0;
    }


    // Open the fasta.gz file and perform double-end sequencing
    int fasta_open_p(
        const string & fastaFileName1, 
        const string & fastaFileName2, 
        const int & kmerLen, 
        const float & matchTd, 
        BloomFilter<string, HashFun1, HashFun2, HashFun3,
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf, 
        const int & threadsNum, 
        const string & prefix
    ) {
        // Output File
        string outputFile1 = prefix + "." +  split(fastaFileName1, "/")[split(fastaFileName1, "/").size()-1];
        if (outputFile1.find(".gz") == string::npos && outputFile1.find(".GZ") == string::npos) {
            outputFile1 += ".gz";
        }
        string outputFile2 = prefix + "." +  split(fastaFileName2, "/")[split(fastaFileName2, "/").size()-1];
        if (outputFile2.find(".gz") == string::npos && outputFile2.find(".GZ") == string::npos) {
            outputFile2 += ".gz";
        }
        // Output file stream
        gzFile gzfpO1 = gzopen(outputFile1.c_str(), "wb");
        gzFile gzfpO2 = gzopen(outputFile2.c_str(), "wb");

        if(!gzfpO1 || !gzfpO2) {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << outputFile1 << "' or '" << outputFile2
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        }

        // Save reads whose matching degree is greater than the threshold
        string outRead1;
        string outRead2;

        // Input file stream
        gzFile gzfpI1 = gzopen(fastaFileName1.c_str(), "rb");
        gzFile gzfpI2 = gzopen(fastaFileName2.c_str(), "rb");

        // Process Pool
        ThreadPool pool(threadsNum);

        // Initialize the thread pool
        pool.init();

        // Opening a file 
        if(!gzfpI1 || !gzfpI2) {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << fastaFileName1 << "' or '" << fastaFileName2
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        } else {
            kseq_t *ks1;
            ks1 = kseq_init(gzfpI1);
            kseq_t *ks2;
            ks2 = kseq_init(gzfpI2);
        
            while(kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
                string readName1 = ks1->name.s;
                readName1 = ">" + readName1;

                string readName2 = ks2->name.s;
                readName2 = ">" + readName2;

                // Determine whether fastq has adapter sequence information
                string readNameLast1;
                if (ks1->comment.s) {
                    readNameLast1 = ' ' + string(ks1->comment.s);
                    readName1 += readNameLast1; // read name plus connector information
                } else {
                    readNameLast1 = "";
                    readName1 += " 1"; // read name plus label
                }

                string readNameLast2;
                if (ks2->comment.s) {
                    readNameLast2 = ' ' + string(ks2->comment.s);
                    readName2 += readNameLast2; // read name plus connector information
                } else {
                    readNameLast2 = "";
                    readName2 += " 2"; // read name plus label
                }

                // The value must be reassigned here, otherwise multithreading will report an error.
                string readSeq1 = ks1->seq.s;
                string readSeq2 = ks2->seq.s;

                // Capitalize
                transform(readSeq1.begin(),readSeq1.end(),readSeq1.begin(),::toupper);
                transform(readSeq2.begin(),readSeq2.end(),readSeq2.begin(),::toupper);
                
                pool.submit(
                    fasta_filter_p, 
                    kmerLen, 
                    matchTd,
                    ref(bf),
                    readName1,
                    readName2,
                    readSeq1,
                    readSeq2,
                    ref(outRead1), 
                    ref(outRead2), 
                    ref(gzfpO1), 
                    ref(gzfpO2)
                );

                // Clear the string and release the memory            
                readName1.clear();
                readNameLast1.clear();
                readSeq1.clear();
                string().swap(readName1);
                string().swap(readNameLast1);
                string().swap(readSeq1);
                readName2.clear();
                readNameLast2.clear();
                readSeq2.clear();
                string().swap(readName2);
                string().swap(readNameLast2);
                string().swap(readSeq2);

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

            // Release memory and close the file
            kseq_destroy(ks1);
            gzclose(gzfpI1);
            kseq_destroy(ks2);
            gzclose(gzfpI2);
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

        // Write the file one last time
        if (outRead1.length() > 0 || outRead2.length() > 0)  {
           gzwrite(gzfpO1, outRead1.c_str(), outRead1.length());
           gzwrite(gzfpO2, outRead2.c_str(), outRead2.length());
        }
        gzclose(gzfpO1);
        gzclose(gzfpO2);

        // Clear a string
        outRead1.clear();
        string().swap(outRead1);
        outRead2.clear();
        string().swap(outRead2);

        // Freeing up memory
        malloc_trim(0);	// 0 is for heap memory
        
        return 0;
    }

    // fasta.gz multithreading function
    int fasta_filter_p(
        const int & kmerLen, 
        const float & matchTd, 
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf, 
        string readName1, 
        string readName2, 
        string readSeq1, 
        string readSeq2, 
        string & outRead1, 
        string & outRead2, 
        gzFile & gzfpO1, 
        gzFile & gzfpO2
    ) {
        // Building kmer index
        vector<string> kmerVector1 = build_kmer_index(readSeq1, kmerLen);
        vector<string> kmerVector2 = build_kmer_index(readSeq2, kmerLen);

        // Vector deduplication
         sort(kmerVector1.begin(), kmerVector1.end());
        kmerVector1.erase(unique(kmerVector1.begin(), kmerVector1.end()), kmerVector1.end());
        sort(kmerVector2.begin(), kmerVector2.end());
        kmerVector2.erase(unique(kmerVector2.begin(), kmerVector2.end()), kmerVector2.end());


        // See if read1 is in the Bloom filter
        int kmerVectorSize1 = kmerVector1.size();
        float minMatchKmerNum1 = kmerVectorSize1*matchTd;
        float maxMisMatchKmerNum1 = kmerVectorSize1*(1-matchTd);
        int matchKmerNum1 = 0;
        int misMatchKmerNum1 = 0;

        for (auto it : kmerVector1) {
            // First determine whether the number of misMatchKmerNum exceeds the maximum value
            // If it exceeds the limit, exit the function directly and judge the situation of read2
            if (misMatchKmerNum1 > maxMisMatchKmerNum1) {
                // Clear Memory
                vector<string>().swap(kmerVector1);
                break;
            }

            if (bf.find(it)) {  // If it exists, the function returns true
                matchKmerNum1++; // The number of kmers matched is increased by 1
            } else {
                misMatchKmerNum1++; // The number of unmatched kmers is increased by 1
            }
            
            // Determine whether the number of MatchKmerNum1 exceeds the minimum value
            // If it exceeds the limit, it will be added directly to outRead1 and the function will be exited.
            if (matchKmerNum1 >= minMatchKmerNum1) {
                // Multithreaded data lock
                std::lock_guard<std::mutex> mtx_locker(mtx);

                // Print progress and save results to map
                cerr << "[" << __func__ << "::" << getTime() << "] pass: " << readName1 << endl;
                outRead1 += readName1 + "\n" + readSeq1 + "\n";
                outRead2 += readName2 + "\n" + readSeq2 + "\n";

                // Clear the string and release the memory            
                readName1.clear();
                readSeq1.clear();
                string().swap(readName1);
                string().swap(readSeq1);
                readName2.clear();
                readSeq2.clear();
                string().swap(readName2);
                string().swap(readSeq2);

                // Check if outRead1 is greater than 10Mb, if so, output and clear outRead1 and outRead2
                if (outRead1.size() > 10 * 1024 * 1024) {
                    gzwrite(gzfpO1, outRead1.c_str(), outRead1.length());
                    gzwrite(gzfpO2, outRead2.c_str(), outRead2.length());

                    // Clear a string
                    outRead1.clear();
                    string().swap(outRead1);
                    outRead2.clear();
                    string().swap(outRead2);
                }

                return 0;
            }
        }

        // If read1 is not in the bloom filter, check read2
        int kmerVectorSize2 = kmerVector2.size();
        float minMatchKmerNum2 = kmerVectorSize2*matchTd;
        float maxMisMatchKmerNum2 = kmerVectorSize2*(1-matchTd);
        int matchKmerNum2 = 0;
        int misMatchKmerNum2 = 0;

        for (auto it : kmerVector2) {
            // First determine whether the number of misMatchKmerNum exceeds the maximum value
            // If it exceeds the limit, exit the function directly without going to the next step
            if (misMatchKmerNum2 > maxMisMatchKmerNum2) {
                // Clear Memory
                vector<string>().swap(kmerVector2);
                return 0;
            }

            if (bf.find(it)) {  // If it exists, the function returns true
                matchKmerNum2++; // The number of kmers matched is increased by 1
            } else {
                misMatchKmerNum2++; // The number of unmatched kmers is increased by 1
            }
            
            // Determine whether the number of MatchKmerNum2 exceeds the minimum value
            // If it exceeds the limit, it will be added directly to outRead2 and the function will be exited.
            if (matchKmerNum2 >= minMatchKmerNum2) {
                // Multithreaded data lock
                std::lock_guard<std::mutex> mtx_locker(mtx);

                // Print progress and save results to map
                cerr << "[" << __func__ << "::" << getTime() << "] pass: " << readName2 << endl;
                outRead1 += readName1 + "\n" + readSeq1 + "\n";
                outRead2 += readName2 + "\n" + readSeq2 + "\n";

                // Clear the string and release the memory            
                readName1.clear();
                readSeq1.clear();
                string().swap(readName1);
                string().swap(readSeq1);
                readName2.clear();
                readSeq2.clear();
                string().swap(readName2);
                string().swap(readSeq2);

                // Check if outRead1 is greater than 10Mb, if so, output and clear outRead1 and outRead2
                if (outRead1.size() > 10 * 1024 * 1024) {
                    gzwrite(gzfpO1, outRead1.c_str(), outRead1.length());
                    gzwrite(gzfpO2, outRead2.c_str(), outRead2.length());

                    // Clear a string
                    outRead1.clear();
                    string().swap(outRead1);
                    outRead2.clear();
                    string().swap(outRead2);
                }

                return 0;
            }
        }

        // Clear Memory
        vector<string>().swap(kmerVector1);
        vector<string>().swap(kmerVector2);

        // Freeing up memory
        malloc_trim(0);	// 0 is for heap memory

        return 0;
    }


    // Open the fasta.gz file and perform single-end sequencing
    int fasta_open_s(
        const string & fastaFileName, 
        const int & kmerLen, 
        const float & matchTd, 
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6, 
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf, 
        const int & threadsNum, 
        const string & prefix
    ) {
        // Output File
        string outputFile = prefix + "." +  split(fastaFileName, "/")[split(fastaFileName, "/").size()-1];
        if (outputFile.find(".gz") == string::npos && outputFile.find(".GZ") == string::npos) {
            outputFile += ".gz";
        }
        // Output file stream
        gzFile gzfpO = gzopen(outputFile.c_str(), "wb");

        if(!gzfpO) {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << outputFile 
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        }

        // Save reads whose matching degree is greater than the threshold
        string outRead;

        // Input file stream
        gzFile gzfpI = gzopen(fastaFileName.c_str(), "rb");

        // Process Pool
        ThreadPool pool(threadsNum);

        // Initialize the thread pool
        pool.init();

        // Opening a file 
        if(!gzfpI) {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << fastaFileName
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        } else {
            kseq_t *ks;
            ks = kseq_init(gzfpI);
        
            while(kseq_read(ks) >= 0) {
                string readName = ks->name.s;
                readName = ">" + readName;

                // Determine whether fastq has adapter sequence information
                string readNameLast;
                if (ks->comment.s) {
                    readNameLast = ' ' + string(ks->comment.s);
                } else {
                    readNameLast = "";
                }

                // Read name plus label and connector information
                readName += readNameLast;

                // The value must be reassigned here, otherwise multithreading will report an error.
                string readSeq = ks->seq.s;

                // Capitalize
                transform(readSeq.begin(),readSeq.end(),readSeq.begin(),::toupper);
                
                pool.submit(
                    fasta_filter_s, 
                    kmerLen, 
                    matchTd,
                    ref(bf),
                    readName,
                    readSeq,
                    ref(outRead),  
                    ref(gzfpO)
                );

                // Clear the string and release the memory            
                readName.clear();
                readNameLast.clear();
                readSeq.clear();
                string().swap(readName);
                string().swap(readNameLast);
                string().swap(readSeq);

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

            // Release memory and close the file
            kseq_destroy(ks);
            gzclose(gzfpI);
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

        // Write the file one last time
        if (outRead.length() > 0) {
           gzwrite(gzfpO, outRead.c_str(), outRead.length());
        }
        gzclose(gzfpO);

        // Clear a string
        outRead.clear();
        string().swap(outRead);

        // Freeing up memory
        malloc_trim(0);	// 0 is for heap memory
        
        return 0;
    }

    // fasta.gz multithreading function
    int fasta_filter_s(
        const int & kmerLen,
        const float & matchTd,
        BloomFilter<string, HashFun1, HashFun2, HashFun3, 
        HashFun4, HashFun5, HashFun6,
        HashFun7, HashFun8, HashFun9, 
        HashFun10, HashFun11, HashFun12> & bf,
        string readName,
        string readSeq,
        string & outRead, 
        gzFile & gzfpO
    ) {
        // Building kmer index
        vector<string> kmerVector = build_kmer_index(readSeq, kmerLen);

        // Vector deduplication
        sort(kmerVector.begin(), kmerVector.end());
        kmerVector.erase(unique(kmerVector.begin(), kmerVector.end()), kmerVector.end());

        // See if read is in the Bloom filter
        int kmerVectorSize = kmerVector.size();
        float minMatchKmerNum = kmerVectorSize*matchTd;
        float maxMisMatchKmerNum = kmerVectorSize*(1-matchTd);
        int matchKmerNum = 0;
        int misMatchKmerNum = 0;

        for (auto it : kmerVector) {
            // First determine whether the number of misMatchKmerNum exceeds the maximum value
            // If it exceeds the limit, exit the function directly without going to the next step
            if (misMatchKmerNum > maxMisMatchKmerNum) {
                // Clear Memory
                vector<string>().swap(kmerVector);
                return 0;
            }

            if (bf.find(it)) {  // If it exists, the function returns true
                matchKmerNum++; // The number of kmers matched is increased by 1
            } else {
                misMatchKmerNum++; // The number of unmatched kmers is increased by 1
            }
            
            // Determine whether the number of MatchKmerNum exceeds the minimum value
            // If it exceeds the limit, it will be added directly to outRead and the function will be exited.
            if (matchKmerNum >= minMatchKmerNum) {
                // Multithreaded data lock
                std::lock_guard<std::mutex> mtx_locker(mtx);

                // Print progress and save results to map
                cerr << "[" << __func__ << "::" << getTime() << "] pass: " << readName << endl;
                outRead += readName + "\n" + readSeq + "\n";

                // Clear the string and release the memory            
                readName.clear();
                readSeq.clear();
                string().swap(readName);
                string().swap(readSeq);

                // Check if outRead1 is greater than 10Mb, if so, output and clear outRead
                if (outRead.size() > 10 * 1024 * 1024) {
                    gzwrite(gzfpO, outRead.c_str(), outRead.length());

                    // Clear a string
                    outRead.clear();
                    string().swap(outRead);
                }

                return 0;
            }
        }

        // Clear Memory
        vector<string>().swap(kmerVector);

        // Freeing up memory
        malloc_trim(0);	// 0 is for heap memory

        return 0;
    }
}

#endif