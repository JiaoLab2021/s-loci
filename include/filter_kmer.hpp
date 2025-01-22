#ifndef filter_kmer_hpp
#define filter_kmer_hpp
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <thread>
#include <malloc.h>
#include "zlib.h"
#include <regex>
#include <mutex>
#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "ThreadPool.hpp"
#include "kmer_bit_filter.hpp"
#include "save.hpp"

using namespace std;

namespace filter_kmer {
    struct readInfoStructure {
        string readName1;
	    string readName2;
	    string readSeq1;
	    string readSeq2;
	    string readQual1;
	    string readQual2;

        uint32_t chromosomeId;
    };
    
    /*#############################################################################################*/
    pair<unordered_map<uint64_t, unordered_map<uint32_t, uint32_t>>, unordered_map<uint32_t, uint32_t>> build_base_index_run(
        const string & inputFileName,  
        const unsigned int & kmerLen, 
        const int & threadsNum, 
        uint32_t & readId, 
        unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> > & kmerIdMap, 
        unordered_map<uint32_t, uint32_t> & refkmerNumMap
    );

    int read_open(
        const vector<string> & inputReadVec, 
        const int & kmerLen, 
        const float & matchTd, 
        const unordered_map<uint64_t, unordered_map<uint32_t, uint32_t>> & kmerIdMap, 
        const unordered_map<uint32_t, uint32_t> & refkmerNumMap, 
        const int & threadsNum, 
        const string & prefix, 
        const bool & calRefCovBool
    );

    vector<string> read_filter(
        const int & kmerLen,
        const float & matchTd, 
        const unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> > & kmerIdMap, 
        const unordered_map<uint32_t, uint32_t> & refkmerNumMap, 
        vector<string> readNameVec, 
        vector<string> readSeqVec, 
        vector<string> readInfoVec, 
        const bool & calRefCovBool
    );
    /*#############################################################################################*/

    // Build base file index
    pair<unordered_map<uint64_t, unordered_map<uint32_t, uint32_t>>, unordered_map<uint32_t, uint32_t>> build_base_index(
        const vector<string> & inputReadVec, 
        const unsigned int & kmerLen, 
        const int & threadsNum
    ) {
        if (inputReadVec.size() == 0) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error: -r.\n";
            exit(1);
        }

        // Record the kmer of the read  map<kmerScore, map<readId, kmerNum> >
        unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> > kmerIdMap;

        // Record the number of kmers for each read to determine coverage <readId, kmerNum (number)>
        unordered_map<uint32_t, uint32_t> refkmerNumMap;

        // Record the number of reads
        uint32_t readId = 0;

        for (const auto& inputRead : inputReadVec) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Indexed files: " << inputRead << ".\n";
            build_base_index_run(
                inputRead,  
                kmerLen, 
                threadsNum, 
                readId, 
                kmerIdMap, 
                refkmerNumMap
            );
        }

        malloc_trim(0);	// 0 is for heap memory

        return make_pair(kmerIdMap, refkmerNumMap);
    }


    /**
	 * Build the index of the base file
	 *
	 * @param inputFileName    Build the index of the base file
	 * @param kmerLen          DNA sequence
	 * @param threadsNum       Number of threads
	 * @param readId           Chromosome index number
     * @param kmerIdMap        Record the kmer of the read map<kmerScore, map<readId, kmerNum> >
     * @param refkmerNumMap    Record the number of kmers for each read to determine coverage <readId, kmerNum (number)>
	 * 
	 * @return 0
	*/
    pair<unordered_map<uint64_t, unordered_map<uint32_t, uint32_t>>, unordered_map<uint32_t, uint32_t>> build_base_index_run(
        const string & inputFileName,  
        const unsigned int & kmerLen, 
        const int & threadsNum, 
        uint32_t & readId, 
        unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> > & kmerIdMap, 
        unordered_map<uint32_t, uint32_t> & refkmerNumMap
    ) {
        // open base file
        gzFile gzfp = gzopen(inputFileName.c_str(), "rb");

        // Process Pool
        ThreadPool pool(threadsNum);

        // Initialize the thread pool
        pool.init();

        // Saving multi-threaded results
        // vector<map<kmerScore, map<readId, kmerNum> > >
        vector<future<unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> > > > kmerIdMapFutureVec;

        // Opening a file
        if(!gzfp) {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'"
                << inputFileName 
                << "': No such file or directory." 
                << endl;
            exit(1);
        } else {
            kseq_t *ks;
            ks = kseq_init(gzfp);
        
            while( kseq_read(ks) >= 0 )
            {
                // ks->name.s The name is recorded
                // ks->seq.s The sequence is recorded
                // Capitalize
                string sequence = ks->seq.s;
                string name = ks->name.s;
                transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

                uint32_t chromosomeId = readId;

                // Multithreaded task submission
                kmerIdMapFutureVec.push_back(
                    pool.submit(
                        kmerBit::kmer_sketch1, 
                        chromosomeId, 
                        sequence,
                        kmerLen
                    )
                );

                // The number of kmers in the sequence
                refkmerNumMap[chromosomeId] = ks->seq.l - kmerLen + 1;

                // Clear the string and release the memory            
                sequence.clear();
                name.clear();
                string().swap(sequence);
                string().swap(name);

                // Update read number
                readId++;

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
            gzclose(gzfp);
        }

        // Multithreaded result saving
        for (auto&& kmerIdMapFuture : kmerIdMapFutureVec) {  // vector<future<unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> > > >
            unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> > kmerIdMapTmp = move(kmerIdMapFuture.get());

            for (const auto& [key1, value1] : kmerIdMapTmp) {  // unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> >
                auto& kmerIdMapKey1 = kmerIdMap[key1];  // Avoid duplicate searches
                
                for (const auto& [key2, value2] : value1) {  // unordered_map<uint32_t, uint32_t>
                    auto& kmerIdMapValue2 = kmerIdMapKey1[key2];  // Avoid duplicate searches

                    kmerIdMapValue2 += value2;  // Directly perform addition operations
                }
            }
        }
        kmerIdMapFutureVec.clear();
        vector<future<unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> > > >().swap(kmerIdMapFutureVec);
        malloc_trim(0);	// 0 is for heap memory

        // Shutdown thread pool
        pool.shutdown();

        return make_pair(kmerIdMap, refkmerNumMap);
    }


    // Filtering Sequencing Files
    int build_read_index(
        const vector<string> & inputReadVec, 
        const int & kmerLen,
        const float & matchTd,
        const unordered_map<uint64_t, unordered_map<uint32_t, uint32_t>> & kmerIdMap, 
        const unordered_map<uint32_t, uint32_t> & refkmerNumMap,
        const int & threadsNum,
        const string & prefix, 
        const bool & calRefCovBool
    ) {
        // print log
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Filtering.\n";

        if(inputReadVec.size() == 0) {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Error: please submit sequencing files."
                << endl;
            exit(1);
        }

        for (const auto& inputRead : inputReadVec) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Filtered file: " << inputRead << ".\n";
        }

        read_open(
            inputReadVec, 
            kmerLen, 
            matchTd, 
            kmerIdMap, 
            refkmerNumMap, 
            threadsNum, 
            prefix, 
            calRefCovBool
        );

        return 0;
    }


    // Open the file you want to filter
    int read_open(
        const vector<string> & inputReadVec, 
        const int & kmerLen, 
        const float & matchTd, 
        const unordered_map<uint64_t, unordered_map<uint32_t, uint32_t>> & kmerIdMap, 
        const unordered_map<uint32_t, uint32_t> & refkmerNumMap, 
        const int & threadsNum, 
        const string & prefix, 
        const bool & calRefCovBool
    ) {
        // Number of input files
        const int& FILE_NUM = inputReadVec.size();

        // Input file stream
        vector<gzFile> gzfpIVec;  // vector<inputFileName>

        // Output file stream
        constexpr uint64_t CACHE_SIZE = 1024 * 1024 * 10;
        vector<SAVE::SAVE*> SaveClassVec;  // vector<outputFile>

        for (auto inputFileName : inputReadVec) {
            // Input file stream
            gzfpIVec.push_back(gzopen(inputFileName.c_str(), "rb"));

            // Output file stream
            string outputFile = prefix + "." +  split(inputFileName, "/")[split(inputFileName, "/").size()-1];
            if (outputFile.find(".gz") == string::npos && outputFile.find(".GZ") == string::npos) {
                outputFile += ".gz";
            }
            
            SaveClassVec.push_back(new SAVE::SAVE(outputFile));
        }

        // Process Pool
        ThreadPool pool(threadsNum);
        const int MAX_THREADS_NUM = threadsNum * 100;  // Limit the length of the task queue

        // Initialize the thread pool
        pool.init();

        // Save multi-threaded results to multiple files
        vector<future<vector<string> > > passReadInfoFutureVec;

        // Save the read information that passes the threshold
        vector<stringstream> outStreamVec(FILE_NUM);  // Use stringstream instead of string concatenation   vector<future<vector<string> > >

        // Opening a file
        vector<kseq_t*> ksVec;
        for (size_t i = 0; i < FILE_NUM; i++) {
            if(!gzfpIVec[i]) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'"
                    << inputReadVec[i]
                    << "': No such file or directory." 
                    << endl;
                exit(1);
            }

            ksVec.push_back(kseq_init(gzfpIVec[i]));
        }
        
        while(kseq_read(ksVec[0]) >= 0) {
            for (size_t i = 1; i < FILE_NUM; i++) {
                if (kseq_read(ksVec[i]) < 0) {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "'"
                        << "'"
                        << inputReadVec[i]
                        << "': Incomplete file." 
                        << endl;
                    exit(1);
                }
            }

            vector<string> readNameVec;
            vector<string> readSeqVec;
            vector<string> readInfoVec;

            for (size_t i = 0; i < FILE_NUM; i++) {
                string readName = ksVec[i]->name.s;

                // Determine whether there is junction sequence information
                string readNameLast = ksVec[i]->comment.s ? ' ' + string(ksVec[i]->comment.s) : "";

                // Read name plus label and connector information
                readName += readNameLast;

                // Capitalize
                string readSeq = ksVec[i]->seq.s;
                transform(readSeq.begin(),readSeq.end(),readSeq.begin(),::toupper);

                string readInfo;

                // Determine whether it is fasta or fastq
                if (ksVec[i]->qual.s) {
                    string readQual = ksVec[i]->qual.s;
                    readInfo = "@" + readName + "\n" + readSeq + "\n+\n" + readQual + "\n";
                } else {
                    readInfo = ">" + readName + "\n" + readSeq + "\n";
                }


                readNameVec.push_back(readName);
                readSeqVec.push_back(readSeq);
                readInfoVec.push_back(readInfo);
            }

            // Multithreaded task submission
            passReadInfoFutureVec.push_back(
                pool.submit(
                    read_filter, 
                    ref(kmerLen), 
                    ref(matchTd),
                    ref(kmerIdMap),
                    ref(refkmerNumMap), 
                    readNameVec, 
                    readSeqVec, 
                    readInfoVec, 
                    ref(calRefCovBool)
                )
            );

            // Write once every MAX_THREADS_NUM tasks to prevent loading data into memory all at once
            if (passReadInfoFutureVec.size() >= MAX_THREADS_NUM) {
                // Multithreaded result saving
                for (auto&& passReadInfoFuture : passReadInfoFutureVec) {  // vector<future<vector<string> > >
                    vector<string> readInfoOutVec = move(passReadInfoFuture.get());  // future<vector<string> >

                    // If it is an empty vector, skip it.
                    if (readInfoOutVec.size() == 0) {
                        continue;
                    }

                    for (size_t i = 0; i < FILE_NUM; i++) {  // string
                        outStreamVec[i] << readInfoOutVec[i];
                    }
                }
                passReadInfoFutureVec.clear();
                vector<future<vector<string> > >().swap(passReadInfoFutureVec);
            }

            if (outStreamVec[0].tellp() >= CACHE_SIZE)  // The cache size is MAX_THREADS_NUM
            {
                for (size_t i = 0; i < FILE_NUM; i++)
                {
                    string outTxt = outStreamVec[i].str();
                    SaveClassVec[i]->save(outTxt);
                    // Clear stringstream
                    outStreamVec[i].str(string());
                    outStreamVec[i].clear();
                }
            }
        }

        // Release memory and close the file
        for (size_t i = 0; i < FILE_NUM; i++) {
            kseq_destroy(ksVec[i]);
            gzclose(gzfpIVec[i]);
        }

        // Multithreaded result saving
        for (auto&& passReadInfoFuture : passReadInfoFutureVec) {  // vector<future<vector<string> > >
            vector<string> readInfoOutVec = move(passReadInfoFuture.get());  // future<vector<string> >

            // If it is an empty vector, skip it.
            if (readInfoOutVec.size() == 0) {
                continue;
            }

            for (size_t i = 0; i < FILE_NUM; i++) {
                outStreamVec[i] << readInfoOutVec[i];
            }
            
        }
        passReadInfoFutureVec.clear();
        vector<future<vector<string> > >().swap(passReadInfoFutureVec);

        if (outStreamVec[0].tellp() >= 0) {  // Write it last time
            for (size_t i = 0; i < FILE_NUM; i++) {
                string outTxt = outStreamVec[i].str();
                SaveClassVec[i]->save(outTxt);
                // clear stringstream
                outStreamVec[i].str(string());
                outStreamVec[i].clear();

                // Deleting a pointer
                if (SaveClassVec[i] != nullptr) {
                    delete SaveClassVec[i];
                    SaveClassVec[i] = nullptr;
                }
            }
            vector<SAVE::SAVE*>().swap(SaveClassVec);
        }

        malloc_trim(0);	// 0 is for heap memory

        // Shutdown thread pool
        pool.shutdown();
        
        return 0;
    }

    // Filtering Files
    vector<string> read_filter(
        const int & kmerLen,
        const float & matchTd, 
        const unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> > & kmerIdMap, 
        const unordered_map<uint32_t, uint32_t> & refkmerNumMap, 
        vector<string> readNameVec, 
        vector<string> readSeqVec, 
        vector<string> readInfoVec, 
        const bool & calRefCovBool
    ) {
        for (size_t i = 0; i < readSeqVec.size(); i++) {
            string readName = readNameVec[i];
            string readSeq = readSeqVec[i];
        
            // Record the kmer of the read  map<kmerScore, map<readId, kmerNum> >
            unordered_map<uint64_t, unordered_map<uint32_t, uint32_t> > kmerNameMap = kmerBit::kmer_sketch1(0, readSeq, kmerLen);

            // The number of kmers in the read
            uint32_t readKmerNum = 0;

            // Record the number of matches for read
            // <readId, matchNum>
            unordered_map<uint32_t, uint32_t> kmerMatchNumMap;
            kmerMatchNumMap.reserve(kmerNameMap.size());  // Preallocate memory

            // Record the number of kmers that match
            uint32_t kmerMatchNum = 0;

            // Check whether the read exceeds the threshold
            for (const auto& [key1, value1] : kmerNameMap) {  // map<kmerScore, map<readId, kmerNum> >
                const uint32_t kmerNum = value1.at(0);

                // Check if the kmer of the read exists in the total kmer hash table
                auto findIter1 = kmerIdMap.find(key1);  //  map<kmerScore, map<readId, kmerNum> >

                if (findIter1 != kmerIdMap.end()) {
                    for (const auto& [key2, value2] : findIter1->second) {  //  map<readId, kmerNum>
                        kmerMatchNumMap[key2] += value2;
                    }

                    // The number of matching kmers is superimposed
                    kmerMatchNum += kmerNum;
                }

                // Overlay of kmer counts of reads
                readKmerNum += kmerNum; 
            }

            // Before judging whether it is greater than the threshold, first check whether the read is a repeated sequence. If so, skip the read
            if (kmerNameMap.size()/(float)readKmerNum < 0.3f) {
                // Return after clearing
                readInfoVec.clear();
                vector<string>().swap(readInfoVec);

                return readInfoVec;
            }

            // Determine whether it is greater than the threshold
            for (auto& [readId, matchNum] : kmerMatchNumMap) {  // <readId, matchNum>
                // First look at the number of kmers corresponding to the reference genome
                uint32_t refkmerNum = 0;
                auto findIter1 = refkmerNumMap.find(readId);  // <readId, kmerNum (number)>
                if (findIter1 != refkmerNumMap.end()) {  // If you have
                    refkmerNum = findIter1->second;
                } else  {  // If there is no read, it means there is a problem with the code and an error is reported.
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "keyError: " << readId << " not in refkmerNumMap.\n";
                    exit(1);
                }
                
                // Determine the proportion of a sequence in the reference genome that is covered. If the proportion of the reference genome is not calculated, the value is assigned to 0
                float refRatio = calRefCovBool ? matchNum/(float)refkmerNum : 0.0f;
                // Determine the proportion of sequencing reads covered
                float readRatio = kmerMatchNum/(float)readKmerNum;
                
                // Judgment threshold
                if (refRatio >= matchTd || readRatio >= matchTd) {
                    // Print progress and return results
                    cerr << "[" << __func__ << "::" << getTime() << "] " 
                        << "pass: " 
                        << readName 
                        << " refgenome:" << refRatio 
                        << " read:" << readRatio 
                        << endl;

                    return readInfoVec;
                }
            }
        }

        // Return after clearing
        readInfoVec.clear();
        vector<string>().swap(readInfoVec);

        return readInfoVec;
    }
}

#endif