#ifndef genotype_hpp
#define genotype_hpp
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <thread>
#include "zlib.h"
#include <regex>
#include <getopt.h>
#include <tuple>
#include "kmer_bit_genotype.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "ThreadPool.hpp"

using namespace std;

// kseq.h open fiile
// KSEQ_INIT(gzFile, gzread)

namespace genotype{
    // Open fastq/a.gz file
    void read_open(
        const vector<string>& inputReadVec, 
        const int& kmerLen, 
        const unordered_map<uint64_t, vector<uint64_t>> & kmerBaseMap, 
        unordered_map<uint64_t, map<uint32_t, uint32_t>> & HapkmerMatMap, 
        unordered_map<uint64_t, uint64_t>& readAllKmerMap, 
        const int& threadsNum
    );

    // read multithreaded function
    tuple<unordered_map<uint64_t, map<uint32_t, uint32_t> >, unordered_map<uint64_t, uint64_t> > read_kmer(
        const int & kmerLen, 
        vector<string> readSeqVec, 
        const unordered_map<uint64_t, vector<uint64_t>> & kmerBaseMap
    );

    // Calculate site_cov
    string site_cov(
        int kmerLen,  
        const unordered_map<uint64_t, tuple<string, uint32_t>> & chromosomeIdLenMap, 
        const unordered_map<uint64_t, map<uint32_t, uint32_t>> & HapkmerMatMap
    );

    /** Reading fasta files 
     * @param inputFileName   
     * @param kmerLen         
     * @param kmerBaseMap     minimizer output hash table
     *          unordered_map<uint64_t, vector<uint64_t>>
     *                        (uint64_t) -> kMer<<8 | kmerSpan
     *                        (uint64_t) -> chromosomeId<<32 | lastPos<<1 | strand
    **/
    tuple<unordered_map<uint64_t, vector<uint64_t> >, unordered_map<uint64_t, tuple<string, uint32_t> > > build_base_index(
        const string & inputFileName, 
        const int & kmerLen
    ) {
        // print log
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Indexed file: " << inputFileName << ".\n";

        // Record all kmer results
        unordered_map<uint64_t, vector<uint64_t> > kmerBaseMap;  // unordered_mapmap<hash, vector<chr/end/strand>>

        // Initialize chromosome id
        uint64_t chromosomeId = 0;
        unordered_map<uint64_t, tuple<string, uint32_t> > chromosomeIdLenMap;  // map<ID, tuple<chr, length> >

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
                // ks->seq.l The sequence length is recorded
                string chromosome = ks->name.s;
                
                // Capitalize
                string sequence = ks->seq.s;
                uint64_t length = ks->seq.l;
                transform(sequence.begin(),sequence.end(),sequence.begin(),::toupper);

                // Building kmer index
                kmerBitGenotype::kmer_sketch(chromosomeId, sequence, kmerLen, kmerBaseMap);

                // Construct chromosome number
                chromosomeIdLenMap[chromosomeId] = make_pair(chromosome, length);

                // Update chromosome ID
                chromosomeId++;

                // Clear the string and release the memory
                chromosome.clear();
                sequence.clear();
                string().swap(chromosome);
                string().swap(sequence);
                
            }

            // Release memory and close the file
            kseq_destroy(ks);
            gzclose(gzfp);
        }
        
        return make_pair(kmerBaseMap, chromosomeIdLenMap);
    }

    // Genotype
    tuple<unordered_map<uint64_t, map<uint32_t, uint32_t> >, unordered_map<uint64_t, uint64_t> > genotype(
        vector<string> inputReadVec,  
        int kmerLen, 
        const unordered_map<uint64_t, vector<uint64_t>> & kmerBaseMap, 
        const unordered_map<uint64_t, tuple<string, uint32_t>> & chromosomeIdLenMap, 
        int threadsNum
    ) {
        // print log
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Genotyping.\n";

        for (const auto& inputRead : inputReadVec) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Sequencing file: " << inputRead << ".\n";
        }

        // Record the chromosome number of the haplotype corresponding to the kmer matched by the sequencing file, and the number of times END corresponds to
        unordered_map<uint64_t, map<uint32_t, uint32_t> > HapkmerMatMap; // unordered_map<chrId, map<end, frequency> >
        unordered_map<uint64_t, uint64_t> readAllKmerMap;  // map<kmerHash, frequence>

        // readBuild index
        read_open(
            inputReadVec,
            kmerLen, 
            kmerBaseMap, 
            HapkmerMatMap, 
            readAllKmerMap, 
            threadsNum
        );
        
        return make_tuple(HapkmerMatMap, readAllKmerMap);
    }


    // open file
    void read_open(
        const vector<string>& inputReadVec, 
        const int& kmerLen, 
        const unordered_map<uint64_t, vector<uint64_t>> & kmerBaseMap, 
        unordered_map<uint64_t, map<uint32_t, uint32_t>> & HapkmerMatMap, 
        unordered_map<uint64_t, uint64_t>& readAllKmerMap, 
        const int& threadsNum
    ) {
        // readBuild index
        const int& FILE_NUM = inputReadVec.size();

        // Input file stream
        vector<gzFile> gzfpIVec;  // vector<inputFileName>
        for (auto inputFileName : inputReadVec) {
            // Input file stream
            gzfpIVec.push_back(gzopen(inputFileName.c_str(), "rb"));
        }

        // Process Pool
        ThreadPool pool(threadsNum);
        const int MAX_THREADS_NUM = threadsNum * 100;  // Limit the length of the task queue

        // Initialize the thread pool
        pool.init();

        // Saving multi-threaded results
        // unordered_map<uint64_t, map<uint32_t, uint32_t> > kmerMatMapTmp;  // map<chromosomeId, map<refEnd, number> >
        // unordered_map<uint64_t, uint64_t> readAllKmerMap;  // map<kmerHash, frequence>
        vector<future<tuple<unordered_map<uint64_t, map<uint32_t, uint32_t> >, unordered_map<uint64_t, uint64_t> > > > kmerMatMapFutureVec;  // vector<future<tuple<map<chromosomeId, map<refEnd, number> >, map<kmerHash, frequence> > > >

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

            vector<string> readSeqVec;

            for (size_t i = 0; i < FILE_NUM; i++) {
                // Capitalize
                string readSeq = ksVec[i]->seq.s;
                transform(readSeq.begin(),readSeq.end(),readSeq.begin(),::toupper);

                readSeqVec.push_back(readSeq);
            }

            // Check if the task queue exceeds the threshold, and wait if it exceeds the threshold to prevent loading the data into memory all at once
            int maxRetries = 1200; // 10 minutes
            int retryCount = 0;
            while (pool.get_queue() >= MAX_THREADS_NUM) {
                if (retryCount >= maxRetries) {
                    break;
                }
                // Check every 0.5 seconds
                sleep(0.5);
                retryCount++;
            }

            kmerMatMapFutureVec.push_back(
                pool.submit(
                    read_kmer, 
                    kmerLen, 
                    readSeqVec,
                    ref(kmerBaseMap)
                )
            );

            // Write once every MAX_THREADS_NUM tasks to prevent loading data into memory all at once
            if (kmerMatMapFutureVec.size() >= MAX_THREADS_NUM) {
                // Multithreaded result saving
                for (auto&& kmerMatMapFuture : kmerMatMapFutureVec) {  // vector<future<map<chromosomeId, map<refEnd, number> > > >
                    // Results returned by multiple threads
                    unordered_map<uint64_t, map<uint32_t, uint32_t> > kmerMatMapTmp;  // future<map<chromosomeId, map<refEnd, number> > >
                    unordered_map<uint64_t, uint64_t> readAllKmerMapTmp;  // map<kmerHash, frequence>
                    tie(kmerMatMapTmp, readAllKmerMapTmp) = move(kmerMatMapFuture.get());

                    // Add the result to the total kmerMatMap
                    for (auto iter1 = kmerMatMapTmp.begin(); iter1 != kmerMatMapTmp.end(); ++iter1) {
                        for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2) {
                            HapkmerMatMap[iter1->first][iter2->first] += iter2->second;
                        }
                    }
                    for (auto iter1 = readAllKmerMapTmp.begin(); iter1 != readAllKmerMapTmp.end(); ++iter1) {
                        readAllKmerMap[iter1->first] += iter1->second;
                    }
                }
                kmerMatMapFutureVec.clear();
                vector<future<tuple<unordered_map<uint64_t, map<uint32_t, uint32_t> >, unordered_map<uint64_t, uint64_t> > > >().swap(kmerMatMapFutureVec);

                malloc_trim(0);	// 0 is for heap memory
            }
        }

        // Release memory and close the file
        for (size_t i = 0; i < ksVec.size(); i++) {
            kseq_destroy(ksVec[i]);
            gzclose(gzfpIVec[i]);
        }

        // Last saved
        if (kmerMatMapFutureVec.size() > 0) {
            // Multithreaded result saving
            for (auto&& kmerMatMapFuture : kmerMatMapFutureVec) {  // vector<future<map<chromosomeId, map<refEnd, number> > > >
                // Results returned by multiple threads
                unordered_map<uint64_t, map<uint32_t, uint32_t> > kmerMatMapTmp;  // future<map<chromosomeId, map<refEnd, number> > >
                unordered_map<uint64_t, uint64_t> readAllKmerMapTmp;  // map<kmerHash, frequence>
                tie(kmerMatMapTmp, readAllKmerMapTmp) = move(kmerMatMapFuture.get());

                // Add the result to the total kmerMatMap
                for (auto iter1 = kmerMatMapTmp.begin(); iter1 != kmerMatMapTmp.end(); ++iter1) {
                    for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2) {
                        HapkmerMatMap[iter1->first][iter2->first] += iter2->second;
                    }
                }
                for (auto iter1 = readAllKmerMapTmp.begin(); iter1 != readAllKmerMapTmp.end(); ++iter1) {
                    readAllKmerMap[iter1->first] += iter1->second;
                }
            }
            kmerMatMapFutureVec.clear();
            vector<future<tuple<unordered_map<uint64_t, map<uint32_t, uint32_t> >, unordered_map<uint64_t, uint64_t> > > >().swap(kmerMatMapFutureVec);

            malloc_trim(0);	// 0 is for heap memory
        }

        // Shutdown thread pool
        pool.shutdown();
    }

    // fastq/a.gz multithreaded function
    tuple<unordered_map<uint64_t, map<uint32_t, uint32_t> >, unordered_map<uint64_t, uint64_t> > read_kmer(
        const int& kmerLen, 
        vector<string> readSeqVec, 
        const unordered_map<uint64_t, vector<uint64_t> >& kmerBaseMap
    ) {
        // Temporary map, save the matched reference chromosome number
        unordered_map<uint64_t, map<uint32_t, uint32_t> > kmerMatMapTmp;  // map<chromosomeId, map<refEnd, number> >
        unordered_map<uint64_t, uint64_t> readAllKmerMapTmp;  // map<kmerHash, frequence>

        // Looping over a list of sequences
        for (size_t i = 0; i < readSeqVec.size(); i++) {
            string readSeq = readSeqVec[i];

            // Building kmer index
            unordered_map<uint64_t, vector<uint64_t> > readKmerMap;  // map<kmerHash, vec<chrId<<32 | end | strand>>
            kmerBitGenotype::kmer_sketch(0, readSeq, kmerLen, readKmerMap);
            
            // Count how many kmers in the read are in kmerBaseMap
            for(auto iter1 = readKmerMap.begin(); iter1 != readKmerMap.end(); ++iter1) {
                auto find1Iter = kmerBaseMap.find(iter1->first);
                if (find1Iter != kmerBaseMap.end()) {
                    // Record the kmer hash value of the read
                    readAllKmerMapTmp[iter1->first]++;

                    // END of record matching
                    for(auto iter2 = find1Iter->second.begin(); iter2 != find1Iter->second.end(); ++iter2) {
                        uint64_t chromosomeId = *iter2>>32;
                        uint32_t refEnd = ((*iter2)-(chromosomeId<<32))>>1;
                        kmerMatMapTmp[chromosomeId][refEnd]++;
                    }
                }
            }
        }

        return make_tuple(kmerMatMapTmp, readAllKmerMapTmp);
    }


    // Calculate Jaccard_score
    string Jaccard_score(
        const int& kmerLen, 
        const unordered_map<uint64_t, tuple<string, uint32_t>> & chromosomeIdLenMap,
        const unordered_map<uint64_t, map<uint32_t, uint32_t>> & HapkmerMatMap
    ) {
        // Calculate Jaccard_score
        string outTxt = "#kmer-" + to_string(kmerLen) + "\n#genotype\tJaccard_score\n";
        for (const auto& it : HapkmerMatMap) {
            auto findIter = chromosomeIdLenMap.find(it.first);
            if (findIter != chromosomeIdLenMap.end()) {
                // Divide the number of matching kmers by the total number of kmers in the haplotype
                float score = (float)it.second.size() / ((float)get<1>(findIter->second)-kmerLen+1);

                outTxt += get<0>(findIter->second) + "\t" + to_string(score) + "\n";
            }
        }
        
        return outTxt;
    }


    // tuple<score, chromosome1, chromosome2>
    int cal_pro(
        unordered_map<uint64_t, tuple<string, uint32_t> > chromosomeIdLenMap, 
        const string& jaccardTxt, 
        vector<string> outHapVec, 
        const string& outputFilePrefix
    ) {
        // Preserve genotype coverage
        map<string, float> covMap;  // map<chr, Jaccard_score>

        // Comment lines for jaccardTxt
        string headTxt;

        // Haplotype index
        for (auto iter : chromosomeIdLenMap) {  // map<ID, tuple<chr, length> >
            string chromosome = get<0>(iter.second);
            covMap[chromosome] = 0.0;
        }

        // Split jaccardTxt
        vector<string> jaccardVec = split(jaccardTxt, "\n");
        for (int i = 0; i < jaccardVec.size(); i++) {
            // Skip Comment Lines
            if (jaccardVec[i].find("#") != string::npos) {
                headTxt += jaccardVec[i] + "\n";
            } else {
                vector<string> jaccardVecVec = split(strip(jaccardVec[i], '\n'), "\t");  //  genotype   Jaccard_score
                if (jaccardVecVec.size() < 2) {  // Skip empty lines
                    continue;
                }
                string chromosome = jaccardVecVec[0];
                float jaccard = stof(jaccardVecVec[1]);
                covMap[chromosome] = jaccard;
            }
        }
        
        // Initial probability matrix
        map<string, float> InitialMatrix;  // map<chr, 1/Jaccard_score>
        for (const auto& iter : covMap) {  // map<chr, Jaccard_score>
            InitialMatrix[iter.first] = 1/(float)covMap.size();
        }

        // correctMatrix, homozygous and heterozygous probabilities of SRNase haplotypes
        float hoRatio = 0.02;
        float heRatio = 1 - hoRatio;
        map<string, float> correctMatrix = {
            {"ho", hoRatio}, 
            {"he", heRatio}
        };

        // Transition probability matrix
        map<string, map<string, float>> transitionMatrix;  // map<chr1, map<chr2, pro>>
        for (auto iter1 : InitialMatrix) {  // map<chr, 1/Jaccard_score>
            string genotype1 = iter1.first;  // Name of state 1

            for (auto iter2 : InitialMatrix) {  // map<chr, 1/Jaccard_score>
                string genotype2 = iter2.first;  // State 2 name

                if (genotype1 == genotype2) {  // Homozygous conversion
                    transitionMatrix[genotype1][genotype2] = correctMatrix["ho"];
                } else  {  // Heterozygous conversion
                    transitionMatrix[genotype1][genotype2] = correctMatrix["he"];
                }
            }
        }

        // Output matrix, 0->true 1->false
        map<string, unordered_map<int, float> > outMatrix;
        for (auto iter : covMap)  // map<chr, Jaccard_score>
        {
            outMatrix[iter.first][0] = iter.second;
            outMatrix[iter.first][1] = 1 - iter.second;
        }

        // Output result matrix
        double outputMatrix[covMap.size()][covMap.size()];
        int iIndex = 0;
        int jIndex = 0;
        map<int, string> genotypeIndexMap;
        for (auto iter1 : covMap) {  // map<chr, Jaccard_score>
            string genotype1 = iter1.first;  // Name of state 1
            for (auto iter2 : covMap) {  // map<chr, Jaccard_score>
                string genotype2 = iter2.first;  // State 2 name
                outputMatrix[iIndex][jIndex] = transitionMatrix[genotype1][genotype2] * 
                                               outMatrix[genotype1][0] * 
                                               outMatrix[genotype2][0];
                jIndex++;
            }
            jIndex = 0;
            genotypeIndexMap[iIndex] = genotype1; // Recording the order of chromosomes
            iIndex++;
        }
        iIndex = 0;
        jIndex = 0;
        
        // The maximum probability of backtracking
        float maxIndexRow = -1.1;
        float maxIndexCol = -1.1;
        float maxScore = 0.0;

        for (size_t i = 0; i < covMap.size(); i++) {
            for (size_t j = 0; j < covMap.size(); j++) {
                if (outputMatrix[i][j] > maxScore) {
                    maxScore = outputMatrix[i][j];
                    maxIndexRow = i;
                    maxIndexCol = j;
                }
            }
        }

        // Output
        ofstream outputFileJaccard;
        outputFileJaccard.open(outputFilePrefix + ".hap", ios::out);
        outputFileJaccard << "#Haplotype-Jaccard 1/2: " + 
                        genotypeIndexMap[maxIndexRow] + "/" +
                        genotypeIndexMap[maxIndexCol] + 
                        "\n" + 
                        "#Haplotype-hmm 1/2: " + join(outHapVec, "/") + 
                        "\n";
        outputFileJaccard << headTxt;
        // Transform the coverage map to sort in ascending order
        map<float, vector<string>> covConvertMap;
        for (auto iter : covMap) {
            covConvertMap[iter.second].push_back(iter.first);
        }
        for (auto iter1 : covConvertMap) {
            for (auto iter2 : iter1.second) {
                outputFileJaccard << iter2 << "\t" << to_string(iter1.first) << endl;
            }
        }
        
        // Close File
        outputFileJaccard.close();

        return 0;
    }


    // Calculate site_cov
    int site_cov(
        const int& kmerLen, 
        const unordered_map<uint64_t, tuple<string, uint32_t> >& chromosomeIdLenMap, 
        const unordered_map<uint64_t, map<uint32_t, uint32_t> >& HapkmerMatMap, 
        const string& outputFilePrefix
    ) {
        // Calculate site coverage
        map<string, map<uint32_t, uint32_t> > sitrCovMap; // map<chr, map<end, cov>>
        for (const auto& it1 : chromosomeIdLenMap) {  // map<ID, tuple<chr, length> >
            // Chromosome number, end and chromosome length information corresponding to kmer
            uint64_t chromosomeId = it1.first;
            string chromosome = get<0>(it1.second);
            uint32_t length = get<1>(it1.second);

            // Check if there is a haplotype match, if not, skip
            auto findIter1 = HapkmerMatMap.find(chromosomeId);  // unordered_map<chrId, map<end, frequency> >
            if (findIter1 != HapkmerMatMap.end()) {
                map<uint32_t, uint32_t> endFeqMap = findIter1->second; // map<end, frequency>
                for (const auto& it2 : endFeqMap) {
                    uint32_t refEnd = it2.first;
                    uint32_t frequency = it2.second;

                    // Go forward kmer bases from the end position
                    for (uint32_t i = refEnd-kmerLen+1; i < refEnd+1; i++) {
                        sitrCovMap[chromosome][i] += frequency;
                    }
                }
            }

            // Fill the uncovered areas with 0
            for (uint32_t i = 1; i < length+1; i++) {
                auto findIter1 = sitrCovMap[chromosome].find(i);
                if (findIter1 == sitrCovMap[chromosome].end()) {
                    sitrCovMap[chromosome][i] = 0;
                }
            }
        }

        // Output coverage information
        string outTxt = "#kmer-" + to_string(kmerLen) + "\n#genotype\tlocation\tcov\n";
        for (auto it1 : sitrCovMap) {
            for (auto it2 : it1.second) {
                outTxt += it1.first + "\t" + to_string(it2.first) + " " + to_string(it2.second) + "\n";
            }
        }

        // Output
        ofstream outputFileCov;
        outputFileCov.open(outputFilePrefix + ".cov", ios::out);
        outputFileCov << outTxt;
        // Close File
        outputFileCov.close();
        
        return 0;
    }
}

#endif