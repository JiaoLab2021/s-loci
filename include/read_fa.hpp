#ifndef filter_hpp
#define filter_hpp
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
#include "kmer.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include "ThreadPool.hpp"
#include "bitmap.hpp"
#include "bloomfilter.hpp"

using namespace std;

// kseq.h Opening a file
KSEQ_INIT(gzFile, gzread)

namespace filter{
    vector<vector<string>> fasta_open(string inputFileName, int kmerLen, int minimizerK, int minimizerW, int threads_num) {
        vector<vector<string>> allMinimizerVectorBase;
        string fastaSeq{};

        // Reading files in a loop
        string information;
        ifstream inputFile;
        inputFile.open(inputFileName, ios::in);

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Collect kmer and minimizer.\n";

        // File open failed
        if(!inputFile.is_open()) {
            cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << inputFileName 
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        } else {
            while(getline(inputFile,information)) {
                if(information.empty()) {
                    continue;
                } else {
                    // '>' is an identifier
					string::size_type idx = information.find(">");

                    // If '>' is encountered and fastaSeq is not empty
                    if (idx != string::npos) {
                        if (fastaSeq.length() > 0) {
                            vector<string> kmerVectorBase = build_kmer_index(fastaSeq, kmerLen);
                            int kmerVectorBaseSize = kmerVectorBase.size();

                            vector<string> minimizerVectorBase;
                            // Enable multithreading
                            thread threads[threads_num];
                            for (int i = 0; i < threads_num; i++) {
                                int indexLeft = (kmerVectorBaseSize/threads_num) * i;
                                int indexRight = (kmerVectorBaseSize/threads_num) * (i + 1);
                                if (i == (threads_num - 1)) {  // The last index of the thread is the length of the vector
                                    indexRight = kmerVectorBaseSize;
                                }
                                
                                threads[i] = thread(
                                    build_minimizer_index,
                                    ref(kmerVectorBase), 
                                    ref(minimizerVectorBase), 
                                    minimizerK, 
                                    minimizerW, 
                                    indexLeft, 
                                    indexRight
                                );
                            }
                            // Thread blocking
                            for (auto& t: threads) {
                                t.join();
                            }

                            allMinimizerVectorBase.push_back(minimizerVectorBase);

                            // Clear the string and prepare for the next loop
                            fastaSeq.clear();
                        }
                    } else {
                        transform(information.begin(),information.end(),information.begin(),::toupper);
                        fastaSeq += information;
                    }
                }
            }
            // If the end of the file is reached, an index is constructed for the string.
            vector<string> kmerVectorBase = build_kmer_index(fastaSeq, kmerLen);
            int kmerVectorBaseSize = kmerVectorBase.size();

            vector<string> minimizerVectorBase;
            // Enable multithreading
            thread threads[threads_num];
            for (int i = 0; i < threads_num; i++) {
                int indexLeft = (kmerVectorBaseSize/threads_num) * i;
                int indexRight = min((kmerVectorBaseSize/threads_num) * (i + 1), kmerVectorBaseSize);
                if (i == threads_num - 1) {  // The last index of the thread is the length of the vector
                    indexRight = kmerVectorBaseSize;
                }

                threads[i] = thread(
                    build_minimizer_index,
                    ref(kmerVectorBase), 
                    ref(minimizerVectorBase), 
                    minimizerK, 
                    minimizerW, 
                    indexLeft, 
                    indexRight
                );
            }
            // Thread blocking
            for (auto& t: threads) {
                t.join();
            }

            allMinimizerVectorBase.push_back(minimizerVectorBase);
            
            // Clear a string
            fastaSeq.clear();
        }
        // Close File
        inputFile.close();

        // Vector deduplication
        for (int i = 0; i < allMinimizerVectorBase.size(); i++) {
            sort(allMinimizerVectorBase[i].begin(), allMinimizerVectorBase[i].end());
            allMinimizerVectorBase[i].erase(
                unique(
                    allMinimizerVectorBase[i].begin(), 
                    allMinimizerVectorBase[i].end()
                ), 
                allMinimizerVectorBase[i].end()
            );
        }
        return allMinimizerVectorBase;
    }
}

#endif