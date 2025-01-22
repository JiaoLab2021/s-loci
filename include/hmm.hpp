#ifndef HMM_ME_hpp
#define HMM_ME_hpp
#include <numeric>
#include <iostream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <vector>
#include <future>
#include "tuple"
#include "include/get_time.hpp"

using namespace std;


class HMM {
private:
    bool debug_ = false;

    // Haplotype number
    uint32_t hapNum_;
    unordered_map<uint64_t, uint64_t> readAllKmerMap_;  // unordered_mapmap<hash, frequency>, read matches all haplotype kmers
    unordered_map<uint64_t, vector<uint64_t> > kmerBaseMap_;  // All kmer information of haplotype unordered_map<kmerHash, vector<chr/end/strand>>
    unordered_map<uint64_t, tuple<string, uint32_t> > chromosomeIdLenMap_;  // Chromosome number and length information of all haplotypes map<chrId, tuple<chr, length> >

    // Initial Probability
    double initialProbabilities_;

    // Observations
    vector<uint64_t> observationsVec_;

    // Transition probability
    long double recombProb_;
    long double noRecombProb_;

    // Emission probability
    unordered_map<uint64_t, unordered_map<uint64_t, uint16_t> > emissionMatrix_;  // map<kmerHash, map<chrId, 0/1> >

    // Score calculated by the forward algorithm
    unordered_map<uint64_t, unordered_map<uint64_t, long double> > forwardMatrix_; // map<kmerHash, map<chrId, score> >

    // The frequency of haplotypes on all kmers
    map<uint32_t, vector<uint64_t> > freHapVecMap_;  // map<frequence, vector<hapID> >

    int initial_states();  // Calculate initial probability
    int observable_states();  // Calculate initial probability
    int transition_probabilities();  // Calculating transition probabilities
    int emission_states();  // Calculating the emission probability
public:
    /*
    * @author zezhen du
    * @date 2023/4/20
    * @version v1.0
    * @brief init
    * 
    * @param readAllKmerMap           map<kmerHash, frequence>
    * @param kmerBaseMap              unordered_map<hash, vector<chr/end/strand>>
    * @param chromosomeIdLenMap       map<ID, tuple<chr, length> >
    * @param noRecombProb             Non-recombination rate
    * @param debug                    Whether to debug the code
    */
    HMM(
        const unordered_map<uint64_t, uint64_t>& readAllKmerMap, 
        const unordered_map<uint64_t, vector<uint64_t> >& kmerBaseMap, 
        const unordered_map<uint64_t, tuple<string, uint32_t> >& chromosomeIdLenMap, 
        long double noRecombProb = 0.99, 
        bool debug = false
    ) {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Hidden Markov Model.\n";

        readAllKmerMap_ = readAllKmerMap;
        kmerBaseMap_ = kmerBaseMap;
        chromosomeIdLenMap_ = chromosomeIdLenMap;
        hapNum_ = chromosomeIdLenMap_.size();
        noRecombProb_ = noRecombProb;
        debug_ = debug;

        // Calculate initial probability
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculating initial probability.\n";
        initial_states();
        // Calculate observations
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculating observation matrix.\n";
        observable_states();
        // Calculating transition probabilities
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculating transition probability matrix.\n";
        transition_probabilities();
        // Calculating the emission probability
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculating emission probability matrix.\n";
        emission_states();
    };
    
    int forward();  // Forward Algorithm
    int get_path();  // Get Path
    vector<string> get_result(const int& needHapNum_);  // Get the specified number of haplotypes
};


/**
 * @author zezhen du
 * @date 2023/4/20
 * @version v1.0
 * @brief initial_states
 * 
 * @param nodeDistence   hapNum haplotype number
 * 
 * @return 0
**/
int HMM::initial_states() {
    initialProbabilities_ = 1 / (double) hapNum_;
    
    return 0;
}


/**
 * @author zezhen du
 * @date 2023/4/20
 * @version v1.0
 * @brief Observable states
 * 
 * @param readAllKmerMap       unordered_mapmap<hash, frequency>, read matches all haplotype kmers
 * 
 * @return 0
**/
int HMM::observable_states() {
    for (const auto& iter1 : readAllKmerMap_) {  // unordered_mapmap<hash, frequence>
        observationsVec_.push_back(iter1.first);
    }
    
    return 0;
}



/**
 * @author zezhen du
 * @date 2023/4/20
 * @version v1.0
 * @brief Transition probabilities
 * 
 * @param nodeDistence   hapNum haplotype number
 * @param noRecombProb   Probability of no recombination
 * 
 * @return 0
**/
int HMM::transition_probabilities() {
    // Transition probability
    recombProb_ = (1.0L - noRecombProb_) * (1.0L / ((long double) hapNum_ - 1));  // (1 - noRecombProb)/(1/(hapNum-1))

    return 0;
}



/**
 * @author zezhen du
 * @date 2023/4/20
 * @version v1.0
 * @brief emission_states
 * 
 * @param kmerBaseMap         All kmer information of haplotype unordered_map<kmerHash, vector<chr/end/strand>>
 * @param chromosomeIdLenMap  Chromosome number and length information for all haplotypes   map<chrId, tuple<chr, length> >
 * @param observationsVec     All matched kmer vectors of haplotypes <kmerHash>, observation matrix
 * 
 * @return 0
**/
int HMM::emission_states() {
    // Calculate the frequency of each haplotype on this kmer
    for (const auto& kmerHash : observationsVec_) {  // vector<kmerHash>
        // On which haplotypes does this kmer contain
        unordered_map<uint64_t, uint8_t> hapKmerMap;  // map<chrID, 1>
        // Find out whether the kmer has appeared in the total kmer hash table of haplotypes
        auto findIter1 = kmerBaseMap_.find(kmerHash);
        if (findIter1 != kmerBaseMap_.end()) {
            // Variable Binding
            const auto& kmerBaseVec = kmerBaseMap_.at(kmerHash);
            for (const auto& hapInfo : kmerBaseVec) {  // vector<chr/end/strand>
                uint64_t chromosomeId = hapInfo>>32;
                hapKmerMap[chromosomeId] = 1;
            }
        }
        
        // Variable Binding
        auto& emissionMatrixTmp = emissionMatrix_[kmerHash];

        for (const auto& [chromosomeId, chrInfoTup] : chromosomeIdLenMap_) {  // map<chrId, tuple<chr, length> >
            // Whether the kmer haplotype matches 1 matches 0 does not match
            uint32_t HapkmerMatchNum = 0;

            // Check whether the haplotype matches the kmer in hapKmerMap
            const auto& findIter2 = hapKmerMap.find(chromosomeId);
            if (findIter2 != hapKmerMap.end()) {
                HapkmerMatchNum = 1;
            }

            // Frequency assignment of this haplotype under this kmer
            emissionMatrixTmp[chromosomeId] = HapkmerMatchNum;
        }
    }

    // Print the emission probability matrix
    if (debug_) {
        for (const auto& [kmerHash, hapInfoMap] : emissionMatrix_) {  // map<kmerHash, map<chrId, 0/1> >
            if (kmerHash == emissionMatrix_.begin()->first) {
                cerr << "\t\t";
                for (const auto& [chrID, kmerState] : hapInfoMap) {
                    cerr << "\t" << chrID;
                }
                cerr << endl;
            }

            cerr << kmerHash << "\t";
            for (const auto& [chrID, kmerState] : hapInfoMap) {
                cerr << kmerState << "\t";
            }
            cerr << endl;
        }
    }
    
    
    return 0;
}


/**
 * @author zezhen du
 * @date 2023/4/20
 * @version v1.0
 * @brief forward
 * 
 * @param initialProbabilities_  initialProbabilities_     Initial probability
 * @param observationsVec_       All kmers matched to the haplotype vector<hash>, the observation matrix
 * @param recombProbTup          make_tuple(recombProb, noRecombProb)
 * @param emissionMatrix         map<kmerHash, map<chrId, frequence> >
 * 
 * @return 0
**/
int HMM::forward() {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Calculating probability using the forward algorithm.\n";

    // Observation Matrix
    const auto& firstKmerHash = observationsVec_[0];  // The first kmer
    uint64_t preKmerHash= 0.0;
    for (const auto& kmerHash : observationsVec_) {  // vector<hash>
        const auto& emissionMatrixTmp = emissionMatrix_.at(kmerHash);  // Emission probability matrix  map<chrId, frequence>

        auto& forwardMatrixTmp = forwardMatrix_[kmerHash];
        
        for (const auto& [chromosomeId, frequence] : emissionMatrixTmp) {
            // The frequency is 1, the sequencing is correct
            double sequenceStatus = frequence == 1 ? 0.99 : 0.01;
            

            // If it is the first node
            if (kmerHash == firstKmerHash) {
                forwardMatrixTmp[chromosomeId] = initialProbabilities_ * sequenceStatus;
                if (debug_) {
                    cerr << kmerHash << " initialProbabilities:" << initialProbabilities_ << " frequence:" << frequence << endl;
                }
            } else {
                const auto& forwardMatrixTmpPre = forwardMatrix_[preKmerHash];  // The result of the forward algorithm at the previous node
                const auto& emissionMatrixTmpPre = emissionMatrix_.at(preKmerHash);  // The emission probability of the previous node  map<chrId, frequence>

                for (const auto& [chromosomeIdPre, frequencePre] : emissionMatrixTmpPre) {
                    if (chromosomeId == chromosomeIdPre) {
                        forwardMatrixTmp[chromosomeId] += forwardMatrixTmpPre.at(chromosomeIdPre) * noRecombProb_ * sequenceStatus;
                        if (debug_) {
                            cerr << kmerHash << " forwardPre:" << forwardMatrixTmpPre.at(chromosomeIdPre) << " np-recombProb:" << noRecombProb_ << " sequenceStatus:" << sequenceStatus << " forward:" << forwardMatrixTmp[chromosomeId] << endl;
                        }
                    } else {
                        forwardMatrixTmp[chromosomeId] += forwardMatrixTmpPre.at(chromosomeIdPre) * recombProb_ * sequenceStatus;
                        if (debug_)
                        {
                            cerr << kmerHash << " forwardPre:" << forwardMatrixTmpPre.at(chromosomeIdPre) << " recombPro:" << recombProb_ << " sequenceStatus:" << sequenceStatus << " forward:" << forwardMatrixTmp[chromosomeId] << endl;
                        }
                    }
                }
            }
        }

        // Record the hash value of the previous node
        preKmerHash = kmerHash;
    }
    return 0;
}


/**
 * @author zezhen du
 * @date 2023/4/20
 * @version v1.0
 * @brief forward
 * 
 * @return 0
**/
int HMM::get_path() {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Backtracking to obtain the most probable path.\n";
    // Find the haplotype with the highest probability for each kmer
    map<uint64_t, uint32_t> hapFreMap;  // map<hapID, frequence>
    for (const auto& [kmerHash, hapInfoMap] : forwardMatrix_) {  // map<kmerHash, map<chrId, score> >
        long double maxScore = 0.0;
        uint64_t maxHapID;

        // Find the maximum value and the corresponding haplotype ID
        for (const auto& [chrId, score] : hapInfoMap) {  // map<chrId, score>
            if (score > maxScore) {
                maxScore = score;
                maxHapID = chrId;
            }
        }

        hapFreMap[maxHapID]++;
    }

    // Key-value pair swap
    for (const auto& [hapID, frequence] : hapFreMap) {  // map<hapID, frequence>
        auto& HapVec = freHapVecMap_[frequence];
        HapVec.push_back(hapID);
    }

    // Print log
    if (debug_) {
        for (const auto& [frequence, hapIDVec] : freHapVecMap_) {  // map<frequence, vector<hapID> >
            cerr << frequence << "\t";
            for (const auto& hapID : hapIDVec) {
                cerr << get<0>(chromosomeIdLenMap_.at(hapID)) << "\t";
            }
            cerr << endl;
        }
    }
    
    return 0;
}



/**
 * @author zezhen du
 * @date 2023/4/20
 * @version v1.0
 * @brief forward
 * 
 * @param hapNum                The number of haplotypes to be obtained
 * 
 * @return outHapVec       vector<string>  vector<hap>
**/
vector<string> HMM::get_result( 
    const int& needHapNum_
) {
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Obtaining haplotype combinations.\n";
    // Output
    vector<string> outHapVec;  // vector<hap>
    for (auto it = freHapVecMap_.rbegin(); it != freHapVecMap_.rend(); ++it) {  // map<frequence, vector<hapID> >
        for (const auto& hapID : it->second) {  // vector<hapID>
            outHapVec.push_back(get<0>(chromosomeIdLenMap_.at(hapID)));

            if (outHapVec.size() == needHapNum_) {
                return outHapVec;
            }
        }
    }
    
    return outHapVec;
}

#endif
