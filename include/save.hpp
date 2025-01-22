#ifndef save_hpp
#define save_hpp

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <tuple>
#include <iomanip>
#include "zlib.h"
#include <unordered_map>
#include "strip_split_join.hpp"
#include "get_time.hpp"

using namespace std;

namespace SAVE {
    /**
     * @brief Open vcf file
     * 
     * @param outputFileName
     * 
     * @return
    **/
    class SAVE {
    private:
        // vcf file
        string outputFileName_;

        // Output file stream
        ofstream fpO;
        gzFile gzfpO;
    public:
        SAVE() {}
        SAVE(string aliFileName) {
            outputFileName_ = aliFileName;

            if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos) {
                // Output file stream
                gzfpO = gzopen(outputFileName_.c_str(), "wb");
                if(!gzfpO) {
                    cerr << "[" << __func__ << "::" << getTime() << "] " 
                        << "'" << outputFileName_ << "': No such file or directory." << endl;
                    exit(1);
                }
            } else if (outputFileName_.size() > 0) {
                // Output file stream
                fpO.open(outputFileName_, ios::out);
                if(!fpO) {
                    cerr << "[" << __func__ << "::" << getTime() << "] " 
                        << "'" << outputFileName_ << "': No such file or directory." << endl;
                    exit(1);
                }
            }
        }
        ~SAVE() {
            if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos) {
                // Close File
                gzclose(gzfpO);
            } else if (outputFileName_.size() > 0) {
                // Close File
                fpO.close();
            }
        }
        int save(string & outTxt);
    };


    /**
     * @brief storage
     * 
     * @return int
    **/
    int SAVE::save(
        string & outTxt
    ) {
        outTxt = strip(outTxt, '\n');  // Remove line breaks

        if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos) {
            outTxt += "\n";
            gzwrite(gzfpO, outTxt.c_str(), outTxt.length());
        } else if (outputFileName_.size() > 0) {
            fpO << outTxt << endl;
        } else {
            cout << outTxt << endl;
        }
        
        return 0;
    }
} // namespace SAVE

#endif