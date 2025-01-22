#ifndef bitmap_hpp
#define bitmap_hpp
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <bitset>

using namespace std;

// Bitmap Implementation
class BitMap
{
public:
	// The memory size of the bitmap is related to the data range
	BitMap(size_t range)
		:_bit(range / 32 + 1) {
		// Bitmap Capacity
		_cap = range;
	}

	// Loading a bitmap from memory
	BitMap(string BitMapFileName) {
		// Bitmap Capacity
		_cap = 0;

		// Open Binary File
		gzFile gzfp = gzopen(BitMapFileName.c_str(), "rb");

		if(!gzfp) {  // Failed to open file
			cerr << "[" << __func__ << "::" << getTime() << "] '" << BitMapFileName << "': No such file or directory." << endl;
			exit(1);
		} else {
			string information;
			char line[100];

			while (gzgets(gzfp, line, 100)) {
				information = line;

				if(information.empty()) {
					continue;
				} else {
					// The first line is the bitmap capacity
					string::size_type idx = information.find("_cap:");
					if (idx != string::npos) {
						_cap = stol(split(strip(information, '\n'), ":")[1]);
					} else {
						_bit.push_back(stol(strip(information, '\n')));
					}
				}
				
				// Clear char
				memset(line,'\0',sizeof(line));
			}
		}
		// Freeing up memory
		gzclose(gzfp);


		// Check the size of _num. If it is 0, it means the bitmap file is incorrect.
		if (_cap == 0) {
			cerr << "[" << __func__ << "::" << getTime() << "] " 
				 << BitMapFileName 
				 << ": Missing bitmap capacity information (_cap:number)." 
				 << endl;
			exit(1);
		}
	}

	void set(const size_t num) {
		// Calculates the subscript in an array
		int idx = num / 32;
		// Calculate the subscript position of num in the corresponding subscript integer
		int bitIdx = num % 32;
		// cout << idx << " " << bitIdx << endl;
		// Set the corresponding bit position to 1
		_bit[idx] |= 1 << bitIdx;
	}

	bool find(const size_t num) {
		int idx = num / 32;
		int bitIdx = num % 32;
		// cout << idx << " " << bitIdx << endl;
		return (_bit[idx] >> bitIdx) & 1;
	}

	// Set the corresponding bit position to 0
	void reset(const size_t num) {
		int idx = num / 32;
		int bitIdx = num % 32;
		_bit[idx] &= ~(1 << bitIdx);
	}

	// Bitmap Capacity
	size_t get_cap() {
		return _cap;
	}

	// Amount of bitmap storage elements
	size_t get_num() {
		size_t _num = 0;
		for (auto it : _bit) {
			if (it == 0) {  // If it is 0, skip
				continue;
			}
			
			bitset<32> tmpIt(it);
			for (int i = 0; i < tmpIt.size(); i++) {
				if (tmpIt[i] != 0) {
					_num++;
				}
			}
		}
		return _num;
	}

	// Save the bitmap to memory
	void save(string BitMapFileName) {
		string outBitStr;

		// The first line of the file is the bitmap capacity '_num'
		outBitStr += "_cap:" + to_string(_cap) + "\n";
		
		for (auto it : _bit) {
			outBitStr += to_string(it) + "\n";
		}

		// Save to file
		gzFile gzfp = gzopen(BitMapFileName.c_str(), "wb");
		gzwrite(gzfp, outBitStr.c_str(), outBitStr.length());

		// Freeing up memory
		gzclose(gzfp);
	}
private:
	vector<int> _bit;
	size_t _cap;  // Bitmap Capacity
};

#endif /* bitmap_hpp */