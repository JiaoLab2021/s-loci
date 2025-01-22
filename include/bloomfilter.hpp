#ifndef bloomfilter_hpp
#define bloomfilter_hpp
#include <string>


using namespace std;

struct HashFun1 {
	// Calculate each character of the string to get a hash value
	size_t operator()(const string& str) {
		size_t hash = 0;
		for (const auto& ch : str) {
			hash = hash * 131 + ch;
		}
        // cout << "HashFun1: " << str << " -> " << hash << endl;
		return hash;
	}
};
struct HashFun2 {
	size_t operator()(const string& str) {
		size_t hash = 0;
		for (const auto& ch : str) {
			hash = hash * 65599 + ch;
		}
        // cout << "HashFun2: " << str << " -> " << hash << endl;
		return hash;
	}
};
struct HashFun3 {
	size_t operator()(const string& str) {
		size_t hash = 0;
		for (const auto& ch : str) {
			hash = hash * 1313131 + ch;
		}
        // cout << "HashFun3: " << str << " -> " << hash << endl;
		return hash;
	}
};

struct HashFun4 {
    size_t operator()(const string& str) {
        size_t hash = 0;
        size_t magic = 63689;	
        for (const auto& ch : str) {
            hash = hash * magic + ch;
            magic *= 378551;
        }
        // cout << "HashFun4: " << str << " -> " << hash << endl;
        return hash;
    }
};

struct HashFun5 {
    size_t operator()(const string& str) {
        size_t hash = 0;
        size_t magic = 63689;
        long i = 0;	
        for (const auto& ch : str) {
            if ((i & 1) == 0) {
                hash ^= ((hash << 7) ^ ch ^ (hash >> 3));
            } else {
                hash ^= (~((hash << 11) ^ ch ^ (hash >> 5)));
            }
            i++;
        }
        // cout << "HashFun5: " << str << " -> " << hash << endl;
        return hash;
    }
};

struct HashFun6 {
    size_t operator()(const string& str) {
        if(str.length() == 0) // An empty string returns 0
        {
		    return 0;
        }

        size_t hash = 1315423911;
        for (const auto& ch : str)
        {
            hash ^= ((hash << 5) + ch + (hash >> 2));
        }
        // cout << "HashFun6: " << str << " -> " << hash << endl;
        return hash;
    }
};

struct HashFun7
{
    size_t operator()(const string& str)
    {
        if(str.length() == 0) {  // An empty string returns 0
		    return 0;
        }

        size_t hash = 1315423911;
        for (const auto& ch : str) {
            hash = ((hash << 5) ^ (hash >> 27)) ^ ch;
        }
        // cout << "HashFun7: " << str << " -> " << hash << endl;
        return hash;
    }
};

struct HashFun8 {
    size_t operator()(const string& str) {
        if(str.length() == 0) {  // An empty string returns 0
		    return 0;
        }

        size_t hash = 2166136261;
        for (const auto& ch : str) {
            hash *= 16777619;
		    hash ^= ch;
        }
        // cout << "HashFun8: " << str << " -> " << hash << endl;
        return hash;
    }
};

struct HashFun9 {
    size_t operator()(const string& str) {
        if(str.length() == 0) {  // An empty string returns 0
		    return 0;
        }

        size_t hash = 5381;
        for (const auto& ch : str) {
            hash += (hash << 5) + ch;
        }
        // cout << "HashFun9: " << str << " -> " << hash << endl;
        return hash;
    }
};

struct HashFun10 {
    size_t operator()(const string& str) {
        if(str.length() == 0) {  // An empty string returns 0
		    return 0;
        }

        size_t hash = 5381;
        for (const auto& ch : str) {
            hash = hash * 33 ^ ch;
        }
        // cout << "HashFun10: " << str << " -> " << hash << endl;
        return hash;
    }
};

struct HashFun11 {
    size_t operator()(const string& str) {
        static const size_t TotalBits		= sizeof(size_t) * 8;
        static const size_t ThreeQuarters	= (TotalBits  * 3) / 4;
        static const size_t OneEighth		= TotalBits / 8;
        static const size_t HighBits		= ((size_t)-1) << (TotalBits - OneEighth);	
        
        size_t hash = 0;
        size_t magic = 0;
        for (const auto& ch : str) {
            hash = (hash << OneEighth) + ch;
            if ((magic = hash & HighBits) != 0) {
                hash = ((hash ^ (magic >> ThreeQuarters)) & (~HighBits));
            }
        }
        // cout << "HashFun11: " << str << " -> " << hash << endl;
        return hash;
    }
};

struct HashFun12 {
    size_t operator()(const string& str) {
        static const size_t TotalBits		= sizeof(size_t) * 8;
        static const size_t ThreeQuarters	= (TotalBits  * 3) / 4;
        static const size_t OneEighth		= TotalBits / 8;
        static const size_t HighBits		= ((size_t)-1) << (TotalBits - OneEighth);	

        size_t hash = 0;
        size_t magic = 0;
        for (const auto& ch : str) {
            hash = (hash << OneEighth) + ch;
            if ((magic = hash & HighBits) != 0) {
                hash ^= (magic >> ThreeQuarters);
                hash &= ~magic;
            }
        }
        // cout << "HashFun12: " << str << " -> " << hash << endl;
        return hash;
    }
};

template<class T, class HashFun1, class HashFun2, class HashFun3, 
         class HashFun4, class HashFun5, class HashFun6, 
         class HashFun7, class HashFun8, class HashFun9, 
         class HashFun10, class HashFun11, class HashFun12>
class BloomFilter
{
public:
    // Newly constructed bitmap
	BloomFilter(const size_t num)
		:_bit(5 * num)
		, _bitCount(5 * num)
	{}

    // Loading a bitmap from memory
    BloomFilter(string BitMapFileName)
        :_bit(BitMapFileName) {
        _bitCount = _bit.get_cap();
    }

	void set(const T& val) {
        // cout << _bitCount << endl;
		HashFun1 h1;
		HashFun2 h2;
		HashFun3 h3;
        HashFun4 h4;
		HashFun5 h5;
		HashFun6 h6;
        HashFun7 h7;
		HashFun8 h8;
		HashFun9 h9;
        HashFun10 h10;
		HashFun11 h11;
		// HashFun12 h12;
        
        
		size_t idx1 = h1(val) % _bitCount;
		size_t idx2 = h2(val) % _bitCount;
		size_t idx3 = h3(val) % _bitCount;
        size_t idx4 = h4(val) % _bitCount;
		size_t idx5 = h5(val) % _bitCount;
		size_t idx6 = h6(val) % _bitCount;
        size_t idx7 = h7(val) % _bitCount;
		size_t idx8 = h8(val) % _bitCount;
		size_t idx9 = h9(val) % _bitCount;
        size_t idx10 = h10(val) % _bitCount;
		size_t idx11 = h11(val) % _bitCount;
        
        // cout << idx12 << endl;
		// int idx12 = h12(val) % _bitCount;
        
        
        // cout << idx1 << endl;
		_bit.set(idx1);
        // cout << idx2 << endl;
		_bit.set(idx2);
        // cout << idx3 << endl;
		_bit.set(idx3);
        // cout << idx4 << endl;
        _bit.set(idx4);
        // cout << idx5 << endl;
		_bit.set(idx5);
        // cout << idx6 << endl;
		_bit.set(idx6);
        // cout << idx7 << endl;
        _bit.set(idx7);
        // cout << idx8 << endl;
		_bit.set(idx8);
        // cout << idx9 << endl;
		_bit.set(idx9);
        // cout << idx10 << endl;
        _bit.set(idx10);
        // cout << idx11 << endl;
		_bit.set(idx11);
        // _bit.set(idx12);
        
		
	}

	bool find(const T& val) {
		HashFun1 h1;
		HashFun2 h2;
		HashFun3 h3;
        HashFun4 h4;
		HashFun5 h5;
		HashFun6 h6;
        HashFun7 h7;
		HashFun8 h8;
		HashFun9 h9;
        HashFun10 h10;
		HashFun11 h11;
		// HashFun12 h12;
		size_t idx1 = h1(val) % _bitCount;
		size_t idx2 = h2(val) % _bitCount;
		size_t idx3 = h3(val) % _bitCount;
        size_t idx4 = h4(val) % _bitCount;
		size_t idx5 = h5(val) % _bitCount;
		size_t idx6 = h6(val) % _bitCount;
        size_t idx7 = h7(val) % _bitCount;
		size_t idx8 = h8(val) % _bitCount;
		size_t idx9 = h9(val) % _bitCount;
        size_t idx10 = h10(val) % _bitCount;
		size_t idx11 = h11(val) % _bitCount;
		// int idx12 = h12(val) % _bitCount;

		if (!_bit.find(idx1))
			return false;
		if (!_bit.find(idx2))
			return false;
		if (!_bit.find(idx3))
			return false;
        if (!_bit.find(idx4))
			return false;
		if (!_bit.find(idx5))
			return false;
		if (!_bit.find(idx6))
			return false;
        if (!_bit.find(idx7))
			return false;
		if (!_bit.find(idx8))
			return false;
		if (!_bit.find(idx9))
			return false;
        if (!_bit.find(idx10))
			return false;
		if (!_bit.find(idx11))
			return false;
		// if (!_bit.find(idx12))
		// 	return false;

		return true; // Possibly
	}

    // Save bitmap to memory
    void save(string BitMapFileName) {
		_bit.save(BitMapFileName);
	}

    // Amount of bitmap storage elements
    size_t get_num() {
		size_t bitMapNum = _bit.get_num();
        return bitMapNum;
	}

    // Bitmap Capacity
    size_t get_cap() {
        return _bitCount;
	}
private:
	BitMap _bit; // Bitmap
	size_t _bitCount;
};

#endif /* bloomfilter_hpp */