
#pragma once

typedef unsigned int uint32;
typedef unsigned short uint16;
typedef unsigned char byte;
typedef unsigned long long ullong;

typedef vector<byte>     ByteAry;
typedef vector<int>      IntAry;
typedef vector<double>   DblAry;
typedef vector<string>   StrAry;
typedef map<byte,int>    ByteIntMap;
typedef map<int,int>     IntIntMap;
typedef map<byte,byte>   ByteByteMap;
typedef map<byte,double> ByteDblMap;
typedef map<int,double>  IntDblMap;
typedef map<string,int>  StrIntMap;
typedef map<string,string>    StrStrMap;
typedef map<string,double>    StrDblMap;
typedef map<byte,string>      ByteStrMap;
typedef vector< vector<int> > IntIntAry;
typedef map<string, IntAry>   StrIntAryMap;

typedef vector<byte>::iterator     ByteAryIter;
typedef vector<int>::iterator      IntAryIter;
typedef vector<double>::iterator   DblAryIter;
typedef vector<string>::iterator   StrAryIter;	
typedef map<byte,int>::iterator    ByteIntMapIter;
typedef map<int,int>::iterator     IntIntMapIter;
typedef map<byte,byte>::iterator   ByteByteMapIter;
typedef map<byte,double>::iterator ByteDblMapIter;
typedef map<int,double>::iterator  IntDblMapIter;
typedef map<string,int>::iterator  StrIntMapIter;
typedef map<string,string>::iterator  StrStrMapIter;
typedef map<string,double>::iterator  StrDblMapIter;
typedef map<byte,string>::iterator    ByteStrMapIter;
typedef map<string, IntAry>::iterator StrIntAryMapIter;

typedef vector< vector<int> >           IntIntAry;
typedef vector< vector<int> >::iterator IntIntAryIter;

#define TT template<class T>
#define TT2 template<class T1, class T2>
#define TTi template<class T> inline

// ------------------------------------------------------------------------------------- filename utils -------------

// DirSepQ(c) tests whether $c is a directory separator.
inline bool DirSepQ(int c) {
#ifdef WIN32
    return c=='\\' || c=='/' || c==':'; // recognize PC, Unix, Mac
#elif defined(MACOS)
    return c=='\\' || c=='/' || c==':'; // recognize PC, Unix, Mac
#else
    return c=='/';                      // on Unix, \ and : are allowed in pathnames
#endif
}

#define TAB "\t"
#define tab "\t"
static const char* nl = "\n";

// -------------------------------------------------------------------------------------- string functions ---------


// TT string ToString(const T& v, char sep=' ');
string ToString(const IntAry& v, char sep=' ');
string ToString(const DblAry& v, char sep=' ');

string padright(string str, int n, char sym=' ');
string padleft(string str, int n, char sym=' ');
bool GetWords(string& line, StrAry& words, string sep);

inline void ltrim(string& str) { str.erase(0, str.find_first_not_of(" \t\r\n")); }
inline void rtrim(string& str) { str.erase(str.find_last_not_of(" \t\r\n")+1); }
inline void trim(string& str) { ltrim(str); rtrim(str); }
inline void toupper(string& str) { std::transform(str.begin(),str.end(),str.begin(), (int(*)(int))toupper); }
inline void tolower(string& str) { std::transform(str.begin(),str.end(),str.begin(), (int(*)(int))tolower); }

// --------- file/filename utils -----------------------------------------------------------------------------------

#ifdef WIN32
#define DIR_SEPARATOR '\\'
#elif defined(MACOS)
#define DIR_SEPARATOR ':'
#else
#define DIR_SEPARATOR '/'
#endif

string StripExt(string& a);
string GetFileName(string fn);
string GetFilePath(string fn);

bool StringToFile(ostringstream& str, string fn);
bool FileToString(string fn, string& str);
bool DirQ(const char* dir);

void qlog(ostringstream& oss);
inline string getstr(StrAry& list) { ostringstream s; for(int k=0; k<list.size(); k++) s << (k ? "," : "") << list[k]; return s.str(); }

TTi void FillInc(vector<T>& list, size_t n) {
	list.resize(n);
	for(size_t k=0; k<n; k++) list[k] = k;
}
struct Elapsed {
	time_t _start, _end;
	Elapsed(bool startQ=false) { if (startQ) start(); }
	void start() { time(&_start); }
	double get() {  time (&_end); return difftime(_end, _start); }
	ostringstream& operator<<(ostringstream& oss) {
		oss << "elapsed : " << get() << " secs.";
		return oss;
	}
};

// ----------------------------- crap like this just to sort two arrays in parallel? great job STL! -----------------

TT struct SortPair {
	T v; int idx;
	bool operator<(const SortPair<T>& a) const { return v < a.v; }
};
TT struct OrderedSort {
	vector<SortPair<T> > _order;
	// this will sort 'ary' and reorder 'order' in parallel
	// ary.size must be equal to order.size
	OrderedSort(vector<T>& ary, IntAry& order) {	
		_order.resize(ary.size());
		for (size_t k=0; k<ary.size(); k++) {
			_order[k].v = ary[k];
			_order[k].idx = order[k];
		}
		std::sort(_order.begin(), _order.end());
		get(ary, order);
	}
	OrderedSort(const vector<T>& ary) {
		_order.resize(ary.size());
		for (size_t k=0; k<ary.size(); k++) {
			_order[k].v = ary[k];
			_order[k].idx = k;
		}
		std::sort(_order.begin(), _order.end());
	}
	void get(vector<T>& ary, IntAry& order) {
		ary.resize(_order.size());
		order.resize(_order.size());
		for (size_t k=0; k<_order.size(); k++) {
			ary[k] = _order[k].v;
			order[k] = _order[k].idx;
		}
	}
};

// --------------------------------------------------------------------------------------- profiler ----------------

struct Profiler {
	vector<Elapsed> timers;
	DblAry times;
	StrAry names;
	int add(string name, bool startQ=false) {
		int id = timers.size();
		timers.push_back(Elapsed());
		times.push_back(0);
		names.push_back(name);
		if (startQ) timers[id].start();
		return id;
	}
	void start(int id) { timers[id].start(); }
	void stop(int id) {
		double diff = timers[id].get();
		times[id] += diff;
	}
	string str() {
		ostringstream oss;
		oss << endl << "Profiler ::" << endl;
		for (int k=0; k<timers.size(); k++) {
			double p = times[0] > DBL_EPSILON ? 100 * times[k] / times[0] : 0;
			oss << names[k] << " -- " << times[k] << " secs. -- " << p << "%" << endl;
		}
		return oss.str();
	}
};
bool GetFileSize(const string& fn, ullong& filesize);
inline size_t GetPairIndex(int i, int j) { return i + j*(j+1)/2; }
void str_replace(string& in, const string& from, const string& to);

// ---------------------------------------------------------------------- CmdArgs -------------------------------

struct CmdArgs {
	StrAry _names;
	StrAry _values;
	StrIntMap _name2idx;
	StrIntMapIter _it;
	bool _goodQ;
	CmdArgs(int argc, char** argv) {
		_goodQ = false;
		for(int k=1; k<argc; k++) {
			string s(argv[k]);
			size_t found = s.find_first_of('=');
			if (found==string::npos) { cout << "arguments must be in the form 'argument=value'" << endl; return; }
			string name(s.substr(0, found));
			string value(s.substr(found+1));
		//	cout << _names.size() << ". " << name << "===" << value << endl;
			_name2idx[name] = _names.size();
			_names.push_back(name);
			_values.push_back(value);
		}
		_goodQ = true;
	}
	TT bool GetValue(string name, T dflt, T& value) {
		_it = _name2idx.find(name);
		if (_it==_name2idx.end()) { value = dflt; return false; }
		istringstream iss(_values[_it->second]);
		if (iss >> value) return true;
		value = dflt;
		return false;
	}
};

// ---------------------------------------------------------------------- FileLinesIter -------------------------------

struct FileLinesIter {
	ifstream _ifs;
	string _line;
	Elapsed _ela;
	bool _echoQ;
	int _polling;
	ullong _nline;
	ullong _fsize, _bytes, _pos;
	void close() { _ifs.close(); }
	bool open(string fn, bool echoQ=false) {
		_echoQ = echoQ; _fsize = 0;
		ullong fsize;
		GetFileSize(fn, fsize); _fsize = fsize;
		_ifs.open(fn.c_str(), ifstream::in);
		if (!_ifs.good()) { cout << endl << "Can't open " << fn << " for reading" << endl; return false; }
		_nline = _bytes = _pos = 0;
		_polling = 10; // print progress every 10 seconds
		if (echoQ) _ela.start();
		return true;
	}
	bool getline() {
		if (_echoQ && _nline%100==0) {
			double ela = _ela.get();
			if (ela>_polling) {
				_pos += _bytes;
				double speed = _bytes / (1024 * ela);
				double prc = 100 * (double(_pos) / double(_fsize));
				cout << "read " << std::fixed << std::setprecision(2) << prc << "% at " << speed << " kb/sec" << endl;
				_bytes = 0; _ela.start();
			}
		}
	//	if (!std::getline(_ifs, _line, '\n')) return false;
		if (!std::getline(_ifs, _line)) return false;
		_bytes += _line.size();
		_nline++;
		return true;
	}
};

// -------------------------------------------------------------------- CRC ------------------------------------------

extern const uint32 crc32_table[256];

inline int CRC(const char *data, int n, int crc=0) throw() {
	const uint32 *t = crc32_table;
	for (int i=0; i<n; i++) crc = (crc << 8) ^ t[((int) (crc >> 24) ^ *data++ ) & 0xff];
	return crc;
}
inline int CRC(int crc, char b) throw() {
	return (crc << 8) ^ crc32_table[((int) (crc >> 24) ^ b ) & 0xff];
}
// Use CRC32 to return a hash value of N bytes in memory
inline int CRCHash(const void *data, int n, int nbuckets) throw() {
	return uint32(CRC((const char*)data,n)) % nbuckets;
}
// Use CRC32 to return a hash value of an object
TTi int TCRCHash(const T &t, int nbuckets) throw() {
	return uint32(CRC((const char*)&t,sizeof(T))) % nbuckets;
}

// ------------------------------------------------------------------
// searches for a largest element smaller (or equal) to given element
// Example:
// 0 1 2 3 4 5 6 7 8 9
// BinSearch(3) == 3
// BinSearch(3.5) == 3
// BinSearch(-1) == 0
// BinSearch(100) == 9
// simply narrows the search down to 1 element.
// 1 5 5 5 7
// BinSearch(6) == 3
//
TTi int BinSearch(vector<T> &ary, T value) {
	size_t lo=0, hi=ary.size();
	while (hi > lo+1) {
		size_t mid = (size_t)((hi+lo)/2);
		if (ary[mid] > value)
			 hi = mid;
		else lo = mid;
	}
	return lo;
}

// ------------------------------------------------------------------------------------------

/*
struct Elapsed2 {
    struct timeval start, end;
    long mtime, seconds, useconds;
	void Start() { gettimeofday(&start, NULL); }
	void Elapsed() {
#ifdef WIN32
		SYSTEMTIME st;
		GetSystemTime(&st);
#else
		gettimeofday(&end, NULL);
		seconds  = end.tv_sec  - start.tv_sec;
		useconds = end.tv_usec - start.tv_usec;
#endif
	}
	ostringstream& operator<<(ostringstream& oss) {
	    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
		oss << mtime;	// printf("Elapsed time: %ld milliseconds\n", mtime);
		return mtime;
	}
}
*/
