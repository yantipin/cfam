
#pragma once

struct CluID {	// -------------------------------- one cluster, list of items -------------------------------------
	IntAry ids;
	double  entropy;
	double  energy;
	CluID() { entropy = energy = 0.; }
	void operator=(const CluID& clu) { ids = clu.ids; entropy = clu.entropy; energy = clu.energy; }
	void Reset() { entropy = energy = 0.; ids.clear(); }
	void Add(int id) { ids.push_back(id); sort(ids.begin(), ids.end()); }
	void Sort() { if (ids.size() < 2) return; sort(ids.begin(), ids.end()); }
	string str() {
		ostringstream str;
		for(int k=0; k<ids.size(); k++) str << (k?",":"") << ids[k];
		return str.str();
	}
};
struct Clu {  // ---------------------------------- clustering as a list of clusters --------------------------------
	vector<CluID>	clu;
	void Reset() { clu.clear(); }
	bool EmptyCluQ(int nclu) {
		if (nclu >= clu.size()) return true;
		return clu[nclu].ids.size()==0;
	}
	int Total() {
		int len=0;
		for (int k=0; k<clu.size(); k++) len += clu[k].ids.size();
		return len;
	}
	void GetAsList(vector<int> &list) {
		list.clear();
		for (int k=0; k<clu.size(); k++)
			for (int j=0; j<clu[k].ids.size(); j++)
				list.push_back(clu[k].ids[j]);
	}
};
struct CluSettings { // --------------------------------------------------------------------------------------------
	double   A1begin;
	double   A1end;
	int      A1steps;
	double   space_threshold;
	int      ncluA;
	CluSettings() { space_threshold = 0.; A1begin = A1end = 0.75; A1steps = 1; ncluA = 0; }
	CluSettings(double a1beg, double a1end, int steps, int nclu=0) {
		space_threshold=0.; A1begin=a1beg; A1end=a1end; A1steps=steps; ncluA=nclu; }
};

struct CluResults : public Clu {// ---------------------------------------------------------------------------------
	double		A1;				// best A1
	IntAry		columns;		// list of columns being used for this clustering
	IntAry		colsbyent;		// columns sorted by clustRealEnt
	double		whole_entropy;	// sum of all clusters entropy
	DblAry		clustEnt;		// columns normalized clusters entropy
	IntAry		clustEntRanks;	// columns normalized clusters entropy ranks
	DblAry		clustRealEnt;	// columns combinatorial entropy
	DblAry		clustRealEntNG;	// same as previous, only gapps (symbol '-') are skipped
	DblAry		clustZS;		// columns z-score
	IntAry		clustZSRanks;	// columns z-score ranks
	CluResults() { Reset(); }
	CluResults(const CluResults& cr) { *this = cr; }
	void CalcWholeEntropy() { whole_entropy = 0; for(int k=0; k<clu.size(); k++) whole_entropy += clu[k].entropy; }
	void clear(bool clearcolumnsQ=true) {
		if (clearcolumnsQ) columns.clear();
		A1 = whole_entropy = 0.;
		clu.clear();
		clustEnt.clear();
		clustEntRanks.clear();
		clustRealEnt.clear();
		clustRealEntNG.clear();
		clustZS.clear();
		clustZSRanks.clear();
	}
};

// ----------------------------------------------------- clu distance map ---------------------------

struct CluPairInfo {	// only pairs of clusters get merged
	double entropy;
	double energy;
	CluPairInfo() { set(0, 0); }
	CluPairInfo(double ent, double eng) { set(ent, eng); }
	void set(double ent, double eng) { entropy = ent; energy = eng; }
};

#ifdef WIN32 // MSVC
#include <unordered_map>
#elif __APPLE__
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

struct CluIDHash;
struct CluIDEqual { bool operator() (const IntAry &c1, const IntAry &c2) const { return c1 == c2; }};
typedef std::unordered_map<IntAry, CluPairInfo, CluIDHash, CluIDEqual> CluIDMap;
static int bucket_count = 0;

struct CluIDHash {
	CluIDMap* _map;	// need unordered_map bucket count
	CluIDHash() : _map(NULL) { }
	void SetMapInstance(CluIDMap* map) {
#ifdef WIN32
		_map = map;
#endif
	}
	long operator() (const IntAry &a) const {
#ifdef WIN32
		return CRCHash((const void *)&a[0], a.size() * sizeof(int), _map->bucket_count());
#else
	//	return CRCHash((const void *)&a[0], a.size() * sizeof(int), bucket_count);
		return CRCHash((const void *)&a[0], a.size() * sizeof(int), 2097152);	// <-- this bucket count never changes on Windows when processing msa.txt (1000 sequences)
#endif
	}
};

struct DistanceMap { // --------------------------------------------------- clu distance map -------------------
	CluIDMap m;
	CluIDMap::iterator it;
	DistanceMap() {
#ifdef WIN32
		CluIDHash& hash = m.comp._Hashobj;
		hash.SetMapInstance(&m);
#endif
	}
//	// well, this doesn't work
//	void UpdateBucketCount() {
//		const CluIDMap::hasher &hfn = m.hash_function();
//		const CluIDHash &hash = static_cast<const CluIDHash &>(hfn);
//		CluIDHash* hp = (CluIDHash*)&hash;
//		hp->UpdateBucketCount(m.bucket_count());
//	}
	void Init() {
		bucket_count = m.bucket_count();
		m.max_load_factor(0.5);
	}
	bool find(IntAry& clu, double& entropy, double& energy) {
		it = m.find(clu);
		if (it==m.end()) return false;
		entropy = it->second.entropy;
		energy  = it->second.energy;
		return true;
	}
	void add(IntAry& clu, double entropy, double energy) {
	//	bucket_count = m.bucket_count();
		m[clu] = CluPairInfo(entropy, energy);
	}
	void clear() { m.clear(); }
	string str() {
		ostringstream s;
		s << m.bucket_count() << "," << m.max_bucket_count() << "," << m.load_factor();
		return s.str();
	}
};

