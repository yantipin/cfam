
#pragma once

struct Sequence {
	string name;
	string seq;
	void FromFasta(string& fasta) {
		size_t to=0, from=0;
		while(string::npos!=(to=fasta.find_first_of('\n', from))) {
			string line(fasta.substr(from, to-from)); from=to+1;
			if (!line.size()) continue;
			if (line[0]=='>') { name = line.substr(1); continue; }
			seq += line;
		}
	}
};

// -------------------------------------------------------------- protein region name (protID/from-to) ---

struct ProtRgnName {	// EGFR_HUMAN/1-300
	string protID;
	int from;
	int to;
	ProtRgnName() : from(0), to(0) {}
	ProtRgnName(string str) { parse(str); }
	string str() { ostringstream s; s << protID << '/' << from << '-' << to; return s.str();}
	bool validQ() { return protID.size() && from>0 && to>from; }
	bool parse(string name) {
		protID = ""; from = to = 0;
		size_t p = name.find_first_of('/');
		if (p==string::npos) { protID = name; return false; }
		protID = name.substr(0, p);
		size_t n = name.find_first_of('-', p);
		if (n==string::npos) return false;
		istringstream s(name.substr(p+1, n-p+1)); s >> from;
		s.str(name.substr(n+1)); s >> to;
		return from>0 && to>from;
	}
};

// -------------------------------------------------------------- Hmmer ----------------------------------

struct ParseHmm {
	istringstream _iss;
	StrAry _domains;
	DblAry _cvalues;
	IntAry _hfrom;
	IntAry _hto;
	string _str;
	StrDblMap _bestdomains;		// only domains with lowest c-value here
	StrDblMapIter _it;
	bool DomainsQ() { return 0!=_bestdomains.size(); }
	bool BestDomainQ(string domain) { return _bestdomains.find(domain) != _bestdomains.end(); }
	bool GetDomainCValue(string domain, double& cvalue) {
		_it = _bestdomains.find(domain);
		if (_it==_bestdomains.end()) return false;
		cvalue = _it->second;
		return true;
	}
	bool ParseDomainHits(string fn) {
		string domain, lastdomain;
		double cval, mincval = 1000;
		ostringstream rgn;
		int hfrom, hto;
		FileLinesIter it;
		if (!it.open(fn, true)) return false;
		StrStrMap proteins;
		StrDblMap cvalues;
		while(it.getline()) {
			int l = it._line.length();
			if (l<2) continue;
			if (it._line[0]=='#') continue;
			_iss.str(it._line);
			_iss >> domain;	SkipWords(10);
			_iss >> cval;	SkipWords( 5);
			_iss >> hfrom;
			_iss >> hto;
			_domains.push_back(domain);
			_cvalues.push_back(cval);
			_hfrom.push_back(hfrom);
			_hto.push_back(hto);
			rgn.str(""); rgn << hfrom << '-' << hto;
			if (domain!=lastdomain) {
			//	cout << "adding " << domain << "/" << hfrom << '-' << hto << " " << cval << endl;
				proteins[domain] = rgn.str();
				cvalues[domain] = cval;
				mincval = cval;
				lastdomain = domain;
				continue;
			}
			if (cval < mincval) {
			//	cout << "updating " << domain << "/" << hfrom << '-' << hto << " " << cval << endl;
				mincval = cval;
				proteins[domain] = rgn.str();
				cvalues[domain] = cval;
			}
		}
		_bestdomains.clear();
		for (StrStrMapIter it=proteins.begin(); it!=proteins.end(); it++) {
			string name(it->first + '/' + it->second);
			_bestdomains[name] = cvalues[it->first];
		}
		return true;
	}
	void SkipWords(int n) { for(int k=0; k<n; k++) _iss >> _str; }
};

// ---------------------------------------------------------------------------------------------------------------

static const int lenAA      = 20;
static const int lenAAext   = 21;
static const int lenAAallow = 23;
static const char AA[lenAA]           = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' };
static const char AAext[lenAAext]     = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U' };
static const char AAallow[lenAAallow] = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'X', '-' };

// ---------------------------------------------------------------------------------------------------------------


