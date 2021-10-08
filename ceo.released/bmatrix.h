
#pragma once

#define HIST_ARYLEN		65535
#define GAMMLN_ARYLEN	1024
#define LETTERS			128
#define SYM_BYTE_IGNORE 'X'		// unknown residue in protein sequence
#define SYM_NOT_FOUND   255		// if you want to ask FIS of not observed residue

// ---------------------------------------------------------------------- gammaLN ------------------

static double gammln_cache_ary[GAMMLN_ARYLEN];
static double gammln(double xx) {
	double x, y, tmp, ser;
	static double cof[6] = {
		76.18009172947146, -86.50532032941677, 24.01409824083091,
		-1.231739572450155, 0.12086509738661179e-2, -0.5395239384953e-5 };
	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for (int j=0; j<=5; j++) ser += cof[j] / ++y;
	double val = -tmp + log(2.5066282746310005 * ser / x);
	return val;
}
static double gammln2(double xx) {
	if (xx <= 1.0) return 0.;
	return gammln(xx + 1);
}
static double cached_gammln(double v, bool treat_as_int) {
	if (treat_as_int) {
		if (v < GAMMLN_ARYLEN) return gammln_cache_ary[(int)v];
		else return gammln(v);
	}
	return gammln(v);
}
static void InitGammaLn() { for(int i=0; i<GAMMLN_ARYLEN; i++) gammln_cache_ary[i] = gammln(i); }

// --- BMatrix methods ------------------------------------------------------------------------

struct BMatrix {

	int					_width;
	int					_height;
	int					_maxcharsize;	// code of max character found + 1
	int					_histogram[HIST_ARYLEN];
	vector<ByteAry>		_transposed;	// transposed MSA for faster access
	IntAry				_symhist;
	IntAry				_colseqs;		// number of sequences in each column (non-gaps)
	vector<ByteIntMap>	_colnsym;
	vector<ByteDblMap>	_colfreq;
	vector<IntDblMap>	_chaos_ent;		// cached chaos entropy for each column and each possible clu length, int is within [2.._height]
	DistanceMap			_dmap;
	CluResults			_cr;
	DblAry				_a1list;
	vector<CluResults>	_shots;			// store cluster snapshot for each iteration
	vector<CluResults>	_bestshots;
	DblAry				_freqent;		// frequential entropy (conservation) - property of matrix
	DblAry				_freqent2;
	DblAry				_freqentNG;		// no-gaps
	bool				_prfQ;			// use profiler and print timings?
	IntAry				_row2clu;		// row index -> cluster index in _cr.clu
	byte				_gap;
	byte				_symignore;
	// encode characters so their codes are 0,1,2,... making hash tables much smaller
	ByteAry				_code2sym;		// code --> character
	ByteByteMap			_sym2code;
	bool				_echoQ;

	BMatrix() : _width(0), _height(0), _maxcharsize(0), _prfQ(false), _echoQ(false) {}

	// ---------------------------------------------------------- clusterize -----------------

	bool FromMSA(MSA& msa) {
		_width = msa._width;
		_height = msa._height;
		_gap = msa.gap;
		_transposed.resize(msa._width, ByteAry());
		_maxcharsize = 0;
		ByteByteMapIter it;
		_sym2code.clear();		// re-coding symbols (so their codes start from 0) makes code faster by 12%
		_code2sym.clear();
		for(int k=0;k<_transposed.size();k++) {
			_transposed[k].resize(_height, 0);
			for(int j=0;j<_height;j++) {
				byte v = msa.at(k,j);
				it = _sym2code.find(v);
				byte code = 0;
				if (it==_sym2code.end()) {
					code = _sym2code.size();
					if (code==254) {
						cout << "More than 254 different symbols found in the data." << endl;
						return false;
					}
					_sym2code[v] = code;
					_code2sym.push_back(v);
					if (v==msa.gap) _gap = code;
					if (v==SYM_BYTE_IGNORE) _symignore = code;
				} else
					code = it->second;
				if (code > _maxcharsize) _maxcharsize = code;
				_transposed[k][j] = code;
			}
		}
		_maxcharsize++;
		if (_echoQ) cout << "unique symbols : " << _maxcharsize << endl;
		return true;
	}

	// ---------------------------------------------------------------------------

	// keep sorted while merging
	void Merge(IntAry &a, IntAry &b, IntAry &c) {
		int i=0, j=0, k=0;
		c.resize(a.size() + b.size(), 0);
		while (i<a.size() && j<b.size()) {
			if (a[i]<b[j]) c[k++] = a[i++]; else c[k++] = b[j++];
		}
		while (i<a.size()) c[k++] = a[i++];
		while (j<b.size()) c[k++] = b[j++];
	}
	void MergeClusters(Clu& clu, int c1, int c2, double entropy, double energy) {
		IntAry c;
		Merge(clu.clu[c1].ids, clu.clu[c2].ids, c);
		clu.clu[c1].ids = c;
		clu.clu[c1].entropy = entropy;
		clu.clu[c1].energy = energy;
		clu.clu.erase(clu.clu.begin() + c2);
	}
	inline byte at(int column, int row) { return _transposed[column][row]; }
	inline byte symbol(int column, int row) { return _code2sym[ _transposed[column][row] ]; }

	void Clusterize(string msafn, CluSettings& cs, bool echoQ=true) {
		_echoQ = echoQ;
		if (!_width || !_height || !_transposed.size()) return;
		CacheColumnsFrequency();
		CacheChaosEntropy();
		GreedyClusterize(cs);
		GetColumnsCombinatorialEntropy(_cr);
		GetConservation();
		double elapsed = 0;
		string xmlfn(StripExt(msafn) + ".clu.xml");
		WriteCluXml(xmlfn, msafn, elapsed, cs);
	}
	void FillRowToClu() {
		_row2clu.resize(_height, 0);
		for(int ci=0; ci<_cr.clu.size(); ci++) {
			for(int k=0; k<_cr.clu[ci].ids.size(); k++) {
				_row2clu[ _cr.clu[ci].ids[k] ] = ci;
			}
		}
	}
	void GreedyClusterize(CluSettings &cs) {

		if (!_cr.columns.size()) FillInc(_cr.columns, _width);	// use all columns if not specified
		Elapsed ela(true);
		ostringstream oss;
		if (_echoQ) oss << endl << "Clustering,  A1 min..max, nsteps : " << cs.A1begin << ".." << cs.A1end << ", " << cs.A1steps << " steps" << endl; qlog(oss);

		int pstep = 0;
		_dmap.Init();
		_a1list.clear();
		_shots.clear();
		_bestshots.clear();
		IntAry tmpclu;

		Profiler prf;
		int prfcluID     = _prfQ ? prf.add("clustering", true) : 0;
		int prfEntropyID = _prfQ ? prf.add("entropy") : 0;
		int prfLookupID  = _prfQ ? prf.add("hash lookup") : 0;

		for (int nstep=0; nstep<cs.A1steps; nstep++) {

			_dmap.clear();
			_shots.clear();
			_cr.clear(false);
			_cr.clu.resize(_height);
			// start with N clusters where N = _height with one sequence each
			for(int k=0;k<_height;k++) _cr.clu[k].ids.push_back(k);

			double energy, entropy;
			double A1 = cs.A1begin + nstep * (cs.A1end - cs.A1begin) / (double)cs.A1steps;
			if (_echoQ) oss << "with A1 = " << A1 << " ---------------------------------" << endl; qlog(oss);

			while (_cr.clu.size() > 1) {

				int mc1 = -1, mc2 = -1;
				double best_energy = FLT_MAX, joined_entropy = 0.;
				for (int c1=0; c1<_cr.clu.size()-1; c1++) {
					for (int c2=c1+1; c2<_cr.clu.size(); c2++) {

						if (_prfQ) prf.start(prfLookupID);
						Merge(_cr.clu[c1].ids, _cr.clu[c2].ids, tmpclu);
						bool foundQ = _dmap.find(tmpclu, entropy, energy);
						if (_prfQ) prf.stop(prfLookupID);

						if (!foundQ) {
							if (_prfQ) prf.start(prfEntropyID);
							entropy = GetClustersEntropyDistance(_cr, c1, c2);
							if (_prfQ) prf.stop(prfEntropyID);

						//	oss << _cr.clu.size() << " " << c1 << " " << c2 << " " << entropy << endl; qlog(oss);
							int length = _cr.clu[c1].ids.size() + _cr.clu[c2].ids.size();
							energy = A1 * entropy / (double)(_width * length) + (1. - A1) * log((double)length);
							_dmap.add(tmpclu, entropy, energy);
						}
						if (fabs(energy - best_energy) > DBL_EPSILON && energy < best_energy) {
						//	oss << "better: " << energy << " < " << best_energy << " : clu " << mc1 << " + " << mc2 << endl; qlog(oss);
							best_energy = energy;
							mc1 = c1; mc2 = c2;
							joined_entropy = entropy;
						}
					}
					if (_echoQ && _cr.clu.size() == _height) { oss << '.'; qlog(oss); cout.flush(); }
				}
				MergeClusters(_cr, mc1, mc2, joined_entropy, best_energy);

				_cr.CalcWholeEntropy();
				_shots.push_back(_cr);
				if (_echoQ)
					oss << endl << "merging clu " << mc1 << " + " << mc2 << " energy: " << best_energy
						<< ", entropy: "  << _cr.whole_entropy << ", nclus: " << _cr.clu.size(); qlog(oss);
					//	<< "   " << _dmap.str();

				if (cs.ncluA && _cr.clu.size() == cs.ncluA) break;
			}

			// find clu snapshot with smallest entropy
			int best = 0;
			if (cs.ncluA) best = _shots.size() - 1;
			else {
				for (int j=1; j<_shots.size(); j++) {
					if (_shots[j].whole_entropy < _shots[best].whole_entropy) best = j;
				}
			}
			if (_echoQ) oss << endl << "Best entropy : " << _shots[best].whole_entropy << endl; qlog(oss);

			_bestshots.push_back( _shots[best] );
			_a1list.push_back(A1);
		}

		if (_prfQ) { prf.stop(prfcluID); cout << prf.str(); }

		// find clu snapshot with smallest entropy
		int best = 0;
		for (int j=1; j<_bestshots.size(); j++)
			if (_bestshots[j].whole_entropy < _bestshots[best].whole_entropy) best = j;

		_cr = _bestshots[best];
		if (_echoQ)
			oss << endl << "Best A1 : " << _cr.A1 << ", best entropy : " << _cr.whole_entropy
				<< ", " << ela.get() << " secs." << endl; qlog(oss);

	}

	// ---------------------------------------------------------------------------------- columns entropy ------------

	void GetColumnsCombinatorialEntropy(CluResults& res, bool calc_ranks=false) {

		if (!res.clu.size()) return;
		DblAry clu_entropies;
		double amin = FLT_MAX, amax = -FLT_MAX;
		ByteDblMap tmpfreq;
		res.clustEnt.resize(_width);
		res.clustRealEnt.resize(_width);
		res.clustRealEntNG.resize(_width);	// same as above, only no gaps

		if (calc_ranks) {
			res.clustEntRanks.resize(_width);
			for(int j=0; j<res.clustEntRanks.size(); j++)
				res.clustEntRanks[j] = (int)99.;
		}
		clu_entropies.resize(res.clu.size(), 0);

		double total_entropy = 0.;
		for (IntAryIter it=res.columns.begin(); it<res.columns.end(); it++) {
			int col = *it, whole_len = 0;
			double whole_clu_entropy = 0, whole_clu_entropy_NG = 0.;
			for(int i=0; i<res.clu.size(); i++) {
				int clu_len = res.clu[i].ids.size();
				if (clu_len == 1) continue;
				int clu_len_NG = 0;					// non-gappy length
				for (int lnIdx=0; lnIdx<res.clu[i].ids.size(); lnIdx++) {
					int line = res.clu[i].ids[lnIdx];
					byte symbol = at(col, line);
					// meaning of SYM_BYTE_IGNORE is "data not available" -- different from the meaning of a gap in MSA
					if (symbol == _symignore) continue;
					if (symbol!=_gap) ++clu_len_NG;
					ByteDblMapIter it = tmpfreq.find(symbol);
					if (it==tmpfreq.end()) tmpfreq[symbol] = 1; else tmpfreq[symbol]++;
				}
				double chaos_entropy   = GetChaosEntropy(col, clu_len);
				double chaos_entropyNG = GetChaosEntropy(col, clu_len_NG);
				double columnEntropyNG = 0;
				double columnEntropy   = GetFreqHashEntropyNG(tmpfreq, true, columnEntropyNG);
				double entropy         = columnEntropy   - chaos_entropy;
				double entropyNG       = columnEntropyNG - chaos_entropyNG;
				whole_clu_entropy     += entropy;
				whole_clu_entropy_NG  += entropyNG;

				tmpfreq.clear();
				clu_entropies[i] += entropy;
			}
			double val = - whole_clu_entropy / (double)res.clu.size();
			if (val<amin) amin = val;
			if (val>amax) amax = val;
			res.clustRealEnt[col]   = whole_clu_entropy;		// real entropy
			res.clustRealEntNG[col] = whole_clu_entropy_NG;
			res.clustEnt[col]     = val;						// normalized entropy
			total_entropy += whole_clu_entropy;
		}

		if (calc_ranks) {
			bool thesame = fabs(amax - amin) < DBL_EPSILON;
			for(int col=0; col<res.clustEnt.size(); col++) {
				int rank = 0;
				if (!thesame)
					rank = floor(99. - 99. * fabs(res.clustEnt[col] - amin) / fabs(amax - amin));
				res.clustEntRanks[col] = rank;
			}
		}
		for(int j=0; j<res.clu.size(); j++)
			res.clu[j].entropy = clu_entropies[j];

		res.whole_entropy = total_entropy;
	}
	void GetConservation() {

		float amin=FLT_MAX, amax=-FLT_MAX;
		_freqent.resize(_width, 0);
		_freqent2.resize(_width, 0);
		_freqentNG.resize(_width, 0);	// no gaps
		for(int col=0; col<_width; col++) {
			double func = 0, funcNG = 0.;
			for (ByteDblMapIter it=_colfreq[col].begin(); it!=_colfreq[col].end(); it++) {
			    byte sym = it->first;
			    double freq = it->second;
				double ff = freq * log(freq);
				func -= ff;
				if (sym==_gap) continue;
				funcNG -= ff;
			}
			if (amin>func) amin=func;
			if (amax<func) amax=func;
			_freqent[col] = func;
			_freqentNG[col] = funcNG;
			func = 0;
			for (ByteIntMapIter it=_colnsym[col].begin(); it!=_colnsym[col].end(); it++) {
			    byte sym = it->first;
			    double nsym = it->second;
				func -= nsym * log(nsym);
			}
			_freqent2[col] = func;
		}
	}

	// --------------------------------------------------------------------------------------------------

	void CacheColumnsFrequency() {
		if (!_width || !_height) return;
		if (_echoQ) cout << "caching columns frequencies ..." << endl;
		double height = _height;
		_colseqs.resize(_width, 0);
		_colnsym.resize(_width);
		_colfreq.resize(_width);
		ByteIntMapIter cit;
		ByteDblMapIter fit;
		for(int j=0; j<_width; j++) {
			ByteIntMap& nsym = _colnsym[j];
			ByteDblMap& freq = _colfreq[j];
			for(int i=0; i<_height; i++) {
				byte sym = at(j,i);
				if (sym!=_gap) _colseqs[j]++;
				cit = nsym.find(sym);
				if (cit==nsym.end()) nsym[sym] = 1; else nsym[sym]++;
			}
			for (cit=nsym.begin(); cit!=nsym.end(); cit++) {
				byte sym = cit->first;
			    int  cnt = cit->second;
			    double f = cnt / height;
			    freq[sym] = f;
			}
		}
	}
	void CacheChaosEntropy() {
		if (!_colfreq.size()) return;
		if (_echoQ) cout << "caching chaos entropies ..." << endl;
		_chaos_ent.resize(_width);
		ByteDblMapIter it;
		for(int col=0; col<_width; col++) {
			IntDblMap& hash = _chaos_ent[col];
			for (int j=2; j<=_height; j++) {
				ByteDblMap ww0;
				for(it=_colfreq[col].begin(); it!=_colfreq[col].end(); it++) {
					byte sym  = it->first;
				    double freq = it->second;
					ww0[sym] = freq * j;
				}
				double entropy = GetFreqHashEntropy(ww0, false);
				hash[j] = entropy;
			}
		}
	}
	double GetChaosEntropy(int col, int length) {
		IntDblMap& hash = _chaos_ent[col];
		IntDblMapIter it = hash.find(length);
		if (it==hash.end()) return 0;
		return it->second;
	}
	double GetFreqHashEntropy(ByteDblMap& hash, bool treat_as_int) {
		double entropy=0.;
		ByteDblMapIter it;
		for (it=hash.begin(); it!=hash.end(); it++) {
			byte sym  = it->first;
		    double freq = it->second;
			if (freq > 1.) entropy -= cached_gammln(freq + 1., treat_as_int);
		}
		return entropy;
	}
	double GetFreqHashEntropyNG(ByteDblMap& hash, bool treat_as_int, double& entropyNG) {
		double entropy=0.; entropyNG=0.;
		ByteDblMapIter it;
		for (it=hash.begin(); it!=hash.end(); it++) {
			byte sym  = it->first;
		    double freq = it->second;
			if (freq > 1.) {
				double v = cached_gammln(freq + 1., treat_as_int);
				entropy -= v;
				if (sym==_gap) continue;
				entropyNG -= v;
			}
		}
		return entropy;
	}
	double GetFreqHashEntropy(int* hist, int size) {
		double entropy = 0.;
		for(int k=0; k<size; k++) {
			if (hist[k] > 1)
				entropy -= cached_gammln(hist[k]+1, true);
		}
		return entropy;
	}
	double GetClustersEntropyDistance(CluResults& cr, int c1, int c2) {
		double entropy = 0.;
		for(int col=0; col<cr.columns.size(); col++) {
			entropy += GetClustersColumnEntropyDistance(cr, c1, c2, cr.columns[col]);
		}
		return entropy;
	}
	double GetClustersColumnEntropyDistance(CluResults& res, int c1, int c2, int col) {
		memset(_histogram, 0, _maxcharsize * sizeof(int));
		for(int i=0;i<res.clu[c1].ids.size();i++) {
			int seq = res.clu[c1].ids[i];
			if (_transposed[col][seq] == _symignore) continue;
			_histogram[ _transposed[col][seq] ]++;
		}
		for(int i=0;i<res.clu[c2].ids.size();i++) {
			int seq = res.clu[c2].ids[i];
			if (_transposed[col][seq] == _symignore) continue;
			_histogram[ _transposed[col][seq] ]++;
		}
		int both_clu_length  = res.clu[c1].ids.size() + res.clu[c2].ids.size();
		double chaos_entropy = GetChaosEntropy(col, both_clu_length);
		double columnEntropy = GetFreqHashEntropy(_histogram, _maxcharsize);
		double entropy = columnEntropy - chaos_entropy;
		return entropy;
	}

	bool WriteCluXml(string fn, string msafn, double elapsed, CluSettings &cs) {
		ostringstream text;
		if (!_cr.clu.size() || !_freqent.size() || !_cr.clustRealEnt.size()) return false;

		text << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl <<
				"<msa version=\"1\" rows=\"" << _height << "\" columns=\"" << _width << "\" elapsed=\"" << elapsed << "\">" << endl <<
				"  <source file=\"" << GetFileName(msafn) << "\"/>" << endl <<
				"  <conservation>" << ToString(_freqent) << "</conservation>" << endl <<
				"  <conservationEx>" << ToString(_freqent2) << "</conservationEx>" << endl <<
				"  <conservationNG>" << ToString(_freqentNG) << "</conservationNG>" << endl <<
				"  <clustering a1begin=\"" << cs.A1begin << "\" a1end=\"" << cs.A1end << "\" a1steps=\"" 
						<< cs.A1steps << "\" a1opt=\"" << _cr.A1 << "\" clusters=\"" << _cr.clu.size() << "\">" << endl <<
				"    <specificity>"    << ToString(_cr.clustRealEnt) << "</specificity>" << endl <<
				"    <specificityNG>"  << ToString(_cr.clustRealEntNG) << "</specificityNG>" << endl;
		for(int k=0; k<_cr.clu.size(); k++) {
		text << "    <cluster id=\"" << k << "\" size=\"" << _cr.clu[k].ids.size() << "\">" <<
				ToString(_cr.clu[k].ids) << "</cluster>" << endl;
		}
		text << "  </clustering>" << endl <<
				"</msa>" << endl;

		StringToFile(text, fn);
		return true;
	}
};
