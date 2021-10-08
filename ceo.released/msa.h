
#pragma once

// requires utils.h
// requires bio_utils.h

struct MSA {	// ---------------------------------------------------------------------- MSA ------------------------------------

	StrAry _seqs;
	StrAry _names;
	vector<ProtRgnName> _rgns;
	int _maxnamelen;
	int _namepad;
	int _annotpad;
	int _seqlen;
	byte gap;
	int	_width;
	int	_height;
	int _gapcount;			// total number of gaps
	bool _echoQ;
	bool _forceuppercaseQ;
	bool _replacedotsQ;		// replace dots with gaps (in Pfam raw MSA)
	IntAry _identity;
	Sequence _rf;			// HMM msa does not have ref.seq in it! have to read it separately and align to hmm "#=GC RF" row
	ParseHmm _hmm;
	bool _hmmQ;
	bool _posmapQ;			// create/maintain/save protein position map?
	IntIntAry _posmap;		// requires valid domain-region ID (EGFR_HUMAN/145-345), sequence index -> list of columns -> absolute protein position

	MSA(bool echoQ=true, bool uppercaseQ=false, bool replacedotsQ=false, bool posmapQ=false) : _echoQ(echoQ), _hmmQ(false), _posmapQ(posmapQ) {
		_maxnamelen = _seqlen = _gapcount = 0; gap = '-'; Init(uppercaseQ, replacedotsQ);
	}
	void Init(bool uppercaseQ=false, bool replacedotsQ=false) {
		_replacedotsQ = replacedotsQ;
		_forceuppercaseQ = uppercaseQ;
		_maxnamelen = _seqlen = 0;
		_namepad = 6;				// space between names and sequences
		_annotpad = 3;
		gap = '-';
	}
	bool GetRefSeqRegion(int &from, int &to) {
		bool Q = _rgns.size() && _rgns[0].validQ();
		from = Q ? _rgns[0].from : 0;
		to = Q ? _rgns[0].to : 0;
		return Q;
	}
	inline byte at(int column, int row) { return _seqs[row][column]; }
	void ReadRefSeq(string fn) {
		string seq;
		FileToString(fn, seq);
		_rf.FromFasta(seq);
	}
	// annot is text annotation per each sequence, printed as a column between seq.name and msa
	// annot must be parallel to rows and has same length
	bool Write(string fn, IntAry& columns, IntAry& rows, StrAry& annot) {
		IntAry seqs;
		if (!rows.size()) FillInc(seqs, _height); else seqs = rows;
		bool allQ = columns.size()==0;
		bool annotQ = annot.size()!=0;
		int annotW = 0;
		if (annotQ) {
			if (annot.size()!=seqs.size()) return false;
			for(int k=0; k<annot.size(); k++) if (annot[k].size()>annotW) annotW = annot[k].size();
			annotW += _annotpad;
		}
		ofstream f(fn.c_str());
		if (!f.good()) { cout << endl << "Can't open " << fn << " for writing" << endl; return false; }
		int pad = _maxnamelen + annotW + _namepad;

		for(int k=0; k<seqs.size(); k++) {
			int row = seqs[k];
			string s(padright(_names[row], pad));
			if (annotQ) s += padleft(annot[k], annotW) + padleft(" ", 3);
			f.write(s.c_str(), s.size());
			if (allQ) {
				f.write(_seqs[row].c_str(), _seqs[row].size());
			} else {
				for(int c=0; c<columns.size(); c++)
					f.write(&(_seqs[row][ columns[c] ]), 1);
			}
			f.write(nl, 1);
		}
		f.close();
		return true;
	}
	// ------------------------------------------------------------------------------------------------------------------
	// this will read STOCKHOLM msa format (hmm output), as well as plain text msa
	//
	bool Read(string fn, int limit=0) {
		_seqs.clear();
		_names.clear();
		_rgns.clear();
		_gapcount = 0;
		string line, name;
		if (_echoQ)	cout << "reading " << fn << (_forceuppercaseQ ? " (force upper case)" : "")
						 << (_replacedotsQ ? " (replace . with -)" : "") << " ... " << endl;
		FileLinesIter it;
		if (!it.open(fn)) return false;
		string hmmrca("#=GC RF");	// look for HMM "reference coordinate annotation" line to align reference sequence to

		while(it.getline()) {
			int l = it._line.length();
			if (l<2) continue;
			if (l>8 && it._line[0]=='#' && hmmrca==it._line.substr(0, 7)) { InsertRefSeqUsingHmmRCA(it._line, it._nline); continue; }
			if (limit && _seqs.size()==limit) continue;				// have to read file to the end because 'GC RF' is very last line
			if (it._line[0]=='#' || it._line[0]=='/') continue;

			if (limit && _seqs.size()==limit) continue;		// limit reached but keep reading to parse "GC RF" last line
			if (!AddSequence(it._line, it._nline)) break;
		}
		_height = _seqs.size();
		if (!_height) { cout << "no sequences found" << endl; return false; }
		if (_echoQ) {
			int size = _seqs.size() * _seqs[0].length();
			double dens = 100 * _gapcount / size;
			cout << "msa dimensions : " << _seqs.size() << " seq x " << _seqs[0].length() << " res" << endl
				 << "fraction of gaps : " << std::fixed << std::setprecision(1) << dens << "%" << endl;
		}
		it.close();
		_width = _seqlen;
		return true;
	}
	bool AddSequence(string& line, int nline, bool refseqQ=false) {
		size_t found = line.find_first_of(' ');
		if (found==string::npos) {
			cout << "space not found at line " << nline << endl << line << endl;
			return false;
		}
		string name = line.substr(0, found);
		if (!refseqQ && name==_rf.name) return true;	// do not read reference sequence, it will be added separately to make sure it's first one
		// if hmm domain hit table is loaded, leave only "good" domains in msa, unless this is reference sequence
		if (!refseqQ && _hmmQ && !_hmm.BestDomainQ(name)) return true;

		if (!refseqQ) {
			_names.push_back(name);
			_rgns.push_back(ProtRgnName(name));
		} else {
			_names.insert(_names.begin(), name);
			_rgns.insert(_rgns.begin(), ProtRgnName(name));
		}
		int len = _names.back().size();
		if (len > _maxnamelen) _maxnamelen = len;
		found = line.find_first_not_of(' ', found);
		if (found==string::npos) { cout << "sequence not found in " << endl << line << endl; return false; }

		string seq(line.substr(found));
		if (_forceuppercaseQ)	// uppercase all sequences
			std::transform(seq.begin(),seq.end(),seq.begin(), (int(*)(int))toupper);	// pretty unbelievable that crap like this is now a c++ standard

		for(int i=0; i<seq.size(); i++) {
			if (_replacedotsQ && seq[i]=='.') seq[i]=gap;
			if (seq[i]==gap) _gapcount++;
		}

		if (!refseqQ) _seqs.push_back(seq); else _seqs.insert(_seqs.begin(), seq);

		int last = _seqs.size()-1;
		if (!_seqlen) _seqlen = _seqs[0].size();
		if (_seqlen!=_seqs[last].size()) {
			cout << "Sequence at line " << nline << " (length=" << _seqs[last].size() << ") : "
					<< endl << _seqs[last] << endl
					<< "previous (length=" << _seqs[last-1].size() << ") : "
					<< _seqs[last-1] << endl;
			return false;
		}
		return true;
	}
	// this assumes reference sequence is a first one
	void RemoveRefSeqGaps(IntAry& columns) {
		if (!_seqs.size()) return;
		columns.clear();
		for (int c=0; c<_width; c++) {
			if (_seqs[0][c]==gap) continue;
			columns.push_back(c);
		}
	}
	// threshold is a fraction [0..1] of gaps greater than given to be removed
	void FilterColumnsByGaps(double threshold, IntAry& columns) {
		if (!(threshold>0)) return;
		columns.clear();
		int thresh = _height * threshold;
		for (int c=0; c<_width; c++) {
			int gaps = 0;
			for (int r=0; r<_height; r++) {
				if (gap==at(c,r)) gaps++; else continue;
				if (gaps>thresh) break;
			}
			if (gaps>thresh) continue;
			columns.push_back(c);
		}
	}
	int GetMSAColumnByRefseqProteinPosition(int pos) {
		if (!_rgns.size()) return -1;			// have to know domain boundaries
		return pos - _rgns[0].from;
	}
	// !!! this function assumes indicies of non-gappy positions in reference sequences are stores in "columns" array !!!
	// position is absolute protein position
	// relative domain boundaries are extracted from reference sequence ID --> "ANDR_HUMAN/687-894"
	//
	bool RemoveSeqsWithGapsInPosition(IntAry& columns, IntAry& seqs, int position) {
		if (!_rgns.size()) return false;			// have to know domain boundaries
		int column = position - _rgns[0].from;
		if (column<0 || column>=columns.size()) return false;
		IntAry rows;
		if (!seqs.size()) FillInc(rows, _height); else rows = seqs;
		seqs.clear();	// output
		int height = rows.size();
		for (int r=0; r<height; r++) {
			int row = rows[r];
			if (gap==at(columns[column],row)) continue;
			seqs.push_back(row);
		}
		return true;
	}
	// threshold is a fraction [0..1] of gaps greater than given to be removed using 'columns'
	void FilterSequencesByGaps(double threshold, IntAry& columns, IntAry& seqs) {
		if (!(threshold>0)) return;
		IntAry cols, rows;
		if (!columns.size()) FillInc(cols, _width); else cols = columns;
		if (!seqs.size()) FillInc(rows, _height); else rows = seqs;
		seqs.clear();	// output
		int height = rows.size();
		int width = cols.size();
		int thresh = width * threshold;
		for (int r=0; r<height; r++) {
			int row = rows[r];
			int gaps=0, col=0;
			for (int i=0; i<width; i++) {
				if (gap==at(cols[i],row)) gaps++; else continue;
				if (gaps>thresh) break;
			}
			if (gaps>thresh) continue;
			seqs.push_back(row);
		}
	}
	// threshold is a fraction [0..1], remove second sequence if identity is greater than given
	// identity is measured in 'columns' only, pass empty array to use all the columns
	void FilterSequencesByIdentity(double threshold, IntAry& columns, IntAry& seqs) {
		if (!(threshold>0)) return;
		IntAry cols, rows;
		if (!columns.size()) FillInc(cols, _width); else cols = columns;
		if (!seqs.size()) FillInc(rows, _height); else rows = seqs;
		int height = rows.size();
		int width = cols.size();
		IntAry removedQ;
		removedQ.resize(rows.size(), 0);
		for (int r1=0; r1<height-1; r1++) {
			int row1 = rows[r1];
			if (removedQ[r1]) continue;
			for (int r2=r1+1; r2<height; r2++) {
				int row2 = rows[r2];
				if (removedQ[r2]) continue;
				int ident=0, len1=0, len2=0;
				for (int i=0; i<width; i++) {
					bool gapQ = false;
					if (gap==at(cols[i],row1)) gapQ = true; else len1++;
					if (gap==at(cols[i],row2)) gapQ = true; else len2++;
					if (gapQ) continue;
					if (at(cols[i],row1)==at(cols[i],row2)) ident++;
				}
				int which = r2;	double len = len2;
				if (len2 > len1) { len = len1; which = r1; }			// shorter one is removed
				double lthresh = ident / len;
				if (which && lthresh>threshold) removedQ[which] = 1;	// very first sequence is considered reference one, don't remove it
			}
		}
		seqs.clear();	// output
		for(int k=0; k<rows.size(); k++) {
			if (removedQ[k]) continue;
			seqs.push_back(rows[k]);
		}
	}
	// !!! this function assumes that there are absolutely no gaps in reference sequence !!!
	// computes two types of sequence identity - relative to msa width and relative to each sequence length
	//
	void GetIdentityToRefSeq(IntAry &rows, IntAry &columns, DblAry &ident, DblAry &identLen) {
		IntAry seqs, cols;
		if (!columns.size()) FillInc(cols, _width); else cols = columns;
		if (!rows.size()) FillInc(seqs, _height); else seqs = rows;
		int height = seqs.size();
		int width = cols.size();
		ident.resize(height, 0); ident[0] = 1;
		identLen.resize(height, 0); identLen[0] = 1;
		for (int r=1; r<height; r++) {
			int len=0, idt=0;
			for (int c=0; c<width; c++) {
				if (gap==at(cols[c],seqs[r])) continue;
				len++;
				if (at(cols[c],seqs[r])==at(cols[c],0)) idt++;
			}
			ident[r] = idt/(double)width;
			identLen[r] = idt/(double)len;
		}
	}

	// ------------------------------------------------------------------------------------------ hmm --------------------

	bool InsertRefSeqUsingHmmRCA(string& line, int nline) {
		if (!_rf.seq.size()) return false;	// need reference sequence
		size_t pos = line.find_first_not_of(' ', 7);	// "#=GC RF"
		if (pos==string::npos) return false;
		string seq(line.substr(pos));
		int rfpos = 0;
		string refseq;
		for(pos=0;pos<seq.size();pos++) {
			if (seq[pos]=='.') { refseq += '-'; continue; }
			refseq += _rf.seq[rfpos++];
		}
		string row(_rf.name + "  " + refseq);
		AddSequence(row, nline, true);
		return true;
	}
	bool ReadHmmDomainHits(string fn) { return _hmmQ = _hmm.ParseDomainHits(fn); }

	bool ReorderByHmmDomainCValue(IntAry& rows) {	// rows must be empty or start with 0 -- reference sequence assumed first one
		if (!rows.size()) FillInc(rows, _height);
		if (!_hmm._cvalues.size()) return false;
		IntAry seqs; DblAry cvals;
		if (_rf.seq.size()) { seqs.push_back(0); cvals.push_back(0); }	// so reference sequence stays first one
		double cval;
		for (int k=1; k<rows.size(); k++) {
			if (!_hmm.GetDomainCValue(_names[rows[k]], cval)) continue;
			cvals.push_back(cval);
			seqs.push_back(rows[k]);
		}
		OrderedSort<double> os(cvals, seqs); // seqs will be sorted in parallel to cvals
		rows = seqs;
		return true;
	}

#if 0

	// since number of pairs can be larger than number of bytes in your RAM we gonna cache just some of them
	struct CachedPairsIdent {
		size_t _size;
		size_t _limit;
		IntAry& _ident;
		CachedPairsIdent(IntAry& ident, int size) : _ident(ident) {
			ullong maxi = size-2, maxj = size-1;
			_limit = 268,435,455; // this many integers in 1GB == 1,073,741,823 / 4 bytes
			ullong maxidx = maxi + maxj * (maxj - 1) / 2;
			_size = maxidx < _limit ? maxidx : _limit;
			_ident.resize(_size, 0);
		}
	};

	void GetAllTheSeqsPairsIdentity(IntAry& columns, IntAry& seqs) {

		IntAry cols, rows;
		if (!columns.size()) FillInc(cols, _width); else cols = columns;
		if (!seqs.size()) FillInc(rows, _height); else rows = seqs;
		int width = cols.size();
		int height = rows.size();
		int h = height-1;
		ullong total = PairIndex(height-2, height-1);

		CachedPairsIdent cpi(_identity, _height);

		for(int i=0; i<height-1; i++) {
			for(int j=i+1; j<height; j++) {
				int ni=0, col, ri, rj;
				for(int p=0; p<width; p++) {
					col = cols[p]; ri = rows[i]; rj = rows[j];
					if ('-'!=_seqs[ri][col] && _seqs[ri][col]==_seqs[rj][col]) ni++;
				}
				int idx = GetPairIndex(i,j);
				_identity[idx] = ni / width;
			}
			cout << i << endl;
		}
	}

#endif

};
