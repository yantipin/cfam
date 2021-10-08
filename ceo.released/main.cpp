
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <memory.h>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

using namespace std;

#include "utils.h"
#include "clutypes.h"
#include "bio_utils.h"
#include "msa.h"
#include "bmatrix.h"

// cmd=clu msa=C:\cfam\data\msa\method.tests\msa.100.txt cfrom=0.65 cto=0.95 csteps=6
// msa=C:\cfam\data\msa\hmmer\ident.msa.txt fs=0.5
// cmd=msa quiet=1 limit=1000 fs=0.1 fi=0.95 refseq=C:\cfam\data\msa\hmmer\EGFR.300.fa msa=C:\cfam\data\msa\hmmer\EGFR.300.msa hmmd=C:\cfam\data\msa\hmmer\EGFR.300.domains cfrom=0.65 cto=0.95 csteps=6
// cmd=variants quiet=1 limit=700 fs=0 fi=0 refseq=C:\cfam\data\msa\hmmer\EGFR_HUMAN.712-968.fa msa=C:\cfam\data\msa\hmmer\EGFR_HUMAN.712-968.msa hmmd=C:\cfam\data\msa\hmmer\EGFR_HUMAN.712-968.domains cfrom=0.65 cto=0.95 csteps=6
// cmd=msa quiet=1 limit=700 fs=0.7 fi=0.98 refseq=C:\cfam\db\external\IARC.TP53\MA\MA.hmm\P53_HUMAN.75-309.fa msa=C:\cfam\db\external\IARC.TP53\MA\MA.hmm\P53_HUMAN.75-309.msa hmmd=C:\cfam\db\external\IARC.TP53\MA\MA.hmm\P53_HUMAN.75-309.domains cfrom=0.65 cto=0.95 csteps=6
// C:\cfam\data\msa\CEA\
// cmd=msa quiet=1 limit=700 fs=0.7 fi=0.98 refseq=C:\cfam\db\external\IARC.TP53\MA\MA.hmm\P53_HUMAN.75-309.fa msa=C:\cfam\db\external\IARC.TP53\MA\MA.hmm\P53_HUMAN.75-309.msa hmmd=C:\cfam\db\external\IARC.TP53\MA\MA.hmm\P53_HUMAN.75-309.domains cfrom=0.65 cto=0.95 csteps=6
// limit=700 fs=0.7 fi=0.98 refseq=C:\cfam\data\msa\CEA\seq.fasta.txt msa=C:\cfam\data\msa\CEA\hmm.msa.txt hmmd=C:\cfam\data\msa\CEA\hmm.domains.txt msaout=C:\cfam\data\msa\CEA\msa.flt.txt cfrom=0.65 cto=0.95 csteps=6
//
// cmd=msa fs=0.7 fi=0.98 msa=C:\cfam\data\msa\CEA\msa.mall.txt msaout=C:\cfam\data\msa\CEA\msa.flt.txt msacluout=C:\cfam\data\msa\CEA\msa.clu.flt.txt msahtmlout=C:\cfam\data\msa\CEA\msa.flt.html cfrom=0.65 cto=0.95 csteps=6
//
// --------------------------------------------------------------------------------------------- MsaMatrix -------------------
//

struct MsaMatrix : public BMatrix {
	ostringstream _ss;
	ByteIntMapIter _bit;

	bool WriteFISofAllVariants(string fn, int domain_start=1) {
		if (domain_start<1) domain_start = 1;	// so 0 can be passed safely when region is unknown
		Elapsed ela; ela.start();
		ByteByteMapIter it;

		// we need scores even if msa consists of only one sequence
		_ss.str(""); _ss << "#variant" << TAB << "FIS" << TAB << "vcons" << TAB << "vspec" << TAB << "gaps" << endl;
		if (_height==1) {
			for(int p=0; p<_width; p++) {
				byte AAfrom = at(p, 0);				// code, not actual character
				for(int res=0; res<lenAA; res++) {
					byte symAAto = AA[res];
					it = _sym2code.find(symAAto);
					byte AAto = it==_sym2code.end() ? SYM_NOT_FOUND : it->second;
					if (AAto==AAfrom) continue;
					_ss << _code2sym[AAfrom] << (domain_start + p) << symAAto << TAB;
					if (AAfrom==_gap) _ss << TAB << TAB << TAB;
					else _ss << 0 << TAB << 0 << TAB << 0 << TAB << 0;
					_ss << endl;
				}
			}
			StringToFile(_ss, fn);
			return true;
		}
		if (!_cr.clu.size()) return false;
		double FIS=0, vc=0, vs=0;
		FillRowToClu();
		ByteIntMap clusymcnt;
		int cluidx = _row2clu[0];				// cluster of reference sequence

		for(int p=0; p<_width; p++) {
			byte AAfrom = at(p, 0);				// code, not actual character
			GetColumnCluSymCount(p, cluidx, clusymcnt);
			double fgaps = (_height-_colseqs[p]) / (double)_height;

			for(int res=0; res<lenAA; res++) {
				byte symAAto = AA[res];
				it = _sym2code.find(symAAto);
				byte AAto = it==_sym2code.end() ? SYM_NOT_FOUND : it->second;
				if (AAto==AAfrom) continue;
				bool ok = GetFIS(clusymcnt, AAfrom, p, AAto, FIS, vc, vs);
				_ss << _code2sym[AAfrom] << (domain_start + p) << symAAto << TAB;
				if (ok) _ss << FIS << TAB << vc << TAB << vs << TAB << fgaps; else _ss << TAB << TAB << TAB;
				_ss << endl;
			}
		}
		if (_echoQ) cout << "FIS for all variants computed in : " << ela.get() << " secs." << endl;
		StringToFile(_ss, fn);
		return true;
	}
	// msacolumn is 0-based
	bool GetFIS(ByteIntMap& clusymcnt, byte AAfrom,
				int msacolumn, byte AAto, double &FIS, double &vc, double &vs) {
		FIS = vc = vs = 0;
		if (AAfrom==_gap) return false;

		// variant conservation score
		_bit = _colnsym[msacolumn].find(AAto);
		double dt = _bit==_colnsym[msacolumn].end() ? 0 : _bit->second;
		vc = log( _colnsym[msacolumn][AAfrom] / (dt + 1) );

		// variant specificity score
		_bit = clusymcnt.find(AAto);
		double df = _bit==clusymcnt.end() ? 0 : _bit->second;
		vs = -log( (df + 1) / clusymcnt[AAfrom] );

		FIS = 0.5 * (vc + vs);

		return true;
	}
	// msacolumn is 0-based
	bool GetFIS2(ByteIntMap& clusymcnt, byte AAfrom,
				int msacolumn, byte AAto, double &FIS, double &vc, double &vs) {
		FIS = vc = vs = 0;
		if (AAfrom==_gap) return false;

		// variant conservation score
		_bit = _colnsym[msacolumn].find(AAto);
		double dt = _bit==_colnsym[msacolumn].end() ? 0 : _bit->second;
		double nonAAto = 0;
		for (ByteIntMapIter it=_colnsym[msacolumn].begin(); it!=_colnsym[msacolumn].end(); it++) {
			byte sym = it->first;
			if (sym==AAto) continue;
			nonAAto += it->second;
		}
		vc = log( nonAAto / (dt + 1) );

		// variant specificity score
		_bit = clusymcnt.find(AAto);
		double df = _bit==clusymcnt.end() ? 0 : _bit->second;
		nonAAto = 0;
		for (ByteIntMapIter it=clusymcnt.begin(); it!=clusymcnt.end(); it++) {
			byte sym = it->first;
			if (sym==AAto) continue;
			nonAAto += it->second;
		}
		vs = log( nonAAto / (df + 1) );

		FIS = 0.5 * (vc + vs);

		return true;
	}

	void GetColumnCluSymCount(int column, int cluidx, ByteIntMap& symcnt) {
		symcnt.clear();
		for(int k=0; k<_cr.clu[cluidx].ids.size(); k++) {
			int row = _cr.clu[cluidx].ids[k];
			byte sym = at(column, row);
			if (sym==_gap) continue;
			_bit = symcnt.find(sym);
			if (_bit==symcnt.end()) symcnt[sym] = 1; else symcnt[sym]++;
		}
	}
	void WriteMsaWithClustering(MSA &msa, string fn) {
		if (msa._height==1) return;
		IntAry columns, rows;
		StrAry cluindecies;
		ostringstream s;
		for(int ci=0; ci<_cr.clu.size(); ci++) {
			s.str(""); s << ci;
			for(int r=0; r<_cr.clu[ci].ids.size(); r++) {
				rows.push_back(_cr.clu[ci].ids[r]);
				cluindecies.push_back(s.str());
			}
		}
		msa.Write(fn, columns, rows, cluindecies);
	}
	void AddMapSyms(ByteStrMap& m, string str, const char* syms) {
		int len = strlen(syms);
		for(int k=0; k<len; k++) m[syms[k]] = str;
	}
	// column is optional column index to highlight (given in 0-based msa column index)
	void WriteMsaHtml(MSA &msa, string fn, int column=-1) {
		if (msa._height==1) return;
		IntAry columns, rows;
		DblAry ident, identLen;
		msa.GetIdentityToRefSeq(rows, columns, ident, identLen);
		bool identQ = ident.size()==msa._height;
		StrAry cluindecies;
		ostringstream s;
		ByteStrMap mview;
		ByteStrMapIter mvit;
		AddMapSyms(mview, "gr", ".-");
		AddMapSyms(mview, "c0", "GAIVLM");
		AddMapSyms(mview, "c1", "FYWH");
		AddMapSyms(mview, "c2", "C");
		AddMapSyms(mview, "c3", "P");
		AddMapSyms(mview, "c4", "KR");
		AddMapSyms(mview, "c5", "DE");
		AddMapSyms(mview, "c6", "QN");
		AddMapSyms(mview, "c7", "ST");
		AddMapSyms(mview, "c8", "BZX");
		AddMapSyms(mview, "c9", "?*");
		AddMapSyms(mview, "bl", "U");
		s << "<html>" << nl
		  << "<style type='text/css'>" << nl
		  << ".msa { font-family:Lucida Console,monospace; font-size:16px; line-height:95%; }" << nl
		  << ".hl  { background-color:#e0e0e0;  }" << nl
		  << ".gr  { color:grey;    }" << nl
		  << ".bl  { color:black;   }" << nl
		  << ".m0  { color:blue;    }" << nl
		  << ".m1  { color:red;     }" << nl
		  << ".c0  { color:#33cc00; }" << nl	// mview residue coloring
		  << ".c1  { color:#009900; }" << nl
		  << ".c2  { color:#ffff00; }" << nl
		  << ".c3  { color:#33cc00; }" << nl
		  << ".c4  { color:#cc0000; }" << nl
		  << ".c5  { color:#0033ff; }" << nl
		  << ".c6  { color:#6600cc; }" << nl
		  << ".c7  { color:#0099ff; }" << nl
		  << ".c8  { color:#666666; }" << nl
		  << ".c9  { color:#999999; }" << nl
		  << "</style>" << nl
		  << "<body>" << nl
		  << "<table border=0 cellpadding=0 cellspacing=0>" << nl;
		for(int ci=0; ci<_cr.clu.size(); ci++) {
			for(int r=0; r<_cr.clu[ci].ids.size(); r++) {
				int row = _cr.clu[ci].ids[r];
				ostringstream idt1, idt2;
				if (identQ) {
					idt1 << std::fixed << std::setprecision(1) << 100*ident[row] << '%';
					idt2 << std::fixed << std::setprecision(1) << 100*identLen[row] << '%';
				}
				s << "<tr><td class='msa' nowrap>" << msa._names[row] << "</td><td width='20' class='msa' nowrap>&nbsp;</td>"
				  << "<td class='msa' align='right' nowrap>" << idt1.str() << "</td><td width='20' class='msa' nowrap>&nbsp;</td>"
				  << "<td class='msa' align='right' nowrap>" << idt2.str() << "</td><td width='20' class='msa' nowrap>&nbsp;</td>"
				  << "<td class='msa' align='center' nowrap>" << (ci+1) << "</td><td width='20' class='msa' nowrap>&nbsp;</td>"
				  << "<td class='msa' nowrap>";
				for(int c=0; c<msa._seqs[row].size(); c++) {
					char sym = msa._seqs[row][c];
					mvit = mview.find(sym);
					bool colQ = c==column;
					if (mvit==mview.end()) {
						if (colQ) s << "<span class='hl'>";	// class name to highlight column
						s << sym;
						if (colQ) s << "</span>";
					}
					else {
						s << "<span class='" << mvit->second;
						if (colQ) s << " hl";	// class name to highlight column
						s << "'>" << sym << "</span>";
					}
				}
				s << "</td></tr>" << nl;
			}
			if (ci!=_cr.clu.size()-1)
				s << "<tr><td height='2' colspan='5' bgcolor='#E0E0E0'></td></tr>" << nl;
		}
		s << "</table>" << nl
		  << "</body></html>" << nl;
		StringToFile(s, fn);
	}
};

// ------------------------------------------------------------------------------------------------------- main ------------
//
// params :

int main(int argc, char** argv) {

//	SanityCheck(); return 0;

	cout << "Combinatorial Entropy Optimization (CEO) for multiple sequence alignment, v.1" << endl
		 << "Please cite: Reva, B.A., Antipin, Y.A. and Sander, C. (2007) Genome Biol, 8, R232." << endl
		 << "             Determinants of protein function revealed by combinatorial entropy optimization" << endl << endl;

	if (argc<2) {
		cout << "usage : msa=<filename> msaout=<filename> cfrom=<float> cto=<float> csteps=<integer> limit=<integer> quiet=<0/1>" << endl
			 << "        cmd - clu: perform clustering, otherwise just write filtered MSA" << endl
			 << "        msa - multiple sequence alignment input file in STOCKHOLM (Pfam/Hmmer) or plain text format" << endl
			 << "        cfrom/cto - optimization region [0..1], defaults are 0.75..0.75, recommended 0.65..0.95" << endl
			 << "        csteps - number of iterations over opt.region, default is 1, recommended 6" << endl
			 << "        limit - number of sequences to read counting from the top of msa" << endl
			 << "        msaout - name of output msa (if processed/filtered)" << endl
			 << "        fc - filter columns (incompatible with 'refseq'), fraction of gaps to remove" << endl
			 << "        fs - filter sequences, fraction of gaps to remove" << endl
			 << "        fi - filter sequences by identity, all pairs are considered, second sequence is removed" << endl
			 << "        hmmd - optional Hmmer domain hits table to read q-values from and sort MSA accordingly," << endl
			 << "               only most significant domain hit per protein is included in MSA" << endl
			 << "        refseq - optional reference sequence in FASTA file (since it's not in Hmmer output MSA)" << endl
			 << "                 to make it first one, columns with gaps in refseq are removed" << endl
			 << "        msahtmlout - produce MSA HTML with mview coloring and subfamilies" << endl
			 << endl;
		return 0;
	}

	CmdArgs ca(argc, argv);		// -------------------------------------------- arguments ----------------------------------

	bool ok;
	string msafn, rffn, rfgaps, hmmd, cmd, msaout, msacluout, msahtmlout;
	int steps, quiet, limit, position;
	double start, end, fc, fs, fi;

	ok = ca.GetValue("msaout", string(""), msaout);
	bool msacluQ = ca.GetValue("msacluout", string(""), msacluout);
	bool msahtmlQ = ca.GetValue("msahtmlout", string(""), msahtmlout);

	ok = ca.GetValue("msa", string(""), msafn);
	if (!ok) { cout << "filename not specified" << endl; return 1; }

	string cmdClu("clu"), cmdVars("vars");

	ok = ca.GetValue("cmd", string(""), cmd);	// command: 'msa' or no command - read/write msa, 'clu' - clustering, 'vars - clustering + write FIS for all variants
	ok = ca.GetValue("cfrom", 0.75, start);		// clustering: start of optimization parameter A1
	ok = ca.GetValue("cto", 0.75, end);			// clustering: end of optimization parameter A1
	ok = ca.GetValue("csteps", 1, steps);		// clustering: number of optimization steps to take over A1 region
	ok = ca.GetValue("quiet", 0, quiet);		// minimal console output
	ok = ca.GetValue("limit", 0, limit);		// limit number of sequences to read from msa
	ok = ca.GetValue("fc", 0., fc);				// filter columns by gaps
	ok = ca.GetValue("fs", 0., fs);				// filter sequences by gaps
	ok = ca.GetValue("fi", 0., fi);				// filter sequences by identity
	ok = ca.GetValue("gapspos", 0, position);	// specify protein position -- all sequences with gaps in this position are removed -- which makes specificity analysis position-specific!

	bool hmmdQ = ca.GetValue("hmmd", string(""), hmmd);
	bool rfQ = ca.GetValue("refseq", string(""), rffn);

	ca.GetValue("refseqgaps", string(""), rfgaps);
	bool rfgapsQ = rfgaps=="1";

	bool echoQ = quiet==0;

	MSA msa(echoQ, true, true);	// ------------------------------------------- MSA ----------------------------------------
	if (rfQ) msa.ReadRefSeq(rffn);
	if (hmmdQ) {
		cout << "reading hmm domain hits " << hmmd << endl;
		msa.ReadHmmDomainHits(hmmd);
	}

	if (!msa.Read(msafn, limit)) return 1;

	IntAry columns, rows;
	if (rfQ || rfgapsQ) msa.RemoveRefSeqGaps(columns);
	if (position) msa.RemoveSeqsWithGapsInPosition(columns, rows, position);

	if (!rfQ && fc > 0) msa.FilterColumnsByGaps(fc, columns);		// column filtering only if reference sequence not set
	if (fs > 0) msa.FilterSequencesByGaps(fs, columns, rows);
	if (fi > 0) msa.FilterSequencesByIdentity(fi, columns, rows);

	if (hmmdQ) msa.ReorderByHmmDomainCValue(rows);

	StrAry annot;
	string rfmsafn = msaout.size() ? msaout : StripExt(msafn) + ".flt.msa";
	msa.Write(rfmsafn, columns, rows, annot);

	if (cmd!=cmdClu && cmd!=cmdVars) return 0;

	// ---------------------------------------------------------------------- clustering -----------------------------------

	cout << "clustering : [" << start << ".." << end << "] -- " << steps << " step(s)" << endl;

	MSA msa2;
	if (!msa2.Read(rfmsafn)) return 1;

	MsaMatrix m;
	if (!m.FromMSA(msa2)) return 1;

	InitGammaLn();
	CluSettings cs(start, end, steps);
	m.Clusterize(rfmsafn, cs, echoQ);

	if (msacluQ) m.WriteMsaWithClustering(msa2, msacluout);

	int column = -1;
	if (position) column = msa2.GetMSAColumnByRefseqProteinPosition(position);
	if (msahtmlQ) m.WriteMsaHtml(msa2, msahtmlout, column);

	if (cmd!=cmdVars) return 0;

	// ------------------------------------------------------------ all possible variants in reference sequence ------------

	string fn(StripExt(rfmsafn) + ".variants");
	int domain_start, domain_end;
	msa2.GetRefSeqRegion(domain_start, domain_end);
	m.WriteFISofAllVariants(fn, domain_start);

    return 0;
}
