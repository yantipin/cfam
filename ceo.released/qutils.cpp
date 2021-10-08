
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
#include <sstream>
#include <string>
#include <vector>
#include <map>

using namespace std;

#include "utils.h"
#include "bio_utils.h"
#include "msa.h"

struct ParseUniref {	// ---------------------------------------------------------------------------------------

	bool Parse(string fn) {
		int nseqs=0;
		string fno(StripExt(fn) + ".parsed.fa"), line;
		cout << "reading " << fn << " ... " << endl;
		ifstream fi(fn.c_str());
		ofstream fo(fno.c_str());
		if (!fi.good()) { cout << "Can't open " << fn << " for reading" << endl; return false; }
		if (!fo.good()) { cout << "Can't open " << fno << " for writing" << endl; return false; }
		int nline=0;
		while (getline(fi, line, '\n')) {
			nline++; if (line.size()<1) continue;
			if (line[0]!='>') { fo << line << endl; continue; }
			size_t found = line.find_last_of('=');
			if (found==string::npos) { cout << "= not found at line " << nline << endl; continue; }
			fo << '>' << line.substr(found+1) << endl;
			nseqs++;
		}
		cout << nseqs << " sequences written into " << fno << endl;
		fi.close();
		fo.close();
		return true;
	}
};

// measure average sequence identity in MSA for each sequence
struct FamilySequenceIdentity {
	MSA msa;
	bool Run(string fn) {
		cout << fn << endl;
		msa.Init(true, true);
		if (!msa.Read(fn)) return false;
		Elapsed ela(true);
		IntAry columns, seqs;

		double colthresh = 0.5;
		cout << "filtering columns with " << colthresh << " gaps ... ";
		msa.FilterColumns(colthresh, columns);
		cout << columns.size() << " left" << endl;

		double rowthresh = 0.5;
		cout << "filtering sequences with " << rowthresh << " gaps ... ";
		msa.FilterSequences(colthresh, columns, seqs);
		cout << seqs.size() << " left" << endl;

		cout << "getting all the pairs identity ..." << endl;
		msa.GetAllTheSeqsPairsIdentity(columns, seqs);

		cout << "elapsed : " << ela.get() << endl;

		return true;
	}
};

// -----------------------------------------------------------------------------------------------
//
// params: <command> <filename>
//         uniref-fasta /home/user/uniref90.fa
//         family-ident /home/user/pfam/PF00089.txt

int main(int argc, char** argv) {

	if (argc<3) return 0;
	string cmd(argv[1]);
	string fn(argv[2]);

	if (cmd==string("uniref-fasta")) {
		ParseUniref pu;
		pu.Parse(fn);
	} else if (cmd==string("family-ident")) {
		FamilySequenceIdentity fsi;
		fsi.Run(fn);
	}

    return 0;
}












