/*
Szymon Rusinkiewicz
Princeton University

mesh_cc.cc
Determine the connected components of a mesh, and possibly write
out only selected components of the object.

Does the same thing as "plycomps", part of the plytools package by Greg Turk.
*/

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "TriMesh.h"
#include "TriMesh_algo.h"
#ifdef _WIN32
# include "wingetopt.h"
#else
# include <unistd.h>
#endif
#include <cstdio>
using namespace std;
using namespace trimesh;


#define BIGNUM numeric_limits<int>::max()


// Print out the connected components that are bigger than morethan and
// smaller than lessthan.  The largest min(nprint, total) components are
// printed out, unless morethan == 0 and lessthan != BIGNUM, in which case
// the smallest min(nprint, total) components are printed
void print_comps(const vector<int> &compsizes,
                 int morethan, int lessthan, int total, int nprint)
{
	printf("%lu connected components total.\n",
		(unsigned long) compsizes.size());

	if (compsizes.size() < 1 || total < 1)
		return;
	int numtoprint = min(min(nprint, total), (int)compsizes.size());

	if (morethan == 0 && lessthan != BIGNUM) {
		// Print numtoprint smallest components
		for (int i = 0; i < numtoprint; i++)
			printf(" Component #%d - %d faces\n", i+1, compsizes[compsizes.size()-1-i]);
	} else {
		// Print numtoprint largest components
		for (int i = 0; i < numtoprint; i++)
			printf(" Component #%d - %d faces\n", i+1, compsizes[i]);
	}
	if (numtoprint != (int)compsizes.size())
		printf(" ...\n");
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [options] in.ply [out.ply]\n", myname);
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "	-v		Base connectivity on vertices, not edges\n");
	fprintf(stderr, "    	-a		Print sizes of ALL the components\n");
	fprintf(stderr, "    	-s		Split components into separate files\n");
	fprintf(stderr, "    	-m n		Select components with >= n faces\n");
	fprintf(stderr, "    	-l n		Select components with <= n faces\n");
	fprintf(stderr, "    	-t n		Select the largest n components\n");
	exit(1);
}


int main(int argc, char *argv[])
{
	// Parse command line
	int morethan = 0;
	int lessthan = BIGNUM;
	int total = BIGNUM;
	int nprint = 20;
	bool conn_vert = false;
	bool splitcc = false;
	const char *infilename=NULL, *outfilename=NULL;

	int c;
	while ((c = getopt(argc, argv, "hvasm:l:t:")) != EOF) {
		switch (c) {
			case 'v': conn_vert = true; break;
			case 'a': nprint = BIGNUM; break;
			case 's': splitcc = true; break;
			case 'm': morethan = atoi(optarg); break;
			case 'l': lessthan = atoi(optarg); break;
			case 't': total = atoi(optarg); break;
			default: usage(argv[0]);
		}
	}
	if (argc - optind < 1)
		usage(argv[0]);
	infilename = argv[optind];
	if (argc - optind >= 2)
		outfilename = argv[optind+1];

	// Read input file
	TriMesh *in = TriMesh::read(infilename);
	if (!in) {
		fprintf(stderr, "Couldn't read input %s\n", infilename);
		exit(1);
	}
	bool had_tstrips = !in->tstrips.empty();
	in->need_faces();
	in->tstrips.clear();

	// Find connected components
	vector<int> comps;
	vector<int> compsizes;
	find_comps(in, comps, compsizes, conn_vert);

	// Print out the top few components
	print_comps(compsizes, morethan, lessthan, total, nprint);

	// Exit here if just printing things out, and not saving anything
	if (!outfilename)
		exit(0);

	// Get rid of the junk we don't want...
	if (lessthan != BIGNUM) {
		select_small_comps(in, comps, compsizes, lessthan,
			morethan != 0 ? BIGNUM : total);
		find_comps(in, comps, compsizes, conn_vert);
	}
	if (morethan != 0 || (lessthan == BIGNUM && total != BIGNUM)) {
		select_big_comps(in, comps, compsizes, morethan, total);
		find_comps(in, comps, compsizes, conn_vert);
	}

	if (splitcc) {
		// Split into separate files, if requested
		int ncomps = compsizes.size();
		for (int i = 0; i < ncomps; i++) {
			TriMesh *tmp = new TriMesh(*in);
			select_comp(tmp, comps, i);
			if (had_tstrips)
				tmp->need_tstrips();
			char filename[1024];
			sprintf(filename, "cc%d-%s", i+1, outfilename);
			tmp->write(filename);
			delete tmp;
		}
	} else {
		// Write the requested components to a file
		if (had_tstrips)
			in->need_tstrips();
		in->write(outfilename);
	}
}
