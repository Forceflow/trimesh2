/*
Szymon Rusinkiewicz
Princeton University

mesh_crunch.cc
Replacement for plycrunch - decimate a mesh using the Rossignac-Borrel
method of vertex collapse.
*/

#include "TriMesh.h"
#include <cstdio>
#include <cstdlib>
using namespace std;
using namespace trimesh;


// A hash table for voxels -
//  stores an index (into the new vertex table) for each voxel.
class VHashTable {
public:
	typedef Vec<3,int> voxel_coord_t;
	typedef size_t key_t;
	typedef size_t data_t;
	static const size_t MAGIC1 = 100000007;
	static const size_t MAGIC2 = 161803409;
	static const size_t MAGIC3 = 423606823;
	static const size_t SIZE_FUDGE = 2;
	static const data_t NO_DATA = 0xffffffffu;

private:
	float scale;
	vector<voxel_coord_t> voxel_coords;
	vector<data_t> data;

public:
	VHashTable(size_t maxpoints, float voxelsize_) :
		scale(1.0f / voxelsize_)
	{
		size_t n = SIZE_FUDGE * maxpoints;
		voxel_coords.resize(n);
		// Work around a bizarre compiler bug
		data_t tmp = NO_DATA;
		data.resize(n, tmp);
	}
	data_t &operator[] (const point &p)
	{
		voxel_coord_t c(int(floor(p[0] * scale)),
				int(floor(p[1] * scale)),
				int(floor(p[2] * scale)));
		key_t key = MAGIC1 * c[0] + MAGIC2 * c[1] + MAGIC3 * c[2];
		key %= data.size();

		// Open hashing
		while (1) {
			if (data[key] == NO_DATA) {
				voxel_coords[key] = c;
				break;
			} else if (voxel_coords[key] == c) {
				break;
			}
			key++;
			if (key == data.size())
				key = 0;
		}

		return data[key];
	}
};


// Perform the vertex collapse
void crunch(TriMesh *in, TriMesh *out, float voxelsize)
{
	size_t nv = in->vertices.size();
	// in->flags stores which output vertex we mapped to
	in->flags.clear();
	in->flags.resize(nv);
	out->vertices.reserve(nv);
	// out->flags stores number of vertices that mapped here
	out->flags.reserve(nv);

	// Go through input vertices.  If we haven't encountered the
	// voxel yet, create a new vertex in the output.  Else, merge with
	// other vertices in the same voxel.
	VHashTable hash(nv, voxelsize);
	for (size_t i = 0; i < nv; i++) {
		size_t &ind = hash[in->vertices[i]];
		if (ind >= nv) {
			// Failed lookup - add a new output vertex
			ind = out->vertices.size();
			out->vertices.push_back(in->vertices[i]);
			out->flags.push_back(1);
		} else {
			out->vertices[ind] += in->vertices[i];
			out->flags[ind]++;
		}
		in->flags[i] = ind;
	}

	// Divide by number of points per voxel, to get final vertex coords
	size_t onv = out->vertices.size();
	for (size_t i = 0; i < onv; i++)
		out->vertices[i] /= float(out->flags[i]);

	// Create new faces: only the non-degenerate ones
	in->need_faces();
	size_t nf = in->faces.size();
	out->faces.reserve(min(nf, onv*3));
	for (size_t i = 0; i < nf; i++) {
		size_t ind0 = in->flags[in->faces[i][0]];
		size_t ind1 = in->flags[in->faces[i][1]];
		size_t ind2 = in->flags[in->faces[i][2]];
		if (ind0 == ind1 || ind0 == ind2 || ind1 == ind2)
			continue;
		out->faces.push_back(TriMesh::Face(ind0,ind1,ind2));
	}
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s in.ply out.ply\n", myname);
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc < 3)
		usage(argv[0]);
	TriMesh *in = TriMesh::read(argv[1]);
	if (!in)
		usage(argv[0]);

	bool had_tstrips = !in->tstrips.empty();

	TriMesh *out = new TriMesh;
	float voxelsize = 2.0f * in->feature_size();
	crunch(in, out, voxelsize);

	if (had_tstrips)
		out->need_tstrips();
	out->write(argv[2]);
}
