/*
Szymon Rusinkiewicz
Princeton University

TriMesh_stats.cc
Computation of various statistics on the mesh.
*/

#include "TriMesh.h"
#include <algorithm>
#include <numeric>
using namespace std;


namespace trimesh {

// Compute a variety of statistics.  Takes a type of statistic to compute,
// and what to do with it.
float TriMesh::stat(StatOp op, StatVal val)
{
	vector<float> vals;

	switch (val) {
		case STAT_VALENCE: {
			need_neighbors();
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back(neighbors[i].size());
			break;
		}
		case STAT_FACEAREA: {
			need_faces();
			int nf = faces.size();
			for (int i = 0; i < nf; i++)
				vals.push_back(len(trinorm(i)));
			break;
		}
		case STAT_ANGLE: {
			need_faces();
			int nf = faces.size();
			for (int i = 0; i < nf; i++)
				for (int j = 0; j < 3; j++)
					vals.push_back(cornerangle(i, j));
			break;
		}
		case STAT_DIHEDRAL: {
			need_across_edge();
			int nf = faces.size();
			for (int i = 0; i < nf; i++)
				for (int j = 0; j < 3; j++) {
					if (across_edge[i][j] < 0)
						continue;
					vals.push_back(dihedral(i, j));
				}
			break;
		}
		case STAT_EDGELEN: {
			need_faces();
			int nf = faces.size();
			for (int i = 0; i < nf; i++)
				for (int j = 0; j < 3; j++)
					vals.push_back(dist(
						vertices[faces[i][j]],
						vertices[faces[i][(j+1)%3]]));
			break;
		}
		case STAT_X: {
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back(vertices[i][0]);
			break;
		}
		case STAT_Y: {
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back(vertices[i][1]);
			break;
		}
		case STAT_Z: {
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back(vertices[i][2]);
			break;
		}
		default:
			return 0.0f;
	}

	int n = vals.size();
	if (!n)
		return 0.0f;

	switch (op) {
		case STAT_MIN:
			return *min_element(vals.begin(), vals.end());

		case STAT_MAX:
			return *max_element(vals.begin(), vals.end());

		case STAT_MEANABS:
			for (int i = 0; i < n; i++)
				if (vals[i] < 0.0f)
					vals[i] = -vals[i];
			// Fall through
		case STAT_MEAN:
			return accumulate(vals.begin(), vals.end(), 0.0f) / n;

		case STAT_RMS:
			for (int i = 0; i < n; i++)
				vals[i] *= vals[i];
			return sqrt(accumulate(vals.begin(), vals.end(), 0.0f)
				/ n);

		case STAT_MEDIAN:
			if (n & 1) {
				nth_element(vals.begin(),
					    vals.begin() + n/2,
					    vals.end());
				return vals[n/2];
			} else {
				nth_element(vals.begin(),
					    vals.begin() + n/2 - 1,
					    vals.end());
				float tmp = vals[n/2 - 1];
				nth_element(vals.begin(),
					    vals.begin() + n/2,
					    vals.end());
				return 0.5f * (tmp + vals[n/2]);
			}

		case STAT_STDEV: {
			float mean = accumulate(vals.begin(), vals.end(), 0.0f)
				/ n;
			for (int i = 0; i < n; i++)
				vals[i] = sqr(vals[i] - mean);
			return sqrt(accumulate(vals.begin(), vals.end(), 0.0f)
				/ n);
		}

		case STAT_TOTAL:
			return accumulate(vals.begin(), vals.end(), 0.0f);
	}

	return 0.0f; // Can't happen, I hope.
}


// A characteristic "feature size" for the mesh.  Computed as an approximation
// to the median edge length
float TriMesh::feature_size()
{
	need_faces();
	if (faces.empty())
		return 0.0f;

	int nf = faces.size();
	int nsamp = min(nf / 2, 333);

	vector<float> samples;
	samples.reserve(nsamp * 3);

	for (int i = 0; i < nsamp; i++) {
		// Quick 'n dirty portable random number generator
		static unsigned randq = 0;
		randq = unsigned(1664525) * randq + unsigned(1013904223);

		int ind = randq % nf;
		const point &p0 = vertices[faces[ind][0]];
		const point &p1 = vertices[faces[ind][1]];
		const point &p2 = vertices[faces[ind][2]];
		samples.push_back(dist2(p0,p1));
		samples.push_back(dist2(p1,p2));
		samples.push_back(dist2(p2,p0));
	}
	nth_element(samples.begin(),
		    samples.begin() + samples.size()/2,
		    samples.end());
	return sqrt(samples[samples.size()/2]);
}

}; // namespace trimesh
