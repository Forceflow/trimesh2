/*
Szymon Rusinkiewicz
Princeton University

TriMesh_stats.cc
Computation of various statistics on the mesh.
*/

#include "TriMesh.h"
#include "KDtree.h"
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
				vals.push_back((float) neighbors[i].size());
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
						vertices[faces[i][NEXT_MOD3(j)]]));
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

	// Take absolute value or square
	switch (op) {
		case STAT_MINABS:
		case STAT_MAXABS:
		case STAT_SUMABS:
		case STAT_MEANABS:
			for (int i = 0; i < n; i++) {
				if (vals[i] < 0.0f)
					vals[i] = -vals[i];
			}
			break;

		case STAT_SUMSQR:
		case STAT_RMS:
			for (int i = 0; i < n; i++)
				vals[i] *= vals[i];
			break;

		default:
			break;
	}

	// Now do the computation
	switch (op) {
		case STAT_MIN:
		case STAT_MINABS:
			return *min_element(vals.begin(), vals.end());

		case STAT_MAX:
		case STAT_MAXABS:
			return *max_element(vals.begin(), vals.end());

		case STAT_SUM:
		case STAT_SUMABS:
		case STAT_SUMSQR:
			return accumulate(vals.begin(), vals.end(), 0.0f);

		case STAT_MEAN:
		case STAT_MEANABS:
			return accumulate(vals.begin(), vals.end(), 0.0f) / n;

		case STAT_RMS:
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

		default:
			return 0.0f; // Can't happen, I hope.
	}
}


// Fast computation of a characteristic "feature size" for the mesh.
// Computed as an approximation to the median edge length, or,
// in a point cloud, the median distance from a point to its nearest neighbor.
float TriMesh::feature_size()
{
	const int nsamples = 999;
	const float approx_eps = 0.05f;
	int nv = vertices.size();
	need_faces();
	int nf = faces.size();

	vector<float> samples;
	samples.reserve(nsamples);

	// We want to return consistent results, even though we're sampling,
	// so reset the RNG
	xorshift_rnd(0);

	// Accumulate samples
	if (nf > nsamples / 3) {
		// Big mesh - do sampling
		while (int(samples.size()) < nsamples) {
			int ind = uniform_rnd(nf);
			const point &p0 = vertices[faces[ind][0]];
			const point &p1 = vertices[faces[ind][1]];
			const point &p2 = vertices[faces[ind][2]];
			samples.push_back(dist2(p0,p1));
			samples.push_back(dist2(p1,p2));
			samples.push_back(dist2(p2,p0));
		}
	} else if (nf > 0) {
		// Small mesh - just loop over all faces
		for (int ind = 0; ind < nf; ind++) {
			const point &p0 = vertices[faces[ind][0]];
			const point &p1 = vertices[faces[ind][1]];
			const point &p2 = vertices[faces[ind][2]];
			samples.push_back(dist2(p0,p1));
			samples.push_back(dist2(p1,p2));
			samples.push_back(dist2(p2,p0));
		}
	} else if (nv > nsamples) {
		// Big point cloud - do sampling
		KDtree kd(vertices);
		while (int(samples.size()) < nsamples) {
			int ind = uniform_rnd(nv);
			const point &p = vertices[ind];
			const float *q = kd.closest_to_pt(p, 0.0f, approx_eps);
			samples.push_back(dist2(p, point(q)));
		}
	} else {
		// Small point cloud - just loop over all vertices
		KDtree kd(vertices);
		for (int ind = 0; ind < nv; ind++) {
			const point &p = vertices[ind];
			const float *q = kd.closest_to_pt(p, 0.0f, approx_eps);
			samples.push_back(dist2(p, point(q)));
		}
	}

	// Find median
	nth_element(samples.begin(),
	            samples.begin() + samples.size()/2,
	            samples.end());
	return sqrt(samples[samples.size()/2]);
}

} // namespace trimesh
