/*
Szymon Rusinkiewicz
Princeton University

ICP.cc
Routines for doing ICP.
*/

#include <cstring>
#include <algorithm>
#include "ICP.h"
#include "KDtree.h"
#include "timestamp.h"
#include "lineqn.h"
using namespace std;


#define MAX_ITERS 100
#define MIN_PAIRS 25
#define DESIRED_PAIRS 500
#define DESIRED_PAIRS_EARLY 50
#define DESIRED_PAIRS_FINAL 2000
#define COMPAT_THRESH 0.7f
#define TERM_THRESH 5
#define TERM_HIST 7
#define EIG_THRESH 0.01f
#define dprintf TriMesh::dprintf


// One or both of the following can be #defined
#define USE_GRID_FOR_OVERLAPS
#undef USE_KD_FOR_OVERLAPS


namespace trimesh {

// Quick 'n dirty portable random number generator 
static inline float tinyrnd()
{
	static unsigned trand = 0;
	trand = 1664525u * trand + 1013904223u;
	return (float) trand / 4294967296.0f;
}


// A pair of points, with an associated normal
struct PtPair {
	point p1, p2;
	vec norm;
	PtPair(const point &p1_, const point &p2_, const vec &norm_) :
			p1(p1_), p2(p2_), norm(norm_)
		{}
};


// A class for evaluating compatibility of normals during KDtree searches
class NormCompat : public KDtree::CompatFunc {
private:
	const vec n;
	TriMesh *m;
	bool pointcloud;

public:
	NormCompat(const vec &n_, TriMesh *m_, bool &p_):
		n(n_), m(m_), pointcloud(p_)
		{}
	virtual bool operator () (const float *p) const
	{
		int idx = (const point *)p - (const point *)&(m->vertices[0]);
		if (pointcloud)
			return fabs(n DOT m->normals[idx]) > COMPAT_THRESH;
		else if (m->is_bdy(idx))
			return true;
		else
			return (n DOT m->normals[idx]) > COMPAT_THRESH;
	}
};


// Find the median squared distance between points
static float median_dist2(const vector<PtPair> &pairs)
{
	size_t n = pairs.size();
	if (!n)
		return 0.0f;

	vector<float> distances2;
	distances2.reserve(n);
	for (size_t i = 0; i < n; i++)
		distances2.push_back(dist2(pairs[i].p1, pairs[i].p2));

	size_t pos = n / 2;
	nth_element(distances2.begin(),
		    distances2.begin() + pos,
		    distances2.end());
	return distances2[pos];
}


// A spatial grid datastructure for fast overlap computation
class Grid {
public:
	enum { GRID_SHIFT = 4, GRID_MAX = (1 << GRID_SHIFT) - 1 };
	float xmin, xmax, ymin, ymax, zmin, zmax, scale;
	vector<char> g;
	bool valid(const point &p)
	{
		return p[0] >= xmin && p[1] >= ymin && p[2] >= zmin &&
		       p[0] <= xmax && p[1] <= ymax && p[2] <= zmax;
	}
	int ind(const point &p)
	{
		int x = clamp(int(scale * (p[0] - xmin)), 0, int(GRID_MAX));
		int y = clamp(int(scale * (p[1] - ymin)), 0, int(GRID_MAX));
		int z = clamp(int(scale * (p[2] - zmin)), 0, int(GRID_MAX));
		return (x << (2*GRID_SHIFT)) + (y << GRID_SHIFT) + z;
	}
	bool overlaps(const point &p) { return valid(p) && g[ind(p)]; }
	Grid(const vector<point> &pts);
};


// Compute a Grid from a list of points
Grid::Grid(const vector<point> &pts)
{
	g.resize(1 << 3*GRID_SHIFT);
	if (pts.empty()) {
		xmin = xmax = ymin = ymax = zmin = zmax = scale = 0.0f;
		return;
	}
	xmin = xmax = pts[0][0];
	ymin = ymax = pts[0][1];
	zmin = zmax = pts[0][2];
	for (size_t i = 1; i < pts.size(); i++) {
		if (pts[i][0] < xmin)  xmin = pts[i][0];
		if (pts[i][0] > xmax)  xmax = pts[i][0];
		if (pts[i][1] < ymin)  ymin = pts[i][1];
		if (pts[i][1] > ymax)  ymax = pts[i][1];
		if (pts[i][2] < zmin)  zmin = pts[i][2];
		if (pts[i][2] > zmax)  zmax = pts[i][2];
	}
	scale = 1.0f / max(max(xmax-xmin, ymax-ymin), zmax-zmin);
	scale *= float(1 << GRID_SHIFT);
	for (size_t i = 0; i < pts.size(); i++)
		g[ind(pts[i])] = 1;
}


// Determine which points on s1 and s2 overlap the other, filling in o1 and o2
// Also fills in maxdist, if it is <= 0 on input
void compute_overlaps(TriMesh *s1, TriMesh *s2,
		      const xform &xf1, const xform &xf2,
		      const KDtree *, const KDtree *,
		      vector<float> &o1, vector<float> &o2,
		      float &maxdist, int verbose)
{
	size_t nv1 = s1->vertices.size(), nv2 = s2->vertices.size();

	timestamp t = now();
	Grid g1(s1->vertices);
	Grid g2(s2->vertices);
	xform xf12 = inv(xf2) * xf1;
	xform xf21 = inv(xf1) * xf2;
	if (maxdist <= 0.0f)
		maxdist = min(1.0f / g1.scale, 1.0f / g2.scale);

#ifdef USE_KD_FOR_OVERLAPS
	float maxdist2 = sqr(maxdist);
	bool pointcloud1 = (s1->faces.empty() && s1->grid.empty() && s1->tstrips.empty());
	bool pointcloud2 = (s2->faces.empty() && s2->grid.empty() && s2->tstrips.empty());
#endif

	o1.resize(nv1);
	for (size_t i = 0; i < nv1; i++) {
		o1[i] = 0;
		point p = xf12 * s1->vertices[i];
#ifdef USE_GRID_FOR_OVERLAPS
		if (!g2.overlaps(p))
			continue;
#endif
#ifdef USE_KD_FOR_OVERLAPS
		const float *match = kd2->closest_to_pt(p, maxdist2);
		if (!match)
			continue;
		if (!pointcloud2 &&
		    s2->is_bdy((match - (const float *) &s2->vertices[0][0]) / 3))
			continue;
#endif
		o1[i] = 1;
	}

	o2.resize(nv2);
	for (size_t i = 0; i < nv2; i++) {
		o2[i] = 0;
		point p = xf21 * s2->vertices[i];
#ifdef USE_GRID_FOR_OVERLAPS
		if (!g1.overlaps(p))
			continue;
#endif
#ifdef USE_KD_FOR_OVERLAPS
		const float *match = kd1->closest_to_pt(p, maxdist2);
		if (!match)
			continue;
		if (!pointcloud1 &&
		    s1->is_bdy((match - (const float *) &s1->vertices[0][0]) / 3))
			continue;
#endif
		o2[i] = 1;
	}
	if (verbose > 1) {
		dprintf("Computed overlaps in %.2f msec.\n",
			(now() - t) * 1000.0);
	}
}


// Select a number of points and find correspondences 
static void select_and_match(TriMesh *s1, TriMesh *s2,
			     const xform &xf1, const xform &xf2,
			     const KDtree *kd2, const vector<float> &sampcdf1,
			     float incr, float maxdist, int /* verbose */,
			     vector<PtPair> &pairs, bool flip)
{
	xform xf1r = norm_xf(xf1);
	xform xf2r = norm_xf(xf2);
	xform xf12 = inv(xf2) * xf1;
	xform xf12r = norm_xf(xf12);
	float maxdist2 = sqr(maxdist);

	size_t i = 0;
	float cval = 0.0f;
	while (1) {
		cval += incr * tinyrnd();
		if (cval >= 1.0f)
			break;
		while (sampcdf1[i] <= cval)
			i++;
		cval = sampcdf1[i];

		point p = xf12 * s1->vertices[i];
		vec n = xf12r * s1->normals[i];

		// Do the matching
		bool pointcloud2 = (s2->faces.empty() && s2->tstrips.empty());
		NormCompat nc(n, s2, pointcloud2);

		const float *match = kd2->closest_to_pt(p, maxdist2, &nc);
		if (!match)
			continue;
		int imatch = (match - (const float *) &(s2->vertices[0][0])) / 3;
		if (!pointcloud2 && s2->is_bdy(imatch))
			continue;

		// Project both points into world coords and save 
		if (flip) {
			pairs.push_back(PtPair(xf2  * s2->vertices[imatch],
					       xf1  * s1->vertices[i],
					       xf2r * s2->normals[imatch]));
		} else {
			pairs.push_back(PtPair(xf1  * s1->vertices[i],
					       xf2  * s2->vertices[imatch],
					       xf1r * s1->normals[i]));
		}
	}
}


// Compute ICP alignment matrix, including eigenvector decomposition
static void compute_ICPmatrix(const vector<PtPair> &pairs,
			      float evec[6][6], float eval[6], float b[6],
			      point &centroid, float &scale, float &err)
{
	size_t n = pairs.size();

	centroid = point(0,0,0);
	for (size_t i = 0; i < n; i++)
		centroid += pairs[i].p2;
	centroid /= float(n);

	scale = 0.0f;
	for (size_t i = 0; i < n; i++)
		scale += dist2(pairs[i].p2, centroid);
	scale /= float(n);
	scale = 1.0f / sqrt(scale);

	memset(&evec[0][0], 0, 6*6*sizeof(float));
	memset(&b[0], 0, 6*sizeof(float));

	err = 0.0f;
	for (size_t i = 0; i < n; i++) {
		const point &p1 = pairs[i].p1;
		const point &p2 = pairs[i].p2;
		const vec &n = pairs[i].norm;

		float d = (p1 - p2) DOT n;
		d *= scale;
		vec p2c = p2 - centroid;
		p2c *= scale;
		vec c = p2c CROSS n;

		err += d * d;
		float x[6] = { c[0], c[1], c[2], n[0], n[1], n[2] };
		for (int j = 0; j < 6; j++) {
			b[j] += d * x[j];
			for (int k = 0; k < 6; k++)
				evec[j][k] += x[j] * x[k];
		}
	}

	err /= float(n);
	err = sqrt(err) / scale;
	eigdc<float,6>(evec, eval);
}


// Compute ICP alignment, given matrix computed by compute_ICPmatrix
static void compute_alignxf(float evec[6][6], float eval[6], float b[6],
			    point &centroid, float scale, xform &alignxf)
{
	float einv[6];
	for (int i = 0; i < 6; i++) {
		if (eval[i] < EIG_THRESH * eval[5])
			einv[i] = 0.0f;
		else
			einv[i] = 1.0f / eval[i];
	}
	float x[6];
	eigmult<float,6>(evec, einv, b, x);

	// Interpret results
	float sx = min(max(x[0], -1.0f), 1.0f);
	float sy = min(max(x[1], -1.0f), 1.0f);
	float sz = min(max(x[2], -1.0f), 1.0f);
	float cx = sqrt(1.0f - sx*sx);
	float cy = sqrt(1.0f - sy*sy);
	float cz = sqrt(1.0f - sz*sz);

	alignxf[0]  = cy*cz;
	alignxf[1]  = sx*sy*cz + cx*sz;
	alignxf[2]  = -cx*sy*cz + sx*sz;
	alignxf[3]  = 0;
	alignxf[4]  = -cy*sz;
	alignxf[5]  = -sx*sy*sz + cx*cz;
	alignxf[6]  = cx*sy*sz + sx*cz;
	alignxf[7]  = 0;
	alignxf[8]  = sy;
	alignxf[9]  = -sx*cy;
	alignxf[10] = cx*cy;
	alignxf[11] = 0;
	alignxf[12] = x[3] / scale + centroid[0] - alignxf[0]*centroid[0] -
		      alignxf[4]*centroid[1] - alignxf[8]*centroid[2];
	alignxf[13] = x[4] / scale + centroid[1] - alignxf[1]*centroid[0] -
		      alignxf[5]*centroid[1] - alignxf[9]*centroid[2];
	alignxf[14] = x[5] / scale + centroid[2] - alignxf[2]*centroid[0] -
		      alignxf[6]*centroid[1] - alignxf[10]*centroid[2];
	alignxf[15] = 1;
}


// Compute isotropic or anisotropic scale
void compute_scale(const vector<PtPair> &pairs, xform &alignxf,
		   int verbose, bool do_affine)
{
	int n = pairs.size();

	// Compute COM
	point centroid;
	for (int i = 0; i < n; i++)
		centroid += pairs[i].p1 + pairs[i].p2;
	centroid /= 2.0f * n;
	xform txf = xform::trans(centroid);

	// Compute covariance matrices
	double cov1[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
	double cov2[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
	for (int i = 0; i < n; i++) {
		vec p = pairs[i].p1 - centroid;
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				cov1[j][k] += p[j]*p[k];
		p = pairs[i].p2 - centroid;
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				cov2[j][k] += p[j]*p[k];
	}

	// Compute eigenstuff of cov
	double eval1[3], eval2[3];
	eigdc<double,3>(cov1, eval1);
	eigdc<double,3>(cov2, eval2);

	if (!do_affine) {
		// Just uniform scale
		float s1 = (float)sqrt(eval1[0] + eval1[1] + eval1[2]);
		float s2 = (float)sqrt(eval2[0] + eval2[1] + eval2[2]);
		alignxf = txf * xform::scale(s1, s1, s1) *
			  inv(xform::scale(s2, s2, s2)) * inv(txf);
		return;
	}

	// Compute sqrt of covariance
	double csqrt1[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
	double csqrt2[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
	for (int i = 0; i < 3; i++)
		eval1[i] = sqrt(eval1[i] / n);
	for (int i = 0; i < 3; i++)
		eigmult<double,3>(cov1, eval1, csqrt1[i], csqrt1[i]);
	for (int i = 0; i < 3; i++)
		eval2[i] = sqrt(eval2[i] / n);
	for (int i = 0; i < 3; i++)
		eigmult<double,3>(cov2, eval2, csqrt2[i], csqrt2[i]);

	if (verbose > 1) {
		dprintf("sqrt covariance 1 =");
		for (int j = 0; j < 3; j++) {
			dprintf("\n\t");
			for (int k = 0; k < 3; k++)
				dprintf("%10.3f ", csqrt1[j][k]);
		}
		dprintf("\nsqrt covariance 2 =");
		for (int j = 0; j < 3; j++) {
			dprintf("\n\t");
			for (int k = 0; k < 3; k++)
				dprintf("%10.3f ", csqrt2[j][k]);
		}
		dprintf("\n");
	}


	xform cxf1 = xform(csqrt1[0][0], csqrt1[1][0], csqrt1[2][0], 0,
			   csqrt1[0][1], csqrt1[1][1], csqrt1[2][1], 0,
			   csqrt1[0][2], csqrt1[1][2], csqrt1[2][2], 0,
			   0, 0, 0, 1);
	xform cxf2 = xform(csqrt2[0][0], csqrt2[1][0], csqrt2[2][0], 0,
			   csqrt2[0][1], csqrt2[1][1], csqrt2[2][1], 0,
			   csqrt2[0][2], csqrt2[1][2], csqrt2[2][2], 0,
			   0, 0, 0, 1);

	alignxf = txf * cxf1 * inv(cxf2) * inv(txf);
}


// Do one iteration of ICP
static float ICP_iter(TriMesh *s1, TriMesh *s2, const xform &xf1, xform &xf2,
		      const KDtree *kd1, const KDtree *kd2,
		      const vector<float> &weights1, const vector<float> &weights2,
		      float &maxdist, int verbose,
		      vector<float> &sampcdf1, vector<float> &sampcdf2,
		      float &incr, bool update_cdfs,
		      bool do_scale, bool do_affine)
{
	// Compute pairs
	timestamp t1 = now();
	if (verbose > 1)
		dprintf("maxdist = %f\n", maxdist);
	vector<PtPair> pairs;
	select_and_match(s1, s2, xf1, xf2, kd2, sampcdf1, incr,
			 maxdist, verbose, pairs, false);
	select_and_match(s2, s1, xf2, xf1, kd1, sampcdf2, incr,
			 maxdist, verbose, pairs, true);

	timestamp t2 = now();
	size_t np = pairs.size();
	if (verbose > 1) {
		dprintf("Generated %lu pairs in %.2f msec.\n",
			(unsigned long) np, (t2-t1) * 1000.0);
	}

	// Reject pairs with distance > 2.5 sigma
	float thresh = 13.73818f * median_dist2(pairs);
	if (verbose > 1)
		dprintf("Rejecting pairs > %f\n", sqrt(thresh));
	size_t next = 0;
	for (size_t i = 0; i < np; i++) {
		if (dist2(pairs[i].p1, pairs[i].p2) <= thresh)
			pairs[next++] = pairs[i];
	}
	pairs.erase(pairs.begin() + next, pairs.end());

	timestamp t3 = now();
	if (verbose > 1) {
		dprintf("Rejected %lu pairs in %.2f msec.\n",
			(unsigned long) (np - pairs.size()), (t3-t2) * 1000.0);
	}
	if (pairs.size() < MIN_PAIRS) {
		if (verbose)
			dprintf("Too few point pairs.\n");
		return -1.0f;
	}

	// Update incr and maxdist based on what happened here
	incr *= (float) pairs.size() / DESIRED_PAIRS;
	maxdist = max(2.0f * sqrt(thresh), 0.7f * maxdist);

	// Do the minimization
	float evec[6][6], eval[6], b[6], scale, err;
	point centroid;
	xform alignxf;
	compute_ICPmatrix(pairs, evec, eval, b, centroid, scale, err);
	if (verbose > 1) {
		dprintf("RMS point-to-plane error = %f\n", err);
		for (int i = 0; i < 5; i++)
			if (eval[i] < EIG_THRESH * eval[5])
				dprintf("Small eigenvalue %f (largest is %f)\n", eval[i], eval[5]);
	}
	compute_alignxf(evec, eval, b, centroid, scale, alignxf);
	xf2 = alignxf * xf2;

	if (do_scale || do_affine) {
		for (size_t i = 0; i < pairs.size(); i++)
			pairs[i].p2 = alignxf * pairs[i].p2;
		compute_scale(pairs, alignxf, verbose, do_affine);
		xf2 = alignxf * xf2;
	}

	timestamp t4 = now();
	if (verbose > 1) {
		dprintf("Computed xform in %.2f msec.\n",
			(t4-t3) * 1000.0);
	}

	// Update CDFs, if necessary
	if (!update_cdfs)
		return err;

	float einv[6];
	for (int i = 0; i < 6; i++)
		einv[i] = 1.0f / max(eval[i], EIG_THRESH * eval[5]);
	float Cinv[6][6];
	for (int i = 0; i < 6; i++) {
		float x[6];
		for (int j = 0; j < 6; j++)
			x[j] = (j == i) ? 1.0f : 0.0f;
		eigmult<float,6>(evec, einv, x, x);
		for (int j = 0; j < 6; j++)
			Cinv[i][j] = x[j];
	}

	xform xf1r = norm_xf(xf1);
	size_t n1 = s1->vertices.size();
	for (size_t i = 0; i < n1; i++) {
		sampcdf1[i] = 0.0;
		if (!weights1[i])
			continue;
		point p = xf1 * s1->vertices[i];
		p -= centroid;
		p *= scale;
		vec n = xf1r * s1->normals[i];
		vec c = p CROSS n;
		for (int j = 0; j < 6; j++) {
			float tmp = Cinv[j][0] * c[0] + Cinv[j][1] * c[1] +
				    Cinv[j][2] * c[2] + Cinv[j][3] * n[0] +
				    Cinv[j][4] * n[1] + Cinv[j][5] * n[2];
			if (j < 3)
				sampcdf1[i] += tmp * c[j];
			else
				sampcdf1[i] += tmp * n[j-3];
		}
		sampcdf1[i] *= weights1[i];
	}
	for (size_t i = 1; i < n1; i++)
		sampcdf1[i] += sampcdf1[i-1];
	if (!sampcdf1[n1-1]) {
		if (verbose)
			dprintf("No overlap.\n");
		return -1.0f;
	}
	float cscale = 1.0f / sampcdf1[n1-1];
	for (size_t i = 0; i < n1-1; i++)
		sampcdf1[i] *= cscale;
	sampcdf1[n1-1] = 1.0f;

	xform xf2r = norm_xf(xf2);
	size_t n2 = s2->vertices.size();
	for (size_t i = 0; i < n2; i++) {
		sampcdf2[i] = 0.0;
		if (!weights2[i])
			continue;
		point p = xf2 * s2->vertices[i];
		p -= centroid;
		p *= scale;
		vec n = xf2r * s2->normals[i];
		vec c = p CROSS n;
		for (int j = 0; j < 6; j++) {
			float tmp = Cinv[j][0] * c[0] + Cinv[j][1] * c[1] +
				    Cinv[j][2] * c[2] + Cinv[j][3] * n[0] +
				    Cinv[j][4] * n[1] + Cinv[j][5] * n[2];
			if (j < 3)
				sampcdf2[i] += tmp * c[j];
			else
				sampcdf2[i] += tmp * n[j-3];
		}
		sampcdf2[i] *= weights2[i];
	}
	for (size_t i = 1; i < n2; i++)
		sampcdf2[i] += sampcdf2[i-1];
	cscale = 1.0f / sampcdf2[n2-1];
	if (!sampcdf2[n2-1]) {
		if (verbose)
			dprintf("No overlap.\n");
		return -1.0f;
	}
	for (size_t i = 0; i < n2-1; i++)
		sampcdf2[i] *= cscale;
	sampcdf2[n2-1] = 1.0f;

	timestamp t5 = now();
	if (verbose > 1) {
		dprintf("Updated CDFs in %.2f msec.\n",
			(t5-t4) * 1000.0);
	}

	return err;
}


// Do one iteration of point-to-point ICP (this is done in the early stages
// to assure stability)
static float ICP_p2pt(TriMesh *s1, TriMesh *s2, const xform &xf1, xform &xf2,
		      const KDtree *kd1, const KDtree *kd2,
		      float &maxdist, int verbose,
		      vector<float> &sampcdf1, vector<float> &sampcdf2,
		      float &incr, bool trans_only)
{
	// Compute pairs
	timestamp t1 = now();
	if (verbose > 1)
		dprintf("maxdist = %f\n", maxdist);
	vector<PtPair> pairs;
	select_and_match(s1, s2, xf1, xf2, kd2, sampcdf1, incr,
			 maxdist, verbose, pairs, false);
	select_and_match(s2, s1, xf2, xf1, kd1, sampcdf2, incr,
			 maxdist, verbose, pairs, true);

	timestamp t2 = now();
	size_t np = pairs.size();
	if (verbose > 1) {
		dprintf("Generated %lu pairs in %.2f msec.\n",
			(unsigned long) np, (t2-t1) * 1000.0);
	}

	// Reject pairs with distance > 3 sigma
	float thresh = 19.782984f * median_dist2(pairs);
	if (verbose > 1)
		dprintf("Rejecting pairs > %f\n", sqrt(thresh));
	size_t next = 0;
	for (size_t i = 0; i < np; i++) {
		if (dist2(pairs[i].p1, pairs[i].p2) <= thresh)
			pairs[next++] = pairs[i];
	}
	pairs.erase(pairs.begin() + next, pairs.end());

	timestamp t3 = now();
	if (verbose > 1) {
		dprintf("Rejected %lu pairs in %.2f msec.\n",
			(unsigned long) (np - pairs.size()), (t3-t2) * 1000.0);
	}
	if ((int)pairs.size() < (trans_only ? 1 : MIN_PAIRS)) {
		if (verbose)
			dprintf("Too few point pairs.\n");
		return -1.0f;
	}

	// Update incr and maxdist based on what happened here
	incr *= (float) pairs.size() / DESIRED_PAIRS_EARLY;
	maxdist = max(1.5f * sqrt(thresh), 0.7f * maxdist);

	// Do the minimization
	point centroid1, centroid2;
	for (size_t i = 0; i < pairs.size(); i++) {
		centroid1 += pairs[i].p1;
		centroid2 += pairs[i].p2;
	}
	centroid1 /= (float) pairs.size();
	centroid2 /= (float) pairs.size();

	xform alignxf = xform::trans(centroid1 - centroid2);

	double A[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
	double B[3] = {0,0,0};
	double sum = 0;
	for (size_t i = 0; i < pairs.size(); i++) {
		vec p12 = pairs[i].p1 - pairs[i].p2;
		vec p2c = pairs[i].p2 - centroid2;
		vec c = p2c CROSS p12;
		sum += len2(p12);
		B[0] += c[0]; B[1] += c[1]; B[2] += c[2];
		A[0][0] += sqr(p2c[1]) + sqr(p2c[2]);
		A[0][1] -= p2c[0] * p2c[1];
		A[0][2] -= p2c[0] * p2c[2];
		A[1][1] += sqr(p2c[0]) + sqr(p2c[2]);
		A[1][2] -= p2c[1] * p2c[2];
		A[2][2] += sqr(p2c[0]) + sqr(p2c[1]);
	}
	float err = (float)sqrt(sum / pairs.size());
	if (verbose > 1)
		dprintf("RMS point-to-point error = %f\n", err);

	if (!trans_only) {
		double diag[3];
		ldltdc<double,3>(A, diag);
		ldltsl<double,3>(A, diag, B, B);
		alignxf = xform::trans(centroid1) *
			  xform::rot(B[0], 1, 0, 0) *
			  xform::rot(B[1], 0, 1, 0) *
			  xform::rot(B[2], 0, 0, 1) *
			  xform::trans(-centroid2);
	}

	xf2 = alignxf * xf2;

	timestamp t4 = now();
	if (verbose > 1) {
		dprintf("Computed xform in %.2f msec.\n",
			(t4-t3) * 1000.0);
	}
	return err;
}


// Do ICP.  Aligns mesh s2 to s1, updating xf2 with the new transform.
// Returns alignment error, or -1 on failure
float ICP(TriMesh *s1, TriMesh *s2, const xform &xf1, xform &xf2,
	  const KDtree *kd1, const KDtree *kd2,
	  vector<float> &weights1, vector<float> &weights2,
	  float maxdist /* = 0.0f */, int verbose /* = 0 */,
	  bool do_scale /* = false */, bool do_affine /* = false */)
{
	// Make sure we have everything precomputed
	s1->need_normals();  s2->need_normals();
	if (!s1->faces.empty() || !s1->tstrips.empty()) {
		s1->need_neighbors();
		s1->need_adjacentfaces();
	}
	if (!s2->faces.empty() || !s2->tstrips.empty()) {
		s2->need_neighbors();
		s2->need_adjacentfaces();
	}
	size_t nv1 = s1->vertices.size(), nv2 = s2->vertices.size();

	timestamp t = now();

	if (maxdist <= 0.0f) {
		s1->need_bbox();
		s2->need_bbox();
		maxdist = 0.5f * min(len(s1->bbox.size()), len(s2->bbox.size()));
	}
	// Compute initial CDFs
	vector<float> sampcdf1(nv1), sampcdf2(nv2);
	for (size_t i = 0; i < nv1-1; i++)
		sampcdf1[i] = (float) (i+1) / nv1;
	sampcdf1[nv1-1] = 1.0f;
	for (size_t i = 0; i < nv2-1; i++)
		sampcdf2[i] = (float) (i+1) / nv2;
	sampcdf2[nv2-1] = 1.0f;

	// Do a few p2pt iterations
	float incr = 4.0f / DESIRED_PAIRS_EARLY;
	for (int i = 0; i < 2; i++) {
		if (ICP_p2pt(s1, s2, xf1, xf2, kd1, kd2, maxdist, verbose,
			     sampcdf1, sampcdf2, incr, true) < 0.0f)
			return -1.0f;
	}
	for (int i = 0; i < 5; i++) {
		if (ICP_p2pt(s1, s2, xf1, xf2, kd1, kd2, maxdist, verbose,
			     sampcdf1, sampcdf2, incr, false) < 0.0f)
			return -1.0f;
	}

	// Do a point-to-plane iteration and update CDFs
	if (weights1.size() != nv1 || weights2.size() != nv2)
		compute_overlaps(s1, s2, xf1, xf2, kd1, kd2,
				 weights1, weights2, maxdist, verbose);
	float err = ICP_iter(s1, s2, xf1, xf2, kd1, kd2, weights1, weights2,
			     maxdist, verbose, sampcdf1, sampcdf2,
			     incr, true, false, false);
	if (verbose > 1) {
		timestamp tnow = now();
		dprintf("Time for initial iterations: %.2f msec.\n\n",
		       (tnow-t) * 1000.0);
		t = tnow;
	}
	if (err < 0.0f)
		return err;

	bool rigid_only = true;
	int iters = 0;
	vector<int> err_delta_history(TERM_HIST);
	do {
		float lasterr = err;
		if (verbose > 1)
			dprintf("Using incr = %f\n", incr);
		bool recompute = (iters % 10 == 9);
		if (recompute)
			compute_overlaps(s1, s2, xf1, xf2, kd1, kd2,
					 weights1, weights2, maxdist, verbose);
		err = ICP_iter(s1, s2, xf1, xf2, kd1, kd2, weights1, weights2,
			       maxdist, verbose, sampcdf1, sampcdf2, incr,
			       recompute, do_scale && !rigid_only,
			       do_affine && !rigid_only);
		if (verbose > 1) {
			timestamp tnow = now();
			dprintf("Time for this iteration: %.2f msec.\n\n",
			       (tnow-t) * 1000.0);
			t = tnow;
		}
		if (err < 0.0f)
			return err;

		// Check whether the error's been going up or down lately.
		// Specifically, we break out if error has gone up in
		// TERM_THRESH out of the last TERM_HIST iterations.
		for (int i = 0; i < TERM_HIST - 1; i++)
			err_delta_history[i] = err_delta_history[i+1];
		err_delta_history[TERM_HIST - 1] = (err >= lasterr);
		int nincreases = 0;
		for (int i = 0; i < TERM_HIST; i++)
			nincreases += err_delta_history[i];
		if (nincreases >= TERM_THRESH) {
			if (!rigid_only || (!do_scale && !do_affine))
				break;
			err_delta_history.clear();
			err_delta_history.resize(TERM_HIST);
			rigid_only = false;
		}
	} while (++iters < MAX_ITERS);

	if (verbose > 1)
		dprintf("Did %d iterations\n\n", iters);

	// One final iteration at a higher sampling rate...
	if (verbose > 1)
		dprintf("Last iteration...\n");
	incr *= (float) DESIRED_PAIRS / DESIRED_PAIRS_FINAL;
	if (verbose > 1)
		dprintf("Using incr = %f\n", incr);
	err = ICP_iter(s1, s2, xf1, xf2, kd1, kd2, weights1, weights2,
		       maxdist, verbose, sampcdf1, sampcdf2, incr,
		       false, do_scale, do_affine);
	if (verbose > 1) {
		timestamp tnow = now();
		dprintf("Time for this iteration: %.2f msec.\n\n",
		       (tnow-t) * 1000.0);
		t = tnow;
	}
	return err;
}


// Easier-to-use interface to ICP
float ICP(TriMesh *s1, TriMesh *s2, const xform &xf1, xform &xf2,
	  int verbose /* = 0 */,
	  bool do_scale /* = false */, bool do_affine /* = false */)
{
	KDtree *kd1 = new KDtree(s1->vertices);
	KDtree *kd2 = new KDtree(s2->vertices);
	vector<float> weights1, weights2;
	float icperr = ICP(s1, s2, xf1, xf2, kd1, kd2,
			   weights1, weights2, 0.0f, verbose,
			   do_scale, do_affine);
	delete kd2;
	delete kd1;
	return icperr;
}

}; // namespace trimesh
