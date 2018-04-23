#ifndef ICP_H
#define ICP_H
/*
Szymon Rusinkiewicz
Princeton University

ICP.h
Iterative Closest Point alignment using symmetric point-to-plane minimization
and adaptive outlier rejection.
*/

#include "TriMesh.h"
#include "XForm.h"
#include "KDtree.h"


namespace trimesh {

class Grid;

// Make a grid for accelerating overlap computation
extern Grid *make_grid(TriMesh *mesh);

// Determine which points on s1 and s2 overlap the other, filling in o1 and o2
// Also fills in maxdist, if it is <= 0 on input
extern void compute_overlaps(TriMesh *mesh1, TriMesh *mesh2,
                             const xform &xf1, const xform &xf2,
                             const KDtree *kd1, const KDtree *kd2,
                             const Grid *g1, const Grid *g2,
                             ::std::vector<float> &o1, ::std::vector<float> &o2,
                             float &maxdist, int verbose);

// Transformation for which ICP is solving
enum ICP_xform_type {
	ICP_TRANSLATION, ICP_RIGID, ICP_SIMILARITY, ICP_AFFINE
};


// Do ICP.  Aligns mesh2 to mesh1, updating xf2 with the new transform.
// Returns alignment error, or -1 on failure.
// Pass in 0 for maxdist to figure it out...
// Pass in empty vector for weights to figure it out...
extern float ICP(TriMesh *mesh1, TriMesh *mesh2,
                 const xform &xf1, xform &xf2,
                 const KDtree *kd1, const KDtree *kd2,
                 ::std::vector<float> &weights1, ::std::vector<float> &weights2,
                 float maxdist = 0.0f, int verbose = 0,
                 ICP_xform_type xform_type = ICP_RIGID);

// Easier-to-use interfaces to ICP
extern float ICP(TriMesh *mesh1, TriMesh *mesh2,
                 const xform &xf1, xform &xf2,
                 const KDtree *kd1, const KDtree *kd2,
                 int verbose = 0, ICP_xform_type xform_type = ICP_RIGID);

extern float ICP(TriMesh *mesh1, TriMesh *mesh2,
                 const xform &xf1, xform &xf2,
                 int verbose = 0, ICP_xform_type xform_type = ICP_RIGID);

// Compatibility interfaces to ICP from before we had ICP_xform_type
extern float ICP(TriMesh *mesh1, TriMesh *mesh2,
                 const xform &xf1, xform &xf2,
                 const KDtree *kd1, const KDtree *kd2,
                 ::std::vector<float> &weights1, ::std::vector<float> &weights2,
                 float maxdist, int verbose,
                 bool do_scale, bool do_affine = false);

extern float ICP(TriMesh *mesh1, TriMesh *mesh2, const xform &xf1, xform &xf2,
                 const KDtree *kd1, const KDtree *kd2,
                 int verbose, bool do_scale, bool do_affine = false);

extern float ICP(TriMesh *mesh1, TriMesh *mesh2, const xform &xf1, xform &xf2,
                 int verbose, bool do_scale, bool do_affine = false);

} // namespace trimesh

#endif
