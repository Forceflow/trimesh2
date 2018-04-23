/*
Szymon Rusinkiewicz
Princeton University

hf.cc
Somewhat hackish tool for filling holes in meshes.
Based on initial ideas developed together with Olaf Hall-Holt.
*/

#include "TriMesh.h"
#ifdef _WIN32
# include "wingetopt.h"
#else
# include <unistd.h>
#endif
#include <cstdio>
#include <set>
#include <map>
#include <queue>
#include <iterator>
#include <functional>
using namespace std;
using namespace trimesh;


// An "edge" type and ordering functor
typedef pair<int,int> edge;
class OrderEdge {
public:
	bool operator () (const edge &e1, const edge &e2) const
	{
		if (e1.first < e2.first)
			return true;
		else if (e1.first > e2.first)
			return false;
		else if (e1.second < e2.second)
			return true;
		else
			return false;
	}
};
typedef set<edge,OrderEdge> edgeset;


// A "hole" type - indices of vertices that make up the hole.
typedef vector<int> hole;


// A triple of vertices, together with a quality metric
struct Triple {
	int v1, v2, v3;
	float quality;
	Triple(int _v1, int _v2, int _v3, float _q) :
		v1(_v1), v2(_v2), v3(_v3), quality(_q)
		{}
	bool operator < (const Triple &rhs) const
	{
		return quality < rhs.quality;
	}
};


#define NO_FACE -1
struct FaceStruct {
	int v1, v2, v3;
	int n12, n23, n31;
	FaceStruct(int _v1, int _v2, int _v3,
	           int _n12 = NO_FACE, int _n23 = NO_FACE, int _n31 = NO_FACE) :
			v1(_v1), v2(_v2), v3(_v3),
			n12(_n12), n23(_n23), n31(_n31)
		{}
};

struct VertStruct {
	point p;
	vec norm;
	float &operator [] (int i) { return p[i]; }
	const float &operator [] (int i) const { return p[i]; }
	VertStruct() {}
	VertStruct(const point &_p, const vec &_norm)
	{
		p[0] = _p[0]; p[1] = _p[1]; p[2] = _p[2];
		norm[0] = _norm[0]; norm[1] = _norm[1]; norm[2] = _norm[2];
	}
};


// Find the neighbor field in a FaceStruct that points to i, and change
// it to point to j
void update_neighbor(FaceStruct &f, int i, int j)
{
	if (f.n12 == i)
		f.n12 = j;
	else if (f.n23 == i)
		f.n23 = j;
	else
		f.n31 = j;
}


// Find all the boundary edges
edgeset *find_boundary_edges(const TriMesh *themesh)
{
	printf("Finding boundary edges... "); fflush(stdout);
	edgeset *edges = new edgeset;

	for (size_t f = 0; f < themesh->faces.size(); f++) {
		for (int i = 0; i < 3; i++) {
			int v1 = themesh->faces[f][i];
			int v2 = themesh->faces[f][NEXT_MOD3(i)];

			// Opposite-pointing edges cancel each other
			if (!edges->erase(make_pair(v2,v1)))
				edges->insert(make_pair(v1,v2));
		}
	}

	printf("Done.\n");
	return edges;
}


// Find the initial (before hole-filling) neighbors of all the boundary verts
void find_initial_edge_neighbors(const TriMesh *themesh,
                                 const edgeset *edges,
                                 map< int, set<int> > &initial_edge_neighbors)
{
	size_t nv = themesh->vertices.size(), nf = themesh->faces.size();
	vector<bool> is_edge(nv);
	for (edgeset::const_iterator i = edges->begin(); i != edges->end(); i++) {
		is_edge[i->first] = true;
		is_edge[i->second] = true;
	}

	for (size_t i = 0; i < nf; i++) {
		int v1 = themesh->faces[i][0];
		int v2 = themesh->faces[i][1];
		int v3 = themesh->faces[i][2];
		if (is_edge[v1]) {
			initial_edge_neighbors[v1].insert(v2);
			initial_edge_neighbors[v1].insert(v3);
		}
		if (is_edge[v2]) {
			initial_edge_neighbors[v1].insert(v3);
			initial_edge_neighbors[v1].insert(v1);
		}
		if (is_edge[v3]) {
			initial_edge_neighbors[v1].insert(v1);
			initial_edge_neighbors[v1].insert(v2);
		}
	}
}


// Find a list of holes, given the boundary edges
vector<hole> *find_holes(edgeset *edges)
{
	printf("Finding holes... "); fflush(stdout);
	vector<hole> *holelist = new vector<hole>;
	edgeset badedges;

	while (!edges->empty()) {
		// Find an edge at which to start
		const edge &firstedge = *(edges->begin());
		int firstvert = firstedge.first;
		int lastvert = firstedge.second;
		edges->erase(edges->begin());

		// Add the verts as the first two in a new hole
		holelist->push_back(hole());
		hole &newhole = holelist->back();
		newhole.push_back(firstvert);
		newhole.push_back(lastvert);

		// Follow edges to find the rest of this hole
		while (1) {

			// Find an edge that starts at lastvert
			edgeset::iterator ei = edges->upper_bound(make_pair(lastvert,-1));
			if (ei == edges->end() || ei->first != lastvert) {
				fprintf(stderr, "\nCouldn't find an edge out of vert %d\n", lastvert);
				exit(1);
			}
			int nextvert = ei->second;
			edges->erase(ei);

			// Are we done?
			if (nextvert == firstvert)
				break;

			// Have we encountered this vertex before in this hole?
			// XXX - linear search.  Yuck.
			hole::iterator hi = find(newhole.begin(), newhole.end(),
			                         nextvert);
			if (hi != newhole.end()) {
				// Assuming everything is OK topologically,
				// this could only have been caused if there
				// was a choice of ways to go the last time
				// we encountered this vertex.  Obviously,
				// we chose the wrong way, so find a different
				// way to go.
				edgeset::iterator nei = edges->upper_bound(make_pair(nextvert,-1));
				if (nei == edges->end() || nei->first != nextvert) {
					fprintf(stderr, "\nCouldn't find an edge out of vert %d\n", nextvert);
					exit(1);
				}
				// XXX - for paranoia's sake, we should check
				// that nei->second is not in newhole

				// Put the bad edges into "badedges"
				for (hole::iterator tmp = hi; tmp+1 != newhole.end(); tmp++)
					badedges.insert(make_pair(*tmp, *(tmp+1)));
				badedges.insert(make_pair(lastvert, nextvert));
				newhole.erase(hi+1, newhole.end());

				// Take the new edge, and run with it
				lastvert = nei->second;
				newhole.push_back(lastvert);
				edges->erase(nei);
			} else {
				// All OK.  Add this vert to the hole and go on
				newhole.push_back(nextvert);
				lastvert = nextvert;
			}
		}
		edges->insert(badedges.begin(), badedges.end());
		badedges.clear();
	}

	printf("Done.\n");
	return holelist;
}


// Print out some statistics about the holes
void print_holes(const vector<hole> *holes, bool verbose)
{
	printf("%d holes total, ", (int)holes->size());
	if (!holes->size())
		return;

	vector<int> holesizes;
	for (size_t i = 0; i < holes->size(); i++)
		holesizes.push_back((*holes)[i].size());
	sort(holesizes.begin(), holesizes.end(), greater<int>());

	if (!verbose || holes->size() == 1) {
		printf("Largest is %d\n", holesizes[0]);
		return;
	}

	int nprint = min((int)holes->size(), 10);
	printf("Largest %d are:\n", nprint);
	for (int i = 0; i < nprint; i++)
		printf("    %d\n", holesizes[i]);
}


// Compute a quality metric for a potential triangle with three vertices
inline float quality(const TriMesh *themesh, float meanedgelen,
                     const vector<VertStruct> &newverts,
                     int v1, int v2, int v3,
                     bool hack = false)
{

#define VERT(v) ((size_t(v) < themesh->vertices.size()) ? \
		themesh->vertices[v] : \
		newverts[size_t(v) - themesh->vertices.size()].p)
#define NORM(v) ((size_t(v) < themesh->vertices.size()) ? \
		themesh->normals[v] : \
		newverts[size_t(v) - themesh->vertices.size()].norm)

	if (v1 == v2 || v2 == v3 || v3 == v1)
		return -1000.0f;

	const point &p1 = VERT(v1);
	const point &p2 = VERT(v2);
	const point &p3 = VERT(v3);
	vec side1 = p1 - p2, side2 = p2 - p3, side3 = p3 - p1;

	vec norm = side2 CROSS side1;
	normalize(norm);

	float dot1 = norm DOT NORM(v1);
	float dot2 = norm DOT NORM(v2);
	float dot3 = norm DOT NORM(v3);
	if (dot1 < -0.999f || dot2 < -0.999f || dot3 < -0.999f)
		return -1000.0f;

	float len1 = len(side1);
	float len2 = len(side2);
	float len3 = len(side3);

	float maxedgelen = max(max(len1,len2),len3);
	//float minedgelen = min(min(len1,len2),len3);

	normalize(side1);
	normalize(side2);
	normalize(side3);

	float f = dot1/(1.0f+dot1) + dot2/(1.0f+dot2) + dot3/(1.0f+dot3);

	//float d = meanedgelen/(maxedgelen+meanedgelen);
	float d1 = 1.0f + maxedgelen/meanedgelen;
	//d1 = sqrt(d1);
	//float d = 0;


	float a;
	if (hack) {
		//a = 0.1f*sqr(2.0f+Dot(side1, side2));
		a = 2.0f + (side1 DOT side2);
		//a = sqrt(1.0f/(1.0f - a));
	} else {
		a = sqr(1.0f + (side1 DOT side2)) +
		    sqr(1.0f + (side2 DOT side3)) +
		    sqr(1.0f + (side3 DOT side1));
	}

	return f*sqrt(d1) - a * d1;
}


// Fill the given hole, and add the newly-created triangles to newtris
// XXX - FIXME!  This is O(n^2)
void fill_hole(const TriMesh *themesh, float meanedgelen,
               const vector<VertStruct> &newverts,
               map< int, set<int> > &initial_edge_neighbors,
               const hole &thehole,
               vector<FaceStruct> &newtris)
{
	vector<bool> used(thehole.size(), false);

	priority_queue<Triple> q;
	for (size_t i = 0; i < thehole.size(); i++) {
		int j = (i+1) % thehole.size();
		int k = (j+1) % thehole.size();
		float qual = quality(themesh, meanedgelen, newverts,
		                     thehole[i], thehole[j], thehole[k], true);
		if (initial_edge_neighbors[thehole[i]].find(thehole[k]) !=
		    initial_edge_neighbors[thehole[i]].end())
			qual = -1000.0f;
		q.push(Triple(i, j, k, qual));
	}

	while (!q.empty()) {
		// Take the highest-quality triple off the queue
		const Triple next = q.top();
		q.pop();

		// Ignore triangles referencing already-used verts
		if (!used[next.v1] && !used[next.v3]) {
			used[next.v2] = true;

			// Create the new face and push it onto newtris
			newtris.push_back(FaceStruct(thehole[next.v3], thehole[next.v2], thehole[next.v1]));

			// Find next verts forward and back
			int forw = next.v3;
			do {
				forw++;
				forw %= thehole.size();
			} while (used[forw]);
			if (forw == next.v1)
				return;
			int back = next.v1;
			do {
				back--;
				if (back < 0) back += thehole.size();
			} while (used[back]);

			// Insert potential new triangles
			float q13f = quality(themesh, meanedgelen, newverts,
			                     thehole[next.v1],
			                     thehole[next.v3],
			                     thehole[forw],
			                     true);
			if (initial_edge_neighbors[thehole[next.v1]].find(thehole[forw]) !=
			    initial_edge_neighbors[thehole[next.v1]].end())
				q13f = -2000.0f;
			float qb13 = quality(themesh, meanedgelen, newverts,
			                     thehole[back],
			                     thehole[next.v1],
			                     thehole[next.v3],
			                     true);
			if (initial_edge_neighbors[thehole[back]].find(thehole[next.v3]) !=
			    initial_edge_neighbors[thehole[back]].end())
				qb13 = -2000.0f;
			q.push(Triple(next.v1, next.v3, forw, q13f));
			q.push(Triple(back, next.v1, next.v3, qb13));
		}
	}
}


// Connect up neighbors in a patch of surface
// XXX - O(n^2)
void find_neighbors(vector<FaceStruct> &tris)
{
	int n = tris.size();
	for (int i = 0; i < n; i++) {
		int v1 = tris[i].v1, v2 = tris[i].v2, v3 = tris[i].v3;
		int found = 0;
		for (int j=0; j < n; j++) {
			if (i == j) continue;
			if ((v1 == tris[j].v2 && v2 == tris[j].v1) ||
			    (v1 == tris[j].v3 && v2 == tris[j].v2) ||
			    (v1 == tris[j].v1 && v2 == tris[j].v3)) {
				tris[i].n12 = j;
				found++; if (found == 3) break;
			} else if ((v2 == tris[j].v2 && v3 == tris[j].v1) ||
			           (v2 == tris[j].v3 && v3 == tris[j].v2) ||
			           (v2 == tris[j].v1 && v3 == tris[j].v3)) {
				tris[i].n23 = j;
				found++; if (found == 3) break;
			} else if ((v3 == tris[j].v2 && v1 == tris[j].v1) ||
			           (v3 == tris[j].v3 && v1 == tris[j].v2) ||
			           (v3 == tris[j].v1 && v1 == tris[j].v3)) {
				tris[i].n31 = j;
				found++; if (found == 3) break;
			}
		}
	}
}


// Should we flip (v1,v2,v3) and (v1,v3,v4) to (v1,v2,v4) and (v2,v3,v4)?
bool should_flip(const TriMesh *themesh, float meanedgelen,
                 const vector<VertStruct> &newverts,
                 int v1, int v2, int v3, int v4)
{
	// Step 1 - Sanity check
	if (v2 == v4)
		return false;
	if (v1 == v3) // Ugh.  Hopefully we'll never get here...
		return true;

	// Step 2 - Check the normals of the two resultant polygons.  If it
	// gets much worse (a crease), shortcut out
	const point &p1 = VERT(v1);
	const point &p2 = VERT(v2);
	const point &p3 = VERT(v3);
	const point &p4 = VERT(v4);
	vec norm123 = trinorm(p1, p2, p3);  normalize(norm123);
	vec norm134 = trinorm(p1, p3, p4);  normalize(norm134);
	vec norm124 = trinorm(p1, p2, p4);  normalize(norm124);
	vec norm234 = trinorm(p2, p3, p4);  normalize(norm234);
	float old_dot = max(-0.999f, (norm123 DOT norm134));
	float new_dot = max(-0.999f, (norm124 DOT norm234));
	float old_badness = 1.0f - old_dot/(1.0f+old_dot);
	float new_badness = 1.0f - new_dot/(1.0f+new_dot);
	if (new_badness > 1.5f*old_badness)
		return false;
	if (new_badness < 0.5f*old_badness)
		return true;

	// Step 3 - Evaluate quality()
	float q123 = quality(themesh, meanedgelen, newverts, v3, v2, v1);
	float q134 = quality(themesh, meanedgelen, newverts, v4, v3, v1);
	float q124 = quality(themesh, meanedgelen, newverts, v4, v2, v1);
	float q234 = quality(themesh, meanedgelen, newverts, v4, v3, v2);

	return (q124+q234 > q123+q134);
}

// Improve the triangulation of a hole by performing a bunch of edge flips
void improve_triangulation(const TriMesh *themesh, float meanedgelen,
                           const vector<VertStruct> &newverts,
                           map< int, set<int> > &initial_edge_neighbors,
                           vector<FaceStruct> &tris)
{
	int n = tris.size();
	if (!n)
		return;

	int nv = themesh->vertices.size();
	edgeset edges;
	set<int> edgeverts;
	for (size_t i = 0; i < tris.size(); i++) {
		edges.insert(make_pair(tris[i].v1, tris[i].v2));
		edges.insert(make_pair(tris[i].v1, tris[i].v3));
		edges.insert(make_pair(tris[i].v2, tris[i].v3));
		edges.insert(make_pair(tris[i].v2, tris[i].v1));
		edges.insert(make_pair(tris[i].v3, tris[i].v1));
		edges.insert(make_pair(tris[i].v3, tris[i].v2));
		if (tris[i].v1 < nv)
			edgeverts.insert(tris[i].v1);
		if (tris[i].v2 < nv)
			edgeverts.insert(tris[i].v2);
		if (tris[i].v3 < nv)
			edgeverts.insert(tris[i].v3);
	}
	for (set<int>::const_iterator i = edgeverts.begin();
	     i != edgeverts.end();
	     i++) {
		for (set<int>::const_iterator j = initial_edge_neighbors[*i].begin();
		     j != initial_edge_neighbors[*i].end();
		     j++) {
			edges.insert(make_pair(*i, *j));
			edges.insert(make_pair(*j, *i));
		}
	}

	for (int iter = 0; iter < 2*n; iter++) {
		//int i = iter % n;
		int i = uniform_rnd(n);
		if (tris[i].n12 != NO_FACE) {
			int j = tris[i].n12;
			int v1 = tris[i].v1, v2 = tris[i].v2, v3 = tris[i].v3;
			int v4 = tris[j].v1, v5 = tris[j].v2, v6 = tris[j].v3;
			int match1 = (v4 == v1) ? 1 : (v5 == v1) ? 2 : 3;
			int v0 = (match1 == 1) ? v5 : (match1 == 2) ? v6 : v4;
			if (should_flip(themesh, meanedgelen, newverts,
					v1, v0, v2, v3)) {
				if (match1 == 1 &&
				    edges.find(make_pair(tris[i].v3, tris[j].v2)) == edges.end()) {
					edges.erase(make_pair(tris[i].v1, tris[i].v2));
					edges.erase(make_pair(tris[i].v2, tris[i].v1));
					edges.insert(make_pair(tris[i].v3, tris[j].v2));
					edges.insert(make_pair(tris[j].v2, tris[i].v3));
					tris[i].v2 = v0;
					tris[j].v1 = v3;
					int neighbor = tris[j].n12;
					tris[i].n12 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], j, i);
					neighbor = tris[i].n23;
					tris[j].n31 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], i, j);
					tris[i].n23 = j;
					tris[j].n12 = i;
				} else if (match1 == 2 &&
				           edges.find(make_pair(tris[i].v3, tris[j].v3)) == edges.end()) {
					edges.erase(make_pair(tris[i].v1, tris[i].v2));
					edges.erase(make_pair(tris[i].v2, tris[i].v1));
					edges.insert(make_pair(tris[i].v3, tris[j].v3));
					edges.insert(make_pair(tris[j].v3, tris[i].v3));
					tris[i].v2 = v0;
					tris[j].v2 = v3;
					int neighbor = tris[j].n23;
					tris[i].n12 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], j, i);
					neighbor = tris[i].n23;
					tris[j].n12 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], i, j);
					tris[i].n23 = j;
					tris[j].n23 = i;
				} else if (match1 == 3 &&
				           edges.find(make_pair(tris[i].v3, tris[j].v1)) == edges.end()) {
					edges.erase(make_pair(tris[i].v1, tris[i].v2));
					edges.erase(make_pair(tris[i].v2, tris[i].v1));
					edges.insert(make_pair(tris[i].v3, tris[j].v1));
					edges.insert(make_pair(tris[j].v1, tris[i].v3));
					tris[i].v2 = v0;
					tris[j].v3 = v3;
					int neighbor = tris[j].n31;
					tris[i].n12 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], j, i);
					neighbor = tris[i].n23;
					tris[j].n23 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], i, j);
					tris[i].n23 = j;
					tris[j].n31 = i;
				}
			}
		}
		if (tris[i].n23 != NO_FACE) {
			int j = tris[i].n23;
			int v1 = tris[i].v1, v2 = tris[i].v2, v3 = tris[i].v3;
			int v4 = tris[j].v1, v5 = tris[j].v2, v6 = tris[j].v3;
			int match2 = (v4 == v2) ? 1 : (v5 == v2) ? 2 : 3;
			int v0 = (match2 == 1) ? v5 : (match2 == 2) ? v6 : v4;
			if (should_flip(themesh, meanedgelen, newverts,
					v2, v0, v3, v1)) {
				if (match2 == 1 &&
				    edges.find(make_pair(tris[i].v1, tris[j].v2)) == edges.end()) {
					edges.erase(make_pair(tris[i].v2, tris[i].v3));
					edges.erase(make_pair(tris[i].v3, tris[i].v2));
					edges.insert(make_pair(tris[i].v1, tris[j].v2));
					edges.insert(make_pair(tris[j].v2, tris[i].v1));
					tris[i].v3 = v0;
					tris[j].v1 = v1;
					int neighbor = tris[j].n12;
					tris[i].n23 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], j, i);
					neighbor = tris[i].n31;
					tris[j].n31 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], i, j);
					tris[i].n31 = j;
					tris[j].n12 = i;
				} else if (match2 == 2 &&
				           edges.find(make_pair(tris[i].v1, tris[j].v3)) == edges.end()) {
					edges.erase(make_pair(tris[i].v2, tris[i].v3));
					edges.erase(make_pair(tris[i].v3, tris[i].v2));
					edges.insert(make_pair(tris[i].v1, tris[j].v3));
					edges.insert(make_pair(tris[j].v3, tris[i].v1));
					tris[i].v3 = v0;
					tris[j].v2 = v1;
					int neighbor = tris[j].n23;
					tris[i].n23 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], j, i);
					neighbor = tris[i].n31;
					tris[j].n12 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], i, j);
					tris[i].n31 = j;
					tris[j].n23 = i;
				} else if (match2 == 3 &&
				           edges.find(make_pair(tris[i].v1, tris[j].v1)) == edges.end()) {
					edges.erase(make_pair(tris[i].v2, tris[i].v2));
					edges.erase(make_pair(tris[i].v3, tris[i].v3));
					edges.insert(make_pair(tris[i].v1, tris[j].v1));
					edges.insert(make_pair(tris[j].v1, tris[i].v1));
					tris[i].v3 = v0;
					tris[j].v3 = v1;
					int neighbor = tris[j].n31;
					tris[i].n23 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], j, i);
					neighbor = tris[i].n31;
					tris[j].n23 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], i, j);
					tris[i].n31 = j;
					tris[j].n31 = i;
				}
			}
		}
		if (tris[i].n31 != NO_FACE) {
			int j = tris[i].n31;
			int v1 = tris[i].v1, v2 = tris[i].v2, v3 = tris[i].v3;
			int v4 = tris[j].v1, v5 = tris[j].v2, v6 = tris[j].v3;
			int match3 = (v4 == v3) ? 1 : (v5 == v3) ? 2 : 3;
			int v0 = (match3 == 1) ? v5 : (match3 == 2) ? v6 : v4;
			if (should_flip(themesh, meanedgelen, newverts,
					v3, v0, v1, v2)) {
				if (match3 == 1 &&
				    edges.find(make_pair(tris[i].v2, tris[j].v2)) == edges.end()) {
					edges.erase(make_pair(tris[i].v3, tris[i].v1));
					edges.erase(make_pair(tris[i].v1, tris[i].v3));
					edges.insert(make_pair(tris[i].v2, tris[j].v2));
					edges.insert(make_pair(tris[j].v2, tris[i].v2));
					tris[i].v1 = v0;
					tris[j].v1 = v2;
					int neighbor = tris[j].n12;
					tris[i].n31 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], j, i);
					neighbor = tris[i].n12;
					tris[j].n31 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], i, j);
					tris[i].n12 = j;
					tris[j].n12 = i;
				} else if (match3 == 2 &&
				           edges.find(make_pair(tris[i].v2, tris[j].v3)) == edges.end()) {
					edges.erase(make_pair(tris[i].v3, tris[i].v1));
					edges.erase(make_pair(tris[i].v1, tris[i].v3));
					edges.insert(make_pair(tris[i].v2, tris[j].v3));
					edges.insert(make_pair(tris[j].v3, tris[i].v2));
					tris[i].v1 = v0;
					tris[j].v2 = v2;
					int neighbor = tris[j].n23;
					tris[i].n31 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], j, i);
					neighbor = tris[i].n12;
					tris[j].n12 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], i, j);
					tris[i].n12 = j;
					tris[j].n23 = i;
				} else if (match3 == 3 &&
				           edges.find(make_pair(tris[i].v2, tris[j].v1)) == edges.end()) {
					edges.erase(make_pair(tris[i].v3, tris[i].v1));
					edges.erase(make_pair(tris[i].v1, tris[i].v3));
					edges.insert(make_pair(tris[i].v2, tris[j].v1));
					edges.insert(make_pair(tris[j].v1, tris[i].v2));
					tris[i].v1 = v0;
					tris[j].v3 = v2;
					int neighbor = tris[j].n31;
					tris[i].n31 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], j, i);
					neighbor = tris[i].n12;
					tris[j].n23 = neighbor;
					if (neighbor != NO_FACE)
						update_neighbor(tris[neighbor], i, j);
					tris[i].n12 = j;
					tris[j].n31 = i;
				}
			}
		}
	}
}


// Subdivide triangles, such that no edge is longer than subdiv_len
// The faces in tris must have correct neighbor information (as produced by
// find_neighbors()).  Works by subdividing long edges...
void tesselate_holefill(const TriMesh *themesh,
			float subdiv_len,
			vector<VertStruct> &newverts,
			vector<FaceStruct> &tris)
{
	int oldverts = themesh->vertices.size();

	int s = tris.size();
	for (int i = 0; i < s; i++) {
		const point &v1 = VERT(tris[i].v1);
		const point &v2 = VERT(tris[i].v2);
		const point &v3 = VERT(tris[i].v3);
		float d12 = dist(v1, v2);
		float d23 = dist(v2, v3);
		float d31 = dist(v3, v1);
		float dmax = max(max(d12, d23), d31);

		if (dmax < subdiv_len) {
			continue;
		} else if (tris[i].n12 != NO_FACE && d12 == dmax) {

			// Subdivide edge 1->2
			int j = tris[i].n12;
			int corner = (tris[j].v1 == tris[i].v1 ? tris[j].v2 :
			              tris[j].v2 == tris[i].v1 ? tris[j].v3 :
			              tris[j].v1);
			int neighbor = (tris[j].v1 == tris[i].v1 ? tris[j].n23 :
					tris[j].v2 == tris[i].v1 ? tris[j].n31 :
					tris[j].n12);
			point newpoint = 0.5f * (v1 + v2);
			vec norm1 = trinorm(v1, v2, v3);
			vec norm2 = trinorm(v1, VERT(corner), v2);
			vec newnorm = norm1 + norm2;
			normalize(newnorm);
			int new_ind = newverts.size() + oldverts;
			newverts.push_back(VertStruct(newpoint, newnorm));
			int new_tri1 = tris.size();
			tris.push_back(FaceStruct(new_ind,
			                          tris[i].v2,
			                          tris[i].v3,
			                          tris[i].n12,
			                          tris[i].n23,
			                          i));
			tris[i].v2 = new_ind;
			tris[i].n23 = new_tri1;
			int new_tri2 = tris.size();
			tris.push_back(FaceStruct(new_ind,
			                          corner,
			                          tris[new_tri1].v2,
			                          j,
			                          neighbor,
			                          new_tri1));
			update_neighbor(tris[j], neighbor, new_tri2);
			update_neighbor(tris[new_tri1], j, new_tri2);
			if (tris[j].v1 == tris[new_tri1].v2)
				tris[j].v1 = new_ind;
			else if (tris[j].v2 == tris[new_tri1].v2)
				tris[j].v2 = new_ind;
			else
				tris[j].v3 = new_ind;
			if (tris[new_tri1].n23 != NO_FACE)
				update_neighbor(tris[tris[new_tri1].n23], i, new_tri1);
			if (neighbor != NO_FACE)
				update_neighbor(tris[neighbor], j, new_tri2);

		} else if (tris[i].n23 != NO_FACE && d23 == dmax) {

			// Subdivide edge 2->3
			int j = tris[i].n23;
			int corner = (tris[j].v1 == tris[i].v2 ? tris[j].v2 :
			              tris[j].v2 == tris[i].v2 ? tris[j].v3 :
			              tris[j].v1);
			int neighbor = (tris[j].v1 == tris[i].v2 ? tris[j].n23 :
					tris[j].v2 == tris[i].v2 ? tris[j].n31 :
					tris[j].n12);
			point newpoint = 0.5f * (v2 + v3);
			vec norm1 = trinorm(v1, v2, v3);
			vec norm2 = trinorm(v2, VERT(corner), v3);
			vec newnorm = norm1 + norm2;
			normalize(newnorm);
			int new_ind = newverts.size() + oldverts;
			newverts.push_back(VertStruct(newpoint, newnorm));
			int new_tri1 = tris.size();
			tris.push_back(FaceStruct(new_ind,
			                          tris[i].v3,
			                          tris[i].v1,
			                          tris[i].n23,
			                          tris[i].n31,
			                          i));
			tris[i].v3 = new_ind;
			tris[i].n31 = new_tri1;
			int new_tri2 = tris.size();
			tris.push_back(FaceStruct(new_ind,
			                          corner,
			                          tris[new_tri1].v2,
			                          j,
			                          neighbor,
			                          new_tri1));
			update_neighbor(tris[j], neighbor, new_tri2);
			update_neighbor(tris[new_tri1], j, new_tri2);
			if (tris[j].v1 == tris[new_tri1].v2)
				tris[j].v1 = new_ind;
			else if (tris[j].v2 == tris[new_tri1].v2)
				tris[j].v2 = new_ind;
			else
				tris[j].v3 = new_ind;
			if (tris[new_tri1].n23 != NO_FACE)
				update_neighbor(tris[tris[new_tri1].n23], i, new_tri1);
			if (neighbor != NO_FACE)
				update_neighbor(tris[neighbor], j, new_tri2);

		} else if (tris[i].n31 != NO_FACE && d31 == dmax) {

			// Subdivide edge 3->1
			int j = tris[i].n31;
			int corner = (tris[j].v1 == tris[i].v3 ? tris[j].v2 :
			              tris[j].v2 == tris[i].v3 ? tris[j].v3 :
			              tris[j].v1);
			int neighbor = (tris[j].v1 == tris[i].v3 ? tris[j].n23 :
					tris[j].v2 == tris[i].v3 ? tris[j].n31 :
					tris[j].n12);
			point newpoint = 0.5f * (v3 + v1);
			vec norm1 = trinorm(v1, v2, v3);
			vec norm2 = trinorm(v3, VERT(corner), v1);
			vec newnorm = norm1 + norm2;
			normalize(newnorm);
			int new_ind = newverts.size() + oldverts;
			newverts.push_back(VertStruct(newpoint, newnorm));
			int new_tri1 = tris.size();
			tris.push_back(FaceStruct(new_ind,
			                          tris[i].v1,
			                          tris[i].v2,
			                          tris[i].n31,
			                          tris[i].n12,
			                          i));
			tris[i].v1 = new_ind;
			tris[i].n12 = new_tri1;
			int new_tri2 = tris.size();
			tris.push_back(FaceStruct(new_ind,
			                          corner,
			                          tris[new_tri1].v2,
			                          j,
			                          neighbor,
			                          new_tri1));
			update_neighbor(tris[j], neighbor, new_tri2);
			update_neighbor(tris[new_tri1], j, new_tri2);
			if (tris[j].v1 == tris[new_tri1].v2)
				tris[j].v1 = new_ind;
			else if (tris[j].v2 == tris[new_tri1].v2)
				tris[j].v2 = new_ind;
			else
				tris[j].v3 = new_ind;
			if (tris[new_tri1].n23 != NO_FACE)
				update_neighbor(tris[tris[new_tri1].n23], i, new_tri1);
			if (neighbor != NO_FACE)
				update_neighbor(tris[neighbor], j, new_tri2);

		}
	}
}


// Relax a triangulation of a hole.  Moves vertices firstvert..newverts.size()-1
void relax_triangulation(TriMesh *themesh,
                         float smooth, float csmooth, int iterations,
                         vector<VertStruct> &newverts,
                         int firstvert,
                         const vector<FaceStruct> &tris)
{
	int firstvert_g = firstvert + themesh->vertices.size();
	int nverts = newverts.size() - firstvert;
	if (nverts == 0)
		return;

	vector<vec> m(nverts);
	vector<float> num(nverts);

	// For each iteration...
	for (int iter = 0; iter < iterations; iter++) {

		// Initialize motion
		for (int i = 0; i < nverts; i++) {
			m[i] = vec();
			num[i] = 0;
		}

		// For each triangle...
		for (size_t i = 0; i < tris.size(); i++) {

			const point &v1 = VERT(tris[i].v1);
			const point &v2 = VERT(tris[i].v2);
			const point &v3 = VERT(tris[i].v3);
			vec edge12 = v2 - v1;
			vec edge23 = v3 - v2;
			vec edge31 = v1 - v3;

			// Tug each vertex towards neighbors
			if (tris[i].v1 >= firstvert_g) {
				m[tris[i].v1 - firstvert_g] += edge12 - edge31;
				num[tris[i].v1 - firstvert_g]++;
			}
			if (tris[i].v2 >= firstvert_g) {
				m[tris[i].v2 - firstvert_g] += edge23 - edge12;
				num[tris[i].v2 - firstvert_g]++;
			}
			if (tris[i].v3 >= firstvert_g) {
				m[tris[i].v3 - firstvert_g] += edge31 - edge23;
				num[tris[i].v3 - firstvert_g]++;
			}

			// Curvature...
			vec mynorm = trinorm(v1, v2, v3);
			normalize(mynorm);
			if (tris[i].v3 >= firstvert_g) {
				vec othernorm;
				if (tris[i].n12 == NO_FACE) {
					othernorm = NORM(tris[i].v1) +
					            NORM(tris[i].v2);
				} else {
					int j = tris[i].n12;
					othernorm = trinorm(
						VERT(tris[j].v1),
						VERT(tris[j].v2),
						VERT(tris[j].v3));
				}
				normalize(othernorm);
				vec cross = mynorm CROSS othernorm;
				float alpha = -(edge12 DOT edge31) / len2(edge12);
				float beta = 1.0f - alpha;
				point p = alpha * v2 + beta * v1;
				float l = csmooth * dist(v3, p);
				if ((mynorm DOT othernorm) > 0.0f)
					l *= len(cross);
				if ((cross DOT edge12) < 0.0f)
					l = -l;
				m[tris[i].v3 - firstvert_g] += l*mynorm;
			}
			if (tris[i].v1 >= firstvert_g) {
				vec othernorm;
				if (tris[i].n23 == NO_FACE) {
					othernorm = NORM(tris[i].v2) +
					            NORM(tris[i].v3);
				} else {
					int j = tris[i].n23;
					othernorm = trinorm(
						VERT(tris[j].v1),
						VERT(tris[j].v2),
						VERT(tris[j].v3));
				}
				normalize(othernorm);
				vec cross = mynorm CROSS othernorm;
				float alpha = -(edge23 DOT edge12) / len2(edge23);
				float beta = 1.0f - alpha;
				point p = alpha * v3 + beta * v2;
				float l = csmooth * dist(v1, p);
				if ((mynorm DOT othernorm) > 0.0f)
					l *= len(cross);
				if ((cross DOT edge23) < 0.0f)
					l = -l;
				m[tris[i].v1 - firstvert_g] += l*mynorm;
			}
			if (tris[i].v2 >= firstvert_g) {
				vec othernorm;
				if (tris[i].n31 == NO_FACE) {
					othernorm = NORM(tris[i].v3) +
					            NORM(tris[i].v1);
				} else {
					int j = tris[i].n31;
					othernorm = trinorm(
						VERT(tris[j].v1),
						VERT(tris[j].v2),
						VERT(tris[j].v3));
				}
				normalize(othernorm);
				vec cross = mynorm CROSS othernorm;
				float alpha = -(edge31 DOT edge23) / len2(edge31);
				float beta = 1.0f - alpha;
				point p = alpha * v1 + beta * v3;
				float l = csmooth * dist(v2, p);
				if ((mynorm DOT othernorm) > 0.0f)
					l *= len(cross);
				if ((cross DOT edge31) < 0.0f)
					l = -l;
				m[tris[i].v2 - firstvert_g] += l*mynorm;
			}

		}

		// Now actually go and move all the verts...
		for (int i = 0; i < nverts; i++) {
			if (!num[i])
				continue;
			point &me = VERT(firstvert_g+i);
			float scale = smooth / (float) num[i];
			me += scale * m[i];
		}

	}
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [-l | -m num] in.ply [out.ply]\n", myname);
	fprintf(stderr, "	-l	Just list the sizes of the largest holes\n");
	fprintf(stderr, "	-m num	Only fill holes smaller than num\n");
	fprintf(stderr, "	-f	Just fill the holes, don't edgeflip / relax\n");
	fprintf(stderr, "	-c	Don't consider curvature when relaxing\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	// Parse command-line options
	bool fillonly = false;
	bool listonly = false;
	bool nocurv = false;
	int maxfill = 0;

	int c;
	while ((c = getopt(argc, argv, "hcflm:")) != EOF) {
		switch (c) {
			case 'f': fillonly = true; break;
			case 'l': listonly = true; break;
			case 'c': nocurv = true; break;
			case 'm': maxfill = atoi(optarg); break;
			default: usage(argv[0]);
		}
	}
	if (optind >= argc)
		usage(argv[0]);
	const char *infilename = argv[optind];
	optind++;
	const char *outfilename = NULL;
	if (!listonly) {
		if (optind >= argc)
			usage(argv[0]);
		outfilename = argv[optind];
	}


	TriMesh *themesh = TriMesh::read(infilename);
	if (!themesh) {
		fprintf(stderr, "Couldn't read file %s\n", infilename);
		exit(1);
	}
	themesh->need_faces();
	bool had_tstrips = !themesh->tstrips.empty();
	themesh->tstrips.clear();

	edgeset *edges = find_boundary_edges(themesh);

	map< int, set<int> > initial_edge_neighbors;
	find_initial_edge_neighbors(themesh, edges, initial_edge_neighbors);

	vector<hole> *holes = find_holes(edges);
	delete edges;

	if (listonly) {
		print_holes(holes, true);
		exit(0);
	} else {
		print_holes(holes, false);
	}

	themesh->need_normals();
	float mel = themesh->feature_size();

	vector<FaceStruct> newtris;
	vector<VertStruct> newverts;
	for (size_t i = 0; i < holes->size(); i++) {
		if (maxfill && int((*holes)[i].size()) > maxfill)
			continue;
		printf("\r                                                                             \r");
		printf("Filling hole %d of %d, size %d: ", (int)i+1, (int)holes->size(), (int)(*holes)[i].size()); fflush(stdout);
		int firstnewvert = newverts.size();
		vector<FaceStruct> tmptris;
		printf("f"); fflush(stdout);
		fill_hole(themesh, mel, newverts,
		          initial_edge_neighbors,
		          (*holes)[i], tmptris);
		if (!fillonly) {
			printf("n"); fflush(stdout);
			find_neighbors(tmptris);
			printf("i"); fflush(stdout);
			improve_triangulation(themesh, mel, newverts, initial_edge_neighbors, tmptris);
			int niters = int(ceil(2.0f*log((float)tmptris.size())*M_LOG2Ef - 3.0f));
			niters = min(niters, 4);
			for (int j=0; j < niters; j++) {
				printf("t"); fflush(stdout);
				tesselate_holefill(themesh, 2.0f*mel, newverts, tmptris);
				printf("i"); fflush(stdout);
				improve_triangulation(themesh, mel, newverts, initial_edge_neighbors, tmptris);
				printf("r"); fflush(stdout);
				relax_triangulation(themesh,
				                    0.2f, (nocurv ? 0 : 2), 10,
				                    newverts, firstnewvert,
				                    tmptris);
				//improve_triangulation(themesh, mel, newverts, initial_edge_neighbors, tmptris);
			}

			printf("t"); fflush(stdout);
			tesselate_holefill(themesh, 1.5f*mel, newverts, tmptris);
			tesselate_holefill(themesh, mel, newverts, tmptris);
			printf("i"); fflush(stdout);
			improve_triangulation(themesh, mel, newverts, initial_edge_neighbors, tmptris);
			printf("r"); fflush(stdout);
			relax_triangulation(themesh,
			                    0.2f, (nocurv ? 0 : 2), 50,
			                    newverts, firstnewvert,
			                    tmptris);
		}

		copy(tmptris.begin(), tmptris.end(), back_inserter(newtris));
	}
	printf("\r                                                                             \r");
	printf("Filling holes... Done.\n");

	delete holes;
	themesh->normals.clear();

	themesh->faces.reserve(themesh->faces.size() + newtris.size());
	for (size_t i = 0; i < newtris.size(); i++) {
		themesh->faces.push_back(TriMesh::Face(
			newtris[i].v1, newtris[i].v2, newtris[i].v3));
	}

	themesh->vertices.reserve(themesh->vertices.size() + newverts.size());
	for (size_t i = 0; i < newverts.size(); i++) {
		themesh->vertices.push_back(newverts[i].p);
	}

	size_t nv = themesh->vertices.size();
	if (!themesh->colors.empty())
		themesh->colors.resize(nv);
	if (!themesh->confidences.empty())
		themesh->confidences.resize(nv);

	if (had_tstrips)
		themesh->need_tstrips();
	themesh->write(outfilename);
}
