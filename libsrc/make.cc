/*
Szymon Rusinkiewicz
Princeton University

make.cc
Routines for creating TriMeshes.
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include "noise3d.h"
using namespace std;


namespace trimesh {

#include "shape_data.inc"


// Local utility functions
static inline void mkpoint(TriMesh *mesh, float x, float y, float z)
{
	mesh->vertices.push_back(point(x,y,z));
}

static inline void mkface(TriMesh *mesh, int v1, int v2, int v3)
{
	if (v1 == v2 || v2 == v3 || v3 == v1)
		return;
	mesh->faces.push_back(TriMesh::Face(v1, v2, v3));
}

static inline void mkquad(TriMesh *mesh, int ll, int lr, int ul, int ur)
{
	mkface(mesh, ll, lr, ur);
	mkface(mesh, ll, ur, ul);
}


// Tessellated square (-1..1, -1..1, 0)
TriMesh *make_plane(int tess_x, int tess_y)
{
	if (tess_x < 1)
		tess_x = 1;
	if (tess_y < 1)
		tess_y = tess_x;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve((tess_x+1) * (tess_y+1));
	for (int j = 0; j < tess_y+1; j++) {
		float y = -1.0f + 2.0f * j / tess_y;
		for (int i = 0; i < tess_x+1; i++) {
			float x = -1.0f + 2.0f * i / tess_x;
			mkpoint(mesh, x, y, 0);
		}
	}

	mesh->faces.reserve(2*tess_x*tess_y);
	for (int j = 0; j < tess_y; j++) {
		for (int i = 0; i < tess_x; i++) {
			int ind = i + j * (tess_x+1);
			mkquad(mesh, ind, ind+1, ind+tess_x+1, ind+tess_x+2);
		}
	}

	return mesh;
}


// Gaussian bump of height 1 and width sigma
TriMesh *make_bump(int tess, float sigma)
{
	if (tess < 2)
		tess = 2;

	TriMesh *mesh = make_plane(tess, tess);
	int nv = mesh->vertices.size();
	float escale = -0.5f * sqr(1.0f / sigma);
	for (int i = 0; i < nv; i++) {
		point &p = mesh->vertices[i];
		p[2] = exp(escale * (sqr(p[0]) + sqr(p[1])));
	}
	return mesh;
}


// Sine wave of angular frequency omega
TriMesh *make_wave(int tess, float omega)
{
	if (tess < 2)
		tess = 2;

	TriMesh *mesh = make_plane(tess, tess);
	int nv = mesh->vertices.size();
	float scale = 1.0f / omega;
	for (int i = 0; i < nv; i++) {
		point &p = mesh->vertices[i];
		p[2] = scale * sin(omega * p[0]) * sin(omega * p[1]);
	}
	return mesh;
}


// Fractal landscape
TriMesh *make_frac(int tess)
{
	if (tess < 2)
		tess = 2;

	TriMesh *mesh = make_plane(tess, tess);
	int nv = mesh->vertices.size();
	PerlinNoise3D noise(tess, tess, tess);
	for (int i = 0; i < nv; i++) {
		point &p = mesh->vertices[i];
		p[2] = noise.lookup(p[0], p[1], 0.5f);
	}
	return mesh;
}


// Tessellated cube
TriMesh *make_cube(int tess)
{
	if (tess < 1)
		tess = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(6*sqr(tess)+2);
	for (int j = 0; j < tess+1; j++) {
		float y = 1.0f - 2.0f * j / tess;
		for (int i = 0; i < tess+1; i++) {
			float x = 1.0f - 2.0f * i / tess;
			mkpoint(mesh, x, y, -1);
		}
	}
	for (int j = 1; j < tess; j++) {
		float z = -1.0f + 2.0f * j / tess;
		for (int i = 0; i < tess; i++) {
			float x = -1.0f + 2.0f * i / tess;
			mkpoint(mesh, x, -1, z);
		}
		for (int i = 0; i < tess; i++) {
			float y = -1.0f + 2.0f * i / tess;
			mkpoint(mesh, 1, y, z);
		}
		for (int i = 0; i < tess; i++) {
			float x = 1.0f - 2.0f * i / tess;
			mkpoint(mesh, x, 1, z);
		}
		for (int i = 0; i < tess; i++) {
			float y = 1.0f - 2.0f * i / tess;
			mkpoint(mesh, -1, y, z);
		}
	}
	for (int j = 0; j < tess+1; j++) {
		float y = -1.0f + 2.0f * j / tess;
		for (int i = 0; i < tess+1; i++) {
			float x = -1.0f + 2.0f * i / tess;
			mkpoint(mesh, x, y, 1);
		}
	}

	mesh->faces.reserve(12*sqr(tess));
	for (int j = 0; j < tess; j++) {
		for (int i = 0; i < tess; i++) {
			int ind = i + j * (tess+1);
			mkquad(mesh, ind, ind+tess+1, ind+1, ind+tess+2);
		}
	}

	int topstart = sqr(tess+1) + 4*tess*(tess-1);
	for (int j = 0; j < tess; j++) {
		int next = sqr(tess+1) + 4*tess*(j-1);
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = sqr(tess+1)-1 - i;
				lr = ll - 1;
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = tess*(tess+1) - i*(tess+1);
				lr = ll - (tess+1);
			}
			if (j == tess-1) {
				ul = topstart + tess + i*(tess+1);
				ur = ul + (tess+1);
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = i;
				lr = i + 1;
			}
			if (j == tess-1) {
				ul = topstart + sqr(tess+1)-1 - i;
				ur = ul - 1;
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = tess + i*(tess+1);
				lr = ll + (tess+1);
			}
			if (j == tess-1) {
				ul = topstart + tess*(tess+1) - i*(tess+1);
				ur = ul - (tess+1);
			}
			if (i == tess-1) {
				if (j != 0)
					lr -= 4*tess;
				if (j != tess-1)
					ur -= 4*tess;
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
	}
	for (int j = 0; j < tess; j++) {
		for (int i = 0; i < tess; i++) {
			int ind = topstart + i + j * (tess+1);
			mkquad(mesh, ind, ind+1, ind+tess+1, ind+tess+2);
		}
	}

	return mesh;
}


// Disc with the given tessellation in angle and radius
TriMesh *make_disc(int tess_th, int tess_r)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_r < 1)
		tess_r = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(1+tess_th*tess_r);
	mkpoint(mesh, 0, 0, 0);
	for (int j = 1; j <= tess_r; j++) {
		float r = (float) j / tess_r;
		for (int i = 0; i < tess_th; i++) {
			float th = M_2PIf * i / tess_th;
			mkpoint(mesh, r*cos(th), r*sin(th), 0);
		}
	}

	mesh->faces.reserve(2*tess_th*tess_r-tess_th);
	for (int i = 0; i < tess_th; i++)
		mkface(mesh, 0, i+1, ((i+1)%tess_th)+1);
	for (int j = 1; j < tess_r; j++) {
		int base = 1 + (j-1) * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base+tess_th+i, base+tess_th+i1,
				base+i, base+i1);
		}
	}

	return mesh;
}


// Open cylinder of height 1 and given radius
TriMesh *make_cyl(int tess_th, int tess_h, float r)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_h < 1)
		tess_h = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(tess_th * (tess_h+1));
	for (int j = 0; j <= tess_h; j++) {
		float z = -1.0f + 2.0f * j / tess_h;
		for (int i = 0; i < tess_th; i++) {
			float th = M_2PIf * i / tess_th;
			mkpoint(mesh, r*cos(th), r*sin(th), z);
		}
	}

	mesh->faces.reserve(2*tess_th*tess_h);
	for (int j = 0; j < tess_h; j++) {
		int base = j * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base + i, base + i1,
				base+tess_th+i, base+tess_th+i1);
		}
	}

	return mesh;
}


// Cylinder capped with discs on both ends
TriMesh *make_ccyl(int tess_th, int tess_h, float r)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_h < 1)
		tess_h = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(2+3*tess_th*tess_h-tess_th);

	mkpoint(mesh, 0, 0, -1);
	for (int j = 1; j <= tess_h; j++) {
		float rr = r * j / tess_h;
		for (int i = 0; i < tess_th; i++) {
			float th = M_2PIf * i / tess_th;
			mkpoint(mesh, rr*cos(th), rr*sin(th), -1);
		}
	}
	int side_start = mesh->vertices.size();
	for (int j = 1; j < tess_h; j++) {
		float z = -1.0f + 2.0f * j / tess_h;
		for (int i = 0; i < tess_th; i++) {
			float th = M_2PIf * i / tess_th;
			mkpoint(mesh, r*cos(th), r*sin(th), z);
		}
	}
	int top_start = mesh->vertices.size();
	for (int j = tess_h; j > 0; j--) {
		float rr = r * j / tess_h;
		for (int i = 0; i < tess_th; i++) {
			float th = M_2PIf * i / tess_th;
			mkpoint(mesh, rr*cos(th), rr*sin(th), 1);
		}
	}
	mkpoint(mesh, 0, 0, 1);

	mesh->faces.reserve(6*tess_th*tess_h - 2*tess_th);

	for (int i = 0; i < tess_th; i++)
		mkface(mesh, 0, ((i+1)%tess_th)+1, i+1);
	for (int j = 1; j < tess_h; j++) {
		int base = 1 + (j-1) * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base+tess_th+i1, base+tess_th+i,
				base+i1, base+i);
		}
	}

	for (int j = 0; j < tess_h; j++) {
		int base = side_start - tess_th + j * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base + i, base + i1,
				base+tess_th+i, base+tess_th+i1);
		}
	}

	for (int j = 0; j < tess_h-1; j++) {
		int base = top_start + j * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base+tess_th+i1, base+tess_th+i,
				base+i1, base+i);
		}
	}
	int base = top_start + (tess_h-1)*tess_th;
	for (int i = 0; i < tess_th; i++)
		mkface(mesh, base+i, base+((i+1)%tess_th), base+tess_th);

	return mesh;
}


// Cylinder capped with hemispheres on both ends
TriMesh *make_scyl(int tess_th, int tess_h, float r)
{
	TriMesh *mesh = make_ccyl(tess_th, tess_h, r);

	int k = 0;
	mesh->vertices[k++][2] -= r;

	for (int j = 1; j < tess_h; j++) {
		float rr = r * sin(M_PI_2f * j / tess_h);
		float z = -1.0f - r * cos(M_PI_2f * j / tess_h);
		for (int i = 0; i < tess_th; i++) {
			float th = M_2PIf * i / tess_th;
			mesh->vertices[k++] = point(rr*cos(th), rr*sin(th), z);
		}
	}

	k += (tess_h + 1) * tess_th;

	for (int j = tess_h - 1; j > 0; j--) {
		float rr = r * sin(M_PI_2f * j / tess_h);
		float z = 1.0f + r * cos(M_PI_2f * j / tess_h);
		for (int i = 0; i < tess_th; i++) {
			float th = M_2PIf * i / tess_th;
			mesh->vertices[k++] = point(rr*cos(th), rr*sin(th), z);
		}
	}

	mesh->vertices[k][2] += r;

	return mesh;
}


// Open cone
TriMesh *make_cone(int tess_th, int tess_r, float r)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_r < 1)
		tess_r = 1;

	TriMesh *mesh = make_disc(tess_th, tess_r);
	for (int j = 0; j <= tess_r; j++) {
		float z = 1.0f - 2.0f * j / tess_r;
		if (j == 0) {
			mesh->vertices[0][2] = z;
			continue;
		}
		for (int i = 0; i < tess_th; i++) {
			int ind = 1 + i + (j-1)*tess_th;
			mesh->vertices[ind][0] *= r;
			mesh->vertices[ind][1] *= r;
			mesh->vertices[ind][2] = z;
		}
	}
	return mesh;
}


// Capped cone
TriMesh *make_ccone(int tess_th, int tess_r, float r)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_r < 1)
		tess_r = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(2+2*tess_th*tess_r-tess_th);

	mkpoint(mesh, 0, 0, -1);
	for (int j = 1; j <= tess_r; j++) {
		float rr = r * j / tess_r;
		for (int i = 0; i < tess_th; i++) {
			float th = M_2PIf * i / tess_th;
			mkpoint(mesh, rr*cos(th), rr*sin(th), -1);
		}
	}
	int side_start = mesh->vertices.size();
	for (int j = 1; j < tess_r; j++) {
		float z = -1.0f + 2.0f * j / tess_r;
		float rr = r * (tess_r - j) / tess_r;
		for (int i = 0; i < tess_th; i++) {
			float th = M_2PIf * i / tess_th;
			mkpoint(mesh, rr*cos(th), rr*sin(th), z);
		}
	}
	mkpoint(mesh, 0, 0, 1);

	mesh->faces.reserve(4*tess_th*tess_r - 2*tess_th);

	for (int i = 0; i < tess_th; i++)
		mkface(mesh, 0, ((i+1)%tess_th)+1, i+1);
	for (int j = 1; j < tess_r; j++) {
		int base = 1 + (j-1) * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base+tess_th+i1, base+tess_th+i,
				base+i1, base+i);
		}
	}

	for (int j = 0; j < tess_r-1; j++) {
		int base = side_start - tess_th + j * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base + i, base + i1,
				base+tess_th+i, base+tess_th+i1);
		}
	}

	int base = side_start + (tess_r-2)*tess_th;
	for (int i = 0; i < tess_th; i++)
		mkface(mesh, base+i, base+((i+1)%tess_th), base+tess_th);

	return mesh;
}


// Torus with major radius 1, given minor radius
TriMesh *make_torus(int tess_th, int tess_ph, float r)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_ph < 3)
		tess_ph = 3;

	TriMesh *mesh = make_cyl(tess_ph, tess_th);

	for (int i = 0; i < tess_ph; i++)
		mesh->vertices.pop_back();
	for (size_t i = 0; i < mesh->faces.size(); i++) {
		mesh->faces[i][0] %= mesh->vertices.size();
		mesh->faces[i][1] %= mesh->vertices.size();
		mesh->faces[i][2] %= mesh->vertices.size();
	}

	for (int j = 0; j < tess_th; j++) {
		float th = M_2PIf * j / tess_th;
		vec circlepos(cos(th), sin(th), 0);
		for (int i = 0; i < tess_ph; i++) {
			float ph = M_2PIf * i / tess_ph;
			mesh->vertices[i+j*tess_ph] = circlepos +
			                              cos(ph)*r*circlepos +
			                              sin(ph)*r*vec(0,0,-1);
		}
	}
	return mesh;
}


// Trefoil knot of the given minor radius
TriMesh *make_knot(int tess_th, int tess_ph, float r)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_ph < 3)
		tess_ph = 3;

	TriMesh *mesh = make_torus(tess_th, tess_ph, r);

	for (int j = 0; j < tess_th; j++) {
		float th = M_2PIf * j / tess_th;
		vec pos(2.0f*sin(2.0f*th) + cos(th),
			2.0f*cos(2.0f*th) + sin(th),
			cos(3.0f*th));
		pos /= 3.0f;
		vec vel( 4.0f*cos(2.0f*th) - sin(th),
			-4.0f*sin(2.0f*th) + cos(th),
			-3.0f*sin(3.0f*th));
		normalize(vel);
		for (int i = 0; i < tess_ph; i++) {
			float ph = M_2PIf * i / tess_ph;
			vec u = vel CROSS vec(0,0,1);
			normalize(u);
			vec v = u CROSS vel;
			mesh->vertices[i+j*tess_ph] =
				pos + cos(ph)*r*v + sin(ph)*r*u;
		}
	}
	return mesh;
}


// Klein bottle
TriMesh *make_klein(int tess_th, int tess_ph)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_ph < 3)
		tess_ph = 3;

	TriMesh *mesh = make_torus(tess_th, tess_ph, 0.2f);

	for (int j = 0; j < tess_th; j++) {
		float v = 8.0f * M_PIf * j / tess_th;
		int part = 0;
		for (int i = 0; i < 7; i++)
			if (v > M_PIf)
				v -= M_PIf, part++;
		float flipx = 1.0f;
		if (part > 3)
			part -= 4, flipx = -1.0f;
		if (part == 0 || part == 2)
			v = M_PIf - v;
		for (int i = 0; i < tess_ph; i++) {
			float u = M_2PIf * i / tess_ph;
			// Based on formulae from
			// http://emsh.calarts.edu/~mathart/sw/klein/Klein.html
			if (part == 0)
				mesh->vertices[i+j*tess_ph].set(
					(2.5f + 1.5f * cos(v)) * flipx * cos(u),
					(2.5f + 1.5f * cos(v)) * sin(u),
					-2.5f * sin(v));
			else if (part == 1)
				mesh->vertices[i+j*tess_ph].set(
					(2.5f + 1.5f * cos(v)) * flipx * cos(u),
					(2.5f + 1.5f * cos(v)) * sin(u),
					3.0f * v);
			else if (part == 2)
				mesh->vertices[i+j*tess_ph].set(
					2.0f + (2.0f - flipx * cos(u)) * cos(v),
					sin(u),
					3.0f * M_PIf + (2.0f - flipx * cos(u)) * sin(v));
			else
				mesh->vertices[i+j*tess_ph].set(
					2.0f + 2.0f * cos(v) - flipx * cos(u),
					sin(u),
					3.0f * (M_PIf - v));
		}
	}
	return mesh;
}


// Helix of major radius 1, given minor radius, and number of turns
TriMesh *make_helix(int tess_th, int tess_ph, float turns, float r)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_ph < 3)
		tess_ph = 3;

	TriMesh *mesh = make_cyl(tess_ph, tess_th);

	const float angle = atan(1.0f / M_PIf);
	for (int j = 0; j <= tess_th; j++) {
		float z = (-1.0f + 2.0f * j / tess_th) * turns;
		float th = M_2PIf * j / tess_th * turns;
		vec helixpos(cos(th), sin(th), z);
		vec xdir(cos(th), sin(th), 0);
		vec ydir = xform::rot(angle, xdir) * vec(0, 0, -1);
		for (int i = 0; i < tess_ph; i++) {
			float ph = M_2PIf * i / tess_ph;
			mesh->vertices[i+j*tess_ph] = helixpos +
			                              r*cos(ph)*xdir +
			                              r*sin(ph)*ydir;
		}
	}
	return mesh;
}


// Lat/long tessellated sphere
TriMesh *make_sphere_polar(int tess_ph, int tess_th)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_ph < 3)
		tess_ph = 3;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(2+tess_ph*(tess_th-1));

	mkpoint(mesh, 0, 0, -1);
	for (int j = 1; j < tess_th; j++) {
		float th = M_PIf * j / tess_th;
		float z = -cos(th);
		float r = sin(th);
		for (int i = 0; i < tess_ph; i++) {
			float ph = M_2PIf * i / tess_ph;
			mkpoint(mesh, r*cos(ph), r*sin(ph), z);
		}
	}
	mkpoint(mesh, 0, 0, 1);

	mesh->faces.reserve(2*tess_th*tess_ph - 2*tess_ph);

	for (int i = 0; i < tess_ph; i++)
		mkface(mesh, 0, ((i+1)%tess_ph)+1, i+1);

	for (int j = 0; j < tess_th-2; j++) {
		int base = 1 + j * tess_ph;
		for (int i = 0; i < tess_ph; i++) {
			int i1 = (i+1)%tess_ph;
			mkquad(mesh, base + i, base + i1,
				base+tess_ph+i, base+tess_ph+i1);
		}
	}

	int base = 1 + (tess_th-2)*tess_ph;
	for (int i = 0; i < tess_ph; i++)
		mkface(mesh, base+i, base+((i+1)%tess_ph), base+tess_ph);

	return mesh;
}


// Sphere subdivided nsubdiv times from a Platonic solid of nfaces faces
TriMesh *make_sphere_subdiv(int nfaces, int nsubdiv)
{
	TriMesh *mesh = make_platonic(nfaces);
	if (!mesh)
		return NULL;
	for (size_t i = 0; i < mesh->vertices.size(); i++)
		normalize(mesh->vertices[i]);
	for (int i = 0; i < nsubdiv; i++)
		subdiv(mesh);
	for (size_t i = 0; i < mesh->vertices.size(); i++)
		normalize(mesh->vertices[i]);
	return mesh;
}


// A pre-baked polyhedron (Platonic solid, etc.)
TriMesh *make_fixed_shape(FixedShape shape)
{
	const float *vert_data;
	const int *face_data;
	int nverts, face_data_size;

#define ASSIGN(shape_name, array_basename) \
	case shape_name: \
		vert_data = array_basename ## _verts; \
		face_data = array_basename ## _faces; \
		nverts = sizeof(array_basename ## _verts) / (3 * sizeof(float)); \
		face_data_size = sizeof(array_basename ## _faces) / sizeof(int); \
		break

	switch (shape) {
		ASSIGN(SHAPE_TETRAHEDRON, tetrahedron);
		ASSIGN(SHAPE_CUBE, cube);
		ASSIGN(SHAPE_OCTAHEDRON, octahedron);
		ASSIGN(SHAPE_DODECAHEDRON, dodecahedron);
		ASSIGN(SHAPE_ICOSAHEDRON, icosahedron);
		ASSIGN(SHAPE_TRUNCATED_TETRAHEDRON, truncated_tetrahedron);
		ASSIGN(SHAPE_CUBOCTAHEDRON, cuboctahedron);
		ASSIGN(SHAPE_TRUNCATED_CUBE, truncated_cube);
		ASSIGN(SHAPE_TRUNCATED_OCTAHEDRON, truncated_octahedron);
		ASSIGN(SHAPE_RHOMBICUBOCTAHEDRON, rhombicuboctahedron);
		ASSIGN(SHAPE_TRUNCATED_CUBOCTAHEDRON, truncated_cuboctahedron);
		ASSIGN(SHAPE_ICOSIDODECAHEDRON, icosidodecahedron);
		ASSIGN(SHAPE_TRUNCATED_DODECAHEDRON, truncated_dodecahedron);
		ASSIGN(SHAPE_TRUNCATED_ICOSAHEDRON, truncated_icosahedron);
		ASSIGN(SHAPE_SNUB_CUBE, snub_cube);
		ASSIGN(SHAPE_RHOMBICOSIDODECAHEDRON, rhombicosidodecahedron);
		ASSIGN(SHAPE_TRUNCATED_ICOSIDODECAHEDRON, truncated_icosidodecahedron);
		ASSIGN(SHAPE_SNUB_DODECAHEDRON, snub_dodecahedron);
		ASSIGN(SHAPE_TRIAKIS_TETRAHEDRON, triakis_tetrahedron);
		ASSIGN(SHAPE_RHOMBIC_DODECAHEDRON, rhombic_dodecahedron);
		ASSIGN(SHAPE_TRIAKIS_OCTAHEDRON, triakis_octahedron);
		ASSIGN(SHAPE_TETRAKIS_HEXAHEDRON, tetrakis_hexahedron);
		ASSIGN(SHAPE_DELTOIDAL_ICOSITETRAHEDRON, deltoidal_icositetrahedron);
		ASSIGN(SHAPE_DISDYAKIS_DODECAHEDRON, disdyakis_dodecahedron);
		ASSIGN(SHAPE_RHOMBIC_TRIACONTAHEDRON, rhombic_triacontahedron);
		ASSIGN(SHAPE_TRIAKIS_ICOSAHEDRON, triakis_icosahedron);
		ASSIGN(SHAPE_PENTAKIS_DODECAHEDRON, pentakis_dodecahedron);
		ASSIGN(SHAPE_PENTAGONAL_ICOSITETRAHEDRON, pentagonal_icositetrahedron);
		ASSIGN(SHAPE_DELTOIDAL_HEXECONTAHEDRON, deltoidal_hexecontahedron);
		ASSIGN(SHAPE_DISDYAKIS_TRIACONTAHEDRON, disdyakis_triacontahedron);
		ASSIGN(SHAPE_PENTAGONAL_HEXECONTAHEDRON, pentagonal_hexecontahedron);

		default:
			TriMesh::eprintf("Unknown shape in make_fixed_shape()\n");
			return NULL;
	}
#undef ASSIGN

	// Make the mesh and assign vertices
	TriMesh *mesh = new TriMesh;
	mesh->vertices.resize(nverts);
	memcpy(&mesh->vertices[0][0], &vert_data[0], nverts * 3 * sizeof(float));

	// Scale to fit inside unit sphere
	float maxlen2 = len2(mesh->vertices[0]);
	for (int i = 1; i < nverts; i++) {
		float l2 = len2(mesh->vertices[i]);
		if (l2 > maxlen2)
			maxlen2 = l2;
	}
	float scale = 1.0f / sqrt(maxlen2);
	for (int i = 0; i < nverts; i++)
		mesh->vertices[i] *= scale;

	// Assign faces
	int i = 0;
	while (i < face_data_size) {
		int face_size = face_data[i++];
		if (face_size == 3) {
			mkface(mesh, face_data[i], face_data[i+1],
				face_data[i+2]);
			i += 3;
			continue;
		}
		point center;
		for (int j = 0; j < face_size; j++)
			center += mesh->vertices[face_data[i + j]];
		center *= 1.0f / face_size;
		int center_ind = mesh->vertices.size();
		mesh->vertices.push_back(center);
		for (int j = 0; j < face_size; j++)
			mkface(mesh, center_ind, face_data[i + j],
			       face_data[i + (j+1) % face_size]);
		i += face_size;
	}

	return mesh;
}


// Backwards compatibility - make a Platonic solid with the given
// number of faces
TriMesh *make_platonic(int nfaces)
{
	switch (nfaces) {
	case 4:
		return make_fixed_shape(SHAPE_TETRAHEDRON);
	case 6:
		return make_fixed_shape(SHAPE_CUBE);
	case 8:
		return make_fixed_shape(SHAPE_OCTAHEDRON);
	case 12:
		return make_fixed_shape(SHAPE_DODECAHEDRON);
	case 20:
		return make_fixed_shape(SHAPE_ICOSAHEDRON);
	default:
		TriMesh::eprintf("make_platonic: Number of faces must be 4, 6, 8, 12, or 20\n");
		return NULL;
	}
}


// Surface of revolution obtained by rotating points on the given curve
// around the z axis.
TriMesh *make_surface_of_revolution(int tess_th, const vector<point> &curve_pts)
{
	int n = curve_pts.size();
	if (n < 2) {
		TriMesh::eprintf("Not enough points in make_sor\n");
		return NULL;
	}

	TriMesh *mesh = make_cyl(tess_th, n - 1);

	for (int i = 0; i < tess_th; i++) {
		xform xf = xform::rot(2.0 * M_PI * i / tess_th, 0, 0, 1);
		for (int j = 0; j < n; j++)
			mesh->vertices[i + j * tess_th] = xf * curve_pts[j];
	}

	return mesh;
}


// Make a teapot.  Purists should pass in true for the last parameters.
TriMesh *make_teapot(int tess, bool omit_bottom, bool taller)
{
	// See Frank Crow's 1987 paper
	const float taller_scale = 3.0f / teapot_verts[2];

	if (tess < 1)
		tess = 1;

	vector<vec4> bezier_vals(tess + 1);
	for (int i = 0; i <= tess; i++) {
		float t = (float) i / tess;
		float s = 1.0f - t;
		bezier_vals[i][0] = s * s * s;
		bezier_vals[i][1] = 3.0f * s * s * t;
		bezier_vals[i][2] = 3.0f * s * t * t;
		bezier_vals[i][3] = t * t * t;
	}

	const int patch_size = 4 * 4;
	int npatches = sizeof(teapot_patches) / (patch_size * sizeof(int));
	if (omit_bottom)
		npatches -= 4;

	TriMesh *mesh = new TriMesh;
	for (int patch_ind = 0; patch_ind < npatches; patch_ind++) {
		const int *patch = teapot_patches + patch_ind * patch_size;
		vector<int> patch_grid(sqr(tess+1));
		for (int x = 0; x <= tess; x++) {
		 for (int y = 0; y <= tess; y++) {
			const vec4 &wx = bezier_vals[x];
			const vec4 &wy = bezier_vals[y];
			point p;
			for (int i = 0; i < 16; i++) {
				p += point(teapot_verts + 3 * patch[i])
					* wx[i % 4] * wy[i / 4];
			}
			if (taller)
				p[2] *= taller_scale;
			int nv = mesh->vertices.size();
			patch_grid[x + y * (tess+1)] = nv;
			if ((x == 0 || x == tess || y == 0 || y == tess) &&
			    !(patch_ind == 14 && x == 0 && y == tess)) {
				for (int i = nv - 1; i >= 0; i--) {
					if (dist2(mesh->vertices[i], p) < 1.0e-12f) {
						patch_grid[x + y * (tess+1)] = i;
						break;
					}
				}
			}
			if (patch_grid[x + y * (tess+1)] == nv)
				mesh->vertices.push_back(p);
		 }
		}
		for (int x = 0; x < tess; x++) {
		 for (int y = 0; y < tess; y++) {
			int ll = patch_grid[x   +  y    * (tess+1)];
			int lr = patch_grid[x+1 +  y    * (tess+1)];
			int ul = patch_grid[x   + (y+1) * (tess+1)];
			int ur = patch_grid[x+1 + (y+1) * (tess+1)];
			mkquad(mesh, ll, lr, ul, ur);
		 }
		}
	}
	return mesh;
}

} // namespace trimesh
