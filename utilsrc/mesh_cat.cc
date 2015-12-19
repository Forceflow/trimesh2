/*
Szymon Rusinkiewicz
Princeton University

mesh_cat.cc
Concatenate meshes together
*/

#include "TriMesh.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace std;
using namespace trimesh;


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s infiles... -o outfile\n", myname);
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc < 4)
		usage(argv[0]);

	TriMesh *outmesh = new TriMesh;
	const char *outfile = NULL;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-o") == 0 && i < argc-1) {
			outfile = argv[i+1];
			i++;
			continue;
		}
		TriMesh *m = TriMesh::read(argv[i]);
		if (!m) {
			fprintf(stderr, "Couldn't read file %s\n", argv[i]);
			continue;
		}
		int onv = outmesh->vertices.size();
		outmesh->vertices.insert(outmesh->vertices.end(),
				         m->vertices.begin(),
					 m->vertices.end());

		if (outmesh->colors.empty() && !m->colors.empty())
			outmesh->colors.resize(onv, Color(1,1,1));
		else if (m->colors.empty() && !outmesh->colors.empty())
			m->colors.resize(m->vertices.size(), Color(1,1,1));
		outmesh->colors.insert(outmesh->colors.end(),
				       m->colors.begin(),
				       m->colors.end());

		if (outmesh->confidences.empty() && !m->confidences.empty())
			outmesh->confidences.resize(onv);
		else if (m->confidences.empty() && !outmesh->confidences.empty())
			m->confidences.resize(m->vertices.size());
		outmesh->confidences.insert(outmesh->confidences.end(),
					    m->confidences.begin(),
					    m->confidences.end());

		if (outmesh->normals.empty() && !m->normals.empty()) {
			outmesh->need_normals();
			outmesh->normals.resize(onv);
		} else if (m->normals.empty() && !outmesh->normals.empty())
			m->need_normals();
		outmesh->normals.insert(outmesh->normals.end(),
					m->normals.begin(),
					m->normals.end());

		m->need_faces();
		for (size_t i = 0; i < m->faces.size(); i++) {
			m->faces[i][0] += onv;
			m->faces[i][1] += onv;
			m->faces[i][2] += onv;
		}
		outmesh->faces.insert(outmesh->faces.end(),
				      m->faces.begin(),
				      m->faces.end());
		delete m;
	}
	if (outfile)
		outmesh->write(outfile);
	else
		fprintf(stderr, "No output file specified\n");
}
