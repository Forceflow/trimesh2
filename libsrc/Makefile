# DEBUG = y
MAKERULESDIR = ..
DESTDIR = ../lib.$(UNAME)
INCLUDES = -I../include

include $(MAKERULESDIR)/Makerules

CCFILES =	TriMesh_bounding.cc \
		TriMesh_connectivity.cc \
		TriMesh_curvature.cc \
		TriMesh_io.cc \
		TriMesh_grid.cc \
		TriMesh_normals.cc \
		TriMesh_pointareas.cc \
		TriMesh_stats.cc \
		TriMesh_tstrips.cc \
		GLCamera.cc \
		ICP.cc \
		KDtree.cc \
		conn_comps.cc \
		diffuse.cc \
		edgeflip.cc \
		faceflip.cc \
		filter.cc \
		lmsmooth.cc \
		overlap.cc \
		remove.cc \
		reorder_verts.cc \
		shared.cc \
		subdiv.cc


OFILES = $(addprefix $(OBJDIR)/,$(CCFILES:.cc=.o))
LIB = $(DESTDIR)/libtrimesh.a
default: $(LIB)

$(LIB): $(OFILES)
	$(STATICLIB)

clean:
	-rm -f $(OFILES) $(OBJDIR)/Makedepend $(OBJDIR)/*.d
	-rm -rf $(OBJDIR)/ii_files
	-rmdir $(OBJDIR)

spotless: clean
	-rm -f $(LIB)
	-rmdir $(DESTDIR)

