# DEBUG = y
MAKERULESDIR = ..
DESTDIR = ../bin.$(UNAME)
INCLUDES = -I../include
LIBDIR = -L../lib.$(UNAME)

include $(MAKERULESDIR)/Makerules

VIEWSOURCES =	mesh_view.cc

OTHERSOURCES =	mesh_align.cc \
		mesh_cat.cc \
		mesh_cc.cc \
		mesh_check.cc \
		mesh_crunch.cc \
		mesh_filter.cc \
		mesh_hf.cc \
		mesh_info.cc \
		mesh_make.cc \
		mesh_shade.cc \
		xf.cc

VIEWOFILES = $(addprefix $(OBJDIR)/,$(VIEWSOURCES:.cc=.o))
OTHEROFILES = $(addprefix $(OBJDIR)/,$(OTHERSOURCES:.cc=.o))
OFILES = $(VIEWOFILES) $(OTHEROFILES)

VIEWPROG = $(DESTDIR)/mesh_view$(EXE)
OTHERPROGS = $(addsuffix $(EXE), $(addprefix $(DESTDIR)/, $(OTHERSOURCES:.cc=)))
PROGS = $(OTHERPROGS) $(VIEWPROG)

default: $(PROGS)


LIBS += -ltrimesh
$(VIEWPROG) : $(VIEWOFILES)
	$(LINK) $(GLLIBS)

$(OTHERPROGS) : $(DESTDIR)/%$(EXE) : $(OBJDIR)/%.o
	$(LINK)

$(PROGS) : ../lib.$(UNAME)/libtrimesh.a

clean :
	-rm -f $(OFILES) $(OBJDIR)/Makedepend $(OBJDIR)/*.d
	-rm -rf $(OBJDIR)/ii_files
	-rmdir $(OBJDIR)

spotless : clean
	-rm -f $(PROGS)
	-rmdir $(DESTDIR)

