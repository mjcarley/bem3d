## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(includedir) -DG_LOG_DOMAIN=\"BEM3D\"

bin_SCRIPTS=bem3d-config

BUILT_SOURCES= \
	bem3d-config

bem3d-config: bem3d-config.in

lib_LTLIBRARIES = libbem3d.la

libbem3d_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

libbem3d_la_SOURCES = \
	acoord.c \
	belement.c \
	logging.c \
	files.c \
	bmesh.c \
	elements.c \
	functions.c \
	geometry.c \
	gfunc.c \
	gmsh.c \
	lookup.c \
	qselect.c \
	quadrature.c \
	hayami.c \
	radiation.c \
	shapefunc.c \
	bem3d-private.h \
	bem3d.h \
	bem3dconfig.h \
	wandzura.c \
	xiao.c \
	location.c \
	matrix.c \
	operators.c \
	newman.c \
	edge.c \
	motion.c \
	util.c \
	data.c \
	configure.c \
	polar.c \
	polar.h \
	series.c \
	intsincos.c \
	gradient.c \
	transforms.c \
	qpanel.c \
	cpanel.c \
	skeleton.c \
	fmm.c \
	fmmlib3d_1_2.c \
	fmmlib3d_1_2.h \
	wbfmm-bem3d.c \
	wbfmm-bem3d.h \
	htriquad.c \
	htriquad.h \
	binomials.c \
	binomials.h \
	expansions.c \
	htri-transforms.c \
	trace.c \
	trace.h \
	reduction.c \
	work.c

if USE_FMMLIB3D_1_2
libbem3d_la_SOURCES += \
	../fmmlib3d-1.2/src/cdjseval3d.f \
	../fmmlib3d-1.2/src/d3mtreeplot.f \
	../fmmlib3d-1.2/src/d3tstrcr.f \
	../fmmlib3d-1.2/src/dfft.f \
	../fmmlib3d-1.2/src/h3dterms.f \
	../fmmlib3d-1.2/src/h3dtrans.f \
	../fmmlib3d-1.2/src/h3dtrirouts.f \
	../fmmlib3d-1.2/src/helmrouts3d.f \
	../fmmlib3d-1.2/src/hfmm3dpart.f \
	../fmmlib3d-1.2/src/hfmm3drouts.f \
	../fmmlib3d-1.2/src/hfmm3dtria.f \
	../fmmlib3d-1.2/src/l3dterms.f \
	../fmmlib3d-1.2/src/l3dtrans.f \
	../fmmlib3d-1.2/src/l3dtrirouts.f \
	../fmmlib3d-1.2/src/laprouts3d.f \
	../fmmlib3d-1.2/src/legeexps.f \
	../fmmlib3d-1.2/src/lfmm3dpart.f \
	../fmmlib3d-1.2/src/lfmm3drouts.f \
	../fmmlib3d-1.2/src/lfmm3dtria.f \
	../fmmlib3d-1.2/src/prini.f \
	../fmmlib3d-1.2/src/prinm.f \
	../fmmlib3d-1.2/src/projections.f \
	../fmmlib3d-1.2/src/rotproj.f \
	../fmmlib3d-1.2/src/rotviarecur3.f \
	../fmmlib3d-1.2/src/triagauc.f \
	../fmmlib3d-1.2/src/triahquad.f \
	../fmmlib3d-1.2/src/triasymq.f \
	../fmmlib3d-1.2/src/trilib.f \
	../fmmlib3d-1.2/src/triquadflatlib.f \
	../fmmlib3d-1.2/src/yrecursion.f
endif

include_HEADERS = \
	bem3d.h \
	bem3dconfig.h

CLEANFILES = $(BUILT_SOURCES)
