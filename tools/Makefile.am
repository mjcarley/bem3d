## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)
LDADD = $(top_builddir)/src/.libs/libbem3d.a 

#AM_CFLAGS = @wmpi_cflags@
#AM_LDFLAGS = @wmpi_libs@

DEPS = $(top_builddir)/src/libbem3d.la

EXTRA_DIST = tools.h

bin_PROGRAMS = \
	bem3d-assemble \
	bem3d-solve \
	bem3d-dump \
	gts2bem3d \
	msh2bem3d \
	bem3d2pos \
	bem3d2msh \
	bem3d-plane \
	bem3d-sphere \
	bem3d-field \
	bem3d-aero \
	bem3d-function \
	bem3d-skeleton \
	bem3d-configuration

#bem3d_assemble_CFLAGS = @gqr_cflags@ @wmpi_cflags@
#bem3d_assemble_LDFLAGS = @gqr_libs@ @wmpi_libs@
bem3d_solve_CFLAGS = @sisl_cflags@ @wmpi_cflags@
bem3d_solve_LDFLAGS = @sisl_libs@ @wmpi_libs@
bem3d_aero_CFLAGS = @gqr_cflags@ @sisl_cflags@ @wmpi_cflags@
bem3d_aero_LDFLAGS = @gqr_libs@ @sisl_libs@ @wmpi_libs@
