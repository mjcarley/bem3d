## Process this file with automake to produce Makefile.in

CC = `wmpi-config --cc --mpi`
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)
LDADD = $(top_builddir)/src/.libs/libbem3d.a 
AM_CFLAGS = @sisl_cflags@ `wmpi-config --cflags --mpi`
DEPS = $(top_builddir)/src/libbem3d.la

bin_SCRIPTS=bem3d-pbs-run \
	bem3d-pbs-launch

bin_PROGRAMS = \
	bem3d-assemble.mpi \
	bem3d-solve.mpi

bem3d_assemble_mpi_SOURCES = ../tools/bem3d-assemble.c
bem3d_solve_mpi_SOURCES = ../tools/bem3d-solve.c

bem3d_assemble_mpi_CFLAGS = @sisl_cflags@ `wmpi-config --cflags --mpi`
bem3d_assemble_mpi_LDFLAGS = -static @sisl_libs@ `wmpi-config --libs --mpi`
bem3d_solve_mpi_CFLAGS = @sisl_cflags@ `wmpi-config --cflags --mpi`
bem3d_solve_mpi_LDFLAGS = -static @sisl_libs@ `wmpi-config --libs --mpi`
