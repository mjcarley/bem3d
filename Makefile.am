## Process this file with automake to produce Makefile.in

# AUTOMAKE_OPTIONS = readme-alpha
ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = \
	octave/bem*.m \
	examples/Sphere/*.geo \
	examples/Sphere/*.bem \
	examples/Sphere/*.fn \
	examples/Sphere/*.cfg \
	examples/Sphere/sphere \
	examples/Scattering/*.geo \
	examples/Scattering/*.bem \
	examples/Scattering/*.fn \
	examples/Scattering/*.cfg \
	examples/Scattering/sphere

DIST_SUBDIRS = src tools mpitools doc

SUBDIRS = src tools @mpi_tools@ doc
