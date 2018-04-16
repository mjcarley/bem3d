SHELL = /bin/sh
CC = gcc -g -O3
CFLAGS = `pkg-config --cflags glib-2.0`
LIBS = `pkg-config --libs glib-2.0` -lgfortran -lm

FC = f95 -g -O3

FMMSOURCE = ..
FMMOBJ = \
$(FMMSOURCE)/src/cdjseval3d.o \
$(FMMSOURCE)/src/d3mtreeplot.o \
$(FMMSOURCE)/src/d3tstrcr.o \
$(FMMSOURCE)/src/dfft.o \
$(FMMSOURCE)/src/h3dterms.o \
$(FMMSOURCE)/src/h3dtrans.o \
$(FMMSOURCE)/src/h3dtrirouts.o \
$(FMMSOURCE)/src/helmrouts3d.o \
$(FMMSOURCE)/src/hfmm3dpart.o \
$(FMMSOURCE)/src/hfmm3drouts.o \
$(FMMSOURCE)/src/hfmm3dtria.o \
$(FMMSOURCE)/src/l3dterms.o \
$(FMMSOURCE)/src/l3dtrans.o \
$(FMMSOURCE)/src/l3dtrirouts.o \
$(FMMSOURCE)/src/laprouts3d.o \
$(FMMSOURCE)/src/legeexps.o \
$(FMMSOURCE)/src/lfmm3dpart.o \
$(FMMSOURCE)/src/lfmm3drouts.o \
$(FMMSOURCE)/src/lfmm3dtria.o \
$(FMMSOURCE)/src/prini.o \
$(FMMSOURCE)/src/prinm.o \
$(FMMSOURCE)/src/projections.o \
$(FMMSOURCE)/src/rotproj.o \
$(FMMSOURCE)/src/rotviarecur3.o \
$(FMMSOURCE)/src/triagauc.o \
$(FMMSOURCE)/src/triahquad.o \
$(FMMSOURCE)/src/triasymq.o \
$(FMMSOURCE)/src/trilib.o \
$(FMMSOURCE)/src/triquadflatlib.o \
$(FMMSOURCE)/src/yrecursion.o

OBJ = fmm-wrapper.o fmm.o

all: fmm-wrapper

fmm-wrapper: $(OBJ) $(FMMOBJ)
	$(CC) $(CFLAGS) $(OBJ) $(FMMOBJ) $(LIBS) -o fmm-wrapper

clean:
	rm -f *.c~ *.h~ *.o $(FMMSOURCE)/src/*.o
