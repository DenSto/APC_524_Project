# for serial version:
#USE_MPI = 0
# for parallel version:
#USE_MPI = 1
#Note: be sure to make clean after changing USE_MPI!

# run test:
#RUN_TEST = 1
# only compile:
RUN_TEST = 0

#ifeq ($(USE_MPI),0)
CC = g++
 # OMPFLAGS = -DUSE_MPI=0
#else
MPICC = mpic++
  #OMPFLAGS = -fopenmp -DUSE_MPI=1
OMPFLAGS = -DUSE_MPI=1
#endif

CFLAGS= -g  -Wall 

include Makefile.obj

SRCDIR = $(PWD)

GTEST_DIR = $(SRCDIR)/../googletest
GTEST_INCL = $(GTEST_DIR)/include

LCONFIG_DIR = $(SRCDIR)/../libconfig-1.5
LCONFIG_INCL = $(LCONFIG_DIR)/include
LCONFIG_LIB = $(LCONFIG_DIR)/lib/libconfig++.a

HDF5_DIR = $(SRCDIR)/../hdf5-1.8.12
HDF5PAR_DIR = $(SRCDIR)/../hdf5-1.8.12-parallel

HDF5_INCL = $(HDF5_DIR)/include
HDF5PAR_INCL = $(HDF5PAR_DIR)/include

HDF5_LIB = $(HDF5_DIR)/lib/libhdf5.a
HDF5PAR_LIB = $(HDF5PAR_DIR)/lib/libhdf5.a

LIBS = $(LCONFIG_LIB) $(HDF5_LIB)
LIBS_PAR = $(LCONFIG_LIB) $(HDF5PAR_LIB)

#ifeq ($(USE_MPI),0)
#  LIBS += $(HDF5_LIB)
#else
#  HDF5_INCL = $(HDF5PAR_DIR)/include
#  LIBS += $(HDF5PAR_LIB)
#endif

PROGRAM = EMOOPIC
PROGRAM_PAR = $(PROGRAM)_par

ifeq ($(RUN_TEST),1)
  all: $(PROGRAM) $(PROGRAM_PAR) run_tests
else
  all: $(PROGRAM) $(PROGRAM_PAR)
endif

$(PROGRAM): $(OBJ) $(LIBS)
	$(CC) -o $@ $^ $(CFLAGS) -lm -lz -ldl 

$(PROGRAM_PAR): $(OBJ_PAR) $(LIBS_PAR)
	$(MPICC) -o $@ $^ $(CFLAGS) -lm -lz -ldl $(OMPFLAGS)

%.o: %.cpp $(LIBS)
	$(CC) -c -o $@ $< $(CFLAGS) -I $(GTEST_INCL) -I $(HDF5_INCL)

%_par.o : %.cpp $(LIBS_PAR)
	$(MPICC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS) -I $(GTEST_INCL) -I $(HDF5PAR_INCL)

$(LCONFIG_LIB):
	cd ..; \
	tar xf libconfig-1.5.tar.gz; \
	cd libconfig-1.5; \
	./configure --prefix=$(LCONFIG_DIR); \
	make; \
	make install; \
	cd $(SRCDIR)

$(HDF5_LIB):
	cd ..; \
	[ -d $(HDF5_DIR) ] || tar xf hdf5-1.8.12.tar.gz; \
	cd $(HDF5_DIR); \
	CC=cc ./configure --prefix=$(HDF5_DIR); \
	make; \
	make install; \
	cd $(SRCDIR)

$(HDF5PAR_LIB):
	cd ..; \
	[ -d $(HDF5_DIR) ] || tar xf hdf5-1.8.12.tar.gz; \
	cd $(HDF5_DIR); \
	CC=mpicc ./configure --enable-parallel --prefix=$(HDF5PAR_DIR); \
	make; \
	make install; \
	cd $(SRCDIR)

# adroit doesn't have makedepend...
#DEPFILE	= .depend
#DEPTOKEN = '\# MAKEDEPENDS'
#DEPFLAGS = -Y -I $(HDF5_DIR) -I $(GTEST_DIR) -f $(DEPFILE) -s $(DEPTOKEN) -p $(OUTDIR)/
#
#
#depend:
#	rm -f .depend
#	make $(DEPFILE)
#
#$(DEPFILE):
#	@echo $(DEPTOKEN) > $(DEPFILE)
#	makedepend $(DEPFLAGS) -- $(CFLAGS) -- */*.cpp

clean:
	rm -f $(PROGRAM) $(PROGRAM_PAR) $(OBJ) $(OBJ_PAR) .depend

run_tests:
ifeq ($(USE_MPI),0)
	cd ../test; \
	make; \
	make run_tests; \
	cd ../src;
endif

sinclude .depend

