# for serial version:
USE_MPI = 0
# for parallel version:
#USE_MPI = 1
#Note: be sure to make clean after changing USE_MPI!

# run test:
#RUN_TEST = 1
# only compile:
RUN_TEST = 0


ifeq ($(USE_MPI),0)
  CC = g++
  OMPFLAGS = -DUSE_MPI=0
else
  CC = mpic++
  #OMPFLAGS = -fopenmp -DUSE_MPI=1
  OMPFLAGS = -DUSE_MPI=1
endif

CFLAGS= -g -lm -Wall 
PROGRAM = EMOOPIC

include Makefile.obj

SRCDIR = $(PWD)
LCONFIG_DIR = ../libconfig-1.5
LCONFIG = $(LCONFIG_DIR)/lib/libconfig++.a

LIBS = $(LCONFIG)

ifeq ($(RUN_TEST),1)
  all: $(PROGRAM) run_tests
else
  all: $(PROGRAM)
endif

$(PROGRAM): $(OBJ) $(LIBS)
	$(CC) -o $@ $^ $(CFLAGS) $(OMPFLAGS)

%.o: %.cpp 
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS) -I ../googletest/include

IO/input.o: IO/input.cpp $(LCONFIG)
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS) -I $(LCONFIG_DIR)/include

$(LCONFIG):
	cd ..; \
	tar xf libconfig-1.5.tar.gz; \
	cd libconfig-1.5; \
	./configure --prefix=$(PWD)/$(LCONFIG_DIR); \
	make; \
	make install; \
	cd $(SRCDIR)
depend:
	$(CC) -MM $(CFLAGS) *.cpp > .depend

clean:
	rm -f $(PROGRAM) $(OBJ) .depend


ifeq ($(RUN_TEST),1)
run_tests:
ifeq ($(USE_MPI),0)
	cd ../test; \
	make; \
	make run_tests; \
	cd ../src;
endif
endif

-include .depend

