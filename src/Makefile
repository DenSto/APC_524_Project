USE_MPI = 1 

ifeq ($(USE_MPI),0)
  CC = icpc
  OMPFLAGS = 
else
  CC = mpic++
  #OMPFLAGS = -fopenmp -DUSE_MPI=1
  OMPFLAGS = -DUSE_MPI=1
endif

CFLAGS= -g -lm -Wall  
PROGRAM = EMOOPIC

IOOBJ = IO.o
PUSHOBJ = boris.o 
GRIDOBJ = grid.o
PARTOBJ = particle.o

OBJ = ./IO/$(IOOBJ) \
      ./pusher/$(PUSHOBJ) \
      ./grid/$(GRIDOBJ) \
      ./particles/$(PARTOBJ) \
      main.o 

all: $(PROGRAM)

$(PROGRAM): $(OBJ) 
	$(CC) -o $@ $^ $(CFLAGS) $(OMPFLAGS)

%.o: %.cpp 
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS)

clean:
	rm -f $(PROGRAM) $(OBJ) 
