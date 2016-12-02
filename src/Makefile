USE_MPI = 0

ifeq ($(USE_MPI),1)
  CC = mpic++
  OMPFLAGS = -openmp -DUSE_MPI=1
else
  CC = icpc
  OMPFLAGS = 
endif

CFLAGS= -lm -Wall  
PROGRAM = EMOOPIC

IOOBJ = ./IO/IO.o
GRIDOBJ =
PARTOBJ =

OBJ = $(IOOBJ) $(GRIDOBJ) $(PARTOBJ) main.o 

all: $(PROGRAM)

$(PROGRAM): $(OBJ) 
	$(CC) -o $@ $^ $(CFLAGS) $(OMPFLAGS)

%.o: %.cpp 
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS)

clean:
	rm -f $(PROGRAM) *.o ./IO/*.o 
