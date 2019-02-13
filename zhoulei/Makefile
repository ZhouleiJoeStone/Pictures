include example.mk

CC=mpic++

LDIR =

#OBJ = test_WendlandQuinticC2Kernel.o
OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: main

run: all
	mpirun -np 4 ./main

.PHONY: clean all run

clean:
	rm -f *.o *~ core main Geo*

