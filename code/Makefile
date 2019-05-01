include ../../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

OPT=

all: main

main_test: OPT += -DTEST_RUN
main_test: all

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

main: $(OBJ)
	$(CC) -o $@ $^ $(OPT) $(CFLAGS) $(LIBS_PATH) $(LIBS)


run: main_test
	mpirun -np 2 ./main

.PHONY: clean all run

clean:
	rm -f *.o *~ core main *.vtk

cleanResult:
	rm -f /home/ubuntu/zhoulei/github/Pravite-Miscellany/codes/OpenFPM_results/DEM/benchmark/*.vtk

