CFLAGS = -I -Isrc -Iinclude -Iinclude/extern
DEPS = src/*.c
DEPS_O = *.o
OBJ = main
ADDITIONAL_FLAGS = -g -fopenmp
RUNFLAGS =
OBJDIR = bin
CC = gcc
OMP_NUM_THREADS=4

make:
	$(CC) $(CFLAGS) $(ADDITIONAL_FLAGS) -c $(DEPS)
	$(CC) $(CFLAGS) $(ADDITIONAL_FLAGS) main.c $(DEPS_O) -o $(OBJ) -lm
	$(RUNFLAGS) ./$(OBJ)
.PHONY: clean
clean:
	rm -rf *.o *.vtk
.PHONY: debug
debug:
	$(CC) $(CFLAGS) $(ADDITIONAL_FLAGS) -c -g $(DEPS)

multi:
	mpic++ $(CFLAGS) -c -g $(DEPS)


.PHONY: run
run:
	$(CC) $(CFLAGS)  -g main.cpp $(DEPS_O) -o $(OBJ)
	$(RUNFLAGS) ./$(OBJ)
.PHONY: graphing
graphing:
	gcc -c include/extern/pbPlots.c -std=c99 -O3 -march=native
	gcc -c include/extern/supportLib.c -std=c99 -O3 -march=native