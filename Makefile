CFLAGS = -I -Isrc -Iinclude
DEPS = src/*.c
DEPS_O = *.o
OBJ = main
ADDITIONAL_FLAGS = -g
RUNFLAGS =
OBJDIR = bin
CC = gcc

make:
	$(CC) $(CFLAGS) $(ADDITIONAL_FLAGS) -c $(DEPS)
	$(CC) $(CFLAGS) $(ADDITIONAL_FLAGS) main.c $(DEPS_O) -o $(OBJ)
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