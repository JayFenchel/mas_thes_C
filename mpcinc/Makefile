CC = gcc
FLAGS = -Os -Wall -Wstrict-prototypes -pedantic
OPT = -O3 -funroll-loops
STD = -std=c99

OBJ := $(patsubst %.c,%.o,$(wildcard *.c))
INC := ./include
LIB := ./lib

main: $(OBJ) libmpcformqpx
	$(CC) $(FLAGS) $(OPT) $(STD) -I$(INC) -L./ -o main $(OBJ) -lm #-lmpcformqpx 

libmpcformqpx: $(OBJ)
	$(CC) $(FLAGS) $(OPT) $(STD) -shared -Wl,-soname,libmpcformqpx.so.1 -o libmpcformqpx.so $(OBJ) -lm

$(OBJ): %.o: %.c
	$(CC) $(FLAGS) $(OPT) $(STD) -fPIC -I$(INC) -c $< -o $@

.PHONY: clean

clean:
	rm $(OBJ) #$(LIB)/*
