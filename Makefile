# Compiler
CC = gcc
# Compiler-Optionen
FLAGS = -Os -Wall -Wstrict-prototypes -pedantic
OPT = -O3 -funroll-loops
STD = -std=c99

# Objekt-Dateien
# mit "$(VAR)" greift man auf den Inhalt einer Variablen VAR zu
# "wildcard *.c" erzeugt eine Variable mit einer Liste aller *.c files
# "patsubst %.c,%.o,TEXT" sucht in TEXT nach der pattern %.c und ersetzt sie durch %.o
# "patsubst %.c,%.o,$(wildcard *.c)" holt sich eine Liste aus *.c files und erzeugt daraus eine Liste der dazugehörigen Objekt-Dateien
OBJ := $(patsubst %.c,%.o,$(wildcard *.c))
INC := ./include
LIB := ./lib

# erstes target wird ausgeführt, wenn make ohne Argumente aufgerufen wird
main: $(OBJ) #libmpcformqpx
	$(CC) $(FLAGS) $(OPT) $(STD) -I$(INC) -L./ -o main $(OBJ) -lm #-lmpcformqpx  # "-I" specifies a directory dir to search for  included  makefiles; Nach "$(OBJ)" stehen zu linkende libraries

# erstellt alle *.o Objekt-Dateien zu den *.c Dateien
$(OBJ): %.o: %.c
	$(CC) $(FLAGS) $(OPT) $(STD) -fPIC -I$(INC) -c $< -o $@

# erzeugt eine Bibliothek
libmpcformqpx: $(OBJ)
	$(CC) $(FLAGS) $(OPT) $(STD) -shared -Wl,-soname,$(LIB)/libmpcformqpx.so.1 -o $(LIB)/libmpcformqpx.so $(OBJ) -lm

# löschen alter Objekt-Dateien
.PHONY: clean  # nicht auf Aktualität prüfen, immer ausführen
clean:
	rm -vfr main $(OBJ) $(LIB)/*

# Teil für Tests
test: math-test solve-test former-test
	./tests/math-test # Test für grundlegendere Mathefunktionen und update P
	./tests/solve-test # Test für Teile des Algorithmus
	./tests/former-test # Test für Teile des Algorithmus

former-test: test_former.o hhmpcipm.o hhmpcsolve.o hhmpcmath.o hhmpcusefull.o mpcincmtxops.o
	gcc -o ./tests/former-test hhmpcipm.o hhmpcsolve.o hhmpcmath.o hhmpcusefull.o mpcincmtxops.o ./tests/test_former.o -lcheck -lpthread -lrt -lm

test_former.o: test_former.c $(INC)/hhmpcsolve.h
	gcc $(FLAGS) -c -o ./tests/test_former.o ./tests/test_former.c

test_former.c: ./tests/test_former.test
	checkmk tests/test_former.test >tests/test_former.c


solve-test: test_solve.o hhmpcsolve.o hhmpcmath.o hhmpcusefull.o mpcincmtxops.o
	gcc -o ./tests/solve-test hhmpcsolve.o hhmpcmath.o hhmpcusefull.o mpcincmtxops.o ./tests/test_solve.o -lcheck -lpthread -lrt -lm

test_solve.o: test_solve.c $(INC)/hhmpcsolve.h
	gcc $(FLAGS) -c -o ./tests/test_solve.o ./tests/test_solve.c

test_solve.c: ./tests/test_solve.test
	checkmk tests/test_solve.test >tests/test_solve.c


math-test: test_math.o hhmpcmath.o hhmpcusefull.o mpcincmtxops.o
	gcc -o ./tests/math-test hhmpcmath.o hhmpcusefull.o mpcincmtxops.o ./tests/test_math.o -lcheck -lpthread -lrt -lm

test_math.o: test_math.c $(INC)/hhmpcmath.h
	gcc $(FLAGS) -c -o ./tests/test_math.o ./tests/test_math.c

test_math.c: ./tests/test_math.test
	checkmk tests/test_math.test >tests/test_math.c