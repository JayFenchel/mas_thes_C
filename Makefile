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
