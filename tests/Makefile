CFLAGS=-Wall
LIBS=-lcheck -lpthread -lrt -lm
DIR=../include/

all: sum2ints

sum2ints: main.o hhmpctestfunc.o
	gcc -o sum2ints main.o hhmpctestfunc.o

main.o: main.c $(DIR)hhmpctestfunc.h
	gcc $(CFLAGS) -c main.c

hhmpctestfunc.o: ../hhmpctestfunc.c $(DIR)hhmpctestfunc.h
	gcc $(CFLAGS) -c ../hhmpctestfunc.c

test: sum2ints-test
	./sum2ints-test

sum2ints-test: impl_test.o hhmpctestfunc.o
	gcc -o sum2ints-test hhmpctestfunc.o impl_test.o $(LIBS)

hhmpctestfunc-test.o: impl_test.c $(DIR)hhmpctestfunc.h
	gcc $(CFLAGS) -c impl_test.c

hhmpctestfunc-test.c: impl_test.check
	~/checkmk impl_test.check >impl_test.c