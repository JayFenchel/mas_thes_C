CC=gcc
CFLAGS="-Wall"

debug:clean
	$(CC) $(CFLAGS) -g -o main main.c
stable:clean
	$(CC) $(CFLAGS) -o main main.c
clean:
	rm -vfr *~ main
