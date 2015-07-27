CC=gcc
CFLAGS="-Wall"

debug:clean
	$(CC) $(CFLAGS) -g -o master_thesis_c main.c
stable:clean
	$(CC) $(CFLAGS) -o master_thesis_c main.c
clean:
	rm -vfr *~ master_thesis_c
