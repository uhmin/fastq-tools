CFLAGS= -Wall -O -D_FILE_OFFSET_BITS=64
OBJECT= -O -D_FILE_OFFSET_BITS=64

all: fastqQtrim.exe

fastqQtrim.exe: fastqQtrim.o
	gcc $(CFLAGS) -o $@ $+
	chmod +x *.pl

clean:
	rm fastqQtrim.exe *.o

%.o: %.c
	gcc $(OBJECT) -c -o $@ $+
