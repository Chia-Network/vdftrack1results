CFLAGS += -O3 -g -Wall -Wno-unused-function -pthread
LDFLAGS += -pthread

LOADLIBES += -lflint -lmpfr -lgmp -lm

all: main

main: main.o qfb.o reduce.o

clean:
	rm -f *.o main
