
CC = gcc
OUTPUT_DIR=../bin/

all: ENCprime SeqCount

ENCprime: ENCprime.o
	$(CC) $(OPTS) -o $(OUTPUT_DIR)/ENCprime ENCprime.o -lm
ENCprime.o: ENCprime.c
	$(CC) $(OPTS) -c ENCprime.c

SeqCount: SeqCount.o
	$(CC) $(OPTS) -o $(OUTPUT_DIR)/SeqCount SeqCount.o
SeqCount.o: SeqCount.c
	$(CC) $(OPTS) -c SeqCount.c

clean:
	rm -f *.o $(OUTPUT_DIR)/ENCprime
	rm -f *.o $(OUTPUT_DIR)/SeqCount