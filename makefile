CC=g++
CFLAGS=-g -Wall

p1t1.o: p1t1.cpp
	$(CC) $(CFLAGS) p1t1.cpp -o p1t1.o

run:
	@$(MAKE) && ./p1t1.o

clean:
	rm p1t1.o
