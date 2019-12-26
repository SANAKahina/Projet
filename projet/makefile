CC=gcc
RM=rm -f

#flags
CFLAGS=-O3
LDFLAGS=-Wall -g -lm

#Files
SOURCE=main.c fonction.c 
OBJECT=$(SOURCE:.c=.o)


%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

main: $(OBJECT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJECT) main
