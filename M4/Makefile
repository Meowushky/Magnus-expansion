CC=gcc
CFLAGS=-Wall -Wextra -pedantic -ggdb
LDFLAGS=-lm

all: M4

#цель
M4: M4.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

M4_1: M4_1.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

lua: lua.c
	$(CC) $(CFLAGS) $(LDFLAGS) -llua -o $@ $<
