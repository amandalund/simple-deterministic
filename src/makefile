# User options

COMPILER = gnu

# Program and source code

PROGRAM = transport

HEADERS = header.h

SOURCE = \
main.c \
utils.c \
io.c \
initialize.c \
solvers.c

OBJECTS = $(SOURCE:.c=.o)

# Set flags

CFLAGS = -Wall -g
LDFLAGS = -lm

ifeq ($(COMPILER),gnu)
  CC = gcc
endif

# Targets to build

$(PROGRAM): $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(PROGRAM)
