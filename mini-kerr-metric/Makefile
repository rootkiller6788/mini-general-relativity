# Makefile for mini-kerr-metric
# Kerr Metric: Rotating Black Hole Solutions in General Relativity
#
# MIT 8.962 / Wald Ch.7 / Chandrasekhar Ch.6

CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99
CFLAGS += -Wno-unused-parameter -Wno-unused-but-set-variable
CFLAGS += -Wno-unused-variable -Wno-maybe-uninitialized
CFLAGS += -I./include
LDFLAGS = -lm

# Source files
SRCS = $(wildcard src/*.c)
OBJS = $(SRCS:.c=.o)

# Test files
TEST_SRCS = $(wildcard tests/*.c)
TEST_BINS = $(TEST_SRCS:.c=.test)

# Example files
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(EXAMPLE_SRCS:.c=.exe)

# Library
LIB = libkerr.a

.PHONY: all clean test examples lib

all: lib examples

# Build static library
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

# Compile C sources to objects
src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

# Build and run tests
test: $(TEST_BINS)
	@echo "=== Running Tests ==="
	@for t in $(TEST_BINS); do echo "--- $$t ---"; ./$$t; done
	@echo "=== All Tests Complete ==="

tests/%.test: tests/%.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lkerr $(LDFLAGS) -o $@

# Build examples
examples: $(EXAMPLE_BINS)

examples/%.exe: examples/%.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lkerr $(LDFLAGS) -o $@

# Clean
clean:
	rm -f src/*.o
	rm -f $(LIB)
	rm -f tests/*.test
	rm -f examples/*.exe

# Help
help:
	@echo "make lib       — Build static library libkerr.a"
	@echo "make test      — Build and run all tests"
	@echo "make examples  — Build all example programs"
	@echo "make clean     — Remove build artifacts"

# Quick compile test (just compile, don't link)
compile-test: $(OBJS)
	@echo "All source files compile successfully."
