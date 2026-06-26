# Makefile for mini-schwarzschild
CC       = gcc
CFLAGS   = -Wall -Wextra -std=c11 -O2 -g
CPPFLAGS = -I./include
LDFLAGS  = -lm

SRC_DIR  = src
INC_DIR  = include
TST_DIR  = tests
EX_DIR   = examples
BLD_DIR  = build

SRCS = src/schwarzschild_defs.c src/schwarzschild_tensors.c src/schwarzschild_geodesics.c src/schwarzschild_observables.c src/schwarzschild_numerical.c src/schwarzschild_extensions.c
OBJS = $(patsubst $(SRC_DIR)/%.c, $(BLD_DIR)/%.o, $(SRCS))

.PHONY: all test clean examples lines

all: test

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

$(BLD_DIR)/libschwarzschild.a: $(OBJS)
	ar rcs $@ $^

$(BLD_DIR)/test_schwarzschild: $(TST_DIR)/test_schwarzschild.c $(BLD_DIR)/libschwarzschild.a | $(BLD_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -L$(BLD_DIR) -lschwarzschild $(LDFLAGS) -o $@

test: $(BLD_DIR)/test_schwarzschild
	@echo "=== Running mini-schwarzschild tests ==="
	./$(BLD_DIR)/test_schwarzschild
	@echo "=== Tests complete ==="

examples: $(BLD_DIR)/libschwarzschild.a
	$(CC) $(CFLAGS) $(CPPFLAGS) $(EX_DIR)/mercury_precession.c -L$(BLD_DIR) -lschwarzschild $(LDFLAGS) -o $(BLD_DIR)/mercury_precession
	$(CC) $(CFLAGS) $(CPPFLAGS) $(EX_DIR)/light_deflection.c -L$(BLD_DIR) -lschwarzschild $(LDFLAGS) -o $(BLD_DIR)/light_deflection
	$(CC) $(CFLAGS) $(CPPFLAGS) $(EX_DIR)/black_hole_basics.c -L$(BLD_DIR) -lschwarzschild $(LDFLAGS) -o $(BLD_DIR)/black_hole_basics
	@echo "=== Running examples ==="
	./$(BLD_DIR)/mercury_precession
	./$(BLD_DIR)/light_deflection
	./$(BLD_DIR)/black_hole_basics

clean:
	rm -rf $(BLD_DIR)

lines:
	@echo "include/ lines:"
	@cat include/*.h | wc -l
	@echo "src/*.c lines:"
	@cat src/*.c | wc -l
	@echo "src/*.lean lines:"
	@cat src/*.lean | wc -l
	@echo "--- Total include/ + src/:"
	@(cat include/*.h src/*.c src/*.lean) | wc -l
