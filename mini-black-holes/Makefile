# Makefile for mini-black-holes
# Black hole physics library and tests
# make / make test / make clean

CC       = gcc
CFLAGS   = -Wall -Wextra -std=c11 -O2 -g
CPPFLAGS = -I./include
LDFLAGS  = -lm

SRC_DIR  = src
INC_DIR  = include
TST_DIR  = tests
BLD_DIR  = build

SRCS = src/black_hole_metrics.c src/black_hole_thermodynamics.c src/black_hole_dynamics.c src/black_hole_waves.c src/black_hole_advanced.c
OBJS = $(patsubst $(SRC_DIR)/%.c, $(BLD_DIR)/%.o, $(SRCS))

.PHONY: all test clean

all: test

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

$(BLD_DIR)/libblackholes.a: $(OBJS)
	ar rcs $@ $^

$(BLD_DIR)/test_black_holes: $(TST_DIR)/test_black_holes.c $(BLD_DIR)/libblackholes.a | $(BLD_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -L$(BLD_DIR) -lblackholes $(LDFLAGS) -o $@

test: $(BLD_DIR)/test_black_holes
	@echo "=== Running mini-black-holes tests ==="
	./$(BLD_DIR)/test_black_holes
	@echo "=== Tests complete ==="

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
