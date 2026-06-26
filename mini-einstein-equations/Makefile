# Makefile for mini-einstein-equations
# Reference: MIT 8.962, GR module of mini-pure-physics

CC       := gcc
CFLAGS   := -Wall -Wextra -std=c11 -O2 -g
INC      := -Iinclude
MATH_LIB := -lm

SRC_DIR  := src
INC_DIR  := include
TST_DIR  := tests
EXA_DIR  := examples
BLD_DIR  := build

# Source files
SRCS := $(SRC_DIR)/tensor.c \
        $(SRC_DIR)/metric.c \
        $(SRC_DIR)/curvature.c \
        $(SRC_DIR)/geodesic.c \
        $(SRC_DIR)/stress_energy.c \
        $(SRC_DIR)/einstein.c \
        $(SRC_DIR)/cosmology.c \
        $(SRC_DIR)/linearized.c \
        $(SRC_DIR)/coordinate.c \
        $(SRC_DIR)/bh_thermo.c

OBJS := $(SRCS:$(SRC_DIR)/%.c=$(BLD_DIR)/%.o)

# Test binary
TEST_SRC := $(TST_DIR)/test_all.c
TEST_BIN := $(BLD_DIR)/test_all

# Example binaries
EXAMPLES := $(BLD_DIR)/demo_schwarzschild \
            $(BLD_DIR)/demo_cosmology \
            $(BLD_DIR)/demo_gravitational_waves

.PHONY: all test clean examples

all: $(BLD_DIR) test examples

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

# Object file compilation
$(BLD_DIR)/%.o: $(SRC_DIR)/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

# Test binary
$(TEST_BIN): $(OBJS) $(TEST_SRC) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INC) $(TEST_SRC) $(OBJS) -o $@ $(MATH_LIB)

# Example binaries
$(BLD_DIR)/demo_schwarzschild: $(OBJS) $(EXA_DIR)/demo_schwarzschild.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INC) $(EXA_DIR)/demo_schwarzschild.c $(OBJS) -o $@ $(MATH_LIB)

$(BLD_DIR)/demo_cosmology: $(OBJS) $(EXA_DIR)/demo_cosmology.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INC) $(EXA_DIR)/demo_cosmology.c $(OBJS) -o $@ $(MATH_LIB)

$(BLD_DIR)/demo_gravitational_waves: $(OBJS) $(EXA_DIR)/demo_gravitational_waves.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INC) $(EXA_DIR)/demo_gravitational_waves.c $(OBJS) -o $@ $(MATH_LIB)

# Build and run tests
test: $(TEST_BIN)
	@echo "=== Running mini-einstein-equations tests ==="
	@./$(TEST_BIN)
	@echo "=== Tests complete ==="

# Build all examples
examples: $(EXAMPLES)

# Run all examples
run-examples: examples
	@echo "=== Running examples ==="
	@echo "\n[1/3] Schwarzschild Geodesics"
	@./$(BLD_DIR)/demo_schwarzschild
	@echo "\n[2/3] FLRW Cosmology"
	@./$(BLD_DIR)/demo_cosmology
	@echo "\n[3/3] Gravitational Waves"
	@./$(BLD_DIR)/demo_gravitational_waves

# Line count
count:
	@echo "include/ + src/ lines:"
	@find include src -name "*.h" -o -name "*.c" -o -name "*.lean" | xargs cat | wc -l

# Clean build artifacts
clean:
	rm -rf $(BLD_DIR)

# Full clean including docs
distclean: clean
	rm -f *.o *.exe

# Help
help:
	@echo "Targets:"
	@echo "  all          - Build and run tests + examples"
	@echo "  test         - Build and run test suite"
	@echo "  examples     - Build example programs"
	@echo "  run-examples - Build and run all examples"
	@echo "  count        - Count lines in include/ + src/"
	@echo "  clean        - Remove build artifacts"
