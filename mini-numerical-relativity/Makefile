# mini-numerical-relativity Makefile
# Knowledge: Numerical methods for solving Einstein's equations

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS  = -lm
AR       = ar
ARFLAGS  = rcs

# Source files
SRC = src/nr_grid.c src/nr_utils.c src/nr_adm.c src/nr_bssn.c \
      src/nr_initial.c src/nr_horizon.c src/nr_wave.c src/nr_gauge.c

# Object files
OBJ = $(SRC:.c=.o)

# Test files
TESTS = tests/test_grid tests/test_adm tests/test_initial

# Example files
EXAMPLES = examples/example_schwarzschild examples/example_puncture \
           examples/example_wave_extract

# Library
LIB = libnr.a

.PHONY: all clean test examples count

all: $(LIB) test examples

# Static library
$(LIB): $(OBJ)
	$(AR) $(ARFLAGS) $@ $^

# Object compilation
src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

# Test compilation
tests/test_grid: tests/test_grid.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lnr $(LDFLAGS) -o $@

tests/test_adm: tests/test_adm.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lnr $(LDFLAGS) -o $@

tests/test_initial: tests/test_initial.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lnr $(LDFLAGS) -o $@

# Example compilation
examples/example_schwarzschild: examples/example_schwarzschild.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lnr $(LDFLAGS) -o $@

examples/example_puncture: examples/example_puncture.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lnr $(LDFLAGS) -o $@

examples/example_wave_extract: examples/example_wave_extract.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lnr $(LDFLAGS) -o $@

# Run tests
test: $(TESTS)
	@echo "=== Running Tests ==="
	@for t in $(TESTS); do \
		echo; \
		./$$t || exit 1; \
	done
	@echo; echo "=== All tests passed ==="

# Run examples
examples: $(EXAMPLES)

# Count lines in include/ + src/
count:
	@echo "include/ lines:"
	@find include -name "*.h" -exec cat {} + | wc -l
	@echo "src/ lines:"
	@find src -name "*.c" -exec cat {} + | wc -l
	@echo "---"
	@echo "Total (include + src):"
	@find include src -name "*.h" -o -name "*.c" | xargs cat | wc -l

# Clean
clean:
	rm -f src/*.o
	rm -f $(LIB)
	rm -f $(TESTS)
	rm -f $(EXAMPLES)
