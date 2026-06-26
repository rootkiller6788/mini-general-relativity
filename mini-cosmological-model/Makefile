CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 -Iinclude
LDFLAGS = -lm

SRC_DIR = src
INC_DIR = include
TEST_DIR = tests
EXAMPLE_DIR = examples
BUILD_DIR = build

SRCS = $(SRC_DIR)/cosmo_model.c $(SRC_DIR)/flrw.c $(SRC_DIR)/friedmann.c        $(SRC_DIR)/distances.c $(SRC_DIR)/perturbations.c        $(SRC_DIR)/power_spectrum.c $(SRC_DIR)/transfer.c

OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

TEST_SRCS = $(TEST_DIR)/test_friedmann.c $(TEST_DIR)/test_distances.c
TEST_BINS = $(TEST_SRCS:$(TEST_DIR)/%.c=$(BUILD_DIR)/%)

EXAMPLE_SRCS = $(EXAMPLE_DIR)/age_of_universe.c                $(EXAMPLE_DIR)/distance_redshift.c                $(EXAMPLE_DIR)/cosmic_evolution.c
EXAMPLE_BINS = $(EXAMPLE_SRCS:$(EXAMPLE_DIR)/%.c=$(BUILD_DIR)/examples/%)

.PHONY: all clean test examples lib

all: lib test examples

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/examples

lib: $(BUILD_DIR) $(OBJS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(INC_DIR)/*.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_%: $(TEST_DIR)/test_%.c $(OBJS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/examples/%: $(EXAMPLE_DIR)/%.c $(OBJS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

test: $(TEST_BINS)
	@echo === Running Tests ===
	$(BUILD_DIR)/test_friedmann
	$(BUILD_DIR)/test_distances
	@echo === All Tests Complete ===

examples: $(EXAMPLE_BINS)
	@echo === Examples Built ===

clean:
	rm -rf $(BUILD_DIR)

help:
	@echo mini-cosmological-model
	@echo Targets: all lib test examples clean help
