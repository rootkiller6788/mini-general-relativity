CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -O2 -g
LDFLAGS  = -lm
INCLUDES = -Iinclude

SRC_DIR  = src
TEST_DIR = tests
EXAM_DIR = examples
BUILD_DIR = build

SRCS = $(SRC_DIR)/dg_manifold.c $(SRC_DIR)/dg_tensor.c $(SRC_DIR)/dg_metric.c \
       $(SRC_DIR)/dg_connection.c $(SRC_DIR)/dg_curvature.c \
       $(SRC_DIR)/dg_geodesic.c $(SRC_DIR)/dg_forms.c \
       $(SRC_DIR)/dg_lie.c $(SRC_DIR)/dg_tetrad.c $(SRC_DIR)/dg_gravwave.c

OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

.PHONY: all clean test examples dirs

all: dirs $(OBJS)
	@echo "mini-differential-geometry library compiled"

dirs:
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

test: dirs $(OBJS)
	@if [ -f $(TEST_DIR)/test_core.c ]; then \
		$(CC) $(CFLAGS) $(INCLUDES) $(TEST_DIR)/test_core.c $(OBJS) $(LDFLAGS) -o $(BUILD_DIR)/test_core && \
		$(BUILD_DIR)/test_core; \
	else \
		echo "No tests found"; \
	fi

examples: dirs $(OBJS)
	@for e in $(EXAM_DIR)/*.c; do \
		if [ -f "$$e" ]; then \
			name=$$(basename "$$e" .c); \
			$(CC) $(CFLAGS) $(INCLUDES) "$$e" $(OBJS) $(LDFLAGS) -o $(BUILD_DIR)/example_$$name && \
			echo "Running example_$$name..." && \
			$(BUILD_DIR)/example_$$name; \
		fi; \
	done

clean:
	rm -rf $(BUILD_DIR)

count:
	@echo "include/ lines:"; cat include/*.h | wc -l
	@echo "src/ lines:"; cat src/*.c | wc -l
	@echo "Total:"; cat include/*.h src/*.c | wc -l
