CC       := gcc
CFLAGS   := -std=c11 -Wall -Wextra -O2 -fPIC
LDFLAGS  := -lm

INCDIR   := include
SRCDIR   := src
TSTDIR   := tests
EXADIR   := examples
BINDIR   := bin

SRCS     := $(wildcard $(SRCDIR)/*.c)
OBJS     := $(patsubst $(SRCDIR)/%.c,$(BINDIR)/%.o,$(SRCS))
TEST_SRCS := $(wildcard $(TSTDIR)/*.c)
TEST_BINS := $(patsubst $(TSTDIR)/%.c,$(BINDIR)/test_%,$(TEST_SRCS))
EXA_SRCS := $(wildcard $(EXADIR)/*.c)
EXA_BINS := $(patsubst $(EXADIR)/%.c,$(BINDIR)/example_%,$(EXA_SRCS))

.PHONY: all clean test examples count

all: $(OBJS)

$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/%.o: $(SRCDIR)/%.c | $(BINDIR)
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

$(BINDIR)/test_%: $(TSTDIR)/%.c $(OBJS) | $(BINDIR)
	$(CC) $(CFLAGS) -I$(INCDIR) $< $(OBJS) $(LDFLAGS) -o $@

$(BINDIR)/example_%: $(EXADIR)/%.c $(OBJS) | $(BINDIR)
	$(CC) $(CFLAGS) -I$(INCDIR) $< $(OBJS) $(LDFLAGS) -o $@

test: $(TEST_BINS)
	@passed=0; failed=0; \
	for t in $(TEST_BINS); do \
		echo "[RUN] $$t"; \
		if $$t; then passed=$$((passed+1)); else failed=$$((failed+1)); fi; \
		echo ""; \
	done; \
	echo "=== $$((passed+failed)) tests: $$passed passed, $$failed failed ==="; \
	if [ $$failed -ne 0 ]; then exit 1; fi

examples: $(EXA_BINS)
	@for e in $(EXA_BINS); do \
		echo "[EXE] $$e"; \
		$$e; \
		echo ""; \
	done

clean:
	rm -rf $(BINDIR)

count:
	@echo "Headers:" && wc -l $(INCDIR)/*.h
	@echo "Sources:" && wc -l $(SRCDIR)/*.c
	@echo "Total include+src:" && cat $(INCDIR)/*.h $(SRCDIR)/*.c | wc -l
