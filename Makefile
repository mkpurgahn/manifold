CC = cc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -I src_c/ -O2
LDFLAGS = -lm

SRCS = src_c/manifold_impl.c src_c/manifold_constructors.c src_c/manifold_properties.c src_c/manifold_api.c src_c/manifold_sdf.c src_c/manifold_boolean.c
OBJS = $(SRCS:.c=.o)

TEST_SRCS = src_c/test_manifold.c
TEST_BIN = test_manifold

.PHONY: all clean test

all: $(OBJS) $(TEST_BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(TEST_SRCS) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_SRCS) $(OBJS) $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -f $(OBJS) $(TEST_BIN)
