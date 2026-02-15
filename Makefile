CC = cc
CXX = c++
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -I src_c/ -O2 -fno-strict-aliasing
CXXFLAGS = -std=c++17 -O2 -I build/_deps/clipper2-src/CPP/Clipper2Lib/include
LDFLAGS = -lm -lc++

SRCS = src_c/manifold_impl.c src_c/manifold_constructors.c src_c/manifold_properties.c src_c/manifold_api.c src_c/manifold_sdf.c src_c/manifold_boolean.c src_c/manifold_hull.c src_c/manifold_smooth.c src_c/manifold_csg_tree.c src_c/manifold_lazy_collider.c src_c/manifold_tree2d.c src_c/manifold_tri_dist.c
OBJS = $(SRCS:.c=.o)

CXX_SRCS = src_c/clipper2_c.cpp
CXX_OBJS = $(CXX_SRCS:.cpp=.o)

TEST_SRCS = src_c/test_manifold.c
TEST_BIN = test_manifold

CLIPPER2_LIB = build/_deps/clipper2-build/libClipper2.a

.PHONY: all clean test

all: $(OBJS) $(CXX_OBJS) $(TEST_BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

src_c/clipper2_c.o: src_c/clipper2_c.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BIN): $(TEST_SRCS) $(OBJS) $(CXX_OBJS) $(CLIPPER2_LIB)
	$(CC) $(CFLAGS) -o $@ $(TEST_SRCS) $(OBJS) $(CXX_OBJS) $(CLIPPER2_LIB) $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -f $(OBJS) $(CXX_OBJS) $(TEST_BIN)
