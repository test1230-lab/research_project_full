CXX = g++-16

TARGET = fit

SRC = fit.cpp Chebyshev.cpp Legendre.cpp TailFit.cpp
OBJ = $(SRC:.cpp=.o)
DEP = $(OBJ:.o=.d)

# Safe math optimizations: faster FP without -ffast-math, so IEEE NaN/Inf
# semantics are preserved (only errno setting and FP traps are given up).
MATHFLAGS = -fno-math-errno -fno-trapping-math -ffp-contract=fast

# Link-time optimization (whole-program); the link step recompiles, so these
# flags are needed on both the compile and link commands.
LTOFLAGS = -flto=auto

CXXFLAGS = -O3 -std=c++23 -Wall -march=native -mtune=native -fopenmp \
           $(MATHFLAGS) $(LTOFLAGS) \
           -I/usr/include/eigen3 -MMD -MP

LDFLAGS = -O3 -march=native -mtune=native -fopenmp \
          $(MATHFLAGS) $(LTOFLAGS) -lquadmath

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include $(DEP)

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)

rebuild: clean all

.PHONY: all clean rebuild
