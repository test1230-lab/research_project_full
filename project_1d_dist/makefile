CXX = g++-13

TARGET = dist1d

SRC = dist1d.cpp
OBJ = $(SRC:.cpp=.o)
DEP = $(OBJ:.o=.d)

# compile via .o with -MMD so header edits (mdarray.h, npy.hpp, libInterpolate)
# trigger rebuilds; the old direct compile-and-link rule missed them
CXXFLAGS = -O3 -march=native -mtune=native -fno-math-errno -std=c++23 -fopenmp \
           -Wall -I/usr/include/eigen3 -DBOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT \
           -MMD -MP

LDFLAGS = -O3 -march=native -mtune=native -fopenmp

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
