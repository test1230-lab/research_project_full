CXX = g++

TARGET = main2

SRC = main2.cpp ElectricField.cpp distrib2d.cpp
OBJ = $(SRC:.cpp=.o)
DEP = $(OBJ:.o=.d)

CXXFLAGS = -O3 -march=native -mtune=native -fno-math-errno -ffp-contract=fast \
           -std=c++23 -fopenmp -flto=auto \
           -Wall -MMD -MP

# -flto redoes codegen at link time, so the link line needs the same opt flags
LDFLAGS = -O3 -march=native -mtune=native -fno-math-errno -ffp-contract=fast \
          -fopenmp -flto=auto

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
