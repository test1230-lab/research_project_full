CXX = g++-13

CXXFLAGS = -O3 -std=c++23 -Wall -march=native -fopenmp -I/usr/include/eigen3 -DBOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT

dist1d: dist1d.cpp
	$(CXX) $(CXXFLAGS) dist1d.cpp -o dist1d

clean:
	rm -f dist1d