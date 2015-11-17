CXXFLAGS=--std=c++11 -fPIC -Wall -O3
SOFLAGS=-shared


all: acgtrie/libcpp_acgtrie.so


acgtrie/libcpp_acgtrie.so: acgtrie/cpp_acgtrie.cpp
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^
