all:
	g++ src_TCEM/main.cc -Wall  -std=c++0x -O3 src_TCEM/allocator.cc src_TCEM/utils.cc src_TCEM/itemGraph.cc src_TCEM/anyoption.cc src_TCEM/sfmt/SFMT.c  -o main_TCEM
