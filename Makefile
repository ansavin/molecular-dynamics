all: main.cpp MDSys.cpp
	g++ -std=c++14 -Wall -Wpedantic -Wno-misleading-indentation -O3 -flto -fopenmp -o MD main.cpp MDSys.cpp

debug: main.cpp MDSys.cpp
	g++ -std=c++14 -Wall -Wpedantic -Wno-misleading-indentation -O0 -g -fopenmp -o MD main.cpp MDSys.cpp
    
clean:
	rm -f MD
