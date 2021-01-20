# *****************************************************
CPP = gcc # Place path of gcc here !!
CPPFLAGS = -Wall -O3 -fopenmp
# ****************************************************

sol: sol.o utils.o sprsmatrix.o ordering.o
	$(CPP) $(CPPFLAGS) -o main sol.o utils.o sprsmatrix.o ordering.o

sol.o: sol.cpp
	$(CPP) $(CPPFLAGS) -c sol.cpp

ordering.o: ordering.cpp 
	$(CPP) $(CPPFLAGS) -c ordering.cpp

utils.o: utils.cpp 
	$(CPP) $(CPPFLAGS) -c utils.cpp

sprsmatrix.o: sprsmatrix.cpp
	$(CPP) $(CPPFLAGS) -c sprsmatrix.cpp 
	
clean:
	rm sol.o sprsmatrix.o utils.o ordering.o main
