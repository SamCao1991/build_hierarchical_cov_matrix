Include = -I./libmorton/include/
Lib = -L. -lh2  -llapack -lblas -lgfortran -lcairo -lm 
CC = g++
CCFlags = -m64 -std=c++11 -O3

test_Jian: test_Jian.o
	$(CC) test_Jian.o $(Lib) -o test_Jian2

%.o: %.cpp
	$(CC) $(Include) $(CCFlags) -c $< -o $@

clean:
	rm *.o
