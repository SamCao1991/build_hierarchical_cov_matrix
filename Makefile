Include = -I./h2lib/ -I./libmorton/include/
Lib = -L.
CC = g++
CCFlags = -m64 -std=c++11 -O3

test_Jian: test_Jian.o libh2.a
	$(CC) test_Jian.o $(Lib) -lh2 -o test_Jian

%.o: %.cpp
	$(CC) $(Include) $(CCFlags) -c $< -o $@

clean:
	rm *.o
