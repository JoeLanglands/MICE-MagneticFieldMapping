CC = g++
PYLIBPATH = $(shell python-config --exec-prefix)/lib
LIB = -L$(PYLIBPATH) $(shell python-config --libs) -lboost_python
OPTS = $(shell python-config --include) -O2 -fPIC

default: makefield_cpp.so


makefield_cpp.so: makefield_cpp.o
	$(CC) $(LIB) -Wl,-rpath,$(PYLIBPATH) -shared $< -o $@

makefield_cpp.o: makefield_cpp.cpp makefile
	$(CC) $(OPTS) -c $< -o $@

clean:
	rm -rf *.so *.o *~

.PHONY: default clean
