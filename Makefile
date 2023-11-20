CC=g++
SOURCES=src/main.cpp src/matrix_t.cpp

run: build
	@./bin/main

build:
	@$(CC) $(SOURCES) -o bin/main 

clean:
	@rm -rf *.o 
