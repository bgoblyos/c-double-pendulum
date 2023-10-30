run: build
	bin/main

build:
	gcc src/main.c src/input.c -o bin/main -lm

clean:
	rm bin/*
