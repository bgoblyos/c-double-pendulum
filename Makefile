run: build
	bin/main

build:
	gcc src/main.c -o bin/main -lm

clean:
	rm bin/main
