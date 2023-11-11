release:
	gcc -s -O2 src/main.c src/input.c src/sim.c -o bin/dpsim -lm -Wall -Werror

run: build-debug
	bin/dpsim-debug

build-debug:
	gcc src/main.c src/input.c src/sim.c -o bin/dpsim-debug -lm -O0

clean:
	rm bin/*


optimized:
	echo "Building with completely unnecessary optimizations"
	gcc -lm -Ofast -s -flto -funroll-loops -finline-functions src/main.c src/input.c src/sim.c-o bin/dpsim

pedantic:
	gcc src/main.c src/input.c src/sim.c -o bin/dpsim-debug -lm -std=iso9899:1990 -pedantic -Wall -Werror
