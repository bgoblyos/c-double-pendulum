release:
	gcc -s -O2 src/main.c src/input.c -o bin/dpsim -lm

run: build-debug
	bin/dpsim-debug

build-debug:
	gcc src/main.c src/input.c -o bin/dpsim-debug -lm -O0

clean:
	rm bin/*


optimized:
	echo "Building with completely unnecessary optimizations"
	gcc -lm -Ofast -s -flto -funroll-loops -finline-functions src/main.c src/input.c -o bin/dpsim
