EIGENLIB = -I/usr/include/eigen3

smoothing: main4.cpp io.cpp
	g++ -Wall -ggdb $(EIGENLIB) -o smoothing main4.cpp io.cpp

test-assignment-image: smoothing
	./smoothing -v bunny.obj threeOneFifty_explicit.obj 3 150
	./smoothing -v -i 0.001 bunny.obj threeOneFifty_implicit.obj 3 150

test-basic: smoothing
	./smoothing bunny.obj oneOne-Basic.obj 1 1
	./smoothing bunny.obj threeOne-Basic.obj 3 1
	./smoothing bunny.obj oneFifty-Basic.obj 1 50
	./smoothing bunny.obj threeFifty-Basic.obj 3 50

test-cotangent: smoothing
	./smoothing -c bunny.obj oneOne-Cotangent.obj 1 1
	./smoothing -c bunny.obj threeOne-Cotangent.obj 3 1
	./smoothing -c bunny.obj oneFifty-Cotangent.obj 1 50
	./smoothing -c bunny.obj threeFifty-Cotangent.obj 3 50

test-biharmonic: smoothing
	./smoothing -b bunny.obj oneOne-Biharmonic.obj 1 1
	./smoothing -b bunny.obj threeOne-Biharmonic.obj 3 1
	./smoothing -b bunny.obj oneFifty-Biharmonic.obj 1 50
	./smoothing -b bunny.obj threeFifty-Biharmonic.obj 3 50

test-preserve: smoothing
	./smoothing -v bunny.obj oneOne-Preserve.obj 1 1
	./smoothing -v bunny.obj threeOne-Preserve.obj 3 1
	./smoothing -v bunny.obj oneFifty-Preserve.obj 1 50
	./smoothing -v bunny.obj threeFifty-Preserve.obj 3 50

test-implicit: smoothing
	./smoothing -i 0.001 bunny.obj oneOne-Implicit.obj 1 1
	./smoothing -i 0.001 bunny.obj threeOne-Implicit.obj 3 1
	./smoothing -i 0.001 bunny.obj oneFifty-Implicit.obj 1 50
	./smoothing -i 0.001 bunny.obj threeFifty-Implicit.obj 3 50

test-loop: smoothing
	./smoothing -s 1 bunny.obj oneOne-Division.obj 1 1
	./smoothing -s 1 bunny.obj threeOne-Division.obj 3 1
	./smoothing -s 1 bunny.obj oneFifty-Division.obj 1 50
	./smoothing -s 1 bunny.obj threeFifty-Division.obj 3 50

test-all: test-basic test-cotangent test-biharmonic test-preserve test-implicit test-loop test-assignment-image
	echo "Complete"

clean: smoothing
	rm smoothing

