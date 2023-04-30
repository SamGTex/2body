all: 2body plot

2body: main.cpp | build
	g++ -o 2body main.cpp
	./2body
	
build:
	mkdir -p build

build_anim:
	mkdir -p videos/build

plot: 2body plot.py | build
	python plot.py

animation: 2body | build_anim
	python animation.py
	
clean:
	rm -rf build
	rm 2body
	rm -rf videos/build

.PHONY: all clean plot