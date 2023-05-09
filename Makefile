CONF = config.txt
all: 2body plot

2body: config.txt main.cpp | build
	g++ -o 2body main.cpp
	./2body
	
build:
	mkdir -p build

build_anim:
	mkdir -p videos/build

plot: config.txt 2body plot.py | build
	python plot.py

animation: config.txt 2body plot | build_anim
	python animation.py
	
clean:
	rm -rf build
	rm 2body
	rm -rf videos/build

.PHONY: all clean plot