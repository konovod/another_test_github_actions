ifndef (ARCH)
ARCH=default
endif

ifndef (ARCHFLAGS)
ifeq ($(ARCH),x86)
ARCHFLAGS=-m32 -march=i386 -mtune=generic
else ifeq ($(ARCH),x64)
ARCHFLAGS=-m64 -march=x86-64 -mtune=generic
else ifeq ($(ARCH),native)
ARCHFLAGS=-march=native
else
ARCHFLAGS=
endif
endif

ifeq ($(OS),Windows_NT)
ARCHFLAGS+=-mstackrealign -DDBC_FFT_ENABLE_MINGW_SIMD
OUTPUT=check.exe
endif

ifeq ($(OUTPUT),)
OUTPUT=check
endif

.PHONY: cpp c preprocessed compress clean

cpp:
	g++ -std=c++11 -Wall -Wextra -Wconversion -Wsign-conversion $(ARCHFLAGS) -DDBC_FFT_CACHE_CPU_DETECTION -DUSE_FLOAT128 -DUSE_FIXEDPOINT -fext-numeric-literals -Ofast -s -o $(OUTPUT) check.cpp -lquadmath

c:
	gcc -std=c99 -Wpedantic -Wall -Wextra -Wconversion -Wsign-conversion $(ARCHFLAGS) -DDBC_FFT_CACHE_CPU_DETECTION -DUSE_FLOAT128 -O3 -s -o $(OUTPUT) check.c -lm -lquadmath

preprocessed:
	gcc -P -E -DPREPROCESSED -std=c99 -Wpedantic -Wall -Wextra -Wconversion -Wsign-conversion $(ARCHFLAGS) -o preprocessed.c check.c

compress:
	cp check.c check.cpp check.inc test.inc dbc_fft.h Makefile ./release
	zip -r9 dbc_fft.zip release
	rm ./release/*

clean:
	rm preprocessed.c
	rm dbc_fft.zip
