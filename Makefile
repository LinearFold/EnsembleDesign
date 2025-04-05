CLA=clang++
CXX=g++
CXXFLAGS=-std=c++11 -Ofast -DFINAL_CHECK -DSPECIAL_HP -fpermissive

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(UNAME_S), Linux)
    LIBNAME = libutils_linux64.so
else
    ifeq ($(UNAME_M), arm64)
         LIBNAME = libutils_Mac_M1.so
    else
         LIBNAME = libutils_Mac_x86.so
    endif
endif


LIBDIR = src/Utils/lib
LIB = $(LIBDIR)/$(LIBNAME)

SOURCES = src/ensemble_design.cpp 
DEPS = src/optimizer.cpp src/optimizer.h src/Utils/reader.h src/Utils/network.h src/Utils/codon.h src/Utils/utility_v.h src/Utils/common.h src/Utils/base.h
BIN = bin/EnsembleDesign

all: $(BIN)

$(BIN): $(SOURCES) $(DEPS) $(LIB) | bin
	@echo "Compiling binary $(BIN) using $(LIBNAME)..."
	mkdir -p bin
	cp $(LIB) $(LIBDIR)/libutils.so
	if $(CLA) $(CXXFLAGS) $(SOURCES) -o $(BIN) $(LIBDIR)/libutils.so; then \
		echo "Compiled with clang++; finished"; \
	else \
		if $(CXX) $(CXXFLAGS) $(SOURCES) -o $(BIN) $(LIBDIR)/libutils.so; then \
			echo "Compiled with g++; finished"; \
		else \
			echo "Compilation failed!"; \
			exit 1; \
		fi; \
	fi

bin:
	mkdir -p bin

.PHONY: clean

clean:
	rm -f $(BIN)
