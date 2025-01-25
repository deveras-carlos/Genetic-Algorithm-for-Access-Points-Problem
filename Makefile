PROGRAM_NAME = Cellular_PGA_Access_Point #Cellular Memetic Parallel Genetic Algorithm (CMPGA) for Access Point problem

PLATFORM ?= linux
ifeq ($(PLATFORM),windows64_gcc)

CC=x86_64-w64-mingw32-gcc
CFLAGS += -DWINDOWS64 -static

else

CC = gcc
CFLAGS += -DUNIX

endif


# Início do makefile
PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))


SOURCE   = $(PROJECT_ROOT)/source
INCLUDE  = $(PROJECT_ROOT)/include
OBJECTS  = \
		$(SOURCE)/main.c \
		$(SOURCE)/genetic_algorithm.c
BUILD_DIR= ./

LIBS = -fopenmp -lm

CFLAGS  += $(LIBS)

.PHONY: all program clean

all: program

program:
	@mkdir -p $(BUILD_DIR)
	$(CC) -o $(BUILD_DIR)/$(PROGRAM_NAME) -I $(INCLUDE) $(CFLAGS) $(OBJECTS)
	
clean:
	rm -fr $(BUILD_DIR)/$(PROGRAM_NAME) build/*
