# If Sherpa-config is not in your $PATH, change this variable to point to your
# Sherpa installation directory manually
SHERPA_PREFIX ?= $(shell Sherpa-config --prefix)

SOURCES=Model.C

CC=g++
CXXFLAGS=-g -O2 -std=c++11 -pedantic -Wall -fPIC -DPIC

CXXFLAGS+= -I$(shell $(SHERPA_PREFIX)/bin/Sherpa-config --incdir)
LDFLAGS=$(shell $(SHERPA_PREFIX)/bin/Sherpa-config --ldflags)

OBJECTS=$(SOURCES:.C=.o)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LIB_SUFFIX=dylib
else
	LIB_SUFFIX=so
endif

all: $(SOURCES) $(OBJECTS) 
	$(CC) -shared -o libSherpaSMGGV.$(LIB_SUFFIX) $(OBJECTS) $(LDFLAGS)

.C.o:
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJECTS) libSherpaSMGGV.$(LIB_SUFFIX)
