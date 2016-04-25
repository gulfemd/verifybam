VERIFYBAM_VERSION := "0.1-alpha"
VERIFYBAM_UPDATE := "September 21, 2015"
VERIFYBAM_DEBUG := 1
BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS =  -O3 -g -I htslib -DVERIFYBAM_VERSION=\"$(VERIFYBAM_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DVERIFYBAM_UPDATE=\"$(VERIFYBAM_UPDATE)\" -DVERIFYBAM_DEBUG=$(VERIFYBAM_DEBUG)
LDFLAGS = htslib/libhts.a -lz -lm -lpthread
SOURCES = verifybam.c cmdline.c common.c processbam.c fastatools.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = verifybam
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) *.o *~ \#*

libs:
	make -C htslib

install:
	cp verifybam $(INSTALLPATH)
