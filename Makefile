CHECKBAM_VERSION := "0.1-alpha"
CHECKBAM_UPDATE := "September 21, 2015"
CHECKBAM_DEBUG := 1
BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS =  -O3 -g -I htslib -I vhc -I vhsc -DCHECKBAM_VERSION=\"$(CHECKBAM_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DCHECKBAM_UPDATE=\"$(CHECKBAM_UPDATE)\" -DCHECKBAM_DEBUG=$(CHECKBAM_DEBUG)
LDFLAGS = htslib/libhts.a -lz -lm -lpthread
SOURCES = checkbam.c cmdline.c common.c processbam.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = checkbam
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) *.o *~

libs:
	make -C htslib

install:
	cp checkbam $(INSTALLPATH)
