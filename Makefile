# CC set in external env., machine specific settings in environment.sh
# CFLAGS config per compiler:

SCALE=1.0

CFLAGS+=${CFLAGS_${CC}}

CFLAGS_mpiicc=-DSCALE=${SCALE} -std=c99 -Iinclude -qopenmp -O2 #-g -O0 -ggdb -gdwarf-2 -g3
CFLAGS_mpicc=-DSCALE=${SCALE} -std=c99 -Iinclude -fopenmp -O2 #-g -O0 -ggdb -gdwarf-2 -g3

# Hash table library should be included/embedded for simplicity, isn't yet
LDFLAGS+=-Llib
LDLIBS+=-lm -ltlhash

# FFPMEG configuration, this needs to be set by env
# FFMPEG=${HOME}/tools/bin/ffmpeg

all: sph dat2txt cp2txt
sph: sph.c sph_io.o particle_hashtab.o lib/libtlhash.a
dat2txt: dat2txt.c
lib/libtlhash.a:
	${MAKE} -C lib
dambreak.mp4: plots
	${FFMPEG} ${FFMPEG_FLAGS} dambreak.mp4
plots: $(shell find plot/ -name '*.dat' | sed s/dat/png/g)
plot/%.txt: plot/%.dat
	@./cp2txt -f plot/$*.dat > plot/$*.txt
plot/%.png: plot/%.txt
	@./convert_dat.sh plot/$*.txt plot/$*.png ${SCALE} 2>&1 > /dev/null
.PHONY: clean plots
clean:
	-rm -f sph dat2txt cp2txt *.o
