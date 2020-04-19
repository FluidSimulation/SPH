# Shell script to contain environment settings for various systems
# Usage: 'source environment.sh <system>'
SYSTEM=$1
case $SYSTEM in
idun)
    module load intel
    export CC=mpiicc
    export KMP_AFFINITY=verbose, compact
    export CFLAGS+=" -DWITH_MPIIO"
    export FFMPEG=${HOME}/tools/bin/ffmpeg
    export FFMPEG_FLAGS="-y -r 25 -i plot/%4d.png -b:v 16384k"
    ;;
vilje)
    module load intelcomp/18.0.1 mpt/2.14
    export KMP_AFFINITY=verbose,compact
    export CC=mpicc
    export CFLAGS+=" -DWITH_MPIIO"
    export FFMPEG=${HOME}/tools/bin/ffmpeg
    export FFMPEG_FLAGS="-y -r 25 -i plot/%4d.png -b:v 16384k"
    ;;
archer)
    module switch PrgEnv-cray PrgEnv-intel
    export CC=cc
    export CFLAGS+="-std=c99 -Iinclude -DWITH_MPIIO -qopenmp"
    ;;
local)
    export CC=mpicc
    export FFMPEG=ffmpeg
    export FFMPEG_FLAGS="-y -r 25 -i plot/%4d.png -tune grain -b:v 16384k"
    ;;
macbook)
    export CC=mpicc;
    export OMPI_CC=gcc-7
    export FFMPEG=ffmpeg
    export FFMPEG_FLAGS="-y -r 25 -i plot/%4d.png -tune grain -b:v 16384k -pix_fmt yuv420p"
    ;;
*)
    echo "Environment not predefined for system '${SYSTEM}'"
    ;;
esac
