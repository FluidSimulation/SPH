#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef int64_t int_t;
typedef double real_t;

int
main ( int argc, char **argv )
{
    if ( argc < 2 )
    {
        fprintf ( stderr, "No file name\n" );
        exit ( EXIT_FAILURE );
    }
    FILE *in = fopen ( argv[1], "r" );
    real_t triple[3];
    while ( fread ( triple, sizeof(real_t), 3, in ) && !feof(in) )
    {
        int_t idx = (int_t) triple[0];
        printf ( "%.5ld %e %e\n", idx, triple[1], triple[2] );
    }
    fclose ( in );
}
