#include "sph.h"
#include <sys/stat.h>

void options ( int argc, char **argv );
void generate_output ( FILE *in );


char *filename = NULL;
size_t filesize = 0;
int_t n_field = 0;


bool
    print_all = false;

int
main ( int argc, char **argv )
{
    options ( argc, argv );
    if ( filename != NULL )
    {
        FILE *in = fopen ( filename, "r" );
        generate_output ( in );
        fclose ( in );
        free ( filename );
    }
    exit ( EXIT_SUCCESS );
}


void
generate_output ( FILE *in )
{
    for ( int p=0; p<n_field; p++ )
    {
        particle_t p;
        int readc = fread ( &p, sizeof(particle_t), 1, in );
        printf ( "%ld %e %e\n", p.idx, p.x[0], p.x[1] );
    }
}


void
options ( int argc, char **argv )
{
    int o;
    struct stat file_stat;
    while ( (o = (getopt (argc, argv, "ahf:") )) != -1 )
    {
        switch (o)
        {
            case 'h': printf ( "Help? This is help!\n" ); break;
            case 'a': print_all = true; break;
            case 'f':
                filename = strdup ( optarg );
                if ( stat ( filename, &file_stat ) )
                {
                    fprintf ( stderr, "%s: Could not stat file '%s'\n",
                        argv[0], filename
                    );
                    free ( filename );
                    exit ( EXIT_FAILURE );
                }
                else
                {
                    filesize = file_stat.st_size;
                    n_field = filesize / sizeof(particle_t);
                    if ( filesize % sizeof(particle_t) )
                    {
                        fprintf ( stderr,
                            "%s: Warning, file '%s' is not a multiple "
                            "of particle structure size\n",
                            argv[0], filename
                        );
                    }
                }
                break;
        }
    }
}
