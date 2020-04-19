#include "sph.h"
#include <errno.h>

static particle_t *checkpoint = NULL;

/* Internals of sph.h required for restarting from checkpoint file */
extern int_t min_iteration, max_iteration, checkpoint_frequency;
extern int_t n_field, n_global_field;
extern int_t n_capacity, n_pair_cap;
extern particle_t *list;
extern pair_t *pairs;
extern real_t subdomain[2];

void
collect_checkpoint ( void )
{
    int_t offsets[size];
    int_t n_local_cp = (n_global_field / size)
        + ( ( rank < (n_global_field % size) ) ? 1 : 0 );

    offsets[0] = 0;
    for ( int_t r=1; r<size; r++ )
    {
        offsets[r] = offsets[r-1] +
            (n_global_field / size) + ( ((r-1)<(n_global_field % size)) ? 1 : 0 );
    }

    if ( checkpoint == NULL )
        checkpoint = malloc ( n_local_cp * sizeof(particle_t) );

    for ( int_t pi=0; pi<n_global_field; pi++ )
    {
        /* Do I have particle #pi? */
        particle_t *p = NULL;
        lookup_particle ( pi, &p );
        /* I have it if p != NULL now */

        /* Who sould have particle #pi? */
        int target = size-1;
        for ( int r=0; r<(size-1); r++ )
            if ( pi >= offsets[r] && pi < offsets[r+1] )
                target = r;
        /* It belongs at target */

        /* At this point, the transfer pair is known to sender, and
         * everyone knows the receiver
         */
        if ( p != NULL && target == rank )  // "Send" to myself
        {
            memcpy ( &checkpoint[pi-offsets[rank]], p, sizeof(particle_t) );
        }
        else if ( p != NULL && target != rank ) // Send to target
        {
            MPI_Ssend ( p, sizeof(particle_t), MPI_BYTE,
                target, pi-offsets[target],
                MPI_COMM_WORLD
            );
        }
        else if ( p == NULL && target == rank ) // Recv from anywhere
        {
            MPI_Recv ( &checkpoint[pi-offsets[rank]],
                sizeof(particle_t), MPI_BYTE,
                MPI_ANY_SOURCE, pi-offsets[rank],
                MPI_COMM_WORLD, MPI_STATUS_IGNORE
            );
        }
        /* Rank offset-relative particle position as tag uniquely identifes
         * receive operations, barrier no longer necessary
         * MPI_Barrier ( MPI_COMM_WORLD );
         */
    }
}




/* NB allocation - this routine builds particle pointer lists
 *  on stack, go to dynamic allocation if it overflows
 */
#ifndef WITH_MPIIO
// Stupid POSIX-I/O, synchronizing iterations w. append to file
// Unsure how this interacts with parallel FS, but it avoids
// passing every particle to a single I/O-responsible rank
void
dump_state ( char *filename )
{
    int_t my_particles = n_particles();
    real_t data[my_particles][3];
    particle_t *actual_ptr[my_particles];
    list_particles ( &(actual_ptr[0]) );

    for ( int_t mp=0; mp<my_particles; mp++ )
    {
        particle_t *p = actual_ptr[mp];
        data[mp][0] = (real_t)p->idx;
        data[mp][1] = p->x[0];
        data[mp][2] = p->x[1];
    }
    /* Token ring synchronization for I/O */
    int token = rank, discard;
    FILE *out;
    switch ( rank )
    {
        case 0:
            out = fopen ( filename, "a" );
            fwrite ( data, 3*sizeof(real_t), my_particles, out );
            fclose ( out );
            MPI_Ssend ( &token, 1, MPI_INT, east, 0, MPI_COMM_WORLD );
            MPI_Recv ( &discard, 1, MPI_INT, west, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE
            );
            break;
        default:
            MPI_Recv ( &discard, 1, MPI_INT, west, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE
            );
            out = fopen ( filename, "a" );
            fwrite ( data, 3*sizeof(real_t), my_particles, out );
            fclose ( out );
            MPI_Ssend ( &token, 1, MPI_INT, east, 0, MPI_COMM_WORLD );
            break;
    }
    /* This barrier is probably not necessary,
     * border exchange also forces sync.
     */
    MPI_Barrier ( MPI_COMM_WORLD );
}


void
write_checkpoint ( char *filename )
{
    int_t offsets[size];
    int_t n_local_cp = (n_global_field / size)
        + ( ( rank < (n_global_field % size) ) ? 1 : 0 );

    offsets[0] = 0;
    for ( int_t r=1; r<size; r++ )
    {
        offsets[r] = offsets[r-1] +
            (n_global_field / size) + (((r-1)<(n_global_field % size))?1:0);
    }
//////////////////////////////////
    int token = rank, discard;
    FILE *out;
    switch ( rank )
    {
        case 0:
            out = fopen ( filename, "a" );
            fwrite ( checkpoint, sizeof(particle_t), n_local_cp, out );
            fclose ( out );
            MPI_Ssend ( &token, 1, MPI_INT, east, 0, MPI_COMM_WORLD );
            MPI_Recv ( &discard, 1, MPI_INT, west, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE
            );
            break;
        default:
            MPI_Recv ( &discard, 1, MPI_INT, west, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE
            );
            out = fopen ( filename, "a" );
            fwrite ( checkpoint, sizeof(particle_t), n_local_cp, out );
            fclose ( out );
            MPI_Ssend ( &token, 1, MPI_INT, east, 0, MPI_COMM_WORLD );
            break;
    }
    /* This barrier is probably not necessary,
     * border exchange also forces sync.
     */
    MPI_Barrier ( MPI_COMM_WORLD );
/*
//////////////////////////////////
    MPI_Info info;
    MPI_Info_create ( &info );
    MPI_Info_set ( info, "access_style", "write_once" );
    MPI_File out;
    MPI_File_open (
        MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,
        info, &out
    );
    MPI_Offset my_offset = offsets[rank];
    MPI_File_write_at_all ( out, my_offset*sizeof(particle_t),
        checkpoint, n_local_cp*sizeof(particle_t), MPI_BYTE,
        MPI_STATUS_IGNORE
    );
    MPI_File_close ( &out );
    MPI_Barrier ( MPI_COMM_WORLD );
*/
}
#else
// Fancy-pants MPI-I/O, for when the file system supports it
void
write_checkpoint ( char *filename )
{
    int_t offsets[size];
    int_t n_local_cp = (n_global_field / size)
        + ( ( rank < (n_global_field % size) ) ? 1 : 0 );

    offsets[0] = 0;
    for ( int_t r=1; r<size; r++ )
    {
        offsets[r] = offsets[r-1] +
            (n_global_field / size) + (((r-1)<(n_global_field % size))?1:0);
    }

    MPI_Info info;
    MPI_Info_create ( &info );
    MPI_Info_set ( info, "access_style", "write_once" );
//   MPI_Info_set ( info, "striping_factor", "4" );
    MPI_File out;
    MPI_File_open (
        MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,
        info, &out
    );
    MPI_Offset my_offset = offsets[rank];
    MPI_File_write_at_all ( out, my_offset*sizeof(particle_t),
        checkpoint, n_local_cp*sizeof(particle_t), MPI_BYTE,
        MPI_STATUS_IGNORE
    );
    /*
    MPI_File_write_ordered (
        out, checkpoint, n_local_cp * sizeof(particle_t), MPI_BYTE,
        MPI_STATUS_IGNORE
    );
    */
    MPI_File_close ( &out );
    MPI_Barrier ( MPI_COMM_WORLD );
}


void
dump_state ( char *filename )
{
    MPI_Barrier ( MPI_COMM_WORLD );
    int_t my_particles = n_particles();

    real_t data[my_particles][3];
    particle_t *actual_ptr[my_particles];
    list_particles ( &(actual_ptr[0]) );

    for ( int_t mp=0; mp<my_particles; mp++ )
    {
        particle_t *p = actual_ptr[mp];
        data[mp][0] = (real_t)p->idx;
        data[mp][1] = p->x[0];
        data[mp][2] = p->x[1];
    }

    MPI_Info info;
    MPI_Info_create ( &info );
    MPI_Info_set ( info, "access_style", "write_once" );
    MPI_File out;
    MPI_File_open (
        MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,
        info, &out
    );
    MPI_File_write_ordered (
        out, data, 3*my_particles*sizeof(real_t), MPI_BYTE, MPI_STATUS_IGNORE
    );
    MPI_File_close ( &out );
    MPI_Barrier ( MPI_COMM_WORLD );
}
#endif


void
restart_checkpoint ( int_t iteration )
{
    int_t file_number = iteration / checkpoint_frequency;

    /* Adjust starting point if requested iter. is not a multiple of freq.
     * Skip the iteration represented by the loaded checkpoint
     */
    min_iteration = file_number * checkpoint_frequency + 1;

    /* Generate the assumed checkpoint file name */
    char filename[256];
    memset ( filename, 0, 256*sizeof(char) );
    sprintf ( filename, "plot/%.4ld.dat", file_number );
    printf (
        "Rank %d resumes from '%s', iteration %ld, freq. %ld, maxiter %ld\n",
        rank, filename, min_iteration, checkpoint_frequency, max_iteration
    );

    /* Calibrate maxiter if it hasn't been set properly */
    if ( max_iteration <= min_iteration )
    {
        if ( rank == 0 )
            fprintf ( stderr,
                "Error: (min,max) iterations set to %ld, %ld, aborting.\n",
                min_iteration, max_iteration
            );
        MPI_Barrier ( MPI_COMM_WORLD );
        /* Stop with errno error code for 'invalid argument' */
        MPI_Abort ( MPI_COMM_WORLD, EINVAL );
    }

    /* The following setup is identical to that of initialize(): */

    /* Initialize hash table for field particles */
    particles_init();

    /* Calculate bounds of local subdomain */
    real_t subdomain_size = B / (real_t)size;
    subdomain[0] = rank * subdomain_size;
    subdomain[1] = (rank+1) * subdomain_size;

    if ( subdomain_size <= sqrt((scale_k*H)*(scale_k*H)) )
        fprintf ( stderr, "Rank %d: "
            "Warning, subdomain size smaller than interaction radius at "
            "scale %lf, results will not be correct.\n", rank, SCALE
        );

    /* Population of local subdomains differs from initialize(),
     * read the particle states from file instead
     */
    
    /* Read-only parallel open should not thrash parallel FS */
    FILE *checkpoint = fopen ( filename, "r" );

    /* Fail if the checkpoint file doesn't exist */
    if ( checkpoint == NULL )
    {
        fprintf ( stderr,
            "Error: Rank %d unable to open '%s', aborting\n", rank, filename
        );
        /* Stop with errno code for 'No such file or directory' */
        MPI_Abort ( MPI_COMM_WORLD, ENOENT );
    }

    /* Scan all the particles in the file, copy those within my subdomain */
    particle_t p;
    int items = fread ( &p, sizeof(particle_t), 1, checkpoint );
    n_global_field = items;
    while ( items != 0 )
    {
        /* This requires a hack so as not to drop particles temporarily
         * outside the domain (due to numerical inaccuracy)
         */
        if (
            (p.x[0] >= subdomain[0] && p.x[0] < subdomain[1])
         || (p.x[0] < 0.0 && rank == 0 )
         || (p.x[0] > B && rank == size-1)
        )
        {
            particle_t *keep_local = malloc ( sizeof(particle_t) );
            memcpy ( keep_local, &p, sizeof(particle_t) );
            insert_particle ( keep_local );
        }
        items = fread ( &p, sizeof(particle_t), 1, checkpoint );
        n_global_field += items;
    } 
    fclose ( checkpoint );

    /* Returning to mimic initialize() */
    n_field = n_particles();

    /* Start-allocation for list of local particles and pairs */
    list = (particle_t *) malloc ( n_capacity * sizeof(particle_t) );
    pairs = malloc ( n_pair_cap * sizeof(pair_t) );
}
