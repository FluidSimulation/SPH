#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <tlhash.h>

#include "sph.h"

static tlhash_t *particles;

void
particles_init ( void )
{
    particles = malloc ( sizeof(tlhash_t) );
    tlhash_init ( particles, 4096 );
}


void
particles_finalize ( void )
{
    tlhash_finalize ( particles );
    free ( particles );
}


void
insert_particle ( particle_t *p )
{
    tlhash_insert ( particles, &(p->idx), sizeof(int_t), (void *)p );
}


void
lookup_particle ( int_t index, particle_t **p )
{
    tlhash_lookup ( particles, &index, sizeof(int_t), (void **)p );
}


void
remove_particle ( particle_t *p )
{
    tlhash_remove ( particles, &(p->idx), sizeof(int_t) );
}


void
list_particles ( particle_t **list_point )
{
    tlhash_values ( particles, (void **)list_point );
}


/* NB - stack allocation of particle pointer list */
void
marshal_particles ( particle_t *list_point )
{
    int_t my_particles = n_particles();
    particle_t *actual_ptr[my_particles];
    list_particles ( &(actual_ptr[0]) );
    for ( int_t i=0; i<my_particles; i++ )
        memcpy ( &(list_point[i]), actual_ptr[i], sizeof(particle_t) );
}

void
unmarshal_particles ( particle_t *list, int_t count )
{
    for ( int_t k=0; k<count; k++ )
    {
        particle_t *listed = &(list[k]), *tabbed;
        lookup_particle ( listed->idx, &tabbed );
        memcpy ( tabbed, listed, sizeof(particle_t) );
    }
}

int_t
n_particles ( void )
{
    return (int_t)tlhash_size ( particles );
}
