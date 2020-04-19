#ifdef __APPLE__
#define _GNU_C_SOURCE
#elif __linux__
#define _XOPEN_SOURCE 500
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <getopt.h>
#include <unistd.h>
#include <math.h>

#define MIN_ITERATION_DEFAULT 0
#define MAX_ITERATION_DEFAULT 200000
#define CHECKPOINT_FREQUENCY_DEFAULT 200

/* Problem parameters */

// Scale coefficient selects problem size
#ifndef SCALE
    #define SCALE (1.0)
#endif //SCALE

// Height of dam
#define T (0.6 * SCALE)

// Width of dam
#define L (1.2 * SCALE)

// Size of tank
#define B (3.22 * SCALE)

// Resolution
#define DELTA (0.01)
#define H (0.94*DELTA*1.4142135623)

typedef int64_t int_t;
#define INT_MACRO_MPI MPI_LONG
typedef double real_t;
#define REAL_MACRO_MPI MPI_DOUBLE

#define BID(X,Y) (Y + N_BUCKETS_Y * X)

#define RADIUS (scale_k * H)
#define BUCKET_RADIUS (1*RADIUS)

#define N_BUCKETS_X ((int_t)(ceil(((subdomain[1]-subdomain[0])+2*RADIUS) / BUCKET_RADIUS))) //RADIUS is the maximum distance (from each sides of the subdomain boundaries) where the mirror particles can be located. RADIUS>1.55*H (otherwise replace RADIUS by 1.55*H when computing N_BUCKETS_X)
#define N_BUCKETS_Y ((int_t)(ceil(((1.5*T)+1.55*H) / BUCKET_RADIUS))) //1.55*H is the boundary used when generating virtual particles


#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))

static const real_t
    D_tank = 1.0 * SCALE,
    scale_k = 3.0,
    density = 1e3,
    dt = 1e-4,
    sos = 50.0;

static const int_t
    free_surface = 30;

typedef struct {
    int_t
        idx,
        interactions;
    real_t
        x[2],   // 2D position
        v[2],   // 2D velocity
        mass,   // Mass
        rho,    // Density
        p,      // Pressure
        type,   // Type of particle
        hsml;   // Related to resolution (not fully understood yet)

    // Differentials
    real_t
        indvxdt[2],
        exdvxdt[2],
        dvx[2],
        drhodt;
    // Density correction
    real_t
        avrho,
        w_sum;
#ifdef BUCKET
    int_t bucket_x;
    int_t bucket_y;
    int_t local_idx;
#endif //BUCKET
} particle_t;

#define X(k)        ((list[k]).x[0])
#define Y(k)        ((list[k]).x[1])
#define VX(k)       ((list[(k)]).v[0])
#define VY(k)       ((list[(k)]).v[1])
#define M(k)        ((list[(k)]).mass)
#define RHO(k)      ((list[(k)]).rho)
#define P(k)        ((list[(k)]).p)
#define TYPE(k)     ((list[(k)]).type)
#define HSML(k)     ((list[(k)]).hsml)
#define INTER(k)    ((list[(k)]).interactions)

#define INDVXDT(k,i)    ((list[(k)]).indvxdt[(i)])
#define EXDVXDT(k,i)    ((list[(k)]).exdvxdt[(i)])
#define DVX(k,i)        ((list[(k)]).dvx[(i)])
#define DRHODT(k)       ((list[(k)]).drhodt)

#define AVRHO(k)    ((list[(k)]).avrho)
#define WSUM(k)     ((list[(k)]).w_sum)

// Pairwise interaction
typedef struct {
    int_t i, j;     // Which particles interact?
#ifdef BUCKET
    particle_t* ip;
    particle_t* jp;
#endif //BUCKET
    real_t
        r,          // Distance between particles (Euclid)
        q,          // Distance normalized to H (resolution)
        w,          // TODO: consult paper for meaning of this
        dwdx[2];    // Influence on velocity
} pair_t;

typedef struct bucket_t bucket_t;
/* A bucket is a node in a linked list */
struct bucket_t{
    particle_t *particle;
    bucket_t *next;
};

/* Global state variables, definitions are in sph.c */
extern int size, rank, east, west;
extern int_t n_global_field, n_field;

/* Setup and takedown */
void initialize ( void );
void finalize ( void );

// Hash table interface (in particle_hashtab.c)
void particles_init ( void );
void particles_finalize ( void );

// Insert a pointer to a particle
void insert_particle ( particle_t *p );
// Remove a pointer to a particle
void remove_particle ( particle_t *p );
// Lookup particle pointer by index
void lookup_particle ( int_t index, particle_t **p );
// Start from particle pointer, serialize contents of table (pointers)
void list_particles ( particle_t **list_point );
// Start from particle pointer, serialize particle data
void marshal_particles ( particle_t *list_point );
// Extract actual particles from list to hash tab
void unmarshal_particles ( particle_t *list, int_t count );
// Count the number of local actuals
int_t n_particles ( void );

// I/O and auxiliary stuff
void dump_state ( char *filename );
void resize_list ( int_t required );
void resize_pair_list ( int_t new_cap );
void collect_checkpoint ( void );
void write_checkpoint ( char *filename );
void restart_checkpoint ( int_t iteration );
void options ( int argc, char **argv );
void free_buckets(bucket_t* bucket);
void print_timing(char* full_string, char* short_string, double value);

// Parts of the solver
void generate_virtual_particles ( void );
void create_pairs(int bx, int by, bucket_t** buckets,
                  particle_t* particle, int_t* n_pairs,
                  int_t* interactions);
void create_pairs_old(int bx, int by, bucket_t* buckets,
                      particle_t* particle, int_t* n_pairs);
#ifdef FILL_BUCKETS_LOCK
static void fill_buckets2(bucket_t**, int, omp_lock_t* lock);
#else
static void fill_buckets1(bucket_t**, int);
#endif //FILL_BUCKETS_LOCK


// MPI communication
void border_exchange( void );
void migrate_particles ( void );
