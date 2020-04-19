#ifndef TLHASH_H
#define TLHASH_H
#include <stddef.h>
typedef struct el {
    void *key, *value;
    size_t key_length;
    struct el *next;
} tlhash_element_t;

typedef struct {
    size_t n_buckets, size;
    tlhash_element_t **buckets;
} tlhash_t;

int tlhash_init ( tlhash_t *tab, size_t n_buckets );
int tlhash_finalize ( tlhash_t *tab );
int tlhash_insert ( tlhash_t *tab, void *key, size_t keylen, void *val );
int tlhash_lookup ( tlhash_t *tab, void *key, size_t keylen, void **val );
int tlhash_remove ( tlhash_t *tab, void *key, size_t key_length );
size_t tlhash_size ( tlhash_t *tab );
void tlhash_keys ( tlhash_t *tab, void **keys );
void tlhash_values ( tlhash_t *tab, void **values );

#define TLHASH_SUCCESS 0    /* Success */
#define TLHASH_ENOMEM 1     /* No memory available */
#define TLHASH_ENOENT 2     /* No such table entry */
#define TLHASH_EEXIST 3     /* Table entry already exists */
#endif
