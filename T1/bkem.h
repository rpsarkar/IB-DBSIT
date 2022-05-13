
#ifndef H_BKEM
#define H_BKEM

#include <string.h>
#include <pbc/pbc.h>
#define MAX_N  1           
#define MAX_SET MAX_n
 
typedef struct bkem_global_params_s {
	pairing_t pairing;
	int N;
	
}* bkem_global_params_t;

typedef struct pubkey_s {
    element_t g; // generator
    element_t gg; // generator
}* pubkey_t;
	

typedef struct bkem_system_s {
	pubkey_t PK;
}* bkem_system_t;

void setup_global_system(bkem_global_params_t *gps, const char *params, int n);
void setup(bkem_system_t *sys, bkem_global_params_t gps);
#endif
