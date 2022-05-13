/**
 * @file BKEM.h
 * @brief General construction of the Boneh-Gentry-Waters 
 * broadcast key encapsulation scheme 
 *
 * BKEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * BKEM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with BKEM.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Oliver Guenther
 * mail@oliverguenther.de
 *
 * 
 * BKEM.h
*/

#ifndef H_BKEM
#define H_BKEM

#include <string.h>
#include <pbc/pbc.h>

 
#define MAX_N  1          
#define MAX_SET MAX_n         // Here, Sets are  Users in Channel Sj

 
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
