#include "bkem.h"
#include <time.h>
clock_t t,t1,t0,t2,t4; 
     
int main(int argc, const char *argv[]) {

	FILE *param = fopen("e.param", "r");
	char buf[4096];
	fread(buf, 1, 4096, param);
	bkem_global_params_t gps;
	setup_global_system(&gps, (const char*) buf, (argc > 1) ? atoi(argv[1]) : 1);
	bkem_system_t sys;
	setup(&sys, gps); 

    }
    
