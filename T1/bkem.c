 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "bkem.h"
#include <time.h>
int num_recip=1000;

clock_t t,t1,t0,t2,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14;

void setup_global_system(bkem_global_params_t *gps, const char *pstr, int N) {
    
    bkem_global_params_t params;
    params = pbc_malloc(sizeof(struct bkem_global_params_s));

    params->N = N;
    
    pairing_init_set_str(params->pairing, pstr);

    *gps = params;
}

void setup(bkem_system_t *sys, bkem_global_params_t gps) 
{
    
    
    // 
    bkem_system_t gbs;
    gbs = pbc_malloc(sizeof(struct bkem_system_s));
    gbs->PK = pbc_malloc(sizeof(struct pubkey_s));

    // ---------------------------------Choose random generator P --------------------------------------------
    element_init_G1(gbs->PK->g, gps->pairing);
    element_random(gbs->PK->g);
    
    
    //==========================================  (P_1) ========================================================
    //=========================================================================================================
    element_t pairing_1;
    element_init_GT(pairing_1, gps->pairing);
    t4 = clock();
    pairing_apply(pairing_1,gbs->PK->g,gbs->PK->g,gps->pairing);
    t4 = clock() - t4;
    double time_taken4 = ((double)t4*1000)/(CLOCKS_PER_SEC); // in miliseconds 
    printf(" Cost to compute T1 symmetric pairing operation (P_1) = %f in miliseconds \n", time_taken4); 
    //=========================================================================================================
    //=========================================================================================================
    
    
    //==========================================  (SE) ========================================================
    //=========================================================================================================
    element_t e1,e2;
    element_init_G1(e1, gps->pairing);
    element_init_Zr(e2, gps->pairing);
    t5 = clock();
    element_pow_zn(e1,gbs->PK->g,e2);
    t5 = clock() - t5;
    double time_taken5 = ((double)t5*1000)/(CLOCKS_PER_SEC); // in seconds 
    printf(" Cost to compute exponentiation in source group of T1 symmetric pairing group (SE) = %f in miliseconds \n", time_taken5); 
    //=========================================================================================================
    //=========================================================================================================
    
    
    //==========================================  (T) ========================================================
    //=========================================================================================================
    element_t e3,e4;
    element_init_G1(e3, gps->pairing);
    t6 = clock();
    element_random(e3);
    t6 = clock() - t6;
    double time_taken6 = ((double)t6*1000)/(CLOCKS_PER_SEC); // in miliseconds 
    printf(" Cost to compute a element from source group of T1 symmetric pairing group (T) = %f in miliseconds \n", time_taken6); 
    //=========================================================================================================
    //=========================================================================================================
    
    
    //==========================================  (M^s b ) ========================================================
    //=========================================================================================================
    element_t e5,e6,e7;
    element_init_G1(e5, gps->pairing);
    element_init_G1(e6, gps->pairing);
    element_init_G1(e7, gps->pairing);
    t7 = clock();
    element_mul(e7,e5,e6);
    t7 = clock() - t7;
    double time_taken7 = ((double)t7*1000)/(CLOCKS_PER_SEC); // in seconds 
    printf(" Cost to compute a multiplication in source group of T1 symmtric pairing group (M^s) = %f in miliseconds \n", time_taken7); 
    //=========================================================================================================
    //=========================================================================================================
    
    
    //==========================================  (Z_p) ========================================================
    //=========================================================================================================
    element_t e8;
    element_init_Zr(e8, gps->pairing);
    t8 = clock();
    element_random(e8);
    t8 = clock() - t8;
    long double time_taken8 = (( long double)t8*1000)/(CLOCKS_PER_SEC); // in miliseconds 
    printf(" Cost to compute  a element from Z_p (|Z_p|) = %Lf in miliseconds \n", time_taken8); 
    //=========================================================================================================
    //=========================================================================================================
    
    
    //==========================================  SM ========================================================
    //=========================================================================================================
    //In the addition, scaler multiplication is equivalent to exponentiation in multiplication group
    printf(" Cost to compute  a scaler multiplication in T1 source group (SM) = %f in miliseconds \n", time_taken5); 
    //=========================================================================================================
    //=========================================================================================================
    
    //==========================================  (PA) ========================================================
    //=========================================================================================================   
    element_t e9,e10,e11;
    element_init_G1(e9, gps->pairing);
    element_init_G1(e10, gps->pairing);
    element_init_G1(e11, gps->pairing);
    element_random(e9);
    element_random(e10);
    t9 = clock();
    element_add(e11,e9,e10);
    t9 = clock() - t9;
    long double time_taken9 = (( long double)t9*1000)/(CLOCKS_PER_SEC); // in miliseconds 
    printf(" Cost to compute a point addition in T1 source group (PA) = %Lf in miliseconds \n", time_taken9); 
    //=========================================================================================================
    //=========================================================================================================
    
    //==========================================  (SM*) ========================================================
    //=========================================================================================================   
    element_t e12,e13;
    element_init_GT(e12, gps->pairing);
    element_init_Zr(e13, gps->pairing);
    element_random(e12);
    element_random(e13);
    t10 = clock();
    element_pow_zn(e12,e12,e13);
    t10 = clock() - t10;
    long double time_taken10 = (( long double)t10*1000)/(CLOCKS_PER_SEC); // in miliseconds 
    printf(" Cost to compute a  scalar multiplication operation of elements of the T1 target group (SM*) = %Lf in miliseconds \n", time_taken10); 
    //=========================================================================================================
    //=========================================================================================================
  
  //==========================================  (PA*) ========================================================
    //=========================================================================================================   
    element_t e14,e15,e16;
    element_init_GT(e14, gps->pairing);
    element_init_GT(e15, gps->pairing);
    element_init_GT(e16, gps->pairing);
    element_random(e14);
    element_random(e15);
    t11 = clock();
    element_add(e14,e15,e16);
    t11 = clock() - t11;
    long double time_taken11 = (( long double)t11*1000)/(CLOCKS_PER_SEC); // in miliseconds 
    printf(" Cost to compute a point addition in T1 target group (PA*) = %Lf in miliseconds \n", time_taken11); 
    //=========================================================================================================
    //=========================================================================================================
  
  
    	
     *sys = gbs;
    
   
 }


