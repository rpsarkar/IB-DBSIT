 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "bkem.h"
#include <time.h>

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
    bkem_system_t gbs;
    gbs = pbc_malloc(sizeof(struct bkem_system_s));
    gbs->PK = pbc_malloc(sizeof(struct pubkey_s));

    // ---------------------------------Choose random generator P --------------------------------------------
    element_init_G1(gbs->PK->g, gps->pairing);
    element_random(gbs->PK->g);
    
    element_init_G2(gbs->PK->gg, gps->pairing);
    element_random(gbs->PK->gg);

    
    //==========================================  (P_2) ========================================================
    //=========================================================================================================
    element_t pairing_1;
    element_init_GT(pairing_1, gps->pairing);
    t4 = clock();
    pairing_apply(pairing_1,gbs->PK->g,gbs->PK->gg,gps->pairing);
    t4 = clock() - t4;
    double time_taken4 = ((double)t4*1000)/(CLOCKS_PER_SEC); // in seconds 
    printf(" Cost to compute T3 asymmtric pairing operation (P_2) = %f in milliseconds \n", time_taken4); 
    //=========================================================================================================
    //=========================================================================================================
    
    
    //==========================================  (AE_1) ========================================================
    //=========================================================================================================
    element_t e1,e2;
    element_init_G1(e1, gps->pairing);
    element_init_Zr(e2, gps->pairing);
    t5 = clock();
    element_pow_zn(e1,gbs->PK->g,e2);
    t5 = clock() - t5;
    double time_taken5 = ((double)t5*1000)/(CLOCKS_PER_SEC); // in seconds 
    printf(" Cost to compute exponentiation in 1st source group of T3 asymmtric pairing group (AE_1) = %f in milliseconds \n", time_taken5); 
    //=========================================================================================================
    //=========================================================================================================
    
    
    //==========================================  (AE_2) ========================================================
    //=========================================================================================================
    element_t e3,e4;
    element_init_G2(e3, gps->pairing);
    element_init_Zr(e4, gps->pairing);
    t6 = clock();
    element_pow_zn(e3,gbs->PK->gg,e4);
    t6 = clock() - t6;
    double time_taken6 = ((double)t6*1000)/(CLOCKS_PER_SEC); // in seconds 
    printf(" Cost to compute exponentiation in 2nd source group of T3 asymmtric pairing group (AE_2) = %f in milliseconds \n", time_taken6); 
    //=========================================================================================================
    //=========================================================================================================
    
    
    //==========================================  (AM_2) ========================================================
    //=========================================================================================================
    element_t e5,e6,e7;
    element_init_G2(e5, gps->pairing);
    element_init_G2(e6, gps->pairing);
    element_init_G2(e7, gps->pairing);
    t7 = clock();
    element_mul(e7,e5,e6);
    t7 = clock() - t7;
    double time_taken7 = ((double)t7*1000)/(CLOCKS_PER_SEC); // in seconds 
    printf(" Cost to compute a multiplication of elements of 2nd source group of T3 asymmtric pairing group (AM_2) = %f in milliseconds \n", time_taken7); 
    //=========================================================================================================
    //=========================================================================================================
    
    
    //==========================================  (Z_p) ========================================================
    //=========================================================================================================
    element_t e8;
    element_init_Zr(e8, gps->pairing);
    t8 = clock();
    element_random(e8);
    t8 = clock() - t8;
    long double time_taken8 = (( long double)t8*1000)/(CLOCKS_PER_SEC); // in seconds 
    printf(" Cost to compute a element from Z_p (|Z_p|) = %Lf in milliseconds \n", time_taken8); 
    //=========================================================================================================
    //=========================================================================================================
    
    
    //==========================================  SM_1  =======================================================
    //=========================================================================================================
    //In the addition, scaler multiplication is equivalent to exponentiation in multiplication group
    printf(" Cost to compute  a scaler multiplication in T3 source group (SM_1) = %f in miliseconds \n", time_taken5); 
    //=========================================================================================================
    //=========================================================================================================
    
    //==========================================  (PA_1) ========================================================
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
    printf(" Cost to compute a point addition in T3 source group (PA_1) = %Lf in miliseconds \n", time_taken9); 
    //=========================================================================================================
    //=========================================================================================================
    
    //==========================================  (SM_2) ========================================================
    //=========================================================================================================   
    element_t e12,e13,e14;
    element_init_G2(e12, gps->pairing);
    element_init_Zr(e13, gps->pairing);
    element_init_G2(e14, gps->pairing);
    element_random(e13);
    element_random(e14);
    t10 = clock();
    element_pow_zn(e14,e12,e13);
    t10 = clock() - t10;
    long double time_taken10 = (( long double)t10*1000)/(CLOCKS_PER_SEC); // in miliseconds 
    printf(" Cost to compute a scalar multiplication in the second T3 source group (SM_2) = %Lf in miliseconds \n", time_taken10); 
    //=========================================================================================================
    //=========================================================================================================
    
    //==========================================  (PA_2) ========================================================
    //=========================================================================================================   
    element_t e15,e16,e17;
    element_init_G2(e15, gps->pairing);
    element_init_G2(e16, gps->pairing);
    element_init_G2(e17, gps->pairing);
    element_random(e15);
    element_random(e16);
    t11 = clock();
    element_add(e17,e15,e16);
    t11 = clock() - t11;
    long double time_taken11 = (( long double)t11*1000)/(CLOCKS_PER_SEC); // in miliseconds 
    printf(" Cost to compute a scalar multiplication in the second T3 source group (PA_2) = %Lf in miliseconds \n", time_taken11); 
    //=========================================================================================================
    //=========================================================================================================
    
    //==========================================  (SM_T) ========================================================
    //=========================================================================================================   
    element_t e121,e131,e141;
    element_init_GT(e121, gps->pairing);
    element_init_Zr(e131, gps->pairing);
    element_init_GT(e141, gps->pairing);
    element_random(e131);
    element_random(e141);
    t12 = clock();
    element_pow_zn(e141,e121,e131);
    t12 = clock() - t12;
    long double time_taken12 = (( long double)t12*1000)/(CLOCKS_PER_SEC); // in miliseconds 
    printf(" Cost to compute a a scalar multiplication in the T3 target group (SM_T) = %Lf in miliseconds \n", time_taken12); 
    //=========================================================================================================
    //=========================================================================================================
    
    //==========================================  (PA_T) ========================================================
    //=========================================================================================================   
    element_t e151,e161,e171;
    element_init_GT(e151, gps->pairing);
    element_init_GT(e161, gps->pairing);
    element_init_GT(e171, gps->pairing);
    element_random(e151);
    element_random(e161);
    t13 = clock();
    element_add(e171,e151,e161);
    t13 = clock() - t13;
    long double time_taken13 = (( long double)t13*1000)/(CLOCKS_PER_SEC); // in miliseconds 
    printf(" Cost to compute a point addition in the T3 target group (PA_T) = %Lf in miliseconds \n", time_taken13); 
    //=========================================================================================================
    //=========================================================================================================
    
    
  
    	
     *sys = gbs;
    
   
 }


