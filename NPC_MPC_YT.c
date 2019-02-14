/*
 *
 */


#define S_FUNCTION_NAME  NPC_MPC_YT 
#define S_FUNCTION_LEVEL 2
#define _USE_MATH_DEFINES

// #define Tsampling 25e-6



#include "simstruc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <simstruc.h>

#define MYFILEPR my_file_pointer_file   /*Define variable representing actually file pointer name */
FILE *MYFILEPR;                         /*Declare file pointer to a file */


/*Variables globales*/


float is_aB[2];
float is_dq[2];

float phir_dq[2];

float Tsampling = 100e-6;

float isJ_aB[2];

float ie_aB[2];

float is_abc[3];

float uabc[3];
int a,b,c;

float phir_aB[2];
// float phir_abc[3];

float wr;

float A[2][2];

float G2;

float G1;

float B1[2][2];

float B2[2][3];

float J;

float lambda;

float incr_u;

float Rs;
float Rr;
float Xls;
float Xlr;
float Xm;
float Xs;
float Xr;
float D;
float tau_s;
float tau_r;

float teta;

float TORQUE;

float TaB_dq[2][2];
float Tabc_aB[2][3];

float flux_r;
float refflux;

float isref_dq[2];

float isref_aB[2];

float histerrorflujo[2];

float histflujoref[2];

float histerrorvel[2];

float histTe[2];

float histteta[2];
float histwr[2];
float histwm[2];

float histphir_a[2];
float histphir_B[2];

float Kp_flux=10.;//*1e-1;

float Ki_flux=20.;//*1e-1;

float Kp_w=16.6*1e3;  //16600.;
    
float Ki_w=27.7*1e3;   //27700.;



/*================*
 * Build checking *
 *================*/


/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *   Setup sizes of the various vectors.
 */
static void mdlInitializeSizes(SimStruct *S)
{
//     MYFILEPR=fopen("my_data_file2","w");  /* Opens a file, whose name is my_data_file, for writing */ 
    
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    if (!ssSetNumInputPorts(S, 3)) return; //aqui defino el numero de entradas
    //ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetInputPortWidth(S, 0, 3); // dimension de la entrada 
    ssSetInputPortDirectFeedThrough(S, 0, 1); //no la usamos, define el comportamiento de una entrada
    ssSetInputPortWidth(S, 1, 2); //
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortWidth(S, 2, 1); //
    ssSetInputPortDirectFeedThrough(S, 2, 1);
//     ssSetInputPortWidth(S, 3, 1); //
//     ssSetInputPortDirectFeedThrough(S, 3, 1);


    if (!ssSetNumOutputPorts(S,11)) return; //
    //ssSetOutputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortWidth(S, 2, 1);
    ssSetOutputPortWidth(S, 3, 1);
    ssSetOutputPortWidth(S, 4, 1);
    ssSetOutputPortWidth(S, 5, 1);
    ssSetOutputPortWidth(S, 6, 1);
    ssSetOutputPortWidth(S, 7, 1);
    ssSetOutputPortWidth(S, 8, 1);
    ssSetOutputPortWidth(S, 9, 1);
    ssSetOutputPortWidth(S, 10, 2);
    
    ssSetNumSampleTimes(S, 1);

    /* specify the sim state compliance to be same as a built-in block */
    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);

    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S,
                 SS_OPTION_WORKS_WITH_CODE_REUSE |
                 SS_OPTION_EXCEPTION_FREE_CODE |
                 SS_OPTION_USE_TLC_WITH_ACCELERATOR);
}


/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy that we inherit our sample time from the driving block.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    //ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetSampleTime(S, 0, Tsampling);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S); 
}

#define MDL_START                      /* Change to #undef to remove function */
#if defined(MDL_START)
/* Function: mdlStart ==========================================================
 * Abstract:
 *
 */
static void mdlStart(SimStruct *S) 
{
    Tabc_aB[0][0]=sqrt(2./3.);
    Tabc_aB[0][1]=-0.5*sqrt(2./3.);
    Tabc_aB[0][2]=-0.5*sqrt(2./3.);
    Tabc_aB[1][0]=0.;
    Tabc_aB[1][1]=sqrt(2./3.)*sqrt(3.)/2.;
    Tabc_aB[1][2]=-sqrt(2./3.)*sqrt(3.)/2.;
    
//     Rs=.6837;
//     Rr=.451;
//     Xls=.004125;
//     Xlr=.004125;
//     Xm=0.1486;
//     Xs=Xls+Xm;
//     Xr=Xlr+Xm;
//     D=Xs*Xr-Xm*Xm;
//     tau_s=Xr*D/(Rs*Xr*Xr+Rr*Xm*Xm);
//     tau_r=Xr/Rr;
// //     T_r=Xlr/Rr;
    
//     Rs=.0108;
//     Rr=.0091;
//     Xls=.1493;
//     Xlr=.1104;
//     Xm=2.394;
    Rs=.05761;
    Rr=.04889;
    Xls=.002544;
    Xlr=.001881;
    Xm=.04001;
    Xs=Xls+Xm;
    Xr=Xlr+Xm;
    D=Xs*Xr-Xm*Xm;
    tau_s=Xr*D/(Rs*Xr*Xr+Rr*Xm*Xm);
    tau_r=Xr/Rr;
    
    A[0][0]=1.-Tsampling/tau_s;
    A[0][1]=0;        
    A[1][0]=0;        
    A[1][1]=1.-Tsampling/tau_s;
    
    G1=Tsampling*Xm/D;
    
    G2=Tsampling*5200.*Xr/(3.*D);
    B2[0][0]=G2;
    B2[0][1]=-G2/2.;
    B2[0][2]=-G2/2.;
    B2[1][0]=0;
    B2[1][1]=sqrt(3.)*G2/2.;
    B2[1][2]=-sqrt(3.)*G2/2.;
    
    uabc[0]=0.;
    uabc[1]=0.;
    uabc[2]=0.;
    
    lambda=0.003; 
    
    histerrorflujo[0]=0.;
    histerrorflujo[1]=0.;
    histflujoref[0]=0.;
    histflujoref[1]=0.;
    histerrorvel[0]=0.;
    histerrorvel[1]=0.;
    histTe[0]=0.;
    histTe[1]=0.;
    histteta[0]=0.;
    histteta[1]=0.;
    histwr[0]=0.;
    histwr[1]=0.;
    histwm[0]=0.;
    histwm[1]=0.;
    histphir_a[0]=0.;
    histphir_a[1]=0.;
    histphir_B[0]=0.;
    histphir_B[1]=0.;

}
#endif /*  MDL_START */


/* Function: mdlOutputs =======================================================
 * Abstract:
 *    
 */
static void mdlOutputs(SimStruct *S, int_T tid) //genera una salida cada vez q se llama al bloque
                                                //aqui es donde va el algoritmo que yo quiera ejecutar
{
    //int_T             i;
    InputRealPtrsType pis_abc= ssGetInputPortRealSignalPtrs(S,0);
    
    InputRealPtrsType pphir_aB = ssGetInputPortRealSignalPtrs(S,1);
    
    InputRealPtrsType pwr = ssGetInputPortRealSignalPtrs(S,2);
    
//     InputRealPtrsType pteta = ssGetInputPortRealSignalPtrs(S,3);
    
    real_T            *S11    = ssGetOutputPortRealSignal(S,0);
    
    real_T            *S12    = ssGetOutputPortRealSignal(S,1);
    
    real_T            *S21    = ssGetOutputPortRealSignal(S,2);
    
    real_T            *S22    = ssGetOutputPortRealSignal(S,3);
    
    real_T            *S31    = ssGetOutputPortRealSignal(S,4);
    
    real_T            *S32    = ssGetOutputPortRealSignal(S,5);
    
    real_T            *ua    = ssGetOutputPortRealSignal(S,6);
    
    real_T            *ub    = ssGetOutputPortRealSignal(S,7);
    
    real_T            *uc    = ssGetOutputPortRealSignal(S,8);
    
    real_T            *pJ    = ssGetOutputPortRealSignal(S,9);
    
    real_T            *pphir_dq    = ssGetOutputPortRealSignal(S,10);
    
   
    // Vector de transiciones
    
    int i;
    
   
    float Uk[27][4]={{-1.,-1.,-1.,true},{-1.,-1.,0.,true},{-1.,-1.,1.,true},{-1.,0.,-1.,true},{-1.,0.,0.,true},{-1.,0.,1.,true},
                     {-1.,1.,-1.,true},{-1.,1.,0.,true},{-1.,1.,1.,true},{0.,-1.,-1.,true},{0.,-1.,0.,true},{0.,-1.,1.,true},
                     {0.,0.,-1.,true},{0.,0.,0.,true},{0.,0.,1.,true},{0.,1.,-1.,true},{0.,1.,0.,true},{0.,1.,1.,true},
                     {1.,-1.,-1.,true},{1.,-1.,0.,true},{1.,-1.,1.,true},{1.,0.,-1.,true},{1.,0.,0.,true},{1.,0.,1.,true},
                     {1.,1.,-1.,true},{1.,1.,0.,true},{1.,1.,1.,true}};
   
    float Jmin=INFINITY;
    
    //Medidas
    
    is_abc[0] = *pis_abc[0];
    is_abc[1] = *pis_abc[1];
    is_abc[2] = *pis_abc[2];
    
    phir_aB[0] = *pphir_aB[0];
    phir_aB[1] = *pphir_aB[1];
    
    histwm[0] = *pwr[0];
    
//     teta = *pteta[0];
//     
//     teta = teta + 2*M_PI*50*Tsampling;
    
    // Recalculo modelo
    
   
    // Norma flujo rotor
//     
    flux_r=sqrt(phir_aB[0]*phir_aB[0]+phir_aB[1]*phir_aB[1]);
//     flux_r=sqrt(histphir_B[0]*histphir_B[0]+histphir_a[0]*histphir_a[0]);
// //     
    if(flux_r<0.001){flux_r=0.001;}
    
    //  ABC-DQ
    
    is_dq[0]=2./3.*(is_abc[0]*cos(histteta[1])+is_abc[1]*cos(histteta[1]-2.*M_PI/3.)+is_abc[2]*cos(histteta[1]+2.*M_PI/3.));
    is_dq[1]=-2./3.*(is_abc[0]*sin(histteta[1])+is_abc[1]*sin(histteta[1]-2.*M_PI/3.)+is_abc[2]*sin(histteta[1]+2.*M_PI/3.));
    
//  CALCULO TETA
    
    histwr[0] = Xm*is_dq[1]/(flux_r*tau_r);
    
    B1[0][0]=G1/tau_r;
    B1[0][1]=G1*histwr[0];
    B1[1][0]=-G1*histwr[0];
    B1[1][1]=G1/tau_r;
    
    
    histteta[1]=Tsampling*(histwr[0]+5.*histwm[0])+histteta[0];
    histteta[0]=histteta[1];
    
//  CALCULO IDS* con control
//     refflux=1.9;
//     histerrorflujo[0]=refflux-flux_r; // 
//     
//     histflujoref[0]=Kp_flux*histerrorflujo[0]+(Ki_flux*Tsampling-Kp_flux)*histerrorflujo[1]+histflujoref[1]; 
//     histflujoref[1]=histflujoref[0];
//     histerrorflujo[1]=histerrorflujo[0];
//     
//     isref_dq[0]=histflujoref[0]+refflux/Xm;
// //     
    
    isref_dq[0]=0.9/Xm;

    
// Controlador velocidad
    
    histerrorvel[0]=55.-histwm[0]; // REFERENCIA NOMINAL 62.413
    
    histTe[0]=Kp_w*histerrorvel[0]+(Ki_w*Tsampling-Kp_w)*histerrorvel[1]+histTe[1]; 
    histTe[1]=histTe[0];
    histerrorvel[1]=histerrorvel[0];
    
    
    
//  CALCULO IQS*    
    isref_dq[1]=2./3./5.*Xr/Xm/flux_r*histTe[0];
    

    
    TaB_dq[0][0]=cos(histteta[1]);
    TaB_dq[0][1]=sin(histteta[1]);  
    TaB_dq[1][0]=-sin(histteta[1]);
    TaB_dq[1][1]=cos(histteta[1]);

    
    isref_aB[0]=TaB_dq[0][0]*isref_dq[0]+TaB_dq[1][0]*isref_dq[1];
    isref_aB[1]=TaB_dq[0][1]*isref_dq[0]+TaB_dq[1][1]*isref_dq[1];
    
    is_aB[0]=Tabc_aB[0][0]*is_abc[0]+Tabc_aB[0][1]*is_abc[1]+Tabc_aB[0][2]*is_abc[2];
    is_aB[1]=Tabc_aB[1][0]*is_abc[0]+Tabc_aB[1][1]*is_abc[1]+Tabc_aB[1][2]*is_abc[2];
    
    histphir_a[1]=histphir_a[0]+Tsampling*(Xm/tau_r*is_aB[0]-histphir_a[0]/tau_r-histphir_B[0]*histwr[0]);
    histphir_B[1]=histphir_B[0]+Tsampling*(Xm/tau_r*is_aB[1]-histphir_B[0]/tau_r+histphir_a[0]*histwr[0]);
    histphir_a[0]=histphir_a[1];
    histphir_B[0]=histphir_B[1];
    
    //// MPC
    
    // Conjunto de transiciones posibles
    
    if(uabc[0]==-1)
    {
        for(i=19;i<=27;i++)
        {
            Uk[i-1][3]=false;
        }
    }
    if(uabc[0]==1)
    {
        for(i=1;i<=9;i++)
        {
            Uk[i-1][3]=false;
        }
    }
    if(uabc[1]==-1)
    {
        for(i=7;i<=9;i++)
        {
            Uk[i-1][3]=false;
        }
        for(i=16;i<=18;i++)
        {
            Uk[i-1][3]=false;
        }
        for(i=25;i<=27;i++)
        {
            Uk[i-1][3]=false;
        }  
    }
    if(uabc[1]==1)
    {
        for(i=1;i<=3;i++)
        {
            Uk[i-1][3]=false;
        }
        for(i=10;i<=12;i++)
        {
            Uk[i-1][3]=false;
        }
        for(i=19;i<=21;i++)
        {
            Uk[i-1][3]=false;
        }  
    }
    if(uabc[2]==-1)
    {
        for(i=3;i<=27;i=i+3)
        {
            Uk[i-1][3]=false;
        }
    }
    if(uabc[2]==1)
    {
        for(i=1;i<=27;i=i+3)
        {
            Uk[i-1][3]=false;
        }
    }
    
    // Calculo funcion de coste
    
//     fprintf(MYFILEPR,"estado actual: %f %f %f \n", uabc[0],uabc[1],uabc[2]);
    
//     for(i=1;i<27;i++)
//     {
//         if(Uk[i-1][3] == true)
//         {
//                 fprintf(MYFILEPR,"posible: %f %f %f \n", Uk[i-1][0],Uk[i-1][1],Uk[i-1][2]);
//         }
//     }
    
    for(i=1;i<27;i++)
    {
        if(Uk[i-1][3] == true)
        {
            isJ_aB[0]=A[0][0]*is_aB[0]+
                     B1[0][0]*phir_aB[0]+B1[0][1]*phir_aB[1]+
                     B2[0][0]*Uk[i-1][0]+B2[0][1]*Uk[i-1][1]+B2[0][2]*Uk[i-1][2];
            isJ_aB[1]=A[1][1]*is_aB[1]+
                     B1[1][0]*phir_aB[0]+B1[1][1]*phir_aB[1]+
                     B2[1][0]*Uk[i-1][0]+B2[1][1]*Uk[i-1][1]+B2[1][2]*Uk[i-1][2];
//             isJ_aB[0]=A[0][0]*is_aB[0]+
//                      B1[0][0]*histphir_a[0]+B1[0][1]*histphir_B[0]+
//                      B2[0][0]*Uk[i-1][0]+B2[0][1]*Uk[i-1][1]+B2[0][2]*Uk[i-1][2];
//             isJ_aB[1]=A[1][1]*is_aB[1]+
//                      B1[1][0]*histphir_a[0]+B1[1][1]*histphir_B[0]+
//                      B2[1][0]*Uk[i-1][0]+B2[1][1]*Uk[i-1][1]+B2[1][2]*Uk[i-1][2];

            ie_aB[0]=isref_aB[0]-isJ_aB[0];
            ie_aB[1]=isref_aB[1]-isJ_aB[1];

    //         ie_abc[0]=Tabc_aB[0][0]*ie_aB[0]+Tabc_aB[1][0]*ie_aB[1];
    //         ie_abc[1]=Tabc_aB[0][1]*ie_aB[0]+Tabc_aB[1][1]*ie_aB[1];
    //         ie_abc[2]=Tabc_aB[0][2]*ie_aB[0]+Tabc_aB[1][2]*ie_aB[1];

            incr_u=fabsf(Uk[i-1][0]-uabc[0])+
                   fabsf(Uk[i-1][1]-uabc[1])+
                   fabsf(Uk[i-1][2]-uabc[2]);

            J=ie_aB[0]*ie_aB[0]+ie_aB[1]*ie_aB[1]+lambda*incr_u;

            if(J<Jmin)
            {
                Jmin = J;
                a = Uk[i-1][0];
                b = Uk[i-1][1];
                c = Uk[i-1][2];
//                 fprintf(MYFILEPR,"cost: %f \n", Jmin);
            }
        }
    }
    
   uabc[0]   = a;
   uabc[1]   = b;
   uabc[2]   = c; 
   TORQUE=Xm/0.78/Xr*(flux_r*is_dq[1]);
//    TORQUE=Xm/0.78/Xr*(histphir_a[0]*is_aB[1]-histphir_B[0]*is_aB[0]);
   
   ua[0]   = a;
   ub[0]   = b;
   uc[0]   = c; 
   pJ[0]   = TORQUE;
   pphir_dq[0] = flux_r;
   pphir_dq[1] = histphir_B[1];

           
//    if(uabc[0] == 1){S11[0]=1;S12[0]=1;}
//    if(uabc[0] == 0){S11[0]=0;S12[0]=1;}
//    if(uabc[0] == -1){S11[0]=0;S12[0]=0;}
//    if(uabc[1] == 1){S21[0]=1;S22[0]=1;}
//    if(uabc[1] == 0){S21[0]=0;S22[0]=1;}
//    if(uabc[1] == -1){S21[0]=0;S22[0]=0;}
//    if(uabc[2] == 1){S31[0]=1;S32[0]=1;}
//    if(uabc[2] == 0){S31[0]=0;S32[0]=1;}
//    if(uabc[2] == -1){S31[0]=0;S32[0]=0;}
       
   switch(a) {

       case 1 :
          S11[0]=1;
          S12[0]=1;
          break; 

       case 0  :
          S11[0]=0;
          S12[0]=1;
          break; 

       case -1 :
          S11[0]=0;
          S12[0]=0;
          break; 

       default : 
          S11[0]=100;
          S12[0]=100;
   }  
   switch(b) {

       case 1 :
          S21[0]=1;
          S22[0]=1;
          break; 

       case 0  :
          S21[0]=0;
          S22[0]=1;
          break; 

       case -1 :
          S21[0]=0;
          S22[0]=0;
          break; 

       default : 
          S21[0]=100;
          S22[0]=100;
   
   } 
   switch(c) {
    
       case 1 :
          S31[0]=1;
          S32[0]=1;
          break; 

       case 0  :
          S31[0]=0;
          S32[0]=1;
          break; 

       case -1 :
          S31[0]=0;
          S32[0]=0;
          break; 

       default : 
          S31[0]=100;
          S32[0]=100;
   } 
   
 
    
}


/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S) //la ultima funciona que se ejecuta antes de terminar la simulacion
{
//     fclose(MYFILEPR); /*Closes the file, my_data_file */
}



#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
