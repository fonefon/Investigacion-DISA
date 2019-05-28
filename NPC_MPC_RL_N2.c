#define S_FUNCTION_NAME  NPC_MPC_RL_N2 
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


float iabc[3],i_aB[2],iref[3],iref_aB[2],iref2[3],iref2_aB[2],iJ_aB[2],iJ2_aB[2],ie[2],ie2[2];
float vabc[3],v_aB[2];
float pref;
float Vph_rms;
float ampli;
float vdc,R,L;
float Tsampling = 100e-6;
float u_0[3],incr_u,u_act1,u_act2,u_act3;
float Uk2[27][4];
int a,b,c;
float lambda,delta,J;
// float seno[10000];

// float A,B;

float Tabc_aB[2][3];

int i,cont,k;
float seq[27][4]={{-1.,-1.,-1.,true},{-1.,-1.,0.,true},{-1.,-1.,1.,true},{-1.,0.,-1.,true},{-1.,0.,0.,true},{-1.,0.,1.,true},
                     {-1.,1.,-1.,true},{-1.,1.,0.,true},{-1.,1.,1.,true},{0.,-1.,-1.,true},{0.,-1.,0.,true},{0.,-1.,1.,true},
                     {0.,0.,-1.,true},{0.,0.,0.,true},{0.,0.,1.,true},{0.,1.,-1.,true},{0.,1.,0.,true},{0.,1.,1.,true},
                     {1.,-1.,-1.,true},{1.,-1.,0.,true},{1.,-1.,1.,true},{1.,0.,-1.,true},{1.,0.,0.,true},{1.,0.,1.,true},
                     {1.,1.,-1.,true},{1.,1.,0.,true},{1.,1.,1.,true}};



/*================*
 * Build checking *
 *================*/


/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *   Setup sizes of the various vectors.
 */
static void mdlInitializeSizes(SimStruct *S)
{
    MYFILEPR=fopen("N2","w");  /* Opens a file, whose name is my_data_file, for writing */ 
    
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    if (!ssSetNumInputPorts(S, 2)) return; //aqui defino el numero de entradas
    //ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetInputPortWidth(S, 0, 3); // dimension de la entrada 
    ssSetInputPortDirectFeedThrough(S, 0, 1); //no la usamos, define el comportamiento de una entrada
    ssSetInputPortWidth(S, 1, 3); //
    ssSetInputPortDirectFeedThrough(S, 1, 1);
//     ssSetInputPortWidth(S, 2, 1); //
//     ssSetInputPortDirectFeedThrough(S, 2, 1);
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

    // for(i = 0; i < 10000; i++)
    // {
    //     seno[i] = sin(50 * (2 * M_PI) * i / 10000);
    // }
	
    R=2.;
    L=2e-3;
    vdc=5200.;
	Vph_rms = 3300./sqrt(3.);
	
	pref = (Vph_rms*Vph_rms)/R;
    ampli = Vph_rms/sqrt(R*R + (2.*M_PI*50.*L)*(2.*M_PI*50.*L));
    
    u_0[0]=0.;
    u_0[1]=0.;
    u_0[2]=0.;

    a = 1.;
    b = 0.;
    c = -1.;
    
    lambda=5e-3; 

    cont = 0;
	
	// A = 1-R*Tsampling/L;
	// B = Tsampling*vdc/L/2.;

    for(i=1;i<=27;i++)
    {
        for(k=1;k<=4;k++)
        {
            Uk2[i-1][k-1] = seq[i-1][k-1];
        }
    }
    // for(j=1;j<=27;j++)
    //             {
                    
                    
    //                     fprintf(MYFILEPR,"estados: %f %f %f \n", Uk2[j-1][0],Uk2[j-1][1],Uk2[j-1][2]);
                    
    //             }
    

}
#endif /*  MDL_START */


/* Function: mdlOutputs =======================================================
 * Abstract:
 *    
 */

static void Useq(float a, float b, float c)
{
    int i;

    for(i=1;i<=27;i++)
    {
        Uk2[i-1][3] = true;
    }

    if(a==-1)
    {
        for(i=19;i<=27;i++)
        {
            Uk2[i-1][3]=false;
        }
    }
    if(a==1)
    {
        for(i=1;i<=9;i++)
        {
            Uk2[i-1][3]=false;
        }
    }
    if(b==-1)
    {
        for(i=7;i<=9;i++)
        {
            Uk2[i-1][3]=false;
        }
        for(i=16;i<=18;i++)
        {
            Uk2[i-1][3]=false;
        }
        for(i=25;i<=27;i++)
        {
            Uk2[i-1][3]=false;
        }  
    }
    if(b==1)
    {
        for(i=1;i<=3;i++)
        {
            Uk2[i-1][3]=false;
        }
        for(i=10;i<=12;i++)
        {
            Uk2[i-1][3]=false;
        }
        for(i=19;i<=21;i++)
        {
            Uk2[i-1][3]=false;
        }  
    }
    if(c==-1)
    {
        for(i=3;i<=27;i=i+3)
        {
            Uk2[i-1][3]=false;
        }
    }
    if(c==1)
    {
        for(i=1;i<=27;i=i+3)
        {
            Uk2[i-1][3]=false;
        }
    }
}


static void mdlOutputs(SimStruct *S, int_T tid) //genera una salida cada vez q se llama al bloque
                                                //aqui es donde va el algoritmo que yo quiera ejecutar
{
    //int_T             i;
    InputRealPtrsType piabc = ssGetInputPortRealSignalPtrs(S,0);
    
	InputRealPtrsType pvabc = ssGetInputPortRealSignalPtrs(S,1);
//     
//     InputRealPtrsType pwr = ssGetInputPortRealSignalPtrs(S,2);
    
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
    
    
    // CONSTRUCCIÓN MATRIZ TOTAL DE SECUENCIAS


   
    float Uk[27][4]={{-1.,-1.,-1.,true},{-1.,-1.,0.,true},{-1.,-1.,1.,true},{-1.,0.,-1.,true},{-1.,0.,0.,true},{-1.,0.,1.,true},
                     {-1.,1.,-1.,true},{-1.,1.,0.,true},{-1.,1.,1.,true},{0.,-1.,-1.,true},{0.,-1.,0.,true},{0.,-1.,1.,true},
                     {0.,0.,-1.,true},{0.,0.,0.,true},{0.,0.,1.,true},{0.,1.,-1.,true},{0.,1.,0.,true},{0.,1.,1.,true},
                     {1.,-1.,-1.,true},{1.,-1.,0.,true},{1.,-1.,1.,true},{1.,0.,-1.,true},{1.,0.,0.,true},{1.,0.,1.,true},
                     {1.,1.,-1.,true},{1.,1.,0.,true},{1.,1.,1.,true}};
   
    float Jmin=INFINITY;
    
    // Medidas
    
    iabc[0] = *piabc[0];
    iabc[1] = *piabc[1];
    iabc[2] = *piabc[2];
	
	vabc[0] = *pvabc[0];
    vabc[1] = *pvabc[1];
    vabc[2] = *pvabc[2];
	
	// Transformadas
    
    v_aB[0] = Tabc_aB[0][0]*vabc[0]+Tabc_aB[0][1]*vabc[1]+Tabc_aB[0][2]*vabc[2];
    v_aB[1] = Tabc_aB[1][0]*vabc[0]+Tabc_aB[1][1]*vabc[1]+Tabc_aB[1][2]*vabc[2];

    i_aB[0] = Tabc_aB[0][0]*iabc[0]+Tabc_aB[0][1]*iabc[1]+Tabc_aB[0][2]*iabc[2];
    i_aB[1] = Tabc_aB[1][0]*iabc[0]+Tabc_aB[1][1]*iabc[1]+Tabc_aB[1][2]*iabc[2];
	
	// Calculo corrientes de referencia

    iref[0] = ampli*sin(50. * (2. * M_PI) * cont / 10000.);
    iref[1] = ampli*sin(50. * (2. * M_PI) * cont / 10000. + 2.*M_PI/3.);
    iref[2] = ampli*sin(50. * (2. * M_PI) * cont / 10000. + 4.*M_PI/3.);

    iref2[0] = ampli*sin(50. * (2. * M_PI) * (cont + 1.) / 10000.);
    iref2[1] = ampli*sin(50. * (2. * M_PI) * (cont + 1.) / 10000. + 2.*M_PI/3.);
    iref2[2] = ampli*sin(50. * (2. * M_PI) * (cont + 1.) / 10000. + 4.*M_PI/3.);

    iref_aB[0] = Tabc_aB[0][0]*iref[0]+Tabc_aB[0][1]*iref[1]+Tabc_aB[0][2]*iref[2];
    iref_aB[1] = Tabc_aB[1][0]*iref[0]+Tabc_aB[1][1]*iref[1]+Tabc_aB[1][2]*iref[2];

    iref2_aB[0] = Tabc_aB[0][0]*iref2[0]+Tabc_aB[0][1]*iref2[1]+Tabc_aB[0][2]*iref2[2];
    iref2_aB[1] = Tabc_aB[1][0]*iref2[0]+Tabc_aB[1][1]*iref2[1]+Tabc_aB[1][2]*iref2[2];
	
	// delta = (v_aB[0]*v_aB[0]+v_aB[1]*v_aB[1]);
    // if(delta == 0)
    //     delta = .001;

	
	// iref_aB[0] = pref*v_aB[0]/delta;
    // iref_aB[1] = pref*v_aB[1]/delta;
	
	// iref[0] = Tabc_aB[0][0]*iref_aB[0]+Tabc_aB[1][0]*iref_aB[1];
	// iref[1] = Tabc_aB[0][1]*iref_aB[0]+Tabc_aB[1][1]*iref_aB[1];
	// iref[2] = Tabc_aB[0][2]*iref_aB[0]+Tabc_aB[1][2]*iref_aB[1];
	
	// fprintf(MYFILEPR,"iaref: %f \n", iref[0]);
	// fprintf(MYFILEPR,"ibref: %f \n", iref[1]);
	// fprintf(MYFILEPR,"icref: %f \n", iref[2]);
    // fprintf(MYFILEPR,"ia: %f \n", iabc[0]);
	// fprintf(MYFILEPR,"ib: %f \n", iabc[1]);
	// fprintf(MYFILEPR,"ic: %f \n", iabc[2]);
	
    // fprintf(MYFILEPR," \n");
    // fprintf(MYFILEPR,"ialfa: %f \n", i_aB[0]);
	// fprintf(MYFILEPR,"ibeta: %f \n", i_aB[1]);

    // fprintf(MYFILEPR," \n");
    // fprintf(MYFILEPR,"irefa: %f \n", iref_aB[0]);
	// fprintf(MYFILEPR,"irefb: %f \n", iref_aB[1]);
    // fprintf(MYFILEPR," \n");

    
    //// MPC
    
    // Conjunto de transiciones posibles
    
    if(u_0[0]==-1)
    {
        for(i=19;i<=27;i++)
        {
            Uk[i-1][3]=false;
        }
    }
    if(u_0[0]==1)
    {
        for(i=1;i<=9;i++)
        {
            Uk[i-1][3]=false;
        }
    }
    if(u_0[1]==-1)
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
    if(u_0[1]==1)
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
    if(u_0[2]==-1)
    {
        for(i=3;i<=27;i=i+3)
        {
            Uk[i-1][3]=false;
        }
    }
    if(u_0[2]==1)
    {
        for(i=1;i<=27;i=i+3)
        {
            Uk[i-1][3]=false;
        }
    }
    
    
    // Calculo funcion de coste
    
    
    
    for(i=1;i<=27;i++)
    {
        if(Uk[i-1][3] == true)
        {
            u_act1 = Uk[i-1][0];
            u_act2 = Uk[i-1][1];
            u_act3 = Uk[i-1][2];

            iJ_aB[0] = (1./(R*Tsampling+L))*(L*i_aB[0]+Tsampling*vdc/3.*sqrt(3./2.)*(Tabc_aB[0][0]*Uk[i-1][0]+Tabc_aB[0][1]*Uk[i-1][1]+Tabc_aB[0][2]*Uk[i-1][2]));
            iJ_aB[1] = (1./(R*Tsampling+L))*(L*i_aB[1]+Tsampling*vdc/3.*sqrt(3./2.)*(Tabc_aB[1][0]*Uk[i-1][0]+Tabc_aB[1][1]*Uk[i-1][1]+Tabc_aB[1][2]*Uk[i-1][2]));

            ie[0] = iref_aB[0]-iJ_aB[0];
            ie[1] = iref_aB[1]-iJ_aB[1];

            incr_u=fabsf(Uk[i-1][0]-u_0[0])+
                   fabsf(Uk[i-1][1]-u_0[1])+
                   fabsf(Uk[i-1][2]-u_0[2]);
			
            if((ie[0]*ie[0]+ie[1]*ie[1] + lambda*incr_u) < Jmin)
            {
                Useq(u_act1,u_act2,u_act3);

                fprintf(MYFILEPR,"estado 1: %f %f %f \n", u_act1,u_act2,u_act3);
    
                for(k=1;k<=27;k++)
                {
                    // fprintf(MYFILEPR,"%d",j);
                    if(Uk2[k-1][3] == true)
                    {
                        fprintf(MYFILEPR,"posible: %f %f %f \n", Uk2[k-1][0],Uk2[k-1][1],Uk2[k-1][2]);
                    }
                }

                for(k=1;k<=27;k++) // SEGUNDO HORIZONTE Y MINIMIZACION
                {

                    if(Uk2[k-1][3] == true)
                    {

                        iJ2_aB[0] = (1./(R*Tsampling+L))*(L*iJ_aB[0]+Tsampling*vdc/3.*sqrt(3./2.)*(Tabc_aB[0][0]*Uk2[k-1][0]+Tabc_aB[0][1]*Uk2[k-1][1]+Tabc_aB[0][2]*Uk2[k-1][2]));
                        iJ2_aB[1] = (1./(R*Tsampling+L))*(L*iJ_aB[1]+Tsampling*vdc/3.*sqrt(3./2.)*(Tabc_aB[1][0]*Uk2[k-1][0]+Tabc_aB[1][1]*Uk2[k-1][1]+Tabc_aB[1][2]*Uk2[k-1][2]));

                        ie2[0] = iref2_aB[0]-iJ2_aB[0];
                        ie2[1] = iref2_aB[1]-iJ2_aB[1];

                        incr_u = incr_u+
                            fabsf(Uk2[k-1][0]-Uk[i-1][0])+
                            fabsf(Uk2[k-1][1]-Uk[i-1][1])+
                            fabsf(Uk2[k-1][2]-Uk[i-1][2]);

                        J=ie[0]*ie[0]+ie[1]*ie[1]+ie2[0]*ie2[0]+ie2[1]*ie2[1]+lambda*incr_u;
                    

                        if(J<Jmin)
                        {
                            Jmin = J;
                            a = Uk[i-1][0];
                            b = Uk[i-1][1];
                            c = Uk[i-1][2];
                            
                            fprintf(MYFILEPR,"a: %f \n", Uk[i-1][0]);
                            fprintf(MYFILEPR,"b: %f \n", Uk[i-1][1]);
                            fprintf(MYFILEPR,"c: %f \n", Uk[i-1][2]);
                            fprintf(MYFILEPR,"a2: %f \n", Uk2[k-1][0]);
                            fprintf(MYFILEPR,"b2: %f \n", Uk2[k-1][1]);
                            fprintf(MYFILEPR,"c2: %f \n", Uk2[k-1][2]);
                            fprintf(MYFILEPR,"cost1: %f \n", (ie[0]*ie[0]+ie[1]*ie[1] + lambda*incr_u));
                            fprintf(MYFILEPR,"cost1y2: %f \n", Jmin);

                        }
                    }
                }
            }
            else
            {
                fprintf(MYFILEPR,"cost1mayor: %f \n", (ie[0]*ie[0]+ie[1]*ie[1] + lambda*incr_u));
            }
            

            

            

            
        }
    }
    
   u_0[0]   = a;
   u_0[1]   = b;
   u_0[2]   = c; 
   // TORQUE=Xm/0.78/Xr*(flux_r*is_dq[1]);
   // TORQUE=Xm/0.78/Xr*(histphir_a[0]*is_aB[1]-histphir_B[0]*is_aB[0]);
   
   ua[0]   = a;
   ub[0]   = b;
   uc[0]   = c; 

//    pJ[0] = 2.*sin(50. * (2. * M_PI) * cont / 10000.);
//    pphir_dq[0] = 2.*sin(50. * (2. * M_PI) * cont / 10000. + 2.*M_PI/3.);
//    pphir_dq[1] = 2.*sin(50. * (2. * M_PI) * cont / 10000. + 4.*M_PI/3.);
   cont++;
   if(cont == 10000)
        cont = 0;

        
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
   
 fprintf(MYFILEPR,"S11: %f \n", S11[0]);
 fprintf(MYFILEPR,"S12: %f \n", S12[0]);
 fprintf(MYFILEPR,"S21: %f \n", S21[0]);
 fprintf(MYFILEPR,"S31: %f \n", S22[0]);
 fprintf(MYFILEPR,"S32: %f \n", S31[0]);  
 fprintf(MYFILEPR,"S22: %f \n", S32[0]);
    
}


/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S) //la ultima funciona que se ejecuta antes de terminar la simulacion
{
    fclose(MYFILEPR); /*Closes the file, my_data_file */
}



#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
