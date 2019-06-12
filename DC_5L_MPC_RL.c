#define S_FUNCTION_NAME  DC_5L_MPC_RL
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


float iabc[3],i_aB[2],iref[3],iref_aB[2],iJ_aB[2],iJ_abc[3],ie[2];
float vabc[3],v_aB[2];
float pref;
float Vph_rms;
float ampli;
float vdc,R,L,idc;
float Tsampling = 100e-6;
float u[3],incr_u,Uk[125*3];
int a,b,c;
float *fa, *fb, *fc;
float C1,C2,C3,C4,vd1,vd2,vd3;
float lambda,delta,J;


// float A,B;

float Tabc_aB[2][3];

int i,cont,j;



/*================*
 * Build checking *
 *================*/


/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *   Setup sizes of the various vectors.
 */
static void mdlInitializeSizes(SimStruct *S)
{
    MYFILEPR=fopen("N1","w");  /* Opens a file, whose name is my_data_file, for writing */ 
    
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    if (!ssSetNumInputPorts(S, 3)) return; //aqui defino el numero de entradas
    //ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetInputPortWidth(S, 0, 3); // dimension de la entrada 
    ssSetInputPortDirectFeedThrough(S, 0, 1); //no la usamos, define el comportamiento de una entrada
    ssSetInputPortWidth(S, 1, 3); //
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortWidth(S, 2, 1); //
    ssSetInputPortDirectFeedThrough(S, 2, 1);
//     ssSetInputPortWidth(S, 3, 1); //
//     ssSetInputPortDirectFeedThrough(S, 3, 1);


    if (!ssSetNumOutputPorts(S,12)) return; //
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
    ssSetOutputPortWidth(S, 10, 1);
    ssSetOutputPortWidth(S, 11, 1);
    
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
static void perm(){	
	int temp;
    int numbers=3;
    float a[4];
	int upto = 2, temp2;
	int index = 0,i;
	for( temp2 = 1 ; temp2 <= numbers; temp2++)
    {
		a[temp2]=-2.;
	}
	a[numbers]=-2.-1.;
	temp=numbers;
	while(1)
    {
		if(a[temp]==upto)
        {
		    temp--;
		    if(temp==0)
			    break;
		}
		else
        {
	        a[temp]++;
	        while(temp<numbers)
            {
		        temp++;
		        a[temp]=-2.;
		    }	
		
		
            //printf("(");
	        for( temp2 = 1 ; temp2 <= numbers; temp2++)
            {
			    Uk[index] = a[temp2];
			    index++;
		    }
			//printf(")");
		}
    }
}

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
	
    R=240.;
    L=2e-3;
    vdc=750.;
	Vph_rms = 230.;
	
    ampli = Vph_rms/sqrt(R*R + (2.*M_PI*50.*L)*(2.*M_PI*50.*L));
    
    u[0]=0.;
    u[1]=0.;
    u[2]=0.;

    a = 1.;
    b = 0.;
    c = -1.;
    
    lambda=1e-2; 

    cont = 0;
    perm();
	
	// A = 1-R*Tsampling/L;
	// B = Tsampling*vdc/L/2.;
    

}
#endif /*  MDL_START */


/* Function: mdlOutputs =======================================================
 * Abstract:
 *    
 */

static float *efes(float u)
{
    static float f[4];

    if(u == 2)
        f[0] = 1.;
    else
    {
        f[1] = 0.;
        f[2] = 0.;
        f[3] = 0.;
    }

    if(u == 1)
        f[1] = 1.;
    else
    {
        f[0] = 0.;
        f[2] = 0.;
        f[3] = 0.;
    }

    if(u == 0)
    {
        f[0] = 0.;
        f[1] = 0.;
        f[2] = 0.;
        f[3] = 0.;
    }

    if(u == -1)
        f[2] = 1.;
    else
    {
        f[1] = 0.;
        f[0] = 0.;
        f[3] = 0.;
    }
    if(u == -2)
        f[3] = 1.;
    else
    {
        f[1] = 0.;
        f[2] = 0.;
        f[0] = 0.;
    }
    return f;
}
	


static void mdlOutputs(SimStruct *S, int_T tid) //genera una salida cada vez q se llama al bloque
                                                //aqui es donde va el algoritmo que yo quiera ejecutar
{
    //int_T             i;
    InputRealPtrsType piabc = ssGetInputPortRealSignalPtrs(S,0);
    
	InputRealPtrsType pvabc = ssGetInputPortRealSignalPtrs(S,1);
//     
    InputRealPtrsType pidc = ssGetInputPortRealSignalPtrs(S,2);
    
//     InputRealPtrsType pteta = ssGetInputPortRealSignalPtrs(S,3);
    
    real_T            *Sa1    = ssGetOutputPortRealSignal(S,0);
    
    real_T            *Sa2    = ssGetOutputPortRealSignal(S,1);
    
    real_T            *Sa3    = ssGetOutputPortRealSignal(S,2);
    
    real_T            *Sa4    = ssGetOutputPortRealSignal(S,3);
    
    real_T            *Sb1    = ssGetOutputPortRealSignal(S,4);
    
    real_T            *Sb2    = ssGetOutputPortRealSignal(S,5);
    
    real_T            *Sb3    = ssGetOutputPortRealSignal(S,6);
    
    real_T            *Sb4    = ssGetOutputPortRealSignal(S,7);

    real_T            *Sc1    = ssGetOutputPortRealSignal(S,8);
    
    real_T            *Sc2    = ssGetOutputPortRealSignal(S,9);
    
    real_T            *Sc3    = ssGetOutputPortRealSignal(S,10);
    
    real_T            *Sc4    = ssGetOutputPortRealSignal(S,11);
   
    // Vector de transiciones
    
    
    
   
   
    float Jmin=INFINITY;
    
    // Medidas
    
    iabc[0] = *piabc[0];
    iabc[1] = *piabc[1];
    iabc[2] = *piabc[2];
	
	vabc[0] = *pvabc[0];
    vabc[1] = *pvabc[1];
    vabc[2] = *pvabc[2];

    idc = *pidc[0];
	
	// Transformadas
    
    v_aB[0] = Tabc_aB[0][0]*vabc[0]+Tabc_aB[0][1]*vabc[1]+Tabc_aB[0][2]*vabc[2];
    v_aB[1] = Tabc_aB[1][0]*vabc[0]+Tabc_aB[1][1]*vabc[1]+Tabc_aB[1][2]*vabc[2];

    i_aB[0] = Tabc_aB[0][0]*iabc[0]+Tabc_aB[0][1]*iabc[1]+Tabc_aB[0][2]*iabc[2];
    i_aB[1] = Tabc_aB[1][0]*iabc[0]+Tabc_aB[1][1]*iabc[1]+Tabc_aB[1][2]*iabc[2];
	
	// Calculo corrientes de referencia

    iref[0] = ampli*sin(50. * (2. * M_PI) * cont / 10000.);
    iref[1] = ampli*sin(50. * (2. * M_PI) * cont / 10000. + 2.*M_PI/3.);
    iref[2] = ampli*sin(50. * (2. * M_PI) * cont / 10000. + 4.*M_PI/3.);

    iref_aB[0] = Tabc_aB[0][0]*iref[0]+Tabc_aB[0][1]*iref[1]+Tabc_aB[0][2]*iref[2];
    iref_aB[1] = Tabc_aB[1][0]*iref[0]+Tabc_aB[1][1]*iref[1]+Tabc_aB[1][2]*iref[2];
	
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
    
    
    
    // Calculo funcion de coste
    
//     fprintf(MYFILEPR,"estado actual: %f %f %f \n", uabc[0],uabc[1],uabc[2]);
    
    // for(i=1;i<=27;i++)
    // {
    //     fprintf(MYFILEPR,"%d",i);
    //     // if(Uk[i-1][3] == true)
    //     // {
    //     //         fprintf(MYFILEPR,"posible: %f %f %f \n", Uk[i-1][0],Uk[i-1][1],Uk[i-1][2]);
    //     // }
    // }
    
    for(i=1;i<=125;i++)
    {
        if((fabsf(Uk[(i-1)*3] - u[0]) <= 1) && (fabsf(Uk[(i-1)*3 + 1] - u[1]) <= 1) && (fabsf(Uk[(i-1)*3 + 2] - u[2]) <= 1))
        {

            iJ_aB[0] = (1./(R*Tsampling+L))*(L*i_aB[0]+Tsampling*vdc/6.*sqrt(3./2.)*(Tabc_aB[0][0]*Uk[(i-1)*3]+Tabc_aB[0][1]*Uk[(i-1)*3 + 1]+Tabc_aB[0][2]*Uk[(i-1)*3 + 2]));
            iJ_aB[1] = (1./(R*Tsampling+L))*(L*i_aB[1]+Tsampling*vdc/6.*sqrt(3./2.)*(Tabc_aB[1][0]*Uk[(i-1)*3]+Tabc_aB[1][1]*Uk[(i-1)*3 + 1]+Tabc_aB[1][2]*Uk[(i-1)*3 + 2]));
			
            // ERROR REFERENCIA

            ie[0] = iref_aB[0]-iJ_aB[0];
            ie[1] = iref_aB[1]-iJ_aB[1];

            // COSTE DISPAROS

            incr_u=fabsf(Uk[(i-1)*3] - u[0])+
                   fabsf(Uk[(i-1)*3 + 1] - u[1])+
                   fabsf(Uk[(i-1)*3 + 2] - u[2]);

            // DESBALANCEO CONDENSADORES

            iJ_abc[0] = Tabc_aB[0][0]*iJ_aB[0]+Tabc_aB[1][0]*iJ_aB[1];
	        iJ_abc[1] = Tabc_aB[0][1]*iJ_aB[0]+Tabc_aB[1][1]*iJ_aB[1];
	        iJ_abc[2] = Tabc_aB[0][2]*iJ_aB[0]+Tabc_aB[1][2]*iJ_aB[1];

            fa = efes(Uk[(i-1)*3]);
            fb = efes(Uk[(i-1)*3 + 1]);
            fc = efes(Uk[(i-1)*3 + 2]);

            C1 = idc - *(fa)*iJ_abc[0] - *(fb)*iJ_abc[1] - *(fc)*iJ_abc[2];
            C2 = idc - *(fa)*iJ_abc[0] - *(fb)*iJ_abc[1] - *(fc)*iJ_abc[2] - *(fa + 1)*iJ_abc[0] - *(fb + 1)*iJ_abc[1] - *(fc + 1)*iJ_abc[2];
            C3 = idc + *(fa + 2)*iJ_abc[0] + *(fb + 2)*iJ_abc[1] + *(fc + 2)*iJ_abc[2] + *(fa + 3)*iJ_abc[0] + *(fb + 3)*iJ_abc[1] + *(fc + 3)*iJ_abc[2];
            C4 = idc + *(fa + 3)*iJ_abc[0] + *(fb + 3)*iJ_abc[1] + *(fc + 3)*iJ_abc[2];

            vd1 = C1 - C4;
            vd2 = C2 - C3;
            vd3 = C3 - C4;

            vd1 = fabsf(vd1);
            vd2 = fabsf(vd2);
            vd3 = fabsf(vd3);

            J=ie[0]*ie[0]+ie[1]*ie[1] + lambda*incr_u + vd1+vd2+vd3;

            if(J<Jmin)
            {
                Jmin = J;
                a = Uk[(i-1)*3];
                b = Uk[(i-1)*3 + 1];
                c = Uk[(i-1)*3 + 2];
                 fprintf(MYFILEPR,"cost: %f \n", Jmin);
                 fprintf(MYFILEPR,"a: %d \n", a);
                 fprintf(MYFILEPR,"b: %d \n", b);
                 fprintf(MYFILEPR,"c: %d \n", c); 
            }
        }
    }
    
   u[0]   = a;
   u[1]   = b;
   u[2]   = c; 
   // TORQUE=Xm/0.78/Xr*(flux_r*is_dq[1]);
   // TORQUE=Xm/0.78/Xr*(histphir_a[0]*is_aB[1]-histphir_B[0]*is_aB[0]);
   
  

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

       case 2 :
          Sa1[0]=1;
          Sa2[0]=1;
          Sa3[0]=1;
          Sa4[0]=1;
          break;

       case 1 :
          Sa1[0]=0;
          Sa2[0]=1;
          Sa3[0]=1;
          Sa4[0]=1;
          break; 

       case 0  :
          Sa1[0]=0;
          Sa2[0]=0;
          Sa3[0]=1;
          Sa4[0]=1;
          break; 

       case -1 :
          Sa1[0]=0;
          Sa2[0]=0;
          Sa3[0]=0;
          Sa4[0]=1;
          break;

       case -2 :
          Sa1[0]=0;
          Sa2[0]=0;
          Sa3[0]=0;
          Sa4[0]=0;
          break; 
     

       default : 
          Sa1[0]=100;
          Sa2[0]=100;
          Sa3[0]=100;
          Sa4[0]=100;
   }  
   switch(b) {

       case 2 :
          Sb1[0]=1;
          Sb2[0]=1;
          Sb3[0]=1;
          Sb4[0]=1;
          break;

       case 1 :
          Sb1[0]=0;
          Sb2[0]=1;
          Sb3[0]=1;
          Sb4[0]=1;
          break; 

       case 0  :
          Sb1[0]=0;
          Sb2[0]=0;
          Sb3[0]=1;
          Sb4[0]=1;
          break; 

       case -1 :
          Sb1[0]=0;
          Sb2[0]=0;
          Sb3[0]=0;
          Sb4[0]=1;
          break;

       case -2 :
          Sb1[0]=0;
          Sb2[0]=0;
          Sb3[0]=0;
          Sb4[0]=0;
          break; 
     

       default : 
          Sb1[0]=100;
          Sb2[0]=100;
          Sb3[0]=100;
          Sb4[0]=100;
   }  
   switch(c) {
       case 2 :
          Sc1[0]=1;
          Sc2[0]=1;
          Sc3[0]=1;
          Sc4[0]=1;
          break;

       case 1 :
          Sc1[0]=0;
          Sc2[0]=1;
          Sc3[0]=1;
          Sc4[0]=1;
          break; 

       case 0  :
          Sc1[0]=0;
          Sc2[0]=0;
          Sc3[0]=1;
          Sc4[0]=1;
          break; 

       case -1 :
          Sc1[0]=0;
          Sc2[0]=0;
          Sc3[0]=0;
          Sc4[0]=1;
          break;

       case -2 :
          Sc1[0]=0;
          Sc2[0]=0;
          Sc3[0]=0;
          Sc4[0]=0;
          break; 
     

       default : 
          Sc1[0]=100;
          Sc2[0]=100;
          Sc3[0]=100;
          Sc4[0]=100;
   }  
   
//  fprintf(MYFILEPR,"S11: %f \n", S11[0]);
//  fprintf(MYFILEPR,"S12: %f \n", S12[0]);
//  fprintf(MYFILEPR,"S21: %f \n", S21[0]);
//  fprintf(MYFILEPR,"S31: %f \n", S22[0]);
//  fprintf(MYFILEPR,"S32: %f \n", S31[0]);  
//  fprintf(MYFILEPR,"S22: %f \n", S32[0]);
    
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
