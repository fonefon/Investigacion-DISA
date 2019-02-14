/*
 * File : PYO.c
 * Abstract:
 *       PERTURBAR Y OBSERVAR
 *
 */


#define S_FUNCTION_NAME  pyo // esto debe tener el mismo nombre que en el simulink
#define S_FUNCTION_LEVEL 2

// #define TPWM 0.001 //Tpwm en seg



#include "simstruc.h"
#include <math.h>
#include <stdio.h>


/*Variables globales*/

float Vpv;
float Vpv_ant;

float Ppv;
float Ppv_ant;

float V_MPPT;

float deltaV;

//int flag;
/*================*
 * Build checking *
 *================*/


/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *   Setup sizes of the various vectors.
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    if (!ssSetNumInputPorts(S, 2)) return; //aqui defino el numero de entradas
    //ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetInputPortWidth(S, 0, 1); // dimension de la entrada #
    ssSetInputPortDirectFeedThrough(S, 0, 1); //no la usamos, define el comportamiento de una entrada
    ssSetInputPortWidth(S, 1, 1); //
    ssSetInputPortDirectFeedThrough(S, 1, 1);
//     ssSetInputPortWidth(S, 2, 1); //
//     ssSetInputPortDirectFeedThrough(S, 2, 1);

    if (!ssSetNumOutputPorts(S,1)) return; // se cuelga
    //ssSetOutputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetOutputPortWidth(S, 0, 1);


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
    ssSetSampleTime(S, 0, 0.25);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S); 
}

#define MDL_START                      /* Change to #undef to remove function */
#if defined(MDL_START)
/* Function: mdlStart ==========================================================
 * Abstract:
 *
 */
static void mdlStart(SimStruct *S) //este bloque solo se ejecuta 1 vez
{
   deltaV=1.;
   Ppv_ant=500.0;
   Vpv_ant=60.0;
   V_MPPT=60.0;
//    flag=0;
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
    InputRealPtrsType pVpv = ssGetInputPortRealSignalPtrs(S,0);
    
    InputRealPtrsType pPpv = ssGetInputPortRealSignalPtrs(S,1);
//     
//     InputRealPtrsType pflag = ssGetInputPortRealSignalPtrs(S,2);
    
    real_T            *y1    = ssGetOutputPortRealSignal(S,0);
    

    
    //Medidas
    
    Vpv = *pVpv[0];
    
    Ppv = *pPpv[0];
    
//     flag = *pflag[0];
   
    // Algoritmo

//     if(flag == 1)
//     {
     
        if(Ppv > Ppv_ant)        
        {
            if(Vpv > Vpv_ant)            
            {
                V_MPPT = V_MPPT + deltaV;
            }
            else          
            {
                V_MPPT = V_MPPT - deltaV;
            }
        }
        else
        {
            if(Ppv < Ppv_ant)
            {
                if(Vpv < Vpv_ant)            
                {
                V_MPPT = V_MPPT + deltaV;
                }
                else          
                {
                V_MPPT = V_MPPT - deltaV;
                }
            }
        }

        Ppv_ant = Ppv;

        Vpv_ant = Vpv;

       y1[0] = V_MPPT;
//     }
//     else
//     {
//        y1[0] = 60.0;
//     }
    
  
}


/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S) //la ultima funciona que se ejecuta antes de terminar la simulacion
{
}



#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
