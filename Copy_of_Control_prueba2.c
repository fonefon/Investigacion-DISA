#define S_FUNCTION_NAME  Copy_of_Control_prueba2
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
//#include "Control_c_cabecera.h"
#include "Copy_of_control_3lvl.c" ///////////////////////////////////////////////////////////////////////////////////////////////


static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 6);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;}
    ssSetSFcnParamTunable(S,0,SS_PRM_TUNABLE);
    ssSetSFcnParamTunable(S,1,SS_PRM_TUNABLE);
    ssSetSFcnParamTunable(S,2,SS_PRM_TUNABLE);
	ssSetSFcnParamTunable(S,3,SS_PRM_TUNABLE);
	ssSetSFcnParamTunable(S,4,SS_PRM_TUNABLE);
    ssSetSFcnParamTunable(S,5,SS_PRM_TUNABLE);
    
    ssSetNumContStates(S, 0);//No estados continuos
    ssSetNumDiscStates(S, 45);//No estados discretos

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, 26); //25 entradas en un bus + 1 modo funcionamiento ////////////////////////////////////////////////////////////////
    ssSetInputPortDirectFeedThrough(S, 0, 1);//Las salidas pueden depender directamente de las entradas

    if (!ssSetNumOutputPorts(S, 3)) return;
    ssSetOutputPortWidth(S, 0, 15); //15 salidas en un bus (5 tiempos por rama)
    ssSetOutputPortWidth(S, 1, 29); //29 salidas para la señal Debug
    ssSetOutputPortWidth(S, 2, 17);  //17 Contactores más resistencia/////////////////////////////////////////////////////////////////////////////////////////

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* Specify the sim state compliance to be same as a built-in block */
    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);

    ssSetOptions(S, 0);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, -1);//Sampleo inherente a la simulación (fijado por el paso fijo de la simulación)
    ssSetOffsetTime(S, 0, 0.0);

}

#define MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */
#if defined(MDL_INITIALIZE_CONDITIONS)
  static void mdlInitializeConditions(SimStruct *S)
  {
    //Aquí se inicializan los estados discretos, de haber alguno
  int i;
  real_T *x0    = ssGetRealDiscStates(S);
  int_T   nDStates = ssGetNumDiscStates(S);
  for (i = 0; i < nDStates; i++) {
    *x0++ = 0.0;     // Iniciamos todos los estados discretos a 0
  }      
  }
#endif /* MDL_INITIALIZE_CONDITIONS */
  
  static void mdlOutputs(SimStruct *S, int_T tid)
{

    InputRealPtrsType u_ptr = ssGetInputPortRealSignalPtrs(S,0);
    real_T      *y = ssGetOutputPortSignal(S,0);//Señales de tiempo de disparos (15)
    real_T      *y1= ssGetOutputPortSignal(S,1);//Señales de Debug (18)
    real_T      *y2= ssGetOutputPortSignal(S,2);//Señales de periféricos (6)
    real_T      *x0 = ssGetDiscStates(S);
	/*
    real_T     	Kcond 	= mxGetScalar(ssGetSFcnParam(S,0));
    real_T    	Kicond  = mxGetScalar(ssGetSFcnParam(S,1));
    real_T     	Kpv 	= mxGetScalar(ssGetSFcnParam(S,2));
    real_T    	Kiv 	= mxGetScalar(ssGetSFcnParam(S,3));
	real_T    	Kpc 	= mxGetScalar(ssGetSFcnParam(S,4));
    real_T      Kic 	= mxGetScalar(ssGetSFcnParam(S,5));
	
    real_T DELT;
    real_T *P_DELT = &DELT; //Se requiere crear primero una variable en MATLAB para poder asignarle el puntero
    
    *P_DELT = 0.0001;// PERIODO DE MUESTREO
  */
    // Las entradas de esta función son:
    //      inicio, Enable, Ia,  Ib,  Ic,  Vc1,  Vc2,  Vc3,   Vc4,  Vdcref,  Va,  Vb,  Vc, disp_a, disp_b, disp_c, Debug, Perifericos, DELT,    TIMEZERO,  TIME
//c3ndpc_(Kcond, Kpv, Kiv, Kpc, Kic, Kicond, u_ptr,y,y+5,y+10,y1,y2,P_DELT);
}
  
#define MDL_UPDATE  /* Change to #undef to remove function */
#if defined(MDL_UPDATE)
  /* Function: mdlUpdate ======================================================
   * Abstract:
   *    This function is called once for every major integration time step.
   *    Discrete states are typically updated here, but this function is useful
   *    for performing any tasks that should only take place once per
   *    integration step.
   */
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
	real_T     	Kcond 	= mxGetScalar(ssGetSFcnParam(S,0));
    real_T    	Kicond  = mxGetScalar(ssGetSFcnParam(S,1));
    real_T     	Kpv 	= mxGetScalar(ssGetSFcnParam(S,2));
    real_T    	Kiv 	= mxGetScalar(ssGetSFcnParam(S,3));
	real_T    	Kpc 	= mxGetScalar(ssGetSFcnParam(S,4));
    real_T      Kic 	= mxGetScalar(ssGetSFcnParam(S,5));	
	InputRealPtrsType u_ptr = ssGetInputPortRealSignalPtrs(S,0);
    real_T      *y = ssGetOutputPortSignal(S,0);//Señales de tiempo de disparos (15)
    real_T      *y1= ssGetOutputPortSignal(S,1);//Señales de Debug (18)
    real_T      *y2= ssGetOutputPortSignal(S,2);//Señales de periféricos (6)
    real_T      *x0 = ssGetDiscStates(S);	
      
	real_T DELT;
    real_T *P_DELT = &DELT; //Se requiere crear primero una variable en MATLAB para poder asignarle el puntero
    
    *P_DELT = 0.0001;// PERIODO DE MUESTREO    
      
    c3ndpc_(Kcond, Kpv, Kiv, Kpc, Kic, Kicond, u_ptr,y,y+5,y+10,y1,y2,P_DELT);   // Almacena las salidas en x0
  }
#endif /* MDL_UPDATE */
  
static void mdlTerminate(SimStruct *S)
{
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif