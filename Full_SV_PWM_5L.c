/*
 *
 */


#define S_FUNCTION_NAME  Full_SV_PWM_5L 
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

float ciclos = 100.;
float Tsampling = 100e-6/100.;
float cont;
char aux;

float i_abc[3];
float i_aBY[3];
float iref_aBY[3];
float ie_aBY[3];
float inte_ia,inte_iB;

float L = 2e-3;

float v_abc[3];
float v_aBY[3];

float vc1,vc2,vc3,vc4;
float vd1,vd2,vd3;
float dvd1,dvd2,dvd3;
float term1,term2,term3;
float sign1,sign2,sign3;
float Vdc;
float Vdcref = 800.;
float inte_Vdc;

float fioj[15],dioj[15],d_aBY[15];
float dyj[15];
float u[8];
float u1z,u2z;
float damax,dbmax,dcmax;
float max1,max2,max4,max5;
char flaga,flagb,flagc;
float Total;

float pref;
float qref = 0.;
float p,q;
float inte_p,inte_q;

float Tabc_aBY[3][3];

float histie_a[2],histie_B[2];
float histv_a[2],histv_B[2];
float histu1[2];
float histu2[2];

float kp = 3e-7;
float kpi = 5e-5;
float kq = 3e-7;
float kqi = 5e-5;
float kpdc = .05;
float kidc = 1.;
float k1 = 5e-3;
float k2 = 5e-3;
float k3 = 5e-3;

float ky1 = .7;
float ky2 = .1;
float ky4 = .1;
float ky5 = .7;

int i,imax,j,i1,i2;


/*================*
 * Build checking *
 *================*/


/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *   Setup sizes of the various vectors.
 */
static void mdlInitializeSizes(SimStruct *S)
{
     MYFILEPR=fopen("my_data_file2","w");  /* Opens a file, whose name is my_data_file, for writing */ 
    
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    if (!ssSetNumInputPorts(S, 8)) return; //aqui defino el numero de entradas
    //ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetInputPortWidth(S, 0, 3); // dimension de la entrada 
    ssSetInputPortDirectFeedThrough(S, 0, 1); //no la usamos, define el comportamiento de una entrada
    ssSetInputPortWidth(S, 1, 3); 
    ssSetInputPortDirectFeedThrough(S, 1, 1); 
	ssSetInputPortWidth(S, 2, 1); //
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    ssSetInputPortWidth(S, 3, 1); //
    ssSetInputPortDirectFeedThrough(S, 3, 1);
    ssSetInputPortWidth(S, 4, 1); //
    ssSetInputPortDirectFeedThrough(S, 4, 1);
	ssSetInputPortWidth(S, 5, 1); //
    ssSetInputPortDirectFeedThrough(S, 5, 1);
	ssSetInputPortWidth(S, 6, 1); //
    ssSetInputPortDirectFeedThrough(S, 6, 1);
	ssSetInputPortWidth(S, 7, 3); //
    ssSetInputPortDirectFeedThrough(S, 7, 1);


    if (!ssSetNumOutputPorts(S,22)) return; //
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
	ssSetOutputPortWidth(S, 12, 1);
    ssSetOutputPortWidth(S, 13, 1);
	ssSetOutputPortWidth(S, 14, 1);
	ssSetOutputPortWidth(S, 15, 1);
	ssSetOutputPortWidth(S, 16, 1);
	ssSetOutputPortWidth(S, 17, 3);
	ssSetOutputPortWidth(S, 18, 5);
	ssSetOutputPortWidth(S, 19, 5);
	ssSetOutputPortWidth(S, 20, 5);
	ssSetOutputPortWidth(S, 21, 3);
	
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
    Tabc_aBY[0][0]=sqrt(2./3.);
    Tabc_aBY[0][1]=-0.5*sqrt(2./3.);
    Tabc_aBY[0][2]=-0.5*sqrt(2./3.);
    Tabc_aBY[1][0]=0.;
    Tabc_aBY[1][1]=sqrt(2./3.)*sqrt(3.)/2.;
    Tabc_aBY[1][2]=-sqrt(2./3.)*sqrt(3.)/2.;
	Tabc_aBY[2][0]=1./sqrt(3.);
    Tabc_aBY[2][1]=1./sqrt(3.);
    Tabc_aBY[2][2]=1./sqrt(3.);
	
	cont=0.;

	
    // histie_a[0]=0.;
    // histie_a[1]=0.;
    // histv_a[0]=0.;
    // histv_a[1]=0.;
	// histie_B[0]=0.;
    // histie_B[1]=0.;
    // histv_B[0]=0.;
    // histv_B[1]=0.;
    // histu1[0]=0.;
    // histu2[0]=0.;
	// histu1[1]=0.;
    // histu2[1]=0.;
	
	// inte_ia = 0.;
	// inte_iB = 0.;
	
	inte_p = 0.;
	inte_q = 0.;
	
	inte_Vdc = 0.;
	
	for(i=1;i<=15;i++)
		{
			dioj[i-1] = 0.;
		}

}
#endif /*  MDL_START */


/* Function: mdlOutputs =======================================================
 * Abstract:
 *    
 */
static void mdlOutputs(SimStruct *S, int_T tid) //genera una salida cada vez q se llama al bloque
                                                //aqui es donde va el algoritmo que yo quiera ejecutar
{
	
    InputRealPtrsType pi_abc= ssGetInputPortRealSignalPtrs(S,0);
	InputRealPtrsType pv_abc= ssGetInputPortRealSignalPtrs(S,1);
    
    InputRealPtrsType pvc1 = ssGetInputPortRealSignalPtrs(S,2);
	InputRealPtrsType pvc2 = ssGetInputPortRealSignalPtrs(S,3);
	InputRealPtrsType pvc3 = ssGetInputPortRealSignalPtrs(S,4);
	InputRealPtrsType pvc4 = ssGetInputPortRealSignalPtrs(S,5);
	
	InputRealPtrsType pvdc = ssGetInputPortRealSignalPtrs(S,6);
	
	InputRealPtrsType pvd = ssGetInputPortRealSignalPtrs(S,7);
    
    
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
    
	real_T            *pvd1    = ssGetOutputPortRealSignal(S,12);
	real_T            *pvd2    = ssGetOutputPortRealSignal(S,13);
    real_T            *pvd3    = ssGetOutputPortRealSignal(S,14);
	
	real_T            *ppot    = ssGetOutputPortRealSignal(S,15);
	real_T            *pu      = ssGetOutputPortRealSignal(S,16);
	real_T            *pterm   = ssGetOutputPortRealSignal(S,17);
	
	real_T            *pda  = ssGetOutputPortRealSignal(S,18);
	real_T            *pdb  = ssGetOutputPortRealSignal(S,19);
	real_T            *pdc  = ssGetOutputPortRealSignal(S,20);
	real_T            *psat  = ssGetOutputPortRealSignal(S,21);
   
   
	
	
    
    //Medidas
    
    i_abc[0] = *pi_abc[0];
    i_abc[1] = *pi_abc[1];
    i_abc[2] = *pi_abc[2];
	
	v_abc[0] = *pv_abc[0];
    v_abc[1] = *pv_abc[1];
    v_abc[2] = *pv_abc[2];
    
    vc1 = *pvc1[0];
    vc2 = *pvc2[0];
	vc3 = *pvc3[0];
	vc4 = *pvc4[0];
	
	vd1 = *pvd[0];
	vd2 = *pvd[1];
	vd3 = *pvd[2];
	
	Vdc = *pvdc[0];
	
	if(cont >= (ciclos - 1.))
    {
		
        cont = 0.;

		// Tensión entre condensadores
		
		// vd1 = vc4 - vc1;
		// vd2 = vc3 - vc2;
		// vd3 = vc2 - vc1;
		
		// Transformaciones
		
		i_aBY[0]=Tabc_aBY[0][0]*i_abc[0]+Tabc_aBY[0][1]*i_abc[1]+Tabc_aBY[0][2]*i_abc[2];
		i_aBY[1]=Tabc_aBY[1][0]*i_abc[0]+Tabc_aBY[1][1]*i_abc[1]+Tabc_aBY[1][2]*i_abc[2];
		i_aBY[2]=Tabc_aBY[2][0]*i_abc[0]+Tabc_aBY[2][1]*i_abc[1]+Tabc_aBY[2][2]*i_abc[2];
		
		v_aBY[0]=Tabc_aBY[0][0]*v_abc[0]+Tabc_aBY[0][1]*v_abc[1]+Tabc_aBY[0][2]*v_abc[2];
		v_aBY[1]=Tabc_aBY[1][0]*v_abc[0]+Tabc_aBY[1][1]*v_abc[1]+Tabc_aBY[1][2]*v_abc[2];
		v_aBY[2]=Tabc_aBY[2][0]*v_abc[0]+Tabc_aBY[2][1]*v_abc[1]+Tabc_aBY[2][2]*v_abc[2];
		
		// histv_a[0]=v_aBY[0];
		// histv_B[0]=v_aBY[1];
		
		// for(i=1;i<=5;i++)
		// {
			// d_aBY[i-1]=Tabc_aBY[0][0]*dioj[i-1]+Tabc_aBY[0][1]*dioj[i-1+5]+Tabc_aBY[0][2]*dioj[i-1+10];
			// d_aBY[i-1+5]=Tabc_aBY[1][0]*dioj[i-1]+Tabc_aBY[1][1]*dioj[i-1+5]+Tabc_aBY[1][2]*dioj[i-1+10];
			// d_aBY[i-1+10]=Tabc_aBY[2][0]*dioj[i-1]+Tabc_aBY[2][1]*dioj[i-1+5]+Tabc_aBY[2][2]*dioj[i-1+10];
		// }
		
		// Calculo potencias y uz
		
		p = v_aBY[0]*i_aBY[0]+v_aBY[1]*i_aBY[1];
		q = v_aBY[0]*i_aBY[1]-v_aBY[1]*i_aBY[0];//???
		
		inte_Vdc += 100e-6*(Vdcref*Vdcref-Vdc*Vdc);
		
		pref=kpdc*(Vdcref*Vdcref-Vdc*Vdc)+kidc*inte_Vdc;
		
		// pref = -pref;
		
		u1z=(4./Vdc)*((1.+2.*M_PI*50.*L*q/(v_aBY[0]*v_aBY[0]+v_aBY[1]*v_aBY[1]))*v_aBY[0]+2.*M_PI*50.*L*p*v_aBY[1]/(v_aBY[0]*v_aBY[0]+v_aBY[1]*v_aBY[1]));
		u2z=(4./Vdc)*((1.+2.*M_PI*50.*L*q/(v_aBY[0]*v_aBY[0]+v_aBY[1]*v_aBY[1]))*v_aBY[1]-2.*M_PI*50.*L*p*v_aBY[0]/(v_aBY[0]*v_aBY[0]+v_aBY[1]*v_aBY[1]));		
		
		// Calculo corrientes referencia
		
		// iref_aBY[0]=(v_aBY[0]*pref-v_aBY[1]*qref)/(v_aBY[0]*v_aBY[0]+v_aBY[1]*v_aBY[1]);
		// iref_aBY[1]=(v_aBY[1]*pref+v_aBY[0]*qref)/(v_aBY[0]*v_aBY[0]+v_aBY[1]*v_aBY[1]);
		
		// Errores
		
		// ie_aBY[0]=iref_aBY[0]-i_aBY[0];
		// ie_aBY[1]=iref_aBY[1]-i_aBY[1];
		
		// histie_a[0]=ie_aBY[0];
		// histie_B[0]=ie_aBY[1];
		
		// inte_ia += 100e-6*ie_aBY[0];
		// inte_iB += 100e-6*ie_aBY[1];
		
		inte_p += 100e-6*(p-pref);
		inte_q += 100e-6*(q-qref);
		
		
		
		// Controladores u1 y u2. Calculo del resto de u
		
		u[0]=u1z+kp*v_aBY[0]*(p-pref)+kpi*v_aBY[0]*inte_p-kq*v_aBY[1]*(q-qref)-kqi*v_aBY[1]*inte_q;
		u[1]=u2z+kp*v_aBY[1]*(p-pref)+kpi*v_aBY[1]*inte_p+kq*v_aBY[0]*(q-qref)+kqi*v_aBY[0]*inte_q;
		
		
		
		// histu1[0]=(4./Vdc)*(kp*histie_a[0]+(100e-6*ki-kp)*histie_a[1]+histv_a[0]-histv_a[1])+histu1[1];
		// histu2[0]=(4./Vdc)*(kp*histie_B[0]+(100e-6*ki-kp)*histie_B[1]+histv_B[0]-histv_B[1])+histu2[1];
		// histu1[1]=histu1[0];
		// histu2[1]=histu2[0];
		// histie_a[1]=histie_a[0];
		// histie_B[1]=histie_B[0];
		// histv_a[1]=histv_a[0];
		// histv_B[1]=histv_B[0];
		
		// u[0]=histu1[0];
		// u[1]=histu2[0];
		u[2]=k1*i_aBY[0]*vd1;
		u[3]=k1*i_aBY[1]*vd1;
		u[4]=k2*i_aBY[0]*vd2;
		u[5]=k2*i_aBY[1]*vd2;
		u[6]=k3*i_aBY[0]*vd3;
		u[7]=k3*i_aBY[1]*vd3;
		
		// Cambio de variables. u -> d
		
		d_aBY[0]=(-u[0]+3.*u[2]-u[4]-2.*u[6])/4.;
		d_aBY[5]=(-u[1]+3.*u[3]-u[5]-2.*u[7])/4.;
		d_aBY[1]=-u[2]+u[4]+u[6];
		d_aBY[6]=-u[3]+u[5]+u[7];
		d_aBY[3]=-u[6];
		d_aBY[8]=-u[7];
		d_aBY[4]=(u[0]+u[2]+u[4]+2.*u[6])/4.;
		d_aBY[9]=(u[1]+u[3]+u[5]+2.*u[7])/4.;
		
		d_aBY[10] = ky1;
		d_aBY[11] = ky2;
		d_aBY[13] = ky4;
		d_aBY[14] = ky5;
		
		
		fprintf(MYFILEPR,"p: %f  \n", p);//u[0]);
		fprintf(MYFILEPR,"q: %f  \n", q);//u[1]);
		fprintf(MYFILEPR,"pref: %f  \n", pref);
		fprintf(MYFILEPR,"u1: %f  \n", u1z);//u[0]);
		fprintf(MYFILEPR,"u2: %f  \n", u2z);//u[1]);
		
		
		// Algoritmo para la selección de los valores dyoj
		
		
		
		
		/* for(i=1;i<=5 && i!=3;i++)
		{
			dyj[i-1] = -sqrt(2.)*d_aBY[i-1];
			dyj[i-1+5] = d_aBY[i-1]/sqrt(2.) - d_aBY[i-1+5]*sqrt(3./2.);
			dyj[i-1+10] = d_aBY[i-1]/sqrt(2.) + d_aBY[i-1+5]*sqrt(3./2.);
		}
		
		// j = 4,5
		
		max5 = dyj[5-1];
		i1 = 0;
		if(dyj[5-1+5] > max5){max5 = dyj[5-1+5];i1 = 5;}
		if(dyj[5-1+10] > max5){max5 = dyj[5-1+10];i1 = 10;}
		
		d_aBY[14] = max5;
		
		max4 = dyj[4-1];
		i2 = 0;
		if(dyj[4-1+5] > max4){max4 = dyj[4-1+5];i2 = 5;}
		if(dyj[4-1+10] > max4){max4 = dyj[4-1+10];i2 = 10;}
		
		if(i2 == i1)
		{
			d_aBY[13] = max4;
		}
		else
		{
			d_aBY[13] = ky4;
		}
		
		// j = 1,2
		
		max1 = dyj[1-1];
		i1 = 0;
		if(dyj[1-1+5] > max1){max1 = dyj[1-1+5];i1 = 5;}
		if(dyj[1-1+10] > max1){max1 = dyj[1-1+10];i1 = 10;}
		
		d_aBY[10] = max1;
		
		max2 = dyj[2-1];
		i2 = 0;
		if(dyj[2-1+5] > max2){max2 = dyj[2-1+5];i2 = 5;}
		if(dyj[2-1+10] > max2){max2 = dyj[2-1+10];i2 = 10;}
		
		if(i2 == i1)
		{
			d_aBY[11] = max2;
		}
		else
		{
			d_aBY[11] = ky2;
		} */
		
		// Destransformada
		
		for(i=1;i<=5;i++)
		{
			if(i != 3)
			{
				dioj[i-1]=Tabc_aBY[0][0]*d_aBY[i-1]+Tabc_aBY[1][0]*d_aBY[i-1+5]+Tabc_aBY[2][0]*d_aBY[i-1+10];
				dioj[i-1+5]=Tabc_aBY[0][1]*d_aBY[i-1]+Tabc_aBY[1][1]*d_aBY[i-1+5]+Tabc_aBY[2][1]*d_aBY[i-1+10];
				dioj[i-1+10]=Tabc_aBY[0][2]*d_aBY[i-1]+Tabc_aBY[1][2]*d_aBY[i-1+5]+Tabc_aBY[2][2]*d_aBY[i-1+10];
			}
		}
		
		  for(i=1;i<=15;i++)
			{
     
      
                fprintf(MYFILEPR,"diof %d: %f  \n", i,dioj[i-1]);
      
			}
		
		// Saturacion al 0

		for(i=1;i<=15;i++)
		{
			if(dioj[i-1] < 0.){dioj[i-1] = 0.;}
		}
		
		// Restricciones
		
		if ((dioj[0]+dioj[1]+dioj[3]+dioj[4]) > 1.)
            {
                Total = dioj[0]+dioj[1]+dioj[3]+dioj[4];
                
                dioj[0] = dioj[0]/Total;
                dioj[1] = dioj[1]/Total;
				dioj[3] = dioj[3]/Total;
				dioj[4] = dioj[4]/Total;
            }
        if ((dioj[5]+dioj[6]+dioj[8]+dioj[9]) > 1.)
            {
                Total = dioj[5]+dioj[6]+dioj[8]+dioj[9];
                
                dioj[5] = dioj[5]/Total;
                dioj[6] = dioj[6]/Total;
				dioj[8] = dioj[8]/Total;
				dioj[9] = dioj[9]/Total;
            }
        if ((dioj[10]+dioj[11]+dioj[13]+dioj[14]) > 1.)
            {
                Total = dioj[10]+dioj[11]+dioj[13]+dioj[14];
                
                dioj[10] = dioj[10]/Total;
                dioj[11] = dioj[11]/Total;
				dioj[13] = dioj[13]/Total;
				dioj[14] = dioj[14]/Total;
            }
		// damax = 0.;
		// dbmax = 0.;
		// dcmax = 0.;
		// flaga = 0;
		// flagb = 0;
		// flagc = 0;
		
		// for(i=1;i<=5;i++)
		// {
			// if(i != 3)
			// {
				// if(dioj[i-1] < 1. && dioj[i-1] > 0. && flaga == 1){dioj[i-1] = 0.;}
				// if(dioj[i-1] > 1.)
				// {
					// flaga = 1;
					// if(dioj[i-1] > damax)
					// {
						// damax = dioj[i-1];
						// dioj[i-1] = 1.;
						// imax = i;
						// for(j=1;j<imax;j++){dioj[j-1] = 0.;}
					// }
					// else
					// {
						// dioj[i-1] = 0.;
					// }
				// }
				// if(dioj[i-1+5] < 1. && dioj[i-1+5] > 0. && flagb == 1){dioj[i-1+5] = 0.;}
				// if(dioj[i-1+5] > 1.)
				// {
					// flagb = 1;
					// if(dioj[i-1+5] > dbmax)
					// {
						// dbmax = dioj[i-1+5];
						// dioj[i-1+5] = 1.;
						// imax = i;
						// for(j=1;j<imax;j++){dioj[j-1+5] = 0.;}
					// }
					// else
					// {
						// dioj[i-1+5] = 0.;
					// }
				// }
				// if(dioj[i-1+10] < 1. && dioj[i-1+10] > 0. && flagc == 1){dioj[i-1+10] = 0.;}
				// if(dioj[i-1+10] > 1.)
				// {
					// flagc = 1;
					// if(dioj[i-1+10] > dcmax)
					// {
						// dcmax = dioj[i-1+10];
						// dioj[i-1+10] = 1.;
						// imax = i;
						// for(j=1;j<imax;j++){dioj[j-1+10] = 0.;}
					// }
					// else
					// {
						// dioj[i-1+10] = 0.;
					// }
				// }
			// }
		// }
		
		
		
		
		// i=0;
		// while((dioj[0]+dioj[1]+dioj[3]+dioj[4]) > 1.)
		// {
			// if(dioj[i] > 0.){dioj[i] = dioj[i] - 0.05;}
			// i++;
			// if(i>4){i=0;}
		// }
		
		// i=5;
		// while((dioj[5]+dioj[6]+dioj[8]+dioj[9]) > 1.)
		// {
			// if(dioj[i] > 0.){dioj[i] = dioj[i] - 0.05;}
			// i++;
			// if(i>9){i=5;}
		// }
		
		// i=10;
		// while((dioj[10]+dioj[11]+dioj[13]+dioj[14]) > 1.)
		// {
			// if(dioj[i] > 0.){dioj[i] = dioj[i] - 0.05;}
			// i++;
			// if(i>14){i=10;}
		// }
		
		// Calculo dio3
		
		dioj[2]=1.-dioj[0]-dioj[1]-dioj[3]-dioj[4];
		dioj[7]=1.-dioj[5]-dioj[6]-dioj[8]-dioj[9];
		dioj[12]=1.-dioj[10]-dioj[11]-dioj[13]-dioj[14];
		
		// Saturacion

		for(i=1;i<=15;i++)
		{
			if(dioj[i-1] < 0.){dioj[i-1] = 0.;}
			if(dioj[i-1] > 1.){dioj[i-1] = 1.;}
		}
		  for(i=1;i<=15;i++)
			{
     
      
                fprintf(MYFILEPR,"diof satured %d: %f  \n", i,dioj[i-1]);
      
			}
		
		
		// Representacion d
		
		for(i=1;i<=5;i++)
		{
			pda[i-1] = dioj[i-1];
			pdb[i-1] = dioj[i-1+5];
			pdc[i-1] = dioj[i-1+10];
		}
		
		psat[0] = dioj[0]+dioj[1]+dioj[2]+dioj[3]+dioj[4];
		psat[1] = dioj[0+5]+dioj[1+5]+dioj[2+5]+dioj[3+5]+dioj[4+5];
		psat[2] = dioj[10]+dioj[11]+dioj[12]+dioj[13]+dioj[14];
		
		
	}
	else
	{
		cont = cont + 1.;
	}
	
	if(cont == 100.){cont = 99.;}
	 
    
	//// Modulacion
	
	// RAMA A
	
	aux = 0;
	
	if(cont>=((dioj[0]/2. + dioj[1] + dioj[2] + dioj[3] + dioj[4])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sa1[0] = 0;
		Sa2[0] = 0;
		Sa3[0] = 0;
		Sa4[0] = 0;
    }
	if(cont>=((dioj[0]/2. + dioj[1]/2. + dioj[2] + dioj[3] + dioj[4])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sa1[0] = 0;
		Sa2[0] = 0;
		Sa3[0] = 0;
		Sa4[0] = 1;
    }
	if(cont>=((dioj[0]/2. + dioj[1]/2. + dioj[2]/2. + dioj[3] + dioj[4])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sa1[0] = 0;
		Sa2[0] = 0;
		Sa3[0] = 1;
		Sa4[0] = 1;
    }
	if(cont>=((dioj[0]/2. + dioj[1]/2. + dioj[2]/2. + dioj[3]/2. + dioj[4])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sa1[0] = 0;
		Sa2[0] = 1;
		Sa3[0] = 1;
		Sa4[0] = 1;
    }
	if(cont>=((dioj[0]/2. + dioj[1]/2. + dioj[2]/2. + dioj[3]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sa1[0] = 1;
		Sa2[0] = 1;
		Sa3[0] = 1;
		Sa4[0] = 1;
    }
    if(cont>=((dioj[0]/2. + dioj[1]/2. + dioj[2]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sa1[0] = 0;
		Sa2[0] = 1;
		Sa3[0] = 1;
		Sa4[0] = 1;
    }
    if(cont>=((dioj[0]/2. + dioj[1]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sa1[0] = 0;
		Sa2[0] = 0;
		Sa3[0] = 1;
		Sa4[0] = 1;
    } 
    if(cont>=((dioj[0]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sa1[0] = 0;
		Sa2[0] = 0;
		Sa3[0] = 0;
		Sa4[0] = 1;
    } 
    if (aux == 0)
    {
        Sa1[0] = 0;
		Sa2[0] = 0;
		Sa3[0] = 0;
		Sa4[0] = 0;
    } 
	
	// RAMA B
	
	aux = 0;
	
	if(cont>=((dioj[5]/2. + dioj[6] + dioj[7] + dioj[8] + dioj[9])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sb1[0] = 0;
		Sb2[0] = 0;
		Sb3[0] = 0;
		Sb4[0] = 0;
    }
	if(cont>=((dioj[5]/2. + dioj[6]/2. + dioj[7] + dioj[8] + dioj[9])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sb1[0] = 0;
		Sb2[0] = 0;
		Sb3[0] = 0;
		Sb4[0] = 1;
    }
	if(cont>=((dioj[5]/2. + dioj[6]/2. + dioj[7]/2. + dioj[8] + dioj[9])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sb1[0] = 0;
		Sb2[0] = 0;
		Sb3[0] = 1;
		Sb4[0] = 1;
    }
	if(cont>=((dioj[5]/2. + dioj[6]/2. + dioj[7]/2. + dioj[8]/2. + dioj[9])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sb1[0] = 0;
		Sb2[0] = 1;
		Sb3[0] = 1;
		Sb4[0] = 1;
    }
	if(cont>=((dioj[5]/2. + dioj[6]/2. + dioj[7]/2. + dioj[8]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sb1[0] = 1;
		Sb2[0] = 1;
		Sb3[0] = 1;
		Sb4[0] = 1;
    }
    if(cont>=((dioj[5]/2. + dioj[6]/2. + dioj[7]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sb1[0] = 0;
		Sb2[0] = 1;
		Sb3[0] = 1;
		Sb4[0] = 1;
    }
    if(cont>=((dioj[5]/2. + dioj[6]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sb1[0] = 0;
		Sb2[0] = 0;
		Sb3[0] = 1;
		Sb4[0] = 1;
    } 
    if(cont>=((dioj[5]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sb1[0] = 0;
		Sb2[0] = 0;
		Sb3[0] = 0;
		Sb4[0] = 1;
    } 
    if (aux == 0)
    {
        Sb1[0] = 0;
		Sb2[0] = 0;
		Sb3[0] = 0;
		Sb4[0] = 0;
    } 
	
	// RAMA C
	
	aux = 0;
	
	if(cont>=((dioj[10]/2. + dioj[11] + dioj[12] + dioj[13] + dioj[14])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sc1[0] = 0;
		Sc2[0] = 0;
		Sc3[0] = 0;
		Sc4[0] = 0;
    }
	if(cont>=((dioj[10]/2. + dioj[11]/2. + dioj[12] + dioj[13] + dioj[14])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sc1[0] = 0;
		Sc2[0] = 0;
		Sc3[0] = 0;
		Sc4[0] = 1;
    }
	if(cont>=((dioj[10]/2. + dioj[11]/2. + dioj[12]/2. + dioj[13] + dioj[14])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sc1[0] = 0;
		Sc2[0] = 0;
		Sc3[0] = 1;
		Sc4[0] = 1;
    }
	if(cont>=((dioj[10]/2. + dioj[11]/2. + dioj[12]/2. + dioj[13]/2. + dioj[14])*ciclos) && aux==0)
    {
        aux = 1;
        
        Sc1[0] = 0;
		Sc2[0] = 1;
		Sc3[0] = 1;
		Sc4[0] = 1;
    }
	if(cont>=((dioj[10]/2. + dioj[11]/2. + dioj[12]/2. + dioj[13]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sc1[0] = 1;
		Sc2[0] = 1;
		Sc3[0] = 1;
		Sc4[0] = 1;
    }
    if(cont>=((dioj[10]/2. + dioj[11]/2. + dioj[12]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sc1[0] = 0;
		Sc2[0] = 1;
		Sc3[0] = 1;
		Sc4[0] = 1;
    }
    if(cont>=((dioj[10]/2. + dioj[11]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sc1[0] = 0;
		Sc2[0] = 0;
		Sc3[0] = 1;
		Sc4[0] = 1;
    } 
    if(cont>=((dioj[10]/2.)*ciclos) && aux==0)
    {
        aux = 1;
        
        Sc1[0] = 0;
		Sc2[0] = 0;
		Sc3[0] = 0;
		Sc4[0] = 1;
    } 
    if (aux == 0)
    {
        Sc1[0] = 0;
		Sc2[0] = 0;
		Sc3[0] = 0;
		Sc4[0] = 0;
    } 
	
	pvd1[0] = vd1;
	pvd2[0] = vd2;
	pvd3[0] = vd3;
	ppot[0] = pref;
	pu[0] = p;
	
	// Calculo derevidas del error
		
		dvd1 = -u[2]*i_aBY[0]-u[3]*i_aBY[1];
		dvd2 = -u[4]*i_aBY[0]-u[5]*i_aBY[1];
		dvd3 = -u[6]*i_aBY[0]-u[7]*i_aBY[1];
		
		// Signo vd
		
		if(vd1 > 0.)
			sign1 = 1.;
		else
			sign1 = -1.;
		
		if(vd2 > 0.)
			sign2 = 1.;
		else
			sign2 = -1.;
		
		if(vd3 > 0.)
			sign3 = 1.;
		else
			sign3 = -1.;
		
		// Termino dvd/dt*sign(vd)
		
		term1 = dvd1*sign1;
		term2 = dvd2*sign2;
		term3 = dvd3*sign3;
		
		pterm[0] = dvd3;
		// pterm[0] = term1;
		pterm[1] = term2;
		pterm[2] = term3;
		
		// pterm[0] = i_aBY[0];
		// pterm[1] = i_aBY[1];
		// pterm[2] = i_aBY[2];
	
    
    
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
