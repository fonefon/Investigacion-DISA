/************* Archivos de cabecera ****************/
#include "math.h"
#include "stdio.h"
#include "funciones_extra.c"
#include "SVM_modificado.c"
#include "SVM_modificado_2.c"
#include "DSC_2_4_8_16_32.c"
/************* Para nombrar como PSCAD (Fortran convention) ****************/
typedef double real_T;
typedef long int integer;
typedef long int logical;
/************* Zona de las etiquetas ****************/

/* Constantes */
#define PI		3.14159265
#define raiz1_3 0.57735027
#define raiz1_2 0.70710678
#define raiz1_6 0.40824829
#define raiz2_3 0.81649658
#define raiz2	1.41421356
#define raiz6	2.44948974

/* Filtrados */
#define f_corte	600.0 			// Esta frecuencia se encuentra dada en rad/s
#define w1	2.0*PI*1500.0 	// Esta es la frecuencia de corte del filtro de tensiones de condensadores
#define wcond	2.0*PI*1000.0 // Frecuencia de filtrado de tension suma de condensadores
#define wr	2*PI*50
#define Q	0.
#define Ts 0.0001

/* Controladores */
//#define kis	0.15//1.5//2.5//2.5 //0.05
//#define kps_Pref 0.05//0.01//0.08//0.5 //0.02

#define C_PRECARGA 		 0
#define C_RED			 1
#define C_PRE_TRAFO		 2
#define C_TRAFO_RED		 3
#define C_DC_SOURCE     16 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define C_IN_EL_15		 5 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define C_OUT_EL_15		 6 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define C_BOBINA_CARGA	14 //Abre y cierra el circuito que conecta el filtro de bobinas con la carga ///////////////////////////////////////////
#define C_CARGA_ESTRELL	12 //Abre y cierra el circuito que pone tres resistencias en estrella o en paralelo respectivamente ////////////////////
#define C_CHOPPER_CARGA	10 //Abre y cierra el circuito que conecta la etapa de chopper y una resistencia con las tres resistencias restantes////

/* Modo de funcionamiento */

#define RECTIFICADOR 		1  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define INVERSOR			2  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*************************/
/* Párametros de Control */
/*************************/
				
//	k_balancing 0.000075 			Proporcional balanceo de condensadores
//	kps_Pref  						Proporcional lazo externo de tensión 
//	kis y kis1 (kis1 = kis*Ts)		Integral lazo externo de tensión
//	Ki 								Integral control de corriente
// 	Kp 								Proporcional control de corriente
//void c3ndpc_(real_T k_balancing, real_T kps_Pref, real_T kis, real_T Ki, real_T Kp, real_T* A_max, real_T* V_max, real_T* A_RMS_max, real_T* V_RMS_max, real_T* Vcr_max,int inicio, int Enable, real_T* iar,real_T* ibr,real_T* icr,real_T* vc1r,real_T* vc2r,real_T* vc3r,real_T* vc4r,real_T* Vdcref,real_T* var,real_T* vbr,real_T* vcr,real_T* tiempos_a,real_T* tiempos_b,real_T* tiempos_c,real_T* Debug, real_T* Perifericos,real_T *DELT,real_T *TIME)
//{
void c3ndpc_(real_T k_balancing, real_T kps_Pref, real_T kis, real_T Kp, real_T Ki, real_T ki_balancing, real_T** Entradas, real_T* tiempos_a, real_T* tiempos_b, real_T* tiempos_c,real_T* Debug, real_T* Perifericos, real_T* DELT)
{
	/************* Definición de variables ***************/

	
/*	//Parámetros de Control
	static float k_balancing;	//	Proporcional balanceo de tensiones
	static float kps_Pref;		//	Proporcional lazo externo de tensión
	static float kis;			//	Integral lazo externo de tensión
	static float Ki;			//	Integral control de corriente
	static float Kp;			// 	Proporcional control de corriente	*/
	//static real_T kis1;			//	Integral k1 multiplicada por Ts
	
	//Variables PLL
	static float theta, P_pll, suma_P_pll, PI_pll, PI_pll_ant;
	static float theta_ant_pll, w, w_ant;

	//Variables errores de tensión
	static float vd1k0 = 0.;
	static float vd2k0 = 0.;
	static float vd3k0 = 0.;
		
	// Variables del integrador para p y q (ORS):
	static float perrork1;
	static float qerrork1;
	static float kipq_discr;
	static float IntPk0;
	static float IntPk1;
	static float IntQk0;
	static float IntQk1;
	
	//ERRORES
	static float error;
	
	 //Variables del PI de control de la diferencia de tensiones:
	static float Int_err_xdif;
	static float vd1k0_1;
	static float kbal1;
	static float kbal2;
	static float PI_err_xdif;
	static int caso;
	static float da, db, dc;
	static float alpha1, alpha2, alpha3, alphamin, alphamax, xmax, xmin;
	static float u_amax, u_amin, u_bmax, u_bmin, u_cmax, u_cmin;
	static float f1, f2, f3, f4, f5;
	static float yx0, zx0, xy0, zy0, xz0, yz0;
	static float vec_min[5], vec_min_ant[5],vec_min_PB[5], vec_min_PB_ant[5], Banda, Vmin, x, y, z;
    static float Vmin_ant, epsilon, aux;
	static int jmin, imin, imin_ant;

	// Señales de control:
	static float da1;
	static float da2;
	static float da3;
	static float da4;
	static float db1;
	static float db2;
	static float db3;
	static float db4;
	static float dc1;
	static float dc2;
	static float dc3;
	static float dc4;
	static float disp1[6];
	float* disp = &disp1;
    static float dalfa1;
	static float dalfa2;
	static float dalfa3;
	static float dalfa4;
	static float dbeta1;
	static float dbeta2;
	static float dbeta3;
	static float dbeta4;
	static float u1;
	static float u2;
	static float u3;
	static float u4;
	static float u5;
	static float u6;
	static float u7;
	static float u8;
	static float dgamma1,dgamma2,dgamma3,dgamma4;

	// Tiempos para la modulación:
	static float t1a,t2a,t3a,t4a,t5a;
    static float t1b,t2b,t3b,t4b,t5b;
	static float t1c,t2c,t3c,t4c,t5c;

    /* Coeficientes de filtrado */
	static float c11,c12,c13,c31,c32,c33;
	
	// Filtrado de perror y q error
	static float perror_1, qerror_1, Pdif, Pdif_1, Qdif, Qdif_1;
		
	/* PI de suma de tensiones */
	static float error_Vdc,error_Vdc1,Ierror_Vdc, Perror_Vdc, Perror_Vdc1;
	static float Pcond, Pcond1;
	
	/* PI de diferencia de tensiones */
	static float signog=1.;
	
	/* Tensiones de condensadores y filtrados */
	static float lvc1,lvc1_1,lvc2,lvc2_1,lvc1_f,lvc1_f1,lvc2_f,lvc2_f1;
    //=====================================================================================================================
	// FUJ INICIO 3
	static float lvc3,lvc3_1,lvc4,lvc4_1,lvc3_f,lvc3_f1,lvc4_f,lvc4_f1;
    // FUJ FINAL 3
    // =====================================================================================================================
	static float Vcond,Vcond1,Vcond_f,Vcond_f1,Inv_Vcond;
	static float lvdcref = 490;
	
	/* Tensiones de red y filtradas */
	static float Vr,Vs,Vt;
	static float valfa,vbeta;
	static float Vd, Vq;
	
	/* Corrientes de red y filtradas */
	static float id, iq;
	static float Iu,Iv,Iw;
	static float ialfa_c=0,ialfa_c1=0,ialfa_c2=0,ibeta_c=0,ibeta_c1=0,ibeta_c2=0;
	static float ialfa_cf=0,ialfa_cf1=0,ialfa_cf2=0,ibeta_cf=0,ibeta_cf1=0,ibeta_cf2=0;
	static float k,nr0,nr1,nr2,dr0,dr1,dr2,inv_dr0;
	
	/******* VARIABLES PARA CALCULO RMS *********/
	static int n=0;//Para contabilizar las veces que entra en el bucle RMS
	static float RMS_Ia = 0;//Para mostrar el resultado final
	static float RMS_Ib = 0;
	static float RMS_Ic = 0;
	static float RMS_Va = 0;
	static float RMS_Vb = 0;
	static float RMS_Vc = 0;
	static float V_rms_Ia = 0;//Para ir almacenando las variables
	static float V_rms_Ib = 0;
	static float V_rms_Ic = 0;
	static float V_rms_Va = 0;
	static float V_rms_Vb = 0;
	static float V_rms_Vc = 0;
	/********************************************/
	
	// DSC for PLL Pre-stage 
	extern float wg, wgf, wgf_ant, wg_ant;
	extern float wcut_lpf, T_DSC;
	extern int n_2, n_4, n_8, n_16, n_32;
	extern int n_2_1, n_4_1, n_8_1, n_16_1, n_32_1;
	extern float val_DSC2_in[120], vbe_DSC2_in[120], val_DSC2_out, vbe_DSC2_out;
	extern float val_DSC4_in[60], val_DSC4_out, vbe_DSC4_in[60], vbe_DSC4_out;
	extern float val_DSC8_in[30], val_DSC8_out, vbe_DSC8_in[30], vbe_DSC8_out;
	extern float val_DSC16_in[15], val_DSC16_out, vbe_DSC16_in[15], vbe_DSC16_out;
	extern float val_DSC32_in[8], val_DSC32_out, vbe_DSC32_in[8], vbe_DSC32_out;
	extern float interp_2, interp_4, interp_8, interp_16, interp_32;
	extern float sampleAl_2, sampleAl_4, sampleAl_8, sampleAl_16, sampleAl_32;
	extern float sampleBe_2, sampleBe_4, sampleBe_8, sampleBe_16, sampleBe_32;
	
	//Variablas para el control DPC
	
	//Potencias instantaneas
	static float p;//Potencia activa instantanea en W
	static float q;//Potencia reactiva instantanea en VA
	static float pref = 0.;//Referencia de la potencia activa instantanea
	static float qref = 0.;//Referencia de la potencia reactiva instantanea
	static float perror = 0.;//Error en la potencia activa instantanea
	static float qerror = 0.;//Error en la potencia reactiva instantanea
	static float X_estimado = 0.;
	static float intx = 0.;
	
	//Variables para la estimación de parametros
	
	static float Xest = 0.;//2*M_PI*50*0.0008;
			
	//Variables para el cálculo de las ORS (DPC)
	
	struct VectorAB
	{
		int indiceAB;
		float Calpha;
		float Cbeta;
	};
	
	static float ModVred;
	static float ModAB1 = 0.;
	static float ModAB2 = 0.;
	static struct VectorAB ORS1AB;
	static struct VectorAB ORS2AB;
	static struct VectorAB ORS;
	
	static float valphaP = 0.;
	static float vbetaP = 0.;
	static float K1;
	static float K2;
	static float Kq;
	static float RefAlpha;
	static float RefBeta;
	static float Refd;
	static float Refq;
    //******//
	
	// Variables para control en corriente
    static float ia_ref;
    static float ib_ref;
	static float id_ref;
	static float iq_ref;
	static float iar_ref, ibr_ref, icr_ref;
    static float Total;
	
	//*******//
	static int initialized = 0; //Variable para inicialización
	int i;
	static int estado = 1; //Máquina de estado del sistema (0: ERROR; 1: STANDBY; 2: PRECARGA; 3: MARCHA)
	static int temp1 = 0; //Temporizador 1
	static int temp2 = 0; //Temporizador 2
	static int cercania;
	// para testeo
	static int iteracion = 0;
	static float señales1[6] = {0.,0.,0.3,0.4,0.6,0.2};//t5c
	static float señales2[6] = {0.3,0.,0.2,0.05,0,0.3};//t1c
	static float señales3[6] = {0.3,0.6,0.1,0.15,0,0.1};//t3c
	
	/**************************/
	/* ENTRADAS DE LA FUNCIÓN */
	/**************************/
	float A_max1 = 0;		
	float V_max1 = 0;		
	float A_RMS_max1 = 0;	
	float V_RMS_max1 = 0;	
	float Vcr_max1 = 0;		
	int inicio1 = 0;			
	int Enable1 = 0;
	int rearme1 = 0;
	float e_drv1 = 0;
	float iar1 = 0; 			
	float ibr1 = 0;			
	float icr1 = 0;			
	float vc1r1 = 0;			
	float vc2r1 = 0;			
	float vc3r1 = 0;			
	float vc4r1 = 0;			
	float Vdcref1 = 0;
	float enable_desbalance1 = 0;
	static float bandera_desb = 0;

	int enable_disp = 0;
	int go1 = 0;
	float K6h1 = 0.;
	float ki_balance1 = 0.;
	float Qref1;
	
	float Vrs1 = 0.;
	float Vst1 = 0.;
	float Vtr1 = 0.;
	float var1 = 0;		
	float vbr1 = 0;	
	float vcr1 = 0;	
	
	int modo1 = 0; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/*** PARA EL RESONANTE ***/
	float hwTs, t1PR_d, t2PR_d, t1PR_q, t2PR_q, dc_d, dc_q, error_d, error_q;
	static float error_d_ant_1, error_d_ant_2, error_q_ant_1, error_q_ant_2, dc_ant_d_2, dc_ant_d_1, dc_ant_q_1, dc_ant_q_2;
	float periodo;
	static float theta_ant;
	float t1PR_a, t1PR_b, t2PR_a, t2PR_b, dc_a, dc_b, error_a, error_b;
	static float error_a_ant_1, error_a_ant_2, error_b_ant_1, error_b_ant_2, dc_ant_a_1, dc_ant_a_2, dc_ant_b_1, dc_ant_b_2;
	float a1, b0, b1, b2;
	static float nume1, nume2, nume3, deno1, KSOGI;
	static float in1_LP[4], in2_LP[4], out1_LP[4], out2_LP[4];
	static float in1_LP_1[4], in2_LP_1[4], out1_LP_1[4], out2_LP_1[4];
	static float in1_LP_2[4], in2_LP_2[4], out1_LP_2[4], out2_LP_2[4];
	static float in3_LP[4], in4_LP[4], out3_LP[4], out4_LP[4];
	static float in3_LP_1[4], in4_LP_1[4], out3_LP_1[4], out4_LP_1[4];
	static float in3_LP_2[4], in4_LP_2[4], out3_LP_2[4], out4_LP_2[4];

	static float in1_LP2[4], in2_LP2[4], out1_LP2[4], out2_LP2[4];
	static float in1_LP2_1[4], in2_LP2_1[4], out1_LP2_1[4], out2_LP2_1[4];
	static float in1_LP2_2[4], in2_LP2_2[4], out1_LP2_2[4], out2_LP2_2[4];	
	static float in3_LP2[4], in4_LP2[4], out3_LP2[4], out4_LP2[4];
	static float in3_LP2_1[4], in4_LP2_1[4], out3_LP2_1[4], out4_LP2_1[4];
	static float in3_LP2_2[4], in4_LP2_2[4], out3_LP2_2[4], out4_LP2_2[4];	
	
	/*** PASO BAJO CORRIENTES ***/
	static float ial,ial_1,ibl,ibl_1,ialfa_c_1, ibeta_c_1;
	
	/***  RESONANTE DE VALFA-VBETA ***/
	static float x1alfa, x1beta, x2alfa_p, x2beta_p, x2alfa_n, x2beta_n, x3alfa_p, x3beta_p;
	static float x3alfa_n, x3beta_n, x4alfa_p, x4beta_p, x4alfa_n, x4beta_n;
	static float x2alfa_p_1, x2alfa_n_1, x2beta_p_1, x2beta_n_1;
	float valfa_res, vbeta_res;
	static float valfa_res_1, vbeta_res_1, gamma;
		
	/*** RESONANTE DE IALFA-IBETA ***/
	static float i1alfa, i1beta, i2alfa_p, i2beta_p, i2alfa_n, i2beta_n, i3alfa_p, i3beta_p;
	static float i3alfa_n, i3beta_n, i4alfa_p, i4beta_p, i4alfa_n, i4beta_n;
	static float i2alfa_p_1, i2alfa_n_1, i2beta_p_1, i2beta_n_1;
	float ialfa_res, ibeta_res;
	static float ialfa_res_1, ibeta_res_1;
	static int temp3;
	
	/*** CONTROL PR AL-BE ***/
	static float wres[3], wcut, ki_res, kt, kp_res;
	static float n0[3], n1[3], n2[3], den1, den2;
	static float in1, in2, out1[3], out2[3];
	static float in1_1[3], in2_1[3], out1_1[3], out2_1[3];
	static float in1_2[3], in2_2[3], out1_2[3], out2_2[3];
	
	/*** HARMONIC DETECTION ***/
	static float idh_p[4],idh_n[4],iqh_p[4],iqh_n[4];
	static float ialh_p[4],ibeh_p[4],ialh_n[4],ibeh_n[4],ialh[4],ibeh[4];
	static float Intdh_p[4],Intdh_n[4],Intqh_p[4],Intqh_n[4];
	static float PIdh_p[4],PIdh_n[4],PIqh_p[4],PIqh_n[4];
	
	float* A_max = &A_max1;		// 14
	float* V_max = &V_max1;		// 15
	float* A_RMS_max = &A_RMS_max1;	// 16
	float* V_RMS_max = &V_RMS_max1;	// 17
	float* Vcr_max = &Vcr_max1;		// 18
	float* rearme = &rearme1;		// 19
	int* inicio = &inicio1;			// 12
	int* Enable = &Enable1;			// 11
	float* e_drv = &e_drv1;			// 13
	float* iar = &iar1; 			// 0
	float* ibr = &ibr1;			// 1
	float* icr = &icr1;			// 2
	float* vc1r = &vc1r1;			// 3
	float* vc2r = &vc2r1;			// 4
	float* vc3r = &vc3r1;			// 5
	float* vc4r = &vc4r1;			// 6
	float* Vdcref = &Vdcref1;		// 7
	float* Vrs = &Vrs1;				// 8
	float* Vst = &Vst1;				// 9
	float* Vtr = &Vtr1;				// 10
	float* go = &go1;		// 20
	float* ki_balance = &ki_balance1;	// 21
	float* K6h = &K6h1;	// 22	
	float* Qref = &Qref1;	//23
	float* enable_desbalance = &enable_desbalance1; //24
	int* modo = &modo1; //25   //////////////////////////////////////////////////////////////////////////////////////////////
	
	float* var = &var1;		
	float* vbr = &vbr1;		
	float* vcr = &vcr1;		
	/*
	float tiempos_a //Salida
	float tiempos_b //Salida
	float tiempos_c // Salida
	float Debug; //Salida
	Perifericos; //Salida*/
	
	// Para el cambio de funcionamiento
	
	
	

	
	//Ts=*DELT;//*DELT = 1e-4 (pasos de simulacion) -> 10 KHz
	
		/*****************/
		/* LEER ENTRADAS */
		/*****************/
		
		//void leer_entradas(real_T** Entradas, float* A_max, float* V_max, float* A_RMS_max, float* V_RMS_max, float* Vcr_max, int* inicio, int* Enable, float* iar, float* ibr, float* icr, float* vc1r, float* vc2r, float* vc3r, float* vc4r, float* Vdcref, float* var, float* vbr, float* vcr);
		*iar = *Entradas[0];
		*ibr = *Entradas[1];
		*icr = *Entradas[2];
		*vc1r = *Entradas[3];
		*vc2r = *Entradas[4];
		*vc3r = *Entradas[5];
		*vc4r = *Entradas[6];
		*Vdcref = *Entradas[7];
		*Vrs = *Entradas[8]; // Vrs
		*Vst = *Entradas[9]; // Vst
		*Vtr = *Entradas[10];// Vtr
		*Enable = *Entradas[11];
		*inicio = *Entradas[12];
		*e_drv 	= *Entradas[13];
		*A_max = *Entradas[14];
		*V_max = *Entradas[15];
		*A_RMS_max = *Entradas[16];
		*V_RMS_max = *Entradas[17];
		*Vcr_max = *Entradas[18];
		*rearme = *Entradas[19];
		*go = *Entradas[20];
		*ki_balance = *Entradas[21];
		*K6h = *Entradas[22];		
		*Qref = *Entradas[23];
		*enable_desbalance = *Entradas[24];
		*modo = *Entradas[25]; //////////////////////////////////////////////////////////////////////////////////
		
		
		*var = (*Vrs - *Vtr)/3; // Vr
		*vbr = (*Vst - *Vrs)/3; // Vs
		*vcr = (*Vtr - *Vst)/3;	// Vt
		
		Vr = (*Vrs - *Vtr)/3; // Vr
		Vs = (*Vst - *Vrs)/3; // Vs
		Vt = (*Vtr - *Vst)/3; // Vt
		
		//Vs = -Vt - Vr;
		
		
		/***********************/
		/* Calculo valores RMS */
		/***********************/
			
		n++; //Suma 1 al bucle RMS
		if (n >= 1./(Ts*50))//¿Se ha almacenado ya un periodo completo?
			{
				//CALCULAR Y MOSTRAR VALOR RMS
				RMS_Ia = sqrt(V_rms_Ia/n);
				RMS_Ib = sqrt(V_rms_Ib/n);
				RMS_Ic = sqrt(V_rms_Ic/n);
				RMS_Va = sqrt(V_rms_Va/n);
				RMS_Vb = sqrt(V_rms_Vb/n);
				RMS_Vc = sqrt(V_rms_Vc/n);
					
				n=1;
					
				V_rms_Ia = 0;
				V_rms_Ib = 0;
				V_rms_Ic = 0;
				V_rms_Va = 0;
				V_rms_Vb = 0;
				V_rms_Vc = 0;
			}
		//Ir Almacenando variables
		V_rms_Ia += (*iar)*(*iar);
		V_rms_Ib += (*ibr)*(*ibr);
		V_rms_Ic += (*icr)*(*icr);
		V_rms_Va += (*var)*(*var);
		V_rms_Vb += (*vbr)*(*vbr);
		V_rms_Vc += (*vcr)*(*vcr);
		
		/*********************/
 		/* Corrientes de red */
		/*********************/
		Iu = *iar;
		Iv = *ibr;
		Iw = *icr;
		
		/******************************/
 		/* Calculo valfa vbeta de red */
		/******************************/
			
		valfa = raiz2_3*Vr - raiz1_6*(Vs+Vt);
		vbeta = raiz1_2*(Vs-Vt);
		
		/******************************/
 		/* Calculo ialfa ibeta de red */
		/******************************/
			
		ialfa_c = raiz2_3*Iu - raiz1_6*(Iv+Iw);
		ibeta_c = raiz1_2*(Iv-Iw);
		
		/*******/
		/* DSC */
		/*******/
		
		DSC_2_4_8_16_32(valfa, vbeta); // Save valfa and vbeta fundamental component into val_DSC32_out and vbe_DSC32_out
			//	valfa = val_DSC32_out
			//	vbeta = vbe_DSC32_out
			
		/*********************/
 		/* TRANSFORMACION DQ */
		/*********************/
		/** Id y Iq **/
			
		id = ialfa_c * cos(theta) + ibeta_c * sin(theta);
		iq = -ialfa_c * sin(theta) + ibeta_c * cos(theta);		

		/** Vd y Vq **/
			
		//Vd = valfa * cos(theta) + vbeta * sin(theta);
		Vd = val_DSC32_out * cos(theta) + vbe_DSC32_out * sin(theta);
		//Vq = -valfa * sin(theta) + vbeta * cos(theta);		
		Vq = -val_DSC32_out * sin(theta) + vbe_DSC32_out * cos(theta);		

		/*******/		
		/* PLL */
		/*******/
		
		suma_P_pll += Vq/sqrt(Vd*Vd + Vq*Vq)/10000;
		
		if (Vq/sqrt(Vd*Vd + Vq*Vq) > 20 || Vq/sqrt(Vd*Vd + Vq*Vq) < -20)	// Anti wind-up
		{	suma_P_pll = 0.;
			Vq = signo(Vq)*20.*sqrt(Vd*Vd + Vq*Vq);
		}
		
		wg = 5000*Vq/sqrt(Vd*Vd + Vq*Vq) + 100000*suma_P_pll;
		
		theta += wg/10000;
		
		if (theta > 2.*pi)
			theta -= 2.*pi;
		if (theta < 0)
			theta += 2.*pi;
		
		/****************************************/
		/* Filtrado de las corrientes alfa-beta */
		/****************************************/
		
		i1alfa = ialfa_c - ialfa_res_1;
		i1beta = ibeta_c - ibeta_res_1;
			
		i2alfa_p = i1alfa*cos(theta) + i1beta*sin(theta);
		i2beta_p = -i1alfa*sin(theta) + i1beta*cos(theta);
		i2alfa_n = i1alfa*cos(theta) - i1beta*sin(theta);
		i2beta_n = i1alfa*sin(theta) + i1beta*cos(theta);
			
		gamma = 1.;
			
		i3alfa_p += Ts/(2*gamma)*(i2alfa_p + i2alfa_p_1);
		i3beta_p += Ts/(2*gamma)*(i2beta_p + i2beta_p_1);
		i3alfa_n += Ts/(2*gamma)*(i2alfa_n + i2alfa_n_1);
		i3beta_n += Ts/(2*gamma)*(i2beta_n + i2beta_n_1);
			
		i4alfa_p = i3alfa_p*cos(theta) - i3beta_p*sin(theta);//ialfa_p
		i4beta_p = i3alfa_p*sin(theta) + i3beta_p*cos(theta);//ibeta_p
		i4alfa_n = i3alfa_n*cos(theta) + i3beta_n*sin(theta);
		i4beta_n = -i3alfa_n*sin(theta) + i3beta_n*cos(theta);
			
		ialfa_res = i4alfa_p + i4alfa_n;
		ibeta_res = i4beta_p + i4beta_n;
		
		ialfa_res_1 = ialfa_res;
		ibeta_res_1 = ibeta_res;
			
		i2alfa_p_1 = i2alfa_p;
		i2beta_p_1 = i2beta_p;
		i2alfa_n_1 = i2alfa_n;
		i2beta_n_1 = i2beta_n;
					
		/***********************/
		/* Calculo de potencia */
		/***********************/
		
		p = (valfa*ialfa_c + vbeta*ibeta_c);
		q = (valfa*ibeta_c - vbeta*ialfa_c);
			
		qref = *Qref;
		
		/*********************/
		/* COMPROBAR ERRORES */
		/*********************/
		
	if (estado != 0) //Si está en estado de error, comprueba los errores dentro del estado
	{	
		//Check error devuelve cero en caso de no error
		error = Check_errors(*A_max, *V_max, *iar, *ibr, *icr, *var, *vbr, *vcr); //La variable error devuelve el tipo de error
		if( error != 0)//¿Algún error?									1,2,3 SOBRECORRIENTE instantanea FASE A,B,C respectivamente
		//																4,5,6 SOBRETENSION instantanea FASE A,B,C respectivamente
			estado = 0;
		
		if (error == 0)// Para no pisar el error anterior
		{
			error = Check_errors(*A_RMS_max, *V_RMS_max,RMS_Ia, RMS_Ib,RMS_Ic,RMS_Va,RMS_Vb,RMS_Vc)+7.;
			if( error != 7)//¿Algún error?									8,9,10 		SOBRECORRIENTE en RMS FASE A,B,C respectivamente
			//																11,12,13 	SOBRETENSIÓN en RMS FASE A,B,C respectivamente
				estado = 0;
		}
		
		if (error == 7)
		{
			error = Check_error_vdc(*Vcr_max, *vc1r, *vc2r, *vc3r, *vc4r)+14.;
			if (error != 14)//												15,16,17,18		SOBRETENSIÓN CONDENSADOR VC1,VC2,VC3,VC4 respectivamente
				estado = 0;
		}
		if (*e_drv != 0)// Si hay error de driver
		{
			estado = 0;
			error = 19;
		}
	}	
	if(initialized == 0) //Si initialized = 0 inicializa todas las variables (necesario realizar al menos una vez)
	{
		/* Usamos esta variable lógica para que estos  */
		/* cálculos sólo se realicen una vez y acelere */
		/* el proceso de simulación.  	               */
		
        //Ts=*DELT;//*DELT = 1e-4 (pasos de simulacion) -> 10 KHz
		initialized = 1;
		estado = 1;
		
		Unset_contactor(Perifericos,C_PRECARGA);
		Unset_contactor(Perifericos,C_RED);
		Unset_contactor(Perifericos,C_PRE_TRAFO);
		Unset_contactor(Perifericos,C_TRAFO_RED);
		Unset_contactor(Perifericos,C_DC_SOURCE);     ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_IN_EL_15);      ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_OUT_EL_15);     ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_BOBINA_CARGA);  ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_CARGA_ESTRELL); ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_CHOPPER_CARGA); ////////////////////////////////////////////////////////////////////////////////////////////
		
		lvdcref = 490;
		
		// Iniciación de los tiempos de modulación
		t1a=0;
		t2a=0;
		t3a=0;
		t4a=0;
		t5a=0;
		t1b=0;
		t2b=0;
		t3b=0;
		t4b=0;
		t5b=0;
		t1c=0;
		t2c=0;
		t3c=0;
		t4c=0;
		t5c=0;

		// Iniciación de las variables y parámetros del integrador del ORS
	    perrork1 = 0.;
		qerrork1 = 0;
		IntPk0 = 0.;
		IntPk1 = 0.;
		IntQk0 = 0.;
		IntQk1 = 0.;
		
		error_d_ant_1 = 0.;
		error_q_ant_1 = 0.;
		error_d_ant_2 = 0.;
		error_q_ant_2 = 0.;
		
		dc_ant_d_1 = 0.;
		dc_ant_d_2 = 0.;
		dc_ant_q_1 = 0.;
		dc_ant_q_2 = 0.;
		
		hwTs = 2*PI*50*Ts;
		
		// Iniciación señales de control:
		da1 = 0.;
		da2 = 0.;
		da3 = 0.;
		da4 = 0.;
		db1 = 0.;
		db2 = 0.;
		db3 = 0.;
		db4 = 0.;
		dc1 = 0.;
		dc2 = 0.;
		dc3 = 0.;
		dc4 = 0.;
		dalfa1 = 0.;
		dalfa2 = 0.;
		dalfa3 = 0.;
		dalfa4 = 0.;
		dbeta1 = 0.;
		dbeta2 = 0.;
		dbeta3 = 0.;
		dbeta4 = 0.;
		u1 = 0;
		u2 = 0;
		u3 = 0.;
		u4 = 0.;
		u5 = 0.;
		u6 = 0.;
		u7 = 0.;
		u8 = 0.;

		//Inicialización de parámetros de controladores
		
		//kis1 = kis;//*0.0001;
		
		Ierror_Vdc = 0.;
		error_Vdc1 = 0.;
		Perror_Vdc1 = 0.;
		
		/**********************************************************/
		/* Coeficientes de la versión discretizada de los filtros */
		/**********************************************************/
		
		/* Pasabajas del filtro de tensión de condensadores (cada una)*/
		c11 = w1*Ts;
		c12 = 2.0 - w1*Ts;
		c13 = 1.0/(2.0 + w1*Ts);
		/* Pasabajas del filtro de tensión de condensadores (suma)*/
		c31 = wcond*Ts;
		c32 = 2.0 - wcond*Ts;
		c33 = 1.0/(2.0 + wcond*Ts);
		
		/* Resonante para las corrientes */
		k = 2./(Ts);
		nr0 = k*wr;
		nr1 = 0;
		nr2 = -k*wr;
		dr0 = k*k*Q + k*wr + Q*wr*wr;
		dr1 = 2.*(Q*wr*wr - k*k*Q);
		dr2 = k*k*Q - k*wr + Q*wr*wr;
		inv_dr0 = 1./dr0;
		
	}
	
	/*******************/
	/* ESTADO DE ERROR */
	/*******************/
	
	if (estado == 0)
	{
		temp1 = 0;
		temp2 = 0;
		Unset_contactor(Perifericos,C_PRECARGA);
		Unset_contactor(Perifericos,C_RED);
		Unset_contactor(Perifericos,C_PRE_TRAFO);
		Unset_contactor(Perifericos,C_TRAFO_RED);
		Unset_contactor(Perifericos,C_DC_SOURCE);     ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_IN_EL_15);      ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_OUT_EL_15);     ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_BOBINA_CARGA);  ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_CARGA_ESTRELL); ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_CHOPPER_CARGA); ////////////////////////////////////////////////////////////////////////////////////////////
		
		//enable_disp = 0;
		
				
		/** PONER A 0 TODOS LOS DISPAROS **/
		for (i=0; i<5; i++)
		{
			*(tiempos_a+i) = 0.0;
			*(tiempos_b+i) = 0.0;
			*(tiempos_c+i) = 0.0;
		}
		/*********************/
		/* COMPROBAR ERRORES */
		/*********************/
		if (*rearme == 1)
		{
			estado = 1;
			error = Check_errors(*A_max, *V_max, *iar, *ibr, *icr, *var, *vbr, *vcr); //La variable error devuelve el tipo de error
			if( error != 0)//¿Algún error?									1,2,3 SOBRECORRIENTE instantanea FASE A,B,C respectivamente
			//																4,5,6 SOBRETENSION instantanea FASE A,B,C respectivamente
				estado = 0;
			if (error == 0)// Para no pisar el error anterior
			{
				error = Check_errors(*A_RMS_max, *V_RMS_max,RMS_Ia, RMS_Ib,RMS_Ic,RMS_Va,RMS_Vb,RMS_Vc)+7.;
				if( error != 7)//¿Algún error?									8,9,10 		SOBRECORRIENTE en RMS FASE A,B,C respectivamente
				//																11,12,13 	SOBRETENSIÓN en RMS FASE A,B,C respectivamente
				estado = 0;
			}
			if (error == 7)
			{
				error = Check_error_vdc(*Vcr_max, *vc1r, *vc2r, *vc3r, *vc4r)+14.;
				if (error != 14)//												15,16,17,18		SOBRETENSIÓN CONDENSADOR VC1,VC2,VC3,VC4 respectivamente
					estado = 0;
			}
			if (*e_drv != 0)// Si hay error de driver
			{
				estado = 0;
				error = 19;
			}
		}
	}
	
	/*********************/
	/* ESTADO DE STANDBY */
	/*********************/
	else if(estado == 1)  		// Reseteo inicial de algunas variables (Cuando sistema no habilitado o en estado de Standby)
	{   
		enable_disp = 0;
		estado = 1;
		/*******************************/
		/* Inicialización de variables */
		/*******************************/
		/* Espera habilitación del control */
		
		/*****	DESACTIVAR CONTACTORES *****/
		Unset_contactor(Perifericos,C_PRECARGA);
		Unset_contactor(Perifericos,C_RED);
		Unset_contactor(Perifericos,C_PRE_TRAFO);
		Unset_contactor(Perifericos,C_TRAFO_RED);
		Unset_contactor(Perifericos,C_DC_SOURCE);     ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_IN_EL_15);      ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_OUT_EL_15);     ////////////////////////////////////////////////////////////////////////////////////////////
		Unset_contactor(Perifericos,C_BOBINA_CARGA);  ////////////////////////////////////////////////////////////////////////////////////////////
		Set_contactor(Perifericos,C_CARGA_ESTRELL); ////////////////////////////////////////////////////////////////////////////////////////////
		Set_contactor(Perifericos,C_CHOPPER_CARGA); ////////////////////////////////////////////////////////////////////////////////////////////
		
		temp1 = 0;
		temp2 = 0;
		
		lvdcref = 490;
		
		/***** PONER A 0 TODOS LOS PARÁMETROS INTEGRALES DEL CONTROL *****/
		Pcond = 0;
		Ierror_Vdc = 0.;
		
		IntPk0 = 0.;
		IntQk0 = 0.;
		
		IntPk1 = 0.;
		IntQk1 = 0.;
		
		error_d_ant_1 = 0.;
		error_d_ant_2 = 0.;
		error_q_ant_1 = 0.;
		error_q_ant_2 = 0.;
		
		dc_ant_d_1 = 0.;
		dc_ant_d_2 = 0.;
		dc_ant_q_1 = 0.;
		dc_ant_q_2 = 0.;
		
		hwTs = 2*PI*50*Ts;
		intx = 0.;
		
		bandera_desb = 0.;
		
		for (i=0; i<3;i++)
		{
			out1_2[i] = 0.;
			out1_1[i] = 0.;
			out2_2[i] = 0.;
			out2_1[i] = 0.;		
		}
		for (i=0; i<4;i++)
		{
			out1_LP_2[i] = 0.;
			out1_LP_1[i] = 0.;
			out2_LP_2[i] = 0.;
			out2_LP_1[i] = 0.;
			out3_LP_2[i] = 0.;
			out3_LP_1[i] = 0.;
			out4_LP_2[i] = 0.;
			out4_LP_1[i] = 0.;
			in1_LP_2[i] = 0.;
			in1_LP_1[i] = 0.;
			in2_LP_2[i] = 0.;
			in2_LP_1[i] = 0.;
			in3_LP_2[i] = 0.;
			in3_LP_1[i] = 0.;
			in4_LP_2[i] = 0.;
			in4_LP_1[i] = 0.;	
			
			out1_LP2_2[i] = 0.;
			out1_LP2_1[i] = 0.;
			out2_LP2_2[i] = 0.;
			out2_LP2_1[i] = 0.;
			out3_LP2_2[i] = 0.;
			out3_LP2_1[i] = 0.;
			out4_LP2_2[i] = 0.;
			out4_LP2_1[i] = 0.;			
			in1_LP2_2[i] = 0.;
			in1_LP2_1[i] = 0.;
			in2_LP2_2[i] = 0.;
			in2_LP2_1[i] = 0.;
			in3_LP2_2[i] = 0.;
			in3_LP2_1[i] = 0.;
			in4_LP2_2[i] = 0.;
			in4_LP2_1[i] = 0.;			
			
		}
		for (i=0;i<4;i++)
		{
			Intdh_p[i] = 0.;
			Intqh_p[i] = 0.;
			Intdh_n[i] = 0.;
			Intqh_n[i] = 0.;		
		}
		
		// Mantenemos inicializados los valores de los filtros
		lvc1   = *vc1r;//*1000.0;
		lvc1_1 = *vc1r;//*1000.0;
		lvc2   = *vc2r;//*1000.0;
		lvc2_1 = *vc2r;//*1000.0;
		lvc1_f1 = *vc1r;//*1000.0;
		lvc2_f1 = *vc2r;//*1000.0;
		lvc3   = *vc3r;//*1000.0;
		lvc3_1 = *vc3r;//*1000.0;
		lvc4   = *vc4r;//*1000.0;
		lvc4_1 = *vc4r;//*1000.0;
		lvc3_f1 = *vc3r;//*1000.0;
		lvc4_f1 = *vc4r;//*1000.0;
		
		/******	PONER A 0 TODOS LOS DISPAROS *******/
		for (i=0; i<5; i++)
		{
			*(tiempos_a+i) = 0.0;
			*(tiempos_b+i) = 0.0;
			*(tiempos_c+i) = 0.0;
		}
		if (*Enable == 1)
		{
			if (*inicio == 1)	estado = 2; //Si se está en STANDBY y se acciona el botón de inicio se puede pasar a modo Precarga
			else	estado = 1;
		}
	}
	/**********************/
	/* ESTADO DE PRECARGA */
	/**********************/
	
	else if (estado == 2) //ESTADO DE PRECARGA
	{
		if (*modo == RECTIFICADOR) //////////////////////////////////////////////////////////////////////////////////////////////////
		{
			/****** ACTIVAR CONTACTORES ****/
			Unset_contactor(Perifericos,C_PRECARGA);
			Unset_contactor(Perifericos,C_RED);
			Set_contactor(Perifericos,C_PRE_TRAFO);//Set
			Unset_contactor(Perifericos,C_TRAFO_RED);
			Unset_contactor(Perifericos,C_DC_SOURCE);     ////////////////////////////////////////////////////////////////////////////////////////////
			Unset_contactor(Perifericos,C_IN_EL_15);      ////////////////////////////////////////////////////////////////////////////////////////////
			Unset_contactor(Perifericos,C_OUT_EL_15);     ////////////////////////////////////////////////////////////////////////////////////////////
			Unset_contactor(Perifericos,C_BOBINA_CARGA);  ////////////////////////////////////////////////////////////////////////////////////////////
			Set_contactor(Perifericos,C_CARGA_ESTRELL);   //Se cierra el contactor para dejar las resistencias en paralelo ///////////////////////////
			Set_contactor(Perifericos,C_CHOPPER_CARGA);   //Se cierra el contactor para conectar el chopper a la carga////////////////////////////////
			/******************************/
			
			enable_disp = 0;
			
			if (temp1 >= 40000) //Contamos 4 segundos para conectar el contactor del trafo a red
			{
				Set_contactor(Perifericos, C_TRAFO_RED); // Conectado precarga del trafo y conexión trafo-red
				Unset_contactor(Perifericos,C_PRECARGA);
				Unset_contactor(Perifericos,C_RED);
				Set_contactor(Perifericos,C_PRE_TRAFO);
				Unset_contactor(Perifericos,C_DC_SOURCE);     ////////////////////////////////////////////////////////////////////////////////////////////
				Unset_contactor(Perifericos,C_IN_EL_15);      ////////////////////////////////////////////////////////////////////////////////////////////
				Unset_contactor(Perifericos,C_OUT_EL_15);     ////////////////////////////////////////////////////////////////////////////////////////////
				Unset_contactor(Perifericos,C_BOBINA_CARGA);  ////////////////////////////////////////////////////////////////////////////////////////////
				Set_contactor(Perifericos,C_CARGA_ESTRELL);   ////////////////////////////////////////////////////////////////////////////////////////////
				Set_contactor(Perifericos,C_CHOPPER_CARGA);   ////////////////////////////////////////////////////////////////////////////////////////////
				
				if (temp2 >= 20000) // Tras 2 segundos se conecta el contactor de precarga del sistema y se desconecta el de precarga del trafo
				{
					Set_contactor(Perifericos,C_PRECARGA); // Desconectar contactor precarga trafo-red y conectar Precarga del sistema
					Unset_contactor(Perifericos,C_RED);
					Unset_contactor(Perifericos,C_PRE_TRAFO);
					Set_contactor(Perifericos,C_TRAFO_RED);
					Unset_contactor(Perifericos,C_DC_SOURCE);     ////////////////////////////////////////////////////////////////////////////////////////////
					Unset_contactor(Perifericos,C_IN_EL_15);      ////////////////////////////////////////////////////////////////////////////////////////////
					Unset_contactor(Perifericos,C_OUT_EL_15);     ////////////////////////////////////////////////////////////////////////////////////////////
					Unset_contactor(Perifericos,C_BOBINA_CARGA);  ////////////////////////////////////////////////////////////////////////////////////////////
					Set_contactor(Perifericos,C_CARGA_ESTRELL);   ////////////////////////////////////////////////////////////////////////////////////////////
					Set_contactor(Perifericos,C_CHOPPER_CARGA);   ////////////////////////////////////////////////////////////////////////////////////////////
				}
				else
					temp2++;
			}
			else
				temp1++;
			
			
			Pcond = 0;
			Ierror_Vdc = 0.;
			
			if (*Enable == 0)
				estado = 1;
			
			/****** COMPROBAR NIVEL DE TENSIÓN *****/
			if ((*vc1r)+(*vc2r)+(*vc3r)+(*vc4r)>470) //Si la tensión de los condensadores es superior a 470 
			{
				estado = 3;
							
				temp1 = 0;
				temp2 = 0;
				
				Set_contactor(Perifericos,C_PRECARGA);// Se conecta el contactor de conexión directa del sistema a red
				Set_contactor(Perifericos,C_RED);
				Unset_contactor(Perifericos,C_PRE_TRAFO);
				Set_contactor(Perifericos,C_TRAFO_RED);
				Unset_contactor(Perifericos,C_DC_SOURCE);     ////////////////////////////////////////////////////////////////////////////////////////////
				Unset_contactor(Perifericos,C_IN_EL_15);      ////////////////////////////////////////////////////////////////////////////////////////////
				Unset_contactor(Perifericos,C_OUT_EL_15);     ////////////////////////////////////////////////////////////////////////////////////////////
				Unset_contactor(Perifericos,C_BOBINA_CARGA);  ////////////////////////////////////////////////////////////////////////////////////////////
				Set_contactor(Perifericos,C_CARGA_ESTRELL);   ////////////////////////////////////////////////////////////////////////////////////////////
				Set_contactor(Perifericos,C_CHOPPER_CARGA);   ////////////////////////////////////////////////////////////////////////////////////////////
			
			}
		}

		else if (*modo == INVERSOR)///////////////////////////////////////////////////////////////////////////////////////////////////
		{
			/****** ACTIVAR CONTACTORES ****/
			Unset_contactor(Perifericos,C_PRECARGA);
			Unset_contactor(Perifericos,C_RED);
			Unset_contactor(Perifericos,C_PRE_TRAFO);
			Unset_contactor(Perifericos,C_TRAFO_RED);
			Set_contactor(Perifericos,C_DC_SOURCE);       //Se conecta la fuente DC al dc-link del convertidor ///////////////////////////////////////
			Unset_contactor(Perifericos,C_IN_EL_15);      ////////////////////////////////////////////////////////////////////////////////////////////
			Unset_contactor(Perifericos,C_OUT_EL_15);     ////////////////////////////////////////////////////////////////////////////////////////////
			Unset_contactor(Perifericos,C_BOBINA_CARGA);  ////////////////////////////////////////////////////////////////////////////////////////////
			Unset_contactor(Perifericos,C_CARGA_ESTRELL); ////////////////////////////////////////////////////////////////////////////////////////////
			Unset_contactor(Perifericos,C_CHOPPER_CARGA); ////////////////////////////////////////////////////////////////////////////////////////////
			/******************************/
			
			enable_disp = 0;
			
			if (temp1 >= 30000) //Contamos 3 segundos para conectar el contactor del convertidor a la carga
			{
				Unset_contactor(Perifericos, C_TRAFO_RED); 
				Unset_contactor(Perifericos,C_PRECARGA);
				Unset_contactor(Perifericos,C_RED);
				Unset_contactor(Perifericos,C_PRE_TRAFO);
				Set_contactor(Perifericos,C_DC_SOURCE);     ////////////////////////////////////////////////////////////////////////////////////////////
				Unset_contactor(Perifericos,C_IN_EL_15);      ////////////////////////////////////////////////////////////////////////////////////////////
				Unset_contactor(Perifericos,C_OUT_EL_15);     ////////////////////////////////////////////////////////////////////////////////////////////
				Set_contactor(Perifericos,C_BOBINA_CARGA);  ////////////////////////////////////////////////////////////////////////////////////////////
				Unset_contactor(Perifericos,C_CARGA_ESTRELL); ////////////////////////////////////////////////////////////////////////////////////////////
				Unset_contactor(Perifericos,C_CHOPPER_CARGA); ////////////////////////////////////////////////////////////////////////////////////////////
			}
			else
				temp1++;
			
			
			Pcond = 0;
			Ierror_Vdc = 0.;
			
			if (*Enable == 0)
				estado = 1;
		}
			
	}/************ FIN ESTADO DE PRECARGA ************/
	
	/************ ESTADO DE MARCHA ************/

	else if (estado == 3) //ESTADO DE MARCHA
	{
		/**************************/
		/* CONSULTAR BAJA TENSIÓN */
		/**************************/
		
		if ((*vc1r)+(*vc2r)+(*vc3r)+(*vc4r)<450) // Si la tensión del dc-link es demasiada baja
		{
			estado = 0;
			error = 20; //error de tensión de condensador baja
		}
		//AÑADIR PARÁMETROS DE CONTROL A LA INTERFAZ
		
		if(temp1 > 20000)// Tras dos segundos se desconecta el contactor de precarga 
		{	Unset_contactor(Perifericos,C_PRECARGA); // Se desconecta el trafo de precarga
			Set_contactor(Perifericos,C_RED);
			Unset_contactor(Perifericos,C_PRE_TRAFO);
			Set_contactor(Perifericos,C_TRAFO_RED);
			Unset_contactor(Perifericos,C_DC_SOURCE);     ////////////////////////////////////////////////////////////////////////////////////////////
			Unset_contactor(Perifericos,C_IN_EL_15);      ////////////////////////////////////////////////////////////////////////////////////////////
			Unset_contactor(Perifericos,C_OUT_EL_15);     ////////////////////////////////////////////////////////////////////////////////////////////
			Unset_contactor(Perifericos,C_BOBINA_CARGA);  ////////////////////////////////////////////////////////////////////////////////////////////
			Set_contactor(Perifericos,C_CARGA_ESTRELL); ////////////////////////////////////////////////////////////////////////////////////////////
			Set_contactor(Perifericos,C_CHOPPER_CARGA); ////////////////////////////////////////////////////////////////////////////////////////////
		}
		else
			temp1++;
		
		/********************************************/
		/* VOLVER A ESTADO DE STANDBY SI ENABLE == 0*/
		/********************************************/
		
		if (*Enable == 0)
			estado = 1;
			
		/****	CICLO DE MUESTRO Y CONTROL	****/
					
			/****************************************************************/
			/* Filtrado de las tensiones medidas de los condensadores       */
			/****************************************************************/
			
			lvc1=*(vc1r);
			lvc2=*(vc2r);
			lvc1_f = (c11*(lvc1+lvc1_1) + c12*lvc1_f1)*c13; //Filtro paso bajo sobre vc1 y vc2
			lvc2_f = (c11*(lvc2+lvc2_1) + c12*lvc2_f1)*c13;
			lvc1_1 = lvc1;
			lvc2_1 = lvc2;
			lvc1_f1 = lvc1_f;
			lvc2_f1 = lvc2_f;

			lvc3=*(vc3r);
			lvc4=*(vc4r);
			lvc3_f = (c11*(lvc3+lvc3_1) + c12*lvc3_f1)*c13; //Filtro paso bajo sobre vc3 y vc4
			lvc4_f = (c11*(lvc4+lvc4_1) + c12*lvc4_f1)*c13;
			lvc3_1 = lvc3;
			lvc4_1 = lvc4;
			lvc3_f1 = lvc3_f;
			lvc4_f1 = lvc4_f;

			Vcond = lvc1_f + lvc2_f + lvc3_f + lvc4_f; //Como vemos Vcond es la tensión total filtrada como vc1, vc2, vc3 y vc4

			Inv_Vcond = 1./Vcond;
			Vcond_f = (c31*(Vcond+Vcond1) + c32*Vcond_f1)*c33;
			Vcond_f1 = Vcond_f;
			Vcond1 = Vcond;
			Vcond = Vcond_f;
			
			Vcond = (*vc1r)+(*vc2r)+(*vc3r)+(*vc4r); //Ignore low-pass filter

			/*********************/
			/* COMPROBAR ERRORES */
			/*********************/
			
			//Check error devuelve cero en caso de no error
			if (error == 0)
			{
				error = Check_errors(*A_max, *V_max, Iu, Iv, Iw, Vr, Vs, Vt); //La variable error devuelve el tipo de error
				if( error != 0)//¿Algún error?									1,2,3 SOBRECORRIENTE instantanea FASE A,B,C respectivamente
				//																4,5,6 SOBRETENSION instantanea FASE A,B,C respectivamente
					estado = 0;
			}
			if (error == 0)// Para no pisar el error anterior
			{
				error = Check_errors(*A_RMS_max, *V_RMS_max,RMS_Ia, RMS_Ib,RMS_Ic,RMS_Va,RMS_Vb,RMS_Vc)+7.;
				if( error != 7)//¿Algún error?									8,9,10 		SOBRECORRIENTE en RMS FASE A,B,C respectivamente
				//																11,12,13 	SOBRETENSIÓN en RMS FASE A,B,C respectivamente
					estado = 0;
			}
			
			if (error == 7)
			{
				error = Check_error_vdc(*Vcr_max, *vc1r, *vc2r, *vc3r, *vc4r)+14.;
				if (error != 14)//												15,16,17,18		SOBRETENSIÓN CONDENSADOR VC1,VC2,VC3,VC4 respectivamente
					estado = 0;
			}
			if (*e_drv != 0)// Si hay error de driver
			{	
				estado = 0;
				error = 19;
			}
				
			if (estado == 0)
			{
						for (i=0; i<5; i++) //Poner tiempos a cero
						{
							*(tiempos_a+i) = 0.0;
							*(tiempos_b+i) = 0.0;
							*(tiempos_c+i) = 0.0;
						}
			}
			
			/***************************************************/
			/* Comprobar cercanías de las magnitudes a errores */
			/***************************************************/
			
			cercania = Check_cercania(*Vcr_max, *A_RMS_max, *V_RMS_max, RMS_Ia, RMS_Ib, RMS_Ic, *vc1r, *vc2r, *vc3r, *vc4r);
			
		/************************************/
		/* CONTROL Y GENERACIÓN DE DISPAROS */
		/************************************/
			
		if (*go == 1) // Si está pulsada la tecla de GO
		{		

			/*******************************************/
			/* Rampa de la tension de referencia de CC */
			/*******************************************/
				
			lvdcref = rampa_V_ref(lvdcref,Ts,(float)((*vc1r)+(*vc2r)+(*vc3r)+(*vc4r)),(*Vdcref));		
		
			enable_disp = 1; // Se habilitan los disparos
			
			/***********************************/
 			/* PI de tension de condensadores  */
			/***********************************/
			
			error_Vdc = lvdcref*lvdcref - Vcond*Vcond; // Con este criterio es potencia a inyectar en condensadores
			Ierror_Vdc += (error_Vdc + error_Vdc1)*0.5;
			
			Perror_Vdc = kps_Pref*error_Vdc;
			//Pcond = Perror_Vdc + kis*Ierror_Vdc;// Proporcional más integral
			Pcond = kps_Pref*(Pcond1*(2 - 2*pi*5000*Ts) +  error_Vdc1*2*pi*5000*Ts + error_Vdc*2*pi*5000*Ts)/(2*pi*5000*Ts + 2) + kis*Ierror_Vdc;// -> Integral más paso bajo
//			salida = (salida_ant*(2-Wc*Ts) + entrada_ant*K*Wc*Ts + entrada*K*Wc*Ts)/(Wc*Ts+2);// -> Filtro paso bajo			
			
			error_Vdc1 = error_Vdc;
			Pcond1 = Pcond;
			Perror_Vdc1 = Perror_Vdc;
			
			/********************************************************/
 			/* Cálculo de la potencia activa y reactiva instantánea */
			/********************************************************/
			
			p = (valfa*ialfa_c + vbeta*ibeta_c);
			q = (valfa*ibeta_c - vbeta*ialfa_c);
			
			qref = *Qref;
						
            ModVred = (valfa*valfa + vbeta*vbeta);
			
			/***************************/
			/* Estimación de la bobina */
			/***************************/
			
			intx = intx + Pcond*q*0.0001;
			X_estimado = -0.00000001*intx; // <-- No funciona
			
			/********************************/
			/* Obtención corrientes Id y Iq */
			/********************************/
			
            ia_ref = 1./ModVred*(valfa*Pcond - vbeta*qref);
            ib_ref = 1./ModVred*(vbeta*Pcond + valfa*qref);
						
			id_ref = ia_ref * cos(theta) + ib_ref * sin(theta);
			iq_ref = -ia_ref * sin(theta) + ib_ref * cos(theta);
			
			iar_ref = raiz2_3*ia_ref;
			ibr_ref = -raiz1_6*ia_ref + raiz1_2*ib_ref;
			icr_ref = -raiz1_6*ia_ref - raiz1_2*ib_ref;
			
			/**************************/
			/* Detección de armónicos */
			/**************************/
			kt = 2./Ts;	
			KSOGI = raiz2;
			for (i=0;i<4;i++)
			{
				switch(i){
				case 0:
					wcut = 2.*pi*70;			
					nume1 = kt*kt + KSOGI*wcut*kt + wcut*wcut;
					nume2 = -2.*kt*kt + 2.*wcut*wcut;
					nume3 = kt*kt - KSOGI*wcut*kt + wcut*wcut;
					deno1 = KSOGI*wcut*wcut;					
					idh_p[i] = ialfa_c * cos(5.*theta) + ibeta_c * sin(5.*theta);	// Synchronous frame for 5th harmonic
					iqh_p[i] = -ialfa_c * sin(5.*theta) + ibeta_c * cos(5.*theta);
					idh_n[i] = ialfa_c * cos(-5.*theta) + ibeta_c * sin(-5.*theta);
					iqh_n[i] = -ialfa_c * sin(-5.*theta) + ibeta_c * cos(-5.*theta);
					break;
				case 1:
					wcut = 2.*pi*100;			
					nume1 = kt*kt + KSOGI*wcut*kt + wcut*wcut;
					nume2 = -2.*kt*kt + 2.*wcut*wcut;
					nume3 = kt*kt - KSOGI*wcut*kt + wcut*wcut;
					deno1 = KSOGI*wcut*wcut;					
					idh_p[i] = ialfa_c * cos(7.*theta) + ibeta_c * sin(7.*theta);	// Synchronous frame for 5th harmonic
					iqh_p[i] = -ialfa_c * sin(7.*theta) + ibeta_c * cos(7.*theta);
					idh_n[i] = ialfa_c * cos(-7.*theta) + ibeta_c * sin(-7.*theta);
					iqh_n[i] = -ialfa_c * sin(-7.*theta) + ibeta_c * cos(-7.*theta);
					break;
				case 2:
					idh_p[i] = ialfa_c * cos(11.*theta) + ibeta_c * sin(11.*theta);	// Synchronous frame for 5th harmonic
					iqh_p[i] = -ialfa_c * sin(11.*theta) + ibeta_c * cos(11.*theta);
					idh_n[i] = ialfa_c * cos(-11.*theta) + ibeta_c * sin(-11.*theta);
					iqh_n[i] = -ialfa_c * sin(-11.*theta) + ibeta_c * cos(-11.*theta);
					break;
				case 3:
					idh_p[i] = ialfa_c * cos(13.*theta) + ibeta_c * sin(13.*theta);	// Synchronous frame for 5th harmonic
					iqh_p[i] = -ialfa_c * sin(13.*theta) + ibeta_c * cos(13.*theta);
					idh_n[i] = ialfa_c * cos(-13.*theta) + ibeta_c * sin(-13.*theta);
					iqh_n[i] = -ialfa_c * sin(-13.*theta) + ibeta_c * cos(-13.*theta);
					break;
				}
				// Low-pass filter to retrieve dc components (two cascaded 2nd order LP filter)
				// 1st filter
				in1_LP[i] = idh_p[i];
				in2_LP[i] = iqh_p[i];
				in3_LP[i] = idh_n[i];
				in4_LP[i] = iqh_n[i];
				
				out1_LP[i] = 1./nume1*(deno1*(in1_LP[i] + 2*in1_LP_1[i] + in1_LP_2[i]) - out1_LP_1[i]*nume2 - out1_LP_2[i]*nume3);
				out2_LP[i] = 1./nume1*(deno1*(in2_LP[i] + 2*in2_LP_1[i] + in2_LP_2[i]) - out2_LP_1[i]*nume2 - out2_LP_2[i]*nume3);				
				out3_LP[i] = 1./nume1*(deno1*(in3_LP[i] + 2*in3_LP_1[i] + in3_LP_2[i]) - out3_LP_1[i]*nume2 - out3_LP_2[i]*nume3);
				out4_LP[i] = 1./nume1*(deno1*(in4_LP[i] + 2*in4_LP_1[i] + in4_LP_2[i]) - out4_LP_1[i]*nume2 - out4_LP_2[i]*nume3);	
				
				out1_LP_2[i] = out1_LP_1[i];
				out1_LP_1[i] = out1_LP[i];
				out2_LP_2[i] = out2_LP_1[i];
				out2_LP_1[i] = out2_LP[i];
				out3_LP_2[i] = out3_LP_1[i];
				out3_LP_1[i] = out3_LP[i];
				out4_LP_2[i] = out4_LP_1[i];
				out4_LP_1[i] = out4_LP[i];
				in1_LP_2[i] = in1_LP_1[i];
				in1_LP_1[i] = in1_LP[i];
				in2_LP_2[i] = in2_LP_1[i];
				in2_LP_1[i] = in2_LP[i];
				in3_LP_2[i] = in3_LP_1[i];
				in3_LP_1[i] = in3_LP[i];
				in4_LP_2[i] = in4_LP_1[i];
				in4_LP_1[i] = in4_LP[i];
				
				// 2nd filter				
				in1_LP2[i] = out1_LP[i];
				in2_LP2[i] = out2_LP[i];
				in3_LP2[i] = out3_LP[i];
				in4_LP2[i] = out4_LP[i];
				
				out1_LP2[i] = 1./nume1*(deno1*(in1_LP2[i] + 2*in1_LP2_1[i] + in1_LP2_2[i]) - out1_LP2_1[i]*nume2 - out1_LP2_2[i]*nume3);
				out2_LP2[i] = 1./nume1*(deno1*(in2_LP2[i] + 2*in2_LP2_1[i] + in2_LP2_2[i]) - out2_LP2_1[i]*nume2 - out2_LP2_2[i]*nume3);				
				out3_LP2[i] = 1./nume1*(deno1*(in3_LP2[i] + 2*in3_LP2_1[i] + in3_LP2_2[i]) - out3_LP2_1[i]*nume2 - out3_LP2_2[i]*nume3);
				out4_LP2[i] = 1./nume1*(deno1*(in4_LP2[i] + 2*in4_LP2_1[i] + in4_LP2_2[i]) - out4_LP2_1[i]*nume2 - out4_LP2_2[i]*nume3);	
				
				out1_LP2_2[i] = out1_LP2_1[i];
				out1_LP2_1[i] = out1_LP2[i];
				out2_LP2_2[i] = out2_LP2_1[i];
				out2_LP2_1[i] = out2_LP2[i];
				out3_LP2_2[i] = out3_LP2_1[i];
				out3_LP2_1[i] = out3_LP2[i];
				out4_LP2_2[i] = out4_LP2_1[i];
				out4_LP2_1[i] = out4_LP2[i];
				in1_LP2_2[i] = in1_LP2_1[i];
				in1_LP2_1[i] = in1_LP2[i];
				in2_LP2_2[i] = in2_LP2_1[i];
				in2_LP2_1[i] = in2_LP2[i];
				in3_LP2_2[i] = in3_LP2_1[i];
				in3_LP2_1[i] = in3_LP2[i];
				in4_LP2_2[i] = in4_LP2_1[i];
				in4_LP2_1[i] = in4_LP2[i];				
				
				// PI Controller to generate the references
				Intdh_p[i] += out1_LP[i]*Ts;//out1_LP2[i]*Ts;
				Intqh_p[i] += out2_LP[i]*Ts;//out2_LP2[i]*Ts;
				Intdh_n[i] += out3_LP[i]*Ts;//out3_LP2[i]*Ts;
				Intqh_n[i] += out4_LP[i]*Ts;//out4_LP2[i]*Ts;
				
				PIdh_p[i] = 0.1*out1_LP[i] + 0.5*Intdh_p[i];//0.1*out1_LP2[i] + 0.2*Intdh_p[i];	// 0.05 0.005	// 0.05 0.005
				PIqh_p[i] = 0.1*out2_LP[i] + 0.5*Intqh_p[i];// 0.1*out2_LP2[i] + 0.2*Intqh_p[i];
				PIdh_n[i] = 0.1*out3_LP[i] + 0.5*Intdh_n[i];//0.1*out3_LP2[i] + 0.2*Intdh_n[i];
				PIqh_n[i] = 0.1*out4_LP[i] + 0.5*Intqh_n[i];//0.1*out4_LP2[i] + 0.2*Intqh_n[i];
				
				// Back to stationary frame
				switch(i){
					case 0:
						ialh_p[i] = PIdh_p[i] * cos(-5.*theta) + PIqh_p[i] * sin(-5.*theta);	// Stationary frame for 5th harmonic
						ibeh_p[i] = -PIdh_p[i] * sin(-5.*theta) + PIqh_p[i] * cos(-5.*theta);
						ialh_n[i] = PIdh_n[i] * cos(5.*theta) + PIqh_n[i] * sin(5.*theta);
						ibeh_n[i] = -PIdh_n[i] * sin(5.*theta) + PIqh_n[i] * cos(5.*theta);
						break;
					case 1:
						ialh_p[i] = PIdh_p[i] * cos(-7.*theta) + PIqh_p[i] * sin(-7.*theta);	// Stationary frame for 7th harmonic
						ibeh_p[i] = -PIdh_p[i] * sin(-7.*theta) + PIqh_p[i] * cos(-7.*theta);
						ialh_n[i] = PIdh_n[i] * cos(7.*theta) + PIqh_n[i] * sin(7.*theta);
						ibeh_n[i] = -PIdh_n[i] * sin(7.*theta) + PIqh_n[i] * cos(7.*theta);
						break;
					case 2:
						ialh_p[i] = PIdh_p[i] * cos(-11.*theta) + PIqh_p[i] * sin(-11.*theta);	// Stationary frame for 11th harmonic
						ibeh_p[i] = -PIdh_p[i] * sin(-11.*theta) + PIqh_p[i] * cos(-11.*theta);
						ialh_n[i] = PIdh_n[i] * cos(11.*theta) + PIqh_n[i] * sin(11.*theta);
						ibeh_n[i] = -PIdh_n[i] * sin(11.*theta) + PIqh_n[i] * cos(11.*theta);
						break;
					case 3:
						ialh_p[i] = PIdh_p[i] * cos(-13.*theta) + PIqh_p[i] * sin(-13.*theta);	// Stationary frame for 13th harmonic
						ibeh_p[i] = -PIdh_p[i] * sin(-13.*theta) + PIqh_p[i] * cos(-13.*theta);
						ialh_n[i] = PIdh_n[i] * cos(13.*theta) + PIqh_n[i] * sin(13.*theta);
						ibeh_n[i] = -PIdh_n[i] * sin(13.*theta) + PIqh_n[i] * cos(13.*theta);
						break;
				}
				
				//if (*enable_desbalance == 0){
					ialh[i] = ialh_p[i] + ialh_n[i];
					ibeh[i] = ibeh_p[i] + ibeh_n[i];
				/*}
				else
				{
					ialh[i] = 0.;
					ibeh[i] = 0.;
				}*/
			}
			ia_ref -= ialh[0] + ialh[1];	// 5th and 7th harmonic sustracted from the reference
			ib_ref -= ibeh[0] + ibeh[1];	// 5th and 7th harmonic sustracted from the reference
			
			/****************/
			/* Alfa-beta PR */
			/****************/
			
			//Resonant controller KP + Ki*s/(s^2 + w^2) -> output = Tz/2*1/(1 + Cw^2)[In - InZ-2] - 2*OuZ-1*(-1 + Cw^2)/(1 + Cw^2) - OuZ-2;	// Cw = Tz/2*Wres
			
			kt = 2./Ts;		
			kp_res = Kp;
			in1 = ia_ref - ialfa_c;			
			in2 = ib_ref - ibeta_c;		
			
			//wres[0] = 2.*pi*50*1.;
			//wres[1] = 2.*pi*50*5.;
			//wres[2] = 2.*pi*50*7.;
			wres[0] = wg;
			wres[1] = wg*5.;
			wres[2] = wg*7.;
			
			for (i=0; i<3; i++){
				
				switch(i){
					case 0:
						wcut = 2.*pi*0.1;
						ki_res = Ki;	// Read from parameter input
						break;
					case 1:
						wcut = 2.*pi*0.2;
						ki_res = 50.;	// Fixed for 5th harmonic
						break;
					case 2:
						wcut = 2.*pi*0.2;
						ki_res = 70.;	// Fixed for 7th harmonic
						break;
				}		

				den1 = 2.*ki_res*kt*wcut;
				den2 = den1;

				n0[i] = kt*kt + 2.*kt*wcut + wres[i]*wres[i];
				n1[i] = 2.*kt*kt - 2.*wres[i]*wres[i];
				n2[i] = kt*kt - 2.*kt*wcut + wres[i]*wres[i];//kt*kt + 2.*kt*wcut + wres*wres ;
				
				out1[i] = 1./n0[i]*(den1*(in1 - in1_2[i]) + n1[i]*out1_1[i] - n2[i]*out1_2[i]);//1./n0*(den1*(in1_1 - in1_2) + n1*out1_1 - n2*out1_2);
				out2[i] = 1./n0[i]*(den1*(in2 - in2_2[i]) + n1[i]*out2_1[i] - n2[i]*out2_2[i]);//1./n0*(den1*(in2_1 - in2_2) + n1*out2_1 - n2*out2_2);

				out1_2[i] = out1_1[i];
				out1_1[i] = out1[i];
				out2_2[i] = out2_1[i];
				out2_1[i] = out2[i];
				
				in1_2[i] = in1_1[i];
				in1_1[i] = in1;			
				in2_2[i] = in2_1[i];
				in2_1[i] = in2;
				
				//	Refalfa = Refalfa + 0.01*out1; -> Actualizar control de la siguiente forma
				//	Refbeta = Refbeta + 0.01*out2;
			}
			
			RefAlpha = 2./Vcond*(-kp_res*in1 - out1[0] - out1[1] - out1[2] + valfa);
			RefBeta = 2./Vcond*(-kp_res*in2 - out2[0] - out2[1] - out2[2] + vbeta);
			
			/*******************************/
			/* Control de corrientes en dq */
			/*******************************/
			/*
			IntPk0 = IntPk1 + Ki*(id_ref - id);
            IntQk0 = IntQk1 + Ki*(iq_ref - iq);
			
			IntPk1=IntPk0;
			IntQk1=IntQk0;

			Refd    	= 2./Vcond*(-dc_d*0 -Kp*(id_ref - id) + iq*2*PI*50*0.002 - IntPk0 + Vd); //originalmente -iq	//Inductancias de 2 mH
            Refq     	= 2./Vcond*(-dc_q*0 -Kp*(iq_ref - iq) - id*2*PI*50*0.002 - IntQk0 + Vq);//originalmente +id 	//Inductancias de 2 mH	

			RefAlpha 	= 	Refd * cos(theta) - Refq * sin(theta);
			RefBeta 	=	Refd * sin(theta) + Refq * cos(theta);
			*/
            /****************************/
            /* Control DPC mediante ORS */
			/****************************/
			
			// Sustracción componente fundamental de secuencia positiva a las tensiones de red
			/*
			x1alfa = valfa - valfa_res_1;
			x1beta = vbeta - vbeta_res_1;
			
			x2alfa_p = x1alfa*cos(theta) + x1beta*sin(theta);
			x2beta_p = -x1alfa*sin(theta) + x1beta*cos(theta);
			x2alfa_n = x1alfa*cos(theta) - x1beta*sin(theta);
			x2beta_n = x1alfa*sin(theta) + x1beta*cos(theta);
			
			gamma = 1.;
			
			x3alfa_p += Ts/(2*gamma)*(x2alfa_p + x2alfa_p_1);
			x3beta_p += Ts/(2*gamma)*(x2beta_p + x2beta_p_1);
			x3alfa_n += Ts/(2*gamma)*(x2alfa_n + x2alfa_n_1);
			x3beta_n += Ts/(2*gamma)*(x2beta_n + x2beta_n_1);
			
			x4alfa_p = x3alfa_p*cos(theta) - x3beta_p*sin(theta);//valfa_p
			x4beta_p = x3alfa_p*sin(theta) + x3beta_p*cos(theta);//vbeta_p
			x4alfa_n = x3alfa_n*cos(theta) + x3beta_n*sin(theta);
			x4beta_n = -x3alfa_n*sin(theta) + x3beta_n*cos(theta);
			
			valfa_res = x4alfa_p + x4alfa_n;
			vbeta_res = x4beta_p + x4beta_n;
		
			valfa_res_1 = valfa_res;
			vbeta_res_1 = vbeta_res;
			
			x2alfa_p_1 = x2alfa_p;
			x2beta_p_1 = x2beta_p;
			x2alfa_n_1 = x2alfa_n;
			x2beta_n_1 = x2beta_n;
			
			/////////////////
			
			if (x4alfa_p*x4alfa_p + x4beta_p*x4beta_p < 0.8*ModVred) //Aseguramos que está bien sintonizado la secuencia positiva del 1º armónico
			{	x4alfa_p = valfa;
				x4beta_p = vbeta;}
				
				// Use direct measure 
				//x4alfa_p = valfa;
				//x4beta_p = vbeta;
			
			Xest = 0.002*2*PI*50;//L = 0.002

			p = (valfa*ialfa_c + vbeta*ibeta_c);
			q = (valfa*ibeta_c - vbeta*ialfa_c);
			
			perror = p - Pcond;
			qerror = q - qref;
			
			//  Addition of harmonic references to supress them 
			//perror = p - (Pcond - valfa*(ialh[0] + ialh[1]) - vbeta*(ibeh[0] + ibeh[1]));	
			//qerror = q - (qref - valfa*(ibeh[0] + ibeh[1]) + vbeta*(ialh[0] + ialh[1]));	
			//
			
			ModAB1 = 2./Vcond*(1 + Xest*q/ModVred);

			ORS1AB.Calpha 	= ModAB1*x4alfa_p;//valfa;// Esta valfa y vbeta se refieren a la magnitud del armónico fundamental
			ORS1AB.Cbeta 	= ModAB1*x4beta_p;//vbeta;

			ModAB2 = -2./Vcond*(Xest*p/ModVred);
			
			// J = [0 -1; 1 0]
			ORS2AB.Calpha 	= -ModAB2*x4beta_p;//vbeta;// Igual al anterior
			ORS2AB.Cbeta 	= ModAB2*x4alfa_p;//valfa;
			
			IntPk0 = IntPk1 + Ki*perror;
			IntQk0 = IntQk1 + Ki*qerror;
			
			IntPk1=IntPk0;
			IntQk1=IntQk0;	
			
			RefAlpha = ORS1AB.Calpha + ORS2AB.Calpha + Kp*(perror*valfa - qerror*vbeta) + IntPk0*valfa - IntQk0*vbeta;//valfa
			RefBeta = ORS1AB.Cbeta + ORS2AB.Cbeta + Kp*(perror*vbeta + qerror*valfa) + IntPk0*vbeta + IntQk0*valfa;//vbeta
			*/
			
			/*************************/
			
			u1 = RefAlpha;
			u2 = RefBeta;
			
			/*****************/
			/* Modo Inversor */
			/*****************/
			/*
			temp3++;
			if (temp3 > 1./(50*Ts))
				temp3 = 0;
			
			u1 = (*K6h)*cos(2*pi*50*Ts*temp3);
			u2 = (*K6h)*sin(2*pi*50*Ts*temp3);	
			*/
			/**************************************/
 			/* PI de diferencia de condensadores  */
			/**************************************/			
	
			//En el convertidor los condensadores C1,C2,C3,C4 van de arriba a abajo, siendo C1 el superior y C4 el inferior
			//xdifk0 = ((*vc2r)-(*vc1r))*1000.0;
			
			vd1k0 = (*vc3r + *vc4r - *vc1r - *vc2r);
			/*** EN modo inversor es al revés ***/
			//vd1k0 = (*vc2r + *vc1r - *vc3r - *vc4r);
			
			/** Código para provocar desbalance **/
			if (*enable_desbalance == 1)
			{
				vd1k0 = vd1k0 + 45;// 3 NIVELES
				k_balancing = 0.02;
				ki_balancing = 0.0000001;
			}
			/** Fin código desbalance **/			
			
			Int_err_xdif += (vd1k0 + vd1k0_1)*0.5;
			PI_err_xdif = k_balancing*vd1k0 + ki_balancing*Int_err_xdif;
			kbal1 = (valfa*p - vbeta*q)/(p*p + q*q);
			kbal2 = (vbeta*p + valfa*q)/(p*p + q*q);
			
			vd1k0_1 = vd1k0;
			
			u3 = PI_err_xdif*kbal1;
			u4 = PI_err_xdif*kbal2;		
						
			// Cambio de variables de las señales de control:
			dalfa1 = -0.5*u1 + 0.5*u3;
			dalfa2 = 0.5*u1 + 0.5*u3;

			dbeta1 = -0.5*u2 + 0.5*u4;
			dbeta2 = 0.5*u2 + 0.5*u4;

			// Cambio de coordenadas de las señales de control:
			
			/**************/
			/* MODULACIÓN */
			/**************/
			
			/***			***/
			/* SVM modified ***/
			/***			***/
			
			//SVM_mod(u1, u2, *vc1r, *vc2r, *vc3r, *vc4r, &Vec[0], &dv[0]);
			
			/***   ***/
			/* IECON */
			/***   ***/
			/*
			dgamma1 = 0.75;
			dgamma2 = 0.75;
			
			da1 = raiz2_3*dalfa1 + raiz1_3*dgamma1;
			da2 = 0.;
			da3 = 0.;
			da4 = raiz2_3*dalfa2 + raiz1_3*dgamma2;
			
			db1 = -raiz1_6*dalfa1 + raiz1_2*dbeta1 + raiz1_3*dgamma1;
			db2 = 0.;
			db3 = 0.;
			db4 = -raiz1_6*dalfa2 + raiz1_2*dbeta2 + raiz1_3*dgamma2;
			
			dc1 = -raiz1_6*dalfa1 - raiz1_2*dbeta1 + raiz1_3*dgamma1;
			dc2 = 0.;
			dc3 = 0.;
			dc4 = -raiz1_6*dalfa2 - raiz1_2*dbeta2 + raiz1_3*dgamma2;
			*/
			/***        ***/
			/* VARIANTE 1 */
			/***        ***/
			/*
			// Nivel positivo
			caso = 0;

			db = raiz2*dbeta2/2 - raiz6*dalfa2/2;
			dc = -raiz6*dalfa2/2 - raiz2*dbeta2/2;

			if ((db >= 0) && (db <= 1) && (dc >= 0) && (dc <= 1)) // Caso 1 -> da2 = 0
			{
				caso = 1;
				da4 = 0.;
				db4 = db;
				dc4 = dc;
			}

			da = raiz6*dalfa2/2 - raiz2*dbeta2/2;
			dc = -raiz2*dbeta2;

			if ((da >= 0) && (da <= 1) && (dc >= 0) && (dc <= 1)) // Caso 2 -> db2 = 0
			{
				caso = 2;
				da4 = da;
				db4 = 0.;
				dc4 = dc;
			}

			da = raiz6*dalfa2/2 + raiz2*dbeta2/2;
			db = raiz2*dbeta2;

			if ((da >= 0) && (da <= 1) && (db >= 0) && (db <= 1)) // Caso 3 -> dc2 = 0
			{
				caso = 3;
				da4 = da;
				db4 = db;
				dc4 = 0.;
			}

			if (caso == 0)
			{
				da4 = raiz2_3*dalfa2 + raiz1_3*0.75;
				db4 = -raiz1_6*dalfa2 + raiz1_2*dbeta2 + raiz1_3*0.75;
				dc4 = -raiz1_6*dalfa2 - raiz1_2*dbeta2 + raiz1_3*0.75;
			}

			//Nivel negativo
			caso = 0;

			db = raiz2*dbeta1/2 - raiz6*dalfa1/2;
			dc = -raiz6*dalfa1/2 - raiz2*dbeta1/2;

			if ((db >= 0 ) && (db <= 1) && (dc >= 0) && (dc <= 1))
			{
				caso=1;
				da1 = 0;
				db1 = db;
				dc1 = dc;
			}

			da = (raiz6*dalfa1)/2 - (raiz2*dbeta1)/2;
			dc = -raiz2*dbeta1;

			if ((da >=0 ) && (da <= 1) && (dc >= 0) && (dc <= 1))
			{
				caso=2;
				da1=da;
				db1=0;
				dc1=dc;
			}

			da = (raiz6*dalfa1)/2 + (raiz2*dbeta1)/2;
			db = raiz2*dbeta1;

			if ((da >=0 ) && (da <= 1) && (db >= 0) && (db <= 1))
			{	
				caso=3;
				da1=da;
				db1=db;
				dc1=0;
			}

			if (caso == 0)
			{
				da1 = raiz2_3*dalfa1 + raiz1_3*0.75;			
				db1 = -raiz1_6*dalfa1 + raiz1_2*dbeta1 + raiz1_3*0.75;
				dc1 = -raiz1_6*dalfa1 - raiz1_2*dbeta1 + raiz1_3*0.75;
			}
			
			da2 = 0.;
			da3 = 0.;
			db2 = 0.;
			db3 = 0.;
			dc2 = 0.;
			dc3 = 0.;
			*/
			/***        ***/
			/* VARIANTE 2 */
			/***        ***/
			
			alpha1 = raiz2_3*RefAlpha;
	     	alpha2 = -RefAlpha*raiz1_6 + RefBeta*raiz1_2;
	     	alpha3 = -RefAlpha*raiz1_6 - RefBeta*raiz1_2;

	     	if (alpha1 > alpha2)
	     	{
	     		alphamin = alpha2;
	     		alphamax = alpha1;
	     	}
	     	else
	     	{
	     		alphamax = alpha2;
	     		alphamin = alpha1;
	     	}
	     	if (alpha3 < alphamin)
	     	{
	     		alphamin = alpha3;
	     	}
	     	if (alpha3 > alphamax)
	     	{
	     		alphamax = alpha3;
	     	}

	     	xmax = (1 - alphamax);
	     	xmin = (-1 - alphamin);

	     	u_amax = alpha1 + xmax;
	     	u_amin = alpha1 + xmin;
	     	u_bmax = alpha2 + xmax;
	     	u_bmin = alpha2 + xmin;
	     	u_cmax = alpha3 + xmax;
	     	u_cmin = alpha3 + xmin;

	     	f1 = Iu*abso(u_amax) + Iv*abso(u_bmax) + Iw*abso(u_cmax);
	     	f2 = Iu*abso(u_amin) + Iv*abso(u_bmin) + Iw*abso(u_cmin);

	     	yx0 = alpha2-alpha1;
	     	zx0 = alpha3-alpha1;
	     	f3 = Iv*abso(yx0) + Iw*abso(zx0);		// Si ua corta al cero

	     	xy0 = alpha1-alpha2;
	     	zy0 = alpha3-alpha2;
	     	f4 = Iu*abso(xy0) + Iw*abso(zy0);		// Si ub corta al cero

	     	xz0 = alpha1-alpha3;
	     	yz0 = alpha2-alpha3;
	     	f5 = Iu*abso(xz0) + Iv*abso(yz0);		// Si uc corta al cero

	     	vec_min[0] = -signo(vd1k0)*f1;
	     	vec_min[1] = -signo(vd1k0)*f2;
	    	vec_min[2] = -signo(vd1k0)*f3;
	     	vec_min[3] = -signo(vd1k0)*f4;
	     	vec_min[4] = -signo(vd1k0)*f5;
			
			// Low-pass filter with cost Functions vec_min[]
			aux = 2*pi*50*Ts;
			
				for (i=0;i<5;i++){
					vec_min_PB[i] = (vec_min_PB_ant[i]*(2 - aux) + aux*(vec_min[i] + vec_min_ant[i]))/(2 + aux);
					vec_min_PB_ant[i] = vec_min_PB[i];
					vec_min_ant[i] = vec_min[i];
					// Replace vec_min with low-passed one			
					vec_min[i] = vec_min_PB[i];
				}
			// End Low-pass filter
			
			imin = imin_ant;	// Select last criteria when no 3 lvl was considered
			
			if (imin > 0 && imin < 6) // in case there were no previous case selected
				Vmin = vec_min[imin-1];
		
			Banda = 0.;
			
			if (signo(vd1k0)*vd1k0 >= Banda || imin == 0)// If its outside the band or previous x is not feasible
			{
				if (vec_min[0] <= vec_min[1])
				{
					Vmin = vec_min[0];
					imin = 1;
				}
				else
				{
					Vmin = vec_min[1];
					imin = 2;
				}
				if (vec_min[2] < Vmin && (0>u_amin) && (0<u_amax))
				{
					Vmin = vec_min[2];
					imin = 3;
				}
				if (vec_min[3] < Vmin && (0>u_bmin) && (0<u_bmax))
				{
					Vmin = vec_min[3];
					imin = 4;
				}
				if (vec_min[4] < Vmin && (0>u_cmin) && (0<u_cmax))
				{
					Vmin = vec_min[4];
					imin = 5;
				}
				imin_ant = imin;
			}
	     	// Until here, the old V2 approach is adopted
	     	// This is the new approach for V2 when the unbalance cannot be solved by using the old one
	     	
	  	Vmin_ant = Vmin;
		epsilon = 0.1;
		jmin = -1;	// Initially no 3rd level is added
			if (Vmin >= 0 && (vd1k0*signo(vd1k0) >= Banda+10+200))// || *enable_desbalance == 1))	// Criteria not met -> Requires 3 level computation for one phase
			// f = -signo(vd1k0)*(Iu*abso(u_amax) + Iv*abso(u_bmax) + Iw*abso(u_cmax))
			{
				///// Reconsider value of x and what phase to add	/////
				// xmax
				if (-signo(vd1k0)*(Iu*(1-epsilon) - Iu*abso(u_amax)) + vec_min[0] < Vmin && abso(u_amax) < 1 - epsilon){	// Using xmax and adding 3rd level to phase a improves the cost function
					Vmin = -signo(vd1k0)*(Iu*(1-epsilon) - Iu*abso(u_amax)) + vec_min[0];	// Update cost function
					imin = 1;	// Update value of x
					jmin = 0;	// Phase a selected for the three levels addition
				}
				if(-signo(vd1k0)*(Iv*(1-epsilon) - Iv*abso(u_bmax)) + vec_min[0] < Vmin && abso(u_bmax) < 1 - epsilon){	// Using xmax and adding 3rd level to phase b improves the cost function
					Vmin = -signo(vd1k0)*(Iv*(1-epsilon) - Iv*abso(u_bmax)) + vec_min[0];	// Update cost function
					imin = 1;	// Update value of x
					jmin = 1;
				}
				if(-signo(vd1k0)*(Iw*(1-epsilon) - Iw*abso(u_cmax)) + vec_min[0] < Vmin && abso(u_cmax) < 1 - epsilon){	// Using xmax and adding 3rd level to phase c improves the cost function
					Vmin = -signo(vd1k0)*(Iw*(1-epsilon) - Iw*abso(u_cmax)) + vec_min[0];	// Update cost function
					imin = 1;	// Update value of x
					jmin = 2;
				}
				
				// xmin
				if(-signo(vd1k0)*(Iu*(1-epsilon) - Iu*abso(u_amin)) + vec_min[1] < Vmin && abso(u_amin) < 1 - epsilon){	// Using xmin and adding 3rd level to phase a improves the cost function
					Vmin = -signo(vd1k0)*(Iu*(1-epsilon) - Iu*abso(u_amin)) + vec_min[1];	// Update cost function
					imin = 2;	// Update value of x
					jmin = 0;
				}
				if(-signo(vd1k0)*(Iv*(1-epsilon) - Iv*abso(u_bmin)) + vec_min[1] < Vmin && abso(u_bmin) < 1 - epsilon){	// Using xmin and adding 3rd level to phase b improves the cost function
					Vmin = -signo(vd1k0)*(Iv*(1-epsilon) - Iv*abso(u_bmin)) + vec_min[1];	// Update cost function
					imin = 2;	// Update value of x
					jmin = 1;
				}
				if(-signo(vd1k0)*(Iw*(1-epsilon) - Iw*abso(u_cmin)) + vec_min[1] < Vmin && abso(u_cmin) < 1 - epsilon){	// Using xmin and adding 3rd level to phase c improves the cost function
					Vmin = -signo(vd1k0)*(Iw*(1-epsilon) - Iw*abso(u_cmin)) + vec_min[1];	// Update cost function
					imin = 2;	// Update value of x
					jmin = 2;
				}
				
				// xa
				if(-signo(vd1k0)*(Iu*(1-epsilon)) + vec_min[2] < Vmin && abso(alpha1 - alpha1) < 1 - epsilon && (0>u_amin) && (0<u_amax)){	// Using xa and adding 3rd level to phase a improves the cost function
					Vmin = -signo(vd1k0)*(Iu*(1-epsilon)) + vec_min[2];	// Update cost function
					imin = 3;	// Update value of x
					jmin = 0;
				}
				if(-signo(vd1k0)*(Iv*(1-epsilon) - Iv*abso(yx0)) + vec_min[2] < Vmin && abso(alpha2 - alpha1) < 1 - epsilon && (0>u_amin) && (0<u_amax)){	// Using xa and adding 3rd level to phase b improves the cost function
					Vmin = -signo(vd1k0)*(Iv*(1-epsilon) - Iv*abso(yx0)) + vec_min[2];	// Update cost function
					imin = 3;	// Update value of x
					jmin = 1;
				}
				if(-signo(vd1k0)*(Iw*(1-epsilon) - Iw*abso(zx0)) + vec_min[2] < Vmin && abso(alpha3 - alpha1) < 1 - epsilon && (0>u_amin) && (0<u_amax)){	// Using xa and adding 3rd level to phase c improves the cost function
					Vmin = -signo(vd1k0)*(Iw*(1-epsilon) - Iw*abso(zx0)) + vec_min[2];	// Update cost function
					imin = 3;	// Update value of x
					jmin = 2;
				}			
				
				// xb
				if(-signo(vd1k0)*(Iu*(1-epsilon) - Iu*abso(xy0)) + vec_min[3] < Vmin && abso(alpha1 - alpha2) < 1 - epsilon && (0>u_bmin) && (0<u_bmax)){	// Using xb and adding 3rd level to phase a improves the cost function
					Vmin = -signo(vd1k0)*(Iu*(1-epsilon) - Iu*abso(xy0)) + vec_min[3];	// Update cost function
					imin = 4;	// Update value of x
					jmin = 0;
				}
				if(-signo(vd1k0)*(Iv*(1-epsilon)) + vec_min[3] < Vmin && abso(alpha2 - alpha2) < 1 - epsilon && (0>u_bmin) && (0<u_bmax)){	// Using xb and adding 3rd level to phase b improves the cost function
					Vmin = -signo(vd1k0)*(Iv*(1-epsilon)) + vec_min[3];	// Update cost function
					imin = 4;	// Update value of x
					jmin = 1;
				}
				if(-signo(vd1k0)*(Iw*(1-epsilon) - Iw*abso(zy0)) + vec_min[3] < Vmin && abso(alpha3 - alpha2) < 1 - epsilon &&  (0>u_bmin) && (0<u_bmax)){	// Using xb and adding 3rd level to phase c improves the cost function
					Vmin = -signo(vd1k0)*(Iw*(1-epsilon) - Iw*abso(zy0)) + vec_min[3];	// Update cost function
					imin = 4;	// Update value of x
					jmin = 2;
				}	
				
				// xc
				if(-signo(vd1k0)*(Iu*(1-epsilon) - Iu*abso(xz0)) + vec_min[4] < Vmin && abso(alpha1 - alpha3) < 1 - epsilon &&  (0>u_cmin) && (0<u_cmax)){	// Using xc and adding 3rd level to phase a improves the cost function
					Vmin = -signo(vd1k0)*(Iu*(1-epsilon) - Iu*abso(xz0)) + vec_min[4];	// Update cost function
					imin = 5;	// Update value of x
					jmin = 0;
				}
				if(-signo(vd1k0)*(Iv*(1-epsilon) - Iv*abso(yz0)) + vec_min[4] < Vmin && abso(alpha2 - alpha3) < 1 - epsilon &&  (0>u_cmin) && (0<u_cmax)){	// Using xc and adding 3rd level to phase b improves the cost function
					Vmin = -signo(vd1k0)*(Iv*(1-epsilon) - Iv*abso(yz0)) + vec_min[4];	// Update cost function
					imin = 5;	// Update value of x
					jmin = 1;
				}
				if(-signo(vd1k0)*(Iw*(1-epsilon)) + vec_min[4] < Vmin && abso(alpha3 - alpha3) < 1 - epsilon &&  (0>u_cmin) && (0<u_cmax)){	// Using xc and adding 3rd level to phase c improves the cost function
					Vmin = -signo(vd1k0)*(Iw*(1-epsilon)) + vec_min[4];	// Update cost function
					imin = 5;	// Update value of x
					jmin = 2;
				}					
				
			}

			switch (imin)
            {
                case 1:
                    x = u_amax;
                    y = u_bmax;
                    z = u_cmax;
                    break;
                case 2:
                    x = u_amin;
                    y = u_bmin;
                    z = u_cmin;
                    break;
                case 3:
                    x = 0;
                    y = yx0;
                    z = zx0;
                    break;
                case 4:
                    x = xy0;
                    y = 0;
                    z = zy0;
                    break;
                case 5:
                    x = xz0;
                    y = yz0;
                    z = 0;
                    break;
            }

	     	da4 = maximo(x,0,0);
			da3 = 0.;
			da2 = 0.;
	     	da1 = maximo(-x,0,0);
	     	db4 = maximo(y,0,0);
			db3 = 0.;
			db2 = 0.;
	     	db1 = maximo(-y,0,0);
	     	dc4 = maximo(z,0,0);
			dc3 = 0.;
			dc2 = 0.;
	     	dc1 = maximo(-z,0,0);
			
			if (Vmin != Vmin_ant || jmin != -1){
				switch(jmin)
				{
					case 0:
						da1 = ((1-epsilon) - x)/2;
						da4 = ((1-epsilon) + x)/2;
						break;
					case 1:
						db1 = ((1-epsilon) - y)/2;
						db4 = ((1-epsilon) + y)/2;						
						break;
					case 2:
						dc1 = ((1-epsilon) - z)/2;
						dc4 = ((1-epsilon) + z)/2;							
						break;
				}
			}			
			
			/*** 			***/
			/* VSV MODIFICADO */
			/***			***/
			/*
			SVM_modificado(disp, RefAlpha, RefBeta, (*vc1r), (*vc2r), (*vc3r), (*vc4r));
			
			da1 = disp[0];
			da2 = 0.;
			da3 = 0.;
			da4 = disp[1];
			db1 = disp[2];
			db2 = 0.;
			db3 = 0.;
			db4 = disp[3];
			dc1 = disp[4];
			dc2 = 0.;
			dc3 = 0.;
			dc4 = disp[5];
			*/
			
			/** Corregir desbalance **/
			if (*enable_desbalance == 1)
				vd1k0 = vd1k0 - 45;
			/** Fin código correción desbalance **/	
			
			/***********************/
			/* Saturación inferior */
			/***********************/
            
            if(da1<0.) da1=0.;
			if(da2<0.) da2=0.;
			if(da3<0.) da3=0.;
			if(da4<0.) da4=0.;
			if(db1<0.) db1=0.;
			if(db2<0.) db2=0.;
			if(db3<0.) db3=0.;
			if(db4<0.) db4=0.;
			if(dc1<0.) dc1=0.;
			if(dc2<0.) dc2=0.;
			if(dc3<0.) dc3=0.;
			if(dc4<0.) dc4=0.;
			        
			/*********************************************************************/
            /* Ponderación de los dutis en caso de sobremodulación (suma dx > 1) */
			/*********************************************************************/
			
            if (da1 + da4 > 1)
            {
                Total = da1 + da4;
                
                da1 = da1/Total;
                da4 = da4/Total;
            }
            if (db1 + db4 > 1)
            {
                Total = db1 + db4;
                
                db1 = db1/Total;
                db4 = db4/Total;
            }
            if (dc1 + dc4 > 1)
            {
                Total = dc1 + dc4;
                
                dc1 = dc1/Total;
                dc4 = dc4/Total;
            }
			
			/****************************************/
			/* Cálculo de los tiempos de modulación */
			/****************************************/
		t1a = 0.;
		t2a = 0.;
		t3a = 0.;
		t4a = 0.;
		t5a = 0.;
		
		t1b = 0.;
		t2b = 0.;
		t3b = 0.;
		t4b = 0.;
		t5b = 0.;
		
		t1c = 0.;
		t2c = 0.;
		t3c = 0.;
		t4c = 0.;
		t5c = 0.;
			
			// Parte primera: cálculo de tiempos para la rama a
		if (estado != 0){
			
			t1a = da1;
			t2a = 0.0;
			t3a = 0.0;
			t4a = da4;
			t5a = 1.0 - da1 - da4;
			
			t1b = db1;
			t2b = 0.0;
			t3b = 0.0;
			t4b = db4;
			t5b = 1.0 - db1 - db4;
			
			t1c = dc1;
			t2c = 0.0;
			t3c = 0.0;
			t4c = dc4;
			t5c = 1.0 - dc1 - dc4;	
		}
		else
		{
			for (i=0; i<5; i++) //Poner tiempos a cero
			{
				*(tiempos_a+i) = 0.0;
				*(tiempos_b+i) = 0.0;
				*(tiempos_c+i) = 0.0;
			}
		}
			
			/*
			D:	   12345678
			t1x -> 00001111, En la FPGA : t1x -> tiempos_a[0]
			t2x -> 00011110, En la FPGA : t2x -> tiempos_a[1]
			t3x -> 01111000, En la FPGA : t3x -> tiempos_a[2]
			t4x -> 11110000, En la FPGA : t4x -> tiempos_a[3]
			t5x -> 00111100, En la FPGA : t5x -> tiempos_a[4]
			*/
			
			tiempos_a[0] = t1a; 
			tiempos_a[1] = t2a;
			tiempos_a[2] = t5a;
			tiempos_a[3] = t3a;
			tiempos_a[4] = t4a;
			
			tiempos_b[0] = t1b; 
			tiempos_b[1] = t2b;
			tiempos_b[2] = t5b;
			tiempos_b[3] = t3b;
			tiempos_b[4] = t4b;
			
			tiempos_c[0] = t1c; 
			tiempos_c[1] = t2c;
			tiempos_c[2] = t5c;
			tiempos_c[3] = t3c;
			tiempos_c[4] = t4c;
			
			/*
										________
									____		____
								____				____
							____						____
						____								____
			
			FPGA		|1	|2 	|3 	|4 	|	5 	|4 	|3 	|2 	|1 	|
			Control		|1 	|2 	|5 	|3 	|	4 	|3 	|5 	|2 	|1	|
			*/
			// FUJ FINAL 7
			//===================================================================================================
			
			/******************/
			//  PARA TESTEO   //
			/******************/
		/*
			tiempos_a[0] = señales1[iteracion];
			tiempos_a[1] = 0.2;
			tiempos_a[2] = señales2[iteracion];
			tiempos_a[3] = 0.2;
			tiempos_a[4] = señales3[iteracion];
			
			tiempos_b[0] = señales3[iteracion];
			tiempos_b[1] = 0.2;
			tiempos_b[2] = señales1[iteracion];
			tiempos_b[3] = 0.2;
			tiempos_b[4] = señales2[iteracion];

			tiempos_c[0] = señales2[iteracion];
			tiempos_c[1] = 0.2;
			tiempos_c[2] = señales3[iteracion];
			tiempos_c[3] = 0.2;
			tiempos_c[4] = señales1[iteracion];

			iteracion++;
			if (iteracion == 6)
				iteracion = 0;
		*/
		/////////////////////////////////////////	

		}
		else //Estado de Marcha pero no de Go
		{
			for (i=0; i<5; i++) //Poner tiempos a cero
			{
				*(tiempos_a+i) = 0.0;
				*(tiempos_b+i) = 0.0;
				*(tiempos_c+i) = 0.0;
			}			
		}
	} /************ FIN ESTADO DE MARCHA ************/

			/********************/
			/* Señales de Debug */
			/********************/
			
			/** UTILIZADAS POR LA INTERFAZ **/
            Debug[0] = enable_disp;
			Debug[1] = estado;
			Debug[2] = error;
			Debug[3] = RMS_Ia;
	        Debug[4] = RMS_Ib;
			Debug[5] = RMS_Ic;
			Debug[6] = RMS_Va;

			Debug[7] = RMS_Vb;
			Debug[8] = RMS_Vc;
			Debug[9] = cercania;
			/** FIN UTILIZADAS POR LA INTERFAZ **/
			
			/** SCOPE 4 **/
			Debug[10] = iar_ref;//iar_ref;//valfa;// 
			Debug[11] = ibr_ref;//ibr_ref;//x4alfa_p;// 
			Debug[12] = icr_ref;//icr_ref;//vbeta;// 
			Debug[13] = Pcond;//Pcond;//x4beta_p;//Vq;
			/** SCOPE 2 **/
			Debug[14] = u1;//Pcond;
			Debug[15] = u2;//valfa*ialfa_c + vbeta*ibeta_c;
			Debug[16] = out1[0];//u3;//-signo(vd1k0)*(Iu*(da1+da4)+Iv*(db1+db4)+Iw*(dc1+dc4));//u3;//qref;
			Debug[17] = out1[1];//u4;//vd1k0;//u4;//valfa*ibeta_c - vbeta*ialfa_c;
			/** SCOPE 5 **/
			Debug[18] = perror;//raiz2_3*dalfa1 + raiz1_3*dgamma1;//Vr;//da1;
			Debug[19] = p;//raiz2_3*dalfa2 + raiz1_3*dgamma2;//da1*(-Vcond/2) + da4*(Vcond/2);//da4;
			Debug[20] = q;//raiz2_3*dalfa1 + raiz1_3*dgamma1 + raiz2_3*dalfa2 + raiz1_3*dgamma2;//Vs;//1 - da1 - da4;
			Debug[21] = lvdcref;//db1*(-Vcond/2) + db4*(Vcond/2);//lvdcref;
			
			Debug[22] = idh_p[0];
			Debug[23] = out1_LP[0];
			Debug[24] = out1_LP2[0];
			Debug[25] = PIdh_p[0];
			Debug[26] = ialh_p[0];
			Debug[27] = ialh[0];
			Debug[28] = ia_ref;
			
			//RefAlpha = ORS1AB.Calpha + ORS2AB.Calpha + Kp*(perror*valfa - qerror*vbeta) + IntPk0*valfa - IntQk0*vbeta;//valfa
			//RefBeta = ORS1AB.Cbeta + ORS2AB.Cbeta + Kp*(perror*vbeta + qerror*valfa) + IntPk0*vbeta + IntQk0*valfa;//vbeta
			
}