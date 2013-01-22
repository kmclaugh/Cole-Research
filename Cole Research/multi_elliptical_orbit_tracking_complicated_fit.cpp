/*************************************************************************/
/*    In this program, we have used in the past (document more here):    */
/*         an adaptive step Bulirsh-Stoer Method and                     */ 
/* 		   a 4th order adaptive step Runge Kutta code.                   */
/*    We used "least square fit" method to plot out "epsilon, p          */
/*    semimajor, semiminor, orbit and radius v.s. time" into files       */
/*    1/r=(1/(epsilon*p)-cos(theta-theta0)/p                             */
/*    use "mid-point method' to track theta in terms of t.               */
/*************************************************************************/


#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <errno.h>
#include "nr.h"
#include "nrutil.h"
#include <unistd.h>

#define PI 3.1415926535897932
#define LC 29979245800.0
#define m 9.1091e-28
#define e 4.80298e-10
#define EPS 0.0000005
#define N 6
#define NUSE 7
#define IMAXX 11
#define SHRINK 0.95
#define GROW 1.2
#define TINY 1.0e-30



double FLorx,FLory;
double Radiationx, Radiationy, Radiationz;
double theta;

void derivs(double x,double y[],double dydx[]);
void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,double yscal[], double *hdid, double *hnext, void (*derivs)(double, double [], double []));
void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep, double yout[], void (*derivs)(double, double[], double[]));
void rzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv, int nuse);
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,double hmin, int *nok, int *nbad, void (*derivs)(double, double [], double []),void (*rkqs)(double [], double [], int, double *, double, double, double [],double *, double *, void (*)(double, double [], double [])));
double search(double z);
double myfun(double thetavalue);

double dxsav, *xp, **yp;
int kount, kmax, Fold, Fnew, theta0_old, F1,F2,F3,F4, Nold, Nnew;
int Cmass, Crad, Cplane, CounterPrint;
double **d,*x;
double radius,smajor;
/***********************************************************************************************
 * Note: the reason why I have A1 and amplitude1, A2 and amplitude2, etc., is that it is then
 * easier to turn the waves on and off.  A1, A2, ..., A11 are the input amplitudes.
 * amplitude1, ..., amplitude11 are the amplitudes in the code at each instant of time.
 ***********************************************************************************************/
double A1,  amplitude1,  frequency1,  sm1,  phase1,  ec1,  tp1;
double A2,  amplitude2,  frequency2,  sm2,  phase2,  ec2,  tp2;
double A3,  amplitude3,  frequency3,  sm3,  phase3,  ec3,  tp3;
double A4,  amplitude4,  frequency4,  sm4,  phase4,  ec4,  tp4;
double A5,  amplitude5,  frequency5,  sm5,  phase5,  ec5,  tp5;
double A6,  amplitude6,  frequency6,  sm6,  phase6,  ec6,  tp6;
double A7,  amplitude7,  frequency7,  sm7,  phase7,  ec7,  tp7;
double A8,  amplitude8,  frequency8,  sm8,  phase8,  ec8,  tp8;
double A9,  amplitude9,  frequency9,  sm9,  phase9,  ec9,  tp9;
double A10, amplitude10, frequency10, sm10, phase10, ec10, tp10;
double A11, amplitude11, frequency11, sm11, phase11, ec11, tp11;
double A12, amplitude12, frequency12, sm12, phase12, ec12, tp12;
double A13, amplitude13, frequency13, sm13, phase13, ec13, tp13;
double A14, amplitude14, frequency14, sm14, phase14, ec14, tp14;
double A15, amplitude15, frequency15, sm15, phase15, ec15, tp15;
double A16, amplitude16, frequency16, sm16, phase16, ec16, tp16;
double A17, amplitude17, frequency17, sm17, phase17, ec17, tp17;
double A18, amplitude18, frequency18, sm18, phase18, ec18, tp18;
double A19, amplitude19, frequency19, sm19, phase19, ec19, tp19;
double A20, amplitude20, frequency20, sm20, phase20, ec20, tp20;
double A21, amplitude21, frequency21, sm21, phase21, ec21, tp21;
double A22, amplitude22, frequency22, sm22, phase22, ec22, tp22;
double A23, amplitude23, frequency23, sm23, phase23, ec23, tp23;
double A24, amplitude24, frequency24, sm24, phase24, ec24, tp24;
double A25, amplitude25, frequency25, sm25, phase25, ec25, tp25;
double A26, amplitude26, frequency26, sm26, phase26, ec26, tp26;
double A27, amplitude27, frequency27, sm27, phase27, ec27, tp27;
double A28, amplitude28, frequency28, sm28, phase28, ec28, tp28;
double A29, amplitude29, frequency29, sm29, phase29, ec29, tp29;
double A30, amplitude30, frequency30, sm30, phase30, ec30, tp30;
double A31, amplitude31, frequency31, sm31, phase31, ec31, tp31;
/*******************************************************************************************************************************
 * Mfrequency and Mamplitude are part of the relativistic code, that need to be gone over again in later versions of this code.
 *******************************************************************************************************************************/
double Mfrequency, Mamplitude;
double k=pow(e,2);
double duration, T;
double eps,epsilon;
double epsilonp;			//epsilon*P = smajor*(1-epsilon^2)
double asint,acost,sintheta0,costheta0,theta0temp;

int main ()
{

/***************************************************************************************************************************************************
 * The key controlling parameters for the simulation are as follows:
 *
 * CounterPrint:  A typical value is 28.  This should be an integer 1, 2, ...
 *          It specifies the number of points caluated before one is printed to file.
 *          Specifying, for example, 28, means that after 28 time steps have been computed, then things like radius, energy, semimajor axis,
 *          etc., are written to file.  Having a large value of CounterPrint, keeps the file sizes from being too large.  However, 
 *          when "zooming" in, then one sees the "blocky" look of the plots.  If detail is needed, then CounterPrint shouuld be like 1, 2, 3, small, 
 *          Usually, when running subharmonic resonances, with precision eps=exp(-22), then CounterPrint=28 is fairly good.
 *
 * Cmass:   Set this equal to 0 for nonrelativistic simulation
 *          Set it equal to 1 for relativistic simulation, although I think this code needs updating.
 *
 * Crad:    Set equal to 1 for having radiation reaction turned on (usual case is on).
 *          "          " 0 to turn radiation reaction off.
 *
 * Cplane:  Set equal to 1 to have applied plane wave radiation on.  (Another way to do this would be to make A1=A1=0.)
 *          "          " 0 to turn plane waves off.
 *
 * eps:     This parameter species the precision of the adaptive step algorithm.  A value of exp((double)-22) is usually fine.
 *          exp((double)-20) may work fine sometimes, but it could be marginal.  exp((double)-28) gets beyond what the 16
 *          digits can handle. 
 *
 * duration:The length of time, in electron time, for the simulation to run. (sec)
 *
 * smajor:  Initial starting value of the semimajor axis (cm).  If epsilon=0, then this is just the radius at t=0.
 *
 * epsilon: Eccentricty of orbit at t=0.  Note: 0<=epsilon<1.  When epsilon=0, then the initial orbit is a circle.  When
 *          epsilon gets close to 1, such as 0.9, then the orbit is highly elliptical.  Note that
 * 
 * Rmax:    This is how long we allow the radius (in cm) to become before the program stops.  Above the value of
 *          Rmax, the program stops and assumes that the electron has "ionized".  A default being used is 1 micron = 1e-4 cm.
 *
 * The present code allows for 11 elliptically polarized plane waves.  
 * The amplitudes (A1,...,A31), 
 * the phases (phase1,...,phase31), 
 * the frequencies (use values sm1,...,sm31, as described below, to set these frequencies), 
 * the eccentricities (ec1,...,ec11), 
 * and the times that the waves turn on (tp1,...,tp31),
 * can all be individually controlled via the paramaters about to be described.
 *
 * A1,...,A31:
 *          Amplitude of 1st, 2nd, ..., 31th set of elliptically polarized plane waves, all in (statvolt/cm).
 *          These values are the amplitudes of the Ex vectors of the corresponding plane waves.  
 *          Consider the 3rd plane wave case, for example.  A3 is the amplitude of the Ex vector of this 
 *          3rd elliptical plane wave.  The amplitude of the Ey vector for this 3rd wave is A3*sqrt(1-ec3*ec3).
 *          Thus, when ec3=0, which is when we have a circularly polarized plane wave for the 3rd wave, then
 *          A3 is the amplitude for both the x and y vector components of the electric field for this
 *          3rd plane wave.  Otherwise, the amplitude of Ey is less than that of Ex, for this 3rd wave.
 *          Note, in the future I could alter this formulation, so that the ellipse formed by the Ex, Ey vectors of each
 *          of these elliptically polarized plane waves, could be tilted.  I.e., right now, each ellipse always has its
 *          semimajor axis along the x direction, which is why the amplitude of the Ex vector is greater or equal to the
 *          amplitude of the Ey vector.   To tilt one of these ellipses, I would just add a phase in, and it should be
 *          general enough to be different for each elliptically polarized plane wave.
 *          Finally, note that if ec3=1, then the amplitude of the Ex vector is A3, while that of the Ey vector is 0.  In
 *          that case we have linearly polarized light for the 3rd wave, with the linear polarization along the x direction.
 * 
 * phase1,...,phase31:
 *          Phase of 1st 1st, 2nd, ..., 31th set of elliptically polarized plane waves, all in radians,
 *          in relation to the initial position of the electron.  More specifically, 
 *          when for example, phase4=0, then that means that at t=0, (-e) times the 
 *          electric field of the 4th E&M wave. points along the y direction, since 
 *          this is the direction that the electron starts along.  See my papers and diagrams on this
 *			for a better discussion.  These phases are in radians. Hence, use values like Pi/4, etc.
 *
 * sm1,...,sm31
 *          These pertain to the 1st, 2nd, ..., 31th set of elliptically polarized plane waves.  Units are in cm.
 *          These values provide the angular frequencies of the waves, but in a bit of an odd way.
 *          Consider sm2, for example.  The angular frequency of the 2nd elliptically polarized wave is calculated by
 *          frequency2=sqrt( k/ (m*pow(sm2,3)) ) .  This calulation is done later in the program.
 *          Now, why this is a convenient way of doing this is that if an electron
 *          travelled around the nucleus in an elliptical orbit with a semimajor axis equal to sm2. then frequency2 would be precisely
 *          equal to 2*Pi / Period, where Period is the period of that orbit, or, equivalently, to 2*Pi times the orbital frequency.
 *          So, for example, if you want to insert the second elliptical plane wave with an angular frequency equal to the 2*Pi times the 
 *          orbital frequency of an electron in an elliptical (or circular) orbit with semimajor axis equal to, say, 0.5 Angstroms 
 *          (close to the Bohr orbit), you would simply enter sm2=0.5e-8 down below.
 *
 * ec1,...,ec31
 *          These pertain to the 1st, 2nd, ..., 31th set of elliptically polarized plane waves.  Units are dimensionless.
 *          All values should be greater than or equal to zero, and less than one, so   0 <= ecn < 1 .
 *          These dimensionless parameters are the eccentricities of the corresponding elliptically polarized plane
 *          waves that can be applied here.  
 *
 * tp1,...,tp31
 *          These pertain to the 1st, 2nd, ..., 31th set of elliptically polarized plane waves.  Units are in seconds.
 *          tp1 is the time that the first elliptically polarized plane wave is turned on.
 *          tp2 is the time that the second elliptically polarized plane wave is turned on.
 *          Etc.  (Yes, very clever notation indeed.)
 **************************************************************************************************************************************************/

CounterPrint = 6;          //Counter for # points before printing for radius
Cmass = 0;                 //Leave as 0 (nonrelativistic) for now.  Later will check to make sure 1 (relativistic) works.
Crad = 1;                  //Radiation reaction.  1 to include.  0 not to include.  Default: 1
Cplane = 1;		 	       //Use of plane waves.  1 for plane waves.  0 for no plane waves.  Default: 1
eps = exp((double)-23);	   //Controls the numerical precision.

duration=1e-9;             //How long, in seconds for the electron, for the simulation to run.  Thus, 1e-12 would mean
                           //the simulation will last for 1e-12 sec in terms of the electron's trajectory.  
                           //For full decay, a value like 1e-8 usually works.  For a few minutes of run time, use something
                           //like 3e-12.
/***************************************************************************************************************************************************
 * The folowing input parameters I am writing in an odd way, as one can then copy them readily and paste them into headings in tecplot plots.
 * The present code provides for 11 elliptically polarized plane waves.  This number can be readily changed in later program versions.
 ***************************************************************************************************************************************************/
 
epsilon=0.0;	 // dimensionless, 0<=epsilon<1 , initial eccentricity of orbit		  
smajor=((0.25e-8)* pow(2.0,2.0/3.0)) + 0.03e-8;  //(0.25e-8)*pow(2,2/3);  //(cm)

A1=2000.0;       // (statvolt/cm)  This is the amplitude of the Ex vector of the first plane wave.
A2=0.0; A3=0.0; A4=0.0; A5=0.0;  A6=0.0;
A7=0.0; A8=0.0; A9=0.0; A10=0.0; A11=0.0;
A12=0.0; A13=0.0; A14=0.0; A15=0.0; A16=0.0,

A17=0.0; A18=0.0; A19=0.0; A20=0.0; A21=0.0;
A22=0.0; A23=0.0; A24=0.0; A25=0.0; A26=0.0;
A27=0.0; A28=0.0; A29=0.0; A30=0.0; A31=0.0; 

sm1=0.25e-8;     // (cm) (but used to compute angular frequency of wave via:  frequency1=sqrt( k/ (m*pow(sm1,3)) )  

sm2=0.25e-8; sm3=0.24804e-8; sm4=0.24818e-8; sm5=0.24832e-8;  sm6=0.24846e-8; 
sm7=0.24860e-8; sm8=0.24874e-8; sm9=0.24888e-8; sm10=0.24902e-8; sm11=0.24916e-8; 
sm12=0.24930e-8; sm13=0.24944e-8; sm14=0.24958e-8; sm15=0.24972e-8; sm16=0.24986e-8,

sm17=0.25014e-8; sm18=0.25028e-8; sm19=0.25042e-8; sm20=0.25056e-8; sm21=0.25070e-8;
sm22=0.25084e-8; sm23=0.25098e-8; sm24=0.25112e-8; sm25=0.25126e-8; sm26=0.25140e-8;
sm27=0.25154e-8; sm28=0.25168e-8; sm29=0.25182e-8; sm30=0.25196e-8; sm31=0.25210e-8;

ec1=0.0;         // dimensionless (eccentricity of plane wave #1)
ec2=0.0; ec3=0.0; ec4=0.0; ec5=0.0; ec6=0.0; 
ec7=0.0; ec8=0.0; ec9=0.0; ec10=0.0; ec11=0.0; 
ec12=0.0; ec13=0.0; ec14=0.0; ec15=0.0; ec16=0.0,
ec17=0.0; ec18=0.0; ec19=0.0; ec20=0.0; ec21=0.0;
ec22=0.0; ec23=0.0; ec24=0.0; ec25=0.0; ec26=0.0;
ec27=0.0; ec28=0.0; ec29=0.0; ec30=0.0; ec31=0.0;

phase1=0.0;     // radians, so phase1 should be between -Pi and Pi.  Phase of 1st wave in relation to initial velocity.  See above.
phase2=0.1; phase3=0.8; phase4=0.3; phase5=1.45;  phase6=2.31; 
phase7=3.7; phase8=5.3; phase9=9.1; phase10=2.4; phase11=4.34; 
phase12=1.112; phase13=7.24; phase14=8.3; phase15=2.6; phase16=3.77,
phase17=2.12; phase18=2.43; phase19=3.765; phase20=3.21; phase21=2.79;
phase22=1.99; phase23=1.234; phase24=8.3; phase25=9.4; phase26=10.2;
phase27=3.2345; phase28=3.21; phase29=2.789; phase30=4.56; phase31=6.543;

/* phase2=0.0; phase3=0.0; phase4=0.0; phase5=0.0;  phase6=0.0; 
phase7=0.0; phase8=0.0; phase9=0.0; phase10=0.0; phase11=0.0; 
phase12=0.0; phase13=0.0; phase14=0.0; phase15=0.0; phase16=0.0,
phase17=0.0; phase18=0.0; phase19=0.0; phase20=0.0; phase21=0.0;
phase22=0.0; phase23=0.0; phase24=0.0; phase25=0.0; phase26=0.0;
phase27=0.0; phase28=0.0; phase29=0.0; phase30=0.0; phase31=0.0; */

tp1=0.0;     // (sec) Time at which the first plane wave turns on.
tp2=0.0; tp3=0.0; tp4=0.0; tp5=0.0; tp6=0.0; 
tp7=0.0; tp8=0.0; tp9=0.0; tp10=0.0; tp11=1.9e-12; 
tp12=0.0; tp13=0.0; tp14=0.0; tp15=0.0; tp16=0.0,
tp17=0.0; tp18=0.0; tp19=0.0; tp20=0.0; tp21=0.0;
tp22=0.0; tp23=0.0; tp24=0.0; tp25=0.0; tp26=0.0;
tp27=0.0; tp28=0.0; tp29=0.0; tp30=0.0; tp31=0.0;

/***************************************************************************************************************************************************
 * NOTE: We once created files using a method Lenny started that used the following for designating a file:
 * conditions = "mass0_rad1_plane1_a500_p0_smjr2.0e-8_rcp0.5e-8_eps-20_t1.6e-9_ecc0.0";
 ***************************************************************************************************************************************************/

theta0_old=0.0;
        Fold=0;
        Nold=0;
	amplitude1 = A1;
	radius = smajor*(1+epsilon);	//Starting radius of electron.  Note that the semimajor axis has length smajor at
	                                //t=0, so the distance from the origin (where the nucleus is, is given by smajor*(1+epsilon).
	epsilonp = smajor*(1-pow(epsilon,2));
	T=2.0*PI*sqrt(m)*pow((epsilonp),1.5)/(e*pow((1-pow(epsilon,2)),1.5));
	
	int nbad, nok;
	int Ceps, Cdur,Cpha;
	double h1, hmin, x1, x2, *ystart;
	double Py;
	

	time_t start, end;

	frequency1 =sqrt( k/ (m*pow(sm1,3)) );
    frequency2 =sqrt( k/ (m*pow(sm2,3)) );
    frequency3 =sqrt( k/ (m*pow(sm3,3)) );
    frequency4 =sqrt( k/ (m*pow(sm4,3)) );
    frequency5 =sqrt( k/ (m*pow(sm5,3)) );
    frequency6 =sqrt( k/ (m*pow(sm6,3)) );
    frequency7 =sqrt( k/ (m*pow(sm7,3)) );
    frequency8 =sqrt( k/ (m*pow(sm8,3)) );
    frequency9 =sqrt( k/ (m*pow(sm9,3)) );
    frequency10=sqrt( k/ (m*pow(sm10,3)) );
    frequency11=sqrt( k/ (m*pow(sm11,3)) );
    frequency12=sqrt( k/ (m*pow(sm12,3)) );
    frequency13=sqrt( k/ (m*pow(sm13,3)) );
    frequency14=sqrt( k/ (m*pow(sm14,3)) );
    frequency15=sqrt( k/ (m*pow(sm15,3)) );
    frequency16=sqrt( k/ (m*pow(sm16,3)) );
    frequency17=sqrt( k/ (m*pow(sm17,3)) );
    frequency18=sqrt( k/ (m*pow(sm18,3)) );
    frequency19=sqrt( k/ (m*pow(sm19,3)) );
    frequency20=sqrt( k/ (m*pow(sm20,3)) );
    frequency21=sqrt( k/ (m*pow(sm21,3)) );
    frequency22=sqrt( k/ (m*pow(sm22,3)) );
    frequency23=sqrt( k/ (m*pow(sm23,3)) );
    frequency24=sqrt( k/ (m*pow(sm24,3)) );
    frequency25=sqrt( k/ (m*pow(sm25,3)) );
    frequency26=sqrt( k/ (m*pow(sm26,3)) );
    frequency27=sqrt( k/ (m*pow(sm27,3)) );
    frequency28=sqrt( k/ (m*pow(sm28,3)) );
    frequency29=sqrt( k/ (m*pow(sm29,3)) );
    frequency30=sqrt( k/ (m*pow(sm30,3)) );
    frequency31=sqrt( k/ (m*pow(sm31,3)) );

    
	if (Cmass==1)
	{
		Py=sqrt(m*k/radius)*sqrt(      sqrt(  1+pow( k/(2*LC*LC*m*radius) , 2)   ) - k/(2*LC*LC*m*radius)  );
		printf("The program will run under the 'relativistic' condition.\n\n\n");   
	}	
	if (Cmass==0)
    {

		Py=pow( (m*k*(1-epsilon)/(radius)),0.5);	
		printf("The program will run under the 'nonrelativistic' condition.\n\n\n");
	}


	if (Crad==1)
		printf("The program will include radiation reaction.\n\n\n");   
	if (Crad==0)
		printf("The program will not include radiation reaction.\n\n\n");   
	
	
	if (Cplane==1)
		printf("The program will include plane waves.\n\n\n");   
	if (Cplane==0)
		printf("The program will not include plane waves.\n\n\n");   
	


  /*      -----------------Driven programming for the subroutine of the algorithm--------------      */
  
	ystart=vector(1,N);
	xp=vector(1,200);
	yp=matrix(1,10,1,200);
	x1=0;
	x2=duration;

/*************************************************************************************************************************************
 * In this code, the electron starts at x=radius=a(1+epsilon), y=0, and it moves upward (vx=0, vy=Py/m) at such a speed as to make it
 * travel in an elliptical orbit with eccentricity epsilon, if no other forces except the Coulombic one was to act.
 *************************************************************************************************************************************/

	ystart[1]=0;   //vstart[1]=	(k*m/r)^0.5;	initial momentum in x;
	ystart[2]=Py;	//initial momentum in y;
	ystart[3]=0;	//initial momentum in z;
	ystart[4]=radius;	//initial position in x;
	ystart[5]=0;	//initial position in y;
	ystart[6]=0;	//initial position in z;
	
	start = time(NULL);      
	printf("\nPlease wait, program is now running...\n");
	
	h1=0.1;
	hmin=0.0;
	kmax=100;
	dxsav=(x2-x1)/20.0;
	
    odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);

	printf("\nProgram executed successfully! Thank you for your patience.\n");
	end = time(NULL);      
	printf("\nThe duration of the program running was: %.0f seconds\n",difftime(end,start));
	
	free_matrix(yp,1,10,1,200);
	free_vector(xp,1,200);
	free_vector(ystart,1,N);
}


  /*      -----------------Subroutine for definition of the differential equations--------------      */
  
  void derivs(double x,double y[],double dydx[])
  {
	  double r3=pow( (y[4]*y[4]+y[5]*y[5]+y[6]*y[6]), 1.5);
	  double r5=pow( (y[4]*y[4]+y[5]*y[5]+y[6]*y[6]), 2.5);
	  double p2=y[1]*y[1]+y[2]*y[2]+y[3]*y[3];
	  double K=p2/(m*m*LC*LC);

	  if (Cplane==1)
	  {
		  theta = search(x/T - (int)(x/T));
		    
		    // The following lines are my concocted way of turning on a plane wave at a certain time.  Make this much better!
		//***************************************************
		//  if (x>=tp2){
		//		 amplitude1=0;
		//		 amplitude2=A2;}
		//  else {
		//		 amplitude1=A1;
		//	     amplitude2=0;
	    //**************************************************
//* In the following code, I have switched things so that one can turn on the different waves at different times.
		  if (x<tp1){
				 amplitude1=0.0;}
		     else {
			     amplitude1=A1;
		     }
		  if (x<tp2){
				 amplitude2=0.0;}
		     else {
			     amplitude2=A2;
		     }
		  if (x<tp3){
				 amplitude3=0.0;}
		     else {
			     amplitude3=A3;
		     }
		  if (x<tp4){
				 amplitude4=0.0;}
		     else {
			     amplitude4=A4;
		     }
		  if (x<tp5){
				 amplitude5=0.0;}
		     else {
			     amplitude5=A5;
		     }
		  if (x<tp6){
				 amplitude6=0.0;}
		     else {
			     amplitude6=A6;
		     }
		  if (x<tp7){
				 amplitude7=0.0;}
		     else {
			     amplitude7=A7;
		     }
		  if (x<tp8){
				 amplitude8=0.0;}
		     else {
			     amplitude8=A8;
		     }
		  if (x<tp9){
				 amplitude9=0.0;}
		     else {
			     amplitude9=A9;
		     }
		  if (x<tp10){
				 amplitude10=0.0;}
		     else {
			     amplitude10=A10;
		     }
		  if (x<tp11){
				 amplitude11=0.0;}
		     else {
			     amplitude11=A11;
		     }
		  if (x<tp12){
				 amplitude12=0.0;}
		     else {
			     amplitude12=A12;
		     }
		  if (x<tp13){
				 amplitude13=0.0;}
		     else {
			     amplitude13=A13;
		     }
		  if (x<tp14){
				 amplitude14=0.0;}
		     else {
			     amplitude14=A14;
		     }
		  if (x<tp15){
				 amplitude15=0.0;}
		     else {
			     amplitude15=A15;
		     }
		  if (x<tp16){
				 amplitude16=0.0;}
		     else {
			     amplitude16=A16;
		     }
		  if (x<tp17){
				 amplitude17=0.0;}
		     else {
			     amplitude17=A17;
		     }
		  if (x<tp18){
				 amplitude18=0.0;}
		     else {
			     amplitude18=A18;
		     }
		  if (x<tp19){
				 amplitude19=0.0;}
		     else {
			     amplitude19=A19;
		     }
		  if (x<tp20){
				 amplitude20=0.0;}
		     else {
			     amplitude20=A20;
		     }
		  if (x<tp21){
				 amplitude21=0.0;}
		     else {
			     amplitude21=A21;
		     }
		  if (x<tp22){
				 amplitude22=0.0;}
		     else {
			     amplitude22=A22;
		     }
		  if (x<tp23){
				 amplitude23=0.0;}
		     else {
			     amplitude23=A23;
		     }
		  if (x<tp24){
				 amplitude24=0.0;}
		     else {
			     amplitude24=A24;
		     }
		  if (x<tp25){
				 amplitude25=0.0;}
		     else {
			     amplitude25=A25;
		     }
		  if (x<tp26){
				 amplitude26=0.0;}
		     else {
			     amplitude26=A26;
		     }
		  if (x<tp27){
				 amplitude27=0.0;}
		     else {
			     amplitude27=A27;
		     }
		  if (x<tp28){
				 amplitude28=0.0;}
		     else {
			     amplitude28=A28;
		     }
		  if (x<tp29){
				 amplitude29=0.0;}
		     else {
			     amplitude29=A29;
		     }
		  if (x<tp30){
				 amplitude30=0.0;}
		     else {
			     amplitude30=A30;
		     }
		  if (x<tp31){
				 amplitude31=0.0;}
		     else {
			     amplitude31=A31;
		     }
// **********************************************************************
/*			Ex= amplitude1*
				( cos(y[6]*frequency1/LC+frequency1*x-PI/2+phase1)
			   + (y[3]/(LC*m))*cos(y[6]*frequency1/LC+frequency1*x-PI/2+phase1) ) ;


			Ey= amplitude1*
				 ( - cos(y[6]*frequency1/LC+frequency1*x+phase1)
				   - (y[3]/(LC*m))*cos(y[6]*frequency1/LC+frequency1*x+phase1) ) ;
//*******************************************************************************************/
		     FLorx = (-e)*(  amplitude1 *(1+(y[3]/(LC*m)))*(cos(y[6]*frequency1 /LC+frequency1*x-PI/2+phase1) )
			               + amplitude2 *(1+(y[3]/(LC*m)))*(cos(y[6]*frequency2 /LC+frequency2*x-PI/2+phase2) )
			               + amplitude3 *(1+(y[3]/(LC*m)))*(cos(y[6]*frequency3 /LC+frequency3*x-PI/2+phase3) )
			               + amplitude4 *(1+(y[3]/(LC*m)))*(cos(y[6]*frequency4 /LC+frequency4*x-PI/2+phase4) )
			               + amplitude5 *(1+(y[3]/(LC*m)))*(cos(y[6]*frequency5 /LC+frequency5*x-PI/2+phase5) )
			               + amplitude6 *(1+(y[3]/(LC*m)))*(cos(y[6]*frequency6 /LC+frequency6*x-PI/2+phase6) )
			               + amplitude7 *(1+(y[3]/(LC*m)))*(cos(y[6]*frequency7 /LC+frequency7*x-PI/2+phase7) )
			               + amplitude8 *(1+(y[3]/(LC*m)))*(cos(y[6]*frequency8 /LC+frequency8*x-PI/2+phase8) )
			               + amplitude9 *(1+(y[3]/(LC*m)))*(cos(y[6]*frequency9 /LC+frequency9*x-PI/2+phase9) )
			               + amplitude10*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency10/LC+frequency10*x-PI/2+phase10))
			               + amplitude11*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency11/LC+frequency11*x-PI/2+phase11))
			               + amplitude12*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency12/LC+frequency12*x-PI/2+phase12))
			               + amplitude13*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency13/LC+frequency13*x-PI/2+phase13))
			               + amplitude14*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency14/LC+frequency14*x-PI/2+phase14))
			               + amplitude15*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency15/LC+frequency15*x-PI/2+phase15))
			               + amplitude16*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency16/LC+frequency16*x-PI/2+phase16))
			               + amplitude17*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency17/LC+frequency17*x-PI/2+phase17))
			               + amplitude18*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency18/LC+frequency18*x-PI/2+phase18))
			               + amplitude19*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency19/LC+frequency19*x-PI/2+phase19))
			               + amplitude20*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency20/LC+frequency20*x-PI/2+phase20))
			               + amplitude21*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency21/LC+frequency21*x-PI/2+phase21))
			               + amplitude22*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency22/LC+frequency22*x-PI/2+phase22))
			               + amplitude23*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency23/LC+frequency23*x-PI/2+phase23))
			               + amplitude24*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency24/LC+frequency24*x-PI/2+phase24))
			               + amplitude25*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency25/LC+frequency25*x-PI/2+phase25))
			               + amplitude26*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency26/LC+frequency26*x-PI/2+phase26))
			               + amplitude27*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency27/LC+frequency27*x-PI/2+phase27))
			               + amplitude28*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency28/LC+frequency28*x-PI/2+phase28))
			               + amplitude29*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency29/LC+frequency29*x-PI/2+phase29))
			               + amplitude30*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency30/LC+frequency30*x-PI/2+phase30))
			               + amplitude31*(1+(y[3]/(LC*m)))*(cos(y[6]*frequency31/LC+frequency31*x-PI/2+phase31))
						   );
		     FLory = (-e)*(  amplitude1*pow(1-ec1*ec1,0.5)   *(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency1/LC  +frequency1*x +phase1) )
                           + amplitude2*pow(1-ec2*ec2,0.5)   *(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency2/LC  +frequency2*x +phase2) )
                           + amplitude3*pow(1-ec3*ec3,0.5)   *(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency3/LC  +frequency3*x +phase3) )
                           + amplitude4*pow(1-ec4*ec4,0.5)   *(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency4/LC  +frequency4*x +phase4) )
                           + amplitude5*pow(1-ec5*ec5,0.5)   *(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency5/LC  +frequency5*x +phase5) )
                           + amplitude6*pow(1-ec6*ec6,0.5)   *(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency6/LC  +frequency6*x +phase6) )
                           + amplitude7*pow(1-ec7*ec7,0.5)   *(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency7/LC  +frequency7*x +phase7) )
                           + amplitude8*pow(1-ec8*ec8,0.5)   *(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency8/LC  +frequency8*x +phase8) )
                           + amplitude9*pow(1-ec9*ec9,0.5)   *(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency9/LC  +frequency9*x +phase9) )
                           + amplitude10*pow(1-ec10*ec10,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency10/LC+frequency10*x+phase10) )
                           + amplitude11*pow(1-ec11*ec11,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency11/LC+frequency11*x+phase11) )
                           + amplitude12*pow(1-ec12*ec12,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency12/LC+frequency12*x+phase12) )
                           + amplitude13*pow(1-ec13*ec13,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency13/LC+frequency13*x+phase13) )
                           + amplitude14*pow(1-ec14*ec14,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency14/LC+frequency14*x+phase14) )
                           + amplitude15*pow(1-ec15*ec15,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency15/LC+frequency15*x+phase15) )
                           + amplitude16*pow(1-ec16*ec16,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency16/LC+frequency16*x+phase16) )
                           + amplitude17*pow(1-ec17*ec17,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency17/LC+frequency17*x+phase17) )
                           + amplitude18*pow(1-ec18*ec18,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency18/LC+frequency18*x+phase18) )
                           + amplitude19*pow(1-ec19*ec19,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency19/LC+frequency19*x+phase19) )
                           + amplitude20*pow(1-ec20*ec20,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency20/LC+frequency20*x+phase20) )
                           + amplitude21*pow(1-ec21*ec21,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency21/LC+frequency21*x+phase21) )
                           + amplitude22*pow(1-ec22*ec22,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency22/LC+frequency22*x+phase22) )
                           + amplitude23*pow(1-ec23*ec23,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency23/LC+frequency23*x+phase23) )
                           + amplitude24*pow(1-ec24*ec24,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency24/LC+frequency24*x+phase24) )
                           + amplitude25*pow(1-ec25*ec25,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency25/LC+frequency25*x+phase25) )
                           + amplitude26*pow(1-ec26*ec26,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency26/LC+frequency26*x+phase26) )
                           + amplitude27*pow(1-ec27*ec27,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency27/LC+frequency27*x+phase27) )
                           + amplitude28*pow(1-ec28*ec28,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency28/LC+frequency28*x+phase28) )
                           + amplitude29*pow(1-ec29*ec29,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency29/LC+frequency29*x+phase29) )
                           + amplitude30*pow(1-ec30*ec30,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency30/LC+frequency30*x+phase30) )
                           + amplitude31*pow(1-ec31*ec31,0.5)*(1+(y[3]/(LC*m)))*( - cos(y[6]*frequency31/LC+frequency31*x+phase31) )
						   );
//******************************************************************************************/
	  }

	  Radiationx=(2.0/3)*(k*k)/(m*LC*LC*LC)*( y[1]/m/r3 -3*y[4]*(y[4]*y[1]/m+y[5]*y[2]/m+y[6]*y[3]/m )/r5  );
	  Radiationy=(2.0/3)*(k*k)/(m*LC*LC*LC)*( y[2]/m/r3 -3*y[5]*(y[4]*y[1]/m+y[5]*y[2]/m+y[6]*y[3]/m )/r5  );
	  Radiationz=(2.0/3)*(k*k)/(m*LC*LC*LC)*( y[3]/m/r3 -3*y[6]*(y[4]*y[1]/m+y[5]*y[2]/m+y[6]*y[3]/m )/r5  );
	  
	  if (Cmass==0)
	  {		  	  
 		  dydx[1]= -k* y[4]/r3 - Crad * Radiationx + Cplane*FLorx;
		  dydx[2]= -k* y[5]/r3 - Crad * Radiationy + Cplane*FLory;
		  dydx[3]= -k* y[6]/r3 - Crad * Radiationz ; 

		  
		  dydx[4]= y[1]/m;
		  dydx[5]= y[2]/m;
		  dydx[6]= y[3]/m;
	  }
	  
	  if (Cmass==1)
	  {
		  dydx[1]= -k* y[4]/r3 - Crad * ( (2.0/3)*(k*k)/(m*LC*LC*LC)*( y[1]/m/r3 -3*y[4]*(y[4]*y[1]/m+y[5]*y[2]/m+y[6]*y[3]/m )/r5  ) ) + Cplane * ( Mamplitude*cos(y[6]*Mfrequency/LC+Mfrequency*x-phase1)       + Mamplitude/LC*y[3]/m*cos(y[6]*Mfrequency/LC+Mfrequency*x-phase1)                 );
		  dydx[2]= -k* y[5]/r3 - Crad * ( (2.0/3)*(k*k)/(m*LC*LC*LC)*( y[2]/m/r3 -3*y[5]*(y[4]*y[1]/m+y[5]*y[2]/m+y[6]*y[3]/m )/r5  ) ) + Cplane * ( Mamplitude*cos(y[6]*Mfrequency/LC+Mfrequency*x+PI/2-phase1)   + Mamplitude/LC*y[3]/m*cos(y[6]*Mfrequency/LC+Mfrequency*x+PI/2-phase1)                  );
		  dydx[3]= -k* y[6]/r3 - Crad * ( (2.0/3)*(k*k)/(m*LC*LC*LC)*( y[3]/m/r3 -3*y[6]*(y[4]*y[1]/m+y[5]*y[2]/m+y[6]*y[3]/m )/r5  ) ) + Cplane * (                                                               - Mamplitude/LC*(y[1]/m*cos(y[6]*Mfrequency/LC+Mfrequency*x-phase1)   -  y[2]/m*cos(y[6]*Mfrequency/LC+Mfrequency*x+PI/2-phase1)    ) ) ;

		  dydx[4]= y[1]/(   m*sqrt(1+  p2/(m*m*LC*LC)  )   );
		  dydx[5]= y[2]/(   m*sqrt(1+  p2/(m*m*LC*LC)  )   );
		  dydx[6]= y[3]/(   m*sqrt(1+  p2/(m*m*LC*LC)  )   );
	  }
}


  /*    -----------------Algorithms for solving the differential equations (from numerical recipes)--------------    */
  

  void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,double yscal[], double *hdid, double *hnext, void (*derivs)(double, double [], double []))
  {
	  int i,j;
	  double xsav,xest,h,errmax,temp;
	  double *ysav,*dysav,*yseq,*yerr;
	  static int nseq[IMAXX+1]={0,2,4, 6, 8, 12, 16, 24, 32, 48, 64, 96};
	  
	  ysav=vector(1,nv);
	  dysav=vector(1,nv);
	  yseq=vector(1,nv);
	  yerr=vector(1,nv);
	  x=vector(1,IMAXX);
	  d=matrix(1,nv,1,NUSE);
	  
	  h=htry;
	  xsav=(*xx);
	  
	  for(i=1; i<=nv;i++) 
	  {
		  ysav[i]=y[i];
		  dysav[i]=dydx[i];
	  }
	  
	  for(;;) 
	  {
		  for (i=1;i<=IMAXX;i++)
		  {
			  mmid(ysav, dysav, nv,xsav,h, nseq[i], yseq,derivs);
			  xest=(temp=h/nseq[i], temp*temp);
			  rzextr(i,xest, yseq, y,yerr,nv, NUSE);
			  
			  if(i>3)
			  {
				  errmax=0.0;
				  for(j=1;j<=nv;j++){
					  if( errmax < fabs(yerr[j]/yscal[j]) )
						  errmax=fabs(yerr[j]/yscal[j]);
				  }
				  
				  errmax /= eps;
				  
				  if(errmax <1.0) 
				  {
					  *xx +=h;
					  *hdid=h;
					  *hnext = i== NUSE? h*SHRINK: i==NUSE-1? h*GROW: ( h* nseq[NUSE -1]) /nseq[i];
					  
					  free_matrix(d,1,nv,1,NUSE);
					  free_vector(x,1,IMAXX);
					  free_vector(yerr,1,nv);
					  free_vector(yseq,1,nv);
					  free_vector(dysav,1,nv);
					  free_vector(ysav,1,nv);
					  
					  return;
				  }
			  }
		  }
		  
		  h *=0.25;
		  for( i=1;i<=( IMAXX - NUSE)/2;i++) h /=2.0;
		  if ((double)(*xx+h) == (*xx)) nrerror("Step size underflow in BSSTEP");
	  }
  }



  void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep, double yout[], void (*derivs)(double, double[], double[]))

/* Modified midpoint step. At xs, input the dependent variable vector y[1..nvar] and its deriva-
   tive vector dydx[1..nvar]. Also input is htot, the total step to be made, and nstep, the
   number of substeps to be used. The output is returned as yout[1..nvar], which need not
   be a distinct array from y; if it is distinct, however, then y and dydx are returned undamaged. */
  
  {
	  int n,i;
	  double x,swap,h2,h,*ym,*yn;
	  
	  ym=vector(1,nvar);
	  yn=vector(1,nvar);
	  h=htot/nstep;  // Stepsize this trip.
	  
	  for (i=1;i<=nvar;i++) 
	  {
		  ym[i]=y[i];
		  yn[i]=y[i]+h*dydx[i]; // First step.
	  }
	  
	  x=xs+h;
	  (*derivs)(x,yn,yout);        // Will use yout for temporary storage of derivatives. 
	  h2=2.0*h;
	  
	  for (n=2;n<=nstep;n++)
	  {       // General step.
		  for (i=1;i<=nvar;i++)
		  {
			  swap=ym[i]+h2*yout[i];
			  ym[i]=yn[i];
			  yn[i]=swap;
		  }
		  
		  x += h;
		  (*derivs)(x,yn,yout);
	  }
	  
	  for (i=1;i<=nvar;i++)          // Last step.
		  yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
	  
	  free_vector(yn,1,nvar);
	  free_vector(ym,1,nvar);
  }



  void rzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv, int nuse)
/* Use polynomial extrapolation to evaluate nv functions at x = 0by fitting a polynomial to a
   sequence of estimates with progressively smaller values x = xest, and corresponding function
   vectors yest[1..nv]. This call is number iest in the sequence of calls. Extrapolated function
   values are output as yz[1..nv], and their estimated error is output as dy[1..nv]. */

  {
	  int m1,k,j;
	  double yy,v,ddy,c,b1,b,*fx;
	  
	  fx=vector(1,nuse);
	  x[iest]=xest;             //Save current independent variable.
	  
	  if (iest ==1)
		  for (j=1;j<=nv;j++) 
		  {
			  yz[j]=yest[j];
			  d[j][1]=yest[j];
			  dy[j]=yest[j];
		  }
		  
		  else
		  {
			  m1=(iest < nuse ? iest:nuse);
			  for (k=1;k<=m1-1;k++)
				  fx[k+1]=x[iest-k]/xest;
			  for (j=1; j<=nv;j++) 
			  {
				  v=d[j][1];
				  d[j][1]=yy=c=yest[j];
				  
				  for (k=2; k<=m1; k++)
				  {
					  b1=fx[k]*v;
					  b=b1-c;
					  
					  if(b)
					  {
						  b=(c-v)/b;
						  ddy=c*b;
						  c=b1*b;
					  }
					  else
						  ddy=v;
					  
					  if (k!=m1) v=d[j][k];
					  d[j][k]=ddy;
					  yy += ddy;
				  }
				  
				  dy[j]=ddy;
				  yz[j]=yy;
			  }
		  }
		  free_vector(fx,1,nuse);
  }



/* User storage for intermediate results. Preset kmax and dxsav in the calling program. If kmax 6 =
   0 results are stored at approximate intervals dxsav in the arrays xp[1..kount], yp[1..nvar]
   [1..kount], wherekount is output by odeint. Defining declarations for these variables, with
   memory allocations xp[1..kmax] and yp[1..nvar][1..kmax] for the arrays, should be in
   the calling program. */


  void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,double hmin, int *nok, int *nbad, void (*derivs)(double, double [], double []),void (*rkqs)(double [], double [], int, double *, double, double, double [],double *, double *, void (*)(double, double [], double [])))

/* Runge-Kutta driver with adaptive stepsize control. Integrate starting values ystart[1..nvar]
   from x1 to x2 with accuracy eps, storing intermediate results in global variables. h1 should
   be set as a guessed first stepsize, hmin as the minimum allowed stepsize (can be zero). On
   output nok and nbad are the number of good and bad (but retried and fixed) steps taken, and
   ystart is replaced by values at the end of the integration interval. derivs is the user-supplied
   routine for calculating the right-hand side derivative, while rkqs is the name of the stepper
   routine to be used. */

  {
	FILE *ofp1,*ofp2,*ofp3,*ofp4,*ofp5,*ofp6,*ofp7,*ofp8,*ofp9,*ofp10,*ofp11,*ofp12,*ofp13,*ofp14,*ofp15,*ofp16;
      
	  //ofp1=fopen("..//radius.dat","w");
	  ofp2=fopen("..//orbit.dat","w");
	  ofp3=fopen("..//smajor.dat","w");
	  //ofp4=fopen("..//sminor.dat","w");
	  ofp5=fopen("..//eccen.dat","w");
	  // ofp6=fopen("..//pvalue.dat","w");	  
	  ofp7=fopen("..//theta0.dat","w");
	  //ofp8=fopen("..//energy.dat","w");
	  //ofp9=fopen("..//Lz.dat","w");
	  //ofp10=fopen("..//period.dat","w");
	  //ofp11=fopen("..//Net-angfreq.dat","w");
      //ofp12=fopen("..//Instantaneous-angfreq.dat","w");
      //ofp13=fopen("..//Approximate-Lz.dat","w");
      //ofp14=fopen("..//Approximate-angfreq.dat","w");
	  //ofp15=fopen("..//Approximate-energy.dat","w");
	  //ofp16=fopen("..//Fv.dat","w");

	  //int nstp;
	  int i,k,j=0,count=0, count1=1;
	  double temp=0.0;
	  double xsav,x,hnext,hdid,h;
	  double *yscal,*y,*dydx;
	  double tempx[30],tempy[30],semimajor,semiminor,eccentricity,pvalue,theta0,costheta,sintheta;
	  double rreverse,energy,Lz,period,angfreq,instant_angfreq,approx_Lz,approx_angfreq;
	  double aa=0.0,bb=0.0,cc=0.0,dd=0.0,ee=0.0,ff=0.0,gg=0.0,hh=0.0,ii=0.0,jj=0.0,kk=0.0,ll=0.0;
	  double parameter1,parameter2,parameter3;
	  double radius1,radius_1;
  

	  yscal=vector(1,nvar);
	  y=vector(1,nvar);
	  dydx=vector(1,nvar);
	  
	  x=x1;
	  
	  h=(x2>x1)? fabs(h1):-fabs(h1);
	  *nok = (*nbad) = kount = 0;
	  
	  for (i=1;i<=nvar;i++) y[i]=ystart[i];
	  
	  if (kmax > 0) xsav=x-dxsav*2.0;    //Assures storage of first step.
	  //  for (nstp=1;nstp<=MAXSTP;nstp++) {    //Take at mos t MAXSTP steps.
	  
	  for (;;)
	  {    //Modified as "Take at unlimited steps.
		  (*derivs)(x,y,dydx);
		  
		  for (i=1;i<=nvar;i++)
			  /* Scaling used to monitor accuracy. This general-purpose choice can be modified if need be. */
			  yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		  if (kmax>0)
		  {
			  if (fabs(x-xsav)>fabs(dxsav))
			  {
				  if (kount<kmax-1)
				  {
					  xp[++kount]=x; 
					  for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
					  xsav=x;
				  }
			  }
		  }
		  
		  if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;     //If stepsize can overshoot, decrease.
		  (*bsstep)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		  if (hdid == h) 
		  {
			  ++(*nok);
			  
			  tempx[j]=y[4];
			  tempy[j]=y[5];

			  if (j==28)
			  {
				  for (k=0;k<=j;k++)
				  {
					  costheta=tempx[k]/sqrt(tempx[k]*tempx[k]+tempy[k]*tempy[k]);
					  sintheta=tempy[k]/sqrt(tempx[k]*tempx[k]+tempy[k]*tempy[k]);
					  rreverse=1.0/sqrt(tempx[k]*tempx[k]+tempy[k]*tempy[k]);
					  aa=aa+1.0;
					  bb=bb+costheta;
					  cc=cc+sintheta;
					  dd=dd+costheta;
					  ee=ee+costheta*costheta;
					  ff=ff+costheta*sintheta;
					  gg=gg+sintheta;
					  hh=hh+sintheta*costheta;
					  ii=ii+sintheta*sintheta;
					  jj=jj+rreverse;
					  kk=kk+rreverse*costheta;
					  ll=ll+rreverse*sintheta;
				  }

				  parameter3=((aa*ll-jj*gg)*(aa*ee-bb*dd)-(aa*kk-jj*dd)*(aa*hh-bb*gg))/((aa*ii-cc*gg)*(aa*ee-bb*dd)-(aa*ff-cc*dd)*(aa*hh-bb*gg));
				  parameter2=((aa*kk-jj*dd)-(aa*ff-cc*dd)*parameter3)/(aa*ee-bb*dd);
				  parameter1=(jj-parameter2*bb-parameter3*cc)/aa;

				  pvalue=sqrt(1.0/(parameter2*parameter2+parameter3*parameter3));
				  eccentricity=1.0/(parameter1*pvalue);
                                      // old way (get rid of if new way works) theta0=atan(parameter3/parameter2);
				  semimajor=pvalue*eccentricity/(1.0-eccentricity*eccentricity);
				  semiminor=pvalue*eccentricity/sqrt(1.0-eccentricity*eccentricity);
                                  
				  costheta0=-parameter2*pvalue;
                                  sintheta0=-parameter3*pvalue;
                                  acost=acos(costheta0);
                                  asint=asin(sintheta0);
				  if (sintheta0>=0 && costheta0>=0) Fnew=1;  //1st quatrant
                                  if (sintheta0>=0 && costheta0< 0) Fnew=2;  //2nd quatrant
                                  if (sintheta0 <0 && costheta0< 0) Fnew=3;  //3rd quatrant
                                  if (sintheta0 <0 && costheta0>=0) Fnew=4;  //4th quatrant
				  if (Fnew==1 && Fold==4) 
				      
				      Nnew=Nold+1;
				  
				  else if (Fnew==4 && Fold==1) 

				    {
  				      Nnew=Nold-1;
				    }

				  if (Fnew==1||Fnew==2) 
				     { 
				       theta0temp = acos(costheta0);
				     }
				  else
				     {
				       theta0temp = 2*PI- acos(costheta0);
				     }
				  
				  theta0 = Nnew*2*PI + theta0temp;
				  Nold = Nnew;
				  Fold = Fnew;

				  // The approximate value of energy should be:  energy=-(4.80298e-10)*(4.80298e-10)/(2*semimajor);
                  radius1=pow(y[4]*y[4]+y[5]*y[5]+y[6]*y[6],0.5);
				  if (radius1>1e-4) nrerror("\nRadius has exceeded Rmax=1e-4cm, so program has terminated.\nIf dissatisfied, please submit & file disatisfaction form #999,999.");
				  Lz=y[4]*y[2]-y[5]*y[1];
				  period=2*PI*( pow(m,0.5) )*( pow(semimajor,1.5) )/4.80298e-10;
				  angfreq=2*PI/period;
				  instant_angfreq=Lz/(m*pow(radius1,2));
				  approx_Lz=e*(pow(m*semimajor*(1-eccentricity*eccentricity),0.5));
				  approx_angfreq=e*(pow((1-eccentricity*cos(theta-theta0)),2))/( (pow(m,0.5))*(pow(semimajor*(1-eccentricity*eccentricity),1.5) ) );
				  
				  aa=0.0;	bb=0.0;		cc=0.0;		dd=0.0;
				  ee=0.0;	ff=0.0;		gg=0.0;		hh=0.0;
				  ii=0.0;	jj=0.0;		kk=0.0;		ll=0.0;				  
				  if (count1>0)  
				  {
					  fprintf(ofp3,"%e, %e\n",x,semimajor);
					  //fprintf(ofp4,"%e, %e\n",x,semiminor);
					  fprintf(ofp5,"%e, %e\n",x,eccentricity);
					  //  fprintf(ofp6,"%e, %e\n",x,pvalue);
					  //  fprintf(ofp7,"%e, %e, %e, %e, %e, %e, %e\n",x, theta0,costheta0, acost,sintheta0,asint);
					  fprintf(ofp7,"%e, %e\n",x,theta0);
					  //fprintf(ofp9,"%e, %e\n",x,Lz);
					  //fprintf(ofp10,"%e, %e\n",x,period);
					  //fprintf(ofp11,"%e, %e\n",x,angfreq);
					  //fprintf(ofp12,"%e, %e\n",x,instant_angfreq);
					  //fprintf(ofp13,"%e, %e\n",x,approx_Lz);
                      //fprintf(ofp14,"%e, %e\n",x,approx_angfreq);
					  //fprintf(ofp15,"%e, %e\n",x,-(e*e)/(2*semimajor));
                      
					  radius_1=pow(y[4]*y[4]+y[5]*y[5]+y[6]*y[6],0.5);
					  energy=-(e*e)/radius_1 + ( y[1]*y[1] + y[2]*y[2]  + y[3]*y[3] )/(2*m);
                      //fprintf(ofp8,"%.16e %.16e\n",x,energy);

					  count1=0;
				  }
				  count1++;
				  j=-1;
			  }

			  j++;

                   //I am experimenting below with CounterPrint for printing//
			  if (count==CounterPrint)
			  {
				  //if ((x>6.7e-11) && (x<6.87e-11))
				  fprintf(ofp2,"%.16e, %.16e, %.16e\n", x, y[4],y[5]);
				  // if ((x>6.7e-11) && (x<6.87e-11))
				      ////  radius_1=pow(y[4]*y[4]+y[5]*y[5]+y[6]*y[6],0.5);
				  //  if ((x>6.7e-11) && (x<6.87e-11)) 
					  ////fprintf(ofp1,"%.16e %.16e\n", x, radius_1 );
                  // if ((x>6.7e-11) && (x<6.87e-11)) 
					  //energy=-(e*e)/radius_1 + ( y[1]*y[1] + y[2]*y[2] + y[3]*y[3] )/(2*m);
					  //// energy=-(e*e)/radius_1 + ( y[1]*y[1] + y[2]*y[2]  + y[3]*y[3] )/(2*m);
				  // if ((x>6.7e-11) && (x<6.87e-11)) 
					  //// fprintf(ofp8,"%.16e %.16e\n",x,energy);
                  // if ((x>6.7e-11) && (x<6.87e-11)) 
				  //    fprintf(ofp16,"%.16e %e.16\n",x,  (y[1]/m)*(- Crad * Radiationx + Cplane*FLorx) +
   				  //                                (y[2]/m)*(- Crad * Radiationy + Cplane*FLory)  );
				  count=0;
			  }
			  count++;  
		  }
		  
		  else ++(*nbad);
		  if ((x-x2)*(x2-x1) >= 0.0)
		  {           //Are we done?   
			  for (i=1;i<=nvar;i++) ystart[i]=y[i];
			  if (kmax)
			  {
				  xp[++kount]=x;                //Save final step.
				  for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			  }
			  
			  //fclose(ofp1);
			  fclose(ofp2);
			  fclose(ofp3);
			  //fclose(ofp4);
			  fclose(ofp5);
			  // fclose(ofp6);
			  fclose(ofp7);
			  //fclose(ofp8);
			  //fclose(ofp9);
			  //fclose(ofp10);
			  //fclose(ofp11);
			  //fclose(ofp12);
			  //fclose(ofp13);
			  //fclose(ofp14);
			  //fclose(ofp15);
              //fclose(ofp16);

			  free_vector(dydx,1,nvar);
			  free_vector(y,1,nvar);
			  free_vector(yscal,1,nvar);
			  return;                       //Normal exit.
		  }
		  
		  if (fabs(hnext) <= hmin) 
			  nrerror("Step size too small in odeint");
		  h=hnext;
	  }
	  nrerror("Too many steps in routine odeint");
  }
  
  
  double search (double z)
  {
	  double z0=0.0;
	  double z1=2*PI;
	  double midpoint;
	  
	  midpoint=(z0+z1)/2.0;
	  
	  while (  fabs(myfun(midpoint)-z) > EPS  )
	  {
		  if (myfun(midpoint)<z)
			  z0=midpoint;
		  else
			  z1=midpoint;
		      midpoint=(z0+z1)/2.0;
	  }
	  return midpoint;
  }


  double myfun (double thetavalue)
  {
	  double myresult;
	  
	  myresult=atan( sqrt(1+epsilon)*tan(thetavalue/2.0)/sqrt(1-epsilon) )/PI+1.0/2.0*sqrt(1-epsilon*epsilon)*epsilon*sin(thetavalue)/PI/(1-epsilon*cos(thetavalue));
	  
	  if (myresult<0)
		  myresult=myresult+1;
	  
	  return myresult;
  }
  
