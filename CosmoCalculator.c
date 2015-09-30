/*
 *  CosmoCalculator.c
 *  This program solves numerically Friedmann Equations and provides a file with the data related to the st. of universe
 *
 *  Created by Guadalupe Cañas Herrera and Palmerina González Izquierdo on 28/11/14.
 *  Copyright 2014 __gch24__pgi25. All rights reserved.
 */


/*Include the library required*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#include "TH1.h"

/* Define constants for the whole program*/

#define C 299792.458 //light velocity in km/s
#define MPARSEC 3.08568E19 //transforms km into Mpc
#define GYEAR 3.1536E16 // transforms s into Gyr 31 556 926
#define N 5

/*Variables used during the whole program*/

double omegaM;  //Matter Density
double omegaDarkE; //Dark Energy Density
double omegaR; //Radiation Density
double omegaK; //Curvature Density
double hubbleConstant; //Hubble Constant
double z; //Redshift
double scaleFactor; //Scale factor of the universe
double deceleration_parameter; //Define the deceleration parameter


//double Pi;
//double lroots[N];
//double weight[N];
//double lcoef[N + 1][N + 1] = {{0}};


/* FUNCTIONS PROTOTYPES (declaring functions that we will use in the code) */
double ecuacion1(double x);
double ecuacion2(double x);
double integral(double m, double n, double k, double (*)(double x));
//double simpson(double a, double b, double N, double (*ecuacion)(double x));
//double lege_inte(double (*ecuacion)(double x), double a, double b);
//void lege_coef();
//void lege_roots();


/*Most important function in the program*/
int main()
{
    //Pi = atan2(1, 1) * 4;
    
    //lege_coef();
    //lege_roots();
    
    int k; //int used in the for
	char filename[256]; //required to write in the document
	FILE* output; //create the document for data

    //We ask in the screen the following values
    printf("WELCOME TO THE COSMOLOGICAL CALCULATOR\n\n");
    printf("Please, introduce the data related to the composition of the Universe\n\n");
    printf("Enter the MATTER DENSITY (with decimals): ");
    scanf("%lf",&omegaM);
    printf("Enter the DARK ENERGY DENSITY (with decimals): ");
    scanf("%lf",&omegaDarkE);
    printf("Enter the RADIATION DENSITY (with decimals): ");
    scanf("%lf",&omegaR);
    printf("Enter the HUBBLE CONSTANT (in km/(Mpc.s)): ");
    scanf("%lf",&hubbleConstant);
    printf("Enter the REDSHIFT (with decimals): ");
    scanf("%lf",&z);
   
    
    /*FIRST PART OF THE MAIN FUNCTION*/
    
    //Calculate the value of curvature density
    omegaK=1.0-(omegaM+omegaR+omegaDarkE);
    
    
    //Calculate scaleFactor
    scaleFactor=1.0/(1.0+z); //By definition, R(t_0) is set up to 1
    
    //Calculate the deceleration Parameter
    deceleration_parameter=0.5*omegaM+omegaR-omegaDarkE;
    
    //Calculate the integral of ecuacion1
	double Integral1 = integral(scaleFactor, 1.0, 1000.0, ecuacion1);

    /*
     * Calculus of distances obtained through the integral Integral1 of equation equation1
     *
     */
    
    //Radial Comoving Distance
    double dc=C*(1.0/hubbleConstant)*Integral1;
    
    double dp;//Transverse Comoving Distance
    
    if(omegaK<0){//k=1
        dp=C*(1.0/(hubbleConstant*sqrt(-omegaK)))*sin(sqrt(-omegaK)*Integral1);
    }else if(omegaK>0){//k=-1
         dp=C*(1.0/(hubbleConstant*sqrt(omegaK)))*sinh(sqrt(omegaK)*Integral1);
    }else{//k=0
         dp=C*(1.0/hubbleConstant)*Integral1;
    }
    
   
    double da=dp*scaleFactor;//Angular Diameter distance
    double dl=dp*(1+z);//Luminosity Distance
    
    //Calculate integral of ecuacion2
	double Integral2 = integral(0.0001, 1.0, 1000.0, ecuacion2); //0 is not allowed as initial point
    double Integral2_z = integral(scaleFactor, 1.0, 1000.0, ecuacion2); //0 is not allowed as initial point
    double Integral2_fromz=integral(0.0001, scaleFactor, 1000.0, ecuacion2);
    
    
    /*
     * Calculus of the age of the universe obtained through the Integral2 of equation2
     *
     */
    
    
    double edad=(MPARSEC/(hubbleConstant*GYEAR))*Integral2;
    double edad_z=(MPARSEC/(hubbleConstant*GYEAR))*Integral2_z;
    double edad_fromz =(MPARSEC/(hubbleConstant*GYEAR))*Integral2_fromz;
    
    /*
     * Show the results in the screen
     */
    
    printf("\n \nRESULTS \n");
    printf("\n \nDeceleration Parameter: %lf \n", deceleration_parameter);
    printf("Light-Travel Distance: %lf Mpc\n", dp);
    printf("Comoving Radial Distance: %lf Mpc\n", dc);
    printf("Angular Diameter Distance: %lf Mpc\n", da);
    printf("Luminosity Distance: %lf Mpc\n", dl);
    printf("Light Travel Time from z=%f until now:  %lf Gyr\n", z, edad_z);
    printf("Age of the Universe at z=%f: %lf Gyr\n", z, edad_fromz);
    printf("Age of the Universe at present time: %lf Gyr\n", edad);
    


	//DOCUMENT FOR GRAPHS DATA
    
	if(omegaK<0){ //define the curvature that allows us to distinguish among documents
		k=1;
	}else if(omegaK>0){
		k=-1;
	}else if (omegaK==0){
		k=0;
	}
    
    printf("The Curvature of the Universe is k=%d \n", k);

    
    
    /*SECOND PART OF THE MAIN FUNCTION*/
    
	sprintf (filename, "datos_cosmo_%d.dat", k); //write the name of the document
	output = fopen (filename, "w"); //open the document

    fprintf (output, "%s\t", "d_propia"); //write the names of the columns
	fprintf (output, "%s\t", "d_angular");
	fprintf (output, "%s\t", "d_luminosidad");
	fprintf (output, "%s\t", "factor_escala");
	fprintf (output, "%s\t", "edad    ");
	fprintf (output, "%s\t", "redshift");
	fprintf (output,"\n");
    
	for (z = 0; z <= 100; z += 0.01) //for loop in order to calculate
                                   //the distances and the age of the universe as a function of z
	{
		//Calculate scaleFactor
        scaleFactor=1.0/(1.0+z);
        //Calculate a new scaleFactor in order to view time futher that R(t)=1
		double scaleFactorNew=2.5/(1.0+z);

		double Integral1 = integral(scaleFactor, 1.0, 1000.0, ecuacion1);

		double dp;//distancia propia
		if(omegaK<0){//k=1
			dp=C*(1.0/(hubbleConstant*sqrt(-omegaK)))*sin(sqrt(-omegaK)*Integral1);
		}else if(omegaK>0){//k=-1
			dp=C*(1.0/(hubbleConstant*sqrt(omegaK)))*sinh(sqrt(omegaK)*Integral1);
		}else{//k=0
			dp=C*(1.0/hubbleConstant)*Integral1;
		}
    
		double da=dp*scaleFactor;//distancia angular
		double dl=dp*(1+z);//distancia luminosidad
    
		//Age of universe
		double Integral2 = integral(0.0001, scaleFactorNew, 1000.0, ecuacion2);

		double edad=(MPARSEC/(hubbleConstant*GYEAR))*Integral2;
        
        //write all the data in the document
		fprintf (output, "%lf\t", dp);
		fprintf (output, "%lf\t", da);
		fprintf (output, "%lf\t", dl);
		fprintf (output, "%lf\t", scaleFactorNew);
		fprintf (output, "%lf\t", edad);
		fprintf (output, "%lf\t", z);
		fprintf (output,"\n");
	}

	fclose(output); //close the document
    
    printf("\n\nRemember: the program has just created a .dat file with data \n\n");
    
    return 0;
    
}


/*
 *
 * Functions required in the main function
 *
 */

/*EQ1: ecuacion1 defines the function required to calculate distances*/
double ecuacion1(double x)
{
    double inverse_x=(1.0/x);
    return 1/(x*x*sqrt(omegaDarkE+omegaK*inverse_x*inverse_x+omegaM*inverse_x*inverse_x*inverse_x+omegaR*inverse_x*inverse_x*inverse_x*inverse_x));
}

/*EQ2: ecuacion2 defines the function required to calculate the age of the universe*/
double ecuacion2(double x)
{
    double inverse_x=(1.0/x);
    return 1/(x*sqrt(omegaDarkE+omegaK*inverse_x*inverse_x+omegaM*inverse_x*inverse_x*inverse_x+ pow(omegaR*inverse_x,4)));
}



/*INTEGRAL: calculates required integrals using trapezoidal method*/
double integral(double m, double n, double k, double (*ecuacion)(double x))
{
	int i;
    double h=fabs((m-n)/k);
    double suma=(ecuacion(m)+ecuacion(n))/2;
    
    for (i = 1; i<k; i++)
    {
        suma += ecuacion(m+i*h);
    }
    
    return suma*h;    
}


/*
double simpson(double a, double b, double N, double (*ecuacion)(double x))
{
   
    double h=(b-a)/N;
    double m=N/2;
    double S=0;
    double P=0;
    int i=0;
    
    for (i = 1; i<m; i++)
    {
        double val=ecuacion(a+(2*i-1)*h);
        S=S+val;
    }
    for (i=1; i<m-1; i++)
    {
        double val2=ecuacion(a+(2*i)*h);
        P=P+val2;
    }
    return 2*((h/3)*(ecuacion(a)+ecuacion(b))+((4/3)*h*S)+(2/3)*h*P);
}
*/




//metodo de cuadratura de gauss

/*
void lege_coef()
{
    int n, i;
    lcoef[0][0] = lcoef[1][1] = 1;
    for (n = 2; n <= N; n++) {
        lcoef[n][0] = -(n - 1) * lcoef[n - 2][0] / n;
        for (i = 1; i <= n; i++)
            lcoef[n][i] = ((2 * n - 1) * lcoef[n - 1][i - 1]
                           - (n - 1) * lcoef[n - 2][i] ) / n;
    }
}

double lege_eval(int n, double x)
{
    int i;
    double s = lcoef[n][n];
    for (i = n; i; i--)
        s = s * x + lcoef[n][i - 1];
    return s;
}

double lege_diff(int n, double x)
{
    return n * (x * lege_eval(n, x) - lege_eval(n - 1, x)) / (x * x - 1);
}

void lege_roots()
{
    int i;
    double x, x1;
    for (i = 1; i <= N; i++) {
        x = cos(Pi * (i - .25) / (N + .5));
        do {
            x1 = x;
            x -= lege_eval(N, x) / lege_diff(N, x);
        } while (x != x1);
 
        lroots[i - 1] = x;
        
        x1 = lege_diff(N, x);
        weight[i - 1] = 2 / ((1 - x * x) * x1 * x1);
    }

double lege_inte(double (*ecuacion)(double x), double a, double b)
{
    double c1 = (b - a) / 2, c2 = (b + a) / 2, sum = 0;
    int i;
    for (i = 0; i < N; i++)
        sum += weight[i] * ecuacion(c1 * lroots[i] + c2);
    return c1 * sum;
}
*/

