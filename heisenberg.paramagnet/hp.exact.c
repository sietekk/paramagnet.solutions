/*
Project: Heisenberg Paramagnet
            Exact Solution
Author: Michael Conroy
Class: PHY 471
Professor: Dr. Enjalran
Date: February 2014
*/

#include <stdio.h> //fprintf, etc.
#include <stdlib.h> //I forget why...
#include <math.h> //Exponents
#include <tgmath.h> //Hyperbolic trig

#define MAX_COUNT 1000 //Temp divisions; # individual temps
#define MAX_TEMP 100.0 //Maximum temperature
#define MIN_TEMP 0.1 //Minimum temperature


#define PI 3.141592653589793239 //Check length
#define N 1 //Number of particles
#define MU 1.0 //Magnetic moment
#define B 1.0 //Magnetic field
#define K_b 1.0 //Boltzmann constant

double partitionFunc(double beta)
{
    //Partition fxn value
    double Z;
    double betaMuB = beta*MU*B;
    double base = ((4*PI)/(betaMuB))*sinh(betaMuB);
    Z = pow(base,N);
    return Z;
}

double avgMagnetization(double beta)
{
    //Avg magnetization value
    double M;
    double betaMuB = beta*MU*B;
    M = N*MU*(atanh(betaMuB)-1.0/(betaMuB));
    return M;
}

double heatCapacity(double beta)
{
    //Heat capacity value
    double C;
    double betaMuB = beta*MU*B;
    C = 2.0*N*K_b*
            (betaMuB)*(betaMuB)*
            (atanh(betaMuB)*atanh(betaMuB));
    return C;
}

double avgEnergy(double beta)
{
    //Avg energy
    double E;
    double betaMuB = beta*MU*B;
    E = -(N*MU*B*(atanh(betaMuB)-1.0/(betaMuB)));
    return E;
}

void print_to_file(double data_array[][8])
{
    int i; //Counter
    char file_name[17] = "hp_exact_data.dat";
    FILE *file;
    file = fopen(file_name,"w+");
    fprintf(file,
            "T\t<E>\t<E>/N\t<M>\t<M>/N\tC/N\tX/N\n");
    for (i = 0; i <= MAX_COUNT; i++)
    {
        fprintf(file,
                "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                data_array[i][1], data_array[i][2], data_array[i][3],
                data_array[i][4], data_array[i][5], data_array[i][6],
                data_array[i][7], data_array[i][8]);
    }
    fclose(file);
}

int main()
{
    //Declare variables
    int i, j; //Loop counter
    double beta; //Inverse temperature
    double dT; //Temperature step/interval magnitude
    double Z, avgM, C, avgE; //Cosh term, exponent term, partition fxn
    double data_array[MAX_COUNT][8]; //Holds data for file creation
        /*[1] = temp, [2] = avg energy, [3] = avg energy per site,
            [4] = avg magentization, [5] = avg mag per site,
            [6] = heat capacity, [7] = avg heat capacity (C/N),
            and [8] = partition function */

    //Temperatures
    dT = (double)((MAX_TEMP-MIN_TEMP)/MAX_COUNT);
    printf("dT: %f\n",dT);
    for (i = 1; i <= MAX_COUNT; i++)
    {
        data_array[i][1] = MAX_TEMP-(i*dT);
        printf("Temp %d: %f\n",i,data_array[i][1]);
    }

    //Calculate and store data
    for (j = 0; j <= MAX_COUNT; j++)
    {
        //Calculations
        beta = 1.0/data_array[j][1];
        avgE = avgEnergy(beta);
        avgM = avgMagnetization(beta);
        C = heatCapacity(beta);
        Z = partitionFunc(beta); //Not sure why I added this lolz.

        //Store data
        data_array[j][2] = avgE;
        data_array[j][3] = avgE/N;
        data_array[j][4] = avgM;
        data_array[j][5] = avgM/N;
        data_array[j][6] = C;
        data_array[j][7] = C/N;
        data_array[j][8] = Z;
    }

    //Print data
    print_to_file(data_array);

    //End
    printf("End of program!");

    return 0;
}

