#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>

//---------------Begin Definitions---------------//
	#define turn_angle  2.*M_PI
//Main parameters
	#define N 1000	 //number of Kuramoto oscillators
	#define dt .001 //time step
	#define T 1000 //simulation runtime

//For fixed value of K-s in simulation
	#define K1 1.
	#define K2 2.
//For sweeping of K
	#define K0 0.
	#define dK .02
	#define K_max 1.5 //
//---------------End Definitions---------------//

//Run bools:
	//check values for initial configurations:
	bool check_initial = false;
	//Gaussian distributed frequencies [N(0,1)]
	bool gaussian_frequencies = false;
	//ODE + sweeping K
	bool sweeping_K = false;
	//ODE + constant K-s


void PrintParams();
float * RandUnifPhase();
float * RandUnifFreq();

float * RandGauss();
float PeriodicPosition(float angular_pos);
void EulerStep(float *phase, float *frequencies, int K);



//Functions declarartion:

int main(void)    
{	
	clock_t begin = clock();
	srand(time(NULL));
	sleep(0);

	PrintParams();

	//Declarations
	int i,j,k;

	float *phases;
	float *ang_freqs;
	float sum_phases = 0;
	float sum_ang_freqs = 0;

	phases = RandUnifPhase();
	if(gaussian_frequencies==true){
		ang_freqs= RandGauss();
	}
	else{
		ang_freqs= RandUnifFreq();
	}

//Test distribution for initial configurations:
	if(check_initial==true){
		for ( i = 0; i < N; ++i) {
		    printf( "\nphases[%d] = %f\n", i, phases[i]);
		    sum_phases+=phases[i];
		    printf( "ang_freqs[%d] = %f\n", i, ang_freqs[i]);
		    sum_ang_freqs+=ang_freqs[i];
		}
		float mean_phases, mean_ang_freq, var_phases, var_ang_freq;
		mean_phases = sum_phases/N;
		mean_ang_freq = sum_ang_freqs/N;

		for ( i = 0; i < N; ++i) {
		    var_phases+=(phases[i]-mean_phases)*(phases[i]-mean_phases);
		    var_ang_freq+=(ang_freqs[i]-mean_ang_freq)*(ang_freqs[i]-mean_ang_freq);
		}
		
		printf("\n\n\nMean Phases = %.5f\n",mean_phases/N);
		printf("Variance Phases = %.5f\n\n",var_phases/N);
		printf("Mean ang_freqs = %.5f\n",mean_ang_freq/N);
		printf("Variance ang_freqs = %.5f\n",var_ang_freq/N);
}

	printf("InitialPhase=%.5f,\tInitialFrequency%.5f\n",phases[N-1],ang_freqs[N-1]);
	int T_split = (int)(T/20);
	for(i=0;i<T;i++){
		if(i%T_split==0)
		{
			printf("\nProcess at %d/100 ", 100*(int)(i)/T);
		}
		EulerStep(phases, ang_freqs, K1);
	}
	printf("\n\n");	

	printf("FinalPhase=%.5f,\tFinalFrequency%.5f\n",phases[N-1],ang_freqs[N-1]);



	//List of input parameters
	PrintParams();
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Execution time: %.5f seconds\n\n",time_spent);
  	return 0;
}


//Function implementations

void PrintParams()
{
	printf("\n\n\n//------------------------Input parameters------------------------//\n");
	printf("//\t\tNumber of Oscillators = %d\n",N);
	printf("//\t\tRuntime simulation = %d\n",T);
	printf("//\t\tdt = %.5f\n",dt);
	printf("//\t\tGaussian initial frequencies? = %s\n", gaussian_frequencies ? "True!" : "False!");
	printf("//\t\tCheck initial conditions? = %s\n", check_initial ? "True!" : "False!");

	printf("//----------------------------------------------------------------//\n\n");
}

float * RandUnifPhase( ) {
	/* Generate array uniformly distributed of float random variables*/
    static float r[N];
    int i;
    float max_freq = 5; //if max2pi=false, the resulting random uniform distrib will be in [-max_freq,max_freq]
    srand( (unsigned)time(NULL) );
	    for ( i = 0; i < N; ++i){
	      		r[i] = ((float)rand()/(float)(RAND_MAX)) * turn_angle;
	    	}
    return r;
}
float * RandUnifFreq( ) {
	/* Generate array uniformly distributed of float random variables*/
    static float r[N];
    int i;
    float max_freq = 5; //if max2pi=false, the resulting random uniform distrib will be in [-max_freq,max_freq]
    srand( (unsigned)time(NULL) );
    for ( i = 0; i < N; ++i) {
		r[i] = ((float)rand()/(float)(RAND_MAX)) * max_freq -max_freq*.5;
    }

    

    return r;
}
float * RandGauss(  ) {
	/*https://www.doc.ic.ac.uk/~wl/papers/07/csur07dt.pdf for gaussian array generation*/
	/*Box-Muller transform*/
	srand( (unsigned)time(NULL) );

    static float r[N];
    float U1,U2,sqrt_term,phase;
    int i;
    for ( i = 0; i < N/2; ++i) {
      	U1 = ((float)rand()/(float)(RAND_MAX));
      	U2 = ((float)rand()/(float)(RAND_MAX));
      	
      	sqrt_term = sqrt(-2*log(U1));
      	phase = U2 * turn_angle;
      	r[2*i] = sqrt_term*cos(phase);
      	r[2*i+1] = sqrt_term*sin(phase);
    }
    return r;
}
float PeriodicPosition(float angular_pos){
	if (angular_pos>turn_angle){
		angular_pos-=turn_angle;
	}
	if (angular_pos<0){
		angular_pos+=turn_angle;
	}
	return angular_pos;
}

void EulerStep(float *phases, float *ang_freqs, int K)
{
	int i,j;

	float phase_updated[N] = {0} ; 
	float frequencies_updated[N] = {0} ; 
	float sum_term = 0.;
	float floatN;
	floatN = (float)N;
  	for(i = 0;i < N; i++)
  	{	
  		//Copy initial phase
  		phase_updated[i] = phases[i];
  		//Perform evaluation of additional term
  		for(j = 0;j < N; j++)
  			{
  				sum_term += sin(phases[j])-sin(phases[i]);
  			}
  			sum_term = sum_term*(K/floatN);
  		frequencies_updated[i]=ang_freqs[i] + sum_term;
  		phase_updated[i] +=  frequencies_updated[i]*dt;
  	}

  	for(i = 0; i < N;i++)
  	{
  		phases[i] = PeriodicPosition(phase_updated[i]);
  		ang_freqs[i] = frequencies_updated[i];
  	}	

}
