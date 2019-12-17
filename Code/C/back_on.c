#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <libgen.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>



//---------------Begin Definitions---------------//
	#define turn_angle  2.*M_PI
//Main parameters
	#define N 100	 //number of Kuramoto oscillators
	
	#define dt .001 //time step
	#define T 1 //simulation runtime

//For fixed value of K-s in simulation
	#define K1 1.
	#define K2 2.
//For sweeping of K
	#define K0 0.
	#define dK .02
	#define K_max 1.5 //
	#define PATH_MAX 1000
//---------------End Definitions---------------//

//Run bools:
	//check values for initial configurations:
	bool check_initial = false;
	//Gaussian distributed frequencies [N(0,1)]
	bool gaussian_frequencies = true;
	//ODE + sweeping K
	bool sweeping_K = false;
	//ODE + constant K-s

const char * CreateResultsFolder();
void PrintParams();
float * RandUnifPhase();
float * RandUnifFreq();

float * RandGauss();
float PeriodicPosition(float angular_pos);
void EulerStep(float *phase, float *frequencies, float K);
float complex OrderParam(float *phases);
float complex FreqOrderParam(float *ang_freqs);


float iN=(float)(1/N);

//Functions declarartion:

int main(void)    
{	
	clock_t begin = clock();
	srand(time(NULL));
	sleep(0);
;

	//filename=strcat(filename,".csv");




	const char* results_path = CreateResultsFolder();
	printf("%s", results_path);

	PrintParams();

	//Declarations
	int i,j,k;
	float complex iN = (float)(1/N)+I*0; //inverse of N
	float *phases;
	float *ang_freqs;
	double complex ord_param = 0 + 0 * I;
	double freq_ord_param = 0 + 0 * I;
	

	phases = RandUnifPhase();
	if(gaussian_frequencies==true){
		ang_freqs= RandGauss();
	}
	else{
		ang_freqs= RandUnifFreq();
	}

//Test distribution for initial configurations:
	if(check_initial==true){
		float sum_phases = 0;
		float sum_ang_freqs = 0;
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
	
	//----------------------START SINGLE RUN LOOP----------------------//
	//printf("InitialPhase=%.5f,\tInitialFrequency%.5f\n",phases[N-1],ang_freqs[N-1]);
	ord_param = OrderParam(phases);
	freq_ord_param = FreqOrderParam(ang_freqs);
	printf("Initial OrderParameter = %.3f + %.3fi\n", creal(ord_param),cimag(ord_param));
	printf("Initial FreqOrderParameter = %.3f + %.3fi", creal(freq_ord_param),cimag(freq_ord_param));

	int T_split = (int)(T/100);
	for(i=0;i<T;i++){
		if(i%T_split==0)
		{
			printf("\nProcess at %d/100\n", 100*(int)(i)/T);
			ord_param = OrderParam(phases);
			printf("Order parameter = %.3f + %.3fi\n", creal(ord_param),cimag(ord_param));
			freq_ord_param = OrderParam(ang_freqs);
			printf("Freq Order parameter = %.3f + %.3fi\n", creal(freq_ord_param),cimag(freq_ord_param));


		}
		EulerStep(phases, ang_freqs, 100);
	}

	//----------------------END SINGLE RUN LOOP----------------------//

	printf("\n\n");	
	
	ord_param = OrderParam(phases);
	freq_ord_param = FreqOrderParam(ang_freqs);

	printf("Final Order parameter = %.3f + %.3f\n", creal(ord_param),cimag(ord_param));
	printf("Final FreqOrderParameter = %.3f + %.3fi", creal(freq_ord_param),cimag(freq_ord_param));



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

float complex OrderParam(float *phases)
	{
		int i;
    	double complex ord_param = 0 + 0 * I;
		for(i = 0;i < N; i++)
		{
			ord_param += cexp(I*phases[i]);
			//printf("Inside OrderParam = %.2f+I%.2f\n",creal(ord_param),cimag(ord_param));

		}
		//ord_param = (ord_param)*iN;

		return ord_param;
	}

float complex FreqOrderParam(float *ang_freqs)
	{
		int i;
    	double complex freq_ord_param = 0 + 0 * I;
		for(i = 0;i < N; i++)
		{
			freq_ord_param += cexp(I*ang_freqs[i]);

		}
		//ord_param = (ord_param)*iN;

		return freq_ord_param;
	}


void EulerStep(float *phases, float *ang_freqs, float K)
{
	int i,j;

	float phase_updated[N] = {0} ; 
	float frequencies_updated[N] = {0} ; 
	float sum_term = 0.;
	float floatN;
	floatN = (float)N;
  	for(i = 0;i < N; i++)
  	{	
  		sum_term = 0;
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

const char * CreateResultsFolder(){
	char path_results[PATH_MAX];
   	if (getcwd(path_results, sizeof(path_results)) != NULL) {
       //printf("Current working dir: %s\n", path_results);
   	} else {
       perror("getcwd() error");
    }
   	strcat(path_results, "/results");
    int ret;
    ret = mkdir(path_results,0777); //creates folder
    //printf("\nPATH IS NOW:\n%s\n\n",path_results);

	return path_results;
}
