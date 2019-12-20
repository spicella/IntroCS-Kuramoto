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

//-------------------------Begin Definitions-------------------------//
	#define turn_angle  2.*M_PI
//Main parameters
	#define N 1000	 //Number of Kuramoto oscillators
	
	#define dt .001 //Time step
	#define T 30000 //End of simulation time

//For fixed value of K-s in simulation
	#define K1 1.
	#define K2 .1
//For sweeping of K
	#define K0 0.
	#define dK .1
	#define K_max 5. //
	#define PATH_MAX 1000
//---------------------------End Definitions-------------------------//

//Run bools:
	//check values for initial configurations:
	bool check_initial = false;
	//Gaussian distributed frequencies [N(0,1)]
	bool gaussian_frequencies = true;
	//ODE + sweeping K
	bool sweeping_K = false;
	//ODE + constant K-s


//-------------------------Functions Declaration-------------------------//
void CreateResultsFolder();
void PrintParams(float K_run);
float * RandUnifPhase();
float * RandUnifFreq();
float * RandGauss();
float PeriodicPosition(float angular_pos);
void EulerStep(float *phases, float *ang_freqs, float K, double o_par[]);
void OrderParamMod(float *phases, double o_param[]);
void ClearResultsFile(float K);
void WriteResults(double o_par[], float K, float t_loop);
//----------------------------------------------------------------------//


float iN=(float)(1/N);

//Functions declarartion:

int main(void)    
{	
	clock_t begin_main = clock();
	srand(time(NULL));
	sleep(0);
	CreateResultsFolder();

	int i,j;
	float complex iN = (float)(1/N)+I*0; //inverse of N

	float K_range = K_max - K0;
	int number_k_steps = K_range/dK+1;

	//--------------------START K LOOP--------------------//
	for(j=0;j<number_k_steps;j++)
	{
		printf("RUNNING %d/%d element in K loop\nK=%.2f",j,number_k_steps,j*dK);
		clock_t begin_loop = clock();

		float K_run = K0+j*dK;
		//printf("%s", results_path);
		PrintParams(K_run);

		//Declarations

		float *phases;
		float *ang_freqs;
		double ord_param[2] = {0};		

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
		int T_split = (int)(T/100);
		ClearResultsFile(K_run);
		for(i=0;i<T+1;i++){
			OrderParamMod(phases,ord_param);
			WriteResults(ord_param, K_run,i);
			EulerStep(phases, ang_freqs, K_run,ord_param);
			
			if(i%T_split==0){
				printf("\n\tProcess at %d/100, K=%.2f\n", 100*(int)(i)/T,K_run);
				printf("\tOrderParameter Modulus = %.5f\n",ord_param[0]);
				}

		}
		//----------------------END SINGLE RUN LOOP----------------------//

		//List of input parameters
		PrintParams(K_run);
		clock_t inner_loop_end = clock();
		double loop_time_spent = (double)(inner_loop_end - begin_loop) / CLOCKS_PER_SEC;
		printf("Loop time: %.5f seconds\n\n",loop_time_spent);

	}
	clock_t outer_loop_end = clock();

	double total_time_spent = (double)(outer_loop_end - begin_main) / CLOCKS_PER_SEC;

	printf("Total execution time: %.5f seconds\n\n",total_time_spent);
  	return 0;
}






//-----------------------------------Function implementations-----------------------------------//

void PrintParams(float K_run){
	printf("\n\n\n//------------------------Input parameters------------------------//\n");
	printf("\t-------> Number of Oscillators = %d\n",N);
	printf("\t-------> Runtime simulation = %d\n",T);
	printf("\t-------> Start K = %.2f\n",K0);
	printf("\t-------> K = %.2f\n",K_run);
	printf("\t-------> dK = %.2f\n",dK);
	printf("\t-------> End K = %.2f\n",K_max);
	printf("\t-------> dt = %.5f\n",dt);
	printf("\t-------> Gaussian initial frequencies? = %s\n", gaussian_frequencies ? "True!" : "False!");
	printf("\t-------> Check initial conditions? = %s\n", check_initial ? "True!" : "False!");
	printf("//----------------------------------------------------------------//\n\n\n");
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

float * RandUnifFreq( ){
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

float * RandGauss(  ){
	/*https://www.doc.ic.ac.uk/~wl/papers/07/csur07dt.pdf for gaussian array generation*/
	/*Box-Muller transform*/
	srand( (unsigned)time(NULL) );

    static float r[N];
    float U1,U2,sqrt_term,phase;
    int i;
    for ( i = 0; i < N/2; ++i){
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

void OrderParamMod(float *phases, double o_par[])
{
		int i;
    	double real_ord_param = 0;
    	double imag_ord_param = 0;
    	double iN = 1./((double)N);
		for(i = 0;i < N; i++)
		{
			real_ord_param += iN*cos(phases[i]);
			imag_ord_param += iN*sin(phases[i]);
		}

		o_par[0] = sqrt(real_ord_param*real_ord_param+imag_ord_param*imag_ord_param); //modulus
 		o_par[1] = 2*atan(imag_ord_param/(o_par[0]+real_ord_param));
}

//Frequency order parameter!


void EulerStep(float *phases, float *ang_freqs, float K, double o_par[]){
	int i;
	double iN = 1./((double)N);
  	for(i = 0;i < N; i++)
  	{	
  		ang_freqs[i] +=  K*o_par[0]*sin(o_par[1]-phases[i]);
  		phases[i] += dt*ang_freqs[i];
  		phases[i] = PeriodicPosition(phases[i]);
  	}
}

void CreateResultsFolder(){
	char path_results[PATH_MAX];
   	if (getcwd(path_results, sizeof(path_results)) != NULL) {
       //printf("Current working dir: %s\n", path_results);
   	} else {
       perror("getcwd() error");
    }
   	strcat(path_results, "/results");
    int ret;
    ret = mkdir(path_results,0777); //creates folder
}

void ClearResultsFile(float K){
		/*Remove File with the same name, avoid overwriting*/
		char filename[64];
		FILE *out;
		if(gaussian_frequencies==true){
			sprintf(filename, "results/gfreq_N%d_T%d_dt%.8f_K%.2f.tsv", N,T,dt,K);
		}
		else{
			sprintf(filename, "results/ufreq_N%d_T%d_dt%.8f_K%.2f.tsv", N,T,dt,K);
		}		if (remove(filename) == 0) 
      		printf("Deleted successfully"); 
   		else
      		printf("Unable to delete the file"); 
}

void WriteResults(double o_par[], float K, float t_loop){
		/*Single shot writing of order param and freq order param (both real and complex)*/
		int i;
		char filename[64];
		FILE *out;
		if(gaussian_frequencies==true){
			sprintf(filename, "results/gfreq_N%d_T%d_dt%.8f_K%.2f.tsv", N,T,dt,K);
		}
		else{
			sprintf(filename, "results/ufreq_N%d_T%d_dt%.8f_K%.2f.tsv", N,T,dt,K);
		}
		out = fopen( filename, "a");
		//fprintf(out, "%.20f\t%.20f\t%.20f\t%.20f\n", creal(ord_param)/N,cimag(ord_param)/N,creal(freq_ord_param)/N,cimag(freq_ord_param)/N);
		fprintf(out,"%.5f\t%.20f\t%.20f\n",t_loop*dt, o_par[0],o_par[1]);
		fclose(out);

}
