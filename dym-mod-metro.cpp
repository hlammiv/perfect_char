/*
 * Discrete Yang-Mills
 */

extern "C" {
#include "group.h"
#include "lattice.h"
}

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <random>
#include "timer.h"

#define K 1000      // Decorrelation time
#define N 10000   // Number of samples

typedef unsigned int uint;

void update();
double S_inv(unsigned i, unsigned d);

std::default_random_engine *rnd;

std::vector<group_t> smallgroup;
std::vector<double> beta;
std::vector<int> idbeta;
unsigned long hit = 0, acc = 0;

int main(int argc, char *argv[]) {
	/* Seed the PRNG. */
	int iseed = 0; //time(NULL) + getpid();

	/* Theory parameters come from the command line. */
	const char *groupfilename = argv[1];

	/* Load the multiplication table. */
	load_group(groupfilename);
        if (argc < 5+C) {
                fprintf(stderr, "usage: %s group D Nt Nx betas[as list]\n", argv[0]);
                return 1;
        }

        D = atoi(argv[2]);
        Nt = atoi(argv[3]);
        Nx = atoi(argv[4]);
	V = Nt*pow(Nx,D-1);

	for(uint i=0; i<C; ++i){
		beta.push_back(atof(argv[5+i]));
		if(fabs(beta[i])>1e-10) idbeta.push_back(i);
	}

	if(argc == 6+C) iseed = atoi(argv[5+C]);
        printf("PARAMS(grp,D,Nt,Nx,seed): %s %d %d %d %d\n", groupfilename, D, Nt, Nx, iseed);

	for(uint i = 0; i<C; ++i){ 
		printf("Beta[%d]: %e   ",i,beta[i]);
		if(i%2==1) printf("\n");
	}
	printf("Nonzero Beta\n");
       for(uint i = 0; i<idbeta.size(); ++i){
                printf("Beta[%d]: %e   ",idbeta[i],beta[idbeta[i]]);
                if(i%2==1) printf("\n");
        }


	printf("\nread the group\n"); fflush(stdout);
	printf("id: %i\n",id);
	printf("check stuff: (a*b)*c: %d a*(b*c): %d\n", mult[mult[1][2]][3], mult[1][mult[2][3]]);
	printf("classes: %d\n",C);


	rnd = new std::default_random_engine[V];
	for(int i=0; i<V; ++i) rnd[i].seed(iseed+i);
	std::uniform_int_distribution<> randgrp(0,P-1);


	// select the group elems of S1080 that are close to identity
	double min_retr=0;
        double nn_retr=0;
        for(uint i=0; i<P; ++i) if(ReTr[i] < min_retr) min_retr=ReTr[i];
        for(uint i=0; i<P; ++i) if(ReTr[i] > min_retr && ReTr[i] < nn_retr) nn_retr=ReTr[i];
        for(uint i=0; i<P; ++i) if(ReTr[i] < nn_retr+1e-6 && ReTr[i] > min_retr) smallgroup.push_back(i);
	printf("min_retr: %e nn_retr: %e \n",min_retr,nn_retr);
	printf("small group size: %lu\n", smallgroup.size());

	/* Initialize the gauge field. */
	a = (group_t*) malloc(sizeof(unsigned) * V * D);
	for (unsigned i = 0; i < V*D; i++){
			a[i] = randgrp(rnd[0]);
//			a[i] = rand()%P;
//			a[i] = id;
//			if(i%2==0){a[i] = rand()%P;}else{a[i] = id;}
	}
step(0,0,1); //prime the nn table
printf("init done\n"); fflush(stdout);

//	for (unsigned i = 0; i < V*D; i++)
//                        printf(" %03d", a[i]);
//                printf("\n");


	for (unsigned k = 0; k < K*0; k++) update();
printf("thermo done\n"); fflush(stdout);

  for (unsigned n = 0; n < N; n++) 
  {
    timer tm;
    tm.start("update");
    for (unsigned k = 0; k < K; k++) update();
    tm.stop();
#if 0
    for (unsigned i = 0; i < V*D; i++)
      printf(" %03d", a[i]);
    printf("\n");
    fflush(stdout);
#endif
    // update done -- print plaq
    tm.start("meas");
//    double plaq = 0;
    
    int Ncc=beta.size();
    double rep, imp; getpoly(&rep, &imp);
    // plaq computed is actually 4*action
    double simpleplaq = 0;
    double wloop[Nx][Nt];
    double wccloop[Ncc][Nx][Nt];
//    for(unsigned int i=0; i<Nx; ++i)
//    for(unsigned int j=0; j<Nt; ++j)
//    {
//	    wloop[i][j]=0;
//    }


    for(unsigned int k=0; k<Ncc; ++k)
    for(unsigned int i=0; i<Nx; ++i)
    for(unsigned int j=0; j<Nt; ++j)
    {
            wccloop[k][i][j]=0;
    }


#pragma omp parallel for reduction(+: simpleplaq)
    for(unsigned int i=0; i<V; ++i)
    for(unsigned int d1=0; d1<D; ++d1)
    for(unsigned int d2=d1+1; d2<D; ++d2)
    {
	    //printf("%e \n",ReTr[wilson(i,d1,d2,3,2)]);
	simpleplaq += -ReTr[plaquette(i, d1, d2)];
    }
   
/*    for(unsigned int k=0; k<Nx; ++k)
    for(unsigned int l=0; l<Nt; ++l)
    {
	double lwloop=0;
	#pragma omp parallel for reduction(+: lwloop)
	for(unsigned int i=0; i<V; ++i)
    	for(unsigned int d1=1; d1<D; ++d1)
    	{
    		lwloop+=-ReTr[wilson(i,d1,0,k+1,l+1)];
    	}
    	wloop[k][l]=lwloop/(V*D*(D-1)/2);
    }
*/  
    for(unsigned int k=0; k<Nx; ++k)
    for(unsigned int l=0; l<Nt; ++l)
    for(unsigned int m=0; m<Ncc; ++m)
    {
        double lwloop=0;
        #pragma omp parallel for reduction(+: lwloop)
        for(unsigned int i=0; i<V; ++i)
        for(unsigned int d1=1; d1<D; ++d1)
        {
                lwloop+=-ReChar[m][conclass[wilson(i,d1,0,k+1,l+1)]];
        }
        wccloop[m][k][l]=lwloop/(V*D*(D-1)/2);
    }

    simpleplaq /= V*D*(D-1)/2;
    double cor;
    printf("GMES: %e %e %e %e", 999.0, rep, imp, simpleplaq);
	printf("\n");
    for(unsigned int m=0; m<Ncc; ++m)
	{
    	for(unsigned int k=0; k<Nx; ++k)
    	{
		printf("CC %d WL %d: ",m,k);
    		for(unsigned int l=0; l<Nt; ++l) printf(" %e",wccloop[m][k][l]);
		printf("\n");
    	}
	}
//    printf("\n");
    fflush(stdout);
    printf("ACC: %f\n", ((double)acc)/hit);

    
    printf("CONFIGS: ");
      for (unsigned i = 0; i < V*D; i++)
                        printf(" %03d", a[i]);
                printf("\n");

    tm.stop();
  }

	free(a);
	return 0;
}

void get_staples(unsigned int i, unsigned int d, group_t* st)
{
  unsigned int i1 = step(i, d, 1); int k=0;
  for(unsigned int d1=0; d1<D; ++d1) if(d1 != d)
  {
//printf("(i,d,d1): %d %d %d\n", i, d, d1);
//printf("[123]: %d %d %d\n", a[i1*D+d1], inv[a[step(i,d1,1)*D+d]], inv[a[i*D+d1]]);
    group_t g = a[i1*D+d1];
    g = mult[g][inv[a[step(i,d1,1)*D+d]]];
    g = mult[g][inv[a[i*D+d1]]];
    st[k++] = g;
    g = inv[a[step(i1,d1,-1)*D+d1]];
    int i2 = step(i,d1,-1);
    g = mult[g][inv[a[i2*D+d]]];
    g = mult[g][a[i2*D+d1]];
    st[k++] = g;
  }
}

extern "C" void getpos(unsigned idx, unsigned int *pos);
void update() 
{
  timer tm(false); tm.start("update");
  std::uniform_int_distribution<> randsmallgrp(0,smallgroup.size()-1);
  std::uniform_real_distribution<> rand01(0.,1.);
  for (unsigned d = 0; d < D; d++) 
  for (int parity = 0; parity < 2; ++parity)
  {
    int lhit=0, lacc=0;
#pragma omp parallel for schedule(static), reduction(+: lacc, lhit)
  for (uint i = 0; i < V; i++) 
  {
    uint pos[4+00];
    getpos(i, pos); if( (pos[0]+pos[1]+pos[2]+pos[3])%2 != parity) continue;
    group_t staples[2*(D-1)+00]; 
    get_staples(i, d, staples);
    int nhit = 20;
    for(int h=0; h<nhit; ++h)
    {
      double r=0;
      group_t b = a[i*D+d];
//      group_t bnew = mult[b][smallgroup[randsmallgrp(rnd[i])]];
      group_t bnew = mult[b][smallgroup[smallgroup.size()*rand01(rnd[i])]];
      for(int j=0; j<2*(D-1); ++j) 
      {
        group_t p1 = mult[b][staples[j]];
        group_t	p2 = mult[bnew][staples[j]];
	for(unsigned cidx = 0; cidx<idbeta.size(); ++cidx){
		unsigned cc=idbeta[cidx];
        	r += beta[cc]*(ReChar[cc][conclass[p1]]-ReChar[cc][conclass[p2]]);
        }

      }
      double probrat = exp(r);
//#pragma omp atomic
//      hit++;
      lhit++;
      if(probrat > rand01(rnd[i])) { 
	a[i*D+d] = bnew; 
//#pragma omp atomic
//	acc++; 
	lacc++; 
//	if(hit%(100*V*D)==0) printf("ACC: %f\n", ((double)acc)/hit);
      }
    }
  }

  hit += lhit;
  acc += lacc;
  }

  tm.stop();
}
