/*
 * Discrete Yang-Mills -- multi-character Metropolis MC that also saves SU(3)
 * configurations in NERSC format for Wilson-flow scale setting (w0/t0).
 *
 * Identical Monte Carlo to dym-mod-metro.cpp (multi-character action
 *   S = -sum_plaq sum_c beta_c Re chi_c(U_plaq) ),
 * but after each measurement it writes the configuration as SU(3) matrices in
 * NERSC (IEEE64BIG) byte order, directly readable by tools/wilsonflow/wilson_flow.
 *
 * Requires a group file with the 18-real (3x3 complex) defining-rep matrices
 * APPENDED after the char-table block (P C / ReChar / ImChar / mult / conclass),
 * e.g. groups/S1080ctm (built by generators/build_s1080_matrix_file.py).
 *
 * usage: ./dym-mod-metro-savecfg group D Nt Nx beta_0 ... beta_{C-1} seed \
 *                                [outprefix] [K] [N] [Ntherm]
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
#include <vector>
#include <fstream>
#include "timer.h"

#define K_DEFAULT       200    // decorrelation sweeps between saved configurations
#define N_DEFAULT       1000   // number of configurations to save
#define NTHERM_DEFAULT  0      // thermalization sweeps
#define NHIT            20     // Metropolis hits per link per visit

typedef unsigned int uint;

extern "C" void getpos(unsigned idx, unsigned int *pos);   // in lattice.c
void update();

std::default_random_engine *rnd;
std::vector<group_t> smallgroup;
std::vector<double> beta;
std::vector<int> idbeta;
char *site_parity = NULL;
unsigned long hit = 0, acc = 0;
const char *outprefix = "./";

/* Defining-rep matrix of each group element: data[part + 2*(c + 3*r)]. */
struct su3 { double data[18]; double& operator()(int r, int c, int part) { return data[part+2*(c+3*r)]; } };
std::vector<su3> groupsu3;

void switchend(unsigned char *buffer, int length)
{
  unsigned char *pos = buffer, save[8];
  for(int j=0; j<length; j++) {
    for(int i=0; i<8; i++) save[i] = pos[i];
    for(int i=7; i>-1; i--, pos++) *pos = save[i];
  }
}

/* Current config -> NERSC 4D_SU3_GAUGE_3X3 file (x-fastest, links x,y,z,t).
 * wilson_flow re-validates LINK_TRACE and PLAQUETTE on load, so a wrong index
 * mapping is caught immediately. Our ReTr stores -Re Tr (defining rep = row 1),
 * so the actual normalized observables use -ReTr. */
void save_nersc(int n)
{
  double sumlink = 0.0;
  for (unsigned i=0;i<V;++i) for (unsigned d=0;d<D;++d) sumlink += -ReTr[a[i*D+d]];
  double link_trace = sumlink / (3.0 * D * V);
  double sumplaq = 0.0;
  for (unsigned i=0;i<V;++i) for (unsigned d1=0;d1<D;++d1) for (unsigned d2=d1+1;d2<D;++d2)
    sumplaq += -ReTr[plaquette(i,d1,d2)];
  double plaq = sumplaq / (3.0 * (V*D*(D-1)/2));

  char name[512];
  double bf = (beta.size()>1)?beta[1]:0.0, ba = (beta.size()>3)?beta[3]:0.0;
  snprintf(name, sizeof(name), "%snersc-bf%.3f-ba%.3f-nt%02d-nx%02d-num%04d",
           outprefix, bf, ba, Nt, Nx, n);
  FILE* f = fopen(name, "w");
  if(!f) { fprintf(stderr, "cannot open %s for writing\n", name); abort(); }
  fprintf(f, "BEGIN_HEADER\n");
  fprintf(f, "HDR_VERSION = 1.0\n");
  fprintf(f, "DATATYPE = 4D_SU3_GAUGE_3X3\n");
  fprintf(f, "DIMENSION_1 = %u\n", Nx);
  fprintf(f, "DIMENSION_2 = %u\n", Nx);
  fprintf(f, "DIMENSION_3 = %u\n", Nx);
  fprintf(f, "DIMENSION_4 = %u\n", Nt);
  fprintf(f, "CHECKSUM = ffffffff\n");
  fprintf(f, "LINK_TRACE = %20.15f\n", link_trace);
  fprintf(f, "PLAQUETTE = %20.15f\n", plaq);
  fprintf(f, "ENSEMBLE_ID = discrete-multichar\n");
  fprintf(f, "FLOATING_POINT = IEEE64BIG\n");
  fprintf(f, "END_HEADER\n");

  const int dirs[4] = {1, 2, 3, 0};   /* file order x,y,z,t -> our dir labels */
  su3 m;
  for (int t=0; t<(int)Nt; ++t)
  for (int z=0; z<(int)Nx; ++z)
  for (int y=0; y<(int)Nx; ++y)
  for (int x=0; x<(int)Nx; ++x)
  {
    unsigned idx = t + Nt*(x + Nx*(y + Nx*z));     /* pos = [t,x,y,z] */
    for (int k=0; k<4; ++k) {
      m = groupsu3[a[idx*D + dirs[k]]];
      switchend((unsigned char*)m.data, 18);
      fwrite(m.data, sizeof(double), 18, f);
    }
  }
  fclose(f);
}

/* 2(D-1) plaquette staples touching link (i,d) -- multi-character action. */
void get_staples(unsigned int i, unsigned int d, group_t* st)
{
  unsigned int i1 = step(i, d, 1); int k=0;
  for(unsigned int d1=0; d1<D; ++d1) if(d1 != d)
  {
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

/* One checkerboard Metropolis sweep, multi-character action. */
void update()
{
  std::uniform_real_distribution<> rand01(0., 1.);
  const int nstaple = 2 * (D - 1);
  for (unsigned d = 0; d < D; d++)
  for (int parity = 0; parity < 2; ++parity)
  {
    long lhit = 0, lacc = 0;
#pragma omp parallel for schedule(static) reduction(+: lacc, lhit)
    for (uint i = 0; i < V; i++)
    {
      if (site_parity[i] != parity) continue;
      group_t staples[2 * (D - 1)];  get_staples(i, d, staples);
      for (int h = 0; h < NHIT; ++h)
      {
        double r = 0;
        group_t b = a[i*D+d];
        group_t bnew = mult[b][ smallgroup[ smallgroup.size() * rand01(rnd[i]) ] ];
        for (int j = 0; j < nstaple; ++j)
        {
          group_t p1 = mult[b][staples[j]];
          group_t p2 = mult[bnew][staples[j]];
          for (unsigned cidx = 0; cidx < idbeta.size(); ++cidx) {
            unsigned cc = idbeta[cidx];
            r += beta[cc]*(ReChar[cc][conclass[p1]] - ReChar[cc][conclass[p2]]);
          }
        }
        ++lhit;
        if (exp(r) > rand01(rnd[i])) { a[i*D+d] = bnew; ++lacc; }
      }
    }
    hit += lhit; acc += lacc;
  }
}

int main(int argc, char *argv[])
{
  if (argc < 2) { fprintf(stderr, "usage: %s group D Nt Nx betas[C] seed [outprefix] [K] [N] [Ntherm]\n", argv[0]); return 1; }
  const char *groupfilename = argv[1];
  load_group(groupfilename);                 /* sets P, C, ReChar, ReTr, mult, ... */

  if (argc < 6 + (int)C) {
    fprintf(stderr, "usage: %s group D Nt Nx betas[%d] seed [outprefix] [K] [N] [Ntherm]\n", argv[0], C);
    return 1;
  }
  D  = atoi(argv[2]); Nt = atoi(argv[3]); Nx = atoi(argv[4]);
  for (uint i = 0; i < C; ++i) {
    beta.push_back(atof(argv[5+i]));
    if (fabs(beta[i]) > 1e-10) idbeta.push_back(i);
  }
  int iseed   = atoi(argv[5+C]);
  if (argc > 6+(int)C) outprefix = argv[6+C];
  unsigned K      = (argc > 7+(int)C) ? (unsigned)atoi(argv[7+C]) : K_DEFAULT;
  unsigned N      = (argc > 8+(int)C) ? (unsigned)atoi(argv[8+C]) : N_DEFAULT;
  unsigned NTHERM = (argc > 9+(int)C) ? (unsigned)atoi(argv[9+C]) : NTHERM_DEFAULT;

  V = Nt; for (unsigned d = 1; d < D; ++d) V *= Nx;
  printf("PARAMS(grp,D,Nt,Nx,seed): %s %d %d %d %d\n", groupfilename, D, Nt, Nx, iseed);
  printf("RUN(K,N,Ntherm,NHIT): %u %u %u %d\n", K, N, NTHERM, NHIT);
  printf("classes: %d   nonzero betas: %lu\n", C, idbeta.size());
  printf("saving NERSC configs to: %s\n", outprefix);

  /* Appended defining-rep matrices: skip the char-table header, then read P*18. */
  {
    std::ifstream ing(groupfilename);
    int p_, c_; ing >> p_ >> c_;
    double dtmp; long itmp;
    for (long i = 0; i < 2L*c_*c_; ++i) ing >> dtmp;     /* ReChar, ImChar */
    for (long i = 0; i < (long)p_*p_; ++i) ing >> itmp;  /* mult */
    for (int i = 0; i < p_; ++i) ing >> itmp;            /* conclass */
    groupsu3.resize(p_);
    for (int i = 0; i < p_; ++i) { double *pm = groupsu3[i].data; for (int j = 0; j < 18; ++j, ++pm) ing >> *pm; }
    if (ing.fail()) {
      fprintf(stderr, "group file '%s' has no appended SU(3) matrices; cannot save configs\n", groupfilename);
      abort();
    }
  }

  rnd = new std::default_random_engine[V];
  for (int i = 0; i < V; ++i) rnd[i].seed(iseed + i);
  std::uniform_int_distribution<> randgrp(0, P-1);
  printf("read the group; id: %i  assoc: %d %d\n", id, mult[mult[1][2]][3], mult[1][mult[2][3]]);

  /* Proposal pool: elements at the smallest non-identity ReTr (near the identity). */
  double min_retr = 0, nn_retr = 10;
  for (uint i = 0; i < P; ++i) if (ReTr[i] < min_retr) min_retr = ReTr[i];
  for (uint i = 0; i < P; ++i) if (ReTr[i] > min_retr && ReTr[i] < nn_retr) nn_retr = ReTr[i];
  for (uint i = 0; i < P; ++i) if (ReTr[i] < nn_retr + 1e-6 && ReTr[i] > min_retr) smallgroup.push_back(i);
  printf("small group size: %lu\n", smallgroup.size());

  a = (group_t*) malloc(sizeof(*a) * V * D);
  for (unsigned i = 0; i < V*D; i++) a[i] = randgrp(rnd[0]);
  step(0, 0, 1);                              /* prime the nn table */

  site_parity = (char*) malloc(V);
  for (unsigned i = 0; i < V; ++i) {
    unsigned pos[4]; getpos(i, pos);
    int s = 0; for (unsigned d = 0; d < D; ++d) s += pos[d];
    site_parity[i] = s & 1;
  }
  printf("init done\n"); fflush(stdout);

  for (unsigned k = 0; k < NTHERM; k++) update();
  printf("thermo done\n"); fflush(stdout);

  { /* reweighting log header: action is linear in beta, S=-sum_c beta_c P_c(U) */
    char pf[600]; snprintf(pf, sizeof(pf), "%spchar.dat", outprefix);
    FILE* ph = fopen(pf, "w");
    if (ph) {
      fprintf(ph, "# NCHAR %d\n# BASE_BETA", C);
      for (unsigned c = 0; c < C; ++c) fprintf(ph, " %.10g", beta[c]);
      fprintf(ph, "\n# cols: num  P_0..P_%d  simpleplaq repoly impoly\n", C-1);
      fclose(ph);
    }
  }

  for (unsigned n = 0; n < N; n++)
  {
    timer tm; tm.start("update");
    for (unsigned k = 0; k < K; k++) update();
    tm.stop();

    tm.start("meas");
    double rep, imp; getpoly(&rep, &imp);
    double simpleplaq = 0;
#pragma omp parallel for reduction(+: simpleplaq)
    for (unsigned int i = 0; i < V; ++i)
    for (unsigned int d1 = 0; d1 < D; ++d1)
    for (unsigned int d2 = d1+1; d2 < D; ++d2)
      simpleplaq += -ReTr[plaquette(i, d1, d2)];
    simpleplaq /= V*D*(D-1)/2;

    /* fundamental Wilson-loop grid W(r,t) for the static potential / string tension */
    double wloop[Nx][Nt];
    for (unsigned int k = 0; k < Nx; ++k)
    for (unsigned int l = 0; l < Nt; ++l) {
      double lwloop = 0;
#pragma omp parallel for reduction(+: lwloop)
      for (unsigned int i = 0; i < V; ++i)
      for (unsigned int d1 = 1; d1 < D; ++d1)
        lwloop += -ReTr[wilson(i, d1, 0, k+1, l+1)];
      wloop[k][l] = lwloop / (V*(D-1));
    }

    printf("GMES: %e %e %e %e", 999.0, rep, imp, simpleplaq);
    for (unsigned int k = 0; k < Nx; ++k)
    for (unsigned int l = 0; l < Nt; ++l) printf(" %e", wloop[k][l]);
    printf("\n");
    printf("ACC: %f\n", ((double)acc)/hit);
    fflush(stdout);

    { /* per-config character plaquette sums P_c = sum_plaq Re chi_c(U_p), for reweighting */
      double Pc[CMAX]; for (unsigned c = 0; c < C; ++c) Pc[c] = 0.0;
      for (unsigned i = 0; i < V; ++i)
      for (unsigned d1 = 0; d1 < D; ++d1)
      for (unsigned d2 = d1+1; d2 < D; ++d2) {
        int cl = conclass[plaquette(i, d1, d2)];
        for (unsigned c = 0; c < C; ++c) Pc[c] += ReChar[c][cl];
      }
      char pf[600]; snprintf(pf, sizeof(pf), "%spchar.dat", outprefix);
      FILE* ph = fopen(pf, "a");
      if (ph) {
        fprintf(ph, "%u", n);
        for (unsigned c = 0; c < C; ++c) fprintf(ph, " %.10e", Pc[c]);
        fprintf(ph, " %.10e %.10e %.10e\n", simpleplaq, rep, imp);
        fclose(ph);
      }
    }

    save_nersc(n);
    tm.stop();
  }

  free(a); free(site_parity); delete[] rnd;
  return 0;
}
