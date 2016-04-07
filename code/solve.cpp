/*
 * STRUMPACK -- STRUctured Matrices PACKage, Copyright (c) 2014, The Regents of
 * the University of California, through Lawrence Berkeley National Laboratory
 * (subject to receipt of any required approvals from the U.S. Dept. of Energy).
 * All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Technology Transfer Department at TTD@lbl.gov.
 *
 * NOTICE. This software is owned by the U.S. Department of Energy. As such, the
 * U.S. Government has been granted for itself and others acting on its behalf a
 * paid-up, nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and display publicly.
 * Beginning five (5) years after the date permission to assert copyright is
 * obtained from the U.S. Department of Energy, and subject to any subsequent five
 * (5) year renewals, the U.S. Government is granted for itself and others acting
 * on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
 * Software to reproduce, prepare derivative works, distribute copies to the
 * public, perform publicly and display publicly, and to permit others to do so.
 *
 * Developers: Francois-Henry Rouet, Xiaoye S. Li, Pieter Ghysels
 *             (Lawrence Berkeley National Lab, Computational Research Division).
 *
 */

/* Solution of a linear system with STRUMPACK */
#include "StrumpackDensePackage.hpp"

#define myscalar double
#define myreal double

int main (int argc, char *argv[]) {
  myscalar *A=NULL, *X=NULL, *B=NULL;
  int descA[BLACSCTXTSIZE], descXB[BLACSCTXTSIZE];
  int n;
  int nrhs;
  int nb;
  int locr, locc;
  int i, j, ii, jj;
  int ierr;
  int dummy;
  int myid, np;
  int myrow, mycol, nprow, npcol;
  int ctxt;

  n=80000; /* Size of the problem */
  nrhs=1; /* Number of RHS */
  nb=16;  /* Blocksize for the 2D block-cyclic distribution */

  /* Initialize MPI */
  if((ierr=MPI_Init(&argc,&argv)))
    return 1;
  myid=-1;
  if((ierr=MPI_Comm_rank(MPI_COMM_WORLD,&myid)))
    return 1;
  np=-1;
  if((ierr=MPI_Comm_size(MPI_COMM_WORLD,&np)))
    return 1;

  /* Initialize the BLACS grid */
  nprow=floor(sqrt((float)np));
  npcol=np/nprow;
  blacs_get_(&IZERO,&IZERO,&ctxt);
  blacs_gridinit_(&ctxt,"R",&nprow,&npcol);
  blacs_gridinfo_(&ctxt,&nprow,&npcol,&myrow,&mycol);

  /* A is a dense n x n distributed Toeplitz matrix */
  if(myid<nprow*npcol) {
    locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
    locc=numroc_(&n,&nb,&mycol,&IZERO,&npcol);
    A=new myscalar[locr*locc];
    dummy=std::max(1,locr);
    descinit_(descA,&n,&n,&nb,&nb,&IZERO,&IZERO,&ctxt,&dummy,&ierr);

    for(i=1;i<=locr;i++)
      for(j=1;j<=locc;j++) {
        ii=indxl2g_(&i,&nb,&myrow,&IZERO,&nprow);
        jj=indxl2g_(&j,&nb,&mycol,&IZERO,&npcol);
        // Toeplitz matrix from Quantum Chemistry.
        myreal pi=3.1416, d=0.1;
        A[locr*(j-1)+(i-1)]=ii==jj?std::pow(pi,2)/6.0/std::pow(d,2):std::pow(-1.0,ii-jj)/std::pow((myreal)ii-jj,2)/std::pow(d,2);

      }
  } else {
    descset_(descA,&n,&n,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
  }

  /* Initialize the solver and set parameters */
  StrumpackDensePackage<myscalar,myreal> sdp(MPI_COMM_WORLD);
  sdp.use_HSS=true;
  sdp.levels_HSS=13;
  sdp.min_rand_HSS=150;
  sdp.lim_rand_HSS=5;
  sdp.inc_rand_HSS=10;
  sdp.max_rand_HSS=150;
  sdp.tol_HSS=1e-4;
  sdp.steps_IR=10;
  sdp.tol_IR=1e-4;

  /* Compression */
  sdp.compress(A,descA);

  /* Accuracy checking */
  /* sdp.check_compression(A,descA); */

  /* Factorization */
  sdp.factor(A,descA);

  /* Set the RHS (random vector) and the solution space */
  if(myid<nprow*npcol) {
    locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
    locc=numroc_(&nrhs,&nb,&mycol,&IZERO,&npcol);
    B=new myscalar[locr*locc];
    dummy=std::max(1,locr);
    descinit_(descXB,&n,&nrhs,&nb,&nb,&IZERO,&IZERO,&ctxt,&dummy,&ierr);
    for(i=0;i<locr*locc;i++) {
      B[i]=static_cast<myreal>(rand())/(static_cast<myreal>(RAND_MAX));
    }
    X=new myscalar[locr*locc];
  } else {
    descset_(descXB,&n,&nrhs,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
  }

  /* Triangular solution */
  sdp.solve(X,descXB,B,descXB);

  /* Accuracy checking */
  sdp.check_solution(A,descA,X,descXB,B,descXB);

  /* Iterative refinement */
  sdp.refine(A,descA,X,descXB,B,descXB);

  /* Statistics */
  sdp.print_statistics();

  /* Clean-up */
  delete[] A;
  delete[] B;
  delete[] X;

  /* The end */
  MPI_Finalize();
  return 0;

}
