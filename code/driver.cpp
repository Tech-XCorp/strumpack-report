/* Compressing and solving a system with a BEM matrix.
 * Collaboration with Ana Manic (Colorado State University).
 */

#include <fstream>
#include "StrumpackDensePackage.hpp"

#define myscalar scomplex
#define myreal float

int main (int argc, char *argv[]) {
  /* A simple driver that reads a BEM matrix (and RHS) from Col State.
   * The matrix is split into 64 files (256 for example 2) corresponding
   * to an 8x8 (16x16 for example 2) non-cyclic grid of processes.
   * The RHS (one column) is split into 8 files (16 for example 2).
   *
   * Parameters:
   *  - argv[1] use HSS? true of false.
   *  - argv[2] number of levels for HSS tree.
   *  - argv[3] number of random vectors.
   *  - argv[4] tolerance for compression.
   *
   */

  myscalar *Atmp=NULL, *A=NULL, *Btmp=NULL, *X=NULL, *B=NULL;
  int descA[BLACSCTXTSIZE], descAtmp[BLACSCTXTSIZE], descXB[BLACSCTXTSIZE], descBtmp[BLACSCTXTSIZE];
  int nrhs=1;
  int nb=64;
  int locr, locc;
  int i, j;
  int ierr;
  int dummy;
  int myid, np, id;
  int myrow, mycol, nprow, npcol;
  int myrowA, mycolA;
  int ctxtA, ctxt, ctxttmp, ctxtglob;
  double tstart, tend;
  std::string filename, prefix, locfile;
  std::ifstream fp;
  scomplex stmp;

  /* Hard-coded matrices */
//  // Example 1
//#define nprowA 8
//#define npcolA 8
//  int n=58800;
//  int nrows[nprowA]={7500,7500,7500,7500,7500,7300,7000,7000};
//  int ncols[npcolA]={7500,7500,7500,7500,7500,7300,7000,7000};
//  prefix="/scratch2/scratchdirs/frouet/Ana/example1/";

  // Example 2
#define nprowA 16
#define npcolA 16
  int n=132300;
  int nrows[nprowA]={8500,8500,8500,8500,8500,8500,8500,8500,8300,8000,8000,8000,8000,8000,8000,8000};
  int ncols[nprowA]={8500,8500,8500,8500,8500,8500,8500,8500,8300,8000,8000,8000,8000,8000,8000,8000};
  prefix="/scratch2/scratchdirs/frouet/Ana/example2/";

//  // Example 3
//#define nprowA 8
//#define npcolA 8
//  int n=27648;
//  int nrows[nprowA]={3456,3456,3456,3456,3456,3456,3456,3456};
//  int ncols[npcolA]={3456,3456,3456,3456,3456,3456,3456,3456};
//  prefix="/scratch2/scratchdirs/frouet/Ana/example3/";

  int rowoffset[nprowA];
  int coloffset[npcolA];

  rowoffset[0]=1;
  for(i=1;i<nprowA;i++)
    rowoffset[i]=rowoffset[i-1]+nrows[i-1];
  coloffset[0]=1;
  for(i=1;i<npcolA;i++)
    coloffset[i]=coloffset[i-1]+ncols[i-1];

  /* Initialize MPI */
  if((ierr=MPI_Init(&argc,&argv)))
    return 1;
  myid=-1;
  if((ierr=MPI_Comm_rank(MPI_COMM_WORLD,&myid)))
    return 1;
  np=-1;
  if((ierr=MPI_Comm_size(MPI_COMM_WORLD,&np)))
    return 1;

  if(np<nprowA*npcolA) {
    std::cout << "This requires " << nprowA*npcolA << " processes or more." << std::endl;
    std::cout << "Aborting." << std::endl;
    MPI_Abort(MPI_COMM_WORLD,-1);
  }

  /* Read arguments */
  if(argc<5) {
    if(!myid)
      std::cout << "Usage: ./main use_HSS levels random_vectors tolerance" << std::endl;
    MPI_Finalize();
    return 0;
  }

  /* Initialize a BLACS grid with nprowA*npcolA processes */
  nprow=nprowA;
  npcol=npcolA;
  blacs_get_(&IZERO,&IZERO,&ctxtA);
  blacs_gridinit_(&ctxtA,"R",&nprow,&npcol);
  blacs_gridinfo_(&ctxtA,&nprow,&npcol,&myrowA,&mycolA);

  /* Processes 0..nprow*npcolA read their piece of the matrix */
  MPI_Barrier(MPI_COMM_WORLD);
  tstart=MPI_Wtime();
  if(myid<nprowA*npcolA) {
    locfile="ZZ_"+SSTR(myrowA)+"_"+SSTR(mycolA)+"_"+SSTR(nrows[myrowA])+"_"+SSTR(ncols[mycolA]);
    filename=prefix+locfile;
    std::cout << "Process " << myid << " reading from file " << locfile << std::endl;
    fp.open(filename.c_str(),std::ios::binary);
    if(!fp.is_open()) {
      std::cout << "Could not open file " << filename << std::endl;
      return -1;
    }

    /* First 4 bytes are an integer */
    fp.read((char *)&ierr,4);
    if(fp.fail() || ierr!=nrows[myrowA]*ncols[mycolA]*8) {
      std::cout << "First 8 bytes should be an integer equal to nrows*ncols*8; instead, " << ierr  << std::endl;
      return -2;
    }

    /* Read 8-byte fields */
    Atmp=new myscalar[nrows[myrowA]*ncols[mycolA]];
    for(i=0;i<nrows[myrowA]*ncols[mycolA];i++) {
      fp.read((char *)&stmp,8);
      Atmp[i]=static_cast<myscalar>(stmp);
      if(fp.fail()) {
        std::cout << "Something went wrong while reading..." << std::endl;
        if(fp.eof())
          std :: cout << "Only " << i << " instead of " << nrows[myrowA]*ncols[mycolA] << std::endl;
        return 2;
      }
    }

    /* Last 4 bytes are an integer */
    fp.read((char *)&ierr,4);
    if(fp.fail() || ierr!=nrows[myrowA]*ncols[mycolA]*8) {
      std::cout << "First 8 bytes should be an integer equal to nrows*ncols*8; instead, " << ierr  << std::endl;
      return -2;
    }
    fp.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  tend=MPI_Wtime();
  if(!myid) std::cout << "Done in: " << tend-tstart << "s" << std::endl;

  /* Initialize a context with all the processes */
  blacs_get_(&IZERO,&IZERO,&ctxtglob);
  blacs_gridinit_(&ctxtglob,"R",&IONE,&np);

  /* Initialize the BLACS grid */
  nprow=floor(sqrt((float)np));
  npcol=np/nprow;
  blacs_get_(&IZERO,&IZERO,&ctxt);
  blacs_gridinit_(&ctxt,"R",&nprow,&npcol);
  blacs_gridinfo_(&ctxt,&nprow,&npcol,&myrow,&mycol);

  /* Create A in 2D block-cyclic form by redistributing each piece */
  if(myid<nprow*npcol) {
    locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
    locc=numroc_(&n,&nb,&mycol,&IZERO,&npcol);
    dummy=std::max(1,locr);
    A=new myscalar[locr*locc];
    descinit_(descA,&n,&n,&nb,&nb,&IZERO,&IZERO,&ctxt,&dummy,&ierr);
  } else {
    descset_(descA,&n,&n,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
  }

  /* Redistribute each piece */
  if(!myid) std::cout << "Redistributing..." << std::endl;
  tstart=MPI_Wtime();
  for(i=0;i<nprowA;i++)
    for(j=0;j<npcolA;j++) {
      /* Initialize a grid that contains only the process that owns piece (i,j) */
      id=i*npcolA+j;
      blacs_get_(&IZERO,&IZERO,&ctxttmp);
      blacs_gridmap_(&ctxttmp,&id,&IONE,&IONE,&IONE);
      if(myid==id) {
        /* myid owns the piece of A to be distributed */
        descinit_(descAtmp,&nrows[i],&ncols[j],&nb,&nb,&IZERO,&IZERO,&ctxttmp,&nrows[i],&ierr);
      } else
        descset_(descAtmp,&nrows[i],&ncols[j],&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);

      pgemr2d(nrows[i],ncols[j],Atmp,IONE,IONE,descAtmp,A,rowoffset[i],coloffset[j],descA,ctxtglob);
    }
  if(myid<nprowA*npcolA)
    delete[] Atmp;
  MPI_Barrier(MPI_COMM_WORLD);
  tend=MPI_Wtime();
  if(!myid) std::cout << "Done in " << tend-tstart << "s" << std::endl;

  /* Initialize the solver and set parameters */
  StrumpackDensePackage<myscalar,myreal> sdp(MPI_COMM_WORLD);
  sdp.use_HSS=atoi(argv[1]);
  sdp.levels_HSS=atoi(argv[2]);
  sdp.min_rand_HSS=atoi(argv[3]);
  sdp.max_rand_HSS=atoi(argv[3]);
  sdp.inc_rand_HSS=0;
  sdp.lim_rand_HSS=0;
  sdp.tol_HSS=atof(argv[4]);
  sdp.steps_IR=20;
  sdp.tol_IR=1e-10;

  /* Compression */
  sdp.compress(A,descA);

  /* Accuracy checking */
  //sdp.check_compression(A,descA);

  /* Factorization */
  sdp.factor(A,descA);

  /* Process 0 reads the RHS */
  if(myid==0) {
    Btmp=new myscalar[n];
    for(j=0;j<nprowA;j++) {
      locfile="CC_"+SSTR(j)+"_"+SSTR(nrows[j]);
      filename=prefix+locfile;
      std::cout << "Process " << myid << " reading from file " << locfile << std::endl;
      fp.open(filename.c_str(),std::ios::binary);
      if(!fp.is_open()) {
        std::cout << "Could not open file " << filename << std::endl;
        return -1;
      }

      fp.read((char *)&ierr,4);
      if(fp.fail() || ierr!=nrows[j]*8) {
        std::cout << "First 4 bytes should be an integer equal to nrows*8; instead, " << ierr  << std::endl;
        return -2;
      }

      for(i=0;i<nrows[j];i++) {
        fp.read((char *)&Btmp[rowoffset[j]-1+i],8); // Binary file
        if(fp.fail()) {
          std::cout << "Something went wrong while reading..." << std::endl;
          if(fp.eof())
            std :: cout << "Only " << i << " instead of " << nrows[i] << std::endl;
          return 2;
        }
      }

      fp.read((char *)&ierr,4);
      if(fp.fail() || ierr!=nrows[j]*8) {
        std::cout << "Last 4 bytes should be an integer equal to nrows*8; instead, " << ierr  << std::endl;
        return -2;
      }

      fp.close();
    }
  }

  /* Set the RHS (random vector) and the solution space */
  if(myid<nprow*npcol) {
    locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
    locc=numroc_(&nrhs,&nb,&mycol,&IZERO,&npcol);
    B=new myscalar[locr*locc];
    dummy=std::max(1,locr);
    descinit_(descXB,&n,&nrhs,&nb,&nb,&IZERO,&IZERO,&ctxt,&dummy,&ierr);
    X=new myscalar[locr*locc];
  } else {
    descset_(descXB,&n,&nrhs,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
  }

  /* Redistribute */
  blacs_get_(&IZERO,&IZERO,&ctxttmp);
  blacs_gridmap_(&ctxttmp,&IZERO,&IONE,&IONE,&IONE);
  if(myid==0)
    /* myid owns the piece of A to be distributed */
    descinit_(descBtmp,&n,&IONE,&nb,&nb,&IZERO,&IZERO,&ctxttmp,&n,&ierr);
  else
    descset_(descBtmp,&n,&IONE,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
  pgemr2d(n,IONE,Btmp,IONE,IONE,descBtmp,B,IONE,IONE,descXB,ctxtglob);
  if(myid==0)
    delete[] Btmp;

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
