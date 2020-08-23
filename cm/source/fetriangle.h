/* Functions defined in fetriangle.c referenced from Fortran */
#ifdef VMS
#define ALE_initialise ALE_INITIALISE
#define ALE_flip ALE_FLIP
#define ALE_destroy ALE_DESTROY
#endif
#if defined(unix) || defined(_AIX) || defined (WIN32)
#  define ALE_initialise ale_initialise_
#  define ALE_flip ale_flip_
#  define ALE_destroy ale_destroy_
#endif

extern void ALE_initialise(
  int *NAM,
  int *nb,
  int *NBFM,
  int *NCM,
  int *NELIST,
  int *NEM,
  int *NEELEM,
  int  NEIM,
  int *NHM,
  int *NIM,
  int *NJM,
  int *NKM,
  int *NNM,
  int *NPLIST,
  int *NPM,
  int *NPNE,
  int *NPT,
  int *NVJE,
  int *NVM,
  int *NXI,
  double *XP,
  double *ZA,
  double *zero_tol,
  int *Err_Flag,
  char *err_StringC
  );

extern void ALE_flip(
  int *NAM,
  int *nb,
  int *NBFM,
  int *NCM,
  int *NELIST,
  int *NEM,
  int *NEELEM,
  int  NEIM,
  int *NHM,
  int *NIM,
  int *NJM,
  int *NKM,
  int *NNM,
  int *NPLIST,
  int *NPM,
  int *NPNE,
  int *NPT,
  int *NVJE,
  int *NVM,
  int *NXI,
  double *XP,
  double *ZA,
  int *Err_Flag,
  char *err_StringC
  );

extern void ALE_destroy(
  );
 
