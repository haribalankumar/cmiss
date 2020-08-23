#include <stdio.h>
#include <string.h>
#include "message.h"
#include "fetriangle.h"

extern void ALE_initialise(
  int *NAM,
  int *nb,
  int *NBFM,
  int *NCM,
  int *NELIST,
  int *NEM,
  int *NEELEM,
  int *NHM,
  int *NIM,
  int *NJM,
  int *NKM,
  int *NNM,
  int *NPLIST,
  int *NPM,
  int *NPNE,
  int *NPNODE,
  int *NVJE,
  int *NVM,
  int *NXI,
  double *XP,
  double *ZA,
  double *zero_tol,
  int *err_Flag,
  char *err_StringC
  )
{
  display_message(ERROR_MESSAGE,">>Link with fetriangle.c>ALE_initialise");
  *err_StringC = NULL;
  *err_Flag = 1;
}

extern void ALE_flip(
  int *NAM,
  int *nb,
  int *NBFM,
  int *NCM,
  int *NELIST,
  int *NEM,
  int *NEELEM,
  int *NHM,
  int *NIM,
  int *NJM,
  int *NKM,
  int *NNM,
  int *NPLIST,
  int *NPM,
  int *NPNE,
  int *NPNODE,
  int *NVJE,
  int *NVM,
  int *NXI,
  double *XP,
  double *ZA,
  int *err_Flag,
  char *err_StringC
  )
{
  display_message(ERROR_MESSAGE,">>Link with fetriangle.c>ALE_flip");
  *err_StringC = NULL;
  *err_Flag = 1;
}

extern void ALE_destroy()
{
}
