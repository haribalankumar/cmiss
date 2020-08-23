      SUBROUTINE OPTIMISER_ALLOC(INTWORK_PTR,REALWORK_PTR,ERROR,*)

C#### Subroutine: OPTIMISER_ALLOC
C###  Description:
C###    Allocate memory for the MinSQP and MINOS optimisers.

      IMPLICIT NONE
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'parameters.inc'

!     Parameter List
      INTEGER*4 INTWORK_PTR,REALWORK_PTR
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER USE_OPT
C     Lengths of workspace used by both optimisers
      INTEGER L_PMIN,L_PMAX,L_PAOPTI,L_PBOPTI,L_ROPTI
C     Lengths of the workspace for the MinSQP optimiser
      INTEGER L_ISTATE,L_NEEDCON,L_CONJAC,L_CONTR,L_R,L_RESID,L_RESJAC,
     '  L_C,L_F,L_CJAC,L_FJAC,L_XC,L_INPSOL,L_RNPSOL
C     Lengths of the workspace for the MINOS optimiser
      INTEGER L_KAM,L_HAM,L_HSM,L_AM,L_BLM,L_BUM,L_XNM,L_PIM,L_RCM,
     '  L_RESIDM,L_FGRADM,L_CM,L_CJACM,L_IMINOS,L_RMINOS
C     INTEGER L_IWORK,L_WORK,L_ZM


      CALL ENTERS('OPTIMISER_ALLOC',*9999)

C CPB 12/10/93 - Redoing the integer work array size calculation. NAG
C suggests the following size
C     L_IWORK=3*NOPM+NLCM+2*NCOM+NREM*(NOPM+3))
C but as there are no linear contraints the following calc. can be used.
C     L_IWORK=(3*NOPM+2*NCOM)*USE_NPSOL

C SEN 22/01/03 - changed L_IWORK and L_WORK to 1, to reflect their new
C use for parsing stuff through to the FUNCT2() routine.
C     L_IWORK=1*USE_NPSOL
C     L_WORK=1*USE_NPSOL

C CPB 20/7/93 Setup optimisation arrays by subdividing the 1D WORK
C     arrays
C SEN 13/2/03 Cleaned up and shifted to a separate routine.

C     Lengths used by both optimisers
      USE_OPT   = MAX(USE_NPSOL,USE_MINOS)
      L_PAOPTI  = USE_OPT* NOPM
      L_PBOPTI  = USE_OPT* NOPM
      L_PMIN    = USE_OPT* (NOPM+NCOM)
      L_PMAX    = USE_OPT* (NOPM+NCOM)
      L_ROPTI   = L_PAOPTI+L_PBOPTI+L_PMIN+L_PMAX

C     Lengths of the NPSOL/MinSQP arrays
C     Integer arrays
C     NCOM      = 1 ! No nonlinear constraints anymore
      L_ISTATE  = USE_NPSOL* (NOPM+NCOM)
      L_IWORK   = USE_NPSOL* 2
      L_NEEDCON = USE_NPSOL* NCOM
      L_INPSOL  = L_ISTATE+L_IWORK+L_NEEDCON

C     Real Arrays
      L_CONJAC  = USE_NPSOL* NOPM*NCOM
      L_CONTR   = USE_NPSOL* NCOM
      L_R       = USE_NPSOL* NOPM*NOPM
      L_RESID   = USE_NPSOL* NREM
      L_RESJAC  = USE_NPSOL* NREM*NOPM
      L_WORK    = USE_NPSOL* 1
      L_C       = USE_NPSOL* NCOM
      L_F       = USE_NPSOL* NREM
      L_CJAC    = USE_NPSOL* NCOM*NOPM
      L_FJAC    = USE_NPSOL* NREM*NOPM
      L_XC      = USE_NPSOL* NOPM
      L_RNPSOL  = L_CONJAC+L_CONTR+L_R+L_RESID+L_RESJAC+L_WORK+L_C
     '           +L_F+L_CJAC+L_FJAC+L_XC

C     Lengths of the MINOS arrays
C     Integer arrays
      L_KAM    = USE_MINOS* (NOPM+1)
      L_HAM    = USE_MINOS* NZ_MINOSM*2
      L_HSM    = USE_MINOS* L_ZM
      L_IMINOS = L_KAM+L_HAM+L_HSM
C     Real Arrays
      L_AM     = USE_MINOS* NZ_MINOSM
      L_BLM    = USE_MINOS* (NOPM+NCOM+NLCM)
      L_BUM    = USE_MINOS* (NOPM+NCOM+NLCM)
      L_XNM    = USE_MINOS* (NOPM+NCOM+NLCM)
      L_PIM    = USE_MINOS* (NCOM+NLCM)
      L_RCM    = USE_MINOS* (NOPM+NCOM+NLCM)
      L_RESIDM = USE_MINOS* NREM
      L_FGRADM = USE_MINOS* NOPM
      L_CM     = USE_MINOS* NCOM
      L_CJACM  = USE_MINOS* NZ_MINOSM
      L_ZM     = USE_MINOS* (NCOM+NLCM)*10000
      L_RMINOS = L_AM+L_BLM+L_BUM+L_XNM+L_PIM+L_RCM+L_RESIDM+L_FGRADM
     '          +L_CM+L_CJACM+L_ZM

C     Total lengths of the required workspace
      NIWM = L_INPSOL + L_IMINOS
      NRWM = L_RNPSOL + L_RMINOS + L_ROPTI

C     Integer arrays
C     NPSOL arrays
      OS_ISTATE = 1
      OS_IWORK  = OS_ISTATE + L_ISTATE
      OS_NEEDCON= OS_IWORK  + L_IWORK
C     MINOS arrays
      OS_KAM    = 1         + L_INPSOL
      OS_HAM    = OS_KAM    + L_KAM
      OS_HSM    = OS_HAM    + L_HAM

C     Real arrays
C     Common arrays
      OS_PAOPTI = 1
      OS_PBOPTI = OS_PAOPTI + L_PAOPTI
      OS_PMIN   = OS_PBOPTI + L_PBOPTI
      OS_PMAX   = OS_PMIN   + L_PMIN
C     NPSOL arrays
      OS_CONJAC = 1         + L_ROPTI
      OS_CONTR  = OS_CONJAC + L_CONJAC
      OS_R      = OS_CONTR  + L_CONTR
      OS_RESID  = OS_R      + L_R
      OS_RESJAC = OS_RESID  + L_RESID
      OS_WORK   = OS_RESJAC + L_RESJAC
      OS_C      = OS_WORK   + L_WORK
      OS_F      = OS_C      + L_C
      OS_CJAC   = OS_F      + L_F
      OS_FJAC   = OS_CJAC   + L_CJAC
      OS_XC     = OS_FJAC   + L_FJAC
C     MINOS arrays
      OS_AM     = 1         + L_INPSOL + L_ROPTI
      OS_BLM    = OS_AM     + L_AM
      OS_BUM    = OS_BLM    + L_BLM
      OS_XNM    = OS_BUM    + L_BUM
      OS_PIM    = OS_XNM    + L_XNM
      OS_RCM    = OS_PIM    + L_PIM
      OS_RESIDM = OS_RCM    + L_RCM
      OS_FGRADM = OS_RESIDM + L_RESIDM
      OS_CM     = OS_FGRADM + L_FGRADM
      OS_CJACM  = OS_CM     + L_CM
      OS_ZM     = OS_CJACM  + L_CJACM

      CALL ALLOCATE_MEMORY(NIWM,1,INTTYPE,INTWORK_PTR,MEM_INIT,
     '  ERROR,*9999)
      CALL ALLOCATE_MEMORY(NRWM,1,DPTYPE,REALWORK_PTR,MEM_INIT,
     '  ERROR,*9999)

      CALL EXITS('OPTIMISER_ALLOC')
      RETURN
 9999 CALL ERRORS('OPTIMISER_ALLOC',ERROR)
      CALL EXITS('OPTIMISER_ALLOC')
      RETURN 1
      END
