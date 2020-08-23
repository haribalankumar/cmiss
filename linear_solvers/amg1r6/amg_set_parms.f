
      SUBROUTINE AMG_SET_PARMS(AMG_CYCLE,AMG_NDIMS,AMG_IPARMS,
     &  AMG_RPARMS,MAXIT,N,NX,NZA,TOL,ERROR,*)

C#### Subroutine: AMG_SET_PARMS
C###  Description:
C###    AMG_SET_PARMS sets various parameter values by defining arrays
C###    defined in 'amgptr.cmn'.  These arrays are necessary for calling
C###    AMG.  To understand some of the choices for parameters read 
C###    through the opening comments in the amg1r6.f file.
C###  Written by Travis Austin 02/12/04

      IMPLICIT NONE
C     Parameter List
      INTEGER AMG_NDIMS(7),AMG_IPARMS(14)
      REAL*8 AMG_RPARMS(4)
      INTEGER MAXIT,N,NX,NZA,AMG_CYCLE
      REAL*8 TOL
      CHARACTER ERROR*(*)
C     Local Variables
      INTEGER NDA, NDU, NDF, NDIA, NDJA, NDIG, MATRIX
      INTEGER ISWTCH,IOUT,IPRINT
      INTEGER LEVELX,IFIRST,IFMG,NCYC,EPS,
     &        MADAPT,NRD,NSOLCO,NRU,NWT,NTR
      REAL*8  ECG1,ECG2,EWT2

      CALL ENTERS('AMG_SET_PARMS',*9999)

      NDA  = 4*NZA+5*N
      NDIA = int(2.8*N)
      NDJA = 4*NZA+5*N
      NDU  = int(2.8*N)
      NDF  = int(2.8*N)
      NDIG = int(5.4*N)

      AMG_NDIMS(1) = NDA
      AMG_NDIMS(2) = NDIA
      AMG_NDIMS(3) = NDJA
      AMG_NDIMS(4) = NDU
      AMG_NDIMS(5) = NDF
      AMG_NDIMS(6) = NDIG
      AMG_NDIMS(7) = N

C NOTE: See opening comments in amg1r6.f

      MATRIX  = 12

C  *  MATRIX   -   INTEGER VALUE CONTAINING INFO ABOUT THE MATRIX L.
C
C                  1ST DIGIT OF MATRIX  --  ISYM:
C                    =1: L IS SYMMETRIC;
C                    =2: L IS NOT SYMMETRIC.
C
C                  2ND DIGIT OF MATRIX  --  IROW0:
C                    =1: L HAS ROWSUM ZERO;
C                    =2: L DOES NOT HAVE ROWSUM ZERO.

      IOUT    = 00
      IPRINT  = 12030
      ISWTCH  = 5      ! SETUP ONLY

C     Set Cycle Types and Number of Cycles
C
C       CYCL = 1 : V-CYCLE
C            = 4 : W-CYCLE
      
      IF(MAXIT.LT.10) THEN
         NCYC = AMG_CYCLE*1000  + 40  + MAXIT
      ELSEIF(MAXIT.LT.100) THEN
         NCYC = AMG_CYCLE*10000  + 400 + MAXIT
      ELSEIF(MAXIT.LT.1000) THEN
         NCYC = AMG_CYCLE*100000  + 4000 + MAXIT
      ENDIF

      LEVELX  = 10
      IFIRST  = 13
      IFMG    = 0
      EPS     = TOL
      MADAPT  = 27 
      NRD     = 1131
      NSOLCO  = 2
      NRU     = 1131
      
      ECG1   = 0.01
      ECG2   = 0.25
      EWT2   = 0.35
      NWT    = 2
      NTR    = 0

      AMG_IPARMS(1)  = MATRIX
      AMG_IPARMS(2)  = ISWTCH  
      AMG_IPARMS(3)  = IOUT
      AMG_IPARMS(4)  = IPRINT
      AMG_IPARMS(5)  = LEVELX
      AMG_IPARMS(6)  = IFIRST
      AMG_IPARMS(7)  = IFMG
      AMG_IPARMS(8)  = NCYC
      AMG_IPARMS(9)  = MADAPT
      AMG_IPARMS(10) = NRD
      AMG_IPARMS(11) = NSOLCO
      AMG_IPARMS(12) = NRU
      AMG_IPARMS(13) = NWT
      AMG_IPARMS(14) = NTR

      AMG_RPARMS(1) = EPS
      AMG_RPARMS(2) = ECG1
      AMG_RPARMS(3) = ECG2
      AMG_RPARMS(4) = EWT2     

      CALL EXITS('AMG_SET_PARMS')
      RETURN

 9999 CALL ERRORS('AMG_SET_PARMS',ERROR)
      RETURN 1
      END
