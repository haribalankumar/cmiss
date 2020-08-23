      INTEGER FUNCTION GETNYR_NONR(nc,NPNY,nrc,nrc1,ny,NYNE,NYNP)

C#### Function: GETNYR_NONR
C###  Type: INTEGER
C###  Description:
C###    GETNYR_NONR returns the corresponding variable nyr with a nrc and
C###    nc for a variable ny from nrc1, eg. the equivalent flux
C###    variable. Note that this function is the same as GETNYR,
C###    without regions.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nc,NPNY(0:6,NYM,0:NRCM),nrc,nrc1,ny,
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM)
!     Local Variables
      INTEGER na,ne,nh,nk,np,nv
      CHARACTER ERROR*100

      GETNYR_NONR=0
      IF(NPNY(0,ny,nrc1).EQ.1) THEN !variable is nodally based
        nk=NPNY(1,ny,nrc1)
        nv=NPNY(2,ny,nrc1)
        nh=NPNY(3,ny,nrc1)
        np=NPNY(4,ny,nrc1)
        GETNYR_NONR=NYNP(nk,nv,nh,np,nrc,nc)
      ELSE IF(NPNY(0,ny,nrc1).EQ.2) THEN !variable is element based
        na=NPNY(1,ny,nrc1)
        nh=NPNY(2,ny,nrc1)
        ne=NPNY(4,ny,nrc1)
        GETNYR_NONR=NYNE(na,nh,nrc,nc,ne)
      ELSE
        WRITE(OP_STRING,'('' >>Error in GETNYR_NONR: NPNY(0,'',I5,'
     '    //''','',I1,'') is out of range.'')') ny,nrc1
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

      IF(GETNYR_NONR.EQ.0) THEN
        WRITE(OP_STRING,'('' >>Error in GETNYR: Could not find nyr'
     '    //' for ny ='',I5)') ny
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

 9999 RETURN
      END


