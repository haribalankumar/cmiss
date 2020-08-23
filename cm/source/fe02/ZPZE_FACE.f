      SUBROUTINE ZPZE_FACE(nf,NBHF,nc,NKHF,NPNF,NVHF,nx,SF,ZE,ZP,
     '  ERROR,*)

C#### Subroutine: ZPZE_FACE
C###  Description:
C###    ZPZE_FACE transfers global node parameters ZP(nk,nv,nh,np,nc)
C###    to element face parameters for face fitting

C CS 13/1/2002 hacked up version of ZPZE for face fitting
C PM 14Aug02 : Removed all unnecessary parameters.
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER nf,NBHF(NHM,NCM,NFM),nc,NKHF(NKM,NNM,NHM),
     '  NPNF(NNM,NBFM),NVHF(NNM,NBFM,NHM),nx
      REAL*8 SF(NSM,NBFM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nh,nhx,nk,nk1,nn,np,ns,nv

      CALL ENTERS('ZPZE_FACE',*9999)

      DO nhx=1,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        nb=NBHF(nh,nc,nf)
C news VJ 3Feb2005. Do not calculate ZE if basis function not available 
C for dependent variable on a face eg hydrostatic pressure for mechanics	
        IF(nb.NE.0) THEN
          ns=0
          DO nn=1,NNT(nb)
            np=NPNF(nn,nb)
            nv=NVHF(nn,nb,nh)
            DO nk1=1,NKT(nn,nb)
              nk=NKHF(nk1,nn,nh)
              ns=ns+1
              ZE(ns,nhx)=ZP(nk,nv,nh,np,nc)*SF(ns,nb)
            ENDDO !nk1
          ENDDO !nn
        ENDIF
      ENDDO !nhx
C news VJ 6Feb2005: Added DOP for ZE matrix.
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          nb=NBHF(nh,nc,nf)
          IF(nb.NE.0) THEN
            WRITE(OP_STRING,
     '        '('' ZE(ns,'',I2,''): '',8D11.3,/(12X,8D11.3))')
     '        nhx,(ZE(ns,nhx),ns=1,NST(NBHF(nh,nc,nf))
     &        +NAT(NBHF(nh,nc,nf)))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !nhx
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZPZE_FACE')
      RETURN
 9999 CALL ERRORS('ZPZE_FACE',ERROR)
      CALL EXITS('ZPZE_FACE')
      RETURN 1
      END


