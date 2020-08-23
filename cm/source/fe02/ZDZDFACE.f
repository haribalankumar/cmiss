      SUBROUTINE ZDZDFACE(NDDL,NDLT,nf,NPF,nx,WD,WDL,XID,XIDL,ZD,ZDL,
     '  ERROR,*)

C#### Subroutine: ZDZDFACE
C###  Description:
C###    ZDZDFACE puts data parameters XID,ZD,WD into element(face)
C###    arrays XIDL,ZDL,WDL. NOTE: In face fitting of volume elements,
C###    faces are treated as elements.
C**** Written by Kumar Mithraratne, Aug. 2002.


      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NDDL(NEM,NDEM),NDLT,nf,NPF(9,NFM),nx
      REAL*8 WD(NJM,NDM),WDL(NHM,NDEM),XID(NIM,NDM),XIDL(NIM,NDEM),
     '  ZD(NJM,NDM),ZDL(NHM,NDEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nd,ndf,nh,nhj,nhx,ni1,ni2,nj,njj

      CALL ENTERS('ZDZDFACE',*9999)

      ni1=NPF(1,nf) ! 1st xi direction of the face
      ni2=NPF(3,nf) ! 2nd xi direction of the face

      DO ndf=1,NDLT
        nd=NDDL(nf,ndf)
        XIDL(1,ndf)=XID(ni1,nd)
        XIDL(2,ndf)=XID(ni2,nd)
        DO njj=1,NUM_FIT(0)
          DO nhj=1,NUM_FIT(njj)
            nj=NJ_FIT(nhj,njj)
            nhx=NLH_FIT(nhj,3,njj)
            nh=NH_LOC(nhx,nx)
            ZDL(nh,ndf)=ZD(nj,nd)
            WDL(nh,ndf)=WD(nj,nd)
          ENDDO !nhj
        ENDDO !njj
      ENDDO !nde

      CALL EXITS('ZDZDFACE')
      RETURN
 9999 CALL ERRORS('ZDZDFACE',ERROR)
      CALL EXITS('ZDZDFACE')
      RETURN 1
      END



