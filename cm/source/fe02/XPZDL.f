      SUBROUTINE XPZDL(IP,NBH,NBJ,NPLIST2,nx,WDL,XIP,XIDL,XP,ZDL,
     '  ERROR,*)

C#### Subroutine: XPZDL
C###  Description:
C###    XPZDL puts data parameters XIP,XD,WD into element arrays XIDL,
C###    ZDL,WDL. This routine is an adaption of ZDZDL and is designed
C###    to be used inconjunction with a node group.
C###    If IP=0 then the geometry information is placed in
C###    the element arrays. If IP=1 then the fitting information is
C###    placed in the element arrays. NOTE: The element arrays are NHM
C###    based and hence if IP=0 then NHM must be >= NJM.

      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IP,NBH(NHM),NBJ(NJM),nx,NPLIST2(0:NPM)
      REAL*8 WDL(NHM,NDEM),XIP(NIM,NPM),XIDL(NIM,NDEM),
     '  XP(NKM,NVM,NJM,NPM),ZDL(NHM,NDEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,npe,nh,nhj,nhx,ni,nj,njj,np

      CALL ENTERS('XPZDL',*9999)

C CPB 5/10/95 Changing xDL arrays to be nh based.

      IF(IP.EQ.0) THEN !Geometry information required
        CALL ASSERT(NHM.GE.NJM,'>>NHM must be >= NJM',ERROR,*9999)
        DO npe=1,nplist2(0)
          np=nplist2(npe)
          DO nj=1,NJT
            nb=NBJ(nj)
            DO ni=1,NIT(nb)
              XIDL(ni,npe)=XIP(ni,np)
            ENDDO !ni
            ZDL(nj,npe)=XP(1,1,nj,np)
            WDL(nj,npe)=1 !hack for my problems - all weights 1
          ENDDO !nj
        ENDDO !npe
      ELSE IF(IP.EQ.1) THEN !Fitting information required
        DO npe=1,nplist2(0)
          np=nplist2(npe)
          DO njj=1,NUM_FIT(0)
            DO nhj=1,NUM_FIT(njj)
              nj=NJ_FIT(nhj,njj) !is data variable nj for fit var njj.
              nhx=NLH_FIT(nhj,3,njj) !is nh where fit var njj is stored.
              nh=NH_LOC(nhx,nx)
              nb=NBH(nh)
              DO ni=1,NIT(nb)
                XIDL(ni,npe)=XIP(ni,np)
              ENDDO !ni
              ZDL(nh,npe)=XP(1,1,nj,np)
              WDL(nh,npe)=1 !hack for my problems - all weights 1
            ENDDO !nhj
          ENDDO !njj
        ENDDO !npe
      ELSE
        ERROR='>>Invalid IP number'
        GOTO 9999
      ENDIF

      CALL EXITS('XPZDL')
      RETURN
 9999 CALL ERRORS('XPZDL',ERROR)
      CALL EXITS('XPZDL')
      RETURN 1
      END


