      SUBROUTINE ZDZDL(IP,NBH,NBJ,NDDL,NDLT,ne,nx,WD,WDL,XID,XIDL,ZD,
     '  ZDL,ERROR,*)

C#### Subroutine: ZDZDL
C###  Description:
C###    ZDZDL puts data parameters XID,ZD,WD into element arrays XIDL,
C###    ZDL,WDL. If IP=0 then the geometry information is placed in
C###    the element arrays. If IP=1 then the fitting information is
C###    placed in the element arrays. NOTE: The element arrays are NHM
C###    based and hence if IP=0 then NHM must be >= NJM.

      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IP,NBH(NHM),NBJ(NJM),NDDL(NEM,NDEM),NDLT,ne,nx
      REAL*8 WD(NJM,NDM),WDL(NHM,NDEM),XID(NIM,NDM),XIDL(NIM,NDEM),
     '  ZD(NJM,NDM),ZDL(NHM,NDEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nd,nde,nh,nhj,nhx,ni,nj,njj

      CALL ENTERS('ZDZDL',*9999)

C CPB 5/10/95 Changing xDL arrays to be nh based.

      IF(IP.EQ.0) THEN !Geometry information required
        CALL ASSERT(NHM.GE.NJM,'>>NHM must be >= NJM',ERROR,*9999)
        DO nde=1,NDLT
          nd=NDDL(ne,nde)
          DO nj=1,NJT
            nb=NBJ(nj)
            DO ni=1,NIT(nb)
              XIDL(ni,nde)=XID(ni,nd)
            ENDDO !ni
            ZDL(nj,nde)=ZD(nj,nd)
            WDL(nj,nde)=WD(nj,nd)
          ENDDO !nj
        ENDDO !nde
      ELSE IF(IP.EQ.1) THEN !Fitting information required
        DO nde=1,NDLT
          nd=NDDL(ne,nde)
          DO njj=1,NUM_FIT(0)
            DO nhj=1,NUM_FIT(njj)
              nj=NJ_FIT(nhj,njj) !is data variable nj for fit var njj.
              nhx=NLH_FIT(nhj,3,njj) !is nh where fit var njj is stored.
              nh=NH_LOC(nhx,nx)
              nb=NBH(nh)
              DO ni=1,NIT(nb)
                XIDL(ni,nde)=XID(ni,nd)
              ENDDO !ni
              ZDL(nh,nde)=ZD(nj,nd)
              WDL(nh,nde)=WD(nj,nd)
            ENDDO !nhj
          ENDDO !njj
        ENDDO !nde
      ELSE
        ERROR='>>Invalid IP number'
        GOTO 9999
      ENDIF

      CALL EXITS('ZDZDL')
      RETURN
 9999 CALL ERRORS('ZDZDL',ERROR)
      CALL EXITS('ZDZDL')
      RETURN 1
      END



