      SUBROUTINE UPFGNODES(FROM,TO,FG,FGNK,FGNV,CONST,CP,INDICES,nc,NKH,
     '  NKJ,NPLIST_LOCAL,NUMVALUES,NVHP,NVJP,NYNP,POS,SCALE,XP,YP,ERROR,
     &  *)

C#### Subroutine: UPFGNODE
C###  Description:
C###    This subroutine copies values and derivatives at node points
C###    between FG and XP/YP/CP


      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NUMVALUES,FGNK(2,0:NUMVALUES),FGNV(2,0:NUMVALUES),
     '  INDICES(10,2),nc,NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NPLIST_LOCAL(0:10000),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJP(NJM,NPM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),POS
      REAL*8  CONST,CP(NMM,NPM,NXM),FG(NKM,NVM,NUMVALUES),SCALE,
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
      CHARACTER FROM*(*),TO*(*),ERROR*(*)
!      LOGICAL
!     Local Variables
      INTEGER nh,niy,nj,nk,nk1,nk2,nkc,nm,nonode,np,nr,nv,nx,ny
!      REAL*8
!      LOGICAL
!      CHARACTER

      CALL ENTERS('UPFGNODES',*9999)
      nm=INDICES(1,POS)
      niy=INDICES(2,POS)
      nh=INDICES(3,POS)
      nj=INDICES(4,POS)
      nx=INDICES(5,POS)
      nr=INDICES(6,POS)

      IF(TO(1:2).EQ.'FG')THEN
        IF(FROM(1:8).EQ.'CONSTANT')THEN
          FGNK(POS,0)=NPLIST_LOCAL(0)
          DO nonode=1,NPLIST_LOCAL(0)
            np=NPLIST_LOCAL(nonode)
            FGNK(POS,nonode)=NKM
            FGNV(POS,nonode)=FGNV(1,nonode)
            DO nv=1,FGNV(POS,nonode)
              FG(1,nv,nonode)=SCALE*CONST
              nkc=0
              DO nk=2,FGNK(2,nonode)
                nkc=nkc+1
                FG(nk,nv,nonode)=0.0D0
              ENDDO
            ENDDO
          ENDDO
        ELSEIF(FROM(1:2).EQ.'CP')THEN
          FGNK(POS,0)=NPLIST_LOCAL(0)
          DO nonode=1,NPLIST_LOCAL(0)
            np=NPLIST_LOCAL(nonode)
            FGNK(POS,nonode)=1
            FGNV(POS,nonode)=1
            FG(1,1,nonode)=SCALE*CP(nm,np,nx)
          ENDDO
        ELSEIF(FROM(1:2).EQ.'XP')THEN
          FGNK(POS,0)=NPLIST_LOCAL(0)
          DO nonode=1,NPLIST_LOCAL(0)
            np=NPLIST_LOCAL(nonode)
            IF(INDICES(7,POS).EQ.0)THEN
              nk1=1
              nk2=NKJ(nj,np)
              FGNK(POS,nonode)=NKJ(nj,np)
            ELSE
              nk1=INDICES(7,POS)
              nk2=INDICES(7,POS)
              FGNK(POS,nonode)=1
            ENDIF
            FGNV(POS,nonode)=NVJP(nj,np)
            DO nv=1,NVJP(nj,np)
              nkc=0
              DO nk=nk1,nk2
                nkc=nkc+1
                FG(nkc,nv,nonode)=SCALE*XP(nk,nv,nj,np)
              ENDDO
            ENDDO
          ENDDO
        ELSEIF(FROM(1:2).EQ.'YP')THEN
          FGNK(POS,0)=NPLIST_LOCAL(0)
          DO nonode=1,NPLIST_LOCAL(0)
            np=NPLIST_LOCAL(nonode)
            IF(INDICES(7,POS).EQ.0)THEN
              nk1=1
              nk2=NKH(nh,np,1,nr)
              FGNK(POS,nonode)=NKH(nh,np,1,nr)
            ELSE
              nk1=INDICES(7,POS)
              nk2=INDICES(7,POS)
              FGNK(POS,nonode)=1
            ENDIF
            FGNV(POS,nonode)=NVHP(nh,np,1,nr)
            DO nv=1,NVHP(nh,np,1,nr)
              nkc=0
              DO nk=nk1,nk2
                nkc=nkc+1
                ny=NYNP(nk,nv,nh,np,0,nc,nr)
                FG(nkc,nv,nonode)=SCALE*YP(ny,niy,nx)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSEIF(FROM(1:2).EQ.'FG')THEN
        IF(TO(1:2).EQ.'CP')THEN
          DO nonode=1,NPLIST_LOCAL(0)
            np=NPLIST_LOCAL(nonode)
            FGNK(POS,nonode)=1
            DO nv=1,FGNV(POS,nonode)
              CP(nm,np,nx)=FG(1,nv,nonode)
            ENDDO
          ENDDO
        ELSEIF(TO(1:2).EQ.'XP')THEN
          DO nonode=1,NPLIST_LOCAL(0)
            np=NPLIST_LOCAL(nonode)
            IF(INDICES(7,POS).EQ.0)THEN
              nk1=1
              nk2=NKJ(nj,np)
              FGNK(POS,nonode)=NKJ(nj,np)
            ELSE
              nk1=INDICES(7,POS)
              nk2=INDICES(7,POS)
              FGNK(POS,nonode)=1
            ENDIF
            nkc=0
            DO nk=nk1,nk2
              nkc=nkc+1
              DO nv=1,NVJP(nj,np)
                XP(nk,nv,nj,np)=FG(nkc,nv,nonode)
              ENDDO
            ENDDO
          ENDDO
        ELSEIF(TO(1:2).EQ.'YP')THEN
          DO nonode=1,NPLIST_LOCAL(0)
            np=NPLIST_LOCAL(nonode)
            IF(INDICES(7,POS).EQ.0)THEN
              nk1=1
              nk2=NKH(nh,np,1,nr)
              FGNK(POS,nonode)=NKH(nh,np,1,nr)
            ELSE
              nk1=INDICES(7,POS)
              nk2=INDICES(7,POS)
              FGNK(POS,nonode)=1
            ENDIF
            DO nv=1,NVHP(nh,np,1,nr)
              nkc=0
              DO nk=nk1,nk2
                nkc=nkc+1
                ny=NYNP(nk,nv,nh,np,0,nc,nr)
                YP(ny,niy,nx)=FG(nkc,nv,nonode)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('UPFGNODES')
      RETURN
 9999 CALL ERRORS('UPFGNODES',ERROR)
      CALL EXITS('UPFGNODES')
      RETURN 1
      END


