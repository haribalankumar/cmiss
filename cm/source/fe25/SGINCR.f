      SUBROUTINE SGINCR(INDEX,ISEG,ISINCR,iw,NHP,NKH,NPNODE,nx,
     '  CSEG,FIX,XP,ZP,ERROR,*)

C#### Subroutine: SGINCR
C###  Description:
C###    SGINCR creates increment segment ISINCR.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'trans00.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISINCR,iw,NHP(NPM),NKH(NHM,NPM,NCM),
     '  NPNODE(0:NP_R_M,0:NRM),nx
      REAL*8 XP(NKM,NVM,NJM,NPM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER INDEX_OLD,nc,nh,nhx,nj,nk,nonode,np,nr,nv,ny
      REAL*8 PTS(3,2),DZ,X1(3),X2(3),Z1(3),Z2(3),Z3(3),Z4(3)

      CALL ENTERS('SGINCR',*9999)
      nc=1 !Temporary MPN 12-Nov-94
      nv=1 !temp
      CALL OPEN_SEGMENT(ISINCR,ISEG,iw,'INCR',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)

      ny=0
      DO nr=1,NRT
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            DO nk=1,NKH(nh,np,nc)
              ny=ny+1
              IF(FIX(ny,3)) THEN
                X1(nh)=XP(1,nv,nh,np)
                DO nj=1,NJT
                  IF(nj.NE.nh) X1(nj)=ZP(1,nv,nj,np,nc)
                  X2(nj)=ZP(1,nv,nj,np,nc)
                ENDDO
                CALL XZ(ITYP10(1),X1,Z1)
                CALL XZ(ITYP10(1),X2,Z2)
                IF(nh.EQ.1) THEN
                  DZ=(Z2(1)-Z1(1))/4.0D0
                  Z3(1)=Z2(1)-DZ
                  Z4(1)=Z2(1)-DZ
                  Z3(2)=Z2(2)+DZ
                  Z4(2)=Z2(2)-DZ
                  Z3(3)=Z2(3)+DZ
                  Z4(3)=Z2(3)-DZ
                ELSE IF(nh.EQ.2) THEN
                  DZ=(Z2(2)-Z1(2))/4.0D0
                  Z3(1)=Z2(1)+DZ
                  Z4(1)=Z2(1)-DZ
                  Z3(2)=Z2(2)-DZ
                  Z4(2)=Z2(2)-DZ
                  Z3(3)=Z2(3)+DZ
                  Z4(3)=Z2(3)-DZ
                ELSE IF(nh.EQ.3) THEN
                  DZ=(Z2(3)-Z1(3))/4.0D0
                  Z3(1)=Z2(1)+DZ
                  Z4(1)=Z2(1)-DZ
                  Z3(2)=Z2(2)+DZ
                  Z4(2)=Z2(2)-DZ
                  Z3(3)=Z2(3)-DZ
                  Z4(3)=Z2(3)-DZ
                ENDIF
                IF(iw.EQ.3) THEN
                  CALL ZZ(Z1,Z1,TRANS)
                  CALL ZZ(Z2,Z2,TRANS)
                  CALL ZZ(Z3,Z3,TRANS)
                  CALL ZZ(Z4,Z4,TRANS)
                ENDIF
                PTS(1,1)=Z1(1)
                PTS(1,2)=Z2(1)
                PTS(2,1)=Z1(2)
                PTS(2,2)=Z2(2)
                PTS(3,1)=Z1(3)
                PTS(3,2)=Z2(3)
                CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
                PTS(1,1)=Z2(1)
                PTS(1,2)=Z3(1)
                PTS(2,1)=Z2(2)
                PTS(2,2)=Z3(2)
                PTS(3,1)=Z2(3)
                PTS(3,2)=Z3(3)
                CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
                PTS(1,1)=Z2(1)
                PTS(1,2)=Z4(1)
                PTS(2,1)=Z2(2)
                PTS(2,2)=Z4(2)
                PTS(3,1)=Z2(3)
                PTS(3,2)=Z4(3)
                CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL CLOSE_SEGMENT(ISINCR,iw,ERROR,*9999)

      CALL EXITS('SGINCR')
      RETURN
 9999 CALL ERRORS('SGINCR',ERROR)
      CALL EXITS('SGINCR')
      RETURN 1
      END


