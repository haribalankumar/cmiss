      SUBROUTINE SGREAC(INDEX,ISEG,ISREAC,iw,IY,nc,NHP,NKH,NPL,
     '  NPNODE,nr,nx,NYNR,CSEG,FIX,XP,YP,ZP,ERROR,*)

C#### Subroutine: SGREAC
C###  Description:
C###    SGREAC creates increment segment ISREAC.
C**** X1(nj) are the curvilinear coords of the arrow head
C**** Z1(nj) are the rect. cart. coords of the arrow head
C**** X2(nj) are the curvilinear coords of the arrow tail
C**** Z2(nj) are the rect. cart. coords of the arrow tail

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'reac00.cmn'
      INCLUDE 'trans00.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISREAC,iw,IY,nc,NHP(NPM),NKH(NHM,NPM,NCM),
     '  NPL(5,0:3,NLM),NPNODE(0:NP_R_M,0:NRM),nr,nx,
     '  NYNR(0:NY_R_M,0:NRCM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER INDEX_OLD,n1,n2,nh,nhx,nj,nk,nl,nonode,no_nynr,np,nv,ny
      REAL*8 DET,DXDXI(3,3),DXIDX(3,3),DZ,DZDX(2),DZDXI(2),
     '  LENGTH,PTS(3,2),sin_THETA,cos_THETA,
     '  X1(3),X2(3),Z1(3),Z2(3),Z3(3),Z4(3),Z5(3)

      LOGICAL DERIV1,DERIV2,PLOT_REACTION,SCALAR

      CALL ENTERS('SGREAC',*9999)
      nv=1 !Temporary
      CALL OPEN_SEGMENT(ISREAC,ISEG,iw,'REAC',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)

      no_nynr=0
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
C       For scalar case with linear basis find line numbers
C       adjoining node np
        IF(NHP(np).EQ.1.AND.NKH(NH_LOC(1,nx),np,1).EQ.1) THEN
          SCALAR=.TRUE.
          IF(DOP) THEN
            WRITE(OP_STRING,'('' np='',I4)') np
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          DERIV1=.FALSE.
          DERIV2=.FALSE.
          DO nl=1,NLT
            IF(.NOT.DERIV1.OR..NOT.DERIV2) THEN
              IF(NPL(2,1,nl).eq.np) THEN
                IF(.NOT.DERIV1.AND.NPL(1,0,nl).EQ.1) THEN
                  n1=np        !n1 is 1st node along nl
                  n2=NPL(3,1,nl) !& n2 is 2nd node in Xi(1) direction
                  DZDXI(1)=ZP(1,1,1,n2,nc)-ZP(1,1,1,n1,nc)
                  DO nj=1,NJT
                    DXDXI(nj,1)=XP(1,nv,nj,n2)-XP(1,nv,nj,n1)
                  ENDDO
                  DERIV1=.TRUE.
                ELSE IF(.NOT.DERIV2.AND.NPL(1,0,nl).EQ.2) THEN
                  n1=np        !n1 is 1st node along nl
                  n2=NPL(3,1,nl) !& n2 is 2nd node in Xi(2) direction
                  DZDXI(2)=ZP(1,1,1,n2,nc)-ZP(1,1,1,n1,nc)
                  DO nj=1,NJT
                    DXDXI(nj,2)=XP(1,nv,nj,n2)-XP(1,nv,nj,n1)
                  ENDDO
                  DERIV2=.TRUE.
                ENDIF
                IF(DOP) THEN
                  WRITE(OP_STRING,
     '              '('' nl='',I4,'' n1='',I4,'' n2='',I4)')
     '              nl,n1,n2
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ELSE IF(NPL(3,1,nl).eq.np) THEN
                IF(.NOT.DERIV1.AND.NPL(1,0,nl).EQ.1) THEN
                  n1=NPL(2,1,nl) !n1 is 1st node along nl
                  n2=np        !& n2 is 2nd node in Xi(1) direction
                  DZDXI(1)=ZP(1,1,1,n2,nc)-ZP(1,1,1,n1,nc)
                  DO nj=1,NJT
                    DXDXI(nj,1)=XP(1,nv,nj,n2)-XP(1,nv,nj,n1)
                  ENDDO
                  DERIV1=.TRUE.
                ELSE IF(.NOT.DERIV2.AND.NPL(1,0,nl).EQ.2) THEN
                  n1=NPL(2,1,nl) !n1 is 1st node along nl
                  n2=np        !& n2 is 2nd node in Xi(2) direction
                  DZDXI(2)=ZP(1,1,1,n2,nc)-ZP(1,1,1,n1,nc)
                  DO nj=1,NJT
                    DXDXI(nj,2)=XP(1,nv,nj,n2)-XP(1,nv,nj,n1)
                  ENDDO
                  DERIV2=.TRUE.
                ENDIF
                IF(DOP) THEN
                  WRITE(OP_STRING,
     '              '('' nl='',I4,'' n1='',I4,'' n2='',I4)') nl,n1,n2
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
            ENDIF
          ENDDO !nl
          CALL INVERT(2,DXDXI,DXIDX,DET)
          DO nj=1,NJT
            DZDX(nj)=DZDXI(1)*DXIDX(1,nj)+DZDXI(2)*DXIDX(2,nj)
          ENDDO
        ELSE
          SCALAR=.FALSE.
        ENDIF !nhp=1 & nhk=1

        IF(SCALAR) THEN
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            DO nk=1,NKH(nh,np,nc)
              no_nynr=no_nynr+1
              ny=NYNR(no_nynr,0) !is global variable#
              IF(ITYP6(nr,nx).EQ.1) THEN       !linear
                IF(.NOT.CONSTR.OR.(CONSTR.AND.(FIX(ny,1).
     '            OR.FIX(ny,2)))) THEN
                  PLOT_REACTION=.TRUE.
                ELSE
                  PLOT_REACTION=.FALSE.
                ENDIF
              ELSE IF(ITYP6(nr,nx).EQ.2) THEN  !nonlinear
                IF(.NOT.CONSTR.OR.(CONSTR.AND.(FIX(ny,2).
     '            OR.FIX(ny,3)))) THEN
                  PLOT_REACTION=.TRUE.
                ELSE
                  PLOT_REACTION=.FALSE.
                ENDIF
              ENDIF !ityp6
              IF(nk.EQ.1.AND.PLOT_REACTION) THEN
                IF(ITYP6(nr,nx).EQ.1) THEN       !linear
! Arrow length is scaled by reaction YP(ny,iy) and user defined factor.
! Arrow direction is set by local solution gradient DZDX(nj).
                  DO nj=1,NJT
                    LENGTH=-FACTOR*YP(ny,IY)*DZDX(nj)
                    X1(nj)=XP(1,nv,nj,np)+LENGTH !is arrow tail
                    X2(nj)=XP(1,nv,nj,np) !is arrow head
                  ENDDO
                ELSE IF(ITYP6(nr,nx).EQ.2) THEN  !nonlinear
                  X2(nh)=ZP(1,nv,nh,np,nc)+FACTOR*YP(ny,IY)
                  DO nj=1,NJT
                    X1(nj)=ZP(1,nv,nj,np,nc)
                    IF(nj.NE.nh) X2(nj)=ZP(1,nv,nj,np,nc)
                  ENDDO
                ENDIF !ityp6
                CALL XZ(ITYP10(1),X1,Z1)
                CALL XZ(ITYP10(1),X2,Z2)
                IF(nh.EQ.1) THEN
                    Z3(1)=Z2(1)+FACTOR*( DZDX(1)/5.0D0-DZDX(2)/15.0D0)
                    Z3(2)=Z2(2)+FACTOR*( DZDX(2)/5.0D0+DZDX(1)/15.0D0)
                    Z3(3)=0.0D0
                    Z4(1)=Z2(1)+FACTOR*( DZDX(1)/5.0D0+DZDX(2)/15.0D0)
                    Z4(2)=Z2(2)+FACTOR*( DZDX(2)/5.0D0-DZDX(1)/15.0D0)
                    Z4(3)=0.0D0
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
                ENDIF !nh
                IF(iw.EQ.3) THEN
                  CALL ZZ(Z1,Z1,TRANS) !not be needed now?. AAY
                  CALL ZZ(Z2,Z2,TRANS)
                  CALL ZZ(Z3,Z3,TRANS)
                  CALL ZZ(Z4,Z4,TRANS)
                ENDIF

C ***           Draw reaction vector
                PTS(1,1)=Z1(1)
                PTS(1,2)=Z2(1)
                PTS(2,1)=Z1(2)
                PTS(2,2)=Z2(2)
                PTS(3,1)=Z1(3)
                PTS(3,2)=Z2(3)
                CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)

C ***           Draw arrow head
                PTS(1,1)=Z2(1)
                PTS(1,2)=Z3(1)
                PTS(2,1)=Z2(2)
                PTS(2,2)=Z3(2)
                PTS(3,1)=Z2(3)
                PTS(3,2)=Z3(3)
                CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)
                PTS(1,1)=Z2(1)
                PTS(1,2)=Z4(1)
                PTS(2,1)=Z2(2)
                PTS(2,2)=Z4(2)
                PTS(3,1)=Z2(3)
                PTS(3,2)=Z4(3)
                CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)
              ENDIF !nk=1 & plot_reaction
            ENDDO !nk
          ENDDO !nh

! PJH 15Apr1994 Draw force vectors rather than just components
        ELSE IF(.NOT.SCALAR) THEN !draw single force vector
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            DO nk=1,NKH(nh,np,nc)
              no_nynr=no_nynr+1
              ny=NYNR(no_nynr,0) !is global variable#
              IF(ITYP6(nr,nx).EQ.1) THEN !linear
                IF(.NOT.CONSTR.OR.(CONSTR.AND.(FIX(ny,1).
     '            OR.FIX(ny,2)))) THEN
                  PLOT_REACTION=.TRUE.
                ELSE
                  PLOT_REACTION=.FALSE.
                ENDIF
              ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear
                IF(.NOT.CONSTR.OR.(CONSTR.AND.(FIX(ny,2).
     '            OR.FIX(ny,3)))) THEN
                  PLOT_REACTION=.TRUE.
                ELSE
                  PLOT_REACTION=.FALSE.
                ENDIF
              ENDIF
              IF(nk.EQ.1.AND.PLOT_REACTION) THEN
                IF(ITYP6(nr,nx).EQ.1) THEN       !linear
                  X1(nh)=XP(1,nv,nh,np)
C moved                  X2(nh)=X1(nh)+FACTOR*YP(ny,IY)
                ELSE IF(ITYP6(nr,nx).EQ.2) THEN  !nonlinear
                  X1(nh)=ZP(1,nv,nh,np,1) !head at deformed node pos.
C wrong                  X2(nh)=X1(nh)+FACTOR*ZP(1,nv,nh,np,nc) !tail
                ENDIF
                X2(nh)=X1(nh)+FACTOR*YP(ny,IY) !tail

              ENDIF
            ENDDO !nk
          ENDDO !nh
          CALL XZ(ITYP10(nr),X1,Z1)
          CALL XZ(ITYP10(nr),X2,Z2)

C ***     Draw reaction vector
          PTS(1,1)=Z1(1)
          PTS(1,2)=Z2(1)
          PTS(2,1)=Z1(2)
          PTS(2,2)=Z2(2)
          PTS(3,1)=Z1(3)
          PTS(3,2)=Z2(3)
          CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)

C ***     Draw arrow head
!         Arrow length is scaled by reaction YP(ny,iy)
!         and user defined factor.
          LENGTH=DSQRT((Z2(1)-Z1(1))**2+(Z2(2)-Z1(2))**2)
          IF(LENGTH.GT.1.D-12) THEN
            DZ=LENGTH/4.0D0 !is arrow head size
!           Z3(j) are coords of pt on vector at start of head
            Z3(1)=Z1(1)+(LENGTH-DZ*0.866d0)/LENGTH*(Z2(1)-Z1(1))
            Z3(2)=Z1(2)+(LENGTH-DZ*0.866d0)/LENGTH*(Z2(2)-Z1(2))
            sin_THETA=(Z2(2)-Z1(2))/LENGTH
            cos_THETA=(Z2(1)-Z1(1))/LENGTH
!           Z4(j) & Z5(j) are coords at end of arrow head
            Z4(1)=Z3(1)+0.5d0*DZ*sin_THETA
            Z4(2)=Z3(2)-0.5d0*DZ*cos_THETA
            Z5(1)=Z3(1)-0.5d0*DZ*sin_THETA
            Z5(2)=Z3(2)+0.5d0*DZ*cos_THETA
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Z3(nj)='',3E12.4)')
     '          (Z3(nj),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Z4(nj)='',3E12.4)')
     '          (Z4(nj),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Z5(nj)='',3E12.4)')
     '          (Z5(nj),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            PTS(1,1)=Z2(1)
            PTS(1,2)=Z4(1)
            PTS(2,1)=Z2(2)
            PTS(2,2)=Z4(2)
            PTS(3,1)=Z2(3)
            PTS(3,2)=Z4(3)
            CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)
            PTS(1,1)=Z2(1)
            PTS(1,2)=Z5(1)
            PTS(2,1)=Z2(2)
            PTS(2,2)=Z5(2)
            PTS(3,1)=Z2(3)
            PTS(3,2)=Z5(3)
            CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)

          ELSE
!           WRITE(OP_STRING,'('' Arrow length is zero'')')
!           CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF !length>0
        ENDIF !scalar
      ENDDO !nonode

      CALL CLOSE_SEGMENT(ISREAC,iw,ERROR,*9999)

      CALL EXITS('SGREAC')
      RETURN
 9999 CALL ERRORS('SGREAC',ERROR)
      CALL EXITS('SGREAC')
      RETURN 1
      END

