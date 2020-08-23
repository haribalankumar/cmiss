      SUBROUTINE SHEET1(INDEX,IBT,IDO,INP,IW,NAN,NBJ,ne,NEELEM,
     '  NKJE,NPF,NPNE,NRE,NVJE,NXI,DXI2,DXI3,SE,THETA,XA,XE,XG,XP,
     '  ERROR,*)

C#### Subroutine: SHEET1
C###  Description:
C###    SHEET1 draws fitted sheet angles on constant Xi_1,X plane.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),IW,
     '  NAN(NIM,NAM,NBFM),NBJ(NJM,NEM),ne,NEELEM(0:NE_R_M,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 DXI2,DXI3,SE(NSM,NBFM,NEM),THETA,
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER j,n,n1elem,nb,NE1,NE_ADJACENT,ni,nj,iXI2,iXI3,nXI2,nXI3
      REAL*8 FITTED_ALFA,FITTED_GAMA,
     '  FITTED_A_VECTOR(3),FITTED_B_VECTOR(3),FITTED_C_VECTOR(3),
     '  GAMA_PROJ,CENTRE(2),
     '  DELTA,DTHETA,DX,DY,DZ,F_VECTOR(3),
     '  NORM_LINE(3),POINTS(3,21),PXI,C_xfv(3),
     '  THETA1,THETA2,X(3),XI(3),XI2,XI3,YZ(3),Z(3)

      CALL ENTERS('SHEET1',*9999)
      CALL ASSERT(IW.EQ.13,
     '  ' Incorrect workstation ID: sb 13 for sheets',ERROR,*9999)

      nXI2=NINT(1.0d0/DXI2)
      nXI3=NINT(1.0d0/DXI3)

C     DO XI2=0.d0,1.d0,DXI2
      DO iXI2=0,nXI2
        XI2=DBLE(iXI2)/DBLE(nXI2)

        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)

        XI(1)=0.d0
        XI(2)=XI2
        XI(3)=0.d0 !inner edge of endocardial element
        nj=3 !theta coordinate
        nb=NBJ(nj,ne)
        THETA1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '    XE(1,nj))
        XI(1)=1.d0
        THETA2=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '    XE(1,nj))
!        write(*,'('' theta='',E12.3,'' theta1='',E12.3,'' theta2='',E12.3)')
!     '    THETA,THETA1,THETA2
        IF(THETA1.GT.THETA.AND.THETA.GT.THETA2) THEN
          NE_ADJACENT=NXI(3,1,ne)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Element '',I4,'
     '        //''' adjacent element is '',I4)') ne,NE_ADJACENT
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          DTHETA=THETA1-THETA2
          IF(DABS(DTHETA).GT.1.d-6) THEN
            XI(1)=(THETA1-THETA)/DTHETA
          ELSE
            XI(1)=0.d0
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,'('' inner edge of endo element'',i4,'
     '        //''' xi(1)='',F5.3,'' xi(2)='',F5.3)') ne,XI(1),XI(2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          N=0
!          DO XI3=DXI3,1.d0-DXI3,DXI3
C          DO XI3=DXI3,1.d0,DXI3
          DO iXI3=1,nXI3
            XI3=DBLE(iXI3)/DBLE(nXI3)
            XI(3)=XI3
            N=N+1
            IF(DOP) THEN
              WRITE(OP_STRING,'('' N= '',I2)') N
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Xi:'',3E12.3)') (XI(ni),ni=1,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO j=1,3
              nb=NBJ(j,ne)
              X(j)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '          XE(1,j))
            ENDDO
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Interpolated coords:'',3E12.3)')
     '          (X(j),j=1,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL XZ(ITYP10(1),X,Z)
            CENTRE(1)=Z(1)
            IF(n.EQ.1) THEN
              CENTRE(2)=DSQRT(Z(2)**2+Z(3)**2)
              YZ(2) =Z(2)  !retains y coord
              YZ(3) =Z(3)  !retains z coord
            ELSE IF(n.GT.1) THEN
              DY=Z(2)-YZ(2) !is change in y from last point
              DZ=Z(3)-YZ(3) !is change in z from last point
              IF(DOP) THEN
                WRITE(OP_STRING,'('' dY,dZ:'',2E12.3)') DY,DZ
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CENTRE(2)=CENTRE(2)+DSQRT(DY**2+DZ**2)
              YZ(2) =Z(2)  !retains y coord
              YZ(3) =Z(3)  !retains z coord
            ENDIF

C new MAT_VEC_XI does not pass back
C FITTED_ALFA,FITTED_BETA,FITTED_GAMA,DELTA,PHI any more.
            CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),NRE(ne),
     '        FITTED_A_VECTOR,FITTED_B_VECTOR,FITTED_C_VECTOR,
     '        XE,XG,XI,.TRUE.,ERROR,*9999)

            CALL ASSERT(.FALSE.,'>> Old code - needs updating',
     '        ERROR,*9999)

C            CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),NRE(ne),
C     '        FITTED_A_VECTOR,FITTED_B_VECTOR,FITTED_C_VECTOR,
C     '        FITTED_ALFA,FITTED_BETA,FITTED_GAMA,DELTA,PHI,
C     '        XE,XI,.TRUE.,ERROR,*9999)

!           Ian Le Grice 28-3-92
!           Define line of intersection of myocardial sheet and
!           (v,x)-plane in x,f,v coordinate system

            F_VECTOR(1)= 0.d0
            F_VECTOR(2)= 1.d0
            F_VECTOR(3)= 0.d0

            C_xfv(1)= DCOS(FITTED_ALFA)*DCOS(FITTED_GAMA)*DCOS(DELTA)
     '               -DSIN(FITTED_GAMA)*DSIN(DELTA)
            C_xfv(2)= DSIN(FITTED_ALFA)*DCOS(FITTED_GAMA)
            C_xfv(3)= DCOS(FITTED_ALFA)*DCOS(FITTED_GAMA)*DSIN(DELTA)
     '               +DSIN(FITTED_GAMA)*DCOS(DELTA)

            NORM_LINE(1)= F_VECTOR(2)*C_xfv(3)-F_VECTOR(3)*C_xfv(2)
            NORM_LINE(2)= F_VECTOR(3)*C_xfv(1)-F_VECTOR(1)*C_xfv(3)
            NORM_LINE(3)= F_VECTOR(1)*C_xfv(2)-F_VECTOR(2)*C_xfv(1)

            GAMA_proj=DATAN2(-NORM_LINE(1),NORM_LINE(3))

            WRITE(OP_STRING,'('' C_xfv_1='',E12.3,'
     '        //''' C_xfv_2='',E12.3,'' C_xfv_3='',E12.3,E12.3)')
     '        C_xfv(1),C_xfv(2),C_xfv(3)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

            WRITE(OP_STRING,'('' norm_1='',E12.3,'' norm_2='',E12.3,'
     '        //''' norm_3='',E12.3, '' gama_proj='',E12.3)')
     '        NORM_LINE(1),NORM_LINE(2),NORM_LINE(3),GAMA_PROJ
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!            !H_vector in xfv coords
!            H_xfv(1)=DSIN(DELTA)
!            H_xfv(2)=0.d0
!            H_xfv(3)=DCOS(DELTA)
!            !GAMA vector in xfv coords
!            GAMA_xfv(1)=-DCOS(ALFA)*DSIN(GAMA)*DCOS(DELTA)+DCOS(GAMA)*DSIN(DELTA)
!            GAMA_xfv(3)= DCOS(ALFA)*DSIN(GAMA)*DSIN(DELTA)+DCOS(GAMA)*DCOS(DELTA)
!            ABS_GAMA_xv=DSQRT(GAMA_xfv(1)**2+GAMA_xfv(3)**2)
!            !unit projection of GAMA vector in xv-plane
!            GAMA_xv(1)=GAMA_xfv(1)/ABS_GAMA_xv
!            GAMA_xv(2)=GAMA_xfv(3)/ABS_GAMA_xv
!            GAMA_proj=DATAN2(-GAMA_xv(1),GAMA_xv(2))

            DX=0.02D0*FOCUS
            POINTS(1,1)=CENTRE(1)-DX*DSIN(GAMA_proj)
            POINTS(2,1)=CENTRE(2)+DX*DCOS(GAMA_proj)
            POINTS(1,2)=CENTRE(1)+DX*DSIN(GAMA_proj)
            POINTS(2,2)=CENTRE(2)-DX*DCOS(GAMA_proj)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Centre: '',2E12.3)') CENTRE(1),
     '          CENTRE(2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL POLYLINE(INDEX,IW,2,POINTS,ERROR,*9999)
          ENDDO

          IF(NE_ADJACENT.GT.0) THEN !draw sheet angles on adjacent element
            NE1=NE_ADJACENT          !.. at same Xi(1) and Xi(2) positions
            CALL XPXE(NBJ(1,NE1),NKJE(1,1,1,NE1),NPF(1,1),
     '        NPNE(1,1,NE1),NRE(ne1),NVJE(1,1,1,ne1),
     '        SE(1,1,NE1),XA(1,1,ne1),XE,XP,ERROR,*9999)
            n=0
C           DO XI3=DXI3,1.0d0,DXI3
            DO iXI3=1,nXI3
              XI3=DBLE(iXI3)/DBLE(nXI3)
              XI(3)=XI3
              n=n+1
              DO j=1,3
                nb=NBJ(j,NE1)
                X(j)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '            XE(1,j))
              ENDDO
              CALL XZ(ITYP10(1),X,Z)
              DY=Z(2)-YZ(2) !is change in y from last point
              DZ=Z(3)-YZ(3) !is change in z from last point
              CENTRE(1)=Z(1)
              CENTRE(2)=CENTRE(2)+DSQRT(DY**2+DZ**2)
              YZ(2) =Z(2)  !retains y coord
              YZ(3) =Z(3)  !retains z coord

C new MAT_VEC_XI does not pass back
C FITTED_ALFA,FITTED_BETA,FITTED_GAMA,DELTA,PHI any more.
              CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),NRE(ne),
     '          FITTED_A_VECTOR,FITTED_B_VECTOR,FITTED_C_VECTOR,
     '          XE,XG,XI,.TRUE.,ERROR,*9999)

              CALL ASSERT(.FALSE.,'>> Old code - needs updating',
     '          ERROR,*9999)

C              CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),NRE(ne),
C     '          FITTED_A_VECTOR,FITTED_B_VECTOR,FITTED_C_VECTOR,
C     '          FITTED_ALFA,FITTED_BETA,FITTED_GAMA,DELTA,PHI,
C     '          XE,XI,.TRUE.,ERROR,*9999)

!           Ian Le Grice 28-3-92
!           Define line of intersection of myocardial sheet and
!           (v,x)-plane in x,g,v coordinate system

              F_VECTOR(1)= 0.d0
              F_VECTOR(2)= 1.d0
              F_VECTOR(3)= 0.d0

              C_xfv(1)= DCOS(FITTED_ALFA)*DCOS(FITTED_GAMA)*DCOS(DELTA)
     '                 +DSIN(FITTED_GAMA)*DSIN(DELTA)
              C_xfv(2)= DSIN(FITTED_ALFA)*DCOS(FITTED_GAMA)
              C_xfv(3)=-DCOS(FITTED_ALFA)*DCOS(FITTED_GAMA)*DSIN(DELTA)
     '                 +DSIN(FITTED_GAMA)*DCOS(DELTA)

              NORM_LINE(1)= F_VECTOR(2)*C_xfv(3)-F_VECTOR(3)*C_xfv(2)
              NORM_LINE(2)= F_VECTOR(3)*C_xfv(1)-F_VECTOR(1)*C_xfv(3)
              NORM_LINE(3)= F_VECTOR(1)*C_xfv(2)-F_VECTOR(2)*C_xfv(1)

              GAMA_proj=DATAN2(-NORM_LINE(1),NORM_LINE(3))

              WRITE(OP_STRING,'('' C_xfv_1='',E12.3,'
     '          //''' C_xfv_2='',E12.3,'' C_xfv_3='',E12.3,E12.3)')
     '          C_xfv(1),C_xfv(2),C_xfv(3)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

              WRITE(OP_STRING,'('' norm_1='',E12.3,'' norm_2='','
     '          //'E12.3,'' norm_3='',E12.3,'' gama_proj='',E12.3)')
     '          NORM_LINE(1),NORM_LINE(2),NORM_LINE(3),GAMA_PROJ
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!              !U vector in xfv coords
!              H_xfv(1)=DSIN(DELTA)
!              H_xfv(2)=0.d0
!              H_xfv(3)=DCOS(DELTA)
!              !GAMA vector in xfv coords
!              GAMA_xfv(1)=-DCOS(ALFA)*DSIN(GAMA)*DCOS(DELTA)+DCOS(GAMA)*DSIN(DELTA)
!              GAMA_xfv(3)= DCOS(ALFA)*DSIN(GAMA)*DSIN(DELTA)+DCOS(GAMA)*DCOS(DELTA)
!              ABS_GAMA_xv=DSQRT(GAMA_xfv(1)**2+GAMA_xfv(3)**2)
!              !unit projection of GAMA vector in xv-plane
!              GAMA_xv(1)=GAMA_xfv(1)/ABS_GAMA_xv
!              GAMA_xv(2)=GAMA_xfv(3)/ABS_GAMA_xv
!              GAMA_proj=DATAN2(-GAMA_xv(1),GAMA_xv(2))

              DX=0.02D0*FOCUS
              POINTS(1,1)=CENTRE(1)-DX*DSIN(GAMA_proj)
              POINTS(2,1)=CENTRE(2)+DX*DCOS(GAMA_proj)
              POINTS(1,2)=CENTRE(1)+DX*DSIN(GAMA_proj)
              POINTS(2,2)=CENTRE(2)-DX*DCOS(GAMA_proj)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Centre: '',2E12.3)') CENTRE(1),
     '            CENTRE(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL POLYLINE(INDEX,IW,2,POINTS,ERROR,*9999)
            ENDDO

          ELSE IF(NE_ADJACENT.EQ.0) THEN !ne is septum and need rv free wall
            THETA1=XE(1,3) !is theta at 1st node of endo element
            THETA2=XE(2,3) !is theta at 2nd node of endo element
            DO n1elem=1,NEELEM(0,1)
              NE1=NEELEM(n1elem,1)
              IF(NET(1).GT.30.AND.NE1.LE.30) THEN !restrict to epi TEMPORARY!!!!
                CALL XPXE(NBJ(1,NE1),NKJE(1,1,1,NE1),
     '            NPF(1,1),NPNE(1,1,NE1),NRE(ne1),
     '            NVJE(1,1,1,ne1),SE(1,1,NE1),XA(1,1,ne1),XE,XP,ERROR,
     '            *9999)
                IF(DABS(XE(1,3)-THETA1).LT.1.d-2.AND
     '            .DABS(XE(2,3)-THETA2).LT.1.d-2) THEN
                  WRITE(OP_STRING,'('' Found rv free wall element '','
     '              //'I4)') NE1
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GO TO 20
                ENDIF
              ENDIF
            ENDDO
 20         CONTINUE

            n=0
C           DO XI3=0.d0,1.d0,DXI3
            DO iXI3=0,nXI3
              XI3=DBLE(iXI3)/DBLE(nXI3)
              XI(3)=XI3
              n=n+1
              DO j=1,3
                nb=NBJ(j,NE1)
                X(j)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '            XE(1,j))
              ENDDO
              CALL XZ(ITYP10(1),X,Z)
              DY=Z(2)-YZ(2) !is change in y from last point
              DZ=Z(3)-YZ(3) !is change in z from last point
              CENTRE(1)=Z(1)
              CENTRE(2)=CENTRE(2)+DSQRT(DY**2+DZ**2)
              YZ(2) =Z(2)  !retains y coord
              YZ(3) =Z(3)  !retains z coord

C new MAT_VEC_XI does not pass back
C FITTED_ALFA,FITTED_BETA,FITTED_GAMA,DELTA,PHI any more.
            CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),NRE(ne),
     '        FITTED_A_VECTOR,FITTED_B_VECTOR,FITTED_C_VECTOR,
     '        XE,XG,XI,.TRUE.,ERROR,*9999)

            CALL ASSERT(.FALSE.,'>> Old code - needs updating',
     '        ERROR,*9999)

C              CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),NRE(ne),
C     '          FITTED_A_VECTOR,FITTED_B_VECTOR,FITTED_C_VECTOR,
C     '          FITTED_ALFA,FITTED_BETA,FITTED_GAMA,DELTA,PHI,
C     '          XE,XI,.TRUE.,ERROR,*9999)
!           Ian Le Grice 28-3-92
!           Define line of intersection of myocardial sheet and
!          (v,x)-plane in x,g,v coordinate system

              F_VECTOR(1)= 0.d0
              F_VECTOR(2)= 1.d0
              F_VECTOR(3)= 0.d0

              C_xfv(1)= DCOS(FITTED_ALFA)*DCOS(FITTED_GAMA)*DCOS(DELTA)
     '                  +DSIN(FITTED_GAMA)*DSIN(DELTA)
              C_xfv(2)= DSIN(FITTED_ALFA)*DCOS(FITTED_GAMA)
              C_xfv(3)=-DCOS(FITTED_ALFA)*DCOS(FITTED_GAMA)*DSIN(DELTA)
     '                  +DSIN(FITTED_GAMA)*DCOS(DELTA)

              NORM_LINE(1)= F_VECTOR(2)*C_xfv(3)-F_VECTOR(3)*C_xfv(2)
              NORM_LINE(2)= F_VECTOR(3)*C_xfv(1)-F_VECTOR(1)*C_xfv(3)
              NORM_LINE(3)= F_VECTOR(1)*C_xfv(2)-F_VECTOR(2)*C_xfv(1)

              GAMA_proj=DATAN2(-NORM_LINE(1),NORM_LINE(3))

              WRITE(OP_STRING,'('' C_xfv_1='',E12.3,'
     '          //''' C_xfv_2='',E12.3,'' C_xfv_3='',E12.3,E12.3)')
     '          C_xfv(1),C_xfv(2),C_xfv(3)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

              WRITE(OP_STRING,'('' norm_1='',E12.3,'' norm_2='','
     '          //'E12.3,'' norm_3='',E12.3, '' gama_proj='',E12.3)')
     '          NORM_LINE(1),NORM_LINE(2),NORM_LINE(3),GAMA_PROJ
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

!              !U vector in xfv coords
!              H_xfv(1)=DSIN(DELTA)
!              H_xfv(2)=0.d0
!              H_xfv(3)=DCOS(DELTA)
!              !GAMA vector in xfv coords
!              GAMA_xfv(1)=-DCOS(ALFA)*DSIN(GAMA)*DCOS(DELTA)+DCOS(GAMA)*DSIN(DELTA)
!              GAMA_xfv(3)= DCOS(ALFA)*DSIN(GAMA)*DSIN(DELTA)+DCOS(GAMA)*DCOS(DELTA)
!              ABS_GAMA_xv=DSQRT(GAMA_xfv(1)**2+GAMA_xfv(3)**2)
!              !unit projection of GAMA vector in xv-plane
!              GAMA_xv(1)=GAMA_xfv(1)/ABS_GAMA_xv
!              GAMA_xv(2)=GAMA_xfv(3)/ABS_GAMA_xv
!              GAMA_proj=DATAN2(-GAMA_xv(1),GAMA_xv(2))

              DX=0.02D0*FOCUS
              POINTS(1,1)=CENTRE(1)-DX*DSIN(GAMA_proj)
              POINTS(2,1)=CENTRE(2)+DX*DCOS(GAMA_proj)
              POINTS(1,2)=CENTRE(1)+DX*DSIN(GAMA_proj)
              POINTS(2,2)=CENTRE(2)-DX*DCOS(GAMA_proj)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Centre: '',2E12.3)') CENTRE(1),
     '            CENTRE(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL POLYLINE(INDEX,IW,2,POINTS,ERROR,*9999)
            ENDDO
          ENDIF

        ENDIF
      ENDDO

C LC 24/2/97 archived section :
C     CPB 30/3/93 Not sure about the nj locations for this commented section

      CALL EXITS('SHEET1')
      RETURN
 9999 CALL ERRORS('SHEET1',ERROR)
      CALL EXITS('SHEET1')
      RETURN 1
      END


