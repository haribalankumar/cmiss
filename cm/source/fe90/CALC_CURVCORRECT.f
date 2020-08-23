      SUBROUTINE CALC_CURVCORRECT(IBT,IDO,INP,nb,ne,NPNE,nr,CE,
     '  CURVCORRECT,DXNDS,SE,XE,XP,INTERFACE,ERROR,*)

C#### Subroutine: CALC_CURVCORRECT
C###  Description:
C###    CALC_CURVCORRECT calculates the curvature corrections when
C###    using cubic Hermite basis functions to interpolate the normal
C###    derivative. If INTERFACE is .TRUE. then the element is in
C###    the slave region of an interface and hence needs to have its
C###    normal reversed.

C*** Note: This subroutine assumes the same basis function in each
C***       geometric direction.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,ne,
     '  NPNE(NNM,NBFM),nr
      REAL*8 CE(NMM),CURVCORRECT(2,2,NNM),DXNDS(3,2),SE(NSM,NBFM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL INTERFACE
!     Local Variables
      INTEGER i,j,M(5),ni,nj,nn,nn1,np,ns,nu
      REAL*8 DET,DSDX(3,3),DXDS(3,3),DXDXI(3,2),D2XDXI(3,2),
     '  D2XDXI1DXI2(3),DXNDXI(3,2),NORMS,PXI,XI(2),
     '  XNBAR(3),DXNBAR(3,2),XNN(3)

      DATA M/1,2,3,1,2/

      CALL ENTERS('CALC_CURVCORRECT',*9999)

      DO nn=1,NNT(nb)
        IF(NKT(nn,nb).GT.1) THEN
          IF(NIT(nb).EQ.1) THEN ! cubic Hermite elements
            XI(1)=DBLE(INP(nn,1)-1)
            np=NPNE(nn,nb)
            NORMS=DSQRT(XP(2,1,1,np)**2+XP(2,1,2,np)**2)
            DXDS(1,1)=XP(2,1,1,np)/NORMS ! dX_1/ds
            DXDS(2,1)=XP(2,1,2,np)/NORMS ! dX_2/ds
            XNN(1)=-DXDS(2,1) ! n_1
            XNN(2)=DXDS(1,1) ! n_2
            DXDS(1,2)=XNN(1)
            DXDS(2,2)=XNN(2)
            DET=DXDS(1,1)*DXDS(2,2)-DXDS(1,2)*DXDS(2,1)
            DSDX(1,1)=DXDS(2,2)/DET
            DSDX(1,2)=-DXDS(1,2)/DET
            DSDX(2,1)=-DXDS(2,1)/DET
            DSDX(2,2)=DXDS(1,1)/DET
            DXDXI(1,1)=PXI(IBT,IDO,INP,nb,2,XI,XE(1,1))
            D2XDXI(1,1)=PXI(IBT,IDO,INP,nb,3,XI,XE(1,1))
            DXDXI(2,1)=PXI(IBT,IDO,INP,nb,2,XI,XE(1,2))
            D2XDXI(2,1)=PXI(IBT,IDO,INP,nb,3,XI,XE(1,2))
            NORMS=DSQRT(DXDXI(1,1)**2+DXDXI(2,1)**2)
            DXNDXI(1,1)=((-1.0d0+XNN(1)**2)*D2XDXI(2,1)-XNN(1)*XNN(2)*
     '        D2XDXI(1,1))/NORMS
            DXNDXI(2,1)=((1.0d0-XNN(2)**2)*D2XDXI(1,1)+XNN(1)*XNN(2)*
     '        D2XDXI(2,1))/NORMS
            ns=(nn-1)*2+2
            DXNDS(1,1)=DXNDXI(1,1)/SE(ns,nb)
            DXNDS(2,1)=DXNDXI(2,1)/SE(ns,nb)
            CURVCORRECT(1,1,nn)=DSDX(1,1)*DXNDS(1,1)+DSDX(1,2)*
     '        DXNDS(2,1)
            IF(IGREN(nr).EQ.7.OR.IGREN(nr).EQ.8) THEN ! Gen laplace
              CURVCORRECT(1,1,nn)=CURVCORRECT(1,1,nn)*CE(1)
            ENDIF
            IF(INTERFACE) THEN ! Reverse curvature correction
              CURVCORRECT(1,1,nn)=-CURVCORRECT(1,1,nn)
            ENDIF
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(CALC_CURVCORRECT_1)
              WRITE(OP_STRING,'('' 1 Xi direction: ne='',I5,'', nn='','
     '          //'I1,'', np='',I5)') ne,nn,np
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XI: '',D12.5)') XI(1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XNN:'',3(1X,D12.5))')
     '          (XNN(nj),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXDS: '',D12.5,1X,D12.5,/,7X,'
     '          //'D12.5,1X,D12.5)') DXDS(1,1),DXDS(1,2),DXDS(2,1),
     '          DXDS(2,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DSDX: '',D12.5,1X,D12.5,/,7X,'
     '          //'D12.5,1X,D12.5)') DSDX(1,1),DSDX(1,2),DSDX(2,1),
     '          DSDX(2,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXDXI:'',3(1X,D12.5))')
     '          (DXDXI(nj,1),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' D2XDXI:'',3(1X,D12.5))')
     '          (D2XDXI(nj,1),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NORM OF DXDXI='',D12.5)') NORMS
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXNDXI:'',3(1X,D12.5))')
     '          (DXNDXI(nj,1),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXNDS:'',3(1X,D12.5))')
     '          (DXNDS(nj,1),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' CURVCORRECT: '',D12.5)')
     '          CURVCORRECT(1,1,nn)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(CALC_CURVCORRECT_1)

            ENDIF
          ELSE IF(NIT(nb).EQ.2) THEN ! bicubic Hermite elements
            XI(1)=DBLE(INP(nn,1)-1)
            XI(2)=DBLE(INP(nn,2)-1)
            np=NPNE(nn,nb)
            NORMS=DSQRT(XP(2,1,1,np)**2+XP(2,1,2,np)**2+
     '        XP(2,1,3,np)**2)
            DXDS(1,1)=XP(2,1,1,np)/NORMS
            DXDS(2,1)=XP(2,1,2,np)/NORMS
            DXDS(3,1)=XP(2,1,3,np)/NORMS
            NORMS=DSQRT(XP(3,1,1,np)**2+XP(3,1,2,np)**2+
     '        XP(3,1,3,np)**2)
            DXDS(1,2)=XP(3,1,1,np)/NORMS
            DXDS(2,2)=XP(3,1,2,np)/NORMS
            DXDS(3,2)=XP(3,1,3,np)/NORMS
            XNN(1)=DXDS(2,1)*DXDS(3,2)-DXDS(3,1)*DXDS(2,2)
            XNN(2)=DXDS(3,1)*DXDS(1,2)-DXDS(1,1)*DXDS(3,2)
            XNN(3)=DXDS(1,1)*DXDS(2,2)-DXDS(2,1)*DXDS(1,2)
            NORMS=DSQRT(XNN(1)**2+XNN(2)**2+XNN(3)**2)
            XNN(1)=XNN(1)/NORMS
            XNN(2)=XNN(2)/NORMS
            XNN(3)=XNN(3)/NORMS
            DXDS(1,3)=XNN(1)
            DXDS(2,3)=XNN(2)
            DXDS(3,3)=XNN(3)
            DET=DXDS(1,1)*(DXDS(2,2)*DXDS(3,3)-DXDS(3,2)*DXDS(2,3))
     '        +DXDS(1,2)*(DXDS(2,3)*DXDS(3,1)-DXDS(3,3)*DXDS(2,1))
     '        +DXDS(1,3)*(DXDS(2,1)*DXDS(3,2)-DXDS(3,1)*DXDS(2,2))
            DO ni=1,3
              DO nj=1,3
                DSDX(ni,nj)=(DXDS(M(nj+1),M(ni+1))*DXDS(M(nj+2),
     '            M(ni+2))-DXDS(M(nj+2),M(ni+1))*DXDS(M(nj+1),
     '            M(ni+2)))/DET
              ENDDO !nj
            ENDDO !ni
            DO ni=1,2
              nu=ni*2
              DO nj=1,NJT
                DXDXI(nj,ni)=PXI(IBT,IDO,INP,nb,nu,XI,XE(1,nj))
                D2XDXI(nj,ni)=PXI(IBT,IDO,INP,nb,nu+1,XI,XE(1,nj))
                D2XDXI1DXI2(nj)=PXI(IBT,IDO,INP,nb,6,XI,XE(1,nj))
              ENDDO !nj
            ENDDO !ni
            XNBAR(1)=DXDXI(2,1)*DXDXI(3,2)-DXDXI(3,1)*DXDXI(2,2)
            XNBAR(2)=DXDXI(3,1)*DXDXI(1,2)-DXDXI(1,1)*DXDXI(3,2)
            XNBAR(3)=DXDXI(1,1)*DXDXI(2,2)-DXDXI(2,1)*DXDXI(1,2)
            NORMS=DSQRT(XNBAR(1)**2+XNBAR(2)**2+XNBAR(3)**2)
            DXNBAR(1,1)=D2XDXI(2,1)*DXDXI(3,2)+DXDXI(2,1)*
     '        D2XDXI1DXI2(3)-D2XDXI(3,1)*DXDXI(2,2)-DXDXI(3,1)*
     '        D2XDXI1DXI2(2)
            DXNBAR(1,2)=D2XDXI1DXI2(2)*DXDXI(3,2)+DXDXI(2,1)*
     '        D2XDXI(3,2)-D2XDXI1DXI2(3)*DXDXI(2,2)-DXDXI(3,1)*
     '        D2XDXI(2,2)
            DXNBAR(2,1)=D2XDXI(3,1)*DXDXI(1,2)+DXDXI(3,1)*
     '        D2XDXI1DXI2(1)-D2XDXI(1,1)*DXDXI(3,2)-DXDXI(1,1)*
     '        D2XDXI1DXI2(3)
            DXNBAR(2,2)=D2XDXI1DXI2(3)*DXDXI(1,2)+DXDXI(3,1)*
     '        D2XDXI(1,2)-D2XDXI1DXI2(1)*DXDXI(3,2)-DXDXI(1,1)*
     '        D2XDXI(3,2)
            DXNBAR(3,1)=D2XDXI(1,1)*DXDXI(2,2)+DXDXI(1,1)*
     '        D2XDXI1DXI2(2)-D2XDXI(2,1)*DXDXI(1,2)-DXDXI(2,1)*
     '        D2XDXI1DXI2(1)
            DXNBAR(3,2)=D2XDXI1DXI2(1)*DXDXI(2,2)+DXDXI(1,1)*
     '        D2XDXI(2,2)-D2XDXI1DXI2(2)*DXDXI(1,2)-DXDXI(2,1)*
     '        D2XDXI(1,2)
C            DXNDXI(1,1)=((XNN(1)**2-1.0d0)*DXNBAR(1,1)+XNN(1)*XNN(2)*
C     '        DXNBAR(2,1)+XNN(1)*XNN(3)*DXNBAR(3,1))/NORMS
C            DXNDXI(1,2)=((XNN(1)**2-1.0d0)*DXNBAR(1,2)+XNN(1)*XNN(2)*
C     '        DXNBAR(2,2)+XNN(1)*XNN(3)*DXNBAR(3,2))/NORMS
C            DXNDXI(2,1)=(XNN(1)*XNN(2)*DXNBAR(1,1)+(XNN(2)**2-1.0d0)*
C     '        DXNBAR(2,1)+XNN(2)*XNN(3)*DXNBAR(3,1))/NORMS
C            DXNDXI(2,2)=(XNN(1)*XNN(2)*DXNBAR(1,2)+(XNN(2)**2-1.0d0)*
C     '        DXNBAR(2,2)+XNN(2)*XNN(3)*DXNBAR(3,2))/NORMS
C            DXNDXI(3,1)=(XNN(1)*XNN(3)*DXNBAR(1,1)+XNN(2)*XNN(3)*
C     '        DXNBAR(2,1)+(XNN(3)**2-1.0d0)*DXNBAR(3,1))/NORMS
C            DXNDXI(3,2)=(XNN(1)*XNN(3)*DXNBAR(1,2)+XNN(2)*XNN(3)*
C     '        DXNBAR(2,2)+(XNN(3)**2-1.0d0)*DXNBAR(3,2))/NORMS
            DXNDXI(1,1)=((XNN(1)*DXNBAR(1,1)+XNN(2)*DXNBAR(2,1)+
     '        XNN(3)*DXNBAR(3,1))*XNN(1)-DXNBAR(1,1))/NORMS
            DXNDXI(1,2)=((XNN(1)*DXNBAR(1,2)+XNN(2)*DXNBAR(2,2)+
     '        XNN(3)*DXNBAR(3,2))*XNN(1)-DXNBAR(1,2))/NORMS
            DXNDXI(2,1)=((XNN(1)*DXNBAR(1,1)+XNN(2)*DXNBAR(2,1)+
     '        XNN(3)*DXNBAR(3,1))*XNN(2)-DXNBAR(2,1))/NORMS
            DXNDXI(2,2)=((XNN(1)*DXNBAR(1,2)+XNN(2)*DXNBAR(2,2)+
     '        XNN(3)*DXNBAR(3,2))*XNN(2)-DXNBAR(2,2))/NORMS
            DXNDXI(3,1)=((XNN(1)*DXNBAR(1,1)+XNN(2)*DXNBAR(2,1)+
     '        XNN(3)*DXNBAR(3,1))*XNN(3)-DXNBAR(3,1))/NORMS
            DXNDXI(3,2)=((XNN(1)*DXNBAR(1,2)+XNN(2)*DXNBAR(2,2)+
     '        XNN(3)*DXNBAR(3,2))*XNN(3)-DXNBAR(3,2))/NORMS
            ns=1
            DO nn1=1,nn-1
              ns=ns+NKT(nn1,nb)
            ENDDO !nn1
            DO j=1,2
              DXNDS(1,j)=DXNDXI(1,j)/SE(ns+j,nb)
              DXNDS(2,j)=DXNDXI(2,j)/SE(ns+j,nb)
              DXNDS(3,j)=DXNDXI(3,j)/SE(ns+j,nb)
            ENDDO !j
            DO i=1,2
              DO j=1,2
                CURVCORRECT(i,j,nn)=0.0d0
                DO nj=1,NJT
                  CURVCORRECT(i,j,nn)=CURVCORRECT(i,j,nn)+
     '              DSDX(i,nj)*DXNDS(nj,j)
                ENDDO !nj
              ENDDO !j
            ENDDO !i
            IF(IGREN(nr).EQ.7.OR.IGREN(nr).EQ.8) THEN ! Gen laplace
              DO i=1,2
                DO j=1,2
                  CURVCORRECT(i,j,nn)=CURVCORRECT(i,j,nn)*CE(1)
                ENDDO !j
              ENDDO !i
            ENDIF
            IF(INTERFACE) THEN !Reverse curvature correction
              CURVCORRECT(1,1,nn)=-CURVCORRECT(1,1,nn)
              CURVCORRECT(1,2,nn)=-CURVCORRECT(1,2,nn)
              CURVCORRECT(2,1,nn)=-CURVCORRECT(2,1,nn)
              CURVCORRECT(2,2,nn)=-CURVCORRECT(2,2,nn)
            ENDIF
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(CALC_CURVCORRECT_2)
              WRITE(OP_STRING,'('' 2 Xi directions: ne='',I5,'
     '          //''', nn='',I1,'', np='',I5)') ne,nn,np
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XI:'',2(1X,D12.5))') (XI(ni),ni=1,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XNN:'',3(1X,D12.5))')
     '          (XNN(nj),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,
     '          '('' DXDS:'',3(1X,D12.5),/,6X,3(1X,D12.5),'
     '          //'/,6X,3(1X,D12.5))') DXDS(1,1),DXDS(1,2),DXDS(1,3),
     '          DXDS(2,1),DXDS(2,2),DXDS(2,3),DXDS(3,1),DXDS(3,2),
     '          DXDS(3,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,
     '          '('' DSDX:'',3(1X,D12.5),/,6X,3(1X,D12.5),'
     '          //'/,6X,3(1X,D12.5))') DSDX(1,1),DSDX(1,2),DSDX(1,3),
     '          DSDX(2,1),DSDX(2,2),DSDX(2,3),DSDX(3,1),DSDX(3,2),
     '          DSDX(3,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXDXI(1..,1):'',3(1X,D12.5))')
     '          (DXDXI(nj,1),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXDXI(1..,2):'',3(1X,D12.5))')
     '          (DXDXI(nj,2),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' D2XDXI(1..,1):'',3(1X,D12.5))')
     '          (D2XDXI(nj,1),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' D2XDXI(1..,2):'',3(1X,D12.5))')
     '          (D2XDXI(nj,2),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' D2XDXI1DXI2:'',3(1X,D12.5))')
     '          (D2XDXI1DXI2(nj),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XNBAR:'',3(1X,D12.5))') (XNBAR(nj),
     '          nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NORM OF XNBAR='',D12.5)') NORMS
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXNBAR(1..,1):'',3(1X,D12.5))')
     '          (DXNBAR(nj,1),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXNBAR(1..,2):'',3(1X,D12.5))')
     '          (DXNBAR(nj,2),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXNDXI(1..,1):'',3(1X,D12.5))')
     '          (DXNDXI(nj,1),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXNDXI(1..,2):'',3(1X,D12.5))')
     '          (DXNDXI(nj,2),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXNDS(1..,1):'',3(1X,D12.5))')
     '          (DXNDS(nj,1),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DXNDS(1..,2):'',3(1X,D12.5))')
     '          (DXNDS(nj,2),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' CURVCORRECT:'',2(1X,D12.5),/,13X,'
     '          //'2(1X,D12.5))') CURVCORRECT(1,1,nn),
     '          CURVCORRECT(1,2,nn),CURVCORRECT(2,1,nn),
     '          CURVCORRECT(2,2,nn)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(CALC_CURVCORRECT_2)
            ENDIF
          ELSE
            ERROR='>>Invalid number of xi directions'
            GOTO 9999
          ENDIF
        ELSE
          CURVCORRECT(1,1,nn)=0.0d0
          CURVCORRECT(1,2,nn)=0.0d0
          CURVCORRECT(2,1,nn)=0.0d0
          CURVCORRECT(2,2,nn)=0.0d0
        ENDIF
      ENDDO ! nn

      CALL EXITS('CALC_CURVCORRECT')
      RETURN
9999  CALL ERRORS('CALC_CURVCORRECT',ERROR)
      CALL EXITS('CALC_CURVCORRECT')
      RETURN 1
      END


