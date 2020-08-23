      SUBROUTINE MG_SOLVE(na,NITB,NAQ,NWQ,NXQ,
     '  ISOTROPIC,TYPE,DT,GCHQ,GUQ,PROPQ,YQ,ERROR,*)

C#### Subroutine: MG_SOLVE
C###  Description:
C###    MGOLVE solves multigrid coarse grid equations.
C###    The solution is returned in YQ(nq,2,na).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER na,NITB,NAQ(NQM),NWQ(8,0:NQM),NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 DT,GCHQ(3,NQM),GUQ(3,3,NQM),PROPQ(3,3,4,2,NQM),
     '  YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*),TYPE*(*)
      LOGICAL ISOTROPIC
!     Local Variables
      INTEGER i,ij,ik,j,k,mq,nii,nij,nik,N_grid,nq,nq_middle
      REAL*8 AA,BB,CC,D,DD,dR,dUdX(3),EU,F(3,3),
     '  GCHQ1,GCHQ2,GUQ11,GUQ22,GUQ12,
     '  RR,SUM1,SUM2,U,UQ(-1:1,-1:1,-1:1)

      CALL ENTERS('MG_SOLVE',*9999)

      DO i=-1,1
        DO j=-1,1
          DO k=-1,1
            UQ(i,j,k)=0.0d0
          ENDDO
        ENDDO
      ENDDO

      N_grid=0
      DO nq=1,NQT !Find #interior grid pts on grid na
        IF(NAQ(nq).EQ.0.AND.NWQ(1,nq).EQ.0) THEN
          N_grid=N_grid+1
          nq_middle=nq
        ENDIF
      ENDDO
      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_SOLVE_1)
        WRITE(OP_STRING,'(/'' #interior gps on coarsest grid='','
     '    //'I4,'' nq_middle='',I4)') N_grid,nq_middle
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_SOLVE_1)
      ENDIF

      IF(N_grid.EQ.1) THEN !one point solution
        D=2.0D0/DBLE(2**na)
        DD=D*D
        DO nq=1,NQT !Loop over interior points of grid na
          IF(NAQ(nq).EQ.0.AND.NWQ(1,nq).EQ.0) THEN
            ik=MAX(0,NITB-2) !zero for 1,2D, one for 3D
            ij=MIN(NITB-1,1) !zero for 1D, one for 2,3D
            DO nik=-ik,ik    !UQ has coeffs of local quadratic
              DO nij=-ij,ij  !element about nq
                DO nii=-1,1
                  mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq)))!adj g.p.
                  IF(mq.GT.0) THEN
                    UQ(nii,nij,nik)=YQ(mq,1,na) !u at mq
                  ELSE
                    write(*,'('' Error: outside grid in MG_SOLVE'')')
                  ENDIF
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik
          ENDIF !nq is interior grid na
        ENDDO !nq

        nq=nq_middle

        IF(ISOTROPIC) THEN !2D isotropy
          GUQ11=GUQ(1,1,nq)*DD
          GUQ22=GUQ(2,2,nq)*DD
          GUQ12=GUQ(1,2,nq)*DD
          GCHQ1=GCHQ( 1,nq)*D
          GCHQ2=GCHQ( 2,nq)*D
          AA=4.0d0*GUQ11*(UQ(-1, 0, 0)+UQ( 1, 0, 0))
     '      +4.0d0*GUQ22*(UQ( 0,-1, 0)+UQ( 0, 1, 0))
     '      +2.0d0*GUQ12*(UQ(-1,-1, 0)+UQ( 1, 1, 0)
     '                   -UQ( 1,-1, 0)-UQ(-1, 1, 0))
     '            +GCHQ1*(UQ(-1, 0, 0)-UQ( 1, 0, 0))
     '            +GCHQ2*(UQ( 0,-1, 0)-UQ( 0, 1, 0))

          YQ(nq,2,na)=(YQ(nq,5,na)+AA)/(8.0D0*(GUQ11+GUQ22))
c         YQ(nq,2,na)=0.25D0*(DD*YQ(nq,3,na)
c    '      +UQ(-1,0,0)+UQ(1,0,0)+UQ(0,-1,0)+UQ(0,1,0))

        ELSE IF(.NOT.ISOTROPIC) THEN !general anisotropic
          DO i=1,NITB
            DO j=1,NITB
              SUM1=0.0d0
              DO k=1,NITB
                SUM1=SUM1+PROPQ(i,k,1,1,nq)*GUQ(j,k,nq)
              ENDDO !k
              F(i,j)=SUM1*DD
            ENDDO !j
          ENDDO !i

          dUdX(1)=UQ(1,0,0)-UQ(-1,0,0)
          dUdX(2)=UQ(0,1,0)-UQ(0,-1,0)
          dUdX(3)=UQ(0,0,1)-UQ(0,0,-1)
          EU=0.d0
          DO i=1,NITB
            SUM1=0.d0
            DO j=1,NITB
              DO k=1,NITB
                SUM1=SUM1+PROPQ(i,k,j+1,1,nq)*GUQ(j,k,nq)
              ENDDO !k
            ENDDO !j
            SUM2=0.d0
            DO k=1,NITB
              SUM2=SUM2+PROPQ(i,k,1,1,nq)*GCHQ(k,nq)
            ENDDO !k
            EU=EU+(SUM1-SUM2)*dUdX(i)
          ENDDO !i
          EU=EU*D

          AA=EU+4.0d0*F(1,1)*(UQ(-1, 0, 0)+UQ( 1, 0, 0))
     '         +4.0d0*F(2,2)*(UQ( 0,-1, 0)+UQ( 0, 1, 0))
     '         +2.0d0*F(1,2)*(UQ(-1,-1, 0)+UQ( 1, 1, 0)
     '                       -UQ( 1,-1, 0)-UQ(-1, 1, 0))

          IF(TYPE(1:6).EQ.'STATIC') THEN
            YQ(nq,2,na)=(YQ(nq,5,na)+AA)/(8.0d0*(F(1,1)+F(2,2)))

          ELSE IF(TYPE(1:7).EQ.'DYNAMIC') THEN
C           (YQ(6) holds previous timestep solution)
            YQ(nq,2,na)=(AA+YQ(nq,6,na)/DT+YQ(nq,5,na))
     '        /(8.0d0*(F(1,1)+F(2,2))+1.d0/DT)

          ELSE IF(TYPE(1:6).EQ.'ACTIVE') THEN
            U=UQ(0,0,0)
            BB=8.0d0*(F(1,1)+F(2,2))+1.d0/DT - 0.15d0
            CC=AA+YQ(nq,6,na)/DT+YQ(nq,5,na)
            RR=U*(U*(U-1.15d0)-BB)+CC  !residual
            dR=U*(3.d0*U-2.3d0)-BB     !derivative
            IF(DABS(dR).GT.1.d-8) THEN
              YQ(nq,2,na)=U-RR/DR      !Newton step
            ELSE
              YQ(nq,2,na)=U
            ENDIF
          ENDIF !type

        ENDIF !isotropic/anisotropic

        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_SOLVE_2)
          WRITE(OP_STRING,'(/'' Solution on coarsest grid: YQ('',I6,'
     '      //''',2,'',I1,'')='',D12.4)') nq,na,YQ(nq,2,na)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_SOLVE_2)
        ENDIF

      ELSE !solve coupled system
      ENDIF

      CALL EXITS('MG_SOLVE')
      RETURN
 9999 CALL ERRORS('MG_SOLVE',ERROR)
      CALL EXITS('MG_SOLVE')
      RETURN 1
      END


