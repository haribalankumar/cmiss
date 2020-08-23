      SUBROUTINE MG_RELAX(N_relax,na,niq,NITB,NAQ,NLQ,NWQ,NXQ,
     '  ADAPTIVE,ISOTROPIC,TYPE,DT,GCHQ,GUQ,PROPQ,YQ,ERROR,*)

C#### Subroutine: MG_RELAX
C###  Description:
C###    MG_RELAX performs N_relax relaxations of current solution
C###    YQ(nq,niq,na) on a given multigrid level.

C**** RHS=Tau=YQ(nq,5,na)
C**** DD makes adjustment between element Xi coords and Xi coords
C**** used in biquadratic patch.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER na,niq,N_relax,NITB,NAQ(NQM),NLQ(NQM),NWQ(8,0:NQM),
     '  NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 DT,GCHQ(3,NQM),GUQ(3,3,NQM),PROPQ(3,3,4,2,NQM),
     '  YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*),TYPE*(*)
      LOGICAL ADAPTIVE,ISOTROPIC
!     Local Variables
      INTEGER i,ij,ik,j,k,mq,nii,nij,nik,nq,Relax
      REAL*8 AA,BB,CC,D,DD,dR,dUdX(3),EU,F(3,3),
     '  GCHQ1,GCHQ2,GUQ11,GUQ22,GUQ12,RR,SUM1,SUM2,
     '  U,U1,U2,UQ(-1:1,-1:1,-1:1)

      CALL ENTERS('MG_RELAX',*9999)

C     DD=(0.5d0*dble(2**na)/dble(2**nat(1)))**2 !grid spacing
      D=2.0D0/DBLE(2**na)
      DD=D*D
      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_RELAX_1)
        WRITE(OP_STRING,'('' na='',I2,'' D='',D13.4)') na,D
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_RELAX_1)
      ENDIF

      DO i=-1,1
        DO j=-1,1
          DO k=-1,1
            UQ(i,j,k)=0.0d0
          ENDDO
        ENDDO
      ENDDO

      DO Relax=1,N_relax
        DO nq=1,NQT !Loop over grid points
          IF(NAQ(nq).EQ.0         !nq is active in grid na
     '      .AND.(.NOT.ADAPTIVE
     '            .OR.(ADAPTIVE.AND.NLQ(nq).EQ.na))) THEN

            IF(NWQ(1,nq).EQ.0) THEN !interior g.p.
              ik=MAX(0,NITB-2) !zero for 1,2D, one for 3D
              ij=MIN(NITB-1,1) !zero for 1D, one for 2,3D
              DO nik=-ik,ik    !UQ has coeffs of local quadratic
                DO nij=-ij,ij  !element about nq
                  DO nii=-1,1
                    mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
                    IF(mq.GT.0) THEN !mq is adjacent g.p.
                      UQ(nii,nij,nik)=YQ(mq,niq,na) !u at mq
                    ELSE
                      write(*,'('' Error: outside grid in MG_RELAX'')')
                    ENDIF
                  ENDDO !nii
                ENDDO !nij
              ENDDO !nik

C           Update solution at nq
              IF(ISOTROPIC) THEN !2D isotropy
                GUQ11=GUQ(1,1,nq)*DD
                GUQ22=GUQ(2,2,nq)*DD
                GUQ12=GUQ(1,2,nq)*DD
                GCHQ1=GCHQ( 1,nq)*D
                GCHQ2=GCHQ( 2,nq)*D
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' nq='',I3,'' GUQ11='',D13.4,
     '              '' GUQ22='',D13.4,'' GUQ12='',D13.4,
     '              '' GCHQ1='',D13.4,'' GCHQ2='',D13.4)')
     '              nq,GUQ11,GUQ22,GUQ12,GCHQ1,GCHQ2
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                AA=4.0D0*GUQ11*(UQ(-1, 0, 0)+UQ( 1, 0, 0))
     '            +4.0D0*GUQ22*(UQ( 0,-1, 0)+UQ( 0, 1, 0))
     '            +2.0D0*GUQ12*(UQ(-1,-1, 0)+UQ( 1, 1, 0)
     '                         -UQ( 1,-1, 0)-UQ(-1, 1, 0))
     '                  +GCHQ1*(UQ(-1, 0, 0)-UQ( 1, 0, 0))
     '                  +GCHQ2*(UQ( 0,-1, 0)-UQ( 0, 1, 0))

                YQ(nq,niq,na)=(YQ(nq,5,na)+AA)/(8.0d0*(GUQ11+GUQ22))
c               YQ(nq,niq,na)=0.25D0*(DD*YQ(nq,5,na)
c    '                               +UQ(-1,0,0)+UQ(1,0,0)
c    '                               +UQ(0,-1,0)+UQ(0,1,0))

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
     '               +4.0d0*F(2,2)*(UQ( 0,-1, 0)+UQ( 0, 1, 0))
     '               +2.0d0*F(1,2)*(UQ(-1,-1, 0)+UQ( 1, 1, 0)
     '                             -UQ( 1,-1, 0)-UQ(-1, 1, 0))

                IF(TYPE(1:6).EQ.'STATIC') THEN
                  YQ(nq,niq,na)=(YQ(nq,5,na)+AA)
     '              /(8.0d0*(F(1,1)+F(2,2)))

                ELSE IF(TYPE(1:7).EQ.'DYNAMIC') THEN
C                 (YQ(6) holds previous timestep solution)
                  YQ(nq,niq,na)=(AA+YQ(nq,6,na)/DT+YQ(nq,5,na))
     '              /(8.0d0*(F(1,1)+F(2,2))+1.d0/DT)

                ELSE IF(TYPE(1:6).EQ.'ACTIVE') THEN
                  U=UQ(0,0,0)
                  BB=8.0d0*(F(1,1)+F(2,2))+1.d0/DT - 0.15d0
                  CC=AA+YQ(nq,6,na)/DT+YQ(nq,5,na)
                  RR=U*(U*(U-1.15d0)-BB)+CC  !residual
                  dR=U*(3.d0*U-2.3d0)-BB     !derivative
                  IF(DABS(dR).GT.1.d-8) THEN
                    YQ(nq,niq,na)=U-RR/DR    !Newton step
                  ELSE
                    YQ(nq,niq,na)=U
                  ENDIF
                ENDIF !type

              ENDIF !isotropic/anisotropic

              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_RELAX_2)
                WRITE(OP_STRING,
     '            '('' Solution YQ('',I6,'','',I1,'','',I1,'')='','
     '            //'D12.4,'' with RHS YQ('',I6,'',5,'',I1,'')='','
     '            //'D12.4,'' at Relax='',I2)')
     '            nq,niq,na,YQ(nq,niq,na),nq,na,YQ(nq,5,na),Relax
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_RELAX_2)
              ENDIF

            ELSE !Boundary point
              IF(NWQ(5,nq).EQ.2) THEN !Neumann b.c.
                IF(NWQ(2,nq).GT.0) THEN !2 interior pts defined
                  U1=YQ(NWQ(1,nq),niq,na)   !Pick up values for g.p.
                  U2=YQ(NWQ(2,nq),niq,na)   !as defined in IPGRID
                  YQ(nq,niq,na)=(4.0D0*U1-U2)/3.0D0
                  IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_RELAX_3)
                   WRITE(OP_STRING,'('' nq:'',I5,'' u1,u2,u:'',3D12.5)')
     '                nq,U1,U2,YQ(nq,niq,na)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_RELAX_3)
                  ENDIF
                ELSE                       !only 1 interior pt defined
                  U1=YQ(NWQ(1,nq),niq,na)
                  YQ(nq,niq,na)=U1
                  IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_RELAX_4)
                    WRITE(OP_STRING,'('' nq:'',I5,'' u1,u:'',2D12.5)')
     '                nq,U1,YQ(nq,niq,na)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_RELAX_4)
                  ENDIF
                ENDIF
              ENDIF !Neumann bc
            ENDIF !interior/bdry point

          ENDIF !active grid point

        ENDDO !nq
      ENDDO !Relax

      CALL EXITS('MG_RELAX')
      RETURN
 9999 CALL ERRORS('MG_RELAX',ERROR)
      CALL EXITS('MG_RELAX')
      RETURN 1
      END


