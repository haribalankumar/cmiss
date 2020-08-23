      SUBROUTINE XPES37(NBH,NBJ,ne,NHE,NPNE,nr,nx,CE,CG,CGE,CP,ED,EM,
     &  ER,ES,PG,RG,WG,XE,XG,UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*)

C#### Subroutine: XPES37
C###  Description:
C###    XPES37 calculates element matrices for electrostatic and 
C###    magnetostatic finite element problems.
C**** Created by Martin Buist, September 2003.

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'

!     Parameter List
      INTEGER NBH(NHM,NCM),NBJ(NJM),ne,NHE,NPNE(NNM,NBFM),nr,nx
      REAL*8 CE(NMM),CG(NMM,NGM),CGE(NMM,NGM),CP(NMM,NPM),
     &  ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     &  ES(NHM*NSM,NHM*NSM),PG(NSM,NUM,NGM,NBM),RG(NGM),WG(NGM,NBM),
     &  XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
      LOGICAL UPDATE_MATRIX,UPDATE_VECTOR
!     Local Variables
      INTEGER mh,mhs,mi,ms,nb,NBMH,NBNH,nc,ng,NGTB,nh,nhs,NHST,ni,NITB,
     &  ns,NSTB,NU1(0:3),NXJ1(5)
      REAL*8 DXIX(3,3),GL(3,3),GU(3,3),PGM,PGMSI(3),PGMX,PGMX1,PGMX2,
     &  PGNSI(3),PGNX,PGNX1,PGNX2,PGX,RWG,SUM

      DATA NU1/1,2,4,7/
      DATA NXJ1/1,2,3,1,2/

      CALL ENTERS('XPES37',*9999)

      IF(ITYP3(nr,nx).EQ.1) THEN 
        !electrostatic
        CALL ASSERT(.FALSE.,'>>Electrostatics not implemented',ERROR,
     &    *9999)
      ELSE IF((ITYP3(nr,nx).EQ.2).OR.(ITYP3(nr,nx).EQ.3)) THEN 
        !magnetostatic
        CALL ASSERT(NHE.EQ.NJT,'>>Invalid no. of dependent variables',
     &    ERROR,*9999)
        CALL ASSERT(NJT.EQ.3,'>>Magnetostatic problems must be 3D',
     &    ERROR,*9999)
      ENDIF

      IF(UPDATE_MATRIX.OR.UPDATE_VECTOR) THEN
        nc=1 !ES is dependent variable
        nb=NBH(NH_LOC(1,nx),nc)
        NGTB=NGT(nb)
        NITB=NIT(nb)
        NSTB=NST(nb)
        NHST=NHE*NSTB
        CALL ASSERT(NST(nb).EQ.NNT(nb),
     &    '>>Only tested with linear elements',ERROR,*9999)

        !Transfer material paramaters to Gauss points
        CALL CPCG(1,nb,NPNE,nr,nx,CE,CG,CGE,CP,PG,ERROR,*9999)

        DO ms=1,NHST
          IF(UPDATE_MATRIX) THEN
            DO ns=1,NHST
              EM(ms,ns)=0.0d0
              ED(ms,ns)=0.0d0
              ES(ms,ns)=0.0d0
            ENDDO !ns
          ENDIF !update
          IF(UPDATE_VECTOR) ER(ms)=0.0d0
        ENDDO !ms

        !Loop over Gauss points
        DO ng=1,NGTB
          !Transfer element parameters to Gauss points
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
          CALL XGMG(0,NITB,nb,nr,DXIX,GL,GU,RG(ng),XG,ERROR,*9999)

          !Transformation Jacobian * Gauss weights
          IF(JTYP4.EQ.1)THEN !unsymmetric
            RWG=RG(ng)*WG(ng,nb)
          ELSE IF(JTYP4.EQ.2)THEN !cylindrical symmetry about x-axis
            RWG=RG(ng)*WG(ng,nb)*2.d0*PI*XG(2,1)
          ELSE IF(JTYP4.EQ.3)THEN !cylindrical symmetry about y-axis
            RWG=RG(ng)*WG(ng,nb)*2.d0*PI*XG(1,1)
          ELSE IF(JTYP4.EQ.4)THEN !spherical symmetry
            RWG=RG(ng)*WG(ng,nb)
          ENDIF

          !Calculate element stiffness matrix
          IF(UPDATE_MATRIX) THEN
            IF(ITYP3(nr,nx).EQ.1) THEN
              !Electrostatic

            ELSE IF(ITYP3(nr,nx).EQ.2) THEN 
              !Magnetostatic with no gauge
              mhs=0
              !Loop over dependent vars (Ax,Ay,Az)
              DO mh=1,NHE
                !Loop over element row dofs !(nn)
                DO ms=1,NSTB
                  mhs=mhs+1
                  nhs=0
                  !Loop over dependent vars (Ax,Ay,Az)
                  DO nh=1,NHE
                    !Loop over element col dofs (nn)
                    DO ns=1,NSTB
                      nhs=nhs+1
                      NBMH=NBH(mh,nc)
                      NBNH=NBH(nh,nc)
                      IF(mh.EQ.nh) THEN
                        PGMX1=PGX(NBMH,NXJ1(mh+1),ms,DXIX,
     &                    PG(1,1,ng,NBMH))
                        PGNX1=PGX(NBNH,NXJ1(mh+1),ns,DXIX,
     &                    PG(1,1,ng,NBNH))
                        PGMX2=PGX(NBMH,NXJ1(mh+2),ms,DXIX,
     &                    PG(1,1,ng,NBMH))
                        PGNX2=PGX(NBNH,NXJ1(mh+2),ns,DXIX,
     &                    PG(1,1,ng,NBNH))
                        ES(mhs,nhs)=ES(mhs,nhs)+(PGMX1*PGNX1+
     &                    PGMX2*PGNX2)*RWG*(1.0d0/CG(5,ng))
                      ELSE
                        PGMX=PGX(NBMH,mh,ms,DXIX,PG(1,1,ng,NBMH))
                        PGNX=PGX(NBNH,nh,ns,DXIX,PG(1,1,ng,NBNH))
                        ES(mhs,nhs)=ES(mhs,nhs)-PGMX*PGNX*RWG*
     &                    (1.0d0/CG(5,ng))
                      ENDIF    
                    ENDDO !ns
                  ENDDO !nh
                ENDDO !ms
              ENDDO !mh

            ELSE IF(ITYP3(nr,nx).EQ.3) THEN
              !Magnetostatic with a Coulomb gauge
              mhs=0
              !Loop over dependent vars (Ax,Ay,Az)
              DO mh=1,NHE
                !Loop over element row dofs (nn)
                DO ms=1,NSTB
                  mhs=mhs+1
                  nhs=(mh-1)*NSTB
                  !Loop over element col dofs (nn)
                  DO ns=1,NSTB
                    nhs=nhs+1                    
                    DO ni=1,NITB
                      PGMSI(ni)=PG(ms,NU1(ni),ng,nb)
                      PGNSI(ni)=PG(ns,NU1(ni),ng,nb)
                    ENDDO !ni
                    SUM =0.0d0
                    DO mi=1,NITB
                      DO ni=1,NITB
                        SUM=SUM+PGMSI(mi)*PGNSI(ni)*GU(mi,ni)
                      ENDDO !ni
                    ENDDO !mi
                    ES(mhs,nhs)=ES(mhs,nhs)+SUM*RWG*(1.0d0/CG(5,ng))
                  ENDDO !ns
                ENDDO !ms
              ENDDO !mh
            ENDIF !ityp3
          ENDIF !update_matrix

          !Calculate RHS vector
          IF(UPDATE_VECTOR)THEN
            mhs=0
            !Loop over dependent vars (Ax,Ay,Az)
            DO mh=1,NHE
              !Loop over element row dofs !(nn)
              DO ms=1,NSTB
                mhs=mhs+1
                PGM=PG(ms,1,ng,nb)
                IF(ITYP3(nr,nx).EQ.1) THEN
                  !Electrostatic
                ELSE IF(ITYP3(nr,nx).EQ.2) THEN
                  !Magnetostatic with no gauge
                  ER(mhs)=ER(mhs)+PGM*CG(mh,ng)*CG(4,ng)*RWG
                ELSE IF(ITYP3(nr,nx).EQ.3) THEN
                  !Magnetostatic with a Coulomb gauge 
                  ER(mhs)=ER(mhs)+PGM*CG(mh,ng)*CG(4,ng)*RWG
                ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF !update_vector
        ENDDO !ng

        IF(DOP) THEN
          WRITE(OP_STRING,'(''Element: '',I6)') ne
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO ms=1,NHST
            DO ns=1,NHST
              WRITE(OP_STRING,'(''ES('',I3,'','',I3,''): '',D12.6)')
     &          ms,ns,ES(ms,ns)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !ns
          ENDDO !ms
          DO ms=1,NHST
            WRITE(OP_STRING,'(''ER('',I3,''): '',D12.6)')
     &        ms,ER(ms)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !ms
        ENDIF
      ENDIF !update

      CALL EXITS('XPES37')
      RETURN
 9999 CALL ERRORS('XPES37',ERROR)
      CALL EXITS('XPES37')
      RETURN 1
      END



