      SUBROUTINE IMGRID_TRUEGRID(NENQ,NLATNE,NLATNQ,NLATPNQ,
     &  NLQ,NQNE,NQNLAT,NQS,NWQ,NXQ,AQ,XQ,COMPUTE_NXQ,ERROR,*)

C#### Subroutine: IMGRID_TRUEGRID
C###  Description:
C###    IMGRID_TRUEGRID imports a finite difference grid from a file containing
C###    grid point location and geometry
C***  Greg Sands, July 2003

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER NENQ(0:8,NQM),NLATNE(NEQM+1),NLATNQ(NEQM*NQEM),
     &  NLATPNQ(NQM),NLQ(NQM),NQNE(NEQM,NQEM),NQNLAT(NEQM*NQEM),
     &  NQS(NEQM),NWQ(8,0:NQM,NAM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)
      LOGICAL COMPUTE_NXQ

!     Local Variables
      INTEGER i,ILIST(10),loc,maq(3),ne,nj,nlat,nq,nq1,nq2,NQTEMP,
     &  VMAP(4,8)
      REAL*8 A(3),B(3),C(3),NQNORM_TOT,RLIST(5)

!     Functions
      REAL*8 DISTANCE

      DATA VMAP/4,5,7,4, 6,3,8,6,  3,6,9,3,  5,4,10,5, 
     &          3,9,8,3, 4,7,10,4, 5,10,7,5, 6,8,9,6/


      CALL ENTERS('IMGRID_TRUEGRID',*9999)

      CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maq(1),MAQ_NORMAL_X,
     &  ERROR,*9999)
      CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maq(2),MAQ_NORMAL_Y,
     &  ERROR,*9999)
      CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maq(3),MAQ_NORMAL_Z,
     &  ERROR,*9999)
C NQT: Number of grid points
      DO NQTEMP=1,NQT
        READ(UNIT=IFILE,FMT=*) nq,RLIST
C nq: the numbering point; 
C RLIST: three dimension coordinates 
        XQ(1,nq)=RLIST(2)
        XQ(2,nq)=RLIST(3)
        XQ(3,nq)=RLIST(4)
        IF(COMPUTE_NXQ) THEN
          NXQ(0,0,nq,1)=1
          NXQ(0,1,nq,1)=nq
        ENDIF
        NLQ(nq)=1
        NLATPNQ(nq)=0
        DO nj=1,3
          AQ(maq(nj),nq)=0.0d0
        ENDDO
      ENDDO !nq

      nlat=1
      DO ne=1,NEQM
C***    TrueGrid file has nodes specified counter-clockwise in xi1;xi3;xi2
C***    => adjust to follow CMISS standard. (1,2,5,6,4,3,8,7)
        READ(UNIT=IFILE,FMT=*) ILIST

C***    Added by Greg Sands 
C        NQTEMP=ILIST(5)
C        ILIST(5)=ILIST(7)
C        ILIST(7)=ILIST(6)
C        ILIST(6)=ILIST(8)
C        ILIST(8)=NQTEMP
C        NQTEMP=ILIST(9)
C        ILIST(9)=ILIST(10)
C        ILIST(10)=NQTEMP

C***    Added by JIZ on Oct. 25, 2007 to swith numberings on each element.
C***    for example: 
C***    (1,2,5,6,4,3,8,7)    ---->     (5,6,1,2,8,7,4,3)        
C***    note that the index of ILIST is from 3 to 10 for numberings
        NQTEMP = ILIST(5)
        ILIST(5) = ILIST(3)
        ILIST(3) = NQTEMP
        NQTEMP = ILIST(10)
        ILIST(10) = ILIST(3)
        ILIST(3) = NQTEMP
        NQTEMP = ILIST(7)        
        ILIST(7) = ILIST(3)
        ILIST(3) = NQTEMP
        NQTEMP = ILIST(6)
        ILIST(6) = ILIST(4)
        ILIST(4) = NQTEMP
        NQTEMP = ILIST(9)
        ILIST(9) = ILIST(4)
        ILIST(4) = NQTEMP
        NQTEMP = ILIST(8)
        ILIST(8) = ILIST(4)
        ILIST(4) = NQTEMP
       
        NLATNE(ne)=nlat
        NQS(ne)=1
        DO NQTEMP=1,8
          nq=ILIST(NQTEMP+2)
          NQNLAT(nlat)=nq
          NQNE(ne,NQTEMP)=nq
          IF(NLATPNQ(nq).EQ.0) THEN
            NLATPNQ(nq)=nlat
            NLATNQ(nlat)=0
            NENQ(1,nq)=ne
          ELSE
            loc=NLATPNQ(nq)
            DO WHILE (NLATNQ(loc).NE.0)
              loc=NLATNQ(loc)
            ENDDO
            NLATNQ(loc)=nlat
            NLATNQ(nlat)=0
          ENDIF
          nlat=nlat+1
          DO i=1,3
            nq1=ILIST(VMAP(i,NQTEMP))
            nq2=ILIST(VMAP(i+1,NQTEMP))
            IF(DISTANCE(3,XQ(1,nq),XQ(1,nq1)).GT.ZERO_TOL.AND.
     &         DISTANCE(3,XQ(1,nq),XQ(1,nq2)).GT.ZERO_TOL) THEN
              DO nj=1,3
                A(nj)=XQ(nj,nq1)-XQ(nj,nq)
                B(nj)=XQ(nj,nq2)-XQ(nj,nq)
              ENDDO
              CALL CROSS(A,B,C)
              CALL NORMALISE(3,C,ERROR,*9999)
              DO nj=1,3
                AQ(maq(nj),nq)=AQ(maq(nj),nq)+C(nj)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        IF(COMPUTE_NXQ) THEN
          IF(NXQ( 1,0,ILIST( 3),1).EQ.0) THEN
            NXQ( 1,0,ILIST( 3),1)=1
            NXQ( 1,1,ILIST( 3),1)=ILIST( 4)
          ENDIF
          IF(NXQ(-1,0,ILIST( 4),1).EQ.0) THEN
            NXQ(-1,0,ILIST( 4),1)=1
            NXQ(-1,1,ILIST( 4),1)=ILIST( 3)
          ENDIF
          IF(NXQ( 1,0,ILIST( 5),1).EQ.0) THEN
            NXQ( 1,0,ILIST( 5),1)=1
            NXQ( 1,1,ILIST( 5),1)=ILIST( 6)
          ENDIF
          IF(NXQ(-1,0,ILIST( 6),1).EQ.0) THEN
            NXQ(-1,0,ILIST( 6),1)=1
            NXQ(-1,1,ILIST( 6),1)=ILIST( 5)
          ENDIF
          IF(NXQ( 1,0,ILIST( 7),1).EQ.0) THEN
            NXQ( 1,0,ILIST( 7),1)=1
            NXQ( 1,1,ILIST( 7),1)=ILIST( 8)
          ENDIF
          IF(NXQ(-1,0,ILIST( 8),1).EQ.0) THEN
            NXQ(-1,0,ILIST( 8),1)=1
            NXQ(-1,1,ILIST( 8),1)=ILIST( 7)
          ENDIF
          IF(NXQ( 1,0,ILIST( 9),1).EQ.0) THEN
            NXQ( 1,0,ILIST( 9),1)=1
            NXQ( 1,1,ILIST( 9),1)=ILIST(10)
          ENDIF
          IF(NXQ(-1,0,ILIST(10),1).EQ.0) THEN
            NXQ(-1,0,ILIST(10),1)=1
            NXQ(-1,1,ILIST(10),1)=ILIST(9)
          ENDIF

          IF(NXQ( 2,0,ILIST( 3),1).EQ.0) THEN
            NXQ( 2,0,ILIST( 3),1)=1
            NXQ( 2,1,ILIST( 3),1)=ILIST( 5)
          ENDIF
          IF(NXQ(-2,0,ILIST( 5),1).EQ.0) THEN
            NXQ(-2,0,ILIST( 5),1)=1
            NXQ(-2,1,ILIST( 5),1)=ILIST( 3)
          ENDIF
          IF(NXQ( 2,0,ILIST( 4),1).EQ.0) THEN
            NXQ( 2,0,ILIST( 4),1)=1
            NXQ( 2,1,ILIST( 4),1)=ILIST( 6)
          ENDIF
          IF(NXQ(-2,0,ILIST( 6),1).EQ.0) THEN
            NXQ(-2,0,ILIST( 6),1)=1
            NXQ(-2,1,ILIST( 6),1)=ILIST( 4)
          ENDIF
          IF(NXQ( 2,0,ILIST( 8),1).EQ.0) THEN
            NXQ( 2,0,ILIST( 8),1)=1
            NXQ( 2,1,ILIST( 8),1)=ILIST(10)
          ENDIF
          IF(NXQ(-2,0,ILIST(10),1).EQ.0) THEN
            NXQ(-2,0,ILIST(10),1)=1
            NXQ(-2,1,ILIST(10),1)=ILIST( 8)
          ENDIF
          IF(NXQ( 2,0,ILIST( 7),1).EQ.0) THEN
            NXQ( 2,0,ILIST( 7),1)=1
            NXQ( 2,1,ILIST( 7),1)=ILIST( 9)
          ENDIF
          IF(NXQ(-2,0,ILIST( 9),1).EQ.0) THEN
            NXQ(-2,0,ILIST( 9),1)=1
            NXQ(-2,1,ILIST( 9),1)=ILIST( 7)
          ENDIF

          IF(NXQ( 3,0,ILIST( 3),1).EQ.0) THEN
            NXQ( 3,0,ILIST( 3),1)=1
            NXQ( 3,1,ILIST( 3),1)=ILIST( 7)
          ENDIF
          IF(NXQ(-3,0,ILIST( 7),1).EQ.0) THEN
            NXQ(-3,0,ILIST( 7),1)=1
            NXQ(-3,1,ILIST( 7),1)=ILIST( 3)
          ENDIF
          IF(NXQ( 3,0,ILIST( 4),1).EQ.0) THEN
            NXQ( 3,0,ILIST( 4),1)=1
            NXQ( 3,1,ILIST( 4),1)=ILIST( 8)
          ENDIF
          IF(NXQ(-3,0,ILIST( 8),1).EQ.0) THEN
            NXQ(-3,0,ILIST( 8),1)=1
            NXQ(-3,1,ILIST( 8),1)=ILIST( 4)
          ENDIF
          IF(NXQ( 3,0,ILIST( 5),1).EQ.0) THEN
            NXQ( 3,0,ILIST( 5),1)=1
            NXQ( 3,1,ILIST( 5),1)=ILIST( 9)
          ENDIF
          IF(NXQ(-3,0,ILIST( 9),1).EQ.0) THEN
            NXQ(-3,0,ILIST( 9),1)=1
            NXQ(-3,1,ILIST( 9),1)=ILIST( 5)
          ENDIF
          IF(NXQ( 3,0,ILIST( 6),1).EQ.0) THEN
            NXQ( 3,0,ILIST( 6),1)=1
            NXQ( 3,1,ILIST( 6),1)=ILIST(10)
          ENDIF
          IF(NXQ(-3,0,ILIST(10),1).EQ.0) THEN
            NXQ(-3,0,ILIST(10),1)=1
            NXQ(-3,1,ILIST(10),1)=ILIST( 6)
          ENDIF
        ENDIF !compute_nxq
      ENDDO !ne
      NLATNE(NEQM+1)=nlat

      IF(COMPUTE_NXQ) THEN
        DO nq=1,NQT
          nq1=nq
          DO i=-3,3
            IF(NXQ(i,0,nq,1).EQ.0) THEN
              nq1=NXQ(-i,1,nq1,1)
            ENDIF
          ENDDO
          IF(nq1.EQ.nq) THEN
            NWQ(1,nq,1)=0
            NWQ(2,nq,1)=0
          ELSE
            NWQ(1,nq,1)=nq1
            DO i=-3,3
              IF(NXQ(i,0,nq,1).EQ.0) THEN
                nq1=NXQ(-i,1,nq1,1)
              ENDIF
            ENDDO
            NWQ(2,nq,1)=nq1
          ENDIF
        ENDDO !nq
      ELSE
        DO nq=1,NQT
          NQNORM_TOT=0.0d0
          DO nj=1,3
            NQNORM_TOT=NQNORM_TOT+AQ(maq(nj),nq)**2
          ENDDO
          NQNORM_TOT=DSQRT(NQNORM_TOT)
          IF(NQNORM_TOT.GT.1.0d-6) THEN
            NWQ(1,nq,1)=1
            DO nj=1,3
              AQ(maq(nj),nq)=AQ(maq(nj),nq)/NQNORM_TOT
            ENDDO
          ELSE
            NWQ(1,nq,1)=0
          ENDIF
C          WRITE(*,'(''nq: '',I6,'' nwq: '',I2,'' norm: '',3F8.4,
C     &      '' |norm|: '',e12.6)')
C     &      nq,NWQ(1,nq,1),(AQ(maq(nj),nq),nj=1,3),NQNORM_TOT
        ENDDO
      ENDIF !compute_nxq
      CALL EXITS('IMGRID_TRUEGRID')
      RETURN
 9999 CALL ERRORS('IMGRID_TRUEGRID',ERROR)
      CALL EXITS('IMGRID_TRUEGRID')
      RETURN 1
      END


