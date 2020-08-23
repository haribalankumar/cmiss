      SUBROUTINE CALC_DTAR(maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,
     '  maqp2i,niq_old,niqV,nnq_min,NQXI,NRLIST,NSOL,NWQ,nx,NXQ,AQ,CQ,
     '  T,YQ,ADD,ERROR,*)

C#### Subroutine: CALC_DTAR
C###  Description:
C###    CALC_DTAR handles both the additions and removals from
C###    the list of active grid points when dynamic tracking of the
C###    active region is used.
C***  Created by Martin Buist, August 1998

      IMPLICIT NONE

      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter list
      INTEGER maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,maqp2i,niq_old,
     '  niqV,nnq_min,NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NSOL,NWQ(8,0:NQM),
     '  nx,NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM),T,YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
      LOGICAL ADD
!     Local variables
      INTEGER IJ,IK,mq,mqi,nbri,nii,nij,nik,NITB,niqSAC,nnq,nq,nr,nrr
      REAL*8 IAPP,VTHRESH
      LOGICAL ERROR_FLAG,SAC

      CALL ENTERS('CALC_DTAR',*9999)

      IF(ADD) THEN
        NITB=NQXI(0,1)
        DO nrr=1,NRLIST(0)
          nr=NRLIST(nrr)
          ERROR_FLAG=.FALSE.

C$OMP PARALLEL DO
C$OMP&  PRIVATE(nq,IAPP,SAC,niqSAC)
C$OMP&  SHARED(AQ,CQ,maqp1i,maqp1t0,maqp1t1,maqp2i,maqp2t0,maqp2t1,
C$OMP&  niq_old,niqV,nr,nx,NWQ,T,YQ,ERROR_FLAG,IODI,IOER)
          DO nq=NQR(1,nr),NQR(2,nr)
            IF(.NOT.ERROR_FLAG) THEN

C SGM 18Dec2000 use YQ(nq,niqV,1) instead of YQ(nq,niq_old,1) for adding
C           points as YQ(nq,niqV,1) hasn't yet been updated to new timestep
C           in MARCH8 and can be used as YQ(nq,niq_old,1).
C           YQ(nq,niq_old,1) is now being used for temporary storage
C           in MARCH8.

C             YQ(nq,niq_old,1)=YQ(nq,niqV,1)
              IF(NWQ(4,nq).EQ.0) THEN !Not already active or was active
                IAPP=0.0d0
                IF((AQ(maqp1i,nq).GT.ZERO_TOL).AND.(T.GE.AQ(maqp1t0,nq))
     '            .AND.(T.LT.AQ(maqp1t1,nq))) THEN
                  IAPP=AQ(maqp1i,nq)
                ELSE IF((AQ(maqp2i,nq).GT.ZERO_TOL).AND.(T.GE.
     '            AQ(maqp2t0,nq)).AND.(T.LT.AQ(maqp2t1,nq))) THEN
                  IAPP=AQ(maqp2i,nq)
                ENDIF
                SAC=.FALSE.
                IF(ITYP3(nr,nx).EQ.2) THEN      !FHN
                  IF(DABS(CQ(18,nq)).GT.ZERO_TOL) SAC=.TRUE.
                ELSEIF(ITYP3(nr,nx).EQ.3) THEN  !VCD
                  IF(DABS(CQ(13,nq)).GT.ZERO_TOL) SAC=.TRUE.
                ENDIF
                IF(SAC) THEN
                  CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqSAC,NIQ_SAC,
     '              ERROR,*100)
                  IAPP=IAPP-YQ(nq,niqSAC,1)
                ENDIF

                IF(IAPP.GT.ZERO_TOL) THEN
                  NWQ(4,nq)=1
                  NWQ(5,0)=NWQ(5,0)+1
                  NWQ(5,NWQ(5,0))=nq
                ENDIF
              ENDIF

              GOTO 102
 100          CONTINUE
C$OMP CRITICAL(CALC_DTAR_1)
              ERROR_FLAG=.TRUE.
              WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
              CALL WRITES(IOER,OP_STRING,ERROR,*101)
              WRITE(OP_STRING,'(/'' >>An error occurred - '
     '          //'results may be unreliable!'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101          CONTINUE
C$OMP END CRITICAL(CALC_DTAR_1)
 102          CONTINUE
            ENDIF !.NOT.ERROR_FLAG
          ENDDO !nq
C$OMP END PARALLEL DO
        ENDDO !nr

        IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
        IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D

C$OMP PARALLEL DO
C$OMP&PRIVATE(nbri,nik,nij,nii,nnq,nq,mq,mqi,VTHRESH)
C$OMP&SHARED(CQ,IJ,IK,niqV,NWQ,NXQ)
        DO nnq=1,NWQ(5,0)
          nq=NWQ(5,nnq)
          VTHRESH=CQ(9,nq)+1.0d0 !Potential shift by 1mV
C SGM 18Dec2000 use niqV instead of niq_old
          IF(NWQ(4,nq).EQ.1.AND.YQ(nq,niqV,1).GT.VTHRESH) THEN
            !Activate all neighbouring grid points
            DO nik=-IK,IK
              DO nij=-IJ,IJ
                DO nii=-1,1
                  IF(NXQ(nii,0,nq).LE.1) THEN
                    mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
                    IF(mq.GT.0.AND.NWQ(4,mq).EQ.0) THEN
                      !Not already active
                      NWQ(4,mq)=1
                      NWQ(5,0)=NWQ(5,0)+1
                      NWQ(5,NWQ(5,0))=mq
                    ENDIF
                  ELSE
                    DO nbri=1,NXQ(nii,0,nq) !branches in xi1 from nq
                      mqi=NXQ(nii,nbri,NXQ(nij*2,1,NXQ(nik*3,1,nq)))
                      IF(mqi.GT.0.AND.NWQ(4,mqi).EQ.0) THEN
                        !if point exists & not already active
                        NWQ(4,mqi)=1
                        NWQ(5,0)=NWQ(5,0)+1
                        NWQ(5,NWQ(5,0))=mqi
                        IF(NWQ(6,mqi).GT.0) THEN
                          NWQ(4,NWQ(6,mqi))=1
                          NWQ(5,0)=NWQ(5,0)+1
                          NWQ(5,NWQ(5,0))=NWQ(6,mqi)
                        ENDIF !coupled nodes
                      ENDIF
                    ENDDO !nbri
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            NWQ(4,nq)=2
          ENDIF
        ENDDO
C$OMP END PARALLEL DO

      ELSE ! Removal from DTAR
C$OMP PARALLEL DO
C$OMP&PRIVATE(nnq,nq)
C$OMP&SHARED(niq_old,niqV,nnq_min,NWQ,YQ)
        DO nnq=nnq_min,NWQ(5,0)
          nq=NWQ(5,nnq)

          IF(DABS(YQ(nq,niq_old,1)-YQ(nq,niqV,1)).LT.
     '      DSQRT(LOOSE_TOL)) THEN
            IF(NWQ(4,nq).GT.0) NWQ(4,nq)=NWQ(4,nq)+1
          ENDIF

          !stationary for 100 time steps
          IF(NWQ(4,nq).GE.100) NWQ(4,nq)=-1
        ENDDO
C$OMP END PARALLEL DO

        NWQ(4,0)=NSOL !Total number of active points
! Determine the minimum index of NWQ(,) containing an active
! data point for following iteration
        nnq=nnq_min
        DO WHILE ((NWQ(4,NWQ(5,nnq)).LE.0).AND.(nnq.LE.NQT))
          nnq=nnq+1
        ENDDO
        nnq_min=nnq
      ENDIF

      CALL EXITS('CALC_DTAR')
      RETURN
 9999 CALL ERRORS('CALC_DTAR',ERROR)
      CALL EXITS('CALC_DTAR')
      RETURN 1
      END



