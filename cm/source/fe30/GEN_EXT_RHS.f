      SUBROUTINE GEN_EXT_RHS(niqV,niqBNDRY,NQGP,NQGP_PIVOT,
     '  NQGW,GM,NRLIST,NTIME_INTERP,NTIME_POINTS,
     '  NWQ,nx_ext,nx_trans,PROPQ,RHS,T,TIME_VALUES,YQ,
     '  CQ,FIXQ,ERROR,*)

C#### Subroutine: GEN_EXT_RHS
C###  Description:
C###    GEN_EXT_RHS is used to generate the right hand side vectors
C###    at each time step for the extracellular domain in the implicit
C###    solution of grid activation problems.
C***  Created by Martin Buist, July 1997

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'

!     Parameter list
      INTEGER niqV,niqBNDRY,NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),
     '  NRLIST(0:NRM),NTIME_INTERP(NTIMEVARSM),NTIME_POINTS(NTIMEVARSM),
     '  NWQ(8,0:NQM),nx_ext,nx_trans
      REAL*8 EVTIME_FCN,NQGW(NQGM,NQM),GM(NZ_GM_M),PROPQ(3,3,4,2,NQM),
     '  RHS(NQM),T,TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM),
     '  YQ(NYQM,NIQM,NAM,NXM),CQ(NMM,NQM)
      LOGICAL FIXQ(NYQM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ERRORCODE,i,nq,nr,nrr

      CALL ENTERS('GEN_EXT_RHS',*9999)

      ERRORCODE=0
      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr)
C$OMP   PARALLEL DO
C$OMP&  PRIVATE(ERRORCODE,i,nq)
C$OMP&  SHARED(CQ,FIXQ,niqBNDRY,niqV,nr,NQGP,NQGP_PIVOT,NQGW,
C$OMP&    NTIME_INTERP,NTIME_POINTS,NWQ,nx_ext,nx_trans,PROPQ,RHS,T,
C$OMP&    TIME_VALUES,YQ)
        DO nq=NQR(1,nr),NQR(2,nr)
C SGM 12Jan01 Grid-based FE implementation
C MLT 29Nov02 Grid FV also
          IF(ITYP4(NRLIST(1),nx_trans).EQ.6.OR.
     '       ITYP4(NRLIST(1),nx_trans).EQ.7)THEN !Grid-based FE
C there is no internal/external differentiation in Grid-based FE
            IF(.NOT.(FIXQ(nq,1).OR.FIXQ(nq,3)))THEN
              RHS(nq)=0.0d0
              DO i=1,NQGP(0,nq)
                IF(NQGP(i,nq).GT.0) RHS(nq)=RHS(nq)-
     '            (NQGW(NQGP_PIVOT(i,nq),nq)*
     '            YQ(NQGP(i,nq),niqV,1,nx_trans))
              ENDDO !i
C MLT 25Feb03 Secondary diffusion contribution for grid based FV
              IF (ITYP4(NRLIST(1),nx_trans).EQ.7) THEN
                DO i=1,3 ! loop over the coordinate directions
                  RHS(nq) = RHS(nq) + PROPQ(2,i,3,2,nq)+ 
     '                                PROPQ(2,i,4,2,nq)
                ENDDO
              ENDIF

              IF(NWQ(8,nq).GT.0) THEN
C!!! Does the following equation have to be multiplied by the mass matrix?
C                RHS(nq)=RHS(nq)+EVTIME_FCN(NTIME_INTERP,NTIME_POINTS,
C     '            NWQ(8,nq),T,TIME_VALUES)
C MLT 16Nov02 The extracellular stimulation contribution to RHS should be 
C             multiplied by the mass matrix GM. Note that the units of the 
C             extracellular stimulation are uA/uF, i.e. extracellular
C             applied current density per membrane capacitance per unit volume
                DO i=1,NQGP(0,nq)
                  RHS(nq)=RHS(nq)-GM(NQGP_PIVOT(i,nq)+NQGM*(nq-1))*
     '                            EVTIME_FCN(NTIME_INTERP,NTIME_POINTS,
     '                            NWQ(8,nq),T,TIME_VALUES)
                ENDDO !i
              ENDIF
            ENDIF
            IF(FIXQ(nq,1)) THEN !dependent variable boundary condition
              RHS(nq)=YQ(nq,niqBNDRY,1,nx_ext)
              IF(NWQ(8,nq).GT.0) THEN
                RHS(nq)=EVTIME_FCN(NTIME_INTERP,NTIME_POINTS,
     '            NWQ(8,nq),T,TIME_VALUES)
              ENDIF
            ELSEIF(FIXQ(nq,2)) THEN !flux boundary condition
              RHS(nq)=RHS(nq)+YQ(nq,niqBNDRY,1,nx_ext)
C KAT don't want to include transmembrane potential term twice.
C              DO i=1,NQGP(0,nq)
C                IF(NQGP(i,nq).GT.0) RHS(nq)=RHS(nq)-
C     '            (NQGW(NQGP_PIVOT(i,nq),nq)*
C     '            YQ(NQGP(i,nq),niqV,1,nx_trans))
C              ENDDO !i
            ELSEIF(FIXQ(nq,3)) THEN !analytic boundary condition
              RHS(nq)=-CQ(3,nq)/(CQ(3,nq)+CQ(6,nq))
              RHS(nq)=RHS(nq)*(YQ(nq,niqV,1,nx_trans)-CQ(9,nq))
              IF(NWQ(8,nq).GT.0) THEN
                RHS(nq)=EVTIME_FCN(NTIME_INTERP,NTIME_POINTS,
     '            NWQ(8,nq),T,TIME_VALUES)
              ENDIF
            ENDIF
          ELSE ! not Grid-based FE
            IF(USE_LAT.EQ.0) THEN
              IF(NWQ(1,nq).EQ.0) THEN !internal
                RHS(nq)=0.0d0
                DO i=1,NQGP(0,nq)
                  IF(NQGP(i,nq).GT.0) RHS(nq)=RHS(nq)-
     '              (NQGW(NQGP_PIVOT(i,nq),nq)*
     '              YQ(NQGP(i,nq),niqV,1,nx_trans))
                ENDDO !i
                IF(NWQ(8,nq).GT.0) THEN
                  RHS(nq)=RHS(nq)+EVTIME_FCN(NTIME_INTERP,NTIME_POINTS,
     '              NWQ(8,nq),T,TIME_VALUES)
                ENDIF
              ELSE !boundary grid point for non Grid-based problems
                IF(FIXQ(nq,1)) THEN
                  RHS(nq)=YQ(nq,niqBNDRY,1,nx_ext)
                  IF(NWQ(8,nq).GT.0) THEN
                    RHS(nq)=EVTIME_FCN(NTIME_INTERP,NTIME_POINTS,
     '                NWQ(8,nq),T,TIME_VALUES)
                  ENDIF
                ELSEIF(FIXQ(nq,2)) THEN
                  RHS(nq)=YQ(nq,niqBNDRY,1,nx_ext)
                  IF(NWQ(8,nq).GT.0) THEN
                    RHS(nq)=RHS(nq)+EVTIME_FCN(NTIME_INTERP,
     '                NTIME_POINTS,NWQ(8,nq),T,TIME_VALUES)
                  ENDIF
                ELSEIF(FIXQ(nq,3)) THEN
                  RHS(nq)=-PROPQ(1,1,1,1,nq)/(PROPQ(1,1,1,1,nq)+
     '              PROPQ(1,1,1,2,nq))
                  RHS(nq)=RHS(nq)*(YQ(nq,niqV,1,nx_trans)-CQ(9,nq))
                  IF(NWQ(8,nq).GT.0) THEN
                    RHS(nq)=EVTIME_FCN(NTIME_INTERP,NTIME_POINTS,
     '                NWQ(8,nq),T,TIME_VALUES)
                  ENDIF
                ELSE
                  ERRORCODE=nq
                ENDIF
              ENDIF !internal
            ELSE !lattice method
              IF(NWQ(1,nq).EQ.0) THEN !internal
                RHS(nq)=0.0d0
                DO i=1,NQGP(0,nq)
                  IF(NQGP(i,nq).GT.0) RHS(nq)=RHS(nq)-
     '              (NQGW(i,nq)*
     '              YQ(NQGP(i,nq),niqV,1,nx_trans))
                ENDDO !i
                IF(NWQ(8,nq).GT.0) THEN
                  RHS(nq)=RHS(nq)+EVTIME_FCN(NTIME_INTERP,NTIME_POINTS,
     '              NWQ(8,nq),T,TIME_VALUES)
                ENDIF
              ELSE !boundary grid point for non Grid-based problems
                IF(FIXQ(nq,1)) THEN
                  RHS(nq)=YQ(nq,niqBNDRY,1,nx_ext)
                  IF(NWQ(8,nq).GT.0) THEN
                    RHS(nq)=EVTIME_FCN(NTIME_INTERP,NTIME_POINTS,
     '                NWQ(8,nq),T,TIME_VALUES)
                  ENDIF
                ELSEIF(FIXQ(nq,2)) THEN
                  RHS(nq)=YQ(nq,niqBNDRY,1,nx_ext)
                  IF(NWQ(8,nq).GT.0) THEN
                    RHS(nq)=RHS(nq)+EVTIME_FCN(NTIME_INTERP,
     '                NTIME_POINTS,NWQ(8,nq),T,TIME_VALUES)
                  ENDIF
                ELSEIF(FIXQ(nq,3)) THEN
                  RHS(nq)=-PROPQ(1,1,1,1,nq)/(PROPQ(1,1,1,1,nq)+
     '              PROPQ(1,1,1,2,nq))
                  RHS(nq)=RHS(nq)*(YQ(nq,niqV,1,nx_trans)-CQ(9,nq))
                  IF(NWQ(8,nq).GT.0) THEN
                    RHS(nq)=EVTIME_FCN(NTIME_INTERP,NTIME_POINTS,
     '                NWQ(8,nq),T,TIME_VALUES)
                  ENDIF
                ELSE
                  ERRORCODE=nq
                ENDIF
              ENDIF !internal              
            ENDIF
          ENDIF !grid-based FE
        ENDDO !nq
C$OMP   END PARALLEL DO
      ENDDO !nr

      IF(ERRORCODE.NE.0) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(GEN_EXT_RHS_1)
        WRITE(OP_STRING,'(''>>Error: No boundary condition applied at'
     '    //' point '',I8)') nq
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP   END CRITICAL(GEN_EXT_RHS_1)
      ENDIF

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(GEN_EXT_RHS_2)
        WRITE(OP_STRING,'(''EXT_RHS:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nq=1,NQT
          WRITE(OP_STRING,'(I8,F12.6)') nq,RHS(nq)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nq
CC$OMP   END CRITICAL(GEN_EXT_RHS_2)
      ENDIF

      CALL EXITS('GEN_EXT_RHS')
      RETURN
 9999 CALL ERRORS('GEN_EXT_RHS',ERROR)
      CALL EXITS('GEN_EXT_RHS')
      RETURN 1
      END



