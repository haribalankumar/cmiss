      SUBROUTINE DIRECTSOLUTION(ITER_USER,NBJ,NEELEM,NORD,NPLIST,NPNE,
     &  NVJE,NVJP,NXI,BBM,CE,CW,dPl,dt,ERR_USER,
     &  InitialVolume,MeanVolume,ppl_current,Pmus,time,ttime,
     &  T_interval,undef,XAB,XP,DIAG_OP,FAST,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn' 
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung00.cmn' 
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn' 

      INTEGER ITER_USER,NBJ(NJM,NEM),NEELEM(0:NE_R_M),NORD(5,NE_R_M),
     &  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),CW,dPl,dt,ERR_USER,
     &  InitialVolume,MeanVolume,Pmus,time,ttime,T_interval,undef,
     &  XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),ppl_current
      LOGICAL DIAG_OP,FAST,INLIST
      CHARACTER ERROR*(*)

      INTEGER nb,ne,ne2,noelem,noelem2,Nend,np,np1,np2,Ntmp,nv,
     &  count_dec,count_inc
      REAL*8 CWW,lambda,mean_PVR,totalC,total_resistance,vratio,aa,
     &  bb,cc,Pe,Ppl_ne,mean_Ptp,iter_step,ErrorEstimate,flow_diff,
     &  FlowSum,sum_change
      LOGICAL CONVERGED,RE_SOLVE

      CALL ENTERS('DIRECTSOLUTION',*9999)
C We have already calculated the flow into each terminal lumped
C parameter unit (assumed to be an acinus), so we can calculate flow
C throughout the rest of the tree simply by summation. After summing
C the flows we can use the resistance equation (P0-P1=R1*Q1) to update
C the pressures throughout the tree.

      NTB=0
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).EQ.0)THEN
           nb=NBJ(1,ne)
           np=NPNE(2,nb,ne)
           XAB(2,ne)=XP(1,1,nj_flow,np)
           NTB=NTB+1
        ENDIF
      ENDDO

      CONVERGED = .FALSE.
      RE_SOLVE=.FALSE.
      iter_step=0
      DO WHILE(.NOT.CONVERGED)
        IF(RE_SOLVE) iter_step=0
        iter_step=iter_step+1
        ErrorEstimate=0.d0
        FlowSum=0.d0
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          IF(NXI(1,0,ne).EQ.0)THEN
            nb=NBJ(1,ne)
            np1=NPNE(1,nb,ne)
            np2=NPNE(2,nb,ne)
            IF(INLIST(np2,NPLIST(1),NPLIST(0),Ntmp))THEN
              XP(2,1,nj_flow,np2)=0.d0
            ELSE
C.......Calculate the mean flow into the unit in the time step
              vratio=BBM(1,ne)/undef
              IF(RE_SOLVE)THEN
                XP(2,1,nj_flow,np2)=XAB(2,ne) !Qinit
                XAB(6,np2)=XAB(2,ne)
                XAB(7,np2)=XAB(2,ne)
                XAB(1,ne)=dPl/dt
                XAB(3,ne)=dPl/dt
              ENDIF
              CALL FLOW_AVERAGE(np1,np2,iter_step,BBM(1,ne),
     &          vratio,CE(1,ne),
     &          dPl,dt,XAB(1,ne),XAB(2,ne),XAB,XP,ERROR,*9999)
              flow_diff=XP(2,1,nj_flow,np2)-XP(1,1,nj_flow,np2)
              ErrorEstimate=ErrorEstimate+
     &          DABS(flow_diff)**2.d0
              FlowSum=FlowSum+XP(2,1,nj_flow,np2)**2.d0
            ENDIF
          ENDIF
        ENDDO

        ErrorEstimate=ErrorEstimate/(FlowSum*DBLE(NTB))

        IF(RE_SOLVE) RE_SOLVE=.FALSE.
        IF(iter_step.gt.1.AND.ErrorEstimate.LT.ERR_USER)THEN
          CONVERGED=.TRUE.
        ELSEIF(iter_step.gt.ITER_USER)THEN
          CONVERGED=.TRUE.
          WRITE(OP_STRING,'('' Warning: lower convergence '//
     &      'tolerance and time step - check values, Error='',D10.3)')
     &      ErrorEstimate
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        CALL FLOW_SUMMATION(NBJ,NEELEM,NPNE,NVJE,NVJP,NXI,XP,
     &    ERROR,*9999)
        !ENDIF

C.....Use the known resistances and flows to calculate nodal pressures      

c        CALL BRANCHRESISTANCE(NBJ,NEELEM,NPNE,NVJE,CE,
c     &    XP,ERROR,*9999)
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          np1=NPNE(1,1,ne)
          np2=NPNE(2,1,ne)
          XP(1,1,nj_pressure,np2)=XP(1,1,nj_pressure,np1)
     &      -CE(nm_R,ne)*XP(2,1,nj_flow,np2)
          DO nv=2,NVJP(nj_pressure,np2)
            XP(1,nv,nj_pressure,np2)=XP(1,nv,nj_pressure,np2)
          ENDDO
        ENDDO !noelem
        CALL UPDATEPRESSUREDT(iter_step,NBJ,NEELEM,NPNE,NXI,CE,dPl,
     &    dt,time,ttime,T_interval,XAB,XP,DIAG_OP,.TRUE.,ERROR,*9999)

       ENDDO
C.....Update pressures for next time step      

        XP(2,1,nj_pressure,1)=XP(1,1,nj_pressure,1) !store previous pressure
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          np1=NPNE(1,1,ne)
          np2=NPNE(2,1,ne)
          XP(2,1,nj_pressure,np2)=XP(1,1,nj_pressure,np2)
          DO nv=2,NVJP(nj_pressure,np2)
            XP(2,nv,nj_pressure,np2)=XP(1,nv,nj_pressure,np2)
          ENDDO
        ENDDO !noelem

      DO noelem=NEELEM(0),1,-1
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).gt.0)THEN
          total_resistance=0.d0
          DO noelem2=1,NXI(1,0,ne)
            ne2=NXI(1,noelem2,ne)
            total_resistance=total_resistance+1.d0/CE(nm_Rt,ne2)
          ENDDO
          CE(nm_Rt,ne)=CE(nm_Rt,ne)+1.d0/total_resistance
        ENDIF
      ENDDO !noelem

      CALL TISSUECOMPLIANCE(NEELEM,NPNE,NXI,BBM,CE,CW,dPl,
     &  undef,.FALSE.,ERROR,*9999)

      Nend=0
      Ppl_current = 0.d0
      mean_Ptp=0.d0
      mean_PVR=0.d0
      totalC=0.d0

      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).EQ.0)THEN
          nb=NBJ(1,ne)
          np=NPNE(2,nb,ne)
C update the volume of the lumped parameter unit
          BBM(1,ne)=BBM(1,ne)+dt*XP(1,1,nj_flow,np) !in mm^3
          Pe = CE(nm_Ppl,ne)
          Palv = XP(1,1,nj_pressure,np)
          mean_Ptp=mean_Ptp+Pe
          Ppl_current = Ppl_current - Pe + Palv
          Nend=Nend+1
          IF(XP(1,1,nj_flow,1).GT.0.d0)THEN !only store inspired volume
            BBM(2,ne)=BBM(2,ne)+dt*XP(1,1,nj_flow,np) !in mm^3, flow
          ENDIF
          mean_PVR=mean_PVR+XP(1,1,nj_pressure,1)-Palv
          totalC=totalC+CE(nm_C,ne)
        ENDIF
      ENDDO

      Ppl_current = Ppl_current/DBLE(Nend)
      mean_PVR=mean_PVR/DBLE(Nend)
      mean_Ptp=mean_Ptp/DBLE(Nend)

      IF(.NOT.FAST)THEN
        CALL VOLUMEOFMESH(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,MeanVolume,
     &   BBM,CE,XP,ERROR,*9999)
      ENDIF


      WRITE(OP_STRING,'(F7.3,11(F14.3))')
     &  time, !time through breath (s)
     &  XP(1,1,nj_flow,1)/1.d3, !flow at the inlet (mL/s)
     &  (CE(nm_vol_bel,1)-InitialVolume)/1.d3, !current tidal volume (mL)
     &  CE(nm_Rt,1)*1.d6/98.0665d0, !airway resistance (cmH2O/L.s)
     &  mean_PVR*1.d6/98.0665d0/XP(1,1,nj_flow,1), !total resistance
     &  totalC*98.0665d0/1.d6, !total model compliance
     &  -mean_PVR/98.0665d0, !-(Pmouth - Palv) (cmH2O)
     &  ppl_current/98.0665d0, !Ppl (cmH2O)
     &  mean_Ptp/98.0665d0, !mean Ptp (cmH2O)
     &  CE(nm_vol_bel,1)/1.d3, !total model volume (mL)
     &  Pmus/98.0665d0,dpl*100.d0 !Pmuscle (cmH2O)
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      
      CALL EXITS('DIRECTSOLUTION')
      RETURN
 9999 CALL ERRORS('DIRECTSOLUTION',ERROR)
      CALL EXITS('DIRECTSOLUTION')
      RETURN 1
      END
