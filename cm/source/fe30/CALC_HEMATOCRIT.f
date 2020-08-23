      SUBROUTINE CALC_HEMATOCRIT(nb,NEELEM,NENP,NPNE,NORD,NYNE,CE,ERR2,
     &  YP,ERROR,*)

C#### Subroutine: CALC_HEMATOCRIT
C###  Description:
C###     CALC_HEMATOCRIT calculates the hematocrit (Hd) of the blood
C###     in the pulmonary capillary network. An initial hematocrit
C###     distribution is assumed (in ipinit) to give a flow & pressure
C###     solution. Solutions are passed in and the hematocrit
C###     distribution is recalculated.
C***   Created by KSB, October 2001.
C***
C***   This method of determining Hematocrit from Kamm et al, 2000 based
C***   on method by Levin et al, 1986.
C***   CE(1,ne) is the flux cutoff parameter
C***   CE(2,ne) is the preferential flux parameter

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER nb,NEELEM(0:NE_R_M),NENP(NPM,0:NEPM),NORD(5,NE_R_M),
     '  NPNE(NNM,NBFM,NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM)
      REAL*8 CE(NMM,NEM),ERR2,YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ne,ne2,nn,noelem,noelem2,nonode,np,NPLIST(0:NP_R_M),ny
      REAL*8 G_RATIO,Q_IN,Q_RATIO,QRBC_IN,RBC_flow
      LOGICAL FLOW_BALANCED

      CALL ENTERS('CALC_HEMATOCRIT',*9999)

      NPLIST(0)=0
      DO nonode=1,NPM
        NPLIST(nonode)=0 !initialising
      ENDDO !nonode
      DO noelem=1,NEELEM(0)!NELIST(0)
        ne=NEELEM(noelem) !NELIST(noelem)
        DO nn=1,NNT(nb) !does each node of an element
          np=NPNE(nn,nb,ne) !end node #
          FLOW_BALANCED=.FALSE.
          DO nonode=1,NPLIST(0)
            IF(np.EQ.NPLIST(nonode)) THEN !node already evaluated
              FLOW_BALANCED=.TRUE.
            ENDIF
          ENDDO !nonode
          IF(NENP(np,0).GT.1.AND..NOT.FLOW_BALANCED)THEN !junction node
            Q_IN=0.d0 !initialise flow into a junction
            QRBC_IN=0.d0
            DO noelem2=1,NENP(np,0)
              ne2=NENP(np,noelem2) !global element #
              ny=NYNE(1,1,0,1,ne2)
C...  ny# for flow in ne,assumes flow is na=1
              IF((np.EQ.NPNE(2,nb,ne2).AND.YP(ny,1).GT.0.d0)
     '          .OR.(np.EQ.NPNE(1,nb,ne2).AND.YP(ny,1).LT.0.d0))THEN
                IF(NORD(5,ne).LT.0.d0) THEN !converging symmetric branch
                  Q_IN=Q_IN+2.d0*DABS(YP(ny,1))
                  QRBC_IN=QRBC_IN+2.d0*(DABS(YP(ny,1))*CE(nm_Hd,ne2))
                ELSE
C...  flow into ne
                  Q_IN=Q_IN+DABS(YP(ny,1)) !sums flows into junction
                  QRBC_IN=QRBC_IN+DABS(YP(ny,1))*CE(nm_Hd,ne2)
                ENDIF !NORD
C...  inlet volume RBCs
              ENDIF
            ENDDO !noelem2
            DO noelem2=1,NENP(np,0)
              ne2=NENP(np,noelem2) !global element #
              ny=NYNE(1,1,0,1,ne2) !d.o.f number for flow in ne
              IF((np.EQ.NPNE(1,nb,ne2).AND.YP(ny,1).GT.0.d0)
     '          .OR.(np.EQ.NPNE(2,nb,ne2).AND.YP(ny,1).LT.0.d0))THEN
C...  flow out of ne
                Q_RATIO=DABS(YP(ny,1)/Q_IN)
               IF(Q_RATIO.GT.1.d0) Q_RATIO=1.d0 !fixes rounding
C...  flow out of ne/flow into junction
                IF(Q_RATIO.GT.0.d0.AND.Q_RATIO.LE.CE(1,ne2)) THEN
                  G_RATIO=0.d0
                ELSEIF(Q_RATIO.GT.CE(1,ne2).AND.
     '              Q_RATIO.LE.(1.d0-CE(1,ne2))) THEN
                  G_RATIO=1.d0/(1.d0+((1.d0-(Q_RATIO+CE(1,ne2)))
     '              /(Q_RATIO-CE(1,ne2)))**CE(2,ne2))
                ELSEIF(Q_RATIO.GT.(1.d0-CE(1,ne2))
     '              .AND.Q_RATIO.LE.1.0d0)THEN 
                  G_RATIO=1.d0
                ENDIF
                IF(G_RATIO.LE.DABS(YP(ny,1)/QRBC_IN)) THEN
                  RBC_flow=G_RATIO*QRBC_IN
                ELSE
                  RBC_flow=DABS(YP(ny,1)) !i.e Hd=1
                ENDIF
                !difference between previous Hd and current Hd value
                Hd=RBC_flow/DABS(YP(ny,1))
                IF(Hd.NE.0.d0) ERR2=ERR2+((CE(nm_Hd,ne2)-Hd)/Hd**2.d0)
                CE(nm_Hd,ne2)=Hd  !RBC_flow/DABS(YP(ny,1))!new Hd value
              ENDIF !flow into ne
            ENDDO !noelem2
            NPLIST(0)=NPLIST(0)+1 !adds node just evaluated
            NPLIST(NPLIST(0))=np !to NPLIST (list of nodes evaluated)
          ENDIF !NENP
        ENDDO !nn
C        IF(ITYP3(nr,nx).EQ.3) THEN !capillary flow
C          CALL CALC_R_SEGMENT(nb,ne,NENP,NPNE,CE,XP,ERROR,*9999)
C        ELSE IF(ITYP3(nr,nx).EQ.5) THEN !Poiseuille flow (arteriole/venule)
C       CALL POISEUILLE_R(nb,ne,NPNE,NVJE,CE,XP,ERROR,*9999)
C       ENDIF !resitance recalculated in CALC_HEMODYNAMICS, called from SOLVE11
C... STORING FLOW SOLUTION IN YP FOR EXPORTING.Z
        ny=NYNE(1,1,1,1,ne) !ny#,assumes na=1 for flow,make general
        CE(nm_flow,ne)=YP(ny,1) !stores flow solution for export
      ENDDO !noelem

      CALL EXITS('CALC_HEMATOCRIT')
      RETURN
 9999 CALL ERRORS('CALC_HEMATOCRIT',ERROR)
      CALL EXITS('CALC_HEMATOCRIT')
      RETURN 1
      END


