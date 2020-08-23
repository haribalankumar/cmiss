      SUBROUTINE ASSEMBLE11_DYNAM(IBT,IDO,INP,ISC_GK,ISR_GK,nb,NBH,
     '  NBJ,ne,NENP,NHE,nhs_cap,NKJE,NORD,NPF,NPLIST2,NPNE,nr,
     '  NVJE,nx,NYNE,NYNP,NZ_ESED,nzr,CE,CG,CGE,CP,ED,EM,ER,ES,GK,GR,PG,
     &  RG,SE,STACK_ED,STACK_EM,STACK_ES,WG,XA,XAB,XG,XP,YG,ZE,
     &  ZG,UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*)

C#### Subroutine: ASSEMBLE11_DYNAM
C###  Description:
C###    ASSEMBLE11_DYNAM assembles the global unreduced matrices GK,GD,
C###    for time dependent pulmonary transport problems.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'host00.inc'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GK(NISC_GKM),ISR_GK(NISR_GKM),nb,NBH(NHM,NCM),NBJ(NJM),ne,
     '  NENP(NPM,0:NEPM),NHE,nhs_cap,NKJE(NKM,NNM,NJM),
     '  NORD(5,NE_R_M),NPF(9),NPLIST2(0:NPM),NPNE(NNM,NBFM,NEM),nr,nx,
     '  NVJE(NNM,NBFM,NJM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM),NZ_ESED(0:12),nzr
      REAL*8 CE(NMM),CG(NMM,NGM),CGE(NMM,NGM),CP(NMM,NPM),
     '  GK(NZ_GK_M),GR(NYROWM),
     '  PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM),
     '  STACK_ED(12),STACK_EM(12),STACK_ES(12),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XAB(NORM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM),ZE(NSM,NHM),ZG(NHM,NUM),
     &  ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     &  ES(NHM*NSM,NHM*NSM),RG(NGM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
      LOGICAL UPDATE_MATRIX,UPDATE_VECTOR
!     Local Variables
      INTEGER mh,mhs,mhx,ms,na,ne2,nh,nhs,nhx,nn,noelem2,nonode,np,
     '  np1,np2,ns,ny1,ny2,ny3,nz,nzi
      REAL*8 flow_term
      LOGICAL FLOW_BALANCED

      CALL ENTERS('ASSEMBLE11_DYNAM',*9999)

      IF(ITYP3(nr,nx).LE.2.OR.ITYP3(nr,nx).EQ.5)THEN
        CALL XPXE(NBJ,NKJE,NPF,NPNE(1,1,ne),nr,NVJE,SE,XA(1,1,ne),XE,XP,
     '    ERROR,*9999)
        CALL XPES30(IBT,IDO,INP,NBH,NBJ,ne,NHE,NORD(5,ne),NPNE,nr,nx,
     '    CE,CG,CGE,CP,ED,EM,ER,ES,PG,RG,SE,WG,XE,XG,YG,ZE,ZG,
     '    UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*9999)
C*** Assemble element matrices into global matrices
        nzi=0
        mhs=0
        DO mhx=1,NHE
          mh=NH_LOC(mhx,nx)
c          DO ms=1,NST(nb)
          DO ms=1,2
            mhs=mhs+1
            np1=NPNE(ms,nb,ne)
            ny1=NYNP(1,1,mh,np1,0,1)
            nhs=0
            DO nhx=1,NHE
              nh=NH_LOC(nhx,nx)
c              DO ns=1,NST(nb)
              DO ns=1,2
                nhs=nhs+1
                np2=NPNE(ns,nb,ne)
                ny2=NYNP(1,1,nh,np2,0,1)
                CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '            NZT(1,nx),ISC_GK,ISR_GK,6,ERROR,*9999)
                IF(((mhx.EQ.1).OR.(mhx.EQ.2.AND.nhx.EQ.2))
     '            .AND.nz.NE.0)THEN
                  nzi=nzi+1
                  STACK_ES(nzi)=ES(mhs,nhs)
                  STACK_ED(nzi)=ED(mhs,nhs)
                  STACK_EM(nzi)=EM(mhs,nhs)
                  NZ_ESED(nzi)=nz
                ENDIF !mh etc.
              ENDDO !ns
            ENDDO !nh
            GR(ny1)=GR(ny1)+ER(mhs)
          ENDDO !ms
        ENDDO !mh
        IF(UPDATE_MATRIX) NZ_ESED(0)=nzi

      ELSE !non-XPES30 set-up
        IF(ITYP3(nr,nx).GE.3)THEN !flow
          nb=NBH(NH_LOC(1,nx),1)
          mhs=0 !ES counter
          IF(NAT(nb).EQ.0) THEN !must have 1 auxillary parameter
            WRITE(ERROR,
     &        '(''Must have 1 auxillary element parameter in basis'')')
            GO TO 9999
          ENDIF
          nzr=nzr+1
          np1=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
           DO na=1,NAT(nb)
            nhs_cap=nhs_cap+1
         ! i.e set up P1-P2-R(12)*Q(12)=0
            DO nhx=1,NHE
              nh=NH_LOC(nhx,nx)
              DO ns=1,NST(nb)
c            DO nh=1,NHE
c              DO ns=1,2
                ny2=NYNP(1,1,nh,NPNE(ns,nb,ne),0,1)
                CALL SPARSE(nhs_cap,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '            NZT(1,nx),ISC_GK,ISR_GK,6,ERROR,*9999)
                IF(ns.EQ.1) THEN
                  GK(nz)=1.0d0 !contribution from P1 (ie 1)
                ELSEIF(ns.EQ.2) THEN
                  GK(nz)=-1.0d0 !contribution from P2 (ie -1)
                ENDIF !ns
              ENDDO !ns
              ny3=NYNE(na,nh,1,1,ne) !flow in element
              CALL SPARSE(nhs_cap,ny3,NYT(1,1,nx),nz,NZ_GK_M,
     '          NZT(1,nx),ISC_GK,ISR_GK,6,ERROR,*9999)
              IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.6)THEN !Kamm capillary blood flow
                GK(nz)=-CE(nm_Rseg) !contribution from Q(1-2)
                IF(GK(nz).EQ.0.D0) then
                  WRITE(*,*) "WHY...",ne
                ENDIF
              ELSEIF(ITYP3(nr,nx).EQ.4)THEN !Simple P-R-airflow solution
                 IF(COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6) THEN
                    GK(nz)=-XAB(nej_resis,ne) !contribution from Q(1-2)
                 ELSE
                    GK(nz)=XAB(nej_resis,ne)
                 ENDIF
              ENDIF
              IF(ITYP3(nr,nx).EQ.4)THEN
C               Apply gravity to all nodes store in GR
C... NB/ In perfusion models we don't want to include gravity in the 
C... capillary elements or the largest vessels.
               IF(COMPLIANCE.EQ.3.AND.COUPLE_VIA_LPM.EQ.'Y'.AND.
     &           (XAB(nej_cap,ne).EQ.1.d0.OR.NORD(1,ne).EQ.1).OR.
     &            XAB(nej_cap,ne).EQ.20)THEN!it is a capillary or an inlet/outlet
                  gr(nzr)=0.d0
               ELSE
                 IF(POSITION(1:7).EQ.'UPRIGHT')THEN
                    IF(COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6) THEN
                      gr(nzr)=PULMAT(1)*(XP(1,1,3,np1)-
     '                   XP(1,1,3,np2))*9810.d0*
     &                   gravfact*grav_vect(3)
                    ELSE
                      gr(nzr)=PULMAT(1)*ABS(XP(1,1,3,np1)-
     '                   XP(1,1,3,np2))*9810.d0*gravfact
                    ENDIF
                  ELSE IF(POSITION(1:8).EQ.'INVERTED')THEN
                    IF(COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6) THEN
                      gr(nzr)=PULMAT(1)*(XP(1,1,3,np1)-
     '                   XP(1,1,3,np2))*9810.d0*
     &                   gravfact*grav_vect(3)
                    ELSE
                      gr(nzr)=PULMAT(1)*ABS(XP(1,1,3,np2)-
     '                   XP(1,1,3,np1))*9810.d0*gravfact
                    ENDIF
                  ELSE IF(POSITION(1:5).EQ.'PRONE')THEN
                    IF(COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6) THEN
                      gr(nzr)=PULMAT(1)*(XP(1,1,2,np1)-
     '                   XP(1,1,2,np2))*9810.d0*
     &                   gravfact*grav_vect(2)
                    ELSE
                      gr(nzr)=PULMAT(1)*ABS(XP(1,1,2,np1)-
     '                   XP(1,1,2,np2))*9810.d0*gravfact
                    ENDIF
                  ELSE IF(POSITION(1:6).EQ.'SUPINE')THEN
                    IF(COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6) THEN
                      gr(nzr)=PULMAT(1)*(XP(1,1,2,np1)-
     '                   XP(1,1,2,np2))*9810.d0*
     &                   gravfact*grav_vect(2)
                    ELSE
                       gr(nzr)=PULMAT(1)*ABS(XP(1,1,2,np2)-
     '                   XP(1,1,2,np1))*9810.d0*gravfact
                    ENDIF
                  ELSE
                    gr(nzr)=0.d0
                  ENDIF
                ENDIF !CAPILLARY or INLET
              ENDIF !GRAVITY LOOP
              NZ_ESED(1)=nz !store the nz location for fast updating next iteration
            ENDDO !nhx
          ENDDO !na
          DO nn=1,NNT(nb) !balances at each node of ne
            FLOW_BALANCED=.FALSE.
            np=NPNE(nn,nb,ne) !balanced at end node of each elem
            DO nonode=1,NPLIST2(0) !where flows already balanced
              IF(np.EQ.NPLIST2(nonode))THEN
                FLOW_BALANCED=.TRUE. !flow balance already done
              ENDIF !np
            ENDDO !nonode
            IF((NENP(np,0).GT.1).AND.(.NOT.FLOW_BALANCED))THEN !balance
              nhs_cap=nhs_cap+1
              nzr=nzr+1
              DO na=1,NAT(nb) !# auxillary variables
                DO nhx=1,NH_LOC(0,nx) !# dependent variables
                  nh=NH_LOC(nhx,nx) !dependent variable #
                  DO noelem2=1,NENP(np,0) !elements at junction
                    ne2=NENP(np,noelem2) !global element #
                    ny2=NYNE(na,nh,1,1,ne2) !variable #
                    CALL SPARSE(nhs_cap,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '                NZT(1,nx),ISC_GK,ISR_GK,6,ERROR,*9999)
                    IF(NORD(5,ne2).EQ.0.OR.ne.EQ.ne2)  THEN !capillary
                      flow_term=1.d0   !Q1-Q2-Q3=0
                    ELSEIF(NORD(5,ne2).EQ.1) THEN
 !symmetric arteriole - flow decreases by 2: Q1-2*Q2=0
                      flow_term=2.d0
                    ELSEIF(NORD(5,ne2).EQ.-1) THEN !venule
 !symmetric venule - flow increases by 2: Q1-0.5*Q2=0 
                      flow_term=0.5d0
                   ENDIF
                    IF(np.EQ.NPNE(2,nb,ne2)) THEN !end node
                      GK(nz)=flow_term
                    ELSEIF(np.EQ.NPNE(1,nb,ne2)) THEN !start node
                      GK(nz)=-flow_term
                    ENDIF
                  ENDDO !noelem2
                ENDDO !nh
              ENDDO !na
              NPLIST2(0)=NPLIST2(0)+1
              NPLIST2(NPLIST2(0))=np !nodes where flow balanced
            ENDIF !NENP
          ENDDO !nn
        ENDIF !ITYP3(nr,nx).GE.3
      ENDIF !ITYP3(nr,nx).LE.2

      CALL EXITS('ASSEMBLE11_DYNAM')
      RETURN
 9999 CALL ERRORS('ASSEMBLE11_DYNAM',ERROR)
      CALL EXITS('ASSEMBLE11_DYNAM')
      RETURN 1
      END


