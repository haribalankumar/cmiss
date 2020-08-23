      SUBROUTINE UPPHI(IBT,IDO,INP,ISIZE_PHI,ISIZE_PHIH,
     '  LIST,LD,NBH,NBJ,NDDATA,
     '  NEELEM,NENP,NPF,NHE,NHP,NKH,NKHE,NPNE,NPNODE,NRE,
     '  NRLIST,NVHE,NVHP,NW,NYNE,NYNP,
     '  CURVCORRECT,PHI,PHI_H,SE,XID,YP,ZA,ZE,ZP,STRING,ERROR,*)

C#### Subroutine: UPPHI
C###  Description:
C###    Update PHI matrix with potential solutions
C***  Created By Leo Cheng on 15-JUL-2002

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISIZE_PHI(2),ISIZE_PHIH(2),
     '  LIST(0:NLISTM),LD(NDM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NRE(NEM),
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NW(NEM,3,NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),
     '  PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  SE(NSM,NBFM,NEM),XID(NIM,NDM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER ERR,IBEG,IEND,N3CO,nb,ne,ni,nlist,nolist,
     '  nr,NSAMPLES,nt,nx,SSTART,SEND
      REAL*8 TSTART,TEND,XI(3)
      CHARACTER TYPE*8
      LOGICAL ALL_REGIONS,PHIH

!     Functions
      INTEGER CALC_SAMPLE_FROM_TIME
      REAL*8 RFROMC,PXI
      LOGICAL ABBREV,CBBREV

      CALL ENTERS('UPPHI',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update PHI
C###  Parameter:    <from (nodes/data)[nodes]>
C###    Specify where to location evaluate the potentials from
C###  Parameter:    <electrodes (#s/all)[all]>
C###    Specify the list of electrodes to update
C###  Parameter:      <tstart (#)[0.0]>
C###    Specify the start time for time dependent problems
C###  Parameter:      <tend (#)[0.0]>
C###    Specify the end time for time dependent problems
C###  Parameter:    <(PHI/PHI_H)[PHI]>
C###    Specify which matrix to update
C###  Parameter:    <electrodes (#s/all)[all]>
C###    Specify which electrode numbers update into PHI
C###  Parameter:    <region #[1]>
C###    Specify the region which the nodes/electrodes belong to
C###  Parameter:    <frequency #>
C###    Over-ride the frequency specified in the iptrsf file.
C###    The value is unchanged if not specified.
C###  Description:
C###    Update the PHI matrix by interpolating potential solutions
C###    at electrodes/nodes

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<from (nodes/data)[nodes]>'
        OP_STRING(3)=BLANK(1:15)//'<electrodes (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<tstart #[0.0]>'
        OP_STRING(5)=BLANK(1:15)//'<tend #[0.0]>'
        OP_STRING(6)=BLANK(1:15)//'<(PHI/PHI_H)[PHI]>'
        OP_STRING(7)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(8)=BLANK(1:15)//'<frequency #>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPMESH',ERROR,*9999)
      ELSE
        nx=1 ! should set using classes


C LKC 11-JAN-2011 Over ride the frequency possibly set in the iptrsf file
        IF(CBBREV(CO,'FREQUENCY',3,noco+1,NTCO,N3CO)) THEN
          TRSF_FREQUENCY=RFROMC(CO(N3CO+1))
        ENDIF
 
        
        IF(CBBREV(CO,'PHI_H',4,noco+1,NTCO,N3CO)) THEN
          PHIH=.TRUE.
        ELSE
          PHIH=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)

        CALL ASSERT(NRLIST(0).EQ.1,
     '    'Only implemented for 1 region',ERROR,*9999)

        IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'NODES',2)) THEN
            TYPE='NODES'
          ELSE IF(ABBREV(CO(N3CO+1),'DATA',2)) THEN
            TYPE='DATA'
          ENDIF
        ELSE
          TYPE='NODES'
        ENDIF


        IF(CBBREV(CO,'ELECTRODES',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NLISTM,LIST(0),LIST(1),
     '      ERROR,*9999)

        ELSE !default to all
          IF(TYPE(1:4).EQ.'DATA') THEN

            CALL ASSERT(NDDATA(0,nr).GT.0,'No data points defined',
     '        ERROR,*9999)

            CALL ASSERT(CALL_XI,'>> Calculate XI positions first',
     '        ERROR,*9999)

            LIST(0)=NDDATA(0,nr)
            DO nlist=1,NDDATA(0,nr)
              LIST(nlist)=NDDATA(nlist,nr)
            ENDDO

          ELSE IF(TYPE(1:5).EQ.'NODES') THEN

            LIST(0)=NPNODE(0,nr)
            DO nlist=1,NPNODE(0,nr)
              LIST(nlist)=NPNODE(nlist,nr)
            ENDDO

          ELSE
            ERROR='Update code for electrode types'
            GOTO 9999
          ENDIF
        ENDIF !electrodes

        CALL ASSERT(LIST(0).GT.0,'>> There are no electrodes to update',
     '    ERROR,*9999)

        IF(CBBREV(CO,'TSTART',3,noco+1,NTCO,N3CO)) THEN
          TSTART=RFROMC(CO(N3CO+1))
        ELSE
          TSTART=0.D0
        ENDIF

        IF(CBBREV(CO,'TEND',3,noco+1,NTCO,N3CO)) THEN
          TEND=RFROMC(CO(N3CO+1))
        ELSE
          TEND=0.D0
        ENDIF

        ERR=0
        SSTART=CALC_SAMPLE_FROM_TIME(TSTART,ERR,ERROR)
        IF(ERR.NE.0) GOTO 9999
        SEND=CALC_SAMPLE_FROM_TIME(TEND,ERR,ERROR)
        IF(ERR.NE.0) GOTO 9999

        NSAMPLES=SEND-SSTART+1
        CALL ASSERT(NSAMPLES.EQ.1,
     '    'Only implemented for 1 time step',ERROR,*9999)

        CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '    NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,
     '    NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)



        DO nt=SSTART,SEND
          DO nlist=1,LIST(0)
            nolist=LIST(nlist)
            IF(TYPE(1:4).EQ.'DATA') THEN
              ne=LD(nolist)
              nb=NBH(NH_LOC(1,nx),1,ne)
              DO ni=1,NIT(nb)
                XI(ni)=XID(ni,nolist)
              ENDDO
            ELSEIF(TYPE(1:5).EQ.'NODES') THEN
              CALL CALC_NP_XI(IBT,INP,NBJ,ne,NENP,nolist,NPNE,NRLIST,
     '          XI,ERROR,*9999)
              nb=NBH(NH_LOC(1,nx),1,ne)
            ENDIF

            CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     '        NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),NRE(ne),
     '        NVHE(1,1,1,ne),NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),
     '        SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)

            IF(PHIH) THEN
              PHI_H(nlist,nt)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '          INP(1,1,nb),nb,1,XI,ZE)
            ELSE
              PHI(nlist,nt)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '          INP(1,1,nb),nb,1,XI,ZE)
            ENDIF
          ENDDO !nolist (nd)
        ENDDO !nt

C*** Setting ISIZE_PHI(1) and ISIZE_PHI(2)

        IF(PHIH) THEN
          IF(ISIZE_PHIH(1).NE.LIST(0)) THEN
            ISIZE_PHIH(1)=LIST(0)
            IEND=0
            CALL APPENDC(IEND,
     '        '>> WARNING: Updating # electrode samples in PHI_H to ',
     '        OP_STRING(1))
            CALL APPENDI(IEND,LIST(0),OP_STRING(1))
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          ENDIF
          IF(ISIZE_PHIH(2).LT.SEND) THEN
            ISIZE_PHIH(2)=SEND
            IEND=0
            CALL APPENDC(IEND,
     '        '>> WARNING: Updating # time samples in PHI_H to ',
     '        OP_STRING(1))
            CALL APPENDI(IEND,SEND,OP_STRING(1))
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE
          IF(ISIZE_PHI(1).NE.LIST(0)) THEN
            ISIZE_PHI(1)=LIST(0)
            IEND=0
            CALL APPENDC(IEND,
     '        '>> WARNING: Updating # electrode samples in PHI to ',
     '        OP_STRING(1))
            CALL APPENDI(IEND,LIST(0),OP_STRING(1))
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          ENDIF

          IF(ISIZE_PHI(2).LT.SEND) THEN
            ISIZE_PHI(2)=SEND
            IEND=0
            CALL APPENDC(IEND,
     '        '>> WARNING: Updating # time samples in PHI to ',
     '        OP_STRING(1))
            CALL APPENDI(IEND,SEND,OP_STRING(1))
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          ENDIF

        ENDIF
        EVALUATE_PHI=.TRUE.

      ENDIF

      CALL EXITS('UPPHI')
      RETURN
 9999 CALL ERRORS('UPPHI',ERROR)
      CALL EXITS('UPPHI')
      RETURN 1
      END


