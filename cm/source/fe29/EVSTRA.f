      SUBROUTINE EVSTRA(NBJ,NEELEM,NELIST,NGLIST,NQNE,NQS,NQSCNB,
     '  NQXI,NRLIST,NXLIST,YQS,XIG,STRING,ERROR,*)

C#### Subroutine: EVSTRA
C###  Description:
C###    Evaluates strains - intially used to calculate strains from
C###    extension ratio values given at grid points
C**** Created by David Nickerson, July 2001

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NGLIST(0:NGM),NQNE(NEQM,NQEM),NQS(NEQM),
     '  NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NXLIST(0:NXM)
      REAL*8 YQS(NIQSM,NQM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,IJF,IJS,IKF,IKS,N3CO,nb,ne,ng,ni,ni1,
     '  ni2,ni3,NITB,no_nelist,no_nglist,no_nrlist,nq,nq_nearest,
     '  NQENG(27),nqsc,nr,nx,nxc,YQS_IDX
      REAL*8 DS,DSH,extension_ratio,strain,XI(3)
      CHARACTER FILE*(MXCH),LOCATION_TYPE*100,METHOD*100
      LOGICAL ALL_REGIONS,ABBREV,CBBREV,OPFILE
C     CHARACTER LOCATION*100

      CALL ENTERS('EVSTRA',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate strain<;FILENAME>
C###  Parameter:        <from cell_extension_ratios>
C###    Specify how to calculate the strains.
C###  Parameter:        <index #[1]>
C###    For the case when using cell_extension_ratios, <index> specifies
C###    the position of extension ratio variable in YQS (can be found
C###    via the "fem inquire" command.
C###  Parameter:        <at gauss/grid all/GROUP/#s [grid all]>
C###    Specify where you would like the strain evaluated - choose
C###    either to evaluate the strain at/close to a given local Gauss
C###    point(s) or specify a global grid point number(s).
C###  Parameter:        <nearest/surrounding[nearest]>
C###    When evaluating at a Guass point, return the value either at
C###    the nearest grid point or for all surrounding grid points.
C###    When evaluating at a grid point, return the value either at the
C###    grid point only or at the grid point and all its neighbours.
C###  Parameter:        <element all/GROUP/#[all]>
C###    Specifies the element containing the Gauss point.
C###  Parameter:        <region all/#s[all]>
C###    The region(s).
C###  Parameter:        <class #[1]>
C###    The class.
C###  Description:
C###    Evaluate the strain from a given source. Currently only
C###    implemented for the case when you wish to evaluate strains from
C###    grid point extension ratio values (i.e. YQS(index)). For this
C###    case, you can specify either a local Gauss point number or a
C###    global grid point number. When a Gauss point number is given,
C###    either the strain at the nearest grid point or the strain at
C###    all the grid points surrounding it evaluating and displayed.
C###    Similarly, if a grid point number is specified, either the
C###    strain at that point or that point and all its neighbours is
C###    returned.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<from cell_extension_ratios>'
        OP_STRING(3)=BLANK(1:15)//'<index #[1]>'
        OP_STRING(4)=BLANK(1:15)
     '    //'<at gauss/grid all/GROUP/#s [grid all]>'
        OP_STRING(5)=BLANK(1:15)//'<nearest/surrounding[nearest]>'
        OP_STRING(6)=BLANK(1:15)//'<element all/GROUP/#[all]>'
        OP_STRING(7)=BLANK(1:15)//'<region all/#s[all]>'
        OP_STRING(8)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVSTRA',ERROR,*9999)
      ELSE

        !Open the file if required
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opstrain','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        !Set the method of strain evaluation
        IF(CBBREV(CO,'FROM',4,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'CELL_EXTENSION_RATIOS',4)) THEN
            METHOD='CELL_EXTENSION_RATIOS'
          ELSE
            CO(noco+1)='?'
            GO TO 1
          ENDIF
        ELSE
          METHOD='CELL_EXTENSION_RATIOS'
        ENDIF

        IF(METHOD(1:21).EQ.'CELL_EXTENSION_RATIOS') THEN
          !Get the index into YQS
          IF(CBBREV(CO,'INDEX',3,noco+1,NTCO,N3CO)) THEN
            YQS_IDX = IFROMC(CO(N3CO+1))
            CALL ASSERT(YQS_IDX.LE.NIQSM,'Invalid index',ERROR,*9999)
          ELSE
            YQS_IDX = 1
          ENDIF
        ENDIF

        !Get the location(s) to evaluate the strain
        IF(CBBREV(CO,'AT',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'GAUSS',4)) THEN
C           LOCATION = 'GAUSS'
          ELSE
            CO(noco+1) = '?'
            GO TO 1
          ENDIF
        ELSE
C         LOCATION = 'GAUSS'
        ENDIF

        !Get the list of local Gauss point numbers
        IF(CBBREV(CO,'GAUSS',4,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NGM,NGLIST(0),NGLIST(1),ERROR,*9999)
        ELSE
          !use geometric gauss points
          NGLIST(0)=NGT(NBJ(1,1))
          DO no_nglist=1,NGLIST(0)
            NGLIST(no_nglist)=no_nglist
          ENDDO
        ENDIF

        !Set the location type
        IF(CBBREV(CO,'NEAREST',4,noco+1,NTCO,N3CO)) THEN
          LOCATION_TYPE = 'NEAREST'
        ELSEIF(CBBREV(CO,'SURROUNDING',2,noco+1,NTCO,N3CO)) THEN
          LOCATION_TYPE = 'SURROUNDING'
        ELSE
          LOCATION_TYPE = 'NEAREST'
        ENDIF

        CALL STRING_TRIM(LOCATION_TYPE,IBEG,IEND)
        DO no_nrlist=1,NRLIST(0)
          !the region number
          nr = NRLIST(no_nrlist)
          DO no_nelist=1,NELIST(0)
            !the element number
            ne = NELIST(no_nelist)
            !the grid scheme for element ne
            nqsc = NQS(ne)
            DO no_nglist=1,NGLIST(0)
              ng = NGLIST(no_nglist)

              !Find the grid point which is closest to the current
              !gauss point in xi space - from GEN_GRID_MAP
              DSH=RMAX
              IJF=1
              IKF=1
              IJS=1
              IKS=1
              nb = NQSCNB(nqsc)
              NITB = NIT(nb)
              IF(NITB.GE.2) THEN
                IJF=NQXI(2,nqsc)
              ENDIF
              IF(NITB.EQ.3) THEN
                IKF=NQXI(3,nqsc)
                IKS=1
              ENDIF
              !Loop over all the grid points in the element
              DO ni3=IKS,IKF
                DO ni2=IJS,IJF
                  DO ni1=1,NQXI(1,nqsc)
                    XI(1)=DBLE(ni1-1)/DBLE(NQXI(1,nqsc)-1)
                    IF(IJF.GT.1) THEN
                      XI(2)=DBLE(ni2-1)/DBLE(NQXI(2,nqsc)-1)
                    ELSE
                      XI(2) = 0.0d0
                    ENDIF
                    IF(IKF.GT.1) THEN
                      XI(3)=DBLE(ni3-1)/DBLE(NQXI(3,nqsc)-1)
                    ELSE
                      XI(3) = 0.0d0
                    ENDIF
                    DS=0.0d0
                    DO ni=1,NITB
                      DS=DS+DABS(XIG(ni,ng,NBJ(1,1))-XI(ni))**2
                    ENDDO !ni
                    DS=DSQRT(DS)
                    IF(DS.LT.DSH) THEN
                      DSH=DS
                      !store local grid point number
                      nq_nearest=ni1+((ni2-1)*NQXI(1,nqsc))+((ni3-1)*
     '                  NQXI(1,nqsc)*NQXI(2,nqsc))
                    ENDIF
                  ENDDO !ni1
                ENDDO !ni2
              ENDDO !ni3
              CALL ASSERT(DSH.LT.RMAX,'>>No grid points found',
     '          ERROR,*9999)

              IF(DOP) THEN
                WRITE(OP_STRING,'(''Grid: '',I6,'' Gauss: '',I6)') nq,ng
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

              IF(LOCATION_TYPE(1:11).EQ.'SURROUNDING') THEN
                !Calulate the local grid point numbers which surround
                !thegiven gauss point
                nq = nq_nearest
                IF(NITB.EQ.1) THEN
                  NQENG(1)=nq-1
                  NQENG(2)=nq
                  NQENG(3)=nq+1
                ELSE IF(NITB.EQ.2) THEN
                  NQENG(1)=nq-NQXI(1,nqsc)-1
                  NQENG(2)=nq-NQXI(1,nqsc)
                  NQENG(3)=nq-NQXI(1,nqsc)+1
                  NQENG(4)=nq-1
                  NQENG(5)=nq
                  NQENG(6)=nq+1
                  NQENG(7)=nq+NQXI(1,nqsc)-1
                  NQENG(8)=nq+NQXI(1,nqsc)
                  NQENG(9)=nq+NQXI(1,nqsc)+1
                ELSE IF(NITB.EQ.3) THEN
                  NQENG(1)=nq-NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))-1
                  NQENG(2)=nq-NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))
                  NQENG(3)=nq-NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))+1
                  NQENG(4)=nq-(NQXI(1,nqsc)*NQXI(2,nqsc))-1
                  NQENG(5)=nq-(NQXI(1,nqsc)*NQXI(2,nqsc))
                  NQENG(6)=nq-(NQXI(1,nqsc)*NQXI(2,nqsc))+1
                  NQENG(7)=nq+NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))-1
                  NQENG(8)=nq+NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))
                  NQENG(9)=nq+NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))+1
                  NQENG(10)=nq-NQXI(1,nqsc)-1
                  NQENG(11)=nq-NQXI(1,nqsc)
                  NQENG(12)=nq-NQXI(1,nqsc)+1
                  NQENG(13)=nq-1
                  NQENG(14)=nq
                  NQENG(15)=nq+1
                  NQENG(16)=nq+NQXI(1,nqsc)-1
                  NQENG(17)=nq+NQXI(1,nqsc)
                  NQENG(18)=nq+NQXI(1,nqsc)+1
                  NQENG(19)=nq-NQXI(1,nqsc)+(NQXI(1,nqsc)*
     '              NQXI(2,nqsc))-1
                  NQENG(20)=nq-NQXI(1,nqsc)+(NQXI(1,nqsc)*NQXI(2,nqsc))
                  NQENG(21)=nq-NQXI(1,nqsc)+(NQXI(1,nqsc)*
     '              NQXI(2,nqsc))+1
                  NQENG(22)=nq+(NQXI(1,nqsc)*NQXI(2,nqsc))-1
                  NQENG(23)=nq+(NQXI(1,nqsc)*NQXI(2,nqsc))
                  NQENG(24)=nq+(NQXI(1,nqsc)*NQXI(2,nqsc))+1
                  NQENG(25)=nq+NQXI(1,nqsc)+(NQXI(1,nqsc)*
     '              NQXI(2,nqsc))-1
                  NQENG(26)=nq+NQXI(1,nqsc)+(NQXI(1,nqsc)*NQXI(2,nqsc))
                  NQENG(27)=nq+NQXI(1,nqsc)+(NQXI(1,nqsc)*
     '              NQXI(2,nqsc))+1
                ENDIF
                IF(DOP) THEN
                  DO ni=1,3**NITB
                    WRITE(OP_STRING,'(''Local grid numbers'',I6)')
     '                NQENG(ni)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDDO !ni
                ENDIF
              ENDIF !LOCATION_TYPE

              !Evaluate the strain at the required grid points
              IF(METHOD(1:21).EQ.'CELL_EXTENSION_RATIOS') THEN
                !Evaluate strains from the cellular extension ratios
                IF(LOCATION_TYPE(1:7).EQ.'NEAREST') THEN
                  !Evaluate only at the nearest grid point
                  nq = NQNE(ne,nq_nearest)
                  extension_ratio = YQS(YQS_IDX,nq)
                  strain = 0.5d0 * (extension_ratio**2 - 1.0d0)
                  WRITE(OP_STRING,'(''  Region:             '',I6)') nr
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(''  Element:            '',I6)') ne
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(''  Gauss Point:        '',I6)') ng
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(''  Nearest grid point: '',I6)') nq
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(''    Strain: '',E12.5,'
     '              //'''  Extension ratio: '',E12.5)') strain,
     '              extension_ratio
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'()')
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ELSEIF(LOCATION_TYPE(1:11).EQ.'SURROUNDING') THEN
                  !Evaluate at all the grid points surrounding the Gauss
                  !point
                ENDIF
              ENDIF

            ENDDO !ng
          ENDDO !ne
        ENDDO !nr

        !Close the file
        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF

      ENDIF

      CALL EXITS('EVSTRA')
      RETURN
 9999 CALL ERRORS('EVSTRA',ERROR)
      CALL EXITS('EVSTRA')
      RETURN 1
      END

