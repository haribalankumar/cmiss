      SUBROUTINE IPMAT3_CELL(CELL_ICQS_SPATIAL,CELL_ICQS_VALUE,
     '  CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,IBT,ICQS_SPATIAL,IDO,
     '  IICQS_SPATIAL,IRCQS_SPATIAL,INP,NEELEM,NELIST,nmq,NPNE,NPNODE,
     '  NQET,NQLIST,NQNE,NQS,NQXI,nr,nx,CELL_CP,CELL_RCQS_VALUE,
     '  CELL_YQS_VALUE,RCQS_SPATIAL,XE,YQS,CELL_ICQS_NAMES,
     '  CELL_RCQS_NAMES,CELL_YQS_NAMES,ERROR,*)

C#### Subroutine: IPMAT3_CELL
C###  Description:
C###    IPMAT3_CELL inputs material parameters for cell
C###    problems. First get the variant number used to set the default
C###    values for all grid points. Then overwrite the YQS array with
C###    spatial variance, and set-up the IICQS_SPATIAL, ICQS_SPATIAL,
C###    IRCQS_SPATIAL, and RCQS_SPATIAL arrays to store any spatial
C###    variance for all other parameters.

C#### Variable: CELL_VARIANTS_USED
C###  Type: LOGICAL
C###  Set_up: IPMAT3_CELL
C###  Description:
C###    CELL_VARIANTS_USED is set to .TRUE. if more than one variant
C###    is used in the model, otherwise it is .FALSE.

C#### Variable: CELL_SPATIALLY_VARYING
C###  Type: LOGICAL
C###  Set_up: IPMAT3_CELL
C###  See-Also: *_INIT_GRID
C###  Description:
C###    CELL_SPATIALLY_VARYING is set to .TRUE. if any parameters are
C###    spatially varying, otherwise it is .FALSE. This variable should
C###    always be .TRUE. for biophysical cellular models, as the variant
C###    is a spatially varying parameter.

C *** Created DPN June 1999 (from CPB's Physiome Merge)

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'file01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER CELL_ICQS_SPATIAL(NQIM,NQVM),CELL_ICQS_VALUE(NQIM,NQVM),
     '  CELL_RCQS_SPATIAL(NQRM,NQVM),CELL_YQS_SPATIAL(NIQSM,NQVM),
     '  IBT(3,NIM,NBFM),ICQS_SPATIAL(NQISVM,NQM),
     '  IDO(NKM,NNM,0:NIM,NBFM),
     '  IICQS_SPATIAL(0:NQISVM,NQVM),IRCQS_SPATIAL(0:NQRSVM,NQVM),
     '  INP(NNM,NIM,NBFM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  nmq,NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),
     '  NQLIST(0:NQM),NQNE(NEQM,NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),nr,nx
      REAL*8 CELL_CP(NMQM,NPM),CELL_RCQS_VALUE(NQRM,NQVM),
     '  CELL_YQS_VALUE(NIQSM,NQVM),
     '  RCQS_SPATIAL(NQRSVM,NQM),XE(NSM,NJM),YQS(NIQSM,NQM)
      CHARACTER CELL_ICQS_NAMES(NQIM,NQVM)*(*),CELL_RCQS_NAMES(NQIM,
     '  NQVM)*(*),CELL_YQS_NAMES(NQIM,NQVM)*(*),ERROR*(*)
!     Local Variables
      INTEGER i,ICHAR,INFO,INTDUMMY,n,NOQUES,nq,nqcv,nqq,nqv
      INTEGER OFFSET_ICQS_SPATIAL,OFFSET_RCQS_SPATIAL
      INTEGER PRE_VARIANT,IBEG,IEND,CELL_VAR_TEMP(CELL_NUM_VARIANTS)
      INTEGER region_cell_var_num,variant_loop
c      REAL*8 REALDUMMY(1)
      CHARACTER*50 TITLESTRING
      CHARACTER*5 CHAR5

      CALL ENTERS('IPMAT3_CELL',*9999)

      CALL ASSERT(ITYP2(nr,nx).EQ.9,
     '  '>>Must define a cellular model first',
     '  ERROR,*9999)

C *** Set the cell variant number for each grid point
C *** The variant number is always spatially varying
      CALL ASSERT(NQISVM.GE.1,'NQISVM must be >= 1',ERROR,*9999)
      DO nqv=1,CELL_NUM_VARIANTS
        IICQS_SPATIAL(1,nqv)=CELL_VARIANT_OFFSET
      ENDDO
      IF(IOTYPE.NE.3) THEN
        DO nq=NQR(1,nr),NQR(2,nr)
          ICQS_SPATIAL(1,nq)=1
        ENDDO !nq
      ENDIF

      !initialise
      NOQUES=0
      ICHAR=0
      region_cell_var_num=0
      FORMAT='('' Enter the cell variant for each collocation '
     '  //'point:'')'
      CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,1,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      nqq=NQR(1,nr)
 6100 FORMAT='($,'' Enter collocation point #s/name [EXIT]: '',I5)'
C *** DPN 14-Apr-2004 - need to use NQLIST for all grid prompts
      IF(IOTYPE.EQ.3) THEN
        IF(nqq.LE.NQR(2,nr)) THEN
          !IDATA(1)=nqq
          NQLIST(0)=1
          NQLIST(1)=nqq
        ELSE
          !IDATA(0)=0
          !IDATA(1)=0
          NQLIST(0)=0
          NQLIST(1)=0
        ENDIF
        nqq=nqq+1
      ENDIF
 6200 CDATA(1)='GRIDS' !for use with group input
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,
     '  0,NQT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '  ERROR,*9999)
      IF(NQLIST(1).NE.0) THEN !not default exit
        IF(IOTYPE.NE.3) THEN
          !NQLIST(0)=IDATA(0)
          DO n=1,NQLIST(0)
            !NQLIST(n)=IDATA(n)
            nq=NQLIST(n)
            IF(nq.LT.NQR(1,nr).OR.nq.GT.NQR(2,nr)) THEN
              WRITE(OP_STRING,'('' >>Grid point '',I7,'' is not '
     '          //'in the current region'')') nq
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              GOTO 6200
            ENDIF
          ENDDO !n
        ENDIF !IOTYPE.NE.3
C       Define variant for first grid point in group. Rest of
C       group filled in at the end.
        nq=NQLIST(1)
        IDEFLT(1)=1
        FORMAT='($,'' The cell variant number is [1]: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=ICQS_SPATIAL(CELL_VARIANT_OFFSET,nq)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,1,CELL_NUM_VARIANTS,LDATA,LDEFLT,RDATA,RDEFLT,
     '    RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          ICQS_SPATIAL(CELL_VARIANT_OFFSET,nq)=IDATA(1)
C         Apply to all grid points in the group
          DO n=2,NQLIST(0)
            nq=NQLIST(n)
            ICQS_SPATIAL(CELL_VARIANT_OFFSET,nq)=IDATA(1)
          ENDDO !n
        ENDIF
C NEWS 20-NOV-2004 JHC store the number of cell variants defined in a region 
        region_cell_var_num=region_cell_var_num+1
        CELL_VAR_TEMP(region_cell_var_num)=IDATA(1)
C NEWE
        GOTO 6100 !For more grid points
      ENDIF
C NEWS 25-NOV-2004 JHC check for zero variants specified 
      IF(region_cell_var_num.EQ.0)THEN
        region_cell_var_num=CELL_NUM_VARIANTS
        CELL_VAR_TEMP(region_cell_var_num)=1
      ENDIF
C NEWE
C *** Initialise all the variables with the default values for the
C *** cell variant for the grid point

C *** State variables (and derived values) are stored at each grid point
      CELL_VARIANTS_USED=.FALSE.
      PRE_VARIANT=-1
      DO nq=NQR(1,nr),NQR(2,nr)
        nqcv=ICQS_SPATIAL(1,nq) ! the variant number
C ***   Test for more than one variant used
        IF(PRE_VARIANT.GE.0) THEN
          IF(nqcv.NE.PRE_VARIANT) THEN
            CELL_VARIANTS_USED=.TRUE.
          ENDIF
        ELSE
          PRE_VARIANT=nqcv
        ENDIF
        DO i=1,CELL_NUM_STATE(nqcv)
          YQS(CELL_STATE_OFFSET(nqcv)+i-1,nq)=
     '      CELL_YQS_VALUE(CELL_STATE_OFFSET(nqcv)+i-1,nqcv)
        ENDDO
      ENDDO

C *** Initialise the count of spatially varying integer and
C     real parameters
      IF (IOTYPE.NE.3) THEN
        NQISVT=1 !the variant
        NQRSVT=0
C        DO nqv=1,CELL_NUM_VARIANTS
C NEWS JHC 20-NOV-2004 modified to loop over cell variants which are defined in
C ipmatc file rather than looping over the whole cell variants
        DO variant_loop=1,region_cell_var_num
          nqv=CELL_VAR_TEMP(variant_loop)
          IICQS_SPATIAL(0,nqv)=1
          IRCQS_SPATIAL(0,nqv)=0
        ENDDO
        CELL_SPATIALLY_VARYING=.TRUE.
      ENDIF

C *** Loop through the variants defining any spatially varying variables

C      DO nqv=1,CELL_NUM_VARIANTS
C NEWS JHC 19-NOV-2004 modified to loop over cell variants which are defiend in C ipmatc file rather than looping over the whole cell variants
      DO variant_loop=1,region_cell_var_num
C NEWE
C        WRITE(CHAR5,'(I4)') nqv
C NEWS JHC 19-NOV-2004 search for cell variants which are defined in the region
        nqv=CELL_VAR_TEMP(variant_loop)
C NEWE
        WRITE(CHAR5,'(I4)') nqv 
        CALL STRING_TRIM(CHAR5,IBEG,IEND)
        FORMAT='(/'' Cellular variant: '//CHAR5(IBEG:IEND)//''')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

        ! Initialise offsets for this variant
        OFFSET_ICQS_SPATIAL=1 !due to variants
        OFFSET_RCQS_SPATIAL=0
        INTDUMMY=0
C ***   Don't ever change the order of these!!!!
C       Do state variables - stored at each grid point
        TITLESTRING='State variables:'
        CALL IPMAT3_CELL1(CELL_ICQS_VALUE,IBT,ICQS_SPATIAL,IDO,
     '    IICQS_SPATIAL,INP,IRCQS_SPATIAL,NIQSM,NEELEM,NELIST,nmq,NPNE,
     '    NPNODE,NQET,NQLIST,NQNE,NQS,nqv,NQXI,nr,CELL_NUM_STATE(nqv),
     '    CELL_STATE_OFFSET(nqv)-1,INTDUMMY,CELL_YQS_SPATIAL,3,
     '    CELL_CP,CELL_RCQS_VALUE,RCQS_SPATIAL,XE,YQS,
     '    CELL_YQS_NAMES,TITLESTRING,ERROR,*9999)
C       Do model variables
        TITLESTRING='Model variables:'
        CALL IPMAT3_CELL1(CELL_ICQS_VALUE,IBT,ICQS_SPATIAL,IDO,
     '    IICQS_SPATIAL,INP,IRCQS_SPATIAL,NQIM,NEELEM,NELIST,nmq,NPNE,
     '    NPNODE,NQET,NQLIST,NQNE,NQS,nqv,NQXI,nr,CELL_NUM_MODEL(nqv),
     '    CELL_MODEL_OFFSET(nqv)-1,OFFSET_ICQS_SPATIAL,
     '    CELL_ICQS_SPATIAL,1,CELL_CP,CELL_RCQS_VALUE,
     '    RCQS_SPATIAL,XE,YQS,CELL_ICQS_NAMES,TITLESTRING,ERROR,*9999)
C       Do control variables
        TITLESTRING='Control variables:'
        CALL IPMAT3_CELL1(CELL_ICQS_VALUE,IBT,ICQS_SPATIAL,IDO,
     '    IICQS_SPATIAL,INP,IRCQS_SPATIAL,NQIM,NEELEM,NELIST,nmq,NPNE,
     '    NPNODE,NQET,NQLIST,NQNE,NQS,nqv,NQXI,nr,CELL_NUM_CONTROL(nqv),
     '    CELL_CONTROL_OFFSET(nqv)-1,OFFSET_ICQS_SPATIAL,
     '    CELL_ICQS_SPATIAL,1,CELL_CP,CELL_RCQS_VALUE,
     '    RCQS_SPATIAL,XE,YQS,CELL_ICQS_NAMES,TITLESTRING,ERROR,*9999)
C       Do parameter variables
        TITLESTRING='Parameter variables:'
        CALL IPMAT3_CELL1(CELL_ICQS_VALUE,IBT,ICQS_SPATIAL,IDO,
     '    IICQS_SPATIAL,INP,IRCQS_SPATIAL,NQRM,NEELEM,NELIST,nmq,NPNE,
     '    NPNODE,NQET,NQLIST,NQNE,NQS,nqv,NQXI,nr,
     '    CELL_NUM_PARAMETERS(nqv),CELL_PARAMETERS_OFFSET(nqv)-1,
     '    OFFSET_RCQS_SPATIAL,CELL_RCQS_SPATIAL,2,CELL_CP,
     '    CELL_RCQS_VALUE,RCQS_SPATIAL,XE,YQS,CELL_RCQS_NAMES,
     '    TITLESTRING,ERROR,*9999)
C       Do protocol variables
        TITLESTRING='Protocol variables:'
        CALL IPMAT3_CELL1(CELL_ICQS_VALUE,IBT,ICQS_SPATIAL,IDO,
     '    IICQS_SPATIAL,INP,IRCQS_SPATIAL,NQRM,NEELEM,NELIST,nmq,NPNE,
     '    NPNODE,NQET,NQLIST,NQNE,NQS,nqv,NQXI,nr,
     '    CELL_NUM_PROTOCOL(nqv),CELL_PROTOCOL_OFFSET(nqv)-1,
     '    OFFSET_RCQS_SPATIAL,CELL_RCQS_SPATIAL,2,CELL_CP,
     '    CELL_RCQS_VALUE,RCQS_SPATIAL,XE,YQS,CELL_RCQS_NAMES,
     '    TITLESTRING,ERROR,*9999)
C       Do additional integer input variables
        TITLESTRING='Additional integer input variables:'
        CALL IPMAT3_CELL1(CELL_ICQS_VALUE,IBT,ICQS_SPATIAL,IDO,
     '    IICQS_SPATIAL,INP,IRCQS_SPATIAL,NQIM,NEELEM,NELIST,nmq,NPNE,
     '    NPNODE,NQET,NQLIST,NQNE,NQS,nqv,NQXI,nr,CELL_NUM_AII(nqv),
     '    CELL_AII_OFFSET(nqv)-1,OFFSET_ICQS_SPATIAL,CELL_ICQS_SPATIAL,
     '    1,CELL_CP,CELL_RCQS_VALUE,RCQS_SPATIAL,XE,YQS,
     '    CELL_ICQS_NAMES,TITLESTRING,ERROR,*9999)
C       Do additional real input variables
        TITLESTRING='Additional real input variables:'
        CALL IPMAT3_CELL1(CELL_ICQS_VALUE,IBT,ICQS_SPATIAL,IDO,
     '    IICQS_SPATIAL,INP,IRCQS_SPATIAL,NQRM,NEELEM,NELIST,nmq,NPNE,
     '    NPNODE,NQET,NQLIST,NQNE,NQS,nqv,NQXI,nr,CELL_NUM_ARI(nqv),
     '    CELL_ARI_OFFSET(nqv)-1,OFFSET_RCQS_SPATIAL,CELL_RCQS_SPATIAL,
     '    2,CELL_CP,CELL_RCQS_VALUE,RCQS_SPATIAL,XE,YQS,
     '    CELL_RCQS_NAMES,TITLESTRING,ERROR,*9999)

        ! Keep track of the maximum sizes
        IF (IICQS_SPATIAL(0,nqv).GT.NQISVT) NQISVT=IICQS_SPATIAL(0,nqv)
        IF (IRCQS_SPATIAL(0,nqv).GT.NQRSVT) NQRSVT=IRCQS_SPATIAL(0,nqv)
      ENDDO !variant loop (nqv)

      IF(NQISVT.GT.0.OR.NQRSVT.GT.0) THEN
        CELL_SPATIALLY_VARYING=.TRUE.
      ELSE
        CELL_SPATIALLY_VARYING=.FALSE.
      ENDIF

      CALL EXITS('IPMAT3_CELL')
      RETURN
 9999 CALL ERRORS('IPMAT3_CELL',ERROR)
      CALL EXITS('IPMAT3_CELL')
      RETURN 1
      END



