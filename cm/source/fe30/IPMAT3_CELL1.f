      SUBROUTINE IPMAT3_CELL1(CELL_ICQS_VALUE,IBT,ICQS_SPATIAL,
     '  IDO,IICQS_SPATIAL,INP,IRCQS_SPATIAL,LD_VAR,NEELEM,NELIST,nmq,
     '  NPNE,NPNODE,NQET,NQLIST,NQNE,NQS,nqv,NQXI,nr,NUMVARS,OFFSET,
     '  OFFSET1,SPATIAL,VARTYPE,CELL_CP,CELL_RCQS_VALUE,
     '  RCQS_SPATIAL,XE,YQS,NAMES,TITLESTRING,ERROR,*)

C#### Subroutine: IPMAT3_CELL1
C###  See-Also: IICQS_SPATIAL
C###  See-Also: IRCQS_SPATIAL
C###  See-Also: ICQS_SPATIAL
C###  See-Also: RCQS_SPATIAL
C###  See-Also: CELL_ASSIGN_SPATIAL
C###  Description:
C###  IPMAT3_CELL1 is used to input spatially varying material
C###  parameters for cellular models. The spatially varying parameters
C###  in a cellular model are flagged in the ipcell file and then
C###  the spatial variation is defined here. Spatial variation is now
C###  defined within a cellular variant, rather than forcing a given
C###  parameter index to be spatially varying across all variants. When
C###  prompting for the spatial variation only grid points for the
C###  current variant will have their values set. If the user is free to
C###  set the value over all grid points but only the grid points with
C###  the correct variant will be used in CELL_ASSIGN_SPATIAL - unless
C###  no grid points are given which are the correct variant, in which
C###  case am error is thrown.

C#### Variable: IICQS_SPATIAL(0:NQISVM,NQVM)
C###  Type: INTEGER
C###  Set_up: IPMAT3_CELL1
C###  See-Also: ICQS_SPATIAL
C###  See-Also: CELL_ASSIGN_SPATIAL
C###  Description:
C###  IICQS_SPATIAL(0,nqv) gives the number of spatially varying integer
C###  parameters in the ICQS array for cellular variant nqv.
C###  IICQS_SPATIAL(nqisv,nqv) gives the index in the ICQS array for
C###  spatially varying parameter nqisv and variant nqv. Together,
C###  IICQS_SPATIAL and ICQS_SPATIAL define the spatial variation of
C###  integer cellular material parameters. Note: There should always be
C###  at least one spatially varying integer parameter as the variant
C###  needs to be defined for each grid point.

C#### Variable: ICQS_SPATIAL(NQISVM,NQM)
C###  Type: INTEGER
C###  Set_up: IPMAT3_CELL1
C###  See-Also: IICQS_SPATIAL
C###  See-Also: CELL_ASSIGN_SPATIAL
C###  Description:
C###  ICQS_SPATIAL(nqisv,nq) gives the value if the nqisv-th spatially
C###  varying cellular integer parameter for the grid point nq. Although
C###  this array is of size NQM, it is only initialised for grid
C###  points with the correct variant.  Note the nqisv is not
C###  (necessary) the same index as that used in IICQS_SPATIAL.

C#### Variable: IRCQS_SPATIAL(0:NQRSVM,NQVM)
C###  Type: INTEGER
C###  Set_up: IPMAT3_CELL1
C###  See-Also: RCQS_SPATIAL
C###  See-Also: CELL_ASSIGN_SPATIAL
C###  Description:
C###  IRCQS_SPATIAL(0,nqv) gives the number of spatially varying real
C###  parameters in the RCQS array for cellular variant nqv.
C###  IRCQS_SPATIAL(nqrsv,nqv) gives the index in the RCQS array for
C###  spatially varying parameter nqrsv and variant nqv. Together,
C###  IRCQS_SPATIAL and RCQS_SPATIAL define the spatial variation of
C###  real cellular material parameters.

C#### Variable: RCQS_SPATIAL(NQRSVM,NQM)
C###  Type: INTEGER
C###  Set_up: IPMAT3_CELL1
C###  See-Also: IRCQS_SPATIAL
C###  See-Also: CELL_ASSIGN_SPATIAL
C###  Description:
C###  RCQS_SPATIAL(nqrsv,nq) gives the value if the nqrsv-th spatially
C###  varying cellular real parameter for the grid point nq. Although
C###  this array is of size NQM, it is only initialised for grid
C###  points with the correct variant.  Note the nqrsv is not
C###  (necessary) the same index as that used in IRCQS_SPATIAL.

C#### Comment: SPATIALLY VARYING CELLULAR PARAMETERS
C###  Description:
C###  Check out the see-also's to find out more.
C###  See-Also: IPMAT3_CELL1
C###  See-Also: IPMAT3_CELL
C###  See-Also: CELL_ASSIGN_SPATIAL

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'file01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER CELL_ICQS_VALUE(NQIM,NQVM),IBT(3,NIM,NBFM),
     '  ICQS_SPATIAL(NQISVM,NQM),
     '  IDO(NKM,NNM,0:NIM,NBFM),IICQS_SPATIAL(0:NQISVM,NQVM),
     '  INP(NNM,NIM,NBFM),IRCQS_SPATIAL(0:NQRSVM,NQVM),LD_VAR,
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  nmq,NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),
     '  NQLIST(0:NQM),NQNE(NEQM,NQEM),NQS(NEQM),nqv,NQXI(0:NIM,NQSCM),
     '  nr,NUMVARS,OFFSET,OFFSET1,SPATIAL(LD_VAR,NQVM),VARTYPE
      REAL*8 CELL_CP(NMQM,NPM),CELL_RCQS_VALUE(NQRM,NQVM),
     '  RCQS_SPATIAL(NQRSVM,NQM),XE(NSM,NJM),YQS(NIQSM,NQM)
      CHARACTER NAMES(LD_VAR,NQVM)*(*),TITLESTRING*(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,IBEG1,ICHAR,IEND,IEND1,II,IJ,IK,INFO,n,
     '  n1,nb,ne,neq,nii,nij,nik,niqs,nn,noelem,nonode,NOQUES,np,nq,
     '  nqq,nqsc,ns
      REAL*8 PXI,XI(3)
      LOGICAL ANYSPATIAL,INLIST,LINEARBASIS,VALUESET
      CHARACTER CHAR3*3,CHAR5*5,CHAR12*12

      CALL ENTERS('IPMAT3_CELL1',*9999)

      CALL ASSERT(IOTYPE.NE.3,
     '  '>>IPMAT3_CELL1 needs a bit more work before it will '//
     '  'write ipmatc files correctly',ERROR,*9999)

      !initialise
      NOQUES=0
      ICHAR=0

      CALL STRING_TRIM(TITLESTRING,IBEG,IEND)
      FORMAT='(/'' '//TITLESTRING(IBEG:IEND)//''')'
      CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      DO niqs=1,NUMVARS
        IF(SPATIAL(OFFSET+niqs,nqv).NE.0) THEN
          ANYSPATIAL=.TRUE.
        ELSE
          ANYSPATIAL=.FALSE.
        ENDIF
        IF(ANYSPATIAL) THEN
          ! We want to ensure that some values are set ??
          VALUESET = .FALSE.
C ***     OFFSET1 stores the current position in the ICQS_SPATIAL and
C ***     RCQS_SPATIAL arrays. Need to use separate variables for ints
C ***     and reals in IPMAT3_CELL().
          IF(VARTYPE.EQ.1.OR.VARTYPE.EQ.2) THEN
            OFFSET1=OFFSET1+1
          ENDIF
C ***     Can not define any spatial variance if there are no grid pts.
          CALL ASSERT(NQT.GT.0,'>>No grid points defined',ERROR,*9999)
C ***     Need to set the default value for all grid points, using the
C         variant already defined in IPMAT3_CELL for each grid point
          IF(IOTYPE.NE.3) THEN
            IF(VARTYPE.EQ.1) THEN !ICQS values
              IICQS_SPATIAL(0,nqv) = IICQS_SPATIAL(0,nqv) + 1
              CALL ASSERT(IICQS_SPATIAL(0,nqv).LE.NQISVM,
     '          'Increase NQISVM',ERROR,*9999)
              IICQS_SPATIAL(IICQS_SPATIAL(0,nqv),nqv)=OFFSET+niqs
C NEW JHC 20-NOV-2004 it makes more sense to loop over NQR(1,nr),NQR(2,nr)
C              DO nq=1,NQT !should use NQR(1,nr),NQR(2,nr) ??
              DO nq=NQR(1,nr),NQR(2,nr)
                IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                  ICQS_SPATIAL(OFFSET1,nq)=
     '              CELL_ICQS_VALUE(OFFSET+niqs,nqv)
C KAT 2006-07-11: Don't set other variants to 0 as they may have already
C                 been set up.
C                 ELSE
C                   ICQS_SPATIAL(OFFSET1,nq)=0 ! JHC maybe this aslo needs to be commented out
                ENDIF
              ENDDO
            ELSEIF(VARTYPE.EQ.2) THEN !RCQS values
              IRCQS_SPATIAL(0,nqv) = IRCQS_SPATIAL(0,nqv) + 1
              CALL ASSERT(IRCQS_SPATIAL(0,nqv).LE.NQRSVM,
     '          'Increase NQRSVM',ERROR,*9999)
              IRCQS_SPATIAL(IRCQS_SPATIAL(0,nqv),nqv)=OFFSET+niqs
C NEW JHC 20-NOV-2004 it makes more sense to loop over NQR(1,nr),NQR(2,nr)       
C              DO nq=1,NQT !should use NQR(1,nr),NQR(2,nr) ??
              DO nq=NQR(1,nr),NQR(2,nr)
                IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                  RCQS_SPATIAL(OFFSET1,nq)=
     '              CELL_RCQS_VALUE(OFFSET+niqs,nqv)
C NEW JHC 22-NOV-2004 don't see point of setting these to zero
C NEW GR 22-FEB-2006 Commenting this out because it prevents multiple
C material laws being used with spatially varying parameters.
                !ELSE
                ! RCQS_SPATIAL(OFFSET1,nq) = 0.0d0
                ENDIF
              ENDDO
            ELSEIF(VARTYPE.EQ.3) THEN !YQS values
              !do nothing, default values already set.
            ELSE
              ERROR='>>Invalid variable type'
              GOTO 9999
            ENDIF !VARTYPE
          ENDIF !IOTYPE.NE.3
          IDEFLT(1)=2
          CALL STRING_TRIM(NAMES(OFFSET+niqs,nqv),IBEG,IEND)
          FORMAT='(/'' Variable '//NAMES(OFFSET+niqs,
     '      nqv)(IBEG:IEND)//' is [2]: '''//
     '      '/''   (1) Piecewise constant (defined by elements)  '''//
     '      '/''   (2) Piecewise linear (defined by nodes)       '''//
     '      '/''   (3) Defined by grid points                    '''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) THEN
            IF(SPATIAL(OFFSET+niqs,nqv).LT.0) THEN
              IDATA(1)=2
            ELSE
              IDATA(1)=SPATIAL(OFFSET+niqs,nqv)
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            SPATIAL(OFFSET+niqs,nqv)=IDATA(1)
          ENDIF

          IF(SPATIAL(OFFSET+niqs,nqv).EQ.1) THEN !defined by elements
            WRITE(CHAR3,'(I3)') niqs
            noelem=0
 6100       FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
            IF(IOTYPE.EQ.3) THEN
              noelem=noelem+1
              IF(noelem.LE.NEELEM(0,nr)) THEN
                ne=NEELEM(noelem,nr)
                IDATA(1)=ne
                NELIST(0)=1
                NELIST(1)=ne
              ELSE
                IDATA(0)=0
                IDATA(1)=0
              ENDIF
            ENDIF
 6200       CDATA(1)='ELEMENTS' !for use with group input
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '        0,NET(nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
            IF(IDATA(1).NE.0) THEN !not default exit
              IF(IOTYPE.NE.3) THEN
                NELIST(0)=IDATA(0)
                DO n=1,IDATA(0)
                  NELIST(n)=IDATA(n)
                  ne=IDATA(n)
                  IF(.NOT.INLIST(ne,NEELEM(1,nr),NEELEM(0,nr),N1))
     '              THEN
                    WRITE(OP_STRING,'('' >>Element '',I5,'' is not '
     '                //'in the current region'')') ne
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 6200
                  ENDIF
                ENDDO !n
              ENDIF !IOTYPE.NE.3
C             Define parameter for first element in group. Rest of
C             group filled in at the end.
              ne=NELIST(1)
              IF(VARTYPE.EQ.1) THEN !ICQS parameter
C ***           The default values are given by the variant of the
C ***           first grid point.
                IDEFLT(1)=CELL_ICQS_VALUE(OFFSET+niqs,nqv)
                WRITE(CHAR5,'(I5)') IDEFLT(1)
                FORMAT='($,'' The value is ['//CHAR5(1:5)//
     '            ']: '',I10)'
                IF(IOTYPE.EQ.3) THEN !writing out values
                  nq=NQNE(ne,1)
                  IDATA(1)=ICQS_SPATIAL(OFFSET1,nq)
                ENDIF
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN !reading in
                  DO nqq=1,NQET(NQS(ne))
                    nq=NQNE(ne,nqq)
                    !only want to set the values for grid points with
                    !the correct variant
                    IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                      VALUESET = .TRUE.
                      ICQS_SPATIAL(OFFSET1,nq)=IDATA(1)
                    ENDIF
                  ENDDO
C                 Apply to all elements in the group
                  DO n=2,NELIST(0)
                    ne=NELIST(n)
                    DO nqq=1,NQET(NQS(ne))
                      nq=NQNE(ne,nqq)
                      !only want to set the values for grid points with
                      !the correct variant
                      IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                        VALUESET = .TRUE.
                        ICQS_SPATIAL(OFFSET1,nq)=IDATA(1)
                      ENDIF
                    ENDDO
                  ENDDO !n
                ENDIF
              ELSE IF(VARTYPE.EQ.2) THEN !RCQS values
C ***           The default values are given by the variant of the
C ***           first grid point.
                RDEFLT(1)=CELL_RCQS_VALUE(OFFSET+niqs,nqv)
                WRITE(CHAR12,'(D12.5)') RDEFLT(1)
                FORMAT='($,'' The value is ['//CHAR12(1:12)
     '            //']: '',D12.5)'
                IF(IOTYPE.EQ.3) THEN
                  nq=NQNE(ne,1)
                  RDATA(1)=RCQS_SPATIAL(OFFSET1,nq)
                ENDIF
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  DO nqq=1,NQET(NQS(ne))
                    nq=NQNE(ne,nqq)
                    !only want to set the values for grid points with
                    !the correct variant
                    IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                      VALUESET = .TRUE.
                      RCQS_SPATIAL(OFFSET1,nq)=RDATA(1)
                    ENDIF
                  ENDDO
C                 Apply to all elements in the group
                  DO n=2,NELIST(0)
                    ne=NELIST(n)
                    DO nqq=1,NQET(NQS(ne))
                      nq=NQNE(ne,nqq)
                      !only want to set the values for grid points with
                      !the correct variant
                      IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                        VALUESET = .TRUE.
                        RCQS_SPATIAL(OFFSET1,nq)=RDATA(1)
                      ENDIF
                    ENDDO
                  ENDDO !n
                ENDIF
              ELSE IF(VARTYPE.EQ.3) THEN !YQS values
                RDEFLT(1)=YQS(OFFSET+niqs,1)
                WRITE(CHAR12,'(D12.5)') RDEFLT(1)
                FORMAT='($,'' The value is ['//CHAR12(1:12)
     '            //']: '',D12.5)'
                IF(IOTYPE.EQ.3) THEN
                  ! Need to find a grid point with the right variant
                  DO n=1,NELIST(0)
                    ne=NELIST(n)
                    DO nqq=1,NQET(NQS(ne))
                      nq=NQNE(ne,nqq)
                      IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                        RDATA(1)=YQS(OFFSET+niqs,nq)
                        GOTO 4321
                      ENDIF
                    ENDDO
                  ENDDO
 4321             CONTINUE
                ENDIF
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  DO nqq=1,NQET(NQS(ne))
                    nq=NQNE(ne,nqq)
                    !Since YQS is stored at all grid points we need to
                    !only overwrite the values for the correct variant
                    IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                      VALUESET = .TRUE.
                      YQS(OFFSET+niqs,nq)=RDATA(1)
                    ENDIF
                  ENDDO
C                 Apply to all elements in the group
                  DO n=2,NELIST(0)
                    ne=NELIST(n)
                    DO nqq=1,NQET(NQS(ne))
                      nq=NQNE(ne,nqq)
                      !Since YQS is stored at all grid points we need to
                      !only overwrite the values for the correct variant
                      IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                        VALUESET = .TRUE.
                        YQS(OFFSET+niqs,nq)=RDATA(1)
                      ENDIF
                    ENDDO
                  ENDDO !n
                ENDIF
              ELSE
                ERROR='>>Invalid variable type'
                GOTO 9999
              ENDIF
              GOTO 6100 !For more elements
            ENDIF

          ELSE IF(SPATIAL(OFFSET+niqs,nqv).EQ.2.OR.
     '        SPATIAL(OFFSET+niqs,nqv).LT.0) THEN !Defined by nodes
            nmq=nmq+1
            IF(nmq.GT.NMQT) NMQT=nmq
            CALL ASSERT(NMQT.LE.NMQM,'>>Increase NMQM',ERROR,*9999)
            WRITE(CHAR3,'(I3)') niqs
            FORMAT='($,'' Enter linear basis type number [1]: '',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=IABS(SPATIAL(OFFSET+niqs,nqv))
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '        IONE,1,NBFM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              SPATIAL(OFFSET+niqs,nqv)=-IDATA(1)
            ENDIF
            nb=-SPATIAL(OFFSET+niqs,nqv)
            LINEARBASIS=.FALSE.
            IF(NIT(nb).EQ.1) THEN
              IF((IBT(1,1,nb).EQ.1).AND.(IBT(2,1,nb).EQ.1))
     '          LINEARBASIS=.TRUE.
            ELSE IF(NIT(nb).EQ.2) THEN
              IF((IBT(1,1,nb).EQ.1).AND.(IBT(2,1,nb).EQ.1).AND.
     '          (IBT(1,2,nb).EQ.1).AND.(IBT(2,2,nb).GE.1))
     '          LINEARBASIS=.TRUE.
            ELSE IF(NIT(nb).EQ.3) THEN
              IF((IBT(1,1,nb).EQ.1).AND.(IBT(2,1,nb).EQ.1).AND.
     '          (IBT(1,2,nb).EQ.1).AND.(IBT(2,2,nb).GE.1).AND.
     '          (IBT(1,3,nb).EQ.1).AND.(IBT(2,3,nb).GE.1))
     '          LINEARBASIS=.TRUE.
            ENDIF
            CALL ASSERT(LINEARBASIS,'>>Invalid linear basis',
     '        ERROR,*9999)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(VARTYPE.EQ.1) THEN !ICQS values
                IF(nonode.EQ.1) THEN
                  IDEFLT(1)=CELL_ICQS_VALUE(OFFSET+niqs,nqv)
                ELSE
                  IDEFLT(1)=INT(CELL_CP(nmq,NPNODE(nonode-1,nr)))
                ENDIF
                WRITE(CHAR12,'(I5)') IDEFLT(1)
                WRITE(CHAR5,'(I5)') np
                CALL STRING_TRIM(CHAR5,IBEG1,IEND1)
                FORMAT='($,'' The value at node '
     '            //CHAR5(IBEG1:IEND1)//' is ['//CHAR12(1:5)
     '            //']: '',I10)'
                IF(IOTYPE.EQ.3) IDATA(1)=INT(CELL_CP(nmq,np))
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,-IMAX,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) CELL_CP(nmq,np)=DBLE(IDATA(1))
              ELSE IF(VARTYPE.EQ.2) THEN !RCQS values
                IF(nonode.EQ.1) THEN
                  RDEFLT(1)=CELL_RCQS_VALUE(OFFSET+niqs,nqv)
                ELSE
                  RDEFLT(1)=CELL_CP(nmq,NPNODE(nonode-1,nr))
                ENDIF
                WRITE(CHAR12,'(D12.5)') RDEFLT(1)
                WRITE(CHAR5,'(I5)') np
                CALL STRING_TRIM(CHAR5,IBEG1,IEND1)
                FORMAT='($,'' The value at node '
     '            //CHAR5(IBEG1:IEND1)//' is ['//CHAR12//']: '',D12.5)'
                IF(IOTYPE.EQ.3) RDATA(1)=CELL_CP(nmq,np)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) CELL_CP(nmq,np)=RDATA(1)
              ELSE IF(VARTYPE.EQ.3) THEN !YQS values
                IF(nonode.EQ.1) THEN
                  RDEFLT(1)=YQS(OFFSET+niqs,1)
                ELSE
                  RDEFLT(1)=CELL_CP(nmq,NPNODE(nonode-1,nr))
                ENDIF
                WRITE(CHAR12,'(D12.5)') RDEFLT(1)
                WRITE(CHAR5,'(I5)') np
                CALL STRING_TRIM(CHAR5,IBEG1,IEND1)
                FORMAT='($,'' The value at node '
     '            //CHAR5(IBEG1:IEND1)//' is ['//CHAR12//']: '',D12.5)'
                IF(IOTYPE.EQ.3) RDATA(1)=CELL_CP(nmq,np)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) CELL_CP(nmq,np)=RDATA(1)
              ELSE
                ERROR='>>Invalid variable type'
                GOTO 9999
              ENDIF
            ENDDO !nonode
            IF(IOTYPE.NE.3) THEN
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                ns=0
                DO nn=1,NNT(nb)
                  ns=ns+1
                  XE(ns,1)=CELL_CP(nmq,NPNE(nn,nb,ne))
                ENDDO !nn
                nqsc=NQS(ne)
                II=MAX(1,NQXI(1,nqsc))
                IJ=1
                IK=1
                IF(NQXI(0,nqsc).GT.1) IJ=MAX(1,NQXI(2,nqsc))
                IF(NQXI(0,nqsc).GT.2) IK=MAX(1,NQXI(3,nqsc))
                DO i=1,3
                  XI(i)=0.0d0
                ENDDO !i
C               Loop over the grid points in each element
                DO nik=1,IK
                  DO nij=1,IJ
                    DO nii=1,II
                      neq=nii+((nij-1)*NQXI(1,nqsc))
                      IF(NQXI(0,nqsc).GT.1) neq=neq+((nik-1)*
     '                  NQXI(1,nqsc)*NQXI(2,nqsc))
                      nq=NQNE(ne,neq)
C                     Local xi coordinates of grid point nq in
C                     element ne
                      IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                      IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                      IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)

                      IF(VARTYPE.EQ.1) THEN !ICQS values
                        !only want to set the values for grid points
                        !with the correct variant
                        IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                          VALUESET = .TRUE.
                          ICQS_SPATIAL(OFFSET1,nq)=
     '                      INT(PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                      INP(1,1,nb),nb,1,XI,XE(1,1)))
                        ENDIF
                      ELSE IF(VARTYPE.EQ.2) THEN !RCQS values
                        !only want to set the values for grid points
                        !with the correct variant
                        IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                          VALUESET = .TRUE.
                          RCQS_SPATIAL(OFFSET1,nq)=PXI(IBT(1,1,nb),
     '                      IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,XE(1,1))
                        ENDIF
                      ELSE IF(VARTYPE.EQ.3) THEN !YQS values
                        !Since YQS is stored at all grid points we need
                        !to only overwrite the values for the correct
                        !variant
                        IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                          VALUESET = .TRUE.
                          YQS(OFFSET+niqs,nq)=PXI(IBT(1,1,nb),
     '                      IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,XE(1,1))
                        ENDIF
                      ENDIF
                    ENDDO !ni1
                  ENDDO !ni2
                ENDDO !ni3
              ENDDO !noelem
            ENDIF

          ELSE IF(SPATIAL(OFFSET+niqs,nqv).EQ.3) THEN !defined b grid pt
            WRITE(CHAR3,'(I3)') niqs
            nqq=NQR(1,nr)
 6300       FORMAT='($,'' Enter collocation point #s/name '
     '        //'[EXIT]: '',I5)'
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
 6400       CDATA(1)='GRIDS' !for use with group input
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,
     '        0,NQT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
            IF(NQLIST(1).NE.0) THEN !not default exit
              IF(IOTYPE.NE.3) THEN
                DO n=1,NQLIST(0)
                  nq=NQLIST(n)
                  IF(nq.LT.NQR(1,nr).OR.nq.GT.NQR(2,nr)) THEN
                    WRITE(OP_STRING,'('' >>Grid point '',I7,'' is not '
     '                //'in the current region'')') nq
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    GOTO 6400
                  ENDIF
                ENDDO !n
              ENDIF !IOTYPE.NE.3
C             Define parameter for first grid point in group. Rest of
C             group filled in at the end.
              nq=NQLIST(1)
              IF(VARTYPE.EQ.1) THEN !ICQS values
                IDEFLT(1)=CELL_ICQS_VALUE(OFFSET+niqs,nqv)
                WRITE(CHAR5,'(I5)') IDEFLT(1)
                FORMAT='($,'' The value is ['//CHAR5(1:5)//']: '',I10)'
                IF(IOTYPE.EQ.3) IDATA(1)=ICQS_SPATIAL(OFFSET1,nq)
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  !only want to set the values for grid points with
                  !the correct variant
                  IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                    VALUESET = .TRUE.
                    ICQS_SPATIAL(OFFSET1,nq)=IDATA(1)
                  ENDIF
C                 Apply to all grid points in the group
                  DO n=2,NQLIST(0)
                    nq=NQLIST(n)
                    !only want to set the values for grid points with
                    !the correct variant
                    IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                      VALUESET = .TRUE.
                      ICQS_SPATIAL(OFFSET1,nq)=IDATA(1)
                    ENDIF
                  ENDDO !n
                ENDIF
              ELSE IF(VARTYPE.EQ.2) THEN !RCQS values
                RDEFLT(1)=CELL_RCQS_VALUE(OFFSET+niqs,nqv)
                WRITE(CHAR12,'(D12.5)') RDEFLT(1)
                FORMAT='($,'' The value is ['//CHAR12(1:12)
     '            //']: '',D12.5)'
                IF(IOTYPE.EQ.3) RDATA(1)=RCQS_SPATIAL(OFFSET1,nq)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  !only want to set the values for grid points with
                  !the correct variant
                  IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                    VALUESET = .TRUE.
                    RCQS_SPATIAL(OFFSET1,nq)=RDATA(1)
                  ENDIF
C                 Apply to all elements in the group
                  DO n=2,NQLIST(0)
                    nq=NQLIST(n)
                    !only want to set the values for grid points with
                    !the correct variant
                    IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                      VALUESET = .TRUE.
                      RCQS_SPATIAL(OFFSET1,nq)=RDATA(1)
                    ENDIF
                  ENDDO !n
                ENDIF
              ELSE IF(VARTYPE.EQ.3) THEN !YQS values
                RDEFLT(1)=YQS(OFFSET+niqs,1)
                WRITE(CHAR12,'(D12.5)') RDEFLT(1)
                FORMAT='($,'' The value is ['//CHAR12(1:12)
     '            //']: '',D12.5)'
                IF(IOTYPE.EQ.3) RDATA(1)=YQS(OFFSET+niqs,nq)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  !Since YQS is stored at all grid points we need to
                  !only overwrite the values for the correct variant
                  IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                    VALUESET = .TRUE.
                    YQS(OFFSET+niqs,nq)=RDATA(1)
                  ENDIF
C                 Apply to all elements in the group
                  DO n=2,NQLIST(0)
                    nq=NQLIST(n)
                    !Since YQS is stored at all grid points we need to
                    !only overwrite the values for the correct variant
                    IF (ICQS_SPATIAL(1,nq).EQ.nqv) THEN
                      VALUESET = .TRUE.
                      YQS(OFFSET+niqs,nq)=RDATA(1)
                    ENDIF
                  ENDDO !n
                ENDIF
              ELSE
                ERROR='>>Invalid variable type'
                GOTO 9999
              ENDIF
              GOTO 6300 !For more grid points
            ENDIF
          ENDIF
C *** DPN 06-06-06 - need to make this check a bit smarter so
C ***                that when you exit it passes through the check.
c          CALL ASSERT(VALUESET,
c     '      '>>No grid points with the correct variant were given',
c     '      ERROR,*9999)
        ENDIF ! ANYSPATIAL
      ENDDO !niqs

      CALL EXITS('IPMAT3_CELL1')
      RETURN
 9999 CALL ERRORS('IPMAT3_CELL1',ERROR)
      CALL EXITS('IPMAT3_CELL1')
      RETURN 1
      END



