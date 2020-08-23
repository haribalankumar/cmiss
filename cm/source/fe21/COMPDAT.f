      SUBROUTINE COMPDAT(NDDATA,NDP,NRLIST,WD,Z_CONT,ZD,ZD2,STRING,
     &  ERROR,*)

C#### Subroutine: COMPDAT
C###  Description:
C###    <HTML>
C###    COMPDAT compares different field values in data files. Displays
C###    a variety of metrics to show the differences. The fields
C###    are stored in the same data file.
C###    <BR> Note: Data weights (WD) are not used.
C###    </HTML>

C*** Created By Leo Cheng 3-JUL-1999

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER NDDATA(0:NDM,0:NRM),NDP(NDM),
     '  NRLIST(0:NRM)
      REAL*8 WD(NJM,NDM),Z_CONT(NDM,2,67),ZD(NJM,NDM),ZD2(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER COUNT,COMP_FIELD,ERR,IBEG,IBEG1,IEND,IEND1,MAST_FIELD,
     '  N3CO,nd,ND_MAX,ND_MIN,nj,NJ_COMP,NJTT,NJ_MAST,nr
      REAL*8 COMPARE,COMPARESUM,COMPARESQUAREDSUM,DIFFSQUAREDSUM,
     '  MASTER,MASTERSUM,MASTERSQUAREDSUM,MAXDIFF,MINDIFF,
     '  PRODUCTSUM,RMS,RELRMS,
     '  SI,TMP,tmp1
      LOGICAL ALL_REGIONS,CBBREV,GEOMETRY,OPFILE,RMS_VAR,RELRMS_VAR,
     &  SI_VAR
      CHARACTER ANGLE_TYPE*7,COMPAREFILE*(MXCH),
     &  MASTERFILE*(MXCH),OUTPUTFILE*(MXCH),TITLE*80,TYPE*20,
     &  RMS_VAR_NAME*(MXCH),RELRMS_VAR_NAME*(MXCH),SI_VAR_NAME*(MXCH)

!     Function
      INTEGER IFROMC

      CALL ENTERS('COMPDAT',*9999)

      CALL STRING_TRIM(FILE00,IBEG1,IEND1)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM compare data;FILENAME field
C###  Description:
C###    Compares the field values associated with data points.
C###  Parameter:      <master_field #[1]>
C###    The master field number.
C###  Parameter:      <compare_field #[1]>
C###    The compare field number.
C###  Parameter:      <region (#s/all)[1]>
C###    The region number which the data/fields correspond to.
C###  Parameter:      <rms_variable NAME>
C###    If present, gives the name of the interpreter variable to
C###    set to the value of the RMS difference.
C###  Parameter:      <relrms_variable NAME>
C###    If present, gives the name of the interpreter variable to
C###    set to the value of the relative RMS difference.
C###  Parameter:      <si_variable NAME>
C###    If present, gives the name of the interpreter variable to
C###    set to the value of the similarity index.

        OP_STRING(1)=STRING(1:IEND)//';FILENAME field'
        OP_STRING(2)=BLANK(1:15)//'<master_field #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<compare_field #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<rms_variable NAME>'
        OP_STRING(5)=BLANK(1:15)//'<relrms_variable NAME>'
        OP_STRING(6)=BLANK(1:15)//'<si_variable NAME>'
        OP_STRING(7)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM compare data;FILENAME geometry
C###  Description:
C###    Compares the data geometric coordinates.

C###  Parameter:      <masterfile FILENAME[$current]>
C###    The master data file.
C###  Parameter:      <comparefile FILENAME[$current]>
C###    The compare data file.
C###  Parameter:      <region (#s/all)[1]>
C###    The region number which the data/fields correspond to.

        OP_STRING(1)=STRING(1:IEND)//';FILENAME geometry'
        OP_STRING(2)=BLANK(1:15)//
     '    '<masterfile FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//
     '    '<comparefile FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------


      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','COMPDAT',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          OUTPUTFILE=COQU(noco,1)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL STRING_TRIM(OUTPUTFILE,IBEG1,IEND1)
          CALL OPENF(IOFI,'DISK',OUTPUTFILE(IBEG1:IEND1)//'.opcomp',
     '      'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL ASSERT(NRLIST(0).EQ.1,
     '    '>> Only implemented for 1 region',ERROR,*9999)
        nr=NRLIST(1)



        GEOMETRY=.FALSE. !defaults field is true
        IF(CBBREV(CO,'GEOMETRY',3,noco+1,NTCO,N3CO)) THEN
          GEOMETRY=.TRUE.
        ENDIF




C*** 26-SEP-2000 LKC adding GEOMETRY comparisons

        IF(GEOMETRY) THEN

          IF(CBBREV(CO,'MASTERFILE',2,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            MASTERFILE=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            MASTERFILE=FILE00(IBEG1:IEND1)
          ENDIF
          IF(CBBREV(CO,'COMPAREFILE',2,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            COMPAREFILE=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            COMPAREFILE=FILE00(IBEG1:IEND1)
          ENDIF


C LKC note: there is currently
C  - no check that the data points are compatible
C  - ie NDT, NJT, WD, NDP etc.
C  - Does not handle fibres or sheets.
C
          NJTT=NJT
          ANGLE_TYPE='-' !initialise angle type
          TYPE='-' !initialise type

C*** Read in the master data
          CALL STRING_TRIM(MASTERFILE,IBEG,IEND)
          CALL OPENF(IOFILE1,'DISK',MASTERFILE(IBEG:IEND)//'.ipdata',
     '      'OLD','SEQUEN','FORMATTED',132,ERROR,*9999)
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C          CALL IODATA('READ',ANGLE_TYPE,TYPE,IOFILE1,NDP,NJTT,
C     '      TITLE,WD,ZD,ERROR,*9999)
          CALL IODATA('READ',ANGLE_TYPE,TYPE,IOFILE1,NDP,NJTT,
     &      TITLE,WD,Z_CONT,ZD,ERROR,*9999)
          CALL CLOSEF(IOFILE1,ERROR,*9999)

C*** Read in the comparison data
          CALL STRING_TRIM(COMPAREFILE,IBEG,IEND)
          CALL OPENF(IOFILE2,'DISK',COMPAREFILE(IBEG:IEND)//'.ipdata',
     '      'OLD','SEQUEN','FORMATTED',132,ERROR,*9999)
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C          CALL IODATA('READ',ANGLE_TYPE,TYPE,IOFILE2,NDP,NJTT,
C     '      TITLE,WD,ZD2,ERROR,*9999)
          CALL IODATA('READ',ANGLE_TYPE,TYPE,IOFILE2,NDP,NJTT,
     &      TITLE,WD,Z_CONT,ZD2,ERROR,*9999)
          CALL CLOSEF(IOFILE2,ERROR,*9999)

          ND_MAX=0 !initialise the data points
          ND_MIN=0

          MAXDIFF=-RMAX !initialise max and min displacements
          MINDIFF=RMAX
          COUNT=0
          RMS=0.D0
          DO nd=1,NDT
            COUNT=COUNT+1 !currently COUNT should == NDT
            DIFFSQUAREDSUM=0.D0
            DO nj=1,NJT
              DIFFSQUAREDSUM=DIFFSQUAREDSUM+(ZD(nj,nd)-ZD2(nj,nd))**2
            ENDDO !nj
            TMP=SQRT(DIFFSQUAREDSUM)

            IF(TMP.LT.MINDIFF) THEN ! set max and min differences
              MINDIFF=TMP
              ND_MIN=COUNT
            ENDIF
            IF(TMP.GT.MAXDIFF) THEN
              MAXDIFF=TMP
              ND_MAX=COUNT
            ENDIF

            RMS=RMS+DIFFSQUAREDSUM
          ENDDO !nd
          RMS=SQRT(RMS/DBLE(COUNT))

          WRITE(OP_STRING(1),'('' '')')
          WRITE(OP_STRING(2),'('' Geometric data comparison'')')
          WRITE(OP_STRING(3),'('' '')')
          WRITE(OP_STRING(4),'('' RMS displacement     :'',F12.3)')
     '      RMS
          WRITE(OP_STRING(5),'('' Maximum displacement :'',F12.3)')
     '      MAXDIFF
          WRITE(OP_STRING(6),'(''   at Electrode       :'',I12)')
     '      ND_MAX
          WRITE(OP_STRING(7),'('' Minimum displacement :'',F12.3)')
     '      MINDIFF
          WRITE(OP_STRING(8),'(''   at Electrode       :'',I12)')
     '      ND_MIN
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)



C*** Field comparisons ....
        ELSE

C*** Get field number to compare
          IF(CBBREV(CO,'MASTER_FIELD',3,noco+1,NTCO,N3CO)) THEN
            MAST_FIELD=IFROMC(CO(N3CO+1))
          ELSE
            MAST_FIELD=1
          ENDIF
          IF(CBBREV(CO,'COMPARE_FIELD',3,noco+1,NTCO,N3CO)) THEN
            COMP_FIELD=IFROMC(CO(N3CO+1))
          ELSE
            COMP_FIELD=1
          ENDIF

          CALL ASSERT(MAST_FIELD.LE.NJ_LOC(NJL_FIEL,0,nr)
     '      ,'>> Non-existing MASTER field',ERROR,*9999)
          CALL ASSERT(COMP_FIELD.LE.NJ_LOC(NJL_FIEL,0,nr)
     '      ,'>> Non-existing COMPARE field',ERROR,*9999)

C*** Set nj values corresponding to fields numbers
          NJ_MAST=NJ_LOC(NJL_FIEL,MAST_FIELD,nr)
          NJ_COMP=NJ_LOC(NJL_FIEL,COMP_FIELD,nr)
          CALL ASSERT(NJ_MAST.GT.0,'>> Invalid MAST field',ERROR,*9999)
          CALL ASSERT(NJ_COMP.GT.0,'>> Invalid COMP field',ERROR,*9999)

          IF(CBBREV(CO,'RMS_VARIABLE',3,noco+2,NTCO,N3CO)) THEN
            RMS_VAR=.TRUE.
C           Set the name of the variable to which the value is returned
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            RMS_VAR_NAME=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            RMS_VAR=.FALSE.
          ENDIF

          IF(CBBREV(CO,'RELRMS_VARIABLE',3,noco+2,NTCO,N3CO)) THEN
            RELRMS_VAR=.TRUE.
C           Set the name of the variable to which the value is returned
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            RELRMS_VAR_NAME=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            RELRMS_VAR=.FALSE.
          ENDIF

          IF(CBBREV(CO,'SI_VARIABLE',2,noco+2,NTCO,N3CO)) THEN
            SI_VAR=.TRUE.
C           Set the name of the variable to which the value is returned
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            SI_VAR_NAME=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            SI_VAR=.FALSE.
          ENDIF

          OP_STRING(1)=' '
          WRITE(OP_STRING(2),'('' Comparison metrics :'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)


C*** Calculates metrics ....
          MASTERSUM=0.0d0
          MASTERSQUAREDSUM=0.0d0
          COMPARESUM=0.0d0
          COMPARESQUAREDSUM=0.0d0
          PRODUCTSUM=0.0d0
          DIFFSQUAREDSUM=0.0d0


          CALL ASSERT(NDDATA(0,nr).GT.0,
     '      '>> No data values to compare',ERROR,*9999)


          COUNT=0 !not strictly required by may be useful for subsets latter
          DO nd=1,NDDATA(0,nr)
            COUNT=COUNT+1
            MASTER=ZD(NJ_MAST,nd)
            COMPARE=ZD(NJ_COMP,nd)
            MASTERSUM=MASTERSUM+MASTER
            COMPARESUM=COMPARESUM+COMPARE
            MASTERSQUAREDSUM=MASTERSQUAREDSUM+MASTER**2
            COMPARESQUAREDSUM=COMPARESQUAREDSUM+COMPARE**2
            PRODUCTSUM=PRODUCTSUM+MASTER*COMPARE
            DIFFSQUAREDSUM=DIFFSQUAREDSUM+(MASTER-COMPARE)**2
          ENDDO !nd
          RMS=DSQRT(DIFFSQUAREDSUM/DBLE(COUNT))
          IF(MASTERSQUAREDSUM.GT.ZERO_TOL) THEN
            RELRMS=DSQRT(DIFFSQUAREDSUM/MASTERSQUAREDSUM)
          ELSE
            RELRMS=0.0d0
          ENDIF

          ! WHY ?
cc          CALL ASSERT(DABS(COMPARESUM).GT.ZERO_TOL,
cc     '      '>> No field in Compare FIELD',ERROR,*9999)
cc          CALL ASSERT(DABS(MASTERSUM).GT.ZERO_TOL,
cc     '      '>> No field in Master FIELD',ERROR,*9999)

          IF(NDDATA(0,nr).GT.1) THEN
            tmp1=(MASTERSQUAREDSUM-MASTERSUM**2/DBLE(COUNT))*
     '        (COMPARESQUAREDSUM-COMPARESUM**2/DBLE(COUNT))
            IF(tmp1.GT.ZERO_TOL) THEN
              SI=(PRODUCTSUM-(MASTERSUM*COMPARESUM/DBLE(COUNT)))/
     '          DSQRT(tmp1)
            ELSE
              SI=0.D0
            ENDIF
          ELSE
            SI=0.D0
          ENDIF

C*** Output results
          OP_STRING(1)='  '
          WRITE(OP_STRING(2),'('' Region #         = '',I12)')
     '      nr
          WRITE(OP_STRING(3),'('' Master field #   = '',I12)')
     '      MAST_FIELD
          WRITE(OP_STRING(4),'(''   with nj #      = '',I12)')
     '      NJ_MAST
          WRITE(OP_STRING(5),'('' Compare field #  = '',I12)')
     '      COMP_FIELD
          WRITE(OP_STRING(6),'(''   with nj #      = '',I12)')
     '      NJ_COMP
          WRITE(OP_STRING(7),'('' # Data points    = '',I12)')
     '      NDDATA(0,nr)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          OP_STRING(1)='  '
          WRITE(OP_STRING(2),'('' RMS              = '',D12.5)') RMS
          WRITE(OP_STRING(3),'('' Relative RMS     = '',D12.5)') RELRMS
          WRITE(OP_STRING(4),'('' Similarity Index = '',D12.5)') SI
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          ! If required, set the interpreter variables
          IF(RMS_VAR) THEN
            CALL STRING_TRIM(RMS_VAR_NAME,IBEG1,IEND1)
            CALL SET_USER_DOUBLE(RMS_VAR_NAME(IBEG1:IEND1),RMS,ERR)
            IF(ERR.NE.0) THEN
              ERROR=
     '          'Unable to set user var "'//RMS_VAR_NAME(IBEG1:IEND1)//
     '          '"'
              GOTO 9999
            ENDIF
          ENDIF
          IF(RELRMS_VAR) THEN
            CALL STRING_TRIM(RELRMS_VAR_NAME,IBEG1,IEND1)
            CALL SET_USER_DOUBLE(RELRMS_VAR_NAME(IBEG1:IEND1),RELRMS,
     '        ERR)
            IF(ERR.NE.0) THEN
              ERROR=
     '          'Unable to set user var "'//RELRMS_VAR_NAME(IBEG1:IEND1)
     '          //'"'
              GOTO 9999
            ENDIF
          ENDIF
          IF(SI_VAR) THEN
            CALL STRING_TRIM(SI_VAR_NAME,IBEG1,IEND1)
            CALL SET_USER_DOUBLE(SI_VAR_NAME(IBEG1:IEND1),SI,ERR)
            IF(ERR.NE.0) THEN
              ERROR=
     '          'Unable to set user var "'//SI_VAR_NAME(IBEG1:IEND1)//
     '          '"'
              GOTO 9999
            ENDIF
          ENDIF

        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('COMPDAT')
      RETURN
 9999 CALL ERRORS('COMPDAT',ERROR)
      CALL EXITS('COMPDAT')
      RETURN 1
      END


