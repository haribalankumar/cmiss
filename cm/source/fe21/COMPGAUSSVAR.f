      SUBROUTINE COMPGAUSSVAR(IBT,IDO,INP,NEELEM,NELIST,NGLIST,NQET,
     '  NQLIST,NQNE,NQS,NQSCNB,NQXI,NRLIST,PGNQE,XIG,YG,YQS,ERROR,*)

C#### Subroutine: COMPGAUSSVAR
C###  Description:
C###    <HTML>
C###    COMPGAUSSVAR compares Gauss point array (YG) values to
C###    (interpolated) grid point variables (YQS) (a corresponding
C###    order of the specified the grid point and Gauss point indices
C###    is assumed) and returns (a list of) the maximum absolute
C###    differences to a user defined variable.
C###    </HTML>
C*** Created by Martyn Nash June 2000

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NGLIST(0:NGM),NQET(NQSCM),
     '  NQLIST(0:NQM),NQNE(NEQM,NQEM),NQS(NEQM),NQSCNB(NQSCM),
     '  NQXI(0:NIM,NQSCM),NRLIST(0:NRM)
      REAL*8 PGNQE(NGM,NQEM,NQSCM),XIG(NIM,NGM,NBM),YG(NIYGM,NGM,NEM),
     '  YQS(NIQSM,NQM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER ERR,i,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,N3CO,ne,ne1,
     '  NEXTOPT,ng,nig,nigg,niqs,noelem,noelem1,no_nrlist,nqq,nr,SCHEME
      REAL*8 RETVALUE,SUM
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,NE_EXCLUDE,RETVAR
      CHARACTER COMMAND_TEMP*200,RETTYPE*20,RETVARNAME*20,TYPE*10

      CALL ENTERS('COMPGAUSSVAR',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(CO(1),IBEG1,IEND1)
        COMMAND_TEMP=CO(1)(IBEG1:IEND1)
        DO i=2,noco
          CALL STRING_TRIM(COMMAND_TEMP,IBEG,IEND)
          CALL STRING_TRIM(CO(i),IBEG1,IEND1)
          COMMAND_TEMP=COMMAND_TEMP(IBEG:IEND)//' '//CO(i)(IBEG1:IEND1)
        ENDDO
        CALL STRING_TRIM(COMMAND_TEMP,IBEG,IEND)

C---------------------------------------------------------------------
C#### Command: FEM compare gaussvars gridvars
C###  Parameter:      <yqs #s[1]>
C###    Specifies the list of grid point variable indices.
C###  Parameter:      <yg #s[1]>
C###    Specifies the list of gauss point location indices.
C###  Parameter:      <exclude element (GROUP/#s)[all]>
C###    Specify the element numbers to exclude (eg. infarcted elements).
C###  Parameter:      <return <maxabsdiff[maxabsdiff]> to_variable NAME>
C###    Specify the value to return and the user variable name to use.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description:
C###    Compares Gauss point array (YG) values to (interpolated)
C###    grid point variables (YQS) (a corresponding order of
C###    the specified the grid point and Gauss point indices is assumed)
C###    and returns (a list of) the maximum absolute differences.

        OP_STRING(1)=COMMAND_TEMP(IBEG:IEND)//' gridvars'
        OP_STRING(2)=BLANK(1:15)//'<yqs #s[1]>'
        OP_STRING(3)=BLANK(1:15)//'<yg #s[1]>'
        OP_STRING(4)=BLANK(1:15)//'<exclude element (GROUP/#s)[all]>'
        OP_STRING(5)=BLANK(1:15)//'<return <maxabsdiff[maxabsdiff]>'
     '    //' to_variable NAME>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','COMPGAUSSVAR',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'GRIDVARS',4,noco+1,NTCO,N3CO)) THEN
          TYPE='GRIDVARS'
        ELSE
          noco=3
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        IF(TYPE(1:8).EQ.'GRIDVARS') THEN
C         Set the grid variable indices
          IF(CBBREV(CO,'YQS',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NQM,NQLIST(0),NQLIST(1),ERROR,*9999)
          ELSE
            NQLIST(0)=1
            NQLIST(1)=1
          ENDIF

C         Set the Gauss variable indices
          IF(CBBREV(CO,'YG',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NGM,NGLIST(0),NGLIST(1),ERROR,*9999)
          ELSE
            NGLIST(0)=1
            NGLIST(1)=1
          ENDIF
          CALL ASSERT(NQLIST(0).EQ.NGLIST(0),
     '      '>>Must be the same number of YG and YQS indices',
     '      ERROR,*9999)

C *** DPN 28 March 2001 - Check the variable indices
          DO nigg=1,NQLIST(0)
            niqs=NQLIST(nigg)
            nig=NGLIST(nigg)
            CALL ASSERT(niqs.LE.NIQSM,'>>Need to increase NIQSM',
     '        ERROR,*9999)
            CALL ASSERT(nig.LE.NIYGM,'>>Need to increase NIYGM',
     '        ERROR,*9999)
          ENDDO

C         Set the elements to exclude from the comparison
          IF(CBBREV(CO,'EXCLUDE',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999) !exclude elements in list
          ELSE
            NELIST(0)=0 !include all elements
          ENDIF !exclude

        ENDIF !type=gridvars

C       Set the return value type and variable name
        IF(CBBREV(CO,'RETURN',3,noco+1,NTCO,N3CO)) THEN
          RETVAR=.TRUE.
          NEXTOPT=N3CO+1
C         Set the type of variable to return (only maxabsdiff for now)
          IF(ABBREV(CO(NEXTOPT),'MAXABSDIFF',3)) THEN
            RETTYPE='MAXABSDIFF'
            NEXTOPT=NEXTOPT+1
          ELSE IF(ABBREV(CO(NEXTOPT),'TO_VARIABLE',2)) THEN
            RETTYPE='MAXABSDIFF'
          ELSE IF(NEXTOPT.LE.NTCO) THEN
            noco=3
            CO(noco+1)='?'
            GO TO 1
          ENDIF
C         Set the name of the variable to which the value is returned
          IF(ABBREV(CO(NEXTOPT),'TO_VARIABLE',2)) THEN
            CALL STRING_TRIM(CO(NEXTOPT+1),IBEG1,IEND1)
            RETVARNAME=CO(NEXTOPT+1)(IBEG1:IEND1)
            CALL ASSERT(IBEG1.LE.IEND1,'>>Must supply a return user '
     '        //'variable name',ERROR,*9999)
          ELSE
            noco=3
            CO(noco+1)='?'
            GO TO 1
          ENDIF
        ELSE
          RETVAR=.FALSE.
          RETTYPE='MAXABSDIFF'
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

C       Compute the selected return value
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)

          IF(TYPE(1:8).EQ.'GRIDVARS') THEN

            IF(.NOT.CALL_GMAP) THEN
              CALL GEN_GRID_MAP(IBT,IDO,INP,NQET,NQSCNB,NQXI,PGNQE,XIG,
     '          ERROR,*9999)
              CALL_GMAP=.TRUE.
            ENDIF

            RETVALUE=0.d0
            DO nigg=1,NGLIST(0)
              nig=NGLIST(nigg)
              niqs=NQLIST(nigg)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                NE_EXCLUDE=.FALSE.
                DO noelem1=1,NELIST(0)
                  ne1=NELIST(noelem1)
                  IF(ne1.EQ.ne) NE_EXCLUDE=.TRUE. !exclude this element
                ENDDO !noelem1 (ne1)
                IF(.NOT.NE_EXCLUDE) THEN
                  SCHEME=NQS(ne)
                  DO ng=1,NGT(NQSCNB(SCHEME))
                    SUM=0.0d0
                    DO nqq=1,NQET(SCHEME)
                      SUM=SUM+
     '                  (YQS(niqs,NQNE(ne,nqq))*PGNQE(ng,nqq,SCHEME))
                    ENDDO !nqq
C                   Compute return value (NB. other types could go here)
                    IF(RETTYPE(1:10).EQ.'MAXABSDIFF') THEN
                      RETVALUE=MAX(RETVALUE,DABS(SUM-YG(nig,ng,ne)))
                    ENDIF !rettype
                  ENDDO !ng
                ENDIF !not ne_exclude
              ENDDO !noelem (ne)
            ENDDO !nigg (nig)
          ENDIF !TYPE

        ENDDO !no_nrlist

        IF(RETVAR) THEN
C         Assign the return value to the given user variable
          CALL STRING_TRIM(RETVARNAME,IBEG1,IEND1)
          CALL SET_USER_DOUBLE(RETVARNAME(IBEG1:IEND1),RETVALUE,ERR)
          IF(ERR.NE.0) THEN
            ERROR=
     '        'Unable to set user var "'//RETVARNAME(IBEG2:IEND2)//'"'
            GOTO 9999
          ENDIF
        ELSE !not RETVAR
C         Output the desired value to screen
          IF(RETTYPE(1:10).EQ.'MAXABSDIFF') THEN
            WRITE(OP_STRING,
     '        '('' The maximum absolute difference is: '',F15.7)')
     '        RETVALUE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF !RETVAR

      ENDIF

      CALL EXITS('COMPGAUSSVAR')
      RETURN
 9999 CALL ERRORS('COMPGAUSSVAR',ERROR)
      CALL EXITS('COMPGAUSSVAR')
      RETURN 1
      END


