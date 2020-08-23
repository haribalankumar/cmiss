      SUBROUTINE COMPGRIDVAR(NQLIST,NRLIST,YQS,STRING,ERROR,*)

C#### Subroutine: COMPGRIDVAR
C###  Description:
C###    <HTML>
C###      Performs the specified comparison of specified indices of the
C###      YQS array. Could easily be extended to handle the RCQS and
C###      ICQS arrays as well?? Compares all grid points in the
C###      specified region(s) - would we ever want to specify grid
C###      points/elements??
C###    </HTML>
C*** Created by David Nickerson May 2001

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter List
      INTEGER NQLIST(0:NQM),NRLIST(0:NRM)
      REAL*8 YQS(NIQSM,NQM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER MAX_INDICES
      PARAMETER(MAX_INDICES=10)
      INTEGER FIELD_LIST1(0:MAX_INDICES),FIELD_LIST2(0:MAX_INDICES)
      INTEGER ERR,i,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,N3CO,NEXTOPT,nq,
     '  niqs,niqs1,niqs2,no_nrlist,nr
      REAL*8 RETVALUE
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,RETVAR
      CHARACTER COMMAND_TEMP*200,RETTYPE*20,RETVARNAME*20,ARRAY*10

      CALL ENTERS('COMPGRIDVAR',*9999)

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
C#### Command: FEM compare gridvars yqs
C###  Parameter:      <field1 #s[1]>
C###    Specifies the list of grid point variable indices.
C###  Parameter:      <field2 #s[2]>
C###    Specifies the list of grid point variable indices.
C###  Parameter:      <return <maxabsdiff[maxabsdiff]> to_variable NAME>
C###    Specify the value to return and the user variable name to use.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to compare in.
C###  Description:
C###    Compares <field1> to <field2> for all the grid points in the
C###    specified region, returning the <maxabsdiff> to
C###    the specified variable. Can be extended to enable the use of
C###    any cell/grid array and various comparisons.

        OP_STRING(1)=COMMAND_TEMP(IBEG:IEND)//' yqs'
        OP_STRING(2)=BLANK(1:15)//'<field1 #s[1]>'
        OP_STRING(3)=BLANK(1:15)//'<field2 #s[2]>'
        OP_STRING(4)=BLANK(1:15)//'<return <maxabsdiff[maxabsdiff]>'
     '    //' to_variable NAME>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','COMPGRIDVAR',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'YQS',3,noco+1,NTCO,N3CO)) THEN
          ARRAY='YQS'
        ELSE
          noco=3
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        IF(ARRAY(1:3).EQ.'YQS') THEN
C         Set the grid variable indices for field list 1
          IF(CBBREV(CO,'FIELD1',6,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NQM,NQLIST(0),NQLIST(1),ERROR,*9999)
          ELSE
            NQLIST(0)=1
            NQLIST(1)=1
          ENDIF
          CALL ASSERT(NQLIST(0).LE.MAX_INDICES,
     '      '>>Too many indices (hardcoded)',ERROR,*9999)
          FIELD_LIST1(0) = NQLIST(0)
          DO niqs=1,NQLIST(0)
            niqs1 = NQLIST(niqs)
            CALL ASSERT(niqs1.LE.NIQSM,'>>Invalid field 1 index',
     '        ERROR,*9999)
            FIELD_LIST1(niqs) = niqs1
          ENDDO
C         Set the grid variable indices for field list 2
          IF(CBBREV(CO,'FIELD2',6,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NQM,NQLIST(0),NQLIST(1),ERROR,*9999)
          ELSE
            NQLIST(0)=1
            NQLIST(1)=2
          ENDIF
          CALL ASSERT(NQLIST(0).LE.MAX_INDICES,
     '      '>>Too many indices (hardcoded)',ERROR,*9999)
          FIELD_LIST2(0) = NQLIST(0)
          CALL ASSERT(FIELD_LIST1(0).EQ.FIELD_LIST1(0),
     '      '>>Must have the same number of indices for fields 1 and 2',
     '      ERROR,*9999)
          DO niqs=1,NQLIST(0)
            niqs2 = NQLIST(niqs)
            CALL ASSERT(niqs2.LE.NIQSM,'>>Invalid field 2 index',
     '        ERROR,*9999)
            FIELD_LIST2(niqs) = niqs2
          ENDDO

        ENDIF !array=yqs

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
        RETVALUE=0.d0
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          IF(ARRAY(1:3).EQ.'YQS') THEN
            DO niqs=1,FIELD_LIST1(0)
              niqs1 = FIELD_LIST1(niqs)
              niqs2 = FIELD_LIST2(niqs)
              DO nq=NQR(1,nr),NQR(2,nr)
C               Compute return value (NB. other types could go here)
                IF(RETTYPE(1:10).EQ.'MAXABSDIFF') THEN
                  RETVALUE=MAX(RETVALUE,
     '              DABS(YQS(niqs1,nq)-YQS(niqs2,nq)))
                ENDIF !rettype
              ENDDO !nq
            ENDDO !niqs
          ENDIF !ARRAY
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
     '        '(1P,'' The maximum absolute difference is: '',E10.4)')
     '        RETVALUE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF !RETVAR

      ENDIF

      CALL EXITS('COMPGRIDVAR')
      RETURN
 9999 CALL ERRORS('COMPGRIDVAR',ERROR)
      CALL EXITS('COMPGRIDVAR')
      RETURN 1
      END


