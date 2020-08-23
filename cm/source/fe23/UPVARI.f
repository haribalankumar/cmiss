      SUBROUTINE UPVARI(XP,STRING,ERROR,*)

C#### Subroutine: UPVARI
C###  Description:
C###    UPVARI updates a interpretor variable with a value from a field.
C###    This function should be extended to cope with geometry, fibre,
C###    solution, material, grid point, data point, etc...  values.

C DMAL 25 SEP 2003

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
!      INTEGER 
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER ERR,IFROMC,IBEG,IEND,INTVALUE,N3CO,
     &  nj,nk,np,NPLIST(0:NPM),nr
      REAL*8 REALVALUE
      CHARACTER FROM*16,VARNAME*255
      LOGICAL CBBREV

      CALL ENTERS('UPVARI',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------


C#### Command: FEM update variable
C###  Parameter:      <variable>
C###  Specify the name of the interpreter variable
C###  Parameter:      <from geometry x/y/z>
C###  Specify the location of the value as geometry
C###  Parameter:      <from fibre fibre#>
C###  Specify the location of the value as fibre
C###  Parameter:      <from field field#>
C###  Specify the location of the value as field
C###  Parameter:      <from material*>
C###  Specify the location of the value as material
C###  Parameter:      <from solution*>
C###  Specify the location of the value as solution
C###  Parameter:        <node #[1]>
C###  Specify node number
C###  Parameter:         <value/dx/dy/dxdy/dz/dxdz/dydz/dxdydz[value]>
C###  Specify value or defivative
C###  Parameter:         <region #[1]>
C###  Specify the region number
C###  Parameter:      <from group GROUP_LABEL>
C###  Specify the location of the value as a group
C###  Parameter:         <index #>
C###  Specify the 
C###  Parameter:         <length>
C###  Returns the length of the group

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<variable>'
        OP_STRING(3)=BLANK(1:15)
     &    //'<from geometry x/y/z>'
        OP_STRING(4)=BLANK(1:15)
     &    //'<from fibre fibre#>'
        OP_STRING(5)=BLANK(1:15)
     &    //'<from field field#>'
        OP_STRING(6)=BLANK(1:15)
     &    //'<from material*>'
        OP_STRING(7)=BLANK(1:15)
     &    //'<from solution*>'
        OP_STRING(8)=BLANK(1:15)
        OP_STRING(9)=BLANK(1:17)//'<node #[1]>'
        OP_STRING(10)=BLANK(1:19)//'<value/dx/dy/dxdy/'//
     &    'dz/dxdz/dydz/dxdydz[value]>'
        OP_STRING(11)=BLANK(1:17)//'<region #[1]>'
        OP_STRING(12)=BLANK(1:15)
     &    //'<from group GROUP_LABEL>'
        OP_STRING(13)=BLANK(1:17)//'<index #[1]>'
        OP_STRING(14)=BLANK(1:17)//'<length>'
        OP_STRING(15)=BLANK(1:15)
        OP_STRING(16)=BLANK(1:15)//'* To be implemented'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPVARI',ERROR,*9999)
      ELSE
        
        IF(CBBREV(CO,'REGION',2,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF
        
        VARNAME=CO(noco+1)
        IF(CBBREV(CO,'GEOMETRY',3,noco+1,NTCO,N3CO)) THEN
          FROM='XP'
          IF(CBBREV(CO,'X',1,N3CO+1,N3CO+1,N3CO)) THEN
            nj=NJ_LOC(NJL_GEOM,1,nr)
          ELSEIF(CBBREV(CO,'Y',1,N3CO+1,N3CO+1,N3CO)) THEN
            nj=NJ_LOC(NJL_GEOM,2,nr)
          ELSEIF(CBBREV(CO,'Z',1,N3CO+1,N3CO+1,N3CO)) THEN
            nj=NJ_LOC(NJL_GEOM,3,nr)
          ELSE
            ERROR='No geometry (x/y/z) defined'
            GO TO 9999
          ENDIF
        ELSEIF(CBBREV(CO,'FIBRE',3,noco+1,NTCO,N3CO)) THEN
          FROM='XP'
          nj=NJ_LOC(NJL_FIBR,IFROMC(CO(N3CO+1)),nr)
        ELSEIF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
          FROM='XP'
          nj=NJ_LOC(NJL_FIEL,IFROMC(CO(N3CO+1)),nr)
        ELSEIF(CBBREV(CO,'GROUP',3,noco+1,NTCO,N3CO)) THEN
          FROM='GROUP'
          CDATA(1)='NODES'
          CALL PARSILG(NPLIST,NPM,CDATA(1),CO(N3CO+1),ERROR,*9999)
        ELSE
          ERROR='This source not defined'
          GO TO 9999
        ENDIF
        
        IF(FROM(1:2).EQ.'XP')THEN
          IF(CBBREV(CO,'NODE',2,noco+1,NTCO,N3CO)) THEN
            np=IFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'VALUE',2,noco+1,NTCO,N3CO)) THEN
              nk=1
            ELSEIF(CBBREV(CO,'DX',2,noco+1,NTCO,N3CO)) THEN
              nk=2
            ELSEIF(CBBREV(CO,'DY',2,noco+1,NTCO,N3CO)) THEN
              nk=3
            ELSEIF(CBBREV(CO,'DXDY',4,noco+1,NTCO,N3CO)) THEN
              nk=4
            ELSEIF(CBBREV(CO,'DZ',2,noco+1,NTCO,N3CO)) THEN
              nk=5
            ELSEIF(CBBREV(CO,'DXDZ',4,noco+1,NTCO,N3CO)) THEN
              nk=6
            ELSEIF(CBBREV(CO,'DYDZ',4,noco+1,NTCO,N3CO)) THEN
              nk=7
            ELSEIF(CBBREV(CO,'DXDYDZ',6,noco+1,NTCO,N3CO)) THEN
              nk=8
            ELSE
              nk=1
            ENDIF
          ELSE
            ERROR='No node defined'
            GO TO 9999
          ENDIF
        ELSEIF(FROM(1:5).EQ.'GROUP')THEN
          IF(CBBREV(CO,'INDEX',3,noco+1,NTCO,N3CO)) THEN
            np=IFROMC(CO(N3CO+1))
            IF(np.GT.NPLIST(0))THEN
              ERROR='Index greater than group length'
              GO TO 9999
            ENDIF
          ELSEIF(CBBREV(CO,'LENGTH',3,noco+1,NTCO,N3CO)) THEN
            np=0
          ELSE
            ERROR='No index defined'
            GO TO 9999
          ENDIF
        ENDIF
        
C SETTING VARIABLE
        IF(FROM(1:2).EQ.'XP')THEN       
          REALVALUE=XP(nk,1,nj,np)
          CALL STRING_TRIM(VARNAME,IBEG,IEND)
          CALL SET_USER_DOUBLE(VARNAME(IBEG:IEND),REALVALUE,ERR)
        ELSEIF(FROM(1:5).EQ.'GROUP')THEN
          INTVALUE=NPLIST(np)
          CALL STRING_TRIM(VARNAME,IBEG,IEND)
          CALL SET_USER_INTEGER(VARNAME(IBEG:IEND),INTVALUE,ERR)
        ENDIF
        
        IF(ERR.NE.0) THEN
          ERROR='Unable to set user var "'
     &      //VARNAME(IBEG:IEND)//'"'
          GOTO 9999
        ENDIF

      ENDIF

      CALL EXITS('UPVARI')
      RETURN
 9999 CALL ERRORS('UPVARI',ERROR)
      CALL EXITS('UPVARI')
      RETURN 1
      END


