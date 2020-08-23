      SUBROUTINE UPFIBR(CP,NKH,NKJ,NPLIST,NPNODE,NRLIST,NVHP,
     '  NVJP,NXLIST,NYNP,XP,YP,ZD,STRING,ERROR,*)

C#### Subroutine: UPFIBR
C###  Description:
C###    Updates fibre variables from geometry, fibre, field, solution
C###    or material variables.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NRLIST(0:NRM),NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),
     '  NXLIST(0:NXM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CP(NMM,NPM,NXM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     &  ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!      LOGICAL
!     Local Variables
      INTEGER IBEG,IEND,N3CO,PART2
!      INTEGER*4
!      REAL*8 
      LOGICAL ABBREV,ALL_REGIONS,CBBREV
      CHARACTER OPERATION*16,TYPE*9,UPDATE*16

      CALL ENTERS('UPFIBR',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

        UPDATE='HELPFIBRE'
        CALL UPFG(UPDATE,%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     '    %VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     '    %VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     &    STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPFIBR',ERROR,*9999)
      ELSE

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)

        IF(CBBREV(CO,'SUBSTITUTE',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPERATION'
          OPERATION='SUBSTITUTE'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'ADD',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPERATION'
          OPERATION='ADD'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'SUBTRACT',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPERATION'
          OPERATION='SUBTRACT'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'MULTIPLY',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPERATION'
          OPERATION='MULTIPLY'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'DIVIDE',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPERATION'
          OPERATION='DIVIDE'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'FROM',2,noco+1,NTCO,n3co)) THEN
          IF(ABBREV(CO(n3co+1),'GEOMETRY',1)) THEN
            TYPE='GEOMETRY'
          ELSE IF(ABBREV(CO(n3co+1),'SOLUTION',1)) THEN
            TYPE='SOLUTION'
          ELSE IF(ABBREV(CO(n3co+1),'MATERIAL',1)) THEN
            TYPE='MATERIAL'
          ELSE
            CO(noco+1)='?'
            GO TO 1
          ENDIF
        ELSE
          CO(noco+1)='?'
          GO TO 1
        ENDIF


        IF(TYPE(1:9).EQ.'OPERATION') THEN

          UPDATE='FIBRE'
          CALL UPFG(UPDATE,CP,%VAL(0),%VAL(0),NKH,NKJ,NPLIST,NPNODE,
     '      NRLIST,NVHP,NVJP,NYNP,OPERATION,PART2,%VAL(0),XP,%VAL(0),YP,
     '      ZD,%VAL(0),ERROR,*9999)
     
        ELSE
          
        ENDIF
      ENDIF

      CALL EXITS('UPFIBR')
      RETURN
 9999 CALL ERRORS('UPFIBR',ERROR)
      CALL EXITS('UPFIBR')
      RETURN 1
      END


