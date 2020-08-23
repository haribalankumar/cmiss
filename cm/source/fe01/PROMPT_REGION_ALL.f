      SUBROUTINE PROMPT_REGION_ALL(nr,ERROR,*)

C#### Subroutine: PROMPT_REGION_ALL
C###  Description:
C###    PROMPT_REGION_ALL prompts for region.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER nr
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,IOTYPE_tmp,NOQUES
      LOGICAL FILEIP
      CHARACTER CHAR1*4

      CALL ENTERS('PROMPT_REGION_ALL',*9999)

      WRITE(CHAR1,'(I3)') nr
      IDATA(1)=nr

      IF(IOTYPE.EQ.1) THEN !prompted input
        OP_STRING(1)=' '

C LKC 1-MAY-2000
C        CALL WRITES(IVDU,OP_STRING,ERROR,*9999)
C        OP_STRING(1)=' Region number :'//CHAR1(1:4)
C        CALL WRITES(IVDU,OP_STRING,ERROR,*9999)
C        FORMAT='(/$,'' Region number :'',I3)'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        OP_STRING(1)=' Region number :'//CHAR1(1:4)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        FORMAT='(/$,'' Region number :'',I3)'


        IOTYPE_tmp=3
        NOQUES=0
        ICHAR=999
        CALL GINOUT(IOTYPE_tmp,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &    NRT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     &    ERROR,*9999)
      ELSE
        FORMAT='(/$,'' Region number :'',I3)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,1,NRT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
      ENDIF !iotype

      CALL EXITS('PROMPT_REGION_ALL')
      RETURN
 9999 CALL ERRORS('PROMPT_REGION_ALL',ERROR)
      CALL EXITS('PROMPT_REGION_ALL')
      RETURN 1
      END

      

