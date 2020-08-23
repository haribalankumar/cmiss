      SUBROUTINE DRINCR(ISEG,ISINCR,NHP,NKH,NPNODE,NRLIST,NXLIST,
     '  XP,ZP,CSEG,STRING,FIX,ERROR,*)

C#### Subroutine: DRINCR
C###  Description:
C###    DRINCR increments nodal dependent variable.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ISEG(*),ISINCR(NWM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NXLIST(0:NXM)
      REAL*8 XP(NKM,NVM,NJM,NPM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYLINE,iw,IWK(6),
     '  N3CO,noiw,nr,NTIW,nx,nxc
      LOGICAL ABBREV,ALL_REGIONS,CBBREV

      CALL ENTERS('DRINCR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw increments
C###  Description:
C###    Draw increments in the nodal dependent variable.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <region #[1]>
C###    Specify the region number to use.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[black]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRINCR',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        IF(ABBREV(COQU(noco,1),'S',1)) THEN
          IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
            INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH2',CO(N3CO+1))
          ELSE
            INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH2','BLACK')
          ENDIF

          DO noiw=1,NTIW
            IW=IWK(noiw)
C CPB 10/9/92 Changed workstation update mode
            CALL ACWK(iw,1,ERROR,*9999)
            CALL SGINCR(INDEX,ISEG,ISINCR(iw),iw,NHP(1,nr,nx),
     '        NKH(1,1,1,nr),NPNODE,nx,CSEG,FIX(1,1,nx),XP,ZP,
     '        ERROR,*9999)
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('DRINCR')
      RETURN
 9999 CALL ERRORS('DRINCR',ERROR)
      CALL EXITS('DRINCR')
      RETURN 1
      END


