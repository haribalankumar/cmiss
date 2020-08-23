      SUBROUTINE DEINCR(NHP,NKH,YP,FIX,STRING,ERROR,*)

C#### Subroutine: DEINCR
C###  Description:
C###    DEINCR increments nodal dependent variable.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM)
      REAL*8 YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,ICHAR,IEND,IFROMC,INFO,
     '  N,nc,nh,nhx,nj,nk,NOQUES,np,nr,nx,NYTOT
      CHARACTER STATUS*3
      LOGICAL ABBREV,CALCU,FILEIP,FILIO,GENER,MOUSE

      CALL ENTERS('DEINCR',*9999)
      ICHAR=999

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define increment;p node #
C###  Description:
C###    Increment segments are defined: a thick line joining the
C###    original undeformed position to the deformed position.

        OP_STRING(1)=STRING(1:IEND)//';p node #'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEINCR',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS('P',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)

        nc=1 !Temporary AJP 19-12-91
        nx=1 !Temporary
        nr=1 !Temporary

        FILEIP=.FALSE.
        NOQUES=0

        IF(ABBREV(COQU(noco,1),'P',1)) THEN
          IF(ABBREV(CO(noco+1),'NODE',1)) THEN
            np=IFROMC(CO(noco+2))
 100        FORMAT='($,'' >>Enter coordinate(1..3) & '//
     '        'derivative(0..6) numbers [exit]: '',2I4)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '        IZERO,0,7,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '        INFO,ERROR,*9999)
c           CALL IINOUT(1,IOIP,0,FORMAT,2,IDATA,IZERO,0,7,INFO,
c    '        ERROR,*9999)
            IF(IDATA(1).GT.0) THEN
              nj=IDATA(1)
              nk=IDATA(2)+1
              NYTOT=0
              DO N=1,np-1
                DO nhx=1,NHP(N,nr,nx)
                  nh=NH_LOC(nhx,nx)
                  NYTOT=NYTOT+NKH(nh,N,nc,nr)
                ENDDO
              ENDDO
              DO nhx=1,nj-1
                nh=NH_LOC(nhx,nx)
                NYTOT=NYTOT+NKH(nh,np,nc,nr)
              ENDDO
              NYTOT=NYTOT+nk
              FORMAT='($,'' >>Enter increment [0]: '',E12.4)'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,
     '          INFO,ERROR,*9999)
c             CALL RINOUT(1,IOIP,0,FORMAT,1,RDATA,RZERO,-RMAX,RMAX,INFO,
c    '          ERROR,*9999)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              YP(NYTOT,3,nx)=RDATA(1)
              FIX(NYTOT,1,nx)=.TRUE.
              FIX(NYTOT,3,nx)=.TRUE.
              GO TO 100
            ENDIF
          ELSE
            CO(noco+1)='?'
            GO TO 1
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('DEINCR')
      RETURN
 9999 CALL ERRORS('DEINCR',ERROR)
      CALL EXITS('DEINCR')
      RETURN 1
      END


