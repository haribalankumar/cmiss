      SUBROUTINE IPREFI(IBT,NBJ,NEELEM,NELIST,NKJE,NPF,NPNE,NRLIST,NVJE,
     '  NEERR,PG,RG,SE,WG,XA,XE,XG,XP,ERROR,*)

C#### Subroutine: IPREFI
C###  Description:
C###    IPREFI inputs refinement parameters
C###  See-Also: IPMESH9_DYNAM
C###  See-Example: 1371,1372
C**** Created by Carey Stevens 1 October 1997

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'mesh00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 NEERR(NEM,3),PG(NSM,NUM,NGM,NBM),
     '  RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,no_nrlist,NOQUES,nr,REFINE_TYPE
      LOGICAL FILEIP

      CALL ENTERS('IPREFI',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      DO no_nrlist=1,NRLIST(0)
        nr=NRLIST(no_nrlist)

        FORMAT='(/'' Enter method for defining refinement [1]:'''//
     '    '/''  (1) By distance field'''//
     '    '/''  (2) '''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) THEN
          REFINE_TYPE=1
          IDATA(1)=REFINE_TYPE
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,1,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) REFINE_TYPE=IDATA(1)

        IF(REFINE_TYPE.EQ.1) THEN ! distance field

          RDEFLT(1)=0.0d0
          FORMAT='($,'' Refine at field values below [0.0]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=REFINE_FIELD_VALUE
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) REFINE_FIELD_VALUE=RDATA(1)

          FORMAT='(/'' Enter field sample point [1]:'''//
     '      '/''  (1) Element center'''//
     '      '/''  (2) Element corners'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=REFINE_FIELD_SAMPLE_POINT
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) REFINE_FIELD_SAMPLE_POINT=IDATA(1)

          FORMAT='($,'' The XI direction to refine is [1]: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NRM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) REFINE_DIR=IDATA(1)

          RDEFLT(1)=0.5d0
          FORMAT='($,'' Refine at XI [0.5]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=REFINE_XI
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) REFINE_XI=RDATA(1)

        ENDIF
      ENDDO !no_nrlist

      CALL EXITS('IPREFI')
      RETURN
 9999 CALL ERRORS('IPREFI',ERROR)
      CALL EXITS('IPREFI')
      RETURN 1
      END


