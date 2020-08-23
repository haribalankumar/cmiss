      SUBROUTINE IPBOUN(NBJ,NEELEM,NENP,NNB,NPNE,nr,NXI,NYNP,YP,
     '  ALL_REGIONS,FIX,ERROR,*)

C#### Subroutine: IPBOUN
C###  Description:
C###    IPBOUN inputs block boundary conditions.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh01.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),NNB(4,4,4,NBFM),
     '  NPNE(NNM,NBFM,NEM),nr,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL ALL_REGIONS,FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER ICHAR,iface,INFO,n,nb,nc,ne,ni,nn,NNI(4,-3:3),
     '  noelem,NOQUES,np,nx,ny
      CHARACTER CHAR1*1,CHAR2*1
      LOGICAL FILEIP

C          NNI: (n,-3)   (n,-2)   (n,-1)   (n, 0)
C               (n, 1)   (n, 2)   (n, 3)
      DATA NNI/1,2,3,4, 1,2,5,6, 1,3,5,7, 0,0,0,0,
     '         2,4,6,8, 3,4,7,8, 5,6,7,8/

      CALL ENTERS('IPBOUN',*9999)

      nc=1 !Temporary CPB 22/11/94
      nx=1 !temporary

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)

      FORMAT='('' Enter mesh type [1]:'''//
     '  '/''   (1) Regular finite element mesh'''//
     '  '/''   (2) Collocation grid'''//
     '  '/''   (3) Unused'''//
     '  '/''   (4) Unused'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP34
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP34=IDATA(1)

      IF(KTYP34.EQ.1) THEN !Regular finite element mesh
        FORMAT='('' Enter type of specified boundary condition [1]:'''//
     '    '/''   (1) Dirichlet'''//
     '    '/''   (2) Unused'''//
     '    '/''   (3) Unused'''//
     '    '/''   (4) Unused'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP35
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,1,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP35=IDATA(1)

        nb=MESH1_NB
        DO ni=1,NIT(nb)
          DO iface=1,2
            WRITE(CHAR1,'(I1)') iface
            WRITE(CHAR2,'(I1)') ni
            RDEFLT(1)=1.0D6
            FORMAT='($,'' Enter value on face '//CHAR1//' in Xi('//CHAR2
     '        //') dir.n [no b.c.]: '',D12.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=0.0D0
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3.AND.RDATA(1).LT.0.99D6) THEN
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(iface.EQ.1.AND.NXI(-ni,1,ne).EQ.0) THEN
                  IF(DOP) WRITE(*,'('' element '',I5,'' nodes:'',20I5)')
     '              ne,(NPNE(NNI(n,-ni),nb,ne),n=1,NNT(nb)/2)
                  DO nn=1,NNT(nb)/2
                    np=NPNE(NNI(nn,-ni),nb,ne)
C GMH 8/1/97 Update cmgui link
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    ny=NYNP(1,1,NH_LOC(1,nx),np,0,nc,nr)
                    FIX(ny,1)=.TRUE.
                    YP(ny,1)=RDATA(1)
                  ENDDO
                ELSE IF(iface.EQ.2.AND.NXI(ni,1,ne).EQ.0) THEN
                  IF(DOP) WRITE(*,'('' element '',I5,'' nodes:'',20I5)')
     '              ne,(NPNE(NNI(n,ni),nb,ne),n=1,NNT(nb)/2)
                  DO nn=1,NNT(nb)/2
                    np=NPNE(NNI(nn,ni),nb,ne)
C GMH 8/1/97 Update cmgui link
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    ny=NYNP(1,1,NH_LOC(1,nx),np,0,nc,nr)
                    FIX(ny,1)=.TRUE.
                    YP(ny,1)=RDATA(1)
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO

      ELSE IF(KTYP34.EQ.1) THEN !Collocation grid
      ELSE
      ENDIF

      IF(FILEIP.AND..NOT.ALL_REGIONS) CALL CLOSEF(IFILE,ERROR,*9999)
      CALL EXITS('IPBOUN')
      RETURN
 9999 CALL ERRORS('IPBOUN',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPBOUN')
      RETURN 1
      END


C      SUBROUTINE IPCELL(nr,nx,CELL_INIT,CELL_PARAM,ERROR,*)
C
CC#### Subroutine: IPCELL
CC###  Description:
CC###    IPCELL inputs material parameters for cellular models which
CC###    have been created in the CMGUI cellular modelling environment.
CC###    Define materials must also be called for continuum modelling
CC###    problems.
CC**** Created by Martin Buist, July 1998
C
C      IMPLICIT NONE
C
C      INCLUDE 'cmiss$reference:b00.cmn'
C      INCLUDE 'cmiss$reference:b01.cmn'
C      INCLUDE 'cmiss$reference:cell00.cmn'
C      INCLUDE 'cmiss$reference:deoxs00.cmn'
C      INCLUDE 'cmiss$reference:inout00.cmn'
C      INCLUDE 'cmiss$reference:ityp00.cmn'
C      INCLUDE 'cmiss$reference:ktyp30.cmn'
C      INCLUDE 'cmiss$reference:loc00.inc'
C
C!     Parameter List
C      INTEGER nr,nx
C      REAL*8 CELL_INIT(*),CELL_PARAM(*)
C      CHARACTER ERROR*(*)
C!     Local Variables
C      INTEGER i,ICHAR,INFO,NOQUES
C      CHARACTER NUM1*4,TYPE*4
C      LOGICAL FILEIP
C
C      CALL ENTERS('IPCELL',*9999)
C
C      FILEIP=.FALSE.
C      NOQUES=0
C      ICHAR=999
C
C! Read in the type of ionic current model being used
C      DO i=1,8
C        WRITE(NUM1,'(I4)') i+1000
C        FORMAT='(''  (modl'//NUM1(2:4)//')   '',I1)'
C        IDEFLT(1)=0
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) CMMODEL(i)=IDATA(1)
C      ENDDO
C
C      IF(CMMODEL(1).EQ.2) THEN !LR
C        MECHANICS_MODEL=0
C        MEMBRANE_MODEL=CMMODEL(1)
C        TYPE(1:2)='LR'
C        CALL ASSERT(ITYP3(nr,nx).EQ.6,'>>LR model not set',ERROR,*9999)
C      ELSEIF(CMMODEL(1).EQ.3) THEN !JRW
C        MECHANICS_MODEL=0
C        MEMBRANE_MODEL=CMMODEL(1)
C        CALL ASSERT(ITYP3(nr,nx).EQ.5,'>>JRW model not set',ERROR,*9999)
C        IF(KTYP33.EQ.1) THEN !Standard
C          TYPE(1:3)='JRW'
C        ELSE IF(KTYP33.EQ.2) THEN !Princeton
C          TYPE(1:3)='JRWP'
C        ELSE
C          ERROR='>>Invalid JRW model type'
C          GOTO 9999
C        ENDIF
C      ELSEIF(CMMODEL(1).EQ.4) THEN !N98
C        MECHANICS_MODEL=0
C        MEMBRANE_MODEL=CMMODEL(1)
C        TYPE(1:3)='N98'
C        CALL ASSERT(ITYP3(nr,nx).EQ.8,'>>N98 model not set',ERROR,*9999)
C      ELSEIF(CMMODEL(1).EQ.5) THEN !HH
C        MECHANICS_MODEL=0
C        MEMBRANE_MODEL=CMMODEL(1)
C        TYPE(1:2)='HH'
C        CALL ASSERT(ITYP3(nr,nx).EQ.9,'>>HH model not set',ERROR,*9999)
C      ELSE
C        CALL ASSERT(.FALSE.,'>>Ionic model not implemented',ERROR,*9999)
C      ENDIF
C
C! Read in the parameters for the selected ionic current model
C      FORMAT='(''Number of electrical real parameters:  '',I2)'
C      IDEFLT(1)=0
C      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '  LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C      IF(IOTYPE.NE.3) NUM_PARAMS=IDATA(1)
C      CALL ASSERT(NUM_PARAMS.GT.0,'>>No parameters found',ERROR,*9999)
C
C      DO i=1,NUM_PARAMS
C        WRITE(NUM1,'(I4)') i+1000
C        FORMAT='(''  (eleR'//NUM1(2:4)//') '',E12.5)'
C        RDEFLT(1)=0.0d0
C        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) CELL_PARAM(i)=RDATA(1)
C      ENDDO
C
C      IF(TYPE(1:2).EQ.'LR') THEN
Cc        CELL_PARAM(49)=0.0d0
Ccc        Ko=CELL_PARAM(50)
Cc        Kb=CELL_PARAM(50)
Cc        Nao=CELL_PARAM(51)
Cc        Cao=1.8d0
Cc        CALL DEFINE_LR(CELL_PARAM)
C      ELSEIF(TYPE(1:3).EQ.'JRW') THEN
C        ! Set the stimulus size to zero
C        CELL_PARAM(62)=0.0d0
C        CALL DEFINE_JRW(CELL_PARAM)
C      ELSEIF(TYPE(1:3).EQ.'JRWP') THEN
C        ! Do Nothing
C      ELSEIF(TYPE(1:3).EQ.'N98') THEN
C        ! Set the stimulus size to zero
C        CELL_PARAM(1)=0.0d0
C        CALL DEFINE_NOBLE98(CELL_PARAM)
C      ELSEIF(TYPE(1:2).EQ.'HH') THEN
C        ! Do Nothing
C      ELSE
C        ERROR='>>Invalid cell type'
C        GOTO 9999
C      ENDIF
C
C      DO i=1,30 !include all currents in our model
C        ISWTCH(i)=1.0d0
C      ENDDO
C
C! Read in the initial values for the selected ionic current model
C      FORMAT='(''Number of initial values:  '',I2)'
C      IDEFLT(1)=0
C      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '  LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C      IF(IOTYPE.NE.3) NUM_INITS=IDATA(1)
C      CALL ASSERT(NUM_INITS.GT.0,'>>No initial values found',
C     '  ERROR,*9999)
C
C      DO i=1,NUM_INITS
C        WRITE(NUM1,'(I4)') i+1000
C        FORMAT='('' (y'//NUM1(2:4)//') '',E12.5)'
C        RDEFLT(1)=0.0d0
C        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) CELL_INIT(i)=RDATA(1)
C      ENDDO
C
C      CALL EXITS('IPCELL')
C      RETURN
C 9999 CALL ERRORS('IPCELL',ERROR)
C      CALL EXITS('IPCELL')
C      RETURN 1
C      END


