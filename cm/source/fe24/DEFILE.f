      SUBROUTINE DEFILE(NBH,NBJ,NEELEM,NENP,NHE,NHP,NKH,NNB,
     '  NPL,NPNE,NPNODE,NVHP,NW,NXI,NYNE,NYNP,
     '  CE,DL,XP,YP,ZA,ZP,FIX,STRING,ERROR,*)

C#### Subroutine: DEFILE
C###  Description:
C###    DEFILE defines i/p file for another program.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'head00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NNB(4,4,4,NBFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NW(NEM,3,NXM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),DL(3,NLM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER i,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IDCOUNT,IDTOTAL,
     '  IFROMC,IMEMBR,IPCOUNT,IPTOTAL,ISHELL,ix,iy,iz,
     '  n,N3CO,NAC1(8),NAC2(8),nb,nc,ne,NE1,NENEXT,NENUM,NEQUAD(8),
     '  NESTRT,nh,nhx,nj,nl,nn,noelem,nonode,np,nr,nv,nx,ny
      REAL*8 T,U,V,VAL(12),W
      CHARACTER DATE*25,FILE*(MXCH),FILE_EXT*30,FILENAME*50,PACKAGE*7,
     '  STATUS*3,TYPE*20
      LOGICAL CALCU,CBBREV,CONTINUE,FILIO,GENER,LINES,MOUSE,
     '  PTS,SURFACES
      DATA NAC1/1,2,4,3,5,6,8,7/,NAC2/1,2,3,6,9,8,7,4/
C      LOGICAL ABBREV

      CALL ENTERS('DEFILE',*9999)

      nx=1 ! temporary 22/11/94

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define file;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines file parameters for various other programs.
C###    File parameters are read from or
C###    written to the file FILENAME (with extension dependent on the
C###    type of file being written) in the directory specified by PATH.
C###    This writes a file in the appropriate format for any of the
C###    programs listed.  When using shell or membrane elements, the
C###    option is given of either writing 4 noded 2D linear elements,
C###    or writing 8 noded 2D quadratic elements which are formed by
C###    combining 4 ordinary elements.  For the latter case, the number
C###    of elements written is reduced by one quarter.  It is
C###    recommended that the mesh is therefore refined once in each
C###    direction before using this option.  If present, the material
C###    parameters and some initial conditions are also written.  In
C###    most cases, however, some user manipulation will be required
C###    after the file is written.
C###  Parameter:      <abaqus/beasy/lusas/nastran/nisa/phoenix>
C###    Specify the target program for the file.
C###  Parameter:      <shell=(4/8)[4]>
C###    Specify if shell elements are being used and if they are
C###    linear or quadratic.
C###  Parameter:      <membrane=(4/8)[4]>
C###    Specify if membrane elements are being used and if they are
C###    linear or quadratic.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)
     '    //'<abaqus/beasy/lusas/nastran/nisa/phoenix>'
        OP_STRING(3)=BLANK(1:15)//'<shell=(4/8)[4]>'
        OP_STRING(4)=BLANK(1:15)//'<membrane=(4/8)[4]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define file;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]> iges
C###  Description:
C###    Defines file parameters. File parameters are read from or
C###    written to the file FILENAME (with extension dependent on the
C###    type of file being written) in the directory specified by PATH.
C###    When writing an IGES file, any combination of points, lines
C###    and/or surfaces may be written.  "Points" steps through all
C###    points 1..NPT(nr), "lines" through all lines 1..NLT(nr), and
C###    "surfaces" through all elements 1..NET(nr).  "Lines" will
C###    create cubic splines if appropriate (cubic Hermite basis),
C###    otherwise writing straight lines.  "Surfaces" will interpolate
C###    a cubic spline surface from the derivative information for each
C###    element.  After writing, it will be necessary to convert the
C###    carriage return format from Fortran control to carriage return
C###    control, if transferring the file to another machine.
C###  Parameter:      <points/lines/surfaces>
C###    Specify what information is to be included in the output file.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
     '    //' iges'
        OP_STRING(2)=BLANK(1:15)//'<points/lines/surfaces>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEFILE',ERROR,*9999)
      ELSE

        nc=1 ! Temporary MPN 12-Nov-94
        nv=1 ! Temporary MPN 12-Nov-94

        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        CALL ASSERT(NPT(1).GT.0,'>>no nodes defined',ERROR,*9999)
        CALL ASSERT(NET(1).GT.0,'>>no elements defined',ERROR,*9999)

        IF(CBBREV(CO,'ABAQUS',1,noco+1,NTCO,N3CO)) THEN
          PACKAGE ='ABAQUS'
          FILE_EXT='INP'
          IF(CBBREV(CO,'SHELL',1,noco+2,NTCO,N3CO)) THEN
            ISHELL=IFROMC(CO(N3CO+1))
          ELSE
            ISHELL=4
          ENDIF
          IF(CBBREV(CO,'MEMBRANE',1,noco+2,NTCO,N3CO)) THEN
            IMEMBR=IFROMC(CO(N3CO+1))
          ELSE
            IMEMBR=4
          ENDIF
        ELSE IF(CBBREV(CO,'BEASY',1,noco+1,NTCO,N3CO)) THEN
          PACKAGE ='BEASY'
          FILE_EXT='BEASY'
        ELSE IF(CBBREV(CO,'IGES',1,noco+1,NTCO,N3CO)) THEN
          PACKAGE ='IGES'
          FILE_EXT='IGS'
          IF(CBBREV(CO,'POINTS',1,noco+2,NTCO,N3CO)) THEN
            PTS=.TRUE.
          ELSE
            PTS=.FALSE.
          ENDIF
          IF(CBBREV(CO,'LINES',1,noco+2,NTCO,N3CO)) THEN
            LINES=.TRUE.
          ELSE
            LINES=.FALSE.
          ENDIF
          IF(CBBREV(CO,'SURFACES',1,noco+2,NTCO,N3CO)) THEN
            SURFACES=.TRUE.
          ELSE
            SURFACES=.FALSE.
          ENDIF
        ELSE IF(CBBREV(CO,'LUSAS',1,noco+1,NTCO,N3CO)) THEN
          PACKAGE ='LUSAS'
          FILE_EXT='DAT'
        ELSE IF(CBBREV(CO,'NASTRAN',1,noco+1,NTCO,N3CO)) THEN
          PACKAGE ='NASTRAN'
          FILE_EXT='dat'
          IF(CBBREV(CO,'SHELL',1,noco+2,NTCO,N3CO)) THEN
            ISHELL=IFROMC(CO(N3CO+1))
          ELSE
            ISHELL=4
          ENDIF
        ELSE IF(CBBREV(CO,'NISA',1,noco+1,NTCO,N3CO)) THEN
          PACKAGE ='NISA'
          FILE_EXT='NISA'
        ELSE IF(CBBREV(CO,'PHOENIX',1,noco+1,NTCO,N3CO)) THEN
          PACKAGE ='PHOENIX'
          FILE_EXT='PHOENIX'
        ENDIF

        nr=1 !Needs fixing

        IF(FILIO) THEN
          CALL STRING_TRIM(FILE,IBEG1,IEND1)
          CALL STRING_TRIM(FILE_EXT,IBEG2,IEND2)
          FILENAME=FILE(IBEG1:IEND1)//'.'//FILE_EXT(IBEG2:IEND2)
          CALL STRING_TRIM(FILENAME,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILENAME(IBEG:IEND),
     '      STATUS,'SEQUEN','FORMATTED',132,ERROR,*9999)

          IF(ITYP2(nr,nx).EQ.1) THEN
            IF(NW(NEELEM(1,1),1,nx).EQ.1) THEN !truss or cable element
              TYPE='TRUSS'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ. 2) THEN !batten element
              TYPE='BATTEN'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ. 3) THEN !beam (Kirchhoff) element
              TYPE='BEAM'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ. 4) THEN !link element
              TYPE='LINK'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ. 5) THEN !membrane element
              TYPE='MEMBRANE'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ. 6) THEN !plate (Kirchhoff) element
              TYPE='PLATE'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ. 7) THEN !shell element
              TYPE='SHELL'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ. 8) THEN !shell/fluid element
              TYPE='SHELL-FLUID'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ. 9) THEN !liquid surface element
              TYPE='?'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ.10) THEN !tank bottom element
              TYPE='?'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ.11) THEN !plane stress element
              TYPE='PLANE-STRESS'
            ELSE IF(NW(NEELEM(1,1),1,nx).EQ.12) THEN !plane strain element
              TYPE='PLANE-STRAIN'
            ENDIF
          ENDIF

!------------------------------ ABAQUS-----------------------------------

          IF(PACKAGE(1:6).EQ.'ABAQUS') THEN

            WRITE(IFILE,'(''*HEADING'')')
            CALL STRING_TRIM(HEADING,IBEG,IEND)
            WRITE(IFILE,'(A)') HEADING(IBEG:IEND)
            WRITE(IFILE,'(''ABAQUS input file generated by CMISS'')')
            CALL GET_DATE_TIME(DATE,1,ERROR,*9999)
            CALL STRING_TRIM(DATE,IBEG,IEND)
            WRITE(IFILE,'(A)') DATE(IBEG:IEND)
            WRITE(IFILE,'(''*RESTART,WRITE'')')
            IF(ITYP10(1).EQ.1) THEN           !Rectangular Cartesian Coords
              WRITE(IFILE,'(''*NODE,NSET=ALL'')')
              DO nonode=1,NPNODE(0,nr)   !Put all nodes in NodeSet 'ALL'
                np=NPNODE(nonode,nr)
                WRITE(IFILE,'(I5,3E10.3)')
     '            np,(XP(1,nv,nj,np),nj=1,NJ_LOC(NJL_GEOM,0,nr))
              ENDDO
            ELSE IF(ITYP10(1).EQ.2) THEN      !Cylindrical Polar Coords
              WRITE(IFILE,'(''*NODE,NSET=ALL,SYSTEM=C'')')
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                WRITE(IFILE,'(I5,3E10.3)') np,
     '            XP(1,nv,1,np),XP(1,nv,2,np)*180.0D0/PI,XP(1,nv,3,np)
              ENDDO
            ENDIF

            IF(TYPE(1:5).EQ.'SHELL') THEN
              IF(ISHELL.EQ.8) THEN
                TYPE='S8R5'
              ELSE IF(ISHELL.EQ.4) THEN
                TYPE='S4R5'
              ENDIF
            ELSE IF(TYPE(1:8).EQ.'MEMBRANE') THEN
              IF(IMEMBR.EQ.8) THEN
                TYPE='M3D8R'
              ELSE IF(IMEMBR.EQ.4) THEN
                TYPE='M3D4R'
              ENDIF
            ENDIF
            CALL STRING_TRIM(TYPE,IBEG,IEND)
            WRITE(IFILE,'(''*ELEMENT,TYPE='//TYPE(1:IEND)
     '        //'ELSET=PRINT'')')

            IF(TYPE(1:4).EQ.'S8R5'.OR.TYPE(1:5).EQ.'M3D8R') THEN
                                   ! 4 4-noded cubics -> 1 8-noded quadratic
              CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999) ! Init vars
              ne=NEELEM(1,nr)
              NE1=ne
              NESTRT=0
              DO WHILE(NXI(-1,1,ne).gt.0.and.NE.NE.NESTRT)
                NESTRT=NE1             !Find elem in lower left corner
                ne=NXI(-1,1,ne)
              ENDDO
              DO WHILE(NXI(-2,1,ne).GT.0)
                ne=NXI(-2,1,ne)
              ENDDO
              NESTRT=ne
              NENEXT=NXI(2,1,ne)
              IF(NENEXT.GT.0) THEN
                NENEXT=NXI(2,1,NENEXT)   ! Find next row
              ELSE
                ERROR='>>Must have even # elements in Xi(2) dir.n'
                GO TO 9999
              ENDIF
              NENUM=1             !initializes count of ABAQUS quad elements
              CONTINUE=.TRUE.
              DO WHILE(CONTINUE)
                NEQUAD(1)=NPNE(1,1,ne)
                NEQUAD(5)=NPNE(2,1,ne)
                NEQUAD(8)=NPNE(3,1,ne)
                ne=NXI(2,1,ne)                ! go to elem above
                NEQUAD(4)=NPNE(3,1,ne)
                NEQUAD(7)=NPNE(4,1,ne)
                ne=NXI(1,1,ne)                ! go to the right
                IF(ne.GT.0) THEN
                  NEQUAD(6)=NPNE(2,1,ne)
                  NEQUAD(3)=NPNE(4,1,ne)
                ELSE
                  ERROR='>>Must have even # elements in Xi(1) dir.n'
                  GO TO 9999
                ENDIF
                ne=NXI(-2,1,ne)               ! and down
                NEQUAD(2)=NPNE(2,1,ne)
                WRITE(IFILE,'(16I5)') NENUM,(NEQUAD(nn),nn=1,8)
                NENUM=NENUM+1
                ne=NXI(1,1,ne)            ! to next starting block of 4
                IF(ne.eq.nestrt.or.NE.LE.0) THEN
                  ne=NENEXT                 ! if completed circuit, move
                  IF(ne.EQ.0) THEN
                    CONTINUE=.FALSE.
                  ELSE
                    NENEXT=NXI(2,1,ne)          ! up to begin another block
                    IF(NENEXT.GT.0) THEN
                      NENEXT=NXI(2,1,NENEXT)    ! Find next row
                    ELSE
                      ERROR='>>Must have even # elements in Xi(2) dir.n'
                      GO TO 9999
                    ENDIF
                    NESTRT=ne
                  ENDIF
                ENDIF
              ENDDO

            ELSE IF(TYPE(1:4).EQ.'S4R5'.OR.TYPE(1:5).EQ.'M3D4R') THEN
              ! 4-Noded Quadrilateral
              DO noelem=1,NEELEM(0,1)
                ne=NEELEM(noelem,1)
                nb=NBJ(1,ne)
                WRITE(IFILE,'(16I5)') ne,(NPNE(NAC1(nn),nb,ne),nn=1,
     '            NNT(nb))
              ENDDO
            ENDIF
            IF(TYPE(1:4).EQ.'S8R5'.OR.TYPE(1:4).EQ.'S4R5') THEN
              WRITE(IFILE,'(''*SHELL SECTION'')')
              WRITE(IFILE,'(E8.3)') CE(3,NEELEM(1,1),nx)
            ELSE IF(TYPE(1:5).EQ.'M3D4R'.OR.TYPE(1:5).EQ.'M3D8R') THEN
              WRITE(IFILE,'(''*SOLID SECTION'')')
              WRITE(IFILE,'(E8.3)') CE(3,NEELEM(1,1),nx)
            ENDIF
            WRITE(IFILE,'(''*MATERIAL'')')
            WRITE(IFILE,'(''*ELASTIC'')')
            WRITE(IFILE,'(2E8.3)') (CE(i,NEELEM(1,1),nx),i=1,2)
            IF(CE(7,NEELEM(1,1),nx).GE.1D-6) THEN
              WRITE(IFILE,'(''*DENSITY'')')
              WRITE(IFILE,'(E8.3)') CE(7,NEELEM(1,1),nx)
            ENDIF
            WRITE(IFILE,'(''*BOUNDARY'')')
            WRITE(IFILE,*)
            WRITE(IFILE,'(''*STEP,LINEAR'')')
            IF(ITYP5(nr,nx).EQ.1) THEN
              WRITE(IFILE,'(''*STATIC'')')
            ELSE IF(ITYP5(nr,nx).EQ.2) THEN
              WRITE(IFILE,'(''*DYNAMIC'')')
            ENDIF
            WRITE(IFILE,'(''*CLOAD'')')
            WRITE(IFILE,*)
            WRITE(IFILE,'(''*EL PRINT, ELSET=PRINT'')')
            WRITE(IFILE,'(''S,E,MISES'')')
            WRITE(IFILE,'(''*END STEP'')')


!------------------------------- BEASY ------------------------------------

          ELSE IF(PACKAGE(1:5).EQ.'BEASY') THEN


!-------------------------------- IGES -------------------------------------
! Note: IGES number: 116 is 2D or 3D point
!                    114 is surface ?

          ELSE IF(PACKAGE(1:5).EQ.'IGES') THEN

            CALL ASSERT(ITYP10(1).EQ.1,
     '        '>>IGES requires nodes defined in r.c. coords',
     '        ERROR,*9999)

       !****Start Section
            FORMAT='(A,T73,A,I7)'
            WRITE(IFILE,FORMAT) 'IGES Data file generated by CMISS',
     '        'S',1
            CALL STRING_TRIM(HEADING,IBEG,IEND)
            WRITE(IFILE,FORMAT) HEADING(IBEG:IEND),'S',2

       !****Global Section
            CALL STRING_TRIM(FILENAME,IBEG,IEND)
            WRITE(IFILE,'('',,5HCMISS,'',I2,''H'',A,'','',T73,A,I7)')
     '        (IEND-IBEG+1),FILENAME(IBEG:IEND),'G',1
            WRITE(IFILE,'(''5HCMISS,1H1,16,8,24,8,56,'',T73,A,I7)')
     '        'G',2
            WRITE(IFILE,'(''5HCMISS,1.0,2,2HMM,1,0.01,'',T73,A,I7)')
     '        'G',3
            CALL GET_DATE_TIME(DATE,2,ERROR,*9999)
            CALL STRING_TRIM(DATE,IBEG,IEND)
            WRITE(IFILE,'(I2,''H'',A,'',0.0001,10000.0,'',T73,A,I7)')
     '        (IEND-IBEG+1),DATE(IBEG:IEND),'G',4
            WRITE(IFILE,
     '        '(''19HENGINEERING SCIENCE,19HAUCKLAND UNIVERSITY;'','
     '        //'T73,A,I7)') 'G',5

       !****Directory Entry Section
            IDCOUNT=0
            IPCOUNT=0

            !--Points
            IF(PTS) THEN
              DO np=1,NPT(1)
                IDCOUNT=IDCOUNT+1
                IPCOUNT=IPCOUNT+1
                WRITE(IFILE,'(5I8,24X,A8,A1,I7)')
     '                116,IPCOUNT,1,1,1,'00000000','D',IDCOUNT
                IDCOUNT=IDCOUNT+1
                WRITE(IFILE,'(4I8,24X,A8,I8,A1,I7)')
     '                116,0,1,1,'POINT',np,'D',IDCOUNT
              ENDDO
            ENDIF !points

            !--Lines
            IF(LINES) THEN
              DO ix=0,1  !Write points to give axes
                DO iy=0,1
                  DO iz=0,1
                    IDCOUNT=IDCOUNT+1
                    IPCOUNT=IPCOUNT+1
                    WRITE(IFILE,'(5I8,24X,A8,A1,I7)')
     '                    116,IPCOUNT,1,1,1,'00000000','D',IDCOUNT
                    IDCOUNT=IDCOUNT+1
                    WRITE(IFILE,'(4I8,24X,A8,I8,A1,I7)')
     '                    116,0,1,1,'POINT',np,'D',IDCOUNT
                  ENDDO
                ENDDO
              ENDDO
              IF(NPL(1,1,1).LE.3) THEN  !Not cubic Hermite - use line segments
                DO nl=1,NLT
                  IDCOUNT=IDCOUNT+1
                  IPCOUNT=IPCOUNT+1
                  WRITE(IFILE,'(5I8,24X,A8,A1,I7)')
     '                  110,IPCOUNT,1,1,1,'00000000','D',IDCOUNT
                  IDCOUNT=IDCOUNT+1
                  WRITE(IFILE,'(4I8,24X,A8,I8,A1,I7)')
     '                  110,0,1,1,'LINE',nl,'D',IDCOUNT
                ENDDO
              ELSE  !Cubic Hermite - use cubic splines
                DO nl=1,NLT
                  IDCOUNT=IDCOUNT+1
                  IPCOUNT=IPCOUNT+1
                  WRITE(IFILE,'(5I8,24X,A8,A1,I7)')
     '                  112,IPCOUNT,1,1,1,'00000000','D',IDCOUNT
                  IDCOUNT=IDCOUNT+1
                  WRITE(IFILE,'(4I8,24X,A8,I8,A1,I7)')
     '                  112,0,1,4,'SPLINE',nl,'D',IDCOUNT
                  IPCOUNT=IPCOUNT+3
                ENDDO

              ENDIF
            ENDIF !lines

            !Surfaces
            IF(SURFACES) THEN
              DO noelem=1,NEELEM(0,1)
                ne=NEELEM(noelem,1)
                IDCOUNT=IDCOUNT+1
                IPCOUNT=IPCOUNT+1
                WRITE(IFILE,'(5I8,24X,A8,A1,I7)')
     '            114,IPCOUNT,1,1,1,'00000000','D',IDCOUNT
                IDCOUNT=IDCOUNT+1
                WRITE(IFILE,'(4I8,24X,A8,I8,A1,I7)')
     '            114,0,1,1,'SURFACE',ne,'D',IDCOUNT
              ENDDO
            ENDIF !surfaces

            IDTOTAL=IDCOUNT

       !****Parameter Data Section
            IPCOUNT=0
            IDCOUNT=-1
            !--Points
            IF(PTS) THEN
              DO np=1,NPT(1)
                IPCOUNT=IPCOUNT+1
                IDCOUNT=IDCOUNT+2
                WRITE(IFILE,'(I8,3('','',F8.4),'';'',T65,I8,A1,I7)')
     '            116,(XP(1,nv,nj,np),nj=1,3),IDCOUNT,'P',IPCOUNT
              ENDDO
            ENDIF !points

            !--Lines
            IF(LINES) THEN
              DO ix=0,1 !Write points to give axes
                DO iy=0,1
                  DO iz=0,1
                    IPCOUNT=IPCOUNT+1
                    IDCOUNT=IDCOUNT+2
                    WRITE(IFILE,'(I8,3('','',F10.4),'';'','
     '                //'T65,I8,A1,I7)')116,REAL(ix*100),REAL(iy*100),
     '                REAL(iz*100),IDCOUNT,'P',IPCOUNT
                  ENDDO
                ENDDO
              ENDDO
              IF(NPL(1,1,1).LE.3) THEN  !Not cubic Hermite - use line segments
                DO nl=1,NLT
                  IDCOUNT=IDCOUNT+2
                  IPCOUNT=IPCOUNT+1
                  WRITE(IFILE,'(I8,6('','',F10.4),'';'',T65,I8,A1,I7)')
     '              110,((XP(1,nv,nj,NPL(i,1,nl)),nj=1,3),i=2,3),
     '              IDCOUNT,'P',IPCOUNT
                ENDDO
              ELSE  !Cubic Hermite - use cubic splines
                DO nl=1,NLT
                  IDCOUNT=IDCOUNT+2
                  IPCOUNT=IPCOUNT+1
                  WRITE(IFILE,
     '            '(I8,4('','',I3),2('','',F10.4),'','',T65,I8,A1,I7)')
     '              112,3,0,3,1,0.0,1.0,IDCOUNT,'P',IPCOUNT
                  DO nj=1,3
                    IPCOUNT=IPCOUNT+1
                    T=XP(1,nv,nj,NPL(2,1,nl))      !Form cubic spline
                    U=XP(NPL(4,1,nl),nv,nj,NPL(2,1,nl))*DL(1,nl)
                    V=XP(1,nv,nj,NPL(3,1,nl))
                    W=XP(NPL(5,1,nl),nv,nj,NPL(3,1,nl))*DL(2,nl)
                    VAL(1)=T
                    VAL(2)=U
                    VAL(3)=3*(V-T)-2*U-W
                    VAL(4)=2*(T-V)+U+W
                    IF(nj.EQ.3) THEN
                      WRITE(IFILE,'(3(F14.7,'',''),F14.7,'';'',T65,I8,'
     '                  //'A1,I7)')(VAL(n),n=1,4),IDCOUNT,'P',IPCOUNT
                    ELSE
                      WRITE(IFILE,'(4(F14.7,'',''),T65,I8,A1,I7)')
     '                  (VAL(n),n=1,4),IDCOUNT,'P',IPCOUNT
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF
            ENDIF !lines

            !Surfaces
            IF(SURFACES) THEN
              DO noelem=1,NEELEM(0,1)
                ne=NEELEM(noelem,1)
                IDCOUNT=IDCOUNT+2
                IPCOUNT=IPCOUNT+1
              ENDDO
            ENDIF !surfaces

            IPTOTAL=IPCOUNT

       !****Terminate Section
            WRITE(IFILE,
     '        '(''S      2G      5D'',I7,''P'',I7,T73,''T      1'')')
     '        IDTOTAL,IPTOTAL


!------------------------------- LUCAS ------------------------------------

          ELSE IF(PACKAGE(1:5).EQ.'LUSAS') THEN
            IF(TYPE(1:5).EQ.'TRUSS') THEN
              TYPE='B???'
            ELSE IF(TYPE(1:5).EQ.'PLATE') THEN
              TYPE='QF4'
            ELSE IF(TYPE(1:5).EQ.'SHELL') THEN
              TYPE='QSL8'
            ENDIF
            CALL STRING_TRIM(TYPE,IBEG,IEND)
            WRITE(IFILE,'(1X,A,''    ELEMENT TOPOLOGY'')') TYPE(1:IEND)
            DO noelem=1,NEELEM(0,1)
              ne=NEELEM(noelem,1)
              nb=NBJ(1,ne)
              IF(TYPE(1:3).EQ.'QF4') THEN
                WRITE(IFILE,'(16I7)') ne,(NPNE(NAC1(nn),nb,ne),nn=1,
     '            NNT(nb))
              ELSE IF(TYPE(1:4).EQ.'QSL8') THEN
                WRITE(IFILE,'(16I7)') ne,(NPNE(NAC2(nn),nb,ne),nn=1,8)
              ENDIF
            ENDDO
            WRITE(IFILE,'('' SOLUTION ORDER AUTOMATIC'')')
            WRITE(IFILE,'('' NODE COORDINATES'')')
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              WRITE(IFILE,'(I9,3E16.4)')
     '          np,(XP(1,nv,nj,np),nj=1,NJ_LOC(NJL_GEOM,0,nr))
            ENDDO
            CALL STRING_TRIM(TYPE,IBEG,IEND)
            WRITE(IFILE,
     '        '(1X,A,''    GEOMETRIC PROPERTIES'')') TYPE(1:IEND)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              IF(TYPE(1:3).EQ.'QF4') THEN
                WRITE(IFILE,'(I7,'' 0 0'',E10.3,3(3X,E10.3))')
     '          ne,(CE(3,ne,nx),nn=1,NNT(nb))
              ELSE IF(TYPE(1:4).EQ.'QSL8') THEN
                WRITE(IFILE,'(I7,'' 0 0'',8E13.3)')
     '          ne,(CE(3,ne,nx),nn=1,8)
              ENDIF
            ENDDO
            WRITE(IFILE,'('' MATERIAL PROPERTIES'')')
            DO noelem=1,NEELEM(0,1)
              ne=NEELEM(noelem,1)
              nb=NBJ(1,ne)
              WRITE(IFILE,'(I7,'' 0 0'',5E14.4)')
     '          ne,(CE(i,ne,nx),i=1,2),0.0,0.0,0.0
            ENDDO
            WRITE(IFILE,'('' SUPPORT NODES'')')
            ny=0
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              DO nhx=1,NHP(np,nr,nx)
                nh=NH_LOC(nhx,nx)
                ny=ny+NKH(nh,np,nc,nr)
              ENDDO
              IF(ny.GT.0) THEN
                IF(TYPE(1:3).EQ.'QF4') THEN
                  IF(FIX(ny,1,nx)) WRITE(IFILE,'(I7,'' 0 0 R F F'')')
     '              np
                ELSE IF(TYPE(1:4).EQ.'QSL8') THEN
                  IF(FIX(ny,1,nx))
     '              WRITE(IFILE,'(I7,'' 0 0 R R R F F'')') np
                ENDIF
              ENDIF
            ENDDO
            WRITE(IFILE,'('' LOAD CASE'')')
            WRITE(IFILE,'('' CL'')')
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(FIX(1+(np-1)*NHP(np,nr,nx),2,nx)) THEN
                WRITE(IFILE,'(I7,'' 0 0 '',3E13.4)')
     '            np,(YP(nh+(np-1)*NHP(np,nr,nx),2,nx),
     '            nh=1,NHP(np,nr,nx))
C AJP 22/1/96 Should use nh_loc
              ENDIF
            ENDDO
            WRITE(IFILE,'('' PLOT FILE'')')
            IF(ITYP5(nr,nx).EQ.1) THEN
            ELSE IF(ITYP5(nr,nx).EQ.2) THEN
            ENDIF
            WRITE(IFILE,'('' END'')')


!------------------------------- NASTRAN ------------------------------------

          ELSE IF(PACKAGE(1:7).EQ.'NASTRAN') THEN

            WRITE(IFILE,'(''ID CMISS,NASTRAN'')')
            WRITE(IFILE,'(''SOL 101    $ V66 - Static Analysis'')') !Default soln proc
            WRITE(IFILE,'(''TIME 15'')')
            WRITE(IFILE,'(''CEND'')')
            WRITE(IFILE,
     '        '(''TITLE = MSC/NASTRAN file generated by CMISS'')')
            CALL STRING_TRIM(HEADING,IBEG,IEND)
            WRITE(IFILE,'(''SUBTITLE = '',A)') HEADING(IBEG:IEND)
            WRITE(IFILE,'(''DISP = ALL'')')            ! Following outputs may
            WRITE(IFILE,'(''STRESS = ALL'')')          ! be modified to suit
            WRITE(IFILE,'(''STRAIN = ALL'')')          ! the user
            WRITE(IFILE,'(''FORCE = ALL'')')           !
            WRITE(IFILE,'(''ESE = ALL'')')             !
            WRITE(IFILE,'(''SUBCASE 1'')')
            WRITE(IFILE,'(''   SPC = 1'')')            ! Default subcase
            WRITE(IFILE,'(''   LOAD = 1'')')
            WRITE(IFILE,'(''BEGIN BULK'')')            ! Start of data definition

            IF(ITYP10(1).EQ.2) THEN            !Cylindrical Polar Coords
              WRITE(IFILE,'(''$'')')
              WRITE(IFILE,
     '          '(''$ Cylindrical polar coordinate definition'')')
              WRITE(IFILE,'(''$'')')
              WRITE(IFILE,
     '         '(''CORD2C  '',2I8,6F8.2,''+'')') 1, 0, (0.0,nc=1,5), 1.0
              WRITE(IFILE,'(''        '',3F8.2)') 1.0, 0.0, 0.0
              WRITE(IFILE,'(''$'')')
              WRITE(IFILE,'(''$ Grid point (cyl.pol.) definition'')')
              WRITE(IFILE,'(''$'')')
              DO nonode=1,NPNODE(0,1)
                np=NPNODE(nonode,1)
                WRITE(IFILE,'(''GRID    '',2I8,3E8.3E1)') np,1,
     '            XP(1,nv,1,np),XP(1,nv,2,np)*180.0d0/PI,XP(1,nv,3,np)
              ENDDO
              !Grid Format: Node #, Coord sys, 3 Coords (angles in degrees)
            ELSE IF(ITYP10(1).EQ.1) THEN       !Rectangular Cartesian Coords
              WRITE(IFILE,'(''$'')')
              WRITE(IFILE,'(''$ Grid point (r.c.) definition'')')
              WRITE(IFILE,'(''$'')')
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                WRITE(IFILE,'(''GRID    '',2I8,3E8.3E1)') np,0,
     '            (XP(1,nv,nj,np),nj=1,NJ_LOC(NJL_GEOM,0,nr))
              ENDDO
            ENDIF

            IF(TYPE(1:5).EQ.'SHELL') THEN
              IF(ISHELL.EQ.8) THEN
                TYPE='CQUAD8'
              ELSE IF(ISHELL.EQ.4) THEN
                TYPE='CQUAD4'
              ENDIF
            ELSE
              TYPE='CQUAD4'
C              TYPE='UNKNOWN'
            ENDIF
            WRITE(IFILE,'(''$'')')
            WRITE(IFILE,'(''$ Element definition  --  '//TYPE//''')')
            WRITE(IFILE,'(''$'')')
            IF(TYPE(1:6).EQ.'CQUAD4') THEN         !4-noded shell elements
              DO noelem=1,NEELEM(0,1)
                ne=NEELEM(noelem,1)
                nb=NBJ(1,ne)
                WRITE(IFILE,'(''CQUAD4  '',6I8)')
     '            ne,1,(NPNE(NAC1(nn),nb,ne),nn=1,NNT(nb))
              ENDDO
            ELSE IF(TYPE(1:6).EQ.'CQUAD8') THEN
                                   ! 4 4-noded cubics -> 1 8-noded quadratic
              CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999) ! Init vars
              ne=NEELEM(1,1)
              NE1=ne
              NESTRT=0
              DO WHILE(NXI(-1,1,ne).gt.0.and.NE.NE.NESTRT)
                NESTRT=NE1             !Find elem in lower left corner
                ne=NXI(-1,1,ne)
              ENDDO
              DO WHILE(NXI(-2,1,ne).GT.0)
                ne=NXI(-2,1,ne)
              ENDDO
              NESTRT=ne
              NENEXT=NXI(2,1,ne)
              IF(NENEXT.GT.0) THEN
                NENEXT=NXI(2,1,NENEXT)   ! Find next row
              ELSE
                ERROR='>>Must have even # elements in Xi(2) dir.n'
                GO TO 9999
              ENDIF
              NENUM=1             !initializes count of NASTRAN quad elements
              CONTINUE=.TRUE.
              DO WHILE(CONTINUE)
                NEQUAD(1)=NPNE(1,1,ne)
                NEQUAD(5)=NPNE(2,1,ne)
                NEQUAD(8)=NPNE(3,1,ne)
                ne=NXI(2,1,ne)                ! go to elem above
                NEQUAD(4)=NPNE(3,1,ne)
                NEQUAD(7)=NPNE(4,1,ne)
                ne=NXI(1,1,ne)                ! go to the right
                IF(ne.GT.0) THEN
                  NEQUAD(6)=NPNE(2,1,ne)
                  NEQUAD(3)=NPNE(4,1,ne)
                ELSE
                  ERROR='>>Must have even # elements in Xi(1) dir.n'
                  GO TO 9999
                ENDIF
                ne=NXI(-2,1,ne)               ! and down
                NEQUAD(2)=NPNE(2,1,ne)
                WRITE(IFILE,'(''CQUAD8  '',8I8)')
     '            NENUM,1,(NEQUAD(nn),nn=1,6)
                WRITE(IFILE,'(8X,2I8)') (NEQUAD(nn),nn=7,8)
                NENUM=NENUM+1
                ne=NXI(1,1,ne)            ! to next starting block of 4
                IF(ne.eq.nestrt.or.NE.LE.0) THEN
                  ne=NENEXT                 ! if completed circuit, move
                  IF(ne.EQ.0) THEN
                    CONTINUE=.FALSE.
                  ELSE
                    NENEXT=NXI(2,1,ne)       ! up to begin another block
                    IF(NENEXT.GT.0) THEN
                      NENEXT=NXI(2,1,NENEXT)    ! Find next row
                    ELSE
                      ERROR='>>Must have even # elements in Xi(2) dir.n'
                      GO TO 9999
                    ENDIF
                    NESTRT=ne
                  ENDIF
                ENDIF
              ENDDO
            ELSE IF(TYPE(1:7).EQ.'UNKNOWN') THEN         !not defined
              DO noelem=1,NEELEM(0,1)
                ne=NEELEM(noelem,1)
                nb=NBJ(1,ne)
                WRITE(IFILE,'(''unknown '',6I8)')
     '            ne,1,(NPNE(NAC1(nn),nb,ne),nn=1,NNT(nb))
              ENDDO
            ENDIF

            WRITE(IFILE,'(''$'')')
            WRITE(IFILE,'(''$ Loads definition'')')
            WRITE(IFILE,'(''$'')')
            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '        NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '        YP(1,1,nx),ZA,ZP,ERROR,*9999)
            ny=0
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              DO nhx=1,NHP(np,nr,nx)
                nh=NH_LOC(nhx,nx)
                ny=ny+NKH(nh,np,nc,nr)
              ENDDO
              IF(FIX(ny,2,nx)) WRITE(IFILE,'(''FORCE   '',3I8,4E8.3)')
     '          1,np,0,1.0,
     '          (ZP(1,nv,NH_LOC(nhx,nx),np,nc),nhx=1,NHP(np,nr,nx))
            ENDDO

            WRITE(IFILE,'(''$'')')
            WRITE(IFILE,'(''$ Single Point Constraint definition'')')
            WRITE(IFILE,'(''$'')')
            ny=0
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              DO nhx=1,NHP(np,nr,nx)
                nh=NH_LOC(nhx,nx)
                ny=ny+NKH(nh,np,nc,nr)
                IF(FIX(ny,1,nx)) WRITE(IFILE,'(''SPC     '',3I8,E8.3)')
     '            1,np,nh,ZP(1,nv,nh,np,nc)
              ENDDO
            ENDDO

            WRITE(IFILE,'(''$'')')
            WRITE(IFILE,'(''$ Property definition'')')
            WRITE(IFILE,'(''$'')')
            IF(TYPE(1:5).EQ.'CQUAD') THEN
              WRITE(IFILE,'(''PSHELL  '',2I8,E8.3,I8)')
     '          1,1,CE(3,NEELEM(1,1),nx),1
            ENDIF
            WRITE(IFILE,'(''$'')')
            WRITE(IFILE,'(''$ Material definition'')')
            WRITE(IFILE,'(''$'')')
            WRITE(IFILE,'(''MAT1    '',I8,E8.3,8X,E8.3)')
     '        1,(CE(i,NEELEM(1,1),nx),i=1,2)

            WRITE(IFILE,'(''$'')')
            WRITE(IFILE,'(''$ Parameters'')')
            WRITE(IFILE,'(''$'')')
            WRITE(IFILE,'(''PARAM   POST    0'')')
            WRITE(IFILE,'(''PARAM   AUTOSPC YES'')')

            WRITE(IFILE,'(''ENDDATA'')')

!--------------------------------- NISA ------------------------------------

          ELSE IF(PACKAGE(1:4).EQ.'NISA') THEN
            WRITE(IFILE,'(''** EXECUTIVE'')')
            WRITE(IFILE,'(''*TITLE'')')
            WRITE(IFILE,'(''NISA file generated by CMISS'')')

            WRITE(IFILE,'(''*ELTYPE'')')

            WRITE(IFILE,'(''*RCTABLE'')')

            WRITE(IFILE,'(''*NODES'')')
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(ITYP10(1).EQ.1) THEN      !rect.cart. coords
                WRITE(IFILE,'(I6,4('',''),1X,3(E12.5,'',''),5X,''0'')')
     '            np,(XP(1,nv,nj,np),nj=1,NJ_LOC(NJL_GEOM,0,nr))
              ELSE IF(ITYP10(1).EQ.2) THEN !cylind. polar coords
                WRITE(IFILE,'(I6,4('',''),1X,3(E12.5,'',''),5X,''0'')')
     '            np,XP(1,nv,1,np),XP(1,nv,2,np)*180.0d0/PI,
     '            XP(1,nv,3,np)
              ENDIF
            ENDDO

            WRITE(IFILE,'(''*ELEMENT'')')
            DO noelem=1,NEELEM(0,1)
              ne=NEELEM(noelem,1)
              nb=NBJ(1,ne)
              WRITE(IFILE,'(1X,I6,'','',5X,''M,'',5X,''E,'',5X,''R,'','
     '          //'5X,''0'')')ne
              WRITE(IFILE,'(1X,9(I6,'',''))') (NPNE(NAC1(nn),nb,ne),
     '          nn=1,NNT(nb))
            ENDDO

            WRITE(IFILE,'(''*MATERIAL'')')

            WRITE(IFILE,'(''*LDCASE'')')

            WRITE(IFILE,'(''*SPDISP'')')

            WRITE(IFILE,'(''*PRESSURE'')')

            WRITE(IFILE,'(''*ENDDATA'')')


!-------------------------------- PHOENIX ------------------------------------

          ELSE IF(PACKAGE(1:6).EQ.'PHOENIX') THEN

          ENDIF !package

!---------------------------------------- ------------------------------------

          CALL CLOSEF(IFILE,ERROR,*9999)

        ENDIF !filio

      ENDIF

      CALL EXITS('DEFILE')
      RETURN
 9999 CALL ERRORS('DEFILE',ERROR)
      CALL EXITS('DEFILE')
      RETURN 1
      END


