      SUBROUTINE EXVORO(NBJ,NEELEM,NXI,VOROLINE,VOROLINE_SIZE,ZA,
     '  STRING,ERROR,*)

C#### Subroutine: EXVORO
C###  Description:
C###    Exports the lines of the Voronoi mesh for region nr

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'parameters.inc'
!     Passed Parameters
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),VOROLINE_SIZE,
     '  VOROLINE(2,VOROLINE_SIZE)
      REAL*8 ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local variables
      INTEGER i,IEND,IBEG,IEND1,IBEG1,nbor,ne,nl,IFROMC,
     '  noelem,nonbor,N3CO,nb,VLT,OUT_FREQ,nr
      CHARACTER FILE*100,XYZ_LISTING(3)*38,XYZ_BASIS(3)*5
      LOGICAL CBBREV
      DATA XYZ_LISTING(1)/'   x.  Value index= 1, #Derivatives= 0'/,
     '     XYZ_LISTING(2)/'   y.  Value index= 2, #Derivatives= 0'/,
     '     XYZ_LISTING(3)/'   z.  Value index= 3, #Derivatives= 0'/
      DATA XYZ_BASIS/'   x.','   y.','   z.'/

      CALL ENTERS('EXVORO',*9999)
C#### Command: FEM export voronoi<;FILENAME[default]>
C###  Parameter:      <region #[1]>
C###    Specify the region number.
C###  Parameter:      <output FREQUENCY>
C###    Report progress every FREQUENCY elements.
C###  Description:
C###    Exports the lines of the Voronoi mesh.

C      CALL ASSERT(CALL_RECONNECT,
C     '  'Need to update Delaunay mesh before export',
C     '  ERROR,*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']>'
        OP_STRING(2)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        OP_STRING(3)=BLANK(1:15)//'<output (frequency)>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
      ELSEIF(CO(noco+1).EQ.'??') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']>'
        OP_STRING(2)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        OP_STRING(3)=BLANK(1:15)//'<output (frequency)>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
      ELSE
C     Calculate the elements that belong to each Voronoi line
        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
          WRITE(*,*) 'Defaulting to region 1.  Is this the correct' //
     '      ' region?'
        ENDIF
        IF(CBBREV(CO,'OUTPUT',3,noco+1,NTCO,N3CO)) THEN
          OUT_FREQ=IFROMC(CO(N3CO+1))
        ELSE
          OUT_FREQ=0
        ENDIF
        VLT=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          IF(OUT_FREQ.NE.0) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            IF(MOD(noelem,OUT_FREQ).EQ.0) THEN
              WRITE(OP_STRING,'('' Done '',I6,'' elements '//
     '          'of '',I6)') noelem,NEELEM(0,nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
CC$          call mp_unsetlock()
          ENDIF
          DO nonbor=1,NNT(nb)
            nbor=NXI(0,nonbor,ne)
            IF(nbor.NE.0) THEN
              VLT=VLT+1
              IF(VLT.GT.VOROLINE_SIZE) THEN
                ERROR='Increase VOROLINE_SIZE'
                GOTO 9999
              ENDIF
              VOROLINE(1,VLT)=ne
              VOROLINE(2,VLT)=nbor
            ENDIF
          ENDDO
        ENDDO

C     Have all the line information we require, now write out
C     the exported nodes ZA (simplex centres) and lines VOROLINE.

C       First the exnode file
        CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        CALL STRING_TRIM(FILE,IBEG,IEND)
        CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'_voro.exnode',
     '    'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
        WRITE(IFILE,'('' Group name:'',A30)')
     '    FILE(IBEG:IEND)//'_voro_node'
        WRITE(IFILE,'('' #Fields=1'')')
        WRITE(IFILE,'('' 1) coordinates, coordinate, rectangular'
     '    //' cartesian, #Components='',I1)') NJ_LOC(NJL_GEOM,0,nr)
        DO i=1,NJ_LOC(NJL_GEOM,0,nr)
          WRITE(IFILE,'(A38)') XYZ_LISTING(i)
        ENDDO
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          WRITE(IFILE,'('' Node: '',I6)') ne
          DO i=1,NJ_LOC(NJL_GEOM,0,nr)
            WRITE(IFILE,'(''      '',E16.6)') ZA(1,i,1,ne)
          ENDDO
        ENDDO
        CALL CLOSEF(IFILE,ERROR,*9999)

C       now the exelem file
        CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'_voro.exelem',
     '    'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
        WRITE(IFILE,'('' Group name: '',A30)')
     '    FILE(IBEG:IEND)//'_voro_elem'
        WRITE(IFILE,'('' Shape.  Dimension=1'')')
        WRITE(IFILE,'('' #Scale factor sets= 1'')')
        WRITE(IFILE,'(''   l.Lagrange, #Scale factors= 2'')')
        WRITE(IFILE,'('' #Nodes= 2'')')
        WRITE(IFILE,'('' #Fields=1'')')
        WRITE(IFILE,'('' 1) coordinates, coordinate, rectangular'
     '    //' cartesian, #Components='',I1)') NJ_LOC(NJL_GEOM,0,nr)
        DO i=1,NJ_LOC(NJL_GEOM,0,nr)
          WRITE(IFILE,'(A5,''  l.Lagrange, no modify, '
     '      //'standard node based.'')') XYZ_BASIS(i)
          WRITE(IFILE,'(''     #Nodes= 2'')')
          WRITE(IFILE,'(''      1.  #Values=1'')')
          WRITE(IFILE,'(''       Value indices:     1'')')
          WRITE(IFILE,'(''       Scale factor indices:   1'')')
          WRITE(IFILE,'(''      2.  #Values=1'')')
          WRITE(IFILE,'(''       Value indices:     1'')')
          WRITE(IFILE,'(''       Scale factor indices:   2'')')
        ENDDO
        DO nl=1,VLT
          WRITE(IFILE,'('' Element: '',I6,'' 0 0'')') nl
          WRITE(IFILE,'(''   Nodes:'')')
          WRITE(IFILE,'(''      '',I6,I6)') VOROLINE(1,nl),
     '      VOROLINE(2,nl)
          WRITE(IFILE,'('' Scale factors:'')')
          WRITE(IFILE,'(''      '',E16.6,''   '',E16.6)') 1.0E0,1.0E0
        ENDDO
        CALL CLOSEF(IFILE,ERROR,*9999)
      ENDIF

      CALL EXITS('EXVORO')
      RETURN
 9999 CALL ERRORS('EXVORO',ERROR)
      CALL EXITS('EXVORO')
      RETURN 1
      END


