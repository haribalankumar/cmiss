      SUBROUTINE OPCOOR(nr,ERROR,*)

C#### Subroutine: OPCOOR
C###  Description:
C###    OPCOOR outputs coordinate data.

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER nr
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj
      CHARACTER OPT3(5)*36,OPT4(7)*36

       DATA OPT3(1)/'rectangular cartesian (x,y,z)       '/,
     '      OPT3(2)/'cylindrical polar (r,theta,z)       '/,
     '      OPT3(3)/'spherical polar (r,theta,phi)       '/,
     '      OPT3(4)/'prolate spheroidal (lambda,mu,theta)'/,
     '      OPT3(5)/'oblate  spheroidal (lambda,mu,theta)'/,
     '      OPT4(1)/'unsymmetric                         '/,
     '      OPT4(2)/'cyl-symmetric (about x (or r) axis) '/,
     '      OPT4(3)/'cyl-symmetric (about y (or z) axis) '/,
     '      OPT4(4)/'sph-symmetric                       '/,
     '      OPT4(5)/'mirror symmetry in x                '/,
     '      OPT4(6)/'mirror symmetry in y                '/,
     '      OPT4(7)/'mirror symmetry in x and y          '/

      CALL ENTERS('OPCOOR',*9999)
      WRITE(OP_STRING,'(/''  The coordinates in region '',I1,'
     '  //''' are '',A,/''  The number of coordinates is '',I1,'
     '  //'/''  The geometry is '',A)')
     '  nr,OPT3(ITYP10(nr)),NJT,OPT4(JTYP4)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''  The origin is at (x,y,z)='',3D12.4)')
     '  (0.0d0,nj=1,3)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(JTYP6.EQ.2) THEN
        WRITE(OP_STRING,'('' The dependent var. coord. system '','
     '    //'''is '',A)') OPT3(ITYP11(nr))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF((JTYP2A.NE.1).AND.(JTYP2B.NE.1).AND.(JTYP2C.NE.1)) THEN
        WRITE(OP_STRING,'(''  All mappings are standard'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(JTYP2A.EQ.1) THEN
        WRITE(OP_STRING,
     '    '(''  Versions are mapped to ensure C0 continuity'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(JTYP2B.EQ.1) THEN
        WRITE(OP_STRING,'(''  Lines have non-standard mappings '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(JTYP2C.EQ.1) THEN
        WRITE(OP_STRING,'(''  Degrees of freedom have non-standard'
     '    //'mappings'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF !jtyp2

      IF(NJ_LOC(njl_fibr,0,nr).GT.0) THEN
        WRITE(OP_STRING,'(''  A fibre direction field is defined'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(OP_STRING,'(2X,I1,'' additional geometric var.s are '','
     '  //'''defined'')') NJ_LOC(njl_fibr,0,nr)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(ITYP10(nr).GT.1) THEN
        IF(JTYP10.EQ.1) THEN
          IF(ITYP10(nr).EQ.2.OR.ITYP10(nr).EQ.3) THEN
            FORMAT='(''  Interpolation in the radial dir.n is in r'')'
          ELSE IF(ITYP10(nr).EQ.4) THEN
            FORMAT='(''  Interpolation is in Lambda'')'
          ENDIF
          WRITE(OP_STRING,FORMAT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          IF(ITYP10(nr).EQ.2) THEN
            FORMAT='(''  Interpolation in the radial dir.n is in r^2'')'
          ELSE IF(ITYP10(nr).EQ.3) THEN
            FORMAT='('' Interpolation in the radial dir.n is in r^3'')'
          ELSE IF(ITYP10(nr).EQ.4.AND.JTYP10.EQ.2) THEN
            FORMAT='(''  Interpolation is in Focus^2.Sinh^2(Lambda)'')'
          ELSE IF(ITYP10(nr).EQ.4.AND.JTYP10.EQ.3) THEN
            FORMAT='(''  Interpolation is in Focus^3.Cosh(Lambda).'','
     '        //'''Sinh^2(Lambda)'')'
          ENDIF
          WRITE(OP_STRING,FORMAT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF !ityp10

      CALL EXITS('OPCOOR')
      RETURN
 9999 CALL ERRORS('OPCOOR',ERROR)
      CALL EXITS('OPCOOR')
      RETURN 1
      END


c     SUBROUTINE OPCORN(NPNODE,nr,NVHP,ERROR,*)
c
C#### Subroutine: OPCORN
C###  Description:
C###    OPCORN outputs corner node data (for BE problems).

c     IMPLICIT NONE
c     INCLUDE 'cmiss$reference:b01.cmn'
c     INCLUDE 'cmiss$reference:b14.cmn'
c     INCLUDE 'cmiss$reference:cbdi02.cmn'
c     INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
c     INTEGER NPNODE(0:NP_R_M,0:NRM),nr,NVHP(NHM,NPM,NCM)
c     CHARACTER ERROR*(*)
!     Local Variables
c     INTEGER nc,NCOR,NEDGE,nh,nonode,np

c     CALL ENTERS('OPCORN',*9999)

C!!!!! Needs to be rewritten

c      nh=1 !Default value
c      NCOR=0
c      NEDGE=0
c      DO nonode=1,NPNODE(0,nr)
c        np=NPNODE(nonode,nr)
c        IF(NVHP(nh,np,0).EQ.3)THEN
c          IF(NJT.EQ.2)THEN
c            NCOR=NCOR+1 !A corner node
c          ELSE
c            NEDGE=NEDGE+1
c          ENDIF
c        ELSEIF(NVHP(nh,np,0).GT.3)THEN
c          NCOR=NCOR+1 !A corner node
c        ENDIF
c      ENDDO
c      IF(NJT.EQ.2) THEN
c        WRITE(OP_STRING,'(''  The number of corner nodes is'',I3)')
c     '    NCOR
c        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c      ELSE IF(NJT.EQ.3) THEN
c        WRITE(OP_STRING,'(''  The number of edge nodes ='',I3)')NEDGE
c        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c        WRITE(OP_STRING,'(''  The number of 3d corner nodes ='',I3)')
c     '    NCOR
c        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c      ENDIF
c      IF((NJT.EQ.2).AND.(NCOR.GT.0)) THEN
c        WRITE(OP_STRING,'(''   Corner node numbers'')')
c        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c      ELSE IF((NJT.EQ.3).AND.(NEDGE.GT.0)) THEN
c        WRITE(OP_STRING,'(''   Edge node numbers'')')
c        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c      ENDIF
c      DO nonode=1,NPNODE(0,nr)
c        np=NPNODE(nonode,nr)
c        IF(NVHP(nh,np,0).EQ.3)THEN
c          WRITE(OP_STRING,'(''    Node  '',I3)') np
c          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c          WRITE(OP_STRING,'(''     Shared between elements '','
c     '      //'4(I3,2X))') (NVHP(nh,np,nc),nc=2,NVHP(nh,np,0))
c          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c        ENDIF
c      ENDDO
c      IF((NJT.EQ.3).AND.(NCOR.GT.0)) THEN
c        WRITE(OP_STRING,'(''   Corner node numbers'')')
c        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c        DO nonode=1,NPNODE(0,nr)
c          np=NPNODE(nonode,nr)
c          IF(NVHP(nh,np,0).GT.3)THEN
c            WRITE(OP_STRING,'(''    Node  '',I3)') np
c            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c            WRITE(OP_STRING,'(''     Shared between elements '','
c     '        //'4(I3,2X))') (NVHP(nh,np,nc),nc=2,NVHP(nh,np,0))
c            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c          ENDIF
c        ENDDO
c      ENDIF

c     CALL EXITS('OPCORN')
c     RETURN
c9999 CALL ERRORS('OPCORN',ERROR)
c     CALL EXITS('OPCORN')
c     RETURN 1
c     END


