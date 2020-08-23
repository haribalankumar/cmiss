      SUBROUTINE IPCOOR(nr,ERROR,*)

C#### Subroutine: IPCOOR
C###  Description:
C###    IPCCOR inputs coordinate system.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
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
      INTEGER IBEG,ICHAR,IEND,INFO,MAXOPT,nj,njj1,njj2,NJTOLD,NOQUES,nrr
      CHARACTER CHAR1*1,CHAR2*2
      LOGICAL FILEIP,NONSTANDARD

      CALL ENTERS('IPCOOR',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      WRITE(CHAR2,'(I2)') nr
      CALL STRING_TRIM(CHAR2,IBEG,IEND)
      FORMAT='('' The global coordinates for region '
     '  //CHAR2(IBEG:IEND)//' are [1]: '''//
     '  '/''   (1) rectangular cartesian (x,y,z)'''//
     '  '/''   (2) cylindrical polar (r,theta,z)'''//
     '  '/''   (3) spherical polar (r,theta,phi)'''//
     '  '/''   (4) prolate spheroidal (lambda,mu,theta)'''//
     '  '/''   (5) oblate  spheroidal (lambda,mu,theta)'''//
     '  '/$,''    '',I1)'
      MAXOPT=5
      IF(IOTYPE.EQ.3) IDATA(1)=ITYP10(nr)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,MAXOPT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) ITYP10(nr)=IDATA(1)

      NJTOLD=NJT
      IDEFLT(1)=NJT
      WRITE(CHAR1,'(I1)') NJT
      FORMAT='($,'' Enter the number of global coordinates ['
     '  //CHAR1//']: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NJT
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NJT=IDATA(1)

      IF(NJT.NE.NJTOLD) THEN
        CALL ASSERT(NJ_LOC(NJL_FIBR,0,nr).EQ.0,
     '    '>>Cancel fibres and sheets',ERROR,*9999)
        CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).EQ.0,
     '    '>>Cancel field',ERROR,*9999)
      ENDIF

      DO nj=1,NJT
        NJ_LOC(NJL_GEOM,nj,nr)=NJ
        NJ_TYPE(nj,1)=NJL_GEOM
        NJ_TYPE(nj,2)=NJ
      ENDDO
      NJ_LOC(NJL_GEOM,0,nr)=NJT
      IF(NJT.GT.NJ_LOC(NJL_GEOM,0,0)) NJ_LOC(NJL_GEOM,0,0)=NJT
      NJ_LOC(0,0,nr)=0
      DO njj1=1,3
        DO njj2=1,NJ_LOC(njj1,0,nr)
          nj=NJ_LOC(njj1,njj2,nr)
          IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=NJ
        ENDDO
      ENDDO
      DO nrr=1,NRT
        IF(NJ_LOC(0,0,nrr).GT.NJ_LOC(0,0,0))
     '    NJ_LOC(0,0,0)=NJ_LOC(0,0,nrr)
      ENDDO !nrr

      FORMAT='($,'' Do you want to specify another coord.'
     '  //' system for dependent variables [N]? '',A)'
      IF(IOTYPE.EQ.3) THEN
        IF(JTYP6.EQ.1) THEN
          ADATA(1)='N'
        ELSE IF(JTYP6.EQ.1) THEN
          ADATA(1)='Y'
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        IF(ADATA(1).EQ.'N') THEN
          JTYP6=1
        ELSE
          JTYP6=2
        ENDIF
      ENDIF

      IF(JTYP6.EQ.1) THEN
        ITYP11(nr)=ITYP10(nr)
      ELSE IF(JTYP6.EQ.2) THEN
        FORMAT='('' The global dependent variable coords are [1]: '''//
     '    '/''   (1) rectangular cartesian'''//
     '    '/''   (2) cylindrical polar'''//
     '    '/''   (3) spherical polar'''//
     '    '/''   (4) prolate spheroidal'''//
     '    '/''   (5) oblate  spheroidal'''//
     '    '/$,''    '',I1)'
        MAXOPT=5
        IF(IOTYPE.EQ.3) IDATA(1)=ITYP11(nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,MAXOPT,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITYP11(nr)=IDATA(1)
      ENDIF

      IF(ITYP10(nr).GT.1.AND.NJT.GT.1) THEN
        IF(IOTYPE.EQ.1) THEN
          IF(ITYP10(nr).EQ.2.OR.ITYP10(nr).EQ.3) THEN
            FORMAT='('' Note: the Xi(1) coord must be in the direction'
     '        //' of increasing THETA'')'
            WRITE(OP_STRING,FORMAT)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(ITYP10(nr).EQ.4.OR.ITYP10(nr).EQ.5) THEN
            FORMAT='('' Note: the Xi(1) coord must be in the direction'
     '        //' of decreasing THETA'')'
            WRITE(OP_STRING,FORMAT)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          IF(ITYP10(nr).EQ.3) THEN
            FORMAT='('' ..and the Xi(2) coord must be in the direction'
     '        //' of increasing PHI'')'
            WRITE(OP_STRING,FORMAT)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(ITYP10(nr).EQ.4.OR.ITYP10(nr).EQ.5) THEN
            FORMAT='('' ..and the Xi(2) coord must be in the direction'
     '        //' of increasing MU'')'
            WRITE(OP_STRING,FORMAT)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        IF(ITYP10(nr).LT.4) THEN
          WRITE(CHAR1,'(I1)') ITYP10(nr)
          FORMAT='(/'' Specify whether radial interpolation'//
     '      ' is in [1]:'''//
     '      '/''   (1) r'''//
     '      '/''   (2) r^'//CHAR1(1:1)//''''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=JTYP10
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) JTYP10=IDATA(1)
        ELSE IF(ITYP10(nr).EQ.4) THEN
          FORMAT='(/'' Specify whether interpolation'//
     '      ' is in [1]: '''//
     '      '/''   (1) Lambda'''//
     '      '/''   (2) Focus^2.Sinh^2(Lambda)'''//
     '      '/''   (3) Focus^3.Cosh(Lambda).Sinh^2(Lambda)'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=JTYP10
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) JTYP10=IDATA(1)
        ENDIF
      ENDIF

      IF(NJT.LE.2) THEN
        FORMAT='('' The geometry is [1]:'''//
     '    '/''   (1) unsymmetric'''//
     '    '/''   (2) cylindrically symmetric about x '//
     '    '(or r for cyl.polar)'''//
     '    '/''   (3) cylindrically symmetric about y '//
     '    '(or z for cyl.polar)'''//
     '    '/''   (4) spherically   symmetric'''//
     '    '/''   (5) mirror symmetry in x'''//
     '    '/''   (6) mirror symmetry in y'''//
     '    '/''   (7) mirror symmetry in x and y'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=JTYP4
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,7,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) JTYP4=IDATA(1)
        CALL ASSERT(JTYP4.LT.5,'>>Sorry, these options are currently '
     '    //'not supported',ERROR,*9999)
      ELSE
        JTYP4=1   !unsymmetric for 3D elements
      ENDIF

C KAT 2002-08-21: This is not used for anything.
      FORMAT='($,'' Enter x,y,z origin of coords'
     '  //' relative to region 0 [0,0,0]:'',3E13.5)'
      IF(IOTYPE.EQ.3) THEN
        DO nj=1,3
          RDATA(nj)=0.0d0
        ENDDO
      ENDIF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,3,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

      FORMAT='($,'' Are there any non-standard mappings [N]? '',A)'
      IF(IOTYPE.EQ.3) THEN
        IF(JTYP2A.EQ.1) THEN
          ADATA(1)='Y'
          NONSTANDARD=.TRUE.
        ELSE IF(JTYP2B.EQ.1) THEN
          ADATA(1)='Y'
          NONSTANDARD=.TRUE.
        ELSE IF(JTYP2C.EQ.1) THEN
          ADATA(1)='Y'
          NONSTANDARD=.TRUE.
        ELSE
          ADATA(1)='N'
          NONSTANDARD=.FALSE.
        ENDIF
      ENDIF
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        IF(ADATA(1).EQ.'N') THEN
          NONSTANDARD=.FALSE.
        ELSE IF(ADATA(1).EQ.'Y') THEN
          NONSTANDARD=.TRUE.
        ENDIF
      ENDIF

      IF(NONSTANDARD) THEN
        FORMAT='($,''    in versions to ensure C0 continuity [N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(JTYP2A.EQ.1) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            JTYP2A=1
            ITYP21(nr)=1
          ELSE
            JTYP2A=0
          ENDIF
        ENDIF

        FORMAT='($,''    in lines [N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(JTYP2B.EQ.1) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            JTYP2B=1
          ELSE
            JTYP2B=0
          ENDIF
        ENDIF

        FORMAT='($,''    in degrees of freedom for hanging '
     '    //'nodes [N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(JTYP2C.EQ.1) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            JTYP2C=1
          ELSE
            JTYP2C=0
          ENDIF
        ENDIF
      ENDIF

C CS 28/6/99 split this into JTYP2A,JTYP2B,JTYP2C above
C      IF(NONSTANDARD) THEN
C        FORMAT='('' The mapping is non-standard in [1]:'''//
C     '    '/''   (1) versions to ensure C0 continuity'''//
C     '    '/''   (2) lines'''//
C     '    '/''   (3) degrees of freedom - hanging nodes'''//
C     '    '/$,''    '',I1)'
C        IF(IOTYPE.EQ.3) THEN
C          IDATA(1)=JTYP2
C        ENDIF
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) THEN
C          JTYP2=IDATA(1)
C        ENDIF
C      ENDIF

      CALL EXITS('IPCOOR')
      RETURN
 9999 CALL ERRORS('IPCOOR',ERROR)
      CALL EXITS('IPCOOR')
      RETURN 1
      END


