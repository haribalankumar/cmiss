      SUBROUTINE TREE_GEOM(ELEM,ELEMLIST,ELEMSMOOTH,ESTART,GENERATIONS,
     '  GEOMFRAC,MAX_BRANCH,NBJ,NOELEMS,NXI,ntt,POINTS,WEIGHTS,
     '  XJPOWER,ERROR,*)

C#### Subroutine: TREE_GEOM
C###  Description:
C###    TREE_GEOM defines the areas in rectangular cartesian
C###    coordinates which points will be placed within for trees
C###    grown with IPMESH8.
C***  Created by Martin Buist December 1996

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'

!     Parameter list
      INTEGER GENERATIONS,MAX_BRANCH,NBJ(NJM,NEM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),ntt
      REAL*8 GEOMFRAC(MAX_BRANCH,2,3),WEIGHTS(MAX_BRANCH),
     '  XJPOWER(MAX_BRANCH,3)
      CHARACTER ERROR*(*)
      LOGICAL ELEM,POINTS(MAX_BRANCH)
!     Local variables
      INTEGER ELEMLIST(MAX_BRANCH,4),ESTART,IBEG,ICHAR,IEND,INFO,
     '  nb,ne,nex,ney,nez,ne1,ne2,ni,nj,NOELEMS(0:6),NOQUES,NUMELEM
      CHARACTER CHAR1*10,CHAR2*10
      LOGICAL FILEIP,ELEMSMOOTH

      CALL ENTERS('TREE_GEOM',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      ELEM=.TRUE.
      ELEMSMOOTH=.TRUE.
      ESTART=1
      GENERATIONS=4
      DO nj=1,6
        NOELEMS(nj)=0
      ENDDO !nj

      DO ne=1,MAX_BRANCH
        ELEMLIST(ne,1)=0
        DO nj=1,NJT
          GEOMFRAC(ne,1,nj)=0.0d0
          GEOMFRAC(ne,2,nj)=1.0d0
          XJPOWER(ne,nj)=0.0d0
        ENDDO !nj
        WEIGHTS(ne)=1.0d0
      ENDDO !ne

      WRITE(CHAR2,'(I1)') ntt
      CALL STRING_TRIM(CHAR2,IBEG,IEND)
      IF(IOTYPE.EQ.3) THEN
        IF(ELEM) THEN
          ADATA(1)='Y'
        ELSE
          ADATA(1)='N'
        ENDIF
      ENDIF
      FORMAT='($,'' Grow tree only in elements for tree '
     '  //CHAR2(IBEG:IEND)//' [Y]? '',A)'
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        IF(ADATA(1).EQ.'Y') THEN
          ELEM=.TRUE.
        ELSE
          ELEM=.FALSE.
        ENDIF
      ENDIF

      IF(ELEM) THEN
        IDEFLT(1)=1
        FORMAT='($,'' Enter the starting element number for tree '
     '    //CHAR2(IBEG:IEND)//' [1]:'',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=ESTART
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NEM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ESTART=IDATA(1)

        DO nj=1,NJT
          WRITE(CHAR1,'(I1)') nj
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
          IDEFLT(1)=0
          FORMAT='($,'' Enter # elements in +'//CHAR1(IBEG:IEND)//
     '    ' direction [0]:'',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=NOELEMS(nj)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NOELEMS(nj)=IDATA(1)

          IDEFLT(1)=0
          FORMAT='($,'' Enter # elements in -'//CHAR1(IBEG:IEND)//
     '    ' direction [0]:'',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=NOELEMS(nj+3)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NOELEMS(nj+3)=IDATA(1)
          CALL ASSERT(((NOELEMS(nj+3)+NOELEMS(nj)).LE.10),
     '      '>>ERROR Increase dimensions of GEOMFRAC in ipmesh8',
     '      ERROR,*9999)
        ENDDO !nj

! Find the bottom/left/near element to start from
        ne=ESTART
        nb=NBJ(1,ne)
        DO ni=1,NIT(nb)
          IF(NOELEMS(ni+3).GT.0) THEN
            DO ne1=1,NOELEMS(ni+3)
              CALL ASSERT(NXI(-ni,1,ne).GT.0,
     '          '>>ERROR no adjacent element',ERROR,*9999)
              ne=NXI(-ni,1,ne)
            ENDDO !ne1
          ENDIF
          NOELEMS(ni)=NOELEMS(ni)+NOELEMS(ni+3)
        ENDDO !ni

        ne1=ne
        ne2=ne
        NUMELEM=0
        NOELEMS(0)=(NOELEMS(3)+1)*(NOELEMS(2)+1)*(NOELEMS(1)+1)
        CALL ASSERT(NOELEMS(0).LE.MAX_BRANCH,
     '    '>>ERROR too many elements selected ',ERROR,*9999)
        DO nez=1,NOELEMS(3)+1
          DO ney=1,NOELEMS(2)+1
            DO nex=1,NOELEMS(1)+1
              NUMELEM=NUMELEM+1
              ELEMLIST(NUMELEM,1)=ne
              ELEMLIST(NUMELEM,2)=nex-1
              ELEMLIST(NUMELEM,3)=ney-1
              ELEMLIST(NUMELEM,4)=nez-1
              ne=NXI(1,1,ne)
            ENDDO !nex
            ne1=NXI(2,1,ne1)
            ne=ne1
          ENDDO !ney
          ne2=NXI(3,1,ne2)
          ne=ne2
          ne1=ne2
        ENDDO !nez
      ELSE
        NUMELEM=1
      ENDIF

      DO ne=1,NUMELEM
        WRITE(OP_STRING,'('' Element number '',I4)') ELEMLIST(ne,1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)

        IF(IOTYPE.EQ.3) THEN
          IF(POINTS(ne)) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        FORMAT='($,'' Put points in current element [Y]? '',A)'
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            POINTS(ne)=.TRUE.
          ELSE
            POINTS(ne)=.FALSE.
          ENDIF
        ENDIF
        IF(.NOT.ELEM) POINTS(ne)=.TRUE.
        IF(POINTS(ne)) THEN
          IF(ELEM) THEN
            RDEFLT(1)=1.0d0
            FORMAT='($,'' Enter element weight [1]: '',F8.6)'
            IF(IOTYPE.EQ.3) RDATA(1)=WEIGHTS(ne)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) WEIGHTS(ne)=RDATA(1)
          ENDIF
          DO nj=1,NJT
            RDEFLT(1)=0.0d0
            WRITE(CHAR1,'(I1)') nj
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            FORMAT='($,'' Enter global coordinate '//CHAR1(IBEG:IEND)//
     '        ' power scaling factor [0]:'',F8.6)'
            IF(IOTYPE.EQ.3) RDATA(1)=XJPOWER(ne,nj)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) XJPOWER(ne,nj)=RDATA(1)

            RDEFLT(1)=0.0d0
            FORMAT='($,'' Enter start fraction of '//CHAR1(IBEG:IEND)//
     '        ' coordinate geometry to grow tree in [0]:'',F8.6)'
            IF(IOTYPE.EQ.3) RDATA(1)=GEOMFRAC(ne,1,nj)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) GEOMFRAC(ne,1,nj)=RDATA(1)

            RDEFLT(1)=1.0d0
            FORMAT='($,'' Enter end fraction of '//CHAR1(IBEG:IEND)//
     '        ' coordinate geometry to grow tree in [1]:'',F8.6)'
            IF(IOTYPE.EQ.3) RDATA(1)=GEOMFRAC(ne,2,nj)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) GEOMFRAC(ne,2,nj)=RDATA(1)
          ENDDO !nj
        ENDIF
      ENDDO !ne

      IDEFLT(1)=4
      FORMAT='($,'' Enter the number of generations [4]:'',I2)'
      IF(IOTYPE.EQ.3) IDATA(1)=GENERATIONS
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) GENERATIONS=IDATA(1)

      IF(IOTYPE.EQ.3) THEN
        IF(ELEMSMOOTH) THEN
          ADATA(1)='Y'
        ELSE
          ADATA(1)='N'
        ENDIF
      ENDIF
      FORMAT='($,'' Do you want to use element smoothing [Y]? '',A)'
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        IF(ADATA(1).EQ.'Y') THEN
          ELEMSMOOTH=.TRUE.
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Using element size smoothing '')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE
          ELEMSMOOTH=.FALSE.
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Using one element per branch '')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('TREE_GEOM')
      RETURN
 9999 CALL ERRORS('TREE_GEOM',ERROR)
      CALL EXITS('TREE_GEOM')
      RETURN 1
      END



