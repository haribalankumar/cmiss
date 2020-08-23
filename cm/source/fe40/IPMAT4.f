      SUBROUTINE IPMAT4(NBJ,NEELEM,NELIST,NPNODE,nr,NW,nx,
     '  CE,CP,YG,ERROR,*)

C#### Subroutine: IPMAT4
C###  Description:
C###    IPMAT4 inputs material parameters.

C**** IMT(ie) is index of material type for element type ie.
C****   (1,2,3,4,5 for isotropic/trans.isotropic (wrt Xi_1/2)
C****   /orthotropic/anisotropic(special case))
C**** ILT(ie,nr,nx) is total number of parameters required for ie;
C****   these comprise material,geometric & load params, where:
C**** NMP(ie) is number of material  parameters reqd for ie.
C****   =NMPP(ie,IMT(ie))
C****   +NTPP(ie,IMT(ie)) if thermal effects are included
C**** NGP(ie) is number of geometric parameters reqd for ie
C**** NLP(ie) is number of loading   parameters reqd for ie
C**** ILP(il,ie,nr,nx),il=1,ILT(ie,nr,nx), records whether eqtn
C**** params are:
C****   (1) Constant spatially  - value  in CE(il,ne)
C****   (2) Piecewise constant  - values in CE(il,ne)
C****   (3) Piecewise linear    - values in CP(il,np)
C****   (4) Defined by Gauss pts- values in YG(5,ng,ne)
C**** For nonlinear probs in which material parameters are incremented
C****   the increments are stored in CE(il+ILT(ie,nr,nx),ne)
C**** NMB(il,nx) is number of material basis used for nodal interpol.n
C**** KTYP43 is 0..3 for type of thermal strain effect
C**** KTYP45 is 1..4  for type of beam cross-section

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp40.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'titl40.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NPNODE(0:NP_R_M,0:NRM),nr,NW(NEM,3),nx
      REAL*8 CE(NMM,NEM),CP(NMM,NPM),YG(NIYGM,NGM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,IBEG1,IBEG2,IBEG3,ICHAR,ie,IEND,IEND1,IEND2,IEND3,
     '  ij,il,INFO,n,nb,n1,
     '  ne,ng,NGPP(4),nj,noelem,nonode,NOQUES,np,NTPP(12,5)
      REAL*8 LAST_INCREMENT,LAST_VALUE,SUM
      CHARACTER CHAR1*50,CHAR2*10,CHAR3*10,TITLE*50
      LOGICAL FILEIP,FIRST,INLIST

      DATA NTPP/12*2,12*3,12*3,12*4,12*2/ !thermal strain parameters
      DATA NGPP/2,3,2,3/

      CALL ENTERS('IPMAT4',*9999)
      CALL ASSERT(ITYP2(nr,nx).NE.0,
     '  ' >>Solution type has not been defined',ERROR,*9999)
      CALL ASSERT(ITYP2(nr,nx).EQ.1,' >>Solution type is incorrect',
     '  ERROR,*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FIRST=.TRUE.
      DO ie=1,12
        IF(ETYP(ie).AND.ie.NE.10) THEN
          IF(.NOT.FIRST) THEN
            FORMAT='(1X)'
            CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '        0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
c           CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
          ENDIF
          FIRST=.FALSE.

          CHAR1=TITL42(ie)
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
          IF(ie.LE.2) THEN      !truss or batten elements
            IMT(ie)=1
            FORMAT='('' Truss, cable or batten elements:'')'
            CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          ELSE IF(ie.EQ.3) THEN !beam elements
            IMT(ie)=1
            FORMAT='('' Cross-section type for '
     '        //CHAR1(IBEG:IEND)//' is [1]: '''//
     '        '/''  (1) Rectangular solid'''//
     '        '/''  (2) Rectangular tube'''//
     '        '/''  (3) Ellipsoidal solid'''//
     '        '/''  (4) Ellipsoidal tube'''//
     '        '/$,''   '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP45
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '        ICHAR,IDATA,IONE,1,4,LDATA,LDEFLT,
     '        RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP45=IDATA(1)
            NGP(ie)=NGPP(KTYP45)
          ELSE IF(ie.EQ.4) THEN !link elements
            IMT(ie)=1
          ELSE                  !all other elements
            FORMAT='('' Specify whether material for '
     '        //CHAR1(IBEG:IEND)//' is [1]: '''//
     '        '/''  (1) Isotropic'''//
     '        '/''  (2) Transversely isotropic wrt Xi_1'''//
     '        '/''  (3) Transversely isotropic wrt Xi_2'''//
     '        '/''  (4) Orthotropic '''//
     '        '/''  (5) Anisotropic (special case)'''//
     '        '/$,''   '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=IMT(ie)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,5,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) IMT(ie)=IDATA(1)
            IF(IMT(ie).GT.1) THEN
              CALL ASSERT(NJ_LOC(njl_fibr,0,nr).GT.0,
     '          ' >>Fibre field has not been defined',ERROR,*9999)
            ENDIF
          ENDIF !ie

          NMP(ie)=NMPP(ie,IMT(ie))
          IF(KTYP43.GT.0) THEN !include thermal strain parameters
            NMP(ie)=NMP(ie)+NTPP(ie,IMT(ie))
          ENDIF
          ILT(ie,nr,nx)=NMP(ie)+NGP(ie)+NLP(ie)+1  !last is density
          CALL ASSERT(ILT(ie,nr,nx).LE.NMM,' >>Increase NMM',
     '      ERROR,*9999)
          DO il=1,ILT(ie,nr,nx)
            IF(il.LE.NMPP(ie,IMT(ie))) THEN        !elastic matl params
              TITLE=TITL46(il,IMT(ie))
            ELSE IF(il.LE.NMP(ie)) THEN            !thermal matl params
              TITLE=TITL41(il-NMPP(ie,IMT(ie)),IMT(ie))
            ELSE IF(il.LE.NMP(ie)+NGP(ie)) THEN    !geometric params
              TITLE=TITL47(il-NMP(ie),ie)
            ELSE IF(il.LE.NMP(ie)+NGP(ie)+NLP(ie)) THEN !element loads
              TITLE=TITL48(il-NMP(ie)-NGP(ie),ie)
            ELSE                                   !mass density
              TITLE='Density (kg/m^3)'
            ENDIF
            CHAR1=TITLE
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            FORMAT='(/'' Specify whether '//CHAR1(IBEG:IEND)//
     '        ' is [1]: '''//
     '        '/''  (1) Constant spatially'''//
     '        '/''  (2) Piecewise constant (defined by elements)'''//
     '        '/''  (3) Piecewise linear   (defined by nodes)'''//
     '        '/''  (4) Defined by Gauss points'''//
     '        '/$,''   '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=ILP(il,ie,nr,nx)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ILP(il,ie,nr,nx)=IDATA(1)
            IF(il.EQ.1.AND.IMT(ie).EQ.1) THEN
              FORMAT='('' NB. Approx: Steel E=200 GPa,'//
     '         ' Timber E=10 GPa'')'
              CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '          ERROR,*9999)
c             CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
            ENDIF

            IF(ILP(il,ie,nr,nx).EQ.1) THEN !constant spatially
              FORMAT='($,'' The value is [0.0]: '',E11.4)'
              IF(IOTYPE.EQ.3) RDATA(1)=CE(il,1)
CGMH 30/6/95.  Must be able to enter negative temperatures
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RZERO,
     '            -RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IF(NW(ne,1).EQ.ie) CE(il,ne)=RDATA(1)
                ENDDO
              ENDIF
              IF(KTYP14.GT.0) THEN
                FORMAT='($,'' The increment is [0.0]: '',E11.4)'
                IF(IOTYPE.EQ.3) RDATA(1)=CE(il,1)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RZERO,
     '            -RMAX,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    IF(NW(ne,1).EQ.ie) CE(il+ILT(ie,nr,nx),ne)=RDATA(1)
                  ENDDO
                ENDIF
              ENDIF

            ELSE IF(ILP(il,ie,nr,nx).EQ.2) THEN !defined by elements
C             Initialise with reasonable values
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(NW(ne,1).EQ.ie) THEN
                  CE(il,ne)=0.0D0
                  IF(KTYP14.GT.0) THEN
                    CE(il+ILT(ie,nr,nx),ne)=0.0D0
                  ENDIF
                ENDIF
              ENDDO
C GMH 5/12/95 Copied from fe50.f IPINI5
              LAST_VALUE=0.0D0
              LAST_INCREMENT=0.0D0
 6720         FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
 6760         CDATA(1)='ELEMENTS' !for use with group input
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '          0,NET(nr),
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IDATA(1).NE.0) THEN !not default exit
                NELIST(0)=IDATA(0)
                DO n=1,IDATA(0)
                  NELIST(n)=IDATA(n)
                  ne=IDATA(n)
                  IF(.NOT.INLIST(ne,NEELEM(1,nr),
     '              NEELEM(0,nr),N1)) THEN
                    WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '                //'in the current region'')') ne
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    GOTO 6760
                  ENDIF
                  IF(NW(ne,1).NE.ie) THEN
                    WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '                //'of the correct type'')') ne
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    GOTO 6760
                  ENDIF
                ENDDO !n

C               Get type for first element in group
                ne=NELIST(1) !rest of group filled at end of loop
                RDEFLT(1)=LAST_VALUE
                WRITE(CHAR2,'(E10.3)') RDEFLT(1)
                CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
                FORMAT='($,'' The value is ['//CHAR2(1:IEND2)
     '            //']:'',E11.4)'
                IF(IOTYPE.EQ.3) RDATA(1)=CE(il,ne)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '            -RMAX,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) LAST_VALUE=RDATA(1)
                IF(KTYP14.GT.0) THEN
                  RDEFLT(1)=LAST_INCREMENT
                  WRITE(CHAR2,'(E10.3)') RDEFLT(1)
                  CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
                  FORMAT='($,'' The increment in element is ['
     '              //CHAR2(1:IEND2)//']:'',E11.4)'
                  IF(IOTYPE.EQ.3) RDATA(1)=CE(il+ILT(ie,nr,nx),ne)
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '              -RMAX,RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) LAST_INCREMENT=RDATA(1)
                ENDIF

C               Apply type to all of elements group
                DO n=1,NELIST(0)
                  ne=NELIST(n)
                  CE(il,ne)=LAST_VALUE
                  IF(KTYP14.GT.0) THEN
                    CE(il+ILT(ie,nr,nx),ne)=LAST_INCREMENT
                  ENDIF
                ENDDO !n

                GO TO 6720 !for more elements
              ENDIF !idata(1).ne.0
C              DO noelem=1,NEELEM(0,nr)
C                ne=NEELEM(noelem,nr)
C                IF(NW(ne,1).EQ.ie) THEN
C                  WRITE(CHAR1,'(I5)') ne
C                  CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
C                  IF(ne.EQ.1) THEN
C                    RDEFLT(1)=0.0d0
C                    CHAR2=' 0.000E+00'
C                  ELSE
C                    RDEFLT(1)=CE(il,ne-1)
C                    WRITE(CHAR2,'(E10.3)') RDEFLT(1)
C                  ENDIF
C                  CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
C                  FORMAT='($,'' The value in element '//
C     '              CHAR1(IBEG1:IEND1)//' is ['//CHAR2(1:IEND2)//
C     '              ']:'',E11.4)'
C                  IF(IOTYPE.EQ.3) RDATA(1)=CE(il,ne)
C                  CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '              FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '              IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
C     '              -RMAX,RMAX,INFO,ERROR,*9999)
C                  IF(IOTYPE.NE.3) CE(il,ne)=RDATA(1)
C                  IF(KTYP14.GT.0) THEN
C                    FORMAT='($,'' The increment in element '//
C     '                CHAR1(IBEG1:IEND1)//' is ['//CHAR2(1:IEND2)//
C     '                ']:'',E11.4)'
C                    IF(IOTYPE.EQ.3) RDATA(1)=CE(il,ne)
C                    CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '                FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '                IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
C     '                -RMAX,RMAX,INFO,ERROR,*9999)
C                    IF(IOTYPE.NE.3) CE(il+ILT(ie,nr,nx),ne)=RDATA(1)
C                  ENDIF
C                ENDIF
C              ENDDO

            ELSE IF(ILP(il,ie,nr,nx).EQ.3) THEN !defined by nodes
              FORMAT='($,'' Enter basis type number [1]: '',I1)'
              IF(IOTYPE.EQ.3) IDATA(1)=NMB(il,ie,nx)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '          ICHAR,IDATA,IONE,1,9,LDATA,LDEFLT,
     '          RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) NMB(il,ie,nx)=IDATA(1)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                WRITE(CHAR1,'(I3)') np
                CALL STRING_TRIM(CHAR1,IBEG,IEND)
                IF(np.EQ.1) THEN
                  RDEFLT(1)=0.0d0
                  CHAR2=' 0.000E+00'
                ELSE
                  RDEFLT(1)=CP(il,np-1)
                  WRITE(CHAR2,'(E10.3)') RDEFLT(1)
                ENDIF
                CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
                FORMAT='($,'' The value at node '//CHAR1(IBEG:IEND)//
     '            ' is ['//CHAR2(1:IEND2)//']: '',E11.4)'
                IF(IOTYPE.EQ.3) RDATA(1)=CP(il,np)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '            -RMAX,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) CP(il,np)=RDATA(1)
                IF(KTYP14.GT.0) THEN
                  FORMAT='($,'' The increment at node '
     '              //CHAR1(IBEG:IEND)
     '              //' is ['//CHAR2(1:IEND2)//']: '',E11.4)'
                  IF(IOTYPE.EQ.3) RDATA(1)=CP(il,np)
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '              -RMAX,RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) CP(il+ILT(ie,nr,nx),np)=RDATA(1)
                ENDIF
              ENDDO

            ELSE IF(ILP(il,ie,nr,nx).EQ.4) THEN !defined by Gauss points
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(NW(ne,1).EQ.ie) THEN
                  WRITE(CHAR1,'(I3)') ne
                  CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
                  nb=NBJ(1,ne)
                  DO ng=1,NGT(nb)
                    WRITE(CHAR2,'(I2)') ng
                    CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
                    IF(ne.EQ.1.AND.ng.EQ.1) THEN
                      RDEFLT(1)=0.0d0
                      CHAR3=' 0.000E+00'
                    ELSE IF(ng.EQ.1) THEN
                      RDEFLT(1)=YG(5,NGT(nb),ne-1)
                      WRITE(CHAR3,'(E10.3)') RDEFLT(1)
                    ELSE IF(ng.GT.1) THEN
                      RDEFLT(1)=YG(5,ng-1,ne)
                      WRITE(CHAR3,'(E10.3)') RDEFLT(1)
                    ENDIF
                    CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
                    FORMAT=
     '                '($,'' The value in element '//CHAR1(IBEG1:IEND1)
     '                //' Gauss pt '//CHAR2(IBEG2:IEND2)
     '                //' is ['//CHAR3(1:IEND3)//']:'',E11.4)'
                    IF(IOTYPE.EQ.3) RDATA(1)=YG(5,ng,ne)
                    CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '                -RMAX,RMAX,INFO,ERROR,*9999)
                    IF(IOTYPE.NE.3) YG(5,ng,ne)=RDATA(1)
                  ENDDO
                ENDIF
              ENDDO !noelem
            ENDIF !ILP(il,ie,nr,nx)
          ENDDO !il
        ENDIF !ETYP(ie)
      ENDDO !ie

      IF(ETYP(3).AND.NJT.EQ.3) THEN !3D beams
        DO nj=1,2
          CHAR1=TITL49(nj)
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
          FORMAT='(/'' Specify the 3 '//CHAR1(IBEG:IEND)//''','
     '      //'/'' (Ratios only are needed): '')'
          CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '      0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '      ERROR,*9999)
c         CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(NW(ne,1).EQ.3) THEN !beam element
              WRITE(CHAR1,'(I3)') ne
              CALL STRING_TRIM(CHAR1,IBEG,IEND)
              FORMAT='($,'' Element '//CHAR1(IBEG:IEND)//
     '          ' [0,.. or previous defined values]: '',3(E12.5))'
              IF(IOTYPE.EQ.3) THEN
                DO i=1,3
                  ij=ILT(3,nr,nx)+3*(nj-1)+i
                  RDATA(i)=CE(ij,ne)
                ENDDO
              ENDIF
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RZERO,
     '          -RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                DO i=1,3
                  ij=ILT(3,nr,nx)+3*(nj-1)+i
                  IF(RDATA(i).EQ.0.0d0) THEN
                    CE(ij,ne)=CE(ij,ne-1)
                  ELSE
                    CE(ij,ne)=RDATA(i)
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        DO nj=1,2 !normalize
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(NW(ne,1).EQ.3) THEN
              ij=ILT(3,nr,nx)+3*(nj-1)
              SUM=CE(ij+1,ne)**2+CE(ij+2,ne)**2+CE(ij+3,ne)**2
              DO i=1,3
                CE(ij+i,ne)=CE(ij+i,ne)/DSQRT(SUM)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF(FILEIP) CALL CLOSEF(7,ERROR,*9999)
      CALL EXITS('IPMAT4')
      RETURN
 9999 CALL ERRORS('IPMAT4',ERROR)
      IF(FILEIP) CLOSE(UNIT=7)
      CALL EXITS('IPMAT4')
      RETURN 1
      END


