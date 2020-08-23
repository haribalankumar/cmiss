      SUBROUTINE IPCOUP(CALCU,NEELEM,NP_INTERFACE,nx,NXQ,
     '  TISSUE_REG,TREE_REG,XQ,ERROR,*)

C#### Subroutine: IPCOUP
C###  Description:
C###    IPCOUP defines appropriate coupling between nodes shared between
C###    more than one region.

C**** For coupled aerofoil flow and stress (KTYP90=5):
C****   NP_INTERFACE(0,i),i=1,3 are #node triplets on aerofoil
C****   NP_INTERFACE(n,i),i=1,3 are node #s on aerofoil
C****   i=1 is upper surface aerofoil flow field
C****   i=2 is lower surface aerofoil flow field
C****   i=3 is aerofoil stress field
C**** AJP 22-7-94
C**** For coupled prob involving a potential jump (Bdy layer:KTYP90=6):
C**** Along the bdy layer interface the nodal value at one node can be
C**** related linearly to the nodal value on the other side of the
C**** interface i.e. u2=a*u1+b where a and b are entered at each node.
C****   NP_INTERFACE(n,1),n=1,NP_INTERFACE(0,1) are the nodes at
C****   which the jump is applied.
C****   NP_INTERFACE(n,2),n=1,NP_INTERFACE(0,1) are the corresponding
C****   nodes on the other side of the jump
C****   i.e. u2 is at node np_interface(np,1)
C****   i.e. u1 is at node np_interface(np,2)
C****   Values a and b at node NP_INTERFACE(n,1) are stored in
C****   coupling_list(n,i), i=1,2 where n=1, np_interface(0,1)

      IMPLICIT NONE
      INCLUDE 'aero00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp90.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NP_INTERFACE(0:NPM,0:3),nx,
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     '  TISSUE_REG,TREE_REG
      REAL*8 XQ(NJM,NQM)
      CHARACTER ERROR*(*)
      LOGICAL CALCU
!     Local Variables
      INTEGER NCONMX
      PARAMETER (NCONMX=1000)
      INTEGER COUP_TYPE,I,IBEG,ICHAR,IEND,INFO,
     '  MAXOPT,N,NCI,ncp,ncp1,ncp2,ncp3,ne,njj,noelem,
     '  NOQUES,NO_COUP_TYPE,no_interface,no_tip,nq,TEMP_IDATA(0:IOIMX)
      REAL*8 DIST,CPDST(NCONMX)
      CHARACTER CHAR1*1,CHAR2*2,CHAR3*3,CHAR5*5,FSTRING*1024
      LOGICAL FILEIP,FOUND

      CALL ENTERS('IPCOUP',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(CALCU) THEN
        KTYP90=8
      ELSE
        FORMAT='('' Specify the problem setup  [1]: '''//
     '    '/''   (0) Reset                          '''//
     '    '/''   (1) Coupled saturated-unsaturated  '''//
     '    '/''   (2) Coupled Laplace                '''//
     '    '/''   (3) Coupled tree growth problem    '''//
     '    '/''   (4) Fluid interface stability      '''//
     '    '/''   (5) Coupled aerofoil flow & stress '''//
     '    '/''   (6) Boundary layer interface       '''//
     '    '/''   (7) Coupled pressure on surface    '''//
     '    '/''   (8) Coupled 1D and 2/3D grids      '''//
     '    '/''   (9) Coupled tube flow and mechanics'''//
     '    '/''  (10) Cellular processes             '''//
     '    '/$,''    '',I3)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP90
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,10,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP90=IDATA(1)
      ENDIF


      NOQUES=0
      IF(KTYP90.EQ.1) THEN !Coupled saturated-unsaturated
        MAXOPT=2
!       Define the optimisation parameters for the free surface
        FORMAT='('' Specify whether [1]:'''//
     '    '/''   (1) All interface nodes free           '''//
     '    '/''   (2) Interface node(s) fixed            '''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=1 !NEEDS FIXING
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,MAXOPT,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) IDATA(1)=1 !NEEDS FIXING

        CALL ASSERT(NCT(0,nx).NE.0,
     '    '>>NCT(0,nx) is not set. Check code',
     '    ERROR,*9999)

        CALL_COUP=.TRUE.
        CALL_SOLV=.FALSE.
      ELSE IF(KTYP90.EQ.2) THEN !Coupled Laplace

c cpb 9/3/95 Set up NCT for the coupled region
        NCT(0,nx)=2

C cpb 23/3/95
C        FORMAT='('' Specify the extent of derivative equation'//
C     '    ' generation at interface nodes [1]:'''//
C     '    '/''   (0) Reset                                    '''//
C     '    '/''   (1) Full generation                          '''//
C     '    '/''   (2) Generated at first interface region only '''//
C     '    '/$,''    '',I1)'
C        IF(IOTYPE.EQ.3) IDATA(1)=KTYP91
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,2,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) KTYP91=IDATA(1)

        CALL ASSERT(NCT(0,nx).NE.0,
     '    '>>NCT(0,nx) is not set. Check code',
     '    ERROR,*9999)
        CALL_COUP=.TRUE.
        CALL_SOLV=.FALSE.

      ELSE IF(KTYP90.EQ.3) THEN !Coupled tree growth problem
!       Define region 2 growing tip nodes in each element of region 1
        DO noelem=1,NEELEM(0,1)
          ne=NEELEM(noelem,1)
          CALL ASSERT(ne.LE.20,
     '      'Array dimension for NP_GROWING_TIP needs'
     '      //' increasing',ERROR,*9999)
          WRITE(CHAR5,'(I5)') ne
          CALL STRING_TRIM(CHAR5,IBEG,IEND)
          FORMAT=
     '      '($,'' Enter the number of growing tip nodes in element '
     '      //CHAR5(IBEG:IEND)//' [0]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=NP_GROWING_TIP(0,ne)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,20,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NP_GROWING_TIP(0,ne)=IDATA(1)

          IF(NP_GROWING_TIP(0,ne).GT.0) THEN
            WRITE(CHAR3,'(I3)') NP_GROWING_TIP(0,ne)
            CALL STRING_TRIM(CHAR3,IBEG,IEND)
            FORMAT='($,'' Enter the '//CHAR3(IBEG:IEND)
     '        //' growing tip nodes: '',20I5)'
            IF(IOTYPE.EQ.3) THEN
              DO no_tip=1,NP_GROWING_TIP(0,ne)
                IDATA(no_tip)=NP_GROWING_TIP(no_tip,ne)
              ENDDO
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NP_GROWING_TIP(0,ne),
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NPT(2),
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO no_tip=1,NP_GROWING_TIP(0,ne)
                NP_GROWING_TIP(no_tip,ne)=IDATA(no_tip)
              ENDDO
            ENDIF

            FORMAT='($,'' Enter the '//CHAR3(IBEG:IEND)
     '        //' growing tip flows: '',10E12.3)'
            IF(IOTYPE.EQ.3) THEN
              DO no_tip=1,NP_GROWING_TIP(0,ne)
                RDATA(no_tip)=FLOW_GROWING_TIP(no_tip,ne)
              ENDDO
            ENDIF
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NP_GROWING_TIP(0,ne),
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO no_tip=1,NP_GROWING_TIP(0,ne)
                FLOW_GROWING_TIP(no_tip,ne)=RDATA(no_tip)
              ENDDO
            ENDIF
          ENDIF
        ENDDO

        CALL ASSERT(NCT(0,nx).NE.0,
     '    '>>NCT(0,nx) is not set. Check code',
     '    ERROR,*9999)
        CALL_COUP=.TRUE.
        CALL_SOLV=.FALSE.

      ELSE IF(KTYP90.EQ.4) THEN !Fluid interface stability problem
        FORMAT='($,'' Enter the number of interface nodes [0]: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=NP_INTERFACE(0,1)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,20,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NP_INTERFACE(0,1)=IDATA(1)
        NP_INTERFACE(0,2)=NP_INTERFACE(0,1)
        NP_INTERFACE(0,3)=NP_INTERFACE(0,1)

        IF(NP_INTERFACE(0,1).GT.0) THEN
          DO I=1,3
            WRITE(CHAR1,'(I1)') I
            WRITE(CHAR3,'(I3)') NP_INTERFACE(0,I)
            CALL STRING_TRIM(CHAR3,IBEG,IEND)
            FORMAT='($,'' Enter the '//CHAR3(IBEG:IEND)
     '        //' interface nodes for region '//CHAR1//': '',20I5)'
            IF(IOTYPE.EQ.3) THEN
              DO no_interface=1,NP_INTERFACE(0,I)
                IDATA(no_interface)=NP_INTERFACE(no_interface,I)
              ENDDO
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NP_INTERFACE(0,I),
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO no_interface=1,NP_INTERFACE(0,I)
                NP_INTERFACE(no_interface,I)=IDATA(no_interface)
              ENDDO
            ENDIF
          ENDDO
        ENDIF

        CALL ASSERT(NCT(0,nx).NE.0,
     '    '>>NCT(0,nx) is not set. Check code',
     '    ERROR,*9999)
        CALL_COUP=.TRUE.
        CALL_SOLV=.FALSE.

      ELSE IF(KTYP90.EQ.5) THEN !Coupled aerofoil flow & stress
        FORMAT='($,'' Enter #node triplets on aerofoil [3]: '',I2)'
        IDEFLT(1)=3
        IF(IOTYPE.EQ.3) IDATA(1)=NP_INTERFACE(0,1)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,20,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NP_INTERFACE(0,1)=IDATA(1)
        NP_INTERFACE(0,2)=NP_INTERFACE(0,1)
        NP_INTERFACE(0,3)=NP_INTERFACE(0,1)

        IF(NP_INTERFACE(0,1).GT.0) THEN
          DO I=1,3
            WRITE(CHAR3,'(I3)') NP_INTERFACE(0,I)
            CALL STRING_TRIM(CHAR3,IBEG,IEND)
            IF(I.EQ.1) THEN
              FORMAT='($,'' Enter the '//CHAR3(IBEG:IEND)
     '          //' upper surface flow field nodes: '',20I5)'
            ELSE IF(I.EQ.2) THEN
              FORMAT='($,'' Enter the '//CHAR3(IBEG:IEND)
     '          //' lower surface flow field nodes: '',20I5)'
            ELSE IF(I.EQ.3) THEN
              FORMAT='($,'' Enter the '//CHAR3(IBEG:IEND)
     '          //' aerofoil stress field nodes: '',20I5)'
            ENDIF
            IF(IOTYPE.EQ.3) THEN
              DO no_interface=1,NP_INTERFACE(0,I)
                IDATA(no_interface)=NP_INTERFACE(no_interface,I)
              ENDDO
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NP_INTERFACE(0,I),
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO no_interface=1,NP_INTERFACE(0,I)
                NP_INTERFACE(no_interface,I)=IDATA(no_interface)
              ENDDO
            ENDIF
          ENDDO

          FORMAT='(/$,'' Enter basis# for aerofoil pressure [1]: '''//
     '      ',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=NB_AERO_PRESS
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NB_AERO_PRESS=IDATA(1)

        ENDIF

        CALL ASSERT(NCT(0,nx).NE.0,
     '    '>>NCT(0,nx) is not set. Check code',
     '    ERROR,*9999)
        CALL_COUP=.TRUE.
        CALL_SOLV=.FALSE.

      ELSE IF(KTYP90.EQ.6) THEN !Boundary layer interface
        FORMAT='($,'' Enter #nodes along interface [3]: '',I2)'
        IDEFLT(1)=3
        IF(IOTYPE.EQ.3) IDATA(1)=NP_INTERFACE(0,1)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,50,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NP_INTERFACE(0,1)=IDATA(1)
        DO n=1,NP_INTERFACE(0,1)
          WRITE(CHAR3,'(I3)') n
          CALL STRING_TRIM(CHAR3,IBEG,IEND)
          FORMAT='($,'' Enter the '//CHAR3(IBEG:IEND)
     '      //' node number at which the jump is applied (u2) '',I5)'
          IF(IOTYPE.EQ.3) THEN
            IDATA(1)=NP_INTERFACE(n,1)
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NP_INTERFACE(n,1)=IDATA(1)
          ENDIF
          FORMAT='($,'' Enter the corresponding node number '
     '      //' on the other side of the jump (u1) '',I5)'
          IF(IOTYPE.EQ.3) THEN
            IDATA(1)=NP_INTERFACE(N,2)
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NP_INTERFACE(n,2)=IDATA(1)
          ENDIF
          FORMAT='($,'' Enter the values of a and b at this jump '
     '          //' (u2=a*u1+b) [1.0,0.0]: '',2(E10.4,X))'
          RDEFLT(1)=1.0D0
          RDEFLT(2)=0.0D0
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=COUPLING_LIST(N,1)
            RDATA(2)=COUPLING_LIST(N,2)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            CALL ASSERT(DABS(RDATA(1)).GE.RDELTA,'>>a must be non zero',
     '        ERROR,*9999)
            COUPLING_LIST(N,1)=RDATA(1)
            COUPLING_LIST(N,2)=RDATA(2)
          ENDIF
        ENDDO

        CALL ASSERT(NCT(0,nx).NE.0,
     '    '>>NCT(0,nx) is not set. Check code',
     '    ERROR,*9999)
        CALL_COUP=.TRUE.
        CALL_SOLV=.FALSE.

      ELSE IF(KTYP90.EQ.7) THEN !Coupled pressure on surface
c       Set up NCT for the coupled region
        NCT(0,nx)=2

        CALL ASSERT(NCT(0,nx).NE.0,
     '    '>>NCT(0,nx) is not set. Check code',
     '    ERROR,*9999)
        CALL_COUP=.TRUE.
        CALL_SOLV=.FALSE.

      ELSE IF(KTYP90.EQ.8) THEN !1D to 2/3D grid coupling

! Check that there are grid points and the arrays are allocated
        CALL ASSERT(USE_GRID.EQ.1,'>>Must have USE_GRID=1 in ippara',
     '    ERROR,*9999)
        CALL ASSERT(NQT.GT.0,'>>No grid points defined.',ERROR,*9999)

C Initialising the default values
        DO ncp=1,NCONMX
          CPDST(ncp)=RMAX
        ENDDO
        DO ncp=0,NCONMX
          DO i=1,2
            CPLST(ncp,i)=0
          ENDDO
        ENDDO

! Automatic calculation of coupling
        IF(CALCU) THEN
          DO nq=NQR(1,TREE_REG),NQR(2,TREE_REG) !Loop over 1d points
            IF(NXQ(1,0,nq,1).EQ.0) THEN         !Is it a fibre end
              CPLST(0,1)=CPLST(0,1)+1           !Add it to the list of
              CPLST(CPLST(0,1),1)=nq            !points to couple
            ENDIF
          ENDDO !nq
! Couple grid points specified by the user only
        ELSE
          FORMAT='('' Enter the region number to couple to [1]:'',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=TISSUE_REG
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TISSUE_REG=IDATA(1)
          CALL ASSERT(TISSUE_REG.LE.NRM,
     '      '>>ERROR:Invalid region number entered',ERROR,*9999)

          FORMAT='('' Enter the number of 1-d grid points to be '
     '      //'coupled [1]:'',I2)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=CPLST(0,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) CPLST(0,1)=IDATA(1)
          CALL ASSERT(CPLST(0,1).LE.NCONMX,
     '      '>>Need to increase parameter NCONMX in IPCOUP',ERROR,*9999)

          DO ncp=1,CPLST(0,1)
            WRITE(CHAR3,'(I2)') ncp
            CALL STRING_TRIM(CHAR3,IBEG,IEND)
            FORMAT='('' Enter 1D grid point number '
     '        //CHAR3(IBEG:IEND)//' [1]:'',I5)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=CPLST(ncp,1)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        99999,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) CPLST(ncp,1)=IDATA(1)
          ENDDO !ncp
        ENDIF !calcu

        IF(CPLST(0,1).LE.0) THEN
          WRITE(OP_STRING,'(''WARNING: No 1d points found to couple'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE
! Find the closest 2 grid points in the 2/3d region
          DO nq=NQR(1,TISSUE_REG),NQR(2,TISSUE_REG)
            DO ncp=1,CPLST(0,1)
              DIST=0.0d0
              DO njj=1,NJT
                DIST=DIST+(XQ(njj,nq)-XQ(njj,CPLST(ncp,1)))**2
              ENDDO !njj
              DIST=DSQRT(DIST)
              IF(DIST.LT.CPDST(ncp)) THEN
                CPDST(ncp)=DIST
                CPLST(ncp,2)=nq
              ENDIF
            ENDDO !ncp
          ENDDO !nq

C MLB 13-April-1999 please leave
C
C! Put into NWQ array in 6th position
C          DO nq=1,CPLST(0,1)
C! Check that each grid point found a pair
C            CALL ASSERT(CPLST(nq,2).GT.0,
C     '        'No grid point found for coupling',ERROR,*9999)
C! Make the closest pair of grid points equivalent
C            NWQ(6,CPLST(nq,1),1)=CPLST(nq,2)  !add 1d - 2/3d link
C            NWQ(6,CPLST(nq,2),1)=CPLST(nq,1)  !add 2/3d - 1d link
C! Make the 1d grid point now be internal
C            NWQ(1,CPLST(nq,1),1)=0
C! Use this for last term in local quadratic element
C            NXQ(1,0,CPLST(nq,1),1)=1
C            NXQ(1,1,CPLST(nq,1),1)=CPLST(nq,2)
C! Write out coupling information
C            IF(DOP) THEN
C              WRITE(OP_STRING,'('' Grid point number '',I5,'' was '
C     '          //'coupled to point '',I5)') CPLST(nq,1),CPLST(nq,2)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              WRITE(OP_STRING,'('' Distance between points '',F9.6)')
C     '          CPDST(nq)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF
C          ENDDO
C        ENDIF

          IF(DOP) THEN
            DO nq=1,CPLST(0,1)
              WRITE(OP_STRING,'('' Grid point number '',I5,'' was '
     '          //'coupled to point '',I5)') CPLST(nq,1),CPLST(nq,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Distance between points '',F9.6)')
     '          CPDST(nq)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
        ENDIF

        CALL_COUP=.TRUE.
C        CALL_SOLV=.FALSE.

      ELSE IF(KTYP90.EQ.9) THEN
C Coupled Mechanics and flow in elastic tubes  problem

        FORMAT='($,'' Do you want start from the'//
     '    ' initial conditions [Y]? '',A)'
        ADEFLT(1)='Y'


        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

        IF(ADATA(1).EQ.'Y') THEN
          FLOW_MECH_INIT=.TRUE.
        ELSE
         FLOW_MECH_INIT=.FALSE.
        ENDIF

        FORMAT='($,'' Do you want to recalculate'//
     '    ' xi positions of grid points [Y]? '',A)'
        ADEFLT(1)='Y'

        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

        IF(ADATA(1).EQ.'Y') THEN
          CALC_GRID_XI=.TRUE.
        ELSE
          CALC_GRID_XI=.FALSE.
        ENDIF

        FORMAT='('' Enter method of stress calculation [0]:'''
     '    //'/''   (0) Interpolated from solution'''
     '    //'/''   (1) Interpolated from fiited values'''
     '    //'/$,''    '',I1)'

        IF(IOTYPE.EQ.3) IDATA(1)=0
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,900,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) STRESS_CALC=IDATA(1)

        FORMAT='($,'' Enter the number time steps per'//
     '    'load increment [10]: '',I3)'
        IF(IOTYPE.EQ.3) IDATA(1)=FLOW_STEPS
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) FLOW_STEPS=IDATA(1)
        CALL_COUP=.TRUE.

      ELSEIF(KTYP90.EQ.10)THEN ! Coupling models

        CALL COUP_LIST(FSTRING,NO_COUP_TYPE,0,ERROR,*9999)

        DO COUP_TYPE=1,NO_COUP_TYPE
          CALL COUP_LIST(FSTRING,MAXOPT,COUP_TYPE,ERROR,*9999)
          IF(IOTYPE.EQ.3) THEN
            IF(COUP_CL(0,1,0).EQ.0) THEN
              IDATA(0)=1
              IDATA(1)=0
            ELSE
              IDATA(0)=0
              IDATA(1)=0
              DO ncp1=1,COUP_CL(0,1,0)
                IF(COUP_CL(ncp1,1,0).EQ.COUP_TYPE)THEN
                  FOUND=.FALSE.
                  DO ncp2=1,IDATA(0)
                    IF(COUP_CL(ncp1,2,0).EQ.IDATA(ncp2))THEN
                      FOUND=.TRUE.
                    ENDIF
                  ENDDO
                  IF(.NOT.FOUND)THEN
                    IDATA(0)=IDATA(0)+1
                    IDATA(IDATA(0))=COUP_CL(ncp1,2,0)
                  ENDIF
                ENDIF
              ENDDO
              IF(IDATA(0).EQ.0)THEN
                IDATA(0)=1
                IDATA(1)=0
              ENDIF
            ENDIF
          ENDIF
          CALL STRING_TRIM(FSTRING,IBEG,IEND)
          FORMAT=FSTRING(IBEG:IEND)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      IDATA(0),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,
     '      MAXOPT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO ncp1=0,IDATA(0)
              TEMP_IDATA(ncp1)=IDATA(ncp1)
            ENDDO
            DO ncp1=1,TEMP_IDATA(0)
              FOUND=.FALSE.
              DO ncp2=1,COUP_CL(0,1,0)
                IF((COUP_CL(ncp2,1,0).EQ.COUP_TYPE).AND.
     '            (COUP_CL(ncp2,2,0).EQ.TEMP_IDATA(ncp1))) THEN
                  FOUND=.TRUE.
                ENDIF
              ENDDO
              IF(.NOT.FOUND.AND.TEMP_IDATA(1).NE.0) THEN
                WRITE(CHAR2,'(I2)') TEMP_IDATA(ncp1)
                CALL STRING_TRIM(CHAR2,IBEG,IEND)
                FORMAT='($,'' Specify the number of instances for '
     '            //'coupling model '//CHAR2(IBEG:IEND)//' [1]: '',I2)'
                IDEFLT(1)=1
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) NCI=IDATA(1)

                DO ncp3=1,NCI
                  KTYP3C(nx)=1
                  COUP_CL(0,1,0)=COUP_CL(0,1,0)+1
                  COUP_CL(COUP_CL(0,1,0),1,0)=COUP_TYPE
                  COUP_CL(COUP_CL(0,1,0),2,0)=TEMP_IDATA(ncp1)
                  COUP_CL(COUP_CL(0,1,0),3,0)=ncp3
                ENDDO

              ENDIF
            ENDDO
          ENDIF
          IF(IOTYPE.EQ.3) THEN
            DO ncp1=1,MAXOPT
              NCI=0
              DO ncp2=1,COUP_CL(0,1,0)
                IF((COUP_CL(ncp2,1,0).EQ.COUP_TYPE).AND.
     '            (COUP_CL(ncp2,2,0).EQ.ncp1)) THEN
                  IF(COUP_CL(ncp2,3,0).GT.NCI)THEN
                    NCI=COUP_CL(ncp2,3,0)
                  ENDIF
                ENDIF
              ENDDO
              IF(NCI.GT.0) THEN
                WRITE(CHAR2,'(I2)') ncp1
                CALL STRING_TRIM(CHAR2,IBEG,IEND)
                FORMAT='($,'' Specify the number of instances for '
     '          //'coupling model '//CHAR2(IBEG:IEND)//' [1]: '',I2)'
                IDATA(0)=1
                IDATA(1)=NCI
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,1,99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)
              ENDIF
            ENDDO
          ENDIF
        ENDDO

      ENDIF !End of KTYP90 decision

      CALL EXITS('IPCOUP')
      RETURN
 9999 CALL ERRORS('IPCOUP',ERROR)
      CALL EXITS('IPCOUP')
      RETURN 1
      END

