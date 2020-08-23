      SUBROUTINE IPANA9(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,INP,
     '  NBJ,NDIPOLES,NEELEM,NENP,NKJE,NKH,NPF,NP_INTERFACE,NPNE,NPNODE,
     '  nr,NRE,NVJE,nx,NYNP,CE,DIPOLE_CEN,DIPOLE_DIR,SE,XA,XE,XP,YP,
     '  ERROR,*)

C#### Subroutine: IPANA9
C###  Description:
C###    IPANA9 inputs analytic functions for FE90 equations for
C###    region nr (but sometimes other regions as a side-effect) and
C###    problem nx.

      IMPLICIT NONE
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh03.cmn'
      INCLUDE 'tol00.cmn'


!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NDIPOLES(NRM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM),
     '  NKJE(NKM,NNM,NJM,NEM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER ICHAR,INFO,MAX_SERIES,
     '  nc,ne,nh,nj,nk,nonode,nonp,NOQUES,np,nrr,
     '  NUM_NODES,nv,ny
      REAL*8 ALPHA,BETA,C_0,COORD(3),DATAN_MOD,E_1,E_11,
     '  EXACT,F_1,F_11,FACTOR,G_0,G_1,G_11,GE,GI,GRADF(3),H_1,
     '  HESSIAN(3,3),K,OFFSET,PHI,R,RAD,RAD_1,RAD_2,RADIUS,RDOTN,
     '  RHO_X,RHO_Y,RHO_Z,SIGMAO,SIGMAI,SIGMAE,SIGN,STRENGTH,
     '  SUM_ANAL,SUM_SOL,TANGENT(3,3),THETA,XN_LOCAL(3)
      CHARACTER CHAR5*5
      LOGICAL AVERAGE,FILEIP,FOUND

!     Functions
      REAL*8 ANALY_PHI_M

      CALL ENTERS('IPANA9',*9999)

      nv=1 ! temporary cpb 7/12/94

      NOQUES=0
      FILEIP=.FALSE.
      ICHAR=999

      IF(ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4) THEN !static or quasi-static
        IF(ITYP2(nr,nx).EQ.3) THEN
          IF(NJT.EQ.2) THEN
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=ANAL_CHOICE(nr)
            ENDIF
            FORMAT='('' Specify the analytic solution [1]:'''//
     '        '/''   (1) K(x-y)'''//
     '        '/''   (2) K(x^2-y^2)'''//
     '        '/''   (3) K(x^2+2xy-y^2)'''//
     '        '/''   (4) K(a.r^n.cos(n.t)+b.r^m.sin(m.t)+c.r.cos(t)+'
     '        //'d.r.sin(t)+e)'''//
     '        '/''   (5) Centre dipole (single circle)'''//
     '        '/''   (6) Centre dipole (multiple circles)'''//
     '        '/''  *(7) Eccentric dipole (single circle)'''//
     '        '/''  *(8) Eccentric dipole (multiple circles)'''//
     '        '/$,''   '',I1)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,7,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_CHOICE(nr)=IDATA(1)
            IF(ANAL_CHOICE(nr).GE.1.AND.ANAL_CHOICE(nr).LE.3) THEN
              FORMAT='(/$,'' Enter the value of K [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANAL_K
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_K=RDATA(1)
            ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN
C              CALL ASSERT(CALL_MESH.AND.JTYP14.EQ.2,
C     '          '>>Define a circular mesh first',ERROR,*9999)
              FORMAT='(/$,'' Enter the value of K [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANAL_K
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_K=RDATA(1)
              FORMAT='(/$,'' Enter the value of n [2]: '',I1)'
              IDEFLT(1)=2
              IF(IOTYPE.EQ.3) IDATA(1)=ANAL_N
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          9,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_N=IDATA(1)
              FORMAT='(/$,'' Enter the value of m [2]: '',I1)'
              IDEFLT(1)=2
              IF(IOTYPE.EQ.3) IDATA(1)=ANAL_M
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          9,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_M=IDATA(1)
              FORMAT='(/$,'' Enter the value of a [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANAL_A
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_A=RDATA(1)
              FORMAT='(/$,'' Enter the value of b [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANAL_B
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_B=RDATA(1)
              FORMAT='(/$,'' Enter the value of c [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANAL_C
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_C=RDATA(1)
              FORMAT='(/$,'' Enter the value of d [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANAL_D
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_D=RDATA(1)
              FORMAT='(/$,'' Enter the value of e [0.0]: '',D12.4)'
              RDEFLT(1)=0.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANAL_E
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_E=RDATA(1)
            ELSE IF(ANAL_CHOICE(nr).GE.5.AND.ANAL_CHOICE(nr).LE.6) THEN
C             Dipole solutions
              CALL ASSERT(CALL_MATE,'>>Define materials first',
     '          ERROR,*9999)
              IF(nr.EQ.1) THEN
                IF(NDIPOLES(1).EQ.0) THEN
                  ERROR='>>Define a dipole in region 1 first'
                  GOTO 9999
                ELSE IF(NDIPOLES(1).GT.1) THEN
                  ERROR='>>Must only have one dipole in region 1'
                  GOTO 9999
                ELSE
                  IF(DIPOLE_CEN(1,0,1,1).NE.0.0d0.OR.
     '              DIPOLE_CEN(2,0,1,1).NE.0.0d0) THEN
                    ERROR='>>The dipole must be centred at the origin'
                    GOTO 9999
                  ENDIF
                ENDIF
                ANAL_A=DIPOLE_DIR(1,0,1,1)
                ANAL_B=DIPOLE_DIR(2,0,1,1)
                IF(ITYP3(nr,nx).EQ.2) THEN !generalised Laplace
                  SIGMA(1)=CE(1,NEELEM(1,1))
                ELSE
                  SIGMA(1)=1.0d0
                ENDIF
                IF(ANAL_CHOICE(nr).EQ.6) THEN
                  CALL ASSERT(CALL_MESH.AND.JTYP14.EQ.2,
     '              '>>Define a circular mesh first',ERROR,*9999)
                  CALL ASSERT(NSPHERES.GT.1,
     '              '>>Must have more than one circle',ERROR,*9999)
                  DO nrr=2,NRT
                    IF(ITYP3(nrr,nx).EQ.2) THEN !generalised Laplace
                      SIGMA(nrr)=CE(1,NEELEM(1,nrr))
                    ELSE
                      SIGMA(nrr)=1.0d0
                    ENDIF
                  ENDDO !nrr
C***              Find the coefficients for the analytic solution
                  CALL DIPOLE_SOLVE(ANAL_CHOICE(nr),
     '              DIPOLE_CEN(1,0,1,1),DIPOLE_DIR(1,0,1,1),
     '              ERROR,*9999)
                ENDIF
                IF(ANAL_CHOICE(nr).EQ.5) THEN
                  IDEFLT(1)=NPNODE(1,1)
                ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN
C***              Find the node number to fix
                  FOUND=.FALSE.
                  nonode=1
                  DO WHILE(.NOT.FOUND.AND.nonode.LE.NPNODE(0,NRT))
                    np=NPNODE(nonode,NRT)
                    R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2)
                    IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL) THEN
                      FOUND=.TRUE.
                    ELSE
                      nonode=nonode+1
                    ENDIF
                  ENDDO !nonode (np)
                  IF(FOUND) THEN
                    IDEFLT(1)=np
                  ELSE
                    ERROR='>>Could not find first outer circle node'
                    GOTO 9999
                  ENDIF
                ENDIF
                WRITE(CHAR5,'(I5)') IDEFLT(1)
                FORMAT='(/$,'' Enter the outer circle node number '
     '            //'to fix the potential at ['//CHAR5//']: '',I5)'
                IF(IOTYPE.EQ.3) IDATA(1)=ANAL_FIXEDNODE
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,1,NPT(0),LDATA,LDEFLT,RDATA,RDEFLT,
     '            RMIN,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  IF(ANAL_CHOICE(nr).EQ.5) THEN
                    nrr=1
                  ELSE
                    nrr=NRT
                  ENDIF
                  FOUND=.FALSE.
                  nonode=1
                  DO WHILE(.NOT.FOUND.AND.nonode.LE.NPNODE(0,nrr))
                    IF(IDATA(1).EQ.NPNODE(nonode,nrr)) THEN
                      FOUND=.TRUE.
                    ELSE
                      nonode=nonode+1
                    ENDIF
                  ENDDO
                  CALL ASSERT(FOUND,'>>Node is not in the outer region',
     '              ERROR,*9999)
                  IF(ANAL_CHOICE(nr).EQ.6) THEN
                    R=DSQRT(XP(1,nv,1,IDATA(1))**2+
     '                XP(1,nv,2,IDATA(1))**2)
                    IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                MESH3_RAD(NSPHERES).GT.SPHERE_RAD_TOL) THEN
                      ERROR='>>Node is not on the outer circle'
                      GOTO 9999
                    ENDIF
                  ENDIF
                  ANAL_FIXEDNODE=IDATA(1)
                ENDIF
                IDEFLT(1)=2 !fix value and deriv
                WRITE(CHAR5,'(I1)') IDEFLT(1)
                FORMAT='('' Specify whether ['//CHAR5(1:1)//']:'''//
     '            '/''   (1) Fix value, solve for normal derivative'''//
     '            '/''   (2) Fix value and normal derivative       '''//
     '            '/$,''    '',I1)'
                IF(IOTYPE.EQ.3) IDATA(1)=POT_BC_TYPE
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,1,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) POT_BC_TYPE=IDATA(1)
              ELSE
                CALL ASSERT(NDIPOLES(nr).EQ.0,
     '            '>>Must only have a dipole in the first region',
     '            ERROR,*9999)
              ENDIF
            ENDIF
          ELSE IF(NJT.EQ.3) THEN
            IF(IOTYPE.EQ.3) IDATA(1)=ANAL_CHOICE(nr)
            FORMAT='('' Specify the analytic solution [1]:'''//
     '        '/''   (1) K(x-y)'''//
     '        '/''   (2) K(z)'''//
     '        '/''   (3) K(x^2-y^2)'''//
     '        '/''   (4) K(x^2-z^2)'''//
     '        '/''   (5) K(y^2-z^2)'''//
     '        '/''   (6) K(x^2+y^2-2z^2)'''//
     '        '/''   (7) K(x^2-2y^2+z^2)'''//
     '        '/''   (8) K(-2x^2+y^2+z^2)'''//
     '        '/''   (9) K(x^2+2xy-y^2)'''//
     '        '/''  (10) K(x^2+2xz-z^2)'''//
     '        '/''  (11) K(y^2+2yz-z^2)'''//
     '        '/''  (12) Centre dipole (single sphere)'''//
     '        '/''  (13) Centre dipole (multiple spheres)'''//
     '        '/''  (14) Eccentric dipole (single sphere)'''//
     '        '/''  (15) Eccentric dipole (multiple spheres)'''//
     '        '/$,''   '',I2)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,15,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_CHOICE(nr)=IDATA(1)
            IF(ANAL_CHOICE(nr).GE.1.AND.ANAL_CHOICE(nr).LE.11) THEN
              FORMAT='(/$,'' Enter the value of K [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANAL_K
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_K=RDATA(1)
            ELSE IF(ANAL_CHOICE(nr).GE.12.AND.
     '          ANAL_CHOICE(nr).LE.15) THEN !Centre dipole solns
C GMH 21/6/96 Move to adaptive determination of series limits
              IF(ANAL_CHOICE(nr).EQ.15) THEN
                FORMAT='(/$,'' Enter the # of series [10]: '',I3)'
                IDEFLT(1)=10
                MAX_SERIES=NSPHERECOEFFMX/2-1
                IF(IOTYPE.EQ.3) IDATA(1)=ANAL_N
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,1,MAX_SERIES,LDATA,LDEFLT,RDATA,
     '            RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_N=IDATA(1)
              ENDIF
              CALL ASSERT(CALL_MATE,'>> Define materials first',
     '          ERROR,*9999)
              IF(nr.EQ.1) THEN
                IF(NDIPOLES(1).EQ.0) THEN
                  ERROR='>>Define a dipole in region 1 first'
                  GOTO 9999
                ELSE IF(NDIPOLES(1).GT.1) THEN
                  ERROR='>>Must only have one dipole in region 1'
                  GOTO 9999
                ELSE
                  IF(ANAL_CHOICE(nr).LE.13) THEN
                    IF(DIPOLE_CEN(1,0,1,1).NE.0.0d0.OR.
     '                DIPOLE_CEN(2,0,1,1).NE.0.0d0.OR.
     '                DIPOLE_CEN(3,0,1,1).NE.0.0d0) THEN
                      ERROR='>>The dipole must be centred at the origin'
                      GOTO 9999
                    ENDIF
                  ENDIF
                ENDIF
                ANAL_A=DIPOLE_DIR(1,0,1,1)
                ANAL_B=DIPOLE_DIR(2,0,1,1)
                ANAL_C=DIPOLE_DIR(3,0,1,1)
                IF(ITYP3(nr,nx).EQ.2) THEN !generalised Laplace
                  SIGMA(1)=CE(1,NEELEM(1,1))
                ELSE
                  SIGMA(1)=1.0d0
                ENDIF
                IF(ANAL_CHOICE(nr).EQ.13.OR.ANAL_CHOICE(nr).EQ.15) THEN
                  CALL ASSERT(CALL_MESH.AND.JTYP14.EQ.2,
     '              '>>Define a spherical mesh first',ERROR,*9999)
                  IF(ANAL_CHOICE(nr).NE.15) THEN !we can handle 1 sphere
                    CALL ASSERT(NSPHERES.GT.1,
     '                '>>Must have more than one sphere',ERROR,*9999)
                  ENDIF
                  DO nrr=2,NRT
                    IF(ITYP3(nrr,nx).EQ.2) THEN !generalised Laplace
                      SIGMA(nrr)=CE(1,NEELEM(1,nrr))
                    ELSE
                      SIGMA(nrr)=1.0d0
                    ENDIF
                  ENDDO !nrr
C***              Find the coefficients for the analytic solution
                  CALL DIPOLE_SOLVE(ANAL_CHOICE(nr),
     '              DIPOLE_CEN(1,0,1,1),DIPOLE_DIR(1,0,1,1),
     '              ERROR,*9999)
                ENDIF
                IF(ANAL_CHOICE(nr).EQ.12.OR.ANAL_CHOICE(nr).EQ.14) THEN
                  IDEFLT(1)=NPNODE(1,1)
                ELSE IF(ANAL_CHOICE(nr).GE.13.OR.
     '              ANAL_CHOICE(nr).EQ.14) THEN
C***              Find the node number to fix
                  FOUND=.FALSE.
                  nonode=1
                  DO WHILE(.NOT.FOUND.AND.nonode.LE.NPNODE(0,NRT))
                    np=NPNODE(nonode,NRT)
                    R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
     '                XP(1,nv,3,np)**2)
                    IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL) THEN
                      FOUND=.TRUE.
                    ELSE
                      nonode=nonode+1
                    ENDIF
                  ENDDO !nonode (np)
                  IF(FOUND) THEN
                    IDEFLT(1)=np
                  ELSE
                    ERROR='>>Could not find first outer sphere node'
                    GOTO 9999
                  ENDIF
                ENDIF
                WRITE(CHAR5,'(I5)') IDEFLT(1)
                FORMAT='(/$,'' Enter the outer sphere node number '
     '            //'to fix the potential at ['//CHAR5//']: '',I5)'
                IF(IOTYPE.EQ.3) IDATA(1)=ANAL_FIXEDNODE
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,1,NPT(0),LDATA,LDEFLT,RDATA,RDEFLT,
     '            RMIN,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  IF(ANAL_CHOICE(nr).EQ.12.OR.
     '              ANAL_CHOICE(nr).EQ.14) THEN
                    nrr=1
                  ELSE
                    nrr=NRT
                  ENDIF
                  FOUND=.FALSE.
                  nonode=1
                  DO WHILE(.NOT.FOUND.AND.nonode.LE.NPNODE(0,nrr))
                    IF(IDATA(1).EQ.NPNODE(nonode,nrr)) THEN
                      FOUND=.TRUE.
                    ELSE
                      nonode=nonode+1
                    ENDIF
                  ENDDO
                  CALL ASSERT(FOUND,'>>Node is not in the outer region',
     '              ERROR,*9999)
                  IF(ANAL_CHOICE(nr).EQ.6) THEN
                    R=DSQRT(XP(1,nv,1,IDATA(1))**2+
     '                XP(1,nv,2,IDATA(1))**2+XP(1,nv,3,IDATA(1))**2)
                    IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                MESH3_RAD(NSPHERES).GT.SPHERE_RAD_TOL) THEN
                      ERROR='>>Node is not on the outer sphere'
                      GOTO 9999
                    ENDIF
                  ENDIF
                  ANAL_FIXEDNODE=IDATA(1)
                ENDIF
                IDEFLT(1)=2 !fix value and deriv
                WRITE(CHAR5,'(I1)') IDEFLT(1)
                FORMAT='('' Specify whether ['//CHAR5(1:1)//']:'''//
     '            '/''   (1) Fix value, solve for derivative'''//
     '            '/''   (2) Fix value and derivative       '''//
     '            '/$,''    '',I1)'
                IF(IOTYPE.EQ.3) IDATA(1)=POT_BC_TYPE
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,1,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) POT_BC_TYPE=IDATA(1)
              ELSE
                CALL ASSERT(NDIPOLES(nr).EQ.0,
     '            '>>Must only have a dipole in the first region',
     '            ERROR,*9999)
              ENDIF
            ENDIF !anal_choice
          ENDIF
C GMH 7/11/95 Allowing analytic solution to be displaced by
C             a fixed amount - ie minimising average error
C             over the domain.
          FORMAT='($,'' Do you want to equate the average potential'
     '      //' to the average solution [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='Y'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'Y') THEN
            AVERAGE=.TRUE.
            WRITE(OP_STRING,'('' >>Remember solution must be '
     '        //'in YP'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE
            AVERAGE=.FALSE.
          ENDIF
C***  Calculate strength of dipole solns (if any have been defined)
          IF(NDIPOLES(nr).GT.0) THEN
            STRENGTH=0.0d0
            DO nj=1,NJT
              STRENGTH=STRENGTH+DIPOLE_DIR(nj,0,1,nr)**2
            ENDDO ! nj
            STRENGTH=DSQRT(STRENGTH)
          ELSE
            STRENGTH=1.0d0
          ENDIF

          nc=1 ! variables first

          NUM_NODES=0
          SUM_ANAL=0.0D0
          SUM_SOL=0.0D0
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            DO nk=1,MAX(NKH(NH_LOC(1,nx),np,nc,nr)-KTYP93(nc,nr),1)
              ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nr) !global var number
              IF(NJT.EQ.2) THEN
                IF(nk.EQ.1) THEN
                  IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                    EXACT=ANAL_K*(XP(1,nv,1,np)-XP(1,nv,2,np))
                  ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(x^2-y^2)
                    EXACT=ANAL_K*(XP(1,nv,1,np)**2-XP(1,nv,2,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2+2xy-y^2)
                    EXACT=ANAL_K*(XP(1,nv,1,np)**2+2.0d0*XP(1,nv,1,np)*
     '                XP(1,nv,2,np)-XP(1,nv,2,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(r^n.cos(n.t) ...
                    RADIUS=DSQRT(XP(1,nv,1,np)*XP(1,nv,1,np)+
     '                XP(1,nv,2,np)*XP(1,nv,2,np))
                    THETA=DATAN_MOD(XP(1,nv,1,np),XP(1,nv,2,np))
                    EXACT=ANAL_K*(ANAL_A*RADIUS**ANAL_N*
     '                DCOS(ANAL_N*THETA)+ANAL_B*RADIUS**ANAL_M*
     '                DSIN(ANAL_M*THETA)+ANAL_C*RADIUS*DCOS(THETA)+
     '                ANAL_D*RADIUS*DSIN(THETA)+ANAL_E)
                  ELSE IF(ANAL_CHOICE(nr).GE.5.AND.
     '                ANAL_CHOICE(nr).LE.6) THEN !dipole solutions
                    ne=0
                    DO nonp=1,NENP(np,0)
                      IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                    ENDDO !nonp
                    CALL ASSERT(ne.NE.0,
     '                '>>Could not find an element in the region',
     '                ERROR,*9999)
                    CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                DIPOLE_DIR_NTIME,1,NDIPOLES,np,
     '                NP_INTERFACE,nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                XP(1,1,1,np),
     '                EXACT,XP(1,1,1,np),XP(1,1,2,np),XP(1,1,1,np),
     '                ERROR,*9999)
                  ELSE
                    ERROR='>>Analytic solution type not implemented'
                    GOTO 9999
                  ENDIF
                ELSE
                  IF(ANAL_CHOICE(nr).EQ.4) THEN
C                   Tangent vector for a circle in polar coordinates
                    RADIUS=DSQRT(XP(1,nv,1,np)*XP(1,nv,1,np)+
     '                XP(1,nv,2,np)*XP(1,nv,2,np))
                    THETA=DATAN_MOD(XP(1,nv,1,np),XP(1,nv,2,np))
                    IF(CALL_MESH.AND.JTYP14.EQ.2.AND.NJT.EQ.2) THEN
                      IF(DABS(RADIUS-MESH3_RAD(1))/
     '                  MESH3_RAD(1).LE.SPHERE_RAD_TOL.AND.
     '                  NSPHERES.GT.1.AND.MESH3_BEMTYPE.EQ.2) THEN
C                       First circle in an annular problem therefore
C                       s is anticlockwise
                        TANGENT(1,1)=0.0d0
                        TANGENT(2,1)=1.0d0
                      ELSE
C                       s is clockwise
                        TANGENT(1,1)=0.0d0
                        TANGENT(2,1)=-1.0d0
                      ENDIF
                    ENDIF
                  ELSE
                    CALL GET_TNVECTOR(IBT,IDO,INP,NBJ,NENP,NKJE,np,NPF,
     '                NP_INTERFACE,NPNE,nr,NRE,NVJE,nx,XN_LOCAL,SE,
     '                TANGENT,XA,XE,XP,ERROR,*9999)
                  ENDIF
                  IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                    GRADF(1)=ANAL_K
                    GRADF(2)=-ANAL_K
                    EXACT=GRADF(1)*TANGENT(1,1)+GRADF(2)*TANGENT(2,1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(x^2-y^2)
                    GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(2)=-2.0d0*ANAL_K*XP(1,nv,2,np)
                    EXACT=GRADF(1)*TANGENT(1,1)+GRADF(2)*TANGENT(2,1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2+2xy-y^2)
                    GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,1,np)+XP(1,nv,2,np))
                    GRADF(2)=2.0d0*ANAL_K*(XP(1,nv,1,np)-XP(1,nv,2,np))
                    EXACT=GRADF(1)*TANGENT(1,1)+GRADF(2)*TANGENT(2,1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(r^n.cos(n.t) ...
                    GRADF(1)=ANAL_K*(ANAL_A*ANAL_N*RADIUS**(ANAL_N-1)*
     '                DCOS(ANAL_N*THETA)+ANAL_B*ANAL_M*
     '                RADIUS**(ANAL_M-1)*DSIN(ANAL_M*THETA)+ANAL_C*
     '                DCOS(THETA)+ANAL_D*DSIN(THETA))
                    GRADF(2)=ANAL_K*(-ANAL_A*ANAL_N*RADIUS**(ANAL_N-1)*
     '                DSIN(ANAL_N*THETA)+ANAL_B*ANAL_M*
     '                RADIUS**(ANAL_M-1)*DCOS(ANAL_M*THETA)-ANAL_C*
     '                DSIN(THETA)+ANAL_D*DCOS(THETA))
                    EXACT=GRADF(1)*TANGENT(1,1)+GRADF(2)*TANGENT(2,1)
                  ELSE IF(ANAL_CHOICE(nr).GE.5.AND.
     '                ANAL_CHOICE(nr).LE.6) THEN !Dipole solutions
                    ne=0
                    DO nonp=1,NENP(np,0)
                      IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                    ENDDO !nonp
                    CALL ASSERT(ne.NE.0,
     '                '>>Could not find an element in the region',
     '                ERROR,*9999)
                    CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                DIPOLE_DIR_NTIME,2,NDIPOLES,np,
     '                NP_INTERFACE,nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                XP(1,1,1,np),
     '                EXACT,XP(1,1,1,np),XP(1,1,2,np),XP(1,1,1,np),
     '                ERROR,*9999)
                  ELSE
                    ERROR='>>Analytic solution type not implemented'
                    GOTO 9999
                  ENDIF
                ENDIF ! nk
              ELSE IF(NJT.EQ.3) THEN
                IF(nk.EQ.1) THEN
                  IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                    EXACT=ANAL_K*(XP(1,nv,1,np)-XP(1,nv,2,np))
                  ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(z)
                    EXACT=ANAL_K*XP(1,nv,3,np)
                  ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2-y^2)
                    EXACT=ANAL_K*(XP(1,nv,1,np)**2-XP(1,nv,2,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(x^2-z^2)
                    EXACT=ANAL_K*(XP(1,nv,1,np)**2-XP(1,nv,3,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.5) THEN !K(y^2-z^2)
                    EXACT=ANAL_K*(XP(1,nv,2,np)**2-XP(1,nv,3,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !K(x^2+y^2-2z^2)
                    EXACT=ANAL_K*(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2-
     '                2.0d0*XP(1,nv,3,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.7) THEN !K(x^2-2y^2+z^2)
                    EXACT=ANAL_K*(XP(1,nv,1,np)**2-2.0d0*
     '                XP(1,nv,2,np)**2+XP(1,nv,3,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.8) THEN !K(-2x^2+y^2+z^2)
                    EXACT=ANAL_K*(-2.0d0*XP(1,nv,1,np)**2+
     '                XP(1,nv,2,np)**2+XP(1,nv,3,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN !K(x^2+2xy-y^2)
                    EXACT=ANAL_K*(XP(1,nv,1,np)**2+2.0d0*XP(1,nv,1,np)*
     '                XP(1,nv,2,np)-XP(1,nv,2,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN !K(x^2+2xz-z^2)
                    EXACT=ANAL_K*(XP(1,nv,1,np)**2+2.0d0*XP(1,nv,1,np)*
     '                XP(1,nv,3,np)-XP(1,nv,3,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN !K(y^2+2yz-z^2)
                    EXACT=ANAL_K*(XP(1,nv,2,np)**2+2.0d0*XP(1,nv,2,np)*
     '                XP(1,nv,3,np)-XP(1,nv,3,np)**2)
                  ELSE IF(ANAL_CHOICE(nr).GE.12.AND.
     '                ANAL_CHOICE(nr).LE.15) THEN
C                   Dipole solutions
                    ne=0
                    DO nonp=1,NENP(np,0)
                      IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                    ENDDO !nonp
                    CALL ASSERT(ne.NE.0,
     '                '>>Could not find an element in the region',
     '                ERROR,*9999)
                    CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                DIPOLE_DIR_NTIME,1,NDIPOLES,np,
     '                NP_INTERFACE,nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                XP(1,1,1,np),EXACT,XP(1,1,1,np),
     '                XP(1,1,2,np),XP(1,1,3,np),ERROR,*9999)
                  ELSE
                    ERROR='>>Analytic solution type not implemented'
                    GOTO 9999
                  ENDIF
                ELSE
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    TANGENT(nj,nk-1)=XP(nk,nv,nj,np)
                  ENDDO
                  IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                    GRADF(1)=ANAL_K
                    GRADF(2)=-ANAL_K
                    EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                TANGENT(2,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(z)
                    GRADF(3)=ANAL_K
                    EXACT=GRADF(3)*TANGENT(3,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2-y^2)
                    GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(2)=-2.0d0*ANAL_K*XP(1,nv,2,np)
                    EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                TANGENT(2,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(x^2-z^2)
                    GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(3)=-2.0d0*ANAL_K*XP(1,nv,3,np)
                    EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(3)*
     '                TANGENT(3,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.5) THEN !K(y^2-z^2)
                    GRADF(2)=2.0d0*ANAL_K*XP(1,nv,2,np)
                    GRADF(3)=-2.0d0*ANAL_K*XP(1,nv,3,np)
                    EXACT=GRADF(2)*TANGENT(2,nk-1)+GRADF(3)*
     '                TANGENT(3,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !K(x^2+y^2-2z^2)
                    GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(2)=2.0d0*ANAL_K*XP(1,nv,2,np)
                    GRADF(3)=-4.0d0*ANAL_K*XP(1,nv,3,np)
                    EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                TANGENT(2,nk-1)+GRADF(3)*TANGENT(3,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.7) THEN !K(x^2-2y^2+z^2)
                    GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(2)=-4.0d0*ANAL_K*XP(1,nv,2,np)
                    GRADF(3)=2.0d0*ANAL_K*XP(1,nv,3,np)
                    EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                TANGENT(2,nk-1)+GRADF(3)*TANGENT(3,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.8) THEN !K(-2x^2+y^2+z^2)
                    GRADF(1)=-4.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(2)=2.0d0*ANAL_K*XP(1,nv,2,np)
                    GRADF(3)=2.0d0*ANAL_K*XP(1,nv,3,np)
                    EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                TANGENT(2,nk-1)+GRADF(3)*TANGENT(3,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN !K(x^2+2xy-y^2)
                    GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,1,np)+XP(1,nv,2,np))
                    GRADF(2)=2.0d0*ANAL_K*(XP(1,nv,1,np)-XP(1,nv,2,np))
                    EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                TANGENT(2,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN !K(x^2+2xz-z^2)
                    GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,1,np)+XP(1,nv,3,np))
                    GRADF(3)=2.0d0*ANAL_K*(XP(1,nv,1,np)-XP(1,nv,3,np))
                    EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(3)*
     '                TANGENT(3,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN !K(y^2+2yz-z^2)
                    GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,2,np)+XP(1,nv,3,np))
                    GRADF(2)=2.0d0*ANAL_K*(XP(1,nv,2,np)-XP(1,nv,3,np))
                    EXACT=GRADF(2)*TANGENT(2,nk-1)+GRADF(3)*
     '                TANGENT(3,nk-1)
                  ELSE IF(ANAL_CHOICE(nr).GE.12.AND
     '                .ANAL_CHOICE(nr).LE.15) THEN
C                   Dipole solutions
                    ne=0
                    DO nonp=1,NENP(np,0)
                      IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                    ENDDO !nonp
                    CALL ASSERT(ne.NE.0,
     '                '>>Could not find an element in the region',
     '                ERROR,*9999)
                    CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                DIPOLE_DIR_NTIME,nk,NDIPOLES,np,
     '                NP_INTERFACE,nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                XP(1,1,1,np),EXACT,
     '                XP(1,1,1,np),XP(1,1,2,np),XP(1,1,3,np),
     '                ERROR,*9999)
                  ELSE
                    ERROR='>>Analytic solution type not implemented'
                    GOTO 9999
                  ENDIF
                ENDIF
              ENDIF
              IF((NP_INTERFACE(np,0).GT.1).AND.(nr.NE.
     '          NP_INTERFACE(np,1))) THEN
C               Interface node in slave region
                IF(NJT.EQ.2) THEN
                  IF(nk.EQ.2) THEN
C                   Reverse sign of s derivative at interface
                    EXACT=-EXACT
                  ENDIF
                ELSE IF(NJT.EQ.3) THEN
                  IF(nk.GE.3) THEN
C                   Reverse sign of s2 and s1s2 derivatives at interface
                    EXACT=-EXACT
                  ENDIF
                ENDIF
              ENDIF
              YP(ny,7)=EXACT
C GMH 7/11/95 Add the potential values at the nodes.
              IF(nk.EQ.1) THEN
                NUM_NODES=NUM_NODES+1
                SUM_ANAL=SUM_ANAL+YP(ny,7)
                IF(AVERAGE) THEN
                  SUM_SOL=SUM_SOL+YP(ny,1)
                ENDIF
              ENDIF
            ENDDO ! nk
          ENDDO ! nonode (np)








C*** Derivatives

          nc=2 ! derivatives

          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            DO nk=1,MAX(NKH(NH_LOC(1,nx),np,nc,nr)-KTYP93(nc,nr),1)
              ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nr) !global var number
C cpb 30/10/96 Get tangents and normals from function
C              IF(.NOT.(NJT.EQ.2.AND.(ANAL_CHOICE(nr).GE.5.AND.
C     '          ANAL_CHOICE(nr).LE.6).OR.(NJT.EQ.3.AND.
C     '          (ANAL_CHOICE(nr).GE.12.AND.ANAL_CHOICE(nr).LE.15))))
C     '          THEN !not dipole solutions
C
CC cpb 8/1/96 basing tangent on circular geometry if there is any.
C                IF(CALL_MESH.AND.JTYP14.EQ.2.AND.NJT.EQ.2) THEN
C                  RADIUS=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2)
C                  THETA=DATAN_MOD(XP(1,nv,1,np),XP(1,nv,2,np))
C                  IF(DABS(RADIUS-MESH3_RAD(1))/
C     '              MESH3_RAD(1).LE.SPHERE_RAD_TOL.AND.
C     '              NSPHERES.GT.1.AND.MESH3_BEMTYPE.EQ.2) THEN
CC                   First circle in an annuluar problem therefore
CC                   s is anticlockwise
C                    IF(ANAL_CHOICE(nr).EQ.4) THEN
CC                     Tangent vector for a circle in polar coordinates
C                      TANGENT(1,1)=0.0d0
C                      TANGENT(2,1)=1.0d0
C                    ELSE
C                      TANGENT(1,1)=-XP(1,1,2,np)/RADIUS
C                      TANGENT(2,1)=XP(1,1,1,np)/RADIUS
C                    ENDIF
C                  ELSE
CC                   s direction is clockwise
C                    IF(ANAL_CHOICE(nr).EQ.4) THEN
CC                     Tangent vector for a circle in polar coordinates
C                      TANGENT(1,1)=0.0d0
C                      TANGENT(2,1)=-1.0d0
C                    ELSE
C                      TANGENT(1,1)=XP(1,1,2,np)/RADIUS
C                      TANGENT(2,1)=-XP(1,1,1,np)/RADIUS
C                    ENDIF
C                  ENDIF
C                ELSE
CC                 Note: relies on hermite geometry to create tangents
C                  IF(NKH(NH_LOC(1,nx),np,nc).EQ.1) THEN
CC                   use geometric hermite interp
C                    DO nkk=2,NKJ(1,np)-KTYP93(nc,nr)
C                      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C                        TANGENT(nj,nkk-1)=XP(nkk,nv,nj,np)
C                      ENDDO !nj
C                    ENDDO !nkk
C                  ELSE
C                    DO nkk=2,NKH(NH_LOC(1,nx),np,nc)-KTYP93(nc,nr)
C                      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C                        TANGENT(nj,nkk-1)=XP(nkk,nv,nj,np)
C                      ENDDO !nj
C                    ENDDO !nkk
C                  ENDIF
C                ENDIF
C              ENDIF
              IF(.NOT.(NJT.EQ.2.AND.(ANAL_CHOICE(nr).GE.5.AND.
     '          ANAL_CHOICE(nr).LE.6).OR.(NJT.EQ.3.AND.
     '          (ANAL_CHOICE(nr).GE.12.AND.ANAL_CHOICE(nr).LE.15)))
     '          .OR.(nk.NE.1))
     '          THEN !not dipole solutions
                IF(ANAL_CHOICE(nr).EQ.4.AND.NJT.EQ.2) THEN
                  IF(CALL_MESH.AND.JTYP14.EQ.2.AND.NJT.EQ.2) THEN
                    RADIUS=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2)
                    THETA=DATAN_MOD(XP(1,nv,1,np),XP(1,nv,2,np))
                    IF(DABS(RADIUS-MESH3_RAD(1))/
     '                MESH3_RAD(1).LE.SPHERE_RAD_TOL.AND.
     '                NSPHERES.GT.1.AND.MESH3_BEMTYPE.EQ.2) THEN
C                     First circle in an annuluar problem therefore
C                     s is anticlockwise
C                     Tangent vector for a circle in polar coordinates
                      TANGENT(1,1)=0.0d0
                      TANGENT(2,1)=1.0d0
                    ELSE
C                     s direction is clockwise
C                     Tangent vector for a circle in polar coordinates
                      TANGENT(1,1)=0.0d0
                      TANGENT(2,1)=-1.0d0
                    ENDIF
                  ENDIF
                ELSE
                  CALL GET_TNVECTOR(IBT,IDO,INP,NBJ,NENP,NKJE,np,NPF,
     '              NP_INTERFACE,NPNE,nr,NRE,NVJE,nx,XN_LOCAL,SE,
     '              TANGENT,XA,XE,XP,ERROR,*9999)
                ENDIF
              ENDIF
              IF(NJT.EQ.2) THEN

C cpb 30/10/96 Get tangents and normals from function
C                IF(.NOT.(ANAL_CHOICE(nr).GE.5.AND.
C     '            ANAL_CHOICE(nr).LE.6)) THEN !not dipole solutions
C                  XN_LOCAL(1)=-TANGENT(2,1)
C                  XN_LOCAL(2)=TANGENT(1,1)
C                  VLENGTH=DSQRT(XN_LOCAL(1)*XN_LOCAL(1)+XN_LOCAL(2)*
C     '              XN_LOCAL(2))
C                  IF(VLENGTH.LE.1.0d-6) THEN
C                    ERROR='>>Normal vector length=0'
C                    GOTO 9999
C                  ENDIF
C                  XN_LOCAL(1)=XN_LOCAL(1)/VLENGTH
C                  XN_LOCAL(2)=XN_LOCAL(2)/VLENGTH
CC cpb 28/1/96 Correct sign of exact at the end rather than reverse
CC normal
CC                  IF((NP_INTERFACE(np,0).GT.1).AND.(NP_INTERFACE(np,1)
CC     '              .NE.nr)) THEN
CC                    XN_LOCAL(1)=-XN_LOCAL(1)
CC                    XN_LOCAL(2)=-XN_LOCAL(2)
CC                  ENDIF
C                ENDIF
                IF(nk.EQ.1) THEN
                  IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                    GRADF(1)=ANAL_K
                    GRADF(2)=-ANAL_K
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(2)*XN_LOCAL(2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(x^2-y^2)
                    GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(2)=-2.0d0*ANAL_K*XP(1,nv,2,np)
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(2)*XN_LOCAL(2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2+2xy-y^2)
                    GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,1,np)+XP(1,nv,2,np))
                    GRADF(2)=2.0d0*ANAL_K*(XP(1,nv,1,np)-XP(1,nv,2,np))
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(2)*XN_LOCAL(2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(r^n.cos(n.t) ...
                    GRADF(1)=ANAL_K*(ANAL_A*ANAL_N*RADIUS**(ANAL_N-1)*
     '                DCOS(ANAL_N*THETA)+ANAL_B*ANAL_M*
     '                RADIUS**(ANAL_M-1)*DSIN(ANAL_M*THETA)+ANAL_C*
     '                DCOS(THETA)+ANAL_D*DSIN(THETA))
                    GRADF(2)=ANAL_K*(-ANAL_A*ANAL_N*RADIUS**(ANAL_N-1)*
     '                DSIN(ANAL_N*THETA)+ANAL_B*ANAL_M*
     '                RADIUS**(ANAL_M-1)*DCOS(ANAL_M*THETA)-ANAL_C*
     '                DSIN(THETA)+ANAL_D*DCOS(THETA))
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(2)*XN_LOCAL(2)
                  ELSE IF(ANAL_CHOICE(nr).GE.5.AND.
     '                ANAL_CHOICE(nr).LE.6) THEN !Dipole solutions
                    ne=0
                    DO nonp=1,NENP(np,0)
                      IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                    ENDDO !nonp
                    CALL ASSERT(ne.NE.0,
     '                '>>Could not find an element in the region',
     '                ERROR,*9999)
                    CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                DIPOLE_DIR_NTIME,4,NDIPOLES,np,
     '                NP_INTERFACE,nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                XP(1,1,1,np),
     '                EXACT,XP(1,1,1,np),XP(1,1,2,np),XP(1,1,1,np),
     '                ERROR,*9999)
                  ELSE
                    ERROR='>>Analytic solution type not implemented'
                    GOTO 9999
                  ENDIF
                ELSE ! Assumes Hermite interpolation in this case
C*** Assumes that the domain is a circle in 2d or a sphere 3d.
C*** NOTE: For a region between 2 spheres or circles only the magnitude
C*** of the solution can be given. It is not known if a node lies on
C*** the inner or outer circle/sphere and so the analytic solution
C*** (which requires this information) is determined only up to the
C*** sign. For this reason only magnitudes are compared for a single
C*** region problem.

                  IF(BEMCURVATURECORRECTION) THEN

C cpb 18/1/97 On the slave surface of an interface only the normal
C is reverse and not the s direction. Hence if we are on the slave
C slave surface we need to reverse the s direction for the analytic
C case.
                    IF(NP_INTERFACE(np,0).GT.1.AND.
     '                NP_INTERFACE(np,1).NE.nr) THEN
                      TANGENT(1,1)=-TANGENT(1,1)
                      TANGENT(2,1)=-TANGENT(2,1)
                    ENDIF

                    IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                      EXACT=0.0d0
                    ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(x^2-y^2)
                      HESSIAN(1,1)=2.0d0*ANAL_K
                      HESSIAN(1,2)=0.0d0
                      HESSIAN(2,1)=0.0d0
                      HESSIAN(2,2)=-2.0d0*ANAL_K
                      EXACT=
     '                  TANGENT(1,1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(1,1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,1)*HESSIAN(2,2)*XN_LOCAL(2)
                    ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2+2xy-y^2)
                      HESSIAN(1,1)=2.0d0*ANAL_K
                      HESSIAN(1,2)=2.0d0*ANAL_K
                      HESSIAN(2,1)=2.0d0*ANAL_K
                      HESSIAN(2,2)=-2.0d0*ANAL_K
                      EXACT=
     '                  TANGENT(1,1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(1,1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,1)*HESSIAN(2,2)*XN_LOCAL(2)
                    ELSE IF(ANAL_CHOICE(nr).EQ.5.OR.
     '                  ANAL_CHOICE(nr).EQ.6) THEN !dipole cases
                      ne=0
                      DO nonp=1,NENP(np,0)
                        IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                      ENDDO !nonp
                      CALL ASSERT(ne.NE.0,
     '                  '>>Could not find an element in the region',
     '                  ERROR,*9999)
                      CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                  DIPOLE_DIR_NTIME,5,NDIPOLES,np,NP_INTERFACE,
     '                  nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                  XP(1,1,1,np),EXACT,XP(1,1,1,np),
     '                  XP(1,1,2,np),XP(1,1,1,np),ERROR,*9999)
                    ELSE
                      ERROR='>>Analytic solution type not implemented'
                      GOTO 9999
                    ENDIF
                  ELSE

C cpb 18/1/97 Try and account for sign of arc-length derivative of
C the flux. The only problem that really remains is that of the inner
C circle/surface of a annuluar problem. Check the dot product of the
C normal with the radius vector. If the dot product is negative and we
C are not on an interface then the normal is inward (i.e. we are on
C the inner surface) and hence reverse the sign of the arc-length
C derivative of the flux for those cases where the analytic formuluae
C has be derived for an outward normal.

                    RDOTN=XP(1,nv,1,np)*XN_LOCAL(1)+XP(1,nv,2,np)*
     '                XN_LOCAL(2)
                    IF(RDOTN.LT.0.0d0.AND.NP_INTERFACE(np,1).EQ.nr)
     '                THEN
                      SIGN=-1.0d0
                    ELSE
                      SIGN=1.0d0
                    ENDIF

                    RAD=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2)
                    IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                      GRADF(1)=ANAL_K/RAD
                      GRADF(2)=-ANAL_K/RAD
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,1)+GRADF(2)*
     '                  TANGENT(2,1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(x^2-y^2)
                      GRADF(1)=4.0d0*ANAL_K*XP(1,nv,1,np)/RAD
                      GRADF(2)=-4.0d0*ANAL_K*XP(1,nv,2,np)/RAD
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,1)+GRADF(2)*
     '                  TANGENT(2,1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2+2xy-y^2)
                      GRADF(1)=4.0d0*ANAL_K/RAD*(XP(1,nv,1,np)+
     '                  XP(1,nv,2,np))
                      GRADF(2)=4.0d0*ANAL_K/RAD*(XP(1,nv,1,np)-
     '                  XP(1,nv,2,np))
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,1)+GRADF(2)*
     '                  TANGENT(2,1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(r^n.cos(n...
                      GRADF(1)=ANAL_K*(ANAL_A*ANAL_N*(ANAL_N-1)*
     '                  RADIUS**(ANAL_N-2)*DCOS(ANAL_N*THETA)+ANAL_B*
     '                  ANAL_M*(ANAL_M-1)*RADIUS**(ANAL_M-2)*
     '                  DSIN(ANAL_M*THETA))
                      GRADF(2)=ANAL_K*(-ANAL_A*ANAL_N**2*
     '                  RADIUS**(ANAL_N-2)*DSIN(ANAL_N*THETA)+ANAL_B*
     '                  ANAL_M**2*RADIUS**(ANAL_M-2)*
     '                  DCOS(ANAL_M*THETA)-ANAL_C/RADIUS*DSIN(THETA)+
     '                  ANAL_D/RADIUS*DCOS(THETA))
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,1)+
     '                  GRADF(2)*TANGENT(2,1))
                    ELSE IF(ANAL_CHOICE(nr).GE.5.AND.
     '                  ANAL_CHOICE(nr).LE.6) THEN !Dipole
                      ne=0
                      DO nonp=1,NENP(np,0)
                        IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                      ENDDO !nonp
                      CALL ASSERT(ne.NE.0,
     '                  '>>Could not find an element in the region',
     '                  ERROR,*9999)
                      CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                  DIPOLE_DIR_NTIME,5,NDIPOLES,np,
     '                  NP_INTERFACE,nr,nx,CE(1,ne),DIPOLE_CEN,
     '                  DIPOLE_DIR,XP(1,1,1,np),
     '                  EXACT,XP(1,1,1,np),XP(1,1,2,np),XP(1,1,1,np),
     '                  ERROR,*9999)
                    ELSE
                      ERROR='>>Analytic solution type not implemented'
                      GOTO 9999
                    ENDIF
                  ENDIF
                ENDIF ! nk


              ELSE IF(NJT.EQ.3) THEN
C cpb 30/10/96 Get tangents and normals from function
C                IF(.NOT.(ANAL_CHOICE(nr).GE.12.AND.ANAL_CHOICE(nr).LE.
C     '            15)) THEN !not dipole solutions
C                  XN_LOCAL(1)=TANGENT(2,1)*TANGENT(3,2)-TANGENT(3,1)*
C     '              TANGENT(2,2)
C                  XN_LOCAL(2)=TANGENT(3,1)*TANGENT(1,2)-TANGENT(1,1)*
C     '              TANGENT(3,2)
C                  XN_LOCAL(3)=TANGENT(1,1)*TANGENT(2,2)-TANGENT(2,1)*
C     '              TANGENT(1,2)
C                  VLENGTH=DSQRT(XN_LOCAL(1)*XN_LOCAL(1)+
C     '              XN_LOCAL(2)*XN_LOCAL(2)+XN_LOCAL(3)*XN_LOCAL(3))
C                  IF(VLENGTH.LE.1.0d-6) THEN
C                    ERROR='>>Normal vector length=0'
C                    GOTO 9999
C                  ENDIF
C                  XN_LOCAL(1)=XN_LOCAL(1)/VLENGTH
C                  XN_LOCAL(2)=XN_LOCAL(2)/VLENGTH
C                  XN_LOCAL(3)=XN_LOCAL(3)/VLENGTH
CC cpb 28/1/96 Correct sign of exact at the end rather than reverse
CC normal
CC                  IF((NP_INTERFACE(np,0).GT.1).AND.(NP_INTERFACE(np,1)
CC     '              .NE.nr)) THEN
CC                    XN_LOCAL(1)=-XN_LOCAL(1)
CC                    XN_LOCAL(2)=-XN_LOCAL(2)
CC                    XN_LOCAL(3)=-XN_LOCAL(3)
CC                  ENDIF ! np_interfaces
C                ENDIF
                IF(nk.EQ.1) THEN
                  IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                    GRADF(1)=ANAL_K
                    GRADF(2)=-ANAL_K
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(2)*XN_LOCAL(2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(z)
                    GRADF(3)=ANAL_K
                    EXACT=GRADF(3)*XN_LOCAL(3)
                  ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2-y^2)
                    GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(2)=-2.0d0*ANAL_K*XP(1,nv,2,np)
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(2)*XN_LOCAL(2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(x^2-z^2)
                    GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(3)=-2.0d0*ANAL_K*XP(1,nv,3,np)
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(3)*XN_LOCAL(3)
                  ELSE IF(ANAL_CHOICE(nr).EQ.5) THEN !K(y^2-z^2)
                    GRADF(2)=2.0d0*ANAL_K*XP(1,nv,2,np)
                    GRADF(3)=-2.0d0*ANAL_K*XP(1,nv,3,np)
                    EXACT=GRADF(2)*XN_LOCAL(2)+GRADF(3)*XN_LOCAL(3)
                  ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !K(x^2+y^2-2z^2)
                    GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(2)=2.0d0*ANAL_K*XP(1,nv,2,np)
                    GRADF(3)=-4.0d0*ANAL_K*XP(1,nv,3,np)
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(2)*XN_LOCAL(2)+
     '                GRADF(3)*XN_LOCAL(3)
                  ELSE IF(ANAL_CHOICE(nr).EQ.7) THEN !K(x^2-2y^2+z^2)
                    GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(2)=-4.0d0*ANAL_K*XP(1,nv,2,np)
                    GRADF(3)=2.0d0*ANAL_K*XP(1,nv,3,np)
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(2)*XN_LOCAL(2)+
     '                GRADF(3)*XN_LOCAL(3)
                  ELSE IF(ANAL_CHOICE(nr).EQ.8) THEN !K(-2x^2+y^2+z^2)
                    GRADF(1)=-4.0d0*ANAL_K*XP(1,nv,1,np)
                    GRADF(2)=2.0d0*ANAL_K*XP(1,nv,2,np)
                    GRADF(3)=2.0d0*ANAL_K*XP(1,nv,3,np)
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(2)*XN_LOCAL(2)+
     '                GRADF(3)*XN_LOCAL(3)
                  ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN !K(x^2+2xy-y^2)
                    GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,1,np)+XP(1,nv,2,np))
                    GRADF(2)=2.0d0*ANAL_K*(XP(1,nv,1,np)-XP(1,nv,2,np))
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(2)*XN_LOCAL(2)
                  ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN !K(x^2+2xz-z^2)
                    GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,1,np)+XP(1,nv,3,np))
                    GRADF(3)=2.0d0*ANAL_K*(XP(1,nv,1,np)-XP(1,nv,3,np))
                    EXACT=GRADF(1)*XN_LOCAL(1)+GRADF(3)*XN_LOCAL(3)
                  ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN !K(y^2+2yz-z^2)
                    GRADF(2)=2.0d0*ANAL_K*(XP(1,nv,2,np)+XP(1,nv,3,np))
                    GRADF(3)=2.0d0*ANAL_K*(XP(1,nv,2,np)-XP(1,nv,3,np))
                    EXACT=GRADF(2)*XN_LOCAL(2)+GRADF(3)*XN_LOCAL(3)
                  ELSE IF(ANAL_CHOICE(nr).GE.12.AND
     '                .ANAL_CHOICE(nr).LE.15) THEN
C                   Dipole solutions
                    ne=0
                    DO nonp=1,NENP(np,0)
                      IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                    ENDDO !nonp
                    CALL ASSERT(ne.NE.0,
     '                '>>Could not find an element in the region',
     '                ERROR,*9999)
                    CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                DIPOLE_DIR_NTIME,4,NDIPOLES,np,
     '                NP_INTERFACE,nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                XP(1,1,1,np),EXACT,
     '                XP(1,1,1,np),XP(1,1,2,np),XP(1,1,3,np),
     '                ERROR,*9999)
                    IF((NP_INTERFACE(np,0).GT.1).AND.
     '                (NP_INTERFACE(np,1).NE.nr)) THEN
C                     To get correct sign on the analytic solution at
C                     an interface
                      EXACT=-EXACT
                    ENDIF
                  ELSE
                    ERROR='>>Analytic solution type not implemented'
                    GOTO 9999
                  ENDIF ! idata
                ELSE
                  IF(BEMCURVATURECORRECTION) THEN

C cpb 18/1/97 On the slave surface of an interface only the normal
C is reverse and not the s2 direction. Hence if we are on the slave
C slave surface we need to reverse the s2 direction for the analytic
C case.
                    IF(NP_INTERFACE(np,0).GT.1.AND.
     '                NP_INTERFACE(np,1).NE.nr) THEN
                      TANGENT(1,2)=-TANGENT(1,2)
                      TANGENT(2,2)=-TANGENT(2,2)
                      TANGENT(3,2)=-TANGENT(3,2)
                    ENDIF

                    IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                      EXACT=0.0d0
                    ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(z)
                      EXACT=0.0d0
                    ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2-y^2)
                      HESSIAN(1,1)=2.0d0*ANAL_K
                      HESSIAN(1,2)=0.0d0
                      HESSIAN(1,3)=0.0d0
                      HESSIAN(2,1)=0.0d0
                      HESSIAN(2,2)=-2.0d0*ANAL_K
                      HESSIAN(2,3)=0.0d0
                      HESSIAN(3,1)=0.0d0
                      HESSIAN(3,2)=0.0d0
                      HESSIAN(3,3)=0.0d0
                      EXACT=
     '                  TANGENT(1,nk-1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,1)*XN_LOCAL(1)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,2)*XN_LOCAL(2)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,2)*XN_LOCAL(2)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,3)*XN_LOCAL(3)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,3)*XN_LOCAL(3)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,3)*XN_LOCAL(3)
                    ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(x^2-z^2)
                      HESSIAN(1,1)=2.0d0*ANAL_K
                      HESSIAN(1,2)=0.0d0
                      HESSIAN(1,3)=0.0d0
                      HESSIAN(2,1)=0.0d0
                      HESSIAN(2,2)=0.0d0
                      HESSIAN(2,3)=0.0d0
                      HESSIAN(3,1)=0.0d0
                      HESSIAN(3,2)=0.0d0
                      HESSIAN(3,3)=-2.0d0*ANAL_K
                      EXACT=
     '                  TANGENT(1,nk-1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,1)*XN_LOCAL(1)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,2)*XN_LOCAL(2)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,2)*XN_LOCAL(2)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,3)*XN_LOCAL(3)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,3)*XN_LOCAL(3)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,3)*XN_LOCAL(3)
                    ELSE IF(ANAL_CHOICE(nr).EQ.5) THEN !K(y^2-z^2)
                      HESSIAN(1,1)=0.0d0
                      HESSIAN(1,2)=0.0d0
                      HESSIAN(1,3)=0.0d0
                      HESSIAN(2,1)=0.0d0
                      HESSIAN(2,2)=2.0d0*ANAL_K
                      HESSIAN(2,3)=0.0d0
                      HESSIAN(3,1)=0.0d0
                      HESSIAN(3,2)=0.0d0
                      HESSIAN(3,3)=-2.0d0*ANAL_K
                      EXACT=
     '                  TANGENT(1,nk-1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,1)*XN_LOCAL(1)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,2)*XN_LOCAL(2)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,2)*XN_LOCAL(2)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,3)*XN_LOCAL(3)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,3)*XN_LOCAL(3)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,3)*XN_LOCAL(3)
                    ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !K(x^2+y^2-2z^2)
                      HESSIAN(1,1)=2.0d0*ANAL_K
                      HESSIAN(1,2)=0.0d0
                      HESSIAN(1,3)=0.0d0
                      HESSIAN(2,1)=0.0d0
                      HESSIAN(2,2)=2.0d0*ANAL_K
                      HESSIAN(2,3)=0.0d0
                      HESSIAN(3,1)=0.0d0
                      HESSIAN(3,2)=0.0d0
                      HESSIAN(3,3)=-4.0d0*ANAL_K
                      EXACT=
     '                  TANGENT(1,nk-1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,1)*XN_LOCAL(1)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,2)*XN_LOCAL(2)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,2)*XN_LOCAL(2)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,3)*XN_LOCAL(3)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,3)*XN_LOCAL(3)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,3)*XN_LOCAL(3)
                    ELSE IF(ANAL_CHOICE(nr).EQ.7) THEN !K(x^2-2y^2+z^2)
                      HESSIAN(1,1)=2.0d0*ANAL_K
                      HESSIAN(1,2)=0.0d0
                      HESSIAN(1,3)=0.0d0
                      HESSIAN(2,1)=0.0d0
                      HESSIAN(2,2)=-4.0d0*ANAL_K
                      HESSIAN(2,3)=0.0d0
                      HESSIAN(3,1)=0.0d0
                      HESSIAN(3,2)=0.0d0
                      HESSIAN(3,3)=2.0d0*ANAL_K
                      EXACT=
     '                  TANGENT(1,nk-1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,1)*XN_LOCAL(1)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,2)*XN_LOCAL(2)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,2)*XN_LOCAL(2)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,3)*XN_LOCAL(3)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,3)*XN_LOCAL(3)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,3)*XN_LOCAL(3)
                    ELSE IF(ANAL_CHOICE(nr).EQ.8) THEN !K(-2x^2+y^2+z^2)
                      HESSIAN(1,1)=-4.0d0*ANAL_K
                      HESSIAN(1,2)=0.0d0
                      HESSIAN(1,3)=0.0d0
                      HESSIAN(2,1)=0.0d0
                      HESSIAN(2,2)=2.0d0*ANAL_K
                      HESSIAN(2,3)=0.0d0
                      HESSIAN(3,1)=0.0d0
                      HESSIAN(3,2)=0.0d0
                      HESSIAN(3,3)=2.0d0*ANAL_K
                      EXACT=
     '                  TANGENT(1,nk-1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,1)*XN_LOCAL(1)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,2)*XN_LOCAL(2)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,2)*XN_LOCAL(2)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,3)*XN_LOCAL(3)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,3)*XN_LOCAL(3)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,3)*XN_LOCAL(3)
                    ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN !K(x^2+2xy-y^2)
                      HESSIAN(1,1)=2.0d0*ANAL_K
                      HESSIAN(1,2)=2.0d0*ANAL_K
                      HESSIAN(1,3)=0.0d0
                      HESSIAN(2,1)=2.0d0*ANAL_K
                      HESSIAN(2,2)=-2.0d0*ANAL_K
                      HESSIAN(2,3)=0.0d0
                      HESSIAN(3,1)=0.0d0
                      HESSIAN(3,2)=0.0d0
                      HESSIAN(3,3)=0.0d0
                      EXACT=
     '                  TANGENT(1,nk-1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,1)*XN_LOCAL(1)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,2)*XN_LOCAL(2)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,2)*XN_LOCAL(2)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,3)*XN_LOCAL(3)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,3)*XN_LOCAL(3)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,3)*XN_LOCAL(3)
                    ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN !K(x^2+2xz-z^2)
                      HESSIAN(1,1)=2.0d0*ANAL_K
                      HESSIAN(1,2)=0.0d0
                      HESSIAN(1,3)=2.0d0*ANAL_K
                      HESSIAN(2,1)=0.0d0
                      HESSIAN(2,2)=0.0d0
                      HESSIAN(2,3)=0.0d0
                      HESSIAN(3,1)=2.0d0*ANAL_K
                      HESSIAN(3,2)=0.0d0
                      HESSIAN(3,3)=-2.0d0*ANAL_K
                      EXACT=
     '                  TANGENT(1,nk-1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,1)*XN_LOCAL(1)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,2)*XN_LOCAL(2)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,2)*XN_LOCAL(2)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,3)*XN_LOCAL(3)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,3)*XN_LOCAL(3)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,3)*XN_LOCAL(3)
                    ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN !K(y^2+2yz-z^2)
                      HESSIAN(1,1)=0.0d0
                      HESSIAN(1,2)=0.0d0
                      HESSIAN(1,3)=0.0d0
                      HESSIAN(2,1)=0.0d0
                      HESSIAN(2,2)=2.0d0*ANAL_K
                      HESSIAN(2,3)=2.0d0*ANAL_K
                      HESSIAN(3,1)=0.0d0
                      HESSIAN(3,2)=2.0d0*ANAL_K
                      HESSIAN(3,3)=-2.0d0*ANAL_K
                      EXACT=
     '                  TANGENT(1,nk-1)*HESSIAN(1,1)*XN_LOCAL(1)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,1)*XN_LOCAL(1)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,1)*XN_LOCAL(1)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,2)*XN_LOCAL(2)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,2)*XN_LOCAL(2)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,2)*XN_LOCAL(2)+
     '                  TANGENT(1,nk-1)*HESSIAN(1,3)*XN_LOCAL(3)+
     '                  TANGENT(2,nk-1)*HESSIAN(2,3)*XN_LOCAL(3)+
     '                  TANGENT(3,nk-1)*HESSIAN(3,3)*XN_LOCAL(3)
C cpb 12/12/96 Adding single sphere until proper dipole cases done.
                    ELSE IF(ANAL_CHOICE(nr).EQ.12) THEN
                      EXACT=0.0d0 !no flux
                    ELSE
                      ERROR='>>Analytic solution type not implemented'
                      GOTO 9999
                    ENDIF ! idata
                  ELSE
C cpb 18/1/97 Try and account for sign of arc-length derivative of
C the flux. The only problem that really remains is that of the inner
C circle/surface of a annuluar problem. Check the dot product of the
C normal with the radius vector. If the dot product is negative and we
C are not on an interface then the normal is inward (i.e. we are on
C the inner surface) and hence reverse the sign of the arc-length
C derivative of the flux for those cases where the analytic formuluae
C has be derived for an outward normal.

                    RDOTN=XP(1,nv,1,np)*XN_LOCAL(1)+XP(1,nv,2,np)*
     '                XN_LOCAL(2)+XP(1,nv,3,np)*XN_LOCAL(3)
                    IF(RDOTN.LT.0.0d0.AND.NP_INTERFACE(np,1).EQ.nr)
     '                THEN
                      SIGN=-1.0d0
                    ELSE
                      SIGN=1.0d0
                    ENDIF
                    RAD=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
     '                XP(1,nv,3,np)**2)
                    IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                      GRADF(1)=ANAL_K/RAD
                      GRADF(2)=-ANAL_K/RAD
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(z)
                      GRADF(3)=ANAL_K/RAD
                      EXACT=SIGN*(GRADF(3)*TANGENT(3,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2-y^2)
                      GRADF(1)=4.0d0*ANAL_K*XP(1,nv,1,np)/RAD
                      GRADF(2)=-4.0d0*ANAL_K*XP(1,nv,2,np)/RAD
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(x^2-z^2)
                      GRADF(1)=4.0d0*ANAL_K*XP(1,nv,1,np)/RAD
                      GRADF(3)=-4.0d0*ANAL_K*XP(1,nv,3,np)/RAD
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,nk-1)+GRADF(3)*
     '                  TANGENT(3,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.5) THEN !K(y^2-z^2)
                      GRADF(2)=4.0d0*ANAL_K*XP(1,nv,2,np)/RAD
                      GRADF(3)=-4.0d0*ANAL_K*XP(1,nv,3,np)/RAD
                      EXACT=SIGN*(GRADF(2)*TANGENT(2,nk-1)+GRADF(3)*
     '                  TANGENT(3,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !K(x^2+y^2-2z^2)
                      GRADF(1)=4.0d0*ANAL_K*XP(1,nv,1,np)/RAD
                      GRADF(2)=4.0d0*ANAL_K*XP(1,nv,2,np)/RAD
                      GRADF(3)=-8.0d0*ANAL_K*XP(1,nv,3,np)/RAD
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)+GRADF(3)*TANGENT(3,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.7) THEN !K(x^2-2y^2+z^2)
                      GRADF(1)=4.0d0*ANAL_K*XP(1,nv,1,np)/RAD
                      GRADF(2)=-8.0d0*ANAL_K*XP(1,nv,2,np)/RAD
                      GRADF(3)=4.0d0*ANAL_K*XP(1,nv,3,np)/RAD
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)+GRADF(3)*TANGENT(3,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.8) THEN !K(-2x^2+y^2+z^2)
                      GRADF(1)=-8.0d0*ANAL_K*XP(1,nv,1,np)/RAD
                      GRADF(2)=4.0d0*ANAL_K*XP(1,nv,2,np)/RAD
                      GRADF(3)=4.0d0*ANAL_K*XP(1,nv,3,np)/RAD
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)+GRADF(3)*TANGENT(3,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN !K(x^2+2xy-y^2)
                      GRADF(1)=4.0d0*ANAL_K*(XP(1,nv,1,np)+
     '                  XP(1,nv,2,np))/RAD
                      GRADF(2)=4.0d0*ANAL_K*(XP(1,nv,1,np)-
     '                  XP(1,nv,2,np))/RAD
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN !K(x^2+2xz-z^2)
                      GRADF(1)=4.0d0*ANAL_K*(XP(1,nv,1,np)+
     '                  XP(1,nv,3,np))/RAD
                      GRADF(3)=4.0d0*ANAL_K*(XP(1,nv,1,np)-
     '                  XP(1,nv,3,np))/RAD
                      EXACT=SIGN*(GRADF(1)*TANGENT(1,nk-1)+GRADF(3)*
     '                  TANGENT(3,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN !K(y^2+2yz-z^2)
                      GRADF(2)=4.0d0*ANAL_K*(XP(1,nv,2,np)+
     '                  XP(1,nv,3,np))/RAD
                      GRADF(3)=4.0d0*ANAL_K*(XP(1,nv,2,np)-
     '                  XP(1,nv,3,np))/RAD
                      EXACT=SIGN*(GRADF(2)*TANGENT(2,nk-1)+GRADF(3)*
     '                  TANGENT(3,nk-1))
                    ELSE IF(ANAL_CHOICE(nr).GE.12.AND
     '                  .ANAL_CHOICE(nr).LE.15) THEN
C                   Dipole solutions
                      ne=0
                      DO nonp=1,NENP(np,0)
                        IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                      ENDDO !nonp
                      CALL ASSERT(ne.NE.0,
     '                  '>>Could not find an element in the region',
     '                  ERROR,*9999)
                      CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                  DIPOLE_DIR_NTIME,nk+3,NDIPOLES,np,NP_INTERFACE,
     '                  nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                  XP(1,1,1,np),EXACT,
     '                  XP(1,1,1,np),XP(1,1,2,np),XP(1,1,3,np),ERROR,
     '                  *9999)
                    ELSE
                      ERROR='>>Analytic solution type not implemented'
                      GOTO 9999
                    ENDIF ! idata
                  ENDIF
                ENDIF ! nk
              ENDIF ! njt
C cpb 30/10/96 Now taken care of in the get_tnvector function
C              IF((NP_INTERFACE(np,0).GT.1).AND.(nr.NE.
C     '          NP_INTERFACE(np,1))) THEN
CC               Interface node in slave region
C                IF(NJT.EQ.2) THEN
C                  IF(nk.EQ.1) THEN
CC                   Reverse sign of normal derivative at interface
C                    EXACT=-EXACT
C                  ENDIF ! nk
C                ELSE IF(NJT.EQ.3) THEN
C                  IF(nk.NE.2) THEN
CC                   Reverse sign of normal and s2 derivatives at
CC                   interface
C                    EXACT=-EXACT
C                  ENDIF ! nk
C                ENDIF ! njt
C              ENDIF ! np_interface
              YP(ny,7)=EXACT
C cpb 30/10/96 Now taken care of in the get_tnvector function
C              IF(CALL_MESH.AND.JTYP14.EQ.2.AND.NSPHERES.GT.1
C     '          .AND.MESH3_BEMTYPE.EQ.2) THEN
CC               Circular/sphererical bem annuluar mesh
C                IF(nk.EQ.2) THEN
C                  IF(NJT.EQ.2) THEN
C                    RADIUS=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2)
C                  ELSE IF(NJT.EQ.3) THEN
C                    RADIUS=DSQRT(XP(1,1,1,np)**2+XP(1,1,2,np)**2+
C     '                XP(1,1,3,np)**2)
C                  ENDIF
C                  IF(DABS(RADIUS-MESH3_RAD(1))/
C     '              MESH3_RAD(1).LE.SPHERE_RAD_TOL
C     '              .AND.NSPHERES.GT.1) THEN !on inner annulus surface
CC                                             so reverse sign
C                    YP(ny,7)=-YP(ny,7)
C                  ENDIF
C                ENDIF
C              ENDIF
            ENDDO ! nk
          ENDDO ! nonode (np)
          IF(AVERAGE.AND.NUM_NODES.GT.0) THEN
            nc=1
            OFFSET=(SUM_SOL-SUM_ANAL)/DBLE(NUM_NODES)
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(IPANA9_1)
              WRITE(OP_STRING,'(/'' Potential offset='',D12.4)')
     '          OFFSET
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(IPANA9_1)
            ENDIF
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              DO nk=1,MAX(NKH(NH_LOC(1,nx),np,nc,nr)-KTYP93(nc,nr),1)
                ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nr) !global var number
                IF(nk.EQ.1) THEN
                  YP(ny,7)=YP(ny,7)+OFFSET
                ENDIF
              ENDDO ! nk
            ENDDO ! nonode (np)
          ENDIF


C news AJP 5/10/99
        ELSEIF (ITYP2(nr,nx).EQ.5) THEN !Poisson equation
          IF(ITYP3(nr,nx).EQ.2) THEN !special rhs term
            K=CE(1,NEELEM(1,nr))
            GI=CE(2,NEELEM(1,nr))
            GE=K-GI
            CALL ASSERT(DABS(K).GT.RDELTA,'>>Zero k material value',
     '        ERROR,*9999)
            CALL ASSERT(DABS(GE).GT.RDELTA,'>>Zero s material value',
     '        ERROR,*9999)
            IF(NJT.EQ.2) THEN
              IF(IOTYPE.EQ.3) THEN
                IDATA(1)=ANAL_CHOICE(nr)
              ENDIF
              FORMAT='('' Specify the analytic solution [1]:'''//
     '          '/''   (1) A*r+B*log(r)+C'''//
     '          '/''   (2) Bidomain transfer test'''//
     '          '/''   (3) Double circle problem'''//
     '          '/$,''   '',I1)'
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_CHOICE(nr)=IDATA(1)

              IF(ANAL_CHOICE(nr).EQ.1) THEN
                FORMAT='(/$,'' Enter the value of A [-1.0]: '',D12.4)'
                RDEFLT(1)=-1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_A
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_A=RDATA(1)
                FORMAT='(/$,'' Enter the value of B  [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_B
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_B=RDATA(1)
                FORMAT='(/$,'' Enter the value of C  [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_C
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_C=RDATA(1)
              ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN
C MLB 2-August-2000 adding 2d circles equivalent to Leo's 3d
C   spheres
                CALL_ANAL_TRSF=.TRUE.

              ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN
C CPB 18-JAN-2001 Adding double circle case which corresponds to the the
C double sphere case.
                FORMAT='(/$,'' Enter the value of C_1 [2.0]: '',D12.4)'
                RDEFLT(1)=2.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_C_1
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_C_1=RDATA(1)
                FORMAT='(/$,'' Enter the value of D_1 [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_D_1
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_D_1=RDATA(1)
                FORMAT='(/$,'' Enter the value of H_1 [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_H_1
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_H_1=RDATA(1)
                FORMAT='(/$,'' Enter the value of the reference '
     '            //'potential [0.0]: '',D12.4)'
                RDEFLT(1)=0.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_REF
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_REF=RDATA(1)
                CALL_ANAL_TRSF=.TRUE.
              ENDIF
C LKC 15-MAR-2000 adding 3D analytic spheres
C   refer to LKC/AJP for specifics about each problem
            ELSE IF(NJT.EQ.3) THEN

              IF(IOTYPE.EQ.3) THEN
                IDATA(1)=ANAL_CHOICE(nr)
              ENDIF
C cpb 28/11/00 Adding double sphere test problem
              FORMAT='('' Specify the analytic solution [1]:'''//
     '          '/''   (1) Case 1b'''//
     '          '/''   (2) Case 3a'''//
     '          '/''   (3) Case 2b'''//
     '          '/''   (4) Case 1a'''//
     '          '/''   (5) Case 2a'''//
     '          '/''   (6) Case 3b'''//
     '          '/''   (7) Double sphere problem'''//
     '          '/''  *(8) Unused'''//
     '          '/''  *(9) Unused'''//
     '          '/'' *(10) Unused'''//
     '          '/$,''   '',I1)'

              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,7,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) ANAL_CHOICE(nr)=IDATA(1)


              IF(ANAL_CHOICE(nr).EQ.7) THEN
                FORMAT='(/$,'' Enter the value of C_11 [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_C_11
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_C_11=RDATA(1)
                FORMAT='(/$,'' Enter the value of D_11 [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_D_11
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_D_11=RDATA(1)
                FORMAT='(/$,'' Enter the value of H_11 [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_H_11
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_H_11=RDATA(1)
                FORMAT='(/$,'' Enter the value of C_1 [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_C_1
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_C_1=RDATA(1)
                FORMAT='(/$,'' Enter the value of the reference '
     '            //'potential [0.0]: '',D12.4)'
                RDEFLT(1)=0.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_REF
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,
     '            1,4,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_REF=RDATA(1)
              ENDIF
C LKC setting a special flag so it does not get confused with
C deanal activation
              CALL_ANAL_TRSF=.TRUE.

            ENDIF !njt
          ELSE
            ERROR='>>Equation type not implemented'
            GOTO 9999
          ENDIF !special rhs term


          IF(NJT.EQ.2) THEN
            IF(ANAL_CHOICE(nr).EQ.1) THEN
              nc=1 ! variables first
              NUM_NODES=0
              SUM_ANAL=0.0D0
              SUM_SOL=0.0D0
C Need to ensure that there is no-flux on the source term for this
C analytic solution to work.
C Also need to have a circular domain (in 2d).
C The no-flux condition is satisfied (on a circular domain) if a=-b/R
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                RADIUS=DSQRT(XP(1,1,1,NPNODE(1,nr))**2+
     '            XP(1,1,2,NPNODE(1,nr))**2)
                CALL ASSERT(ANAL_A+ANAL_B/RADIUS.LT.RDELTA,
     '            '>>Invalid a,b and radius',ERROR,*9999)
                DO nk=1,MAX(NKH(NH_LOC(1,nx),np,nc,nr)-KTYP93(nc,nr),1)
                  ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nr) !global var num
                  IF(NJT.EQ.2) THEN
                    IF(nk.EQ.1) THEN
                      IF(ANAL_CHOICE(nr).EQ.1) THEN !A*r+B*log(r)+C
                        EXACT=-GI/(GI+GE)*ANAL_A*RADIUS+
     '                    ANAL_B*DLOG(RADIUS)+
     '                    ANAL_C
                      ELSE
                        ERROR='>>Analytic solution type not implemented'
                        GOTO 9999
                      ENDIF
                    ELSE !nk.gt.1
                      CALL GET_TNVECTOR(IBT,IDO,INP,NBJ,NENP,NKJE,np,
     '                  NPF,NP_INTERFACE,NPNE,nr,NRE,NVJE,nx,XN_LOCAL,
     '                  SE,TANGENT,XA,XE,XP,ERROR,*9999)
                      IF(ANAL_CHOICE(nr).EQ.1) THEN !A*r+B*log(r)+C
                        GRADF(1)=-GI/(GI+GE)*ANAL_A+ANAL_B/RADIUS
                        GRADF(2)=0.0d0
                        EXACT=GRADF(1)*XP(1,1,1,NPNODE(1,nr))/RADIUS
                      ELSE
                        ERROR='>>Analytic solution type not implemented'
                        GOTO 9999
                      ENDIF
                    ENDIF ! nk
                  ELSE IF(NJT.EQ.3) THEN
                    ERROR='>>Analytic solution type not implemented'
                    IF(nk.EQ.1) THEN
                      IF(ANAL_CHOICE(nr).EQ.1) THEN
C*** Ksin(nx)sin(ny)sin(nz)
                        EXACT=ANAL_K*GI/(GI+GE)*
     '                    (DSIN(DBLE(ANAL_N)*XP(1,nv,1,np))*
     '                    DSIN(DBLE(ANAL_N)*XP(1,nv,2,np))*
     '                    DSIN(DBLE(ANAL_N)*XP(1,nv,3,np)))
                      ELSE
                        ERROR='>>Analytic solution type not implemented'
                        GOTO 9999
                      ENDIF
                    ELSE
                      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                        TANGENT(nj,nk-1)=XP(nk,nv,nj,np)
                      ENDDO
                      IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                        GRADF(1)=ANAL_K*GI/(GI+GE)*
     '                    (DCOS(DBLE(ANAL_N)*XP(1,nv,1,np))*
     '                    DSIN(DBLE(ANAL_N)*XP(1,nv,2,np))*
     '                    DSIN(DBLE(ANAL_N)*XP(1,nv,3,np)))
                        GRADF(2)=ANAL_K*GI/(GI+GE)*
     '                    (DSIN(DBLE(ANAL_N)*XP(1,nv,1,np))*
     '                    DCOS(DBLE(ANAL_N)*XP(1,nv,2,np))*
     '                    DSIN(DBLE(ANAL_N)*XP(1,nv,3,np)))
                        GRADF(3)=ANAL_K*GI/(GI+GE)*
     '                    (DSIN(DBLE(ANAL_N)*XP(1,nv,1,np))*
     '                    DSIN(DBLE(ANAL_N)*XP(1,nv,2,np))*
     '                    DCOS(DBLE(ANAL_N)*XP(1,nv,3,np)))
                        EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                    TANGENT(2,nk-1)+GRADF(3)*TANGENT(3,nk-1)
                      ELSE
                        ERROR='>>Analytic solution type not implemented'
                        GOTO 9999
                      ENDIF
                    ENDIF !nk
                  ENDIF !njt
                  YP(ny,7)=EXACT
                ENDDO !nk
              ENDDO !np
              nc=2 ! derivatives
C Must do derivatives.


            ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN
C***          Circles for testing transfer matrices
              nc=1           !variables first
              nv=1           !variables first
              nk=1           !linear basis
              nh=NH_LOC(1,nx)

C***          Calcuate the transmembrane potential
              RADIUS=1.0d0   !the inner radius
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny=NYNP(nk,nv,nh,np,0,nc,nr) !glob var num
                THETA=DACOS(XP(nk,nv,1,np)/RADIUS)
                IF(XP(nk,nv,2,np).LT.ZERO_TOL) THETA=-THETA

                COORD(1)=-5.0d0 !dipole strengths
                COORD(2)=5.0d0

C***            Analytic transmembrane values
                YP(ny,1)=-(((RADIUS**2.0d0)+1)/RADIUS)*
     '            ((COORD(1)*DCOS(THETA))+(COORD(2)*DSIN(THETA)))
                YP(ny,7)=YP(ny,1)
              ENDDO !nonode

C***          Calculate the extracellular solution
              nrr=2
              ANAL_CHOICE(nrr)=ANAL_CHOICE(nr)

              DO nonode=1,NPNODE(0,nrr)
                np=NPNODE(nonode,nrr)
                ny=NYNP(nk,nv,nh,np,0,nc,nrr) !glob var num
                IF(NP_INTERFACE(np,0).EQ.1) THEN
                  RADIUS=3.0d0
                ELSE
                  RADIUS=1.0d0
                ENDIF
                THETA=DACOS(XP(nk,nv,1,np)/RADIUS)
                IF(XP(nk,nv,2,np).LT.ZERO_TOL) THETA=-THETA

                YP(ny,7)=(DCOS(THETA)-DSIN(THETA))*
     '            (RADIUS+(9.0d0/RADIUS))-6.0d0
              ENDDO !nonode

C***          Analytic derivatives
              nc=2 !normal derivatives (fluxes)
              DO nonode=1,NPNODE(0,nrr)
                np=NPNODE(nonode,nrr)
                ny=NYNP(nk,nv,nh,np,0,nc,nrr) !glob var num
                IF(NP_INTERFACE(np,0).EQ.1) THEN
                  RADIUS=3.0d0
                ELSE
                  RADIUS=1.0d0
                ENDIF
                THETA=DACOS(XP(nk,nv,1,np)/RADIUS)
                IF(XP(nk,nv,2,np).LT.ZERO_TOL) THETA=-THETA

                YP(ny,7)=-(DCOS(THETA)-DSIN(THETA))*
     '            (1.0d0-(9.0d0/(RADIUS**2)))*CE(1,NEELEM(1,nrr))
              ENDDO !nonode
            ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN
              IF(NDIPOLES(1).EQ.0) THEN
                ERROR='>>Define a dipole in region 1 first'
                GOTO 9999
              ELSE IF(NDIPOLES(1).GT.1) THEN
                ERROR='>>Must only have one dipole in region 1'
                GOTO 9999
              ELSE
                IF(DIPOLE_CEN(1,0,1,1).NE.0.0d0.OR.
     '            DIPOLE_CEN(2,0,1,1).NE.0.0d0) THEN
                  ERROR='>>The dipole must be centred at the origin'
                  GOTO 9999
                ENDIF
              ENDIF
              RHO_X=DIPOLE_DIR(1,0,1,1)
              RHO_Y=DIPOLE_DIR(2,0,1,1)
              CALL ASSERT(CALL_MESH.AND.JTYP14.EQ.2,
     '          '>>Define a circular mesh first',ERROR,*9999)
              CALL ASSERT(NSPHERES.EQ.2,
     '          '>>Must have only two circles',ERROR,*9999)
              np=NPNODE(1,1)
              RAD_1=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2)
              np=NPNODE(NPNODE(0,2),2)
              RAD_2=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2)

              CALL ASSERT(DABS(RAD_1).GT.ZERO_TOL,
     '          '>>Circle 1 has a zero radius',ERROR,*9999)
              CALL ASSERT(DABS(RAD_2).GT.ZERO_TOL,
     '          '>>Circle 2 has a zero radius',ERROR,*9999)
              CALL ASSERT(DABS(RAD_2-RAD_1).GT.ZERO_TOL,
     '          '>>Circle 2 radius must be > Circle 1 radius',
     '          ERROR,*9999)

              K=CE(1,NEELEM(1,nr)) !setup the material properties
              SIGMAI=CE(2,NEELEM(1,nr))
              SIGMAE=K-SIGMAI
              SIGMAO=CE(1,NEELEM(1,2))

              ALPHA=((SIGMAI+SIGMAE)*RAD_2**2*(RAD_1**2*ANAL_C_1-
     '          ANAL_D_1))/(SIGMAO*(RAD_1**2-RAD_2**2)*ANAL_H_1)

              BETA=(2.0d0*SIGMAI*RAD_2**2)/((SIGMAI+SIGMAE)*(ALPHA*
     '          (RAD_1**2+RAD_2**2)*ANAL_H_1-RAD_2**2*
     '          (ANAL_C_1*RAD_1**2+ANAL_D_1)))

              ANAL_A_1=BETA*RHO_X
              ANAL_B_1=BETA*RHO_Y

              E_1=ALPHA*ANAL_A_1
              F_1=ALPHA*ANAL_B_1

              C_0=ANAL_REF-2.0d0*E_1*ANAL_H_1/RAD_2

              G_1=ANAL_H_1/RAD_2**2

              G_0=C_0

              nc=1

              DO nrr=1,2
                DO nonode=1,NPNODE(0,nrr)
                  np=NPNODE(nonode,nrr)
                  DO nk=1,
     &              MAX(NKH(NH_LOC(1,nx),np,nc,nrr)-KTYP93(nc,nrr),1)
                    ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nrr) !glob var num

                    RADIUS=DSQRT(XP(nk,nv,1,np)**2+XP(nk,nv,2,np)**2)

                    IF(XP(nk,nv,2,np).GT.ZERO_TOL) THEN ! positive theta
                      THETA=DACOS(XP(nk,nv,1,np)/RADIUS)
                    ELSE ! negative theta
                      THETA=-DACOS(XP(nk,nv,1,np)/RADIUS)
                    ENDIF

                    IF(nrr.EQ.1) THEN

                      YP(ny,1)=-(RAD_1**2+RADIUS**2)*(RHO_X*DCOS(THETA)+
     '                  RHO_Y*DSIN(THETA))/(RAD_1**2*RADIUS)
                      YP(ny,7)=(ANAL_A_1*DCOS(THETA)+ANAL_B_1*
     '                  DSIN(THETA))*(ANAL_C_1*RADIUS+
     '                  ANAL_D_1/RADIUS)+C_0-
     '                  SIGMAI*YP(ny,1)/(SIGMAI+SIGMAE)

                    ELSE IF(nrr.EQ.2) THEN

                      YP(ny,7)=(E_1*DCOS(THETA)+F_1*DSIN(THETA))*
     '                  (G_1*RADIUS+ANAL_H_1/RADIUS)+G_0

                    ENDIF
                  ENDDO !nk
                ENDDO !nonode
              ENDDO !nrr

              IF(nr.EQ.1) THEN
                nrr=2
              ELSE
                nrr=1
              ENDIF
              ANAL_CHOICE(nrr)=ANAL_CHOICE(nr)

C*** Derivatives
              nc=2
              DO nrr=1,2
                DO nonode=1,NPNODE(0,nrr)
                  IF(nonode.LE.NPNODE(0,1)) THEN
                    FACTOR=-1.0d0
                  ELSE
                    FACTOR=1.0d0
                  ENDIF
                  np=NPNODE(nonode,nrr)
                  DO nk=1,
     &              MAX(NKH(NH_LOC(1,nx),np,nc,nrr)-KTYP93(nc,nrr),1)
                    ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nrr) !glob var num

                    RADIUS=DSQRT(XP(nk,nv,1,np)**2 +XP(nk,nv,2,np)**2)
                    IF(XP(nk,nv,2,np).GT.ZERO_TOL) THEN ! positive theta
                      THETA=DACOS(XP(nk,nv,1,np)/RADIUS)
                    ELSE ! negative theta
                      THETA=-DACOS(XP(nk,nv,1,np)/RADIUS)
                    ENDIF

                    IF(nrr.EQ.1) THEN

                      YP(ny,7)=(SIGMAI+SIGMAE)*(ANAL_A_1*DCOS(THETA)+
     '                  ANAL_B_1*DSIN(THETA))*(ANAL_C_1-
     '                  ANAL_D_1/RADIUS**2)

                    ELSE IF(nrr.EQ.2) THEN

                      YP(ny,7)=FACTOR*SIGMAO*(E_1*DCOS(THETA)+F_1*
     '                  DSIN(THETA))*(G_1-ANAL_H_1/RADIUS**2)

                    ENDIF

                  ENDDO !nk
                ENDDO !nonode
              ENDDO !nrr
            ENDIF

          ELSE

C*** Spheres for testing transfer matrices

            IF(ANAL_CHOICE(nr).EQ.7) THEN

              IF(NDIPOLES(1).EQ.0) THEN
                ERROR='>>Define a dipole in region 1 first'
                GOTO 9999
              ELSE IF(NDIPOLES(1).GT.1) THEN
                ERROR='>>Must only have one dipole in region 1'
                GOTO 9999
              ELSE
                IF(DIPOLE_CEN(1,0,1,1).NE.0.0d0.OR.
     '            DIPOLE_CEN(2,0,1,1).NE.0.0d0.OR.
     '            DIPOLE_CEN(3,0,1,1).NE.0.0d0) THEN
                  ERROR='>>The dipole must be centred at the origin'
                  GOTO 9999
                ENDIF
              ENDIF
              RHO_X=DIPOLE_DIR(1,0,1,1)
              RHO_Y=DIPOLE_DIR(2,0,1,1)
              RHO_Z=DIPOLE_DIR(3,0,1,1)
              CALL ASSERT(CALL_MESH.AND.JTYP14.EQ.2,
     '          '>>Define a spherical mesh first',ERROR,*9999)
              CALL ASSERT(NSPHERES.EQ.2,
     '          '>>Must have only two spheres',ERROR,*9999)
              np=NPNODE(1,1)
              RAD_1=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
     '          XP(1,nv,3,np)**2)
              np=NPNODE(NPNODE(0,2),2)
              RAD_2=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
     '          XP(1,nv,3,np)**2)

              CALL ASSERT(DABS(RAD_1).GT.ZERO_TOL,
     '          '>>Sphere 1 has a zero radius',ERROR,*9999)
              CALL ASSERT(DABS(RAD_2).GT.ZERO_TOL,
     '          '>>Sphere 2 has a zero radius',ERROR,*9999)
              CALL ASSERT(DABS(RAD_2-RAD_1).GT.ZERO_TOL,
     '          '>>Sphere 2 radius must be > Sphere 1 radius',
     '          ERROR,*9999)

              K=CE(1,NEELEM(1,nr)) !setup the material properties
              SIGMAI=CE(2,NEELEM(1,nr))
              SIGMAE=K-SIGMAI
              SIGMAO=CE(1,NEELEM(2,2))

              ALPHA=((SIGMAI+SIGMAE)*RAD_2**3*(RAD_1**3*ANAL_C_11-
     '          2.0d0*ANAL_D_11))/(2.0d0*SIGMAO*(RAD_1**3-RAD_2**3)*
     '          ANAL_H_11)

              BETA=(-3.0d0*SIGMAI*RAD_2**3)/((SIGMAI+SIGMAE)*(ALPHA*
     '          (2.0d0*RAD_1**3+RAD_2**3)*ANAL_H_11-RAD_2**3*
     '          (ANAL_C_11*RAD_1**3+ANAL_D_11)))

              ANAL_A_11=BETA*RHO_X
              ANAL_B_11=BETA*RHO_Y

              FACTOR=((SIGMAI+SIGMAE)*RAD_1**3*ANAL_C_1-
     '          2.0d0*SIGMAI*RHO_Z)/(2.0d0*((2.0d0*SIGMAI+
     '          2.0d0*SIGMAE+SIGMAO)*RAD_1**3+(SIGMAI+SIGMAE-
     '          SIGMAO)*RAD_2**3))

              ANAL_D_1=3.0d0*(2.0d0*RAD_1**3+RAD_2**3)*FACTOR-
     '          RAD_1**3*ANAL_C_1+3.0d0*SIGMAI*RHO_Z/(SIGMAI+SIGMAE)

              C_0=ANAL_REF-9.0d0*RAD_2*FACTOR

              E_11=ALPHA*ANAL_A_11
              F_11=ALPHA*ANAL_B_11

              G_11=2.0d0*ANAL_H_11/RAD_2**3

              H_1=3.0d0*RAD_2**3*FACTOR

              G_1=2.0d0*H_1/RAD_2**3

              G_0=C_0

              nc=1

              DO nrr=1,2
                DO nonode=1,NPNODE(0,nrr)
                  np=NPNODE(nonode,nrr)
                  DO nk=1,
     &              MAX(NKH(NH_LOC(1,nx),np,nc,nrr)-KTYP93(nc,nrr),1)
                    ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nrr) !glob var num

                    RADIUS=DSQRT(XP(nk,nv,1,np)**2+XP(nk,nv,2,np)**2
     '                +XP(nk,nv,3,np)**2)

                    RAD=XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2
                    IF(DABS(RAD).LT.ZERO_TOL) THEN
                      THETA=0.0d0 !straight above or below
                    ELSE
                      IF(XP(nk,nv,2,np).GT.ZERO_TOL) THEN ! positive theta
                        THETA=DACOS(XP(nk,nv,1,np)/DSQRT(RAD))
                      ELSE ! negative theta
                        THETA=-DACOS(XP(nk,nv,1,np)/DSQRT(RAD))
                      ENDIF
                    ENDIF
                    IF(XP(nk,nv,3,np).GT.ZERO_TOL) THEN
                      PHI=PI/2.0d0-DASIN(XP(nk,nv,3,np)/RADIUS)
                    ELSE
                      PHI=PI/2.0d0+DABS(DASIN(XP(nk,nv,3,np)/RADIUS))
                    ENDIF

                    IF(nrr.EQ.1) THEN

                      YP(ny,1)=(RAD_1**3+2.0d0*RADIUS**3)*(
     '                  RHO_X*DCOS(THETA)*DSIN(PHI)+
     '                  RHO_Y*DSIN(THETA)*DSIN(PHI)+RHO_Z*
     '                  DCOS(PHI))/(RAD_1**3*RADIUS**2)
                      YP(ny,7)=(ANAL_A_11*DCOS(THETA)+ANAL_B_11*
     '                  DSIN(THETA))*(ANAL_C_11*RADIUS+
     '                  ANAL_D_11/RADIUS**2)*DSIN(PHI)+(ANAL_C_1*
     '                  RADIUS+ANAL_D_1/RADIUS**2)*DCOS(PHI)+C_0-
     '                  SIGMAI*YP(ny,1)/(SIGMAI+SIGMAE)

                    ELSE IF(nrr.EQ.2) THEN

                      YP(ny,7)=(E_11*DCOS(THETA)+F_11*DSIN(THETA))*
     '                  (G_11*RADIUS+ANAL_H_11/RADIUS**2)*DSIN(PHI)+
     '                  (G_1*RADIUS+H_1/RADIUS**2)*DCOS(PHI)+G_0

                    ENDIF
                  ENDDO !nk
                ENDDO !nonode
              ENDDO !nrr

              IF(nr.EQ.1) THEN
                nrr=2
              ELSE
                nrr=1
              ENDIF
              ANAL_CHOICE(nrr)=ANAL_CHOICE(nr)

            ELSE

              K=CE(1,NEELEM(1,nr)) !setup the material properties
              GI=CE(2,NEELEM(1,nr))
              GE=K-GI
              radius=1.D0 ! the inner radius
              nc=1 ! variables first
              nv=1 ! variables first



C*** Calcuate the transmembrane potential
              WRITE(OP_STRING,'(
     '          '' >>Calculating transmembrane for region '',I3)') nr
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)


              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nk=1,MAX(NKH(NH_LOC(1,nx),np,nc,nr)-KTYP93(nc,nr),1)
                  ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nr) !glob var num

C*** r
                  radius=DSQRT( XP(nk,nv,1,np)**2 +XP(nk,nv,2,np)**2
     '              +XP(nk,nv,3,np)**2 )
                  COORD(1)=radius
                  rad=XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2
C*** theta
                  IF(DABS(rad).LT.ZERO_TOL) THEN
                    COORD(2)=0.0d0 !straight above or below
                  ELSE
                    IF(XP(nk,nv,2,np).GT.ZERO_TOL) THEN ! positive theta
                      COORD(2)=DACOS(XP(nk,nv,1,np)/
     '                  DSQRT((XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2)) )
                    ELSE ! negative theta
                      COORD(2)=-DACOS(XP(nk,nv,1,np)/
     '                  DSQRT((XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2)) )
                    ENDIF
                  ENDIF
C*** phi
                  IF(XP(nk,nv,3,np).GT.ZERO_TOL) THEN
                    COORD(3)=PI/2.0d0-DASIN(XP(nk,nv,3,np)/radius)
                  ELSE
                    COORD(3)=PI/2.0d0+DABS(DASIN(XP(nk,nv,3,np)/radius))
                  ENDIF


C*** Analytic transmembrane values
                  IF(ANAL_CHOICE(nr).EQ.1) THEN !case1b
                    YP(ny,1)=ANALY_PHI_M(1.D0,2.D0,0.D0,
     '                coord(1),coord(2),coord(3))
                    YP(ny,7)=YP(ny,1)
                  ELSEIF(ANAL_CHOICE(nr).EQ.2) THEN !case3a
                    YP(ny,1)=ANALY_PHI_M(1.D0,2.D0,1.D0,
     '                coord(1),coord(2),coord(3))
                    YP(ny,7)=YP(ny,1)
                  ELSEIF(ANAL_CHOICE(nr).EQ.3) THEN !case2b
                    YP(ny,1)=ANALY_PHI_M(0.D0,0.D0,1.D0,
     '                coord(1),coord(2),coord(3))
                    YP(ny,7)=YP(ny,1)
                  ELSEIF(ANAL_CHOICE(nr).EQ.4) THEN !case1a
                    YP(ny,1)=ANALY_PHI_M(1.D0,2.D0,0.D0,
     '                coord(1),coord(2),coord(3))
                    YP(ny,7)=YP(ny,1)
                  ELSEIF(ANAL_CHOICE(nr).EQ.5) THEN !case2a
                    YP(ny,1)=ANALY_PHI_M(0.D0,0.D0,1.D0,
     '                coord(1),coord(2),coord(3))
                    YP(ny,7)=YP(ny,1)
                  ELSEIF(ANAL_CHOICE(nr).EQ.6) THEN !case3b
                    YP(ny,1)=ANALY_PHI_M(1.D0,2.D0,1.D0,
     '                coord(1),coord(2),coord(3))
                    YP(ny,7)=YP(ny,1)
                  ELSE
                    ERROR='Analytic transmembrane not implemented'
                    GOTO 9999
                  ENDIF
C               WRITE(*,*) 'Setting tranmem value: ',np,ny,YP(ny,1)
                ENDDO !nk
              ENDDO !nonode



C*** Calculate the extracellular analytic solution
              IF(nr.EQ.1) THEN
                nrr=2
              ELSE
                nrr=1
              ENDIF

C LKC - this is so that CHKSOL know that an analytic solution has been
C  set for the other region
              ANAL_CHOICE(nrr)=ANAL_CHOICE(nr)



              WRITE(OP_STRING,'(
     '          '' >>Calculating extracellular for region '',I3)') nrr
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

              DO nonode=1,NPNODE(0,nrr)
                np=NPNODE(nonode,nrr)
                DO nk=1,
     &            MAX(NKH(NH_LOC(1,nx),np,nc,nrr)-KTYP93(nc,nrr),1)
                  ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nrr) !glob var num

C*** r
                  radius=DSQRT( XP(nk,nv,1,np)**2 +XP(nk,nv,2,np)**2
     '              +XP(nk,nv,3,np)**2 )
                  COORD(1)=radius
C*** theta
                  rad=XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2
                  IF(DABS(rad).LT.ZERO_TOL) THEN
                    COORD(2)=0.0d0 !straight above or below
                  ELSE
                    IF(XP(nk,nv,2,np).GT.ZERO_TOL) THEN ! positive theta
                      COORD(2)=DACOS(XP(nk,nv,1,np)/
     '                  DSQRT((XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2)) )
                    ELSE ! negative theta
                      COORD(2)=-DACOS(XP(nk,nv,1,np)/
     '                  DSQRT((XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2)) )
                    ENDIF
                  ENDIF
C*** phi
                  IF(XP(nk,nv,3,np).GT.ZERO_TOL) THEN
                    COORD(3)=PI/2.0d0-DASIN(XP(nk,nv,3,np)/radius)
                  ELSE
                    COORD(3)=PI/2.0d0+DABS(DASIN(XP(nk,nv,3,np)/radius))
                  ENDIF


C*** Analytic values
                  IF(ANAL_CHOICE(nr).EQ.1) THEN !case1b
                    YP(ny,7)=(DCOS(COORD(2))+2.D0*DSIN(COORD(2)))
     '                *(0.19672D0*COORD(1) + 2.6557D0/COORD(1)**2)
     '                *DSIN(COORD(3))
                  ELSEIF(ANAL_CHOICE(nr).EQ.2) THEN !case3a
                    YP(ny,7)=(DCOS(COORD(2)) + 2.D0*DSIN(COORD(2)))
     '                *(0.52174D0*COORD(1) + 2.087D0/COORD(1)**2)
     '                *DSIN(COORD(3))
     '                + ( 2.D0*COORD(1)+ 8.D0/COORD(1)**2 )
     '                *DCOS(COORD(3))
                  ELSEIF(ANAL_CHOICE(nr).EQ.3) THEN !case2b
                    YP(ny,7)=(2.D0*COORD(1)+27.D0/COORD(1)**2)
     '                *DCOS(COORD(3))
     '                -9.D0
                  ELSEIF(ANAL_CHOICE(nr).EQ.4) THEN !case1a
                    YP(ny,7)=(DCOS(COORD(2))+2.D0*DSIN(COORD(2)))
     '                *(0.52174D0*COORD(1) + 2.087D0/COORD(1)**2)
     '                *DSIN(COORD(3))
                  ELSEIF(ANAL_CHOICE(nr).EQ.5) THEN !case2a
                    YP(ny,7)=(2.D0*COORD(1)+8.D0/COORD(1)**2)
     '                *DCOS(COORD(3))
     '                -6.D0
                  ELSEIF(ANAL_CHOICE(nr).EQ.6) THEN !case3b
                    YP(ny,7)=(DCOS(COORD(2))+2.D0*DSIN(COORD(2)))
     '                *(0.19672D0*COORD(1)+2.6557D0/COORD(1)**2)
     '                *DSIN(COORD(3))
     '                +(2.D0*COORD(1) + 27.D0/COORD(1)**2)
     '                *DCOS(COORD(3))
     '                -9.D0
                  ELSE
                    ERROR='Analytic extracellular soln not implemented'
                    GOTO 9999
                  ENDIF
C                WRITE(*,*) 'Setting torso value: ',np,ny,YP(ny,7)

                ENDDO !nk
              ENDDO !nonode

            ENDIF

C*** Derivatives NC=2
            nc=2 !normal derivatives (fluxes)


            IF(ANAL_CHOICE(nr).EQ.7) THEN

              K=CE(1,NEELEM(1,nr)) !setup the material properties
              SIGMAI=CE(2,NEELEM(1,nr))
              SIGMAE=K-SIGMAI
              SIGMAO=CE(1,NEELEM(2,2))
              DO nrr=1,2
                DO nonode=1,NPNODE(0,nrr)
                  IF(nonode.LE.NPNODE(0,1)) THEN
                    FACTOR=-1.0d0
                  ELSE
                    FACTOR=1.0d0
                  ENDIF
                  np=NPNODE(nonode,nrr)
                  DO nk=1,
     &              MAX(NKH(NH_LOC(1,nx),np,nc,nrr)-KTYP93(nc,nrr),1)
                    ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nrr) !glob var num

                    RADIUS=DSQRT(XP(nk,nv,1,np)**2 +XP(nk,nv,2,np)**2
     '                +XP(nk,nv,3,np)**2)

                    RAD=XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2
                    IF(DABS(RAD).LT.ZERO_TOL) THEN
                      THETA=0.0d0 !straight above or below
                    ELSE
                      IF(XP(nk,nv,2,np).GT.ZERO_TOL) THEN ! positive theta
                        THETA=DACOS(XP(nk,nv,1,np)/DSQRT(RAD))
                      ELSE ! negative theta
                        THETA=-DACOS(XP(nk,nv,1,np)/DSQRT(RAD))
                      ENDIF
                    ENDIF
                    IF(XP(nk,nv,3,np).GT.ZERO_TOL) THEN
                      PHI=PI/2.0d0-DASIN(XP(nk,nv,3,np)/RADIUS)
                    ELSE
                      PHI=PI/2.0d0+DABS(DASIN(XP(nk,nv,3,np)/RADIUS))
                    ENDIF

                    IF(nrr.EQ.1) THEN

                      YP(ny,7)=(SIGMAI+SIGMAE)*((ANAL_A_11*DCOS(THETA)+
     '                  ANAL_B_11*DSIN(THETA))*(ANAL_C_11-2.0d0*
     '                  ANAL_D_11/RADIUS**3)*DSIN(PHI)+(ANAL_C_1-
     '                  2.0d0*ANAL_D_1/RADIUS**3)*DCOS(PHI))

                    ELSE IF(nrr.EQ.2) THEN

                      YP(ny,7)=FACTOR*SIGMAO*((E_11*DCOS(THETA)+F_11*
     '                  DSIN(THETA))*(G_11-2.0d0*ANAL_H_11/RADIUS**3)*
     '                  DSIN(PHI)+(G_1-2.0d0*H_1/RADIUS**3)*DCOS(PHI))

                    ENDIF

                  ENDDO !nk
                ENDDO !nonode
              ENDDO !nrr

            ELSE

C*** assumes there are only 2 region - nr is the region where the
C***  transmembrane solution is so nrr is the region which we want to
C***  solve the extracelluar potential for.
              IF(nr.EQ.1) THEN
                nrr=2 !outer
              ELSE
                nrr=1 !outer
              ENDIF

              DO nonode=1,NPNODE(0,nrr)
                np=NPNODE(nonode,nrr)
                DO nk=1,
     &            MAX(NKH(NH_LOC(1,nx),np,nc,nrr)-KTYP93(nc,nrr),1)
                  ny=NYNP(nk,1,NH_LOC(1,nx),np,0,nc,nrr) !glob var num

C*** r
                  radius=DSQRT( XP(nk,nv,1,np)**2 +XP(nk,nv,2,np)**2
     '              +XP(nk,nv,3,np)**2 )
                  COORD(1)=radius
C*** theta
                  rad=XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2
                  IF(DABS(rad).LT.ZERO_TOL) THEN
                    COORD(2)=0.D0 !straight above or below
                  ELSE
                    IF(XP(nk,nv,2,np).GT.ZERO_TOL) THEN ! positive theta
                      COORD(2)=DACOS(XP(nk,nv,1,np)/
     '                  DSQRT((XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2)) )
                    ELSE ! negative theta
                      COORD(2)=-DACOS(XP(nk,nv,1,np)/
     '                  DSQRT((XP(nk,nv,1,np)**2+ XP(nk,nv,2,np)**2)) )
                    ENDIF
                  ENDIF
C*** phi
                  IF(XP(nk,nv,3,np).GT.ZERO_TOL) THEN
                    COORD(3)=PI/2.D0-DASIN(XP(nk,nv,3,np)/radius)
                  ELSE
                    COORD(3)=PI/2.D0+DABS(DASIN(XP(nk,nv,3,np)/radius))
                  ENDIF



C*** Analytic derivatives
                  IF(ANAL_CHOICE(nr).EQ.1) THEN !case1b
                    YP(ny,7)=-(DCOS(COORD(2))+2.D0*DSIN(COORD(2)))
     '                *(0.19672D0 + (-2.D0)*2.6557D0/COORD(1)**3)
     '                *DSIN(COORD(3))
     '                *CE(1,NEELEM(1,nrr))
                  ELSEIF(ANAL_CHOICE(nr).EQ.2) THEN !case3a
                    YP(ny,7)=-(DCOS(COORD(2))+2.D0*DSIN(COORD(2)))
     '                *(0.5217D0 + (-2.D0)*2.087D0/COORD(1)**3)
     '                *DSIN(COORD(3))
     '                -(2.D0+(-2.D0)*8.D0/COORD(1)**3)
     '                *DCOS(COORD(3))
     '                *CE(1,NEELEM(1,nrr))
                  ELSEIF(ANAL_CHOICE(nr).EQ.3) THEN !case2b
                    YP(ny,7)=-(2.D0+(-2.D0*27.D0/COORD(1)**3))
     '                *DCOS(COORD(3))
     '                *CE(1,NEELEM(1,nrr))
                  ELSEIF(ANAL_CHOICE(nr).EQ.4) THEN !case1a
                    YP(ny,7)=-(DCOS(COORD(2))+2.D0*DSIN(COORD(2)))
     '                *(0.52174D0+(-2.D0*2.087D0/COORD(1)**3))
     '                *DSIN(COORD(3))
                  ELSEIF(ANAL_CHOICE(nr).EQ.5) THEN !case2a
                    YP(ny,7)=-(2.D0+(-2.D0)*8.D0/COORD(1)**3)
     '                *DCOS(COORD(3))
     '                *CE(1,NEELEM(1,nrr))
                  ELSEIF(ANAL_CHOICE(nr).EQ.6) THEN !case3b
                    YP(ny,7)=
     '                -(
     '                (DCOS(COORD(2))+2.D0*DSIN(COORD(2)))
     '                *(0.19672D0+(-2.D0)*2.6557D0/COORD(1)**3)
     '                *DSIN(COORD(3))
     '                +(2.D0 + (-2.D0)*27.D0/COORD(1)**3)
     '                *DCOS(COORD(3))
     '                )
     '                *CE(1,NEELEM(1,nrr))
                  ELSE
                    ERROR='Analytic deriviatives not implemented'
                    GOTO 9999
                  ENDIF
C                WRITE(*,*) 'setting flux ',np,ny,YP(ny,7)
                ENDDO !nk
              ENDDO !nonode
            ENDIF
          ENDIF !NJT
        ELSE
          ERROR='>>Equation type not implemented'
          GOTO 9999
        ENDIF !IF Poisson equation


      ELSE
        ERROR='>>Analysis type not implemented'
        GOTO 9999
      ENDIF

      CALL EXITS('IPANA9')
      RETURN
 9999 CALL ERRORS('IPANA9',ERROR)
      CALL EXITS('IPANA9')
      RETURN 1
      END


