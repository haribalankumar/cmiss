      SUBROUTINE IPANA3(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,INP,
     '  NBJ,NDIPOLES,NEELEM,NENP,NKJE,NKH,NPF,
     '  NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVHP,NVJE,nx,NYNP,CE,
     '  DIPOLE_CEN,DIPOLE_DIR,SE,
     '  XA,XE,XP,YP,ERROR,*)

C#### Subroutine: IPANA3
C###  Description:
C###    IPANA3 inputs analytic functions for FE30 equations for
C###    region nr and problem nx.

      IMPLICIT NONE
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh03.cmn'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NDIPOLES(NRM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM),
     '  NKJE(NKM,NNM,NJM,NEM),NKH(NHM,NPM,NCM),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,ne,nh,nj,nj2,nj3,njj,njj2,
     '  MAXOPT,MAX_SERIES,nk,nk1,nk1r,nk2,nk3,NKS(8,8),nks1,nks2,
     '  NKST(8),nonode,nonp,NOQUES,np,nrr,nv,ny,ny1,ny2,ny3,ny4
      REAL*8 EXACT,DATAN_MOD,GRADF(3),NORMAL(3),R,R1,R2,
     '  RADIUS,RM,
     '  S11,S12,S22,SG,TANGENT(3,3),THETA,LAMBDA1,LAMBDA2
      CHARACTER CHAR5*5
      LOGICAL FOUND,FILEIP

      DATA NKST/1,2,2,4,2,4,4,8/
      DATA NKS/1,7*0,
     '         1,2,6*0,
     '         1,3,6*0,
     '         1,2,3,4,4*0,
     '         1,5,6*0,
     '         1,2,5,6,4*0,
     '         1,3,5,7,4*0,
     '         1,2,3,4,5,6,7,8/

      CALL ENTERS('IPANA3',*9999)

      NOQUES=0
      ICHAR=999
      FILEIP=.FALSE.
      IF(ITYP5(nr,nx).EQ.1) THEN
        IF(ITYP2(nr,nx).EQ.3) THEN !Laplaces equation
          IF(NJT.EQ.2) THEN
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=ANAL_CHOICE(nr)
            ENDIF
            IF(ITYP3(nr,nx).EQ.2.AND.NJ_LOC(njl_fibr,0,nr).GT.0)
     '        THEN
              FORMAT='('' Specify the analytic solution [1]:'''//
     '          '/''   (1) K(x-y)'''//
     '          '/''   (2) K(x^2-y^2)'''//
     '          '/''   (3) K(x^2+2xy-y^2)'''//
     '          '/''   (4) K(a.r^n.cos(n.t)+b.r^m.sin(m.t)+c.r.cos(t)+'
     '          //'d.r.sin(t)+e)'''//
     '          '/''   (5) Centre dipole (single circle)'''//
     '          '/''   (6) Centre dipole (mulitple circles)'''//
     '          '/''  *(7) Eccentric dipole (single circle)'''//
     '          '/''  *(8) Eccentric dipole (mulitple circles)'''//
     '          '/''   (9) Anisotropic annulus'''//
     '          '/''  (10) Anisotropic plate, f(z)=Ke^z'''//
     '         '/''  (11) Anisotropic plate, f(z)=a.sin(z)+b.cos(z)'''//
     '        '/$,''   '',I1)'
              MAXOPT=11
            ELSE
            FORMAT='('' Specify the analytic solution [1]:'''//
     '        '/''   (1) K(x-y)'''//
     '        '/''   (2) K(x^2-y^2)'''//
     '        '/''   (3) K(x^2+2xy-y^2)'''//
     '        '/''   (4) K(a.r^n.cos(n.t)+b.r^m.sin(m.t)+c.r.cos(t)+'
     '        //'d.r.sin(t)+e)'''//
     '        '/''   (5) Centre dipole (single circle)'''//
     '        '/''   (6) Centre dipole (mulitple circles)'''//
     '        '/''  *(7) Eccentric dipole (single circle)'''//
     '        '/''  *(8) Eccentric dipole (mulitple circles)'''//
     '        '/$,''   '',I1)'
              MAXOPT=6
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        MAXOPT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
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
                    nv=NVHP(1,np,1)
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
     '            IDEFLT,1,NPNODE(0,0),LDATA,LDEFLT,RDATA,RDEFLT,
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
                    nv=NVHP(1,nonode,1)
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
              ELSE
                CALL ASSERT(NDIPOLES(nr).EQ.0,
     '            '>>Must only have a dipole in the first region',
     '            ERROR,*9999)
              ENDIF
            ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN
              CALL ASSERT(NRT.EQ.1,'>>Define one region only',ERROR,
     '          *9999)
              CALL ASSERT(ITYP10(1).EQ.2,'>>Coordinates must be polar',
     '          ERROR,*9999)
              CALL ASSERT(CALL_MATE,
     '          '>>Define material properities first',ERROR,*9999)
              CALL ASSERT(CALL_MESH,'>>Define mesh first',ERROR,*9999)
              FORMAT='(/$,'' Enter the value of K [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANISO_K
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     &          ERROR,*9999)
              IF(IOTYPE.NE.3) ANISO_K=RDATA(1)
              FORMAT='(/$,'' Enter the value of m [2]: '',I2)'
              IDEFLT(1)=2
              IF(IOTYPE.EQ.3) IDATA(1)=ANISO_M
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) ANISO_M=IDATA(1)
              R1=MESH3_RAD(1)
              R2=MESH3_RAD(NSPHERES)
              RM=DBLE(ANISO_M)
            ELSE IF(ANAL_CHOICE(nr).EQ.10.OR.ANAL_CHOICE(nr).EQ.11) THEN
              CALL ASSERT(ITYP10(1).EQ.1,
     '          '>>Coordinates must be rectangular cartesian',
     '          ERROR,*9999)
              CALL ASSERT(CALL_MATE,
     '          '>>Define material properities first',ERROR,*9999)
              IF(ANAL_CHOICE(nr).EQ.10) THEN
                FORMAT='(/$,'' Enter the value of K [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANISO_K
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) ANISO_K=RDATA(1)
              ELSE
                FORMAT='(/$,'' Enter the value of a [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_A
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_A=RDATA(1)
                FORMAT='(/$,'' Enter the value of b [1.0]: '',D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=ANAL_B
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_B=RDATA(1)
              ENDIF
              FORMAT='(/$,'' Enter the value of L [2.0]: '',D12.4)'
              RDEFLT(1)=2.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANISO_L
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) ANISO_L=RDATA(1)
              FORMAT='(/$,'' Enter the value of H [1.0]: '',D12.4)'
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANISO_H
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) ANISO_H=RDATA(1)
              FORMAT='(/$,'' Enter the value of theta (deg) [0.0]: '','
     '          //'D12.4)'
              RDEFLT(1)=0.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=ANISO_THETA*180.0d0/PI
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) ANISO_THETA=RDATA(1)*PI/180.0d0
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
     '        '/''  (13) Centre dipole (mulitiple spheres)'''//
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
            ELSE IF(ANAL_CHOICE(nr).GE.12.AND.ANAL_CHOICE(nr).LE.15)
     '          THEN !Centre dipole solutions
C GMH 19/8/96 Move to adaptive determination of series limits
              IF(ANAL_CHOICE(nr).EQ.15) THEN
                FORMAT='(/$,'' Enter the # of series [10]: '',I1)'
                IDEFLT(1)=10
                MAX_SERIES=NSPHERECOEFFMX/2
                IF(IOTYPE.EQ.3) IDATA(1)=ANAL_N
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,1,MAX_SERIES,LDATA,LDEFLT,RDATA,
     '            RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) ANAL_N=IDATA(1)
              ENDIF
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
                  IF(ANAL_CHOICE(nr).LE.13) THEN
                    IF(DIPOLE_CEN(1,0,1,1).NE.0.0d0.OR.
     '                DIPOLE_CEN(2,0,1,1).NE.0.0d0.OR.
     '                DIPOLE_CEN(3,0,1,1).NE.0.0d0) THEN
                      ERROR='>>The dipole must be centred at the origin'
                      GOTO 9999
                    ENDIF
                  ENDIF !ANAL_CHOICE
                ENDIF !NDIPOLES
                ANAL_A=DIPOLE_DIR(1,0,1,1)
                ANAL_B=DIPOLE_DIR(2,0,1,1)
                ANAL_C=DIPOLE_DIR(3,0,1,1)
                IF(ITYP3(nr,nx).EQ.2) THEN !generalised Laplace
                  SIGMA(1)=CE(1,NEELEM(1,1))
                ELSE
                  SIGMA(1)=1.0d0
                ENDIF !ITYP3(nr,nx)
                IF(ANAL_CHOICE(nr).EQ.13.OR.ANAL_CHOICE(nr).EQ.15) THEN
                  CALL ASSERT(CALL_MESH.AND.JTYP14.EQ.2,
     '              '>>Define a spherical mesh first',ERROR,*9999)
                  IF(ANAL_CHOICE(nr).NE.15) THEN !we can handle 1 sphere
                    CALL ASSERT(NSPHERES.GT.1,
     '                '>>Must have more than one sphere',ERROR,*9999)
                  ENDIF !ANAL_CHOICE
                  DO nrr=2,NRT
                    IF(ITYP3(nrr,nx).EQ.2) THEN !generalised Laplace
                      SIGMA(nrr)=CE(1,NEELEM(1,nrr))
                    ELSE
                      SIGMA(nrr)=1.0d0
                    ENDIF !ITYP3
                  ENDDO !nrr
C***              Find the coefficients for the analytic solution
                  CALL DIPOLE_SOLVE(ANAL_CHOICE(nr),
     '              DIPOLE_CEN(1,0,1,1),DIPOLE_DIR(1,0,1,1),
     '              ERROR,*9999)
                ENDIF !ANAL_CHOICE
                IF(ANAL_CHOICE(nr).EQ.12.OR.ANAL_CHOICE(nr).EQ.14) THEN
                  IDEFLT(1)=NPNODE(1,1)
                ELSE IF(ANAL_CHOICE(nr).GE.13.OR.
     '              ANAL_CHOICE(nr).EQ.15) THEN
C***              Find the node number to fix
                  FOUND=.FALSE.
                  nonode=1
                  DO WHILE(.NOT.FOUND.AND.nonode.LE.NPNODE(0,NRT))
                    np=NPNODE(nonode,NRT)
                    nv=NVHP(1,nonode,1)
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
                ENDIF !ANAL_CHOICE
                WRITE(CHAR5,'(I5)') IDEFLT(1)
                FORMAT='(/$,'' Enter the outer sphere node number '
     '            //'to fix the potential at ['//CHAR5//']: '',I5)'
                IF(IOTYPE.EQ.3) IDATA(1)=ANAL_FIXEDNODE
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,1,NPNODE(0,0),LDATA,LDEFLT,RDATA,RDEFLT,
     '            RMIN,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  IF(ANAL_CHOICE(nr).EQ.12.OR.
     '              ANAL_CHOICE(nr).EQ.14) THEN
                    nrr=1
                  ELSE
                    nrr=NRT
                  ENDIF !ANAL_CHOICE
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
                    nv=NVHP(1,nonode,1)
                    R=DSQRT(XP(1,nv,1,IDATA(1))**2+
     '                XP(1,nv,2,IDATA(1))**2+XP(1,nv,3,IDATA(1))**2)
                    IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                MESH3_RAD(NSPHERES).GT.SPHERE_RAD_TOL) THEN
                      ERROR='>>Node is not on the outer sphere'
                      GOTO 9999
                    ENDIF
                  ENDIF
                  ANAL_FIXEDNODE=IDATA(1)
                ENDIF !IOTYPE
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
              ENDIF !nr
            ENDIF !ANAL_CHOICE
          ENDIF !NJT

          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO nv=1,NVHP(1,np,1)
              DO nk=1,MAX(NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr),1)
                ny=NYNP(nk,nv,NH_LOC(1,nx),np,0,1,nr) !global var number
                IF(NJT.EQ.2) THEN
                  IF(ANAL_CHOICE(nr).EQ.9) THEN
                    R=XP(1,1,1,np)
                    THETA=XP(1,1,2,np)
                    ne=NENP(np,1)
                    SG=DSQRT(CE(2,ne)/CE(1,ne))
                  ELSE IF(ANAL_CHOICE(nr).EQ.10.OR.
     '                ANAL_CHOICE(nr).EQ.11) THEN
                    ne=NENP(np,1)
                    S11=CE(1,ne)*DCOS(ANISO_THETA)*DCOS(ANISO_THETA)+
     '                CE(2,ne)*DSIN(ANISO_THETA)*DSIN(ANISO_THETA)
                    S12=(CE(1,ne)-CE(2,ne))*DSIN(ANISO_THETA)*
     '                DCOS(ANISO_THETA)
                    S22=CE(1,ne)*DSIN(ANISO_THETA)*DSIN(ANISO_THETA)+
     '                CE(2,ne)*DCOS(ANISO_THETA)*DCOS(ANISO_THETA)
                    LAMBDA1=S12/S22
                    LAMBDA2=DSQRT(S11*S22-S12*S12)/S22
                  ENDIF
                  IF(nk.EQ.1) THEN
                    IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                      EXACT=ANAL_K*(XP(1,nv,1,np)-XP(1,nv,2,np))
                    ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(x^2-y^2)
                      EXACT=ANAL_K*(XP(1,nv,1,np)**2-XP(1,nv,2,np)**2)
                    ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2+2xy-y^2)
                      EXACT=ANAL_K*(XP(1,nv,1,np)**2+2.0d0*
     '                  XP(1,nv,1,np)*XP(1,nv,2,np)-XP(1,nv,2,np)**2)
                    ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(r^n.cos(n.t)..
                      RADIUS=DSQRT(XP(1,nv,1,np)*XP(1,nv,1,np)+
     '                  XP(1,nv,2,np)*XP(1,nv,2,np))
                      THETA=DATAN_MOD(XP(1,nv,1,np),XP(1,nv,2,np))
                      EXACT=ANAL_K*(ANAL_A*RADIUS**ANAL_N*
     '                  DCOS(ANAL_N*THETA)+ANAL_B*RADIUS**ANAL_M*
     '                  DSIN(ANAL_M*THETA)+ANAL_C*RADIUS*DCOS(THETA)+
     '                  ANAL_D*RADIUS*DSIN(THETA)+ANAL_E)
                    ELSE IF(ANAL_CHOICE(nr).GE.5.AND.
     '                  ANAL_CHOICE(nr).LE.6) THEN !dipole solutions
                      ne=0
                      DO nonp=1,NENP(np,0)
                        IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                      ENDDO !nonp
                      CALL ASSERT(ne.NE.0,
     '                  '>>Could not find an element in the region',
     '                  ERROR,*9999)
                      CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                  DIPOLE_DIR_NTIME,1,NDIPOLES,np,NP_INTERFACE,
     '                  nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                  XP(1,1,1,np),
     '                  EXACT,XP(1,1,1,np),XP(1,1,2,np),XP(1,1,1,np),
     '                  ERROR,*9999)
                    ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN
                      EXACT=(ANISO_K*DCOS(RM*THETA)*
     '                  (R2**(-2.0d0*RM/SG)*R**(RM/SG)+R**(-RM/SG)))/
     '                  (R2**(-2.0d0*RM/SG)*R1**(RM/SG)+R1**(-RM/SG))
                    ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN
                      EXACT=2.0d0*ANISO_K*DEXP(XP(1,1,1,np))*
     '                  DEXP(-1.0d0*LAMBDA1*XP(1,1,2,np))*
     '                  DCOS(LAMBDA2*XP(1,1,2,np))
                    ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN
                      EXACT=2.0d0*((ANAL_A*DCOS(XP(1,1,1,np))+ANAL_B*
     '                  DSIN(XP(1,1,1,np)))*DCOS(LAMBDA1*XP(1,1,2,np))+
     '                  (ANAL_A*DSIN(XP(1,1,1,np))-ANAL_B*
     '                  DCOS(XP(1,1,1,np)))*DSIN(LAMBDA1*XP(1,1,2,np)))*
     '                  DCOSH(LAMBDA2*XP(1,1,2,np))
                    ELSE
                      ERROR='>>Analytic solution type not implemented'
                      GOTO 9999
                    ENDIF
                  ELSE IF(nk.EQ.2.OR.nk.EQ.3) THEN
                    IF(ANAL_CHOICE(nr).EQ.4) THEN
C                     Tangent vector for a circle in polar coordinates
                      RADIUS=DSQRT(XP(1,nv,1,np)*XP(1,nv,1,np)+
     '                  XP(1,nv,2,np)*XP(1,nv,2,np))
                      THETA=DATAN_MOD(XP(1,nv,1,np),XP(1,nv,2,np))
                      IF(CALL_MESH.AND.JTYP14.EQ.2.AND.NJT.EQ.2) THEN
C                       s is clockwise
                        TANGENT(1,1)=0.0d0
                        TANGENT(2,1)=-1.0d0
                      ENDIF
                    ELSE
C cpb 28/10/96 Get tangent from function
CC                     Assumes standard Hermite interpolatn in these
CC                     cases. Don't normalise the tangent length and
CC                     don't check to see if it is zero length.  It
CC                     will be zero length at some hermite simplex
CC                     nodes but this will not affect anything.
C                      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C                        TANGENT(nj,nk-1)=XP(nk,nv,nj,np)
C                      ENDDO !nj
                      CALL GET_TNVECTOR(IBT,IDO,INP,NBJ,NENP,NKJE,np,
     '                  NPF,NP_INTERFACE,NPNE,nr,NRE,NVJE,nx,
     '                  NORMAL,SE,TANGENT,XA,XE,XP,ERROR,*9999)
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
                      GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,1,np)+
     '                  XP(1,nv,2,np))
                      GRADF(2)=2.0d0*ANAL_K*(XP(1,nv,1,np)-
     '                  XP(1,nv,2,np))
                      EXACT=GRADF(1)*TANGENT(1,1)+GRADF(2)*TANGENT(2,1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(r^n.cos(n.t)..
                      GRADF(1)=ANAL_K*(ANAL_A*ANAL_N*RADIUS**(ANAL_N-1)*
     '                  DCOS(ANAL_N*THETA)+ANAL_B*ANAL_M*
     '                  RADIUS**(ANAL_M-1)*DSIN(ANAL_M*THETA)+ANAL_C*
     '                  DCOS(THETA)+ANAL_D*DSIN(THETA))
                      GRADF(2)=ANAL_K*(-ANAL_A*ANAL_N*RADIUS**
     '                  (ANAL_N-1)*DSIN(ANAL_N*THETA)+ANAL_B*ANAL_M*
     '                  RADIUS**(ANAL_M-1)*DCOS(ANAL_M*THETA)-ANAL_C*
     '                  DSIN(THETA)+ANAL_D*DCOS(THETA))
                      EXACT=GRADF(1)*TANGENT(1,1)+GRADF(2)*TANGENT(2,1)
                    ELSE IF(ANAL_CHOICE(nr).GE.5.AND.
     '                  ANAL_CHOICE(nr).LE.6) THEN !Dipole solutions
                      ne=0
                      DO nonp=1,NENP(np,0)
                        IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                      ENDDO !nonp
                      CALL ASSERT(ne.NE.0,
     '                  '>>Could not find an element in the region',
     '                  ERROR,*9999)
                      CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                  DIPOLE_DIR_NTIME,2,NDIPOLES,np,NP_INTERFACE,
     '                  nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                  XP(1,1,1,np),
     '                  EXACT,XP(1,1,1,np),XP(1,1,2,np),XP(1,1,1,np),
     '                  ERROR,*9999)
                    ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN
                      IF(nk.EQ.2) THEN
                        EXACT=(-RM*ANISO_K*DSIN(RM*THETA)*
     '                    (R2**(-2.0d0*RM/SG)*R**(RM/SG)+R**(-RM/SG)))/
     '                    (R2**(-2.0d0*RM/SG)*R1**(RM/SG)+R1**(-RM/SG))
                      ELSE IF(nk.EQ.3) THEN
                        EXACT=(RM*ANISO_K*DCOS(RM*THETA)*(R2**(-2.0d0*
     '                    RM/SG)*R**((RM-1.0d0)/SG)+R**
     '                    ((-RM+1.0d0)/SG)))/(SG*(R2**(-2.0d0*RM/SG)*
     '                    R1**(RM/SG)+R1**(-RM/SG)))
                      ENDIF
                    ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN
                      GRADF(1)=2.0d0*ANISO_K*DEXP(XP(1,1,1,np))*
     '                  DEXP(-1.0d0*LAMBDA1*XP(1,1,2,np))*
     '                  DCOS(LAMBDA2*XP(1,1,2,np))
                      GRADF(2)=-2.0d0*ANISO_K*LAMBDA1*
     '                  DEXP(XP(1,1,1,np))*DEXP(-1.0d0*LAMBDA1*
     '                  XP(1,1,2,np))*DCOS(LAMBDA2*XP(1,1,2,np))
     '                  -2.0d0*ANISO_K*LAMBDA2*DEXP(XP(1,1,1,np))*
     '                  DEXP(-1.0d0*LAMBDA1*XP(1,1,2,np))*
     '                  DSIN(LAMBDA2*XP(1,1,2,np))
                      EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN
                      GRADF(1)=2.0d0*((-ANAL_A*DSIN(XP(1,1,1,np))+
     '                  ANAL_B*DCOS(XP(1,1,1,np)))*DCOS(LAMBDA1*
     '                  XP(1,1,2,np))+(ANAL_A*DCOS(XP(1,1,1,np))+ANAL_B*
     '                  DSIN(XP(1,1,1,np)))*DSIN(LAMBDA1*XP(1,1,2,np)))*
     '                  DCOSH(LAMBDA2*XP(1,1,2,np))
                      GRADF(2)=2.0d0*LAMBDA1*((ANAL_A*
     '                  DCOS(XP(1,1,1,np))+ANAL_B*DSIN(XP(1,1,1,np)))*
     '                  (-DSIN(LAMBDA1*XP(1,1,2,np)))+(ANAL_A*
     '                  DSIN(XP(1,1,1,np))-ANAL_B*DCOS(XP(1,1,1,np)))*
     '                  DCOS(LAMBDA1*XP(1,1,2,np)))*
     '                  DCOSH(LAMBDA2*XP(1,1,2,np))+
     '                  2.0d0*LAMBDA2*((ANAL_A*DCOS(XP(1,1,1,np))+
     '                  ANAL_B*DSIN(XP(1,1,1,np)))*DCOS(LAMBDA1*
     '                  XP(1,1,2,np))+(ANAL_A*DSIN(XP(1,1,1,np))-
     '                  ANAL_B*DCOS(XP(1,1,1,np)))*DSIN(LAMBDA1*
     '                  XP(1,1,2,np)))*DSINH(LAMBDA2*XP(1,1,2,np))
                      EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)
                    ELSE
                      ERROR='>>Analytic solution type not implemented'
                      GOTO 9999
                    ENDIF
                  ELSE IF(nk.EQ.4) THEN
                    IF(ANAL_CHOICE(nr).EQ.9) THEN
                      EXACT=(-RM*RM*ANISO_K*DSIN(RM*THETA)*
     '                  (R2**(-2.0d0*RM/SG)*R**((RM-1.0d0)/SG)+
     '                  R**((-RM+1.0d0)/SG)))/(SG*
     '                  (R2**(-2.0d0*RM/SG)*R1**(RM/SG)+R1**(-RM/SG)))
                    ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN
C !!! cpb 28/10/96 Only code for elements aligned with the xy axes for
C !!! the moment.
                      EXACT=-2.0d0*ANISO_K*LAMBDA1*
     '                  DEXP(XP(1,1,1,np))*DEXP(-1.0d0*LAMBDA1*
     '                  XP(1,1,2,np))*DCOS(LAMBDA2*XP(1,1,2,np))
     '                  -2.0d0*ANISO_K*LAMBDA2*DEXP(XP(1,1,1,np))*
     '                  DEXP(-1.0d0*LAMBDA1*XP(1,1,2,np))*
     '                  DSIN(LAMBDA2*XP(1,1,2,np))
                    ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN
C !!! cpb 28/10/96 Only code for elements aligned with the xy axes for
C !!! the moment.
                      EXACT=2.0d0*LAMBDA1*((-ANAL_A*
     '                  DSIN(XP(1,1,1,np))+ANAL_B*DCOS(XP(1,1,1,np)))*
     '                  (-DSIN(LAMBDA1*XP(1,1,2,np)))+(ANAL_A*
     '                  DCOS(XP(1,1,1,np))+ANAL_B*DSIN(XP(1,1,1,np)))*
     '                  DCOS(LAMBDA1*XP(1,1,2,np)))*
     '                  DCOSH(LAMBDA2*XP(1,1,2,np))+
     '                  2.0d0*LAMBDA2*((-ANAL_A*DSIN(XP(1,1,1,np))+
     '                  ANAL_B*DCOS(XP(1,1,1,np)))*DCOS(LAMBDA1*
     '                  XP(1,1,2,np))+(ANAL_A*DCOS(XP(1,1,1,np))-
     '                  ANAL_B*DSIN(XP(1,1,1,np)))*DSIN(LAMBDA1*
     '                  XP(1,1,2,np)))*DSINH(LAMBDA2*XP(1,1,2,np))
                    ENDIF
                  ENDIF !nk
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
     '                  2.0d0*XP(1,nv,3,np)**2)
                    ELSE IF(ANAL_CHOICE(nr).EQ.7) THEN !K(x^2-2y^2+z^2)
                      EXACT=ANAL_K*(XP(1,nv,1,np)**2-2.0d0*
     '                  XP(1,nv,2,np)**2+XP(1,nv,3,np)**2)
                    ELSE IF(ANAL_CHOICE(nr).EQ.8) THEN !K(-2x^2+y^2+z^2)
                      EXACT=ANAL_K*(-2.0d0*XP(1,nv,1,np)**2+
     '                  XP(1,nv,2,np)**2+XP(1,nv,3,np)**2)
                    ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN !K(x^2+2xy-y^2)
                      EXACT=ANAL_K*(XP(1,nv,1,np)**2+2.0d0*
     '                  XP(1,nv,1,np)*XP(1,nv,2,np)-XP(1,nv,2,np)**2)
                    ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN !K(x^2+2xz-z^2)
                      EXACT=ANAL_K*(XP(1,nv,1,np)**2+2.0d0*
     '                  XP(1,nv,1,np)*XP(1,nv,3,np)-XP(1,nv,3,np)**2)
                    ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN !K(y^2+2yz-z^2)
                      EXACT=ANAL_K*(XP(1,nv,2,np)**2+2.0d0*
     '                  XP(1,nv,2,np)*XP(1,nv,3,np)-XP(1,nv,3,np)**2)
                    ELSE IF(ANAL_CHOICE(nr).GE.12.AND.
     '                  ANAL_CHOICE(nr).LE.15) THEN
C                     Dipole solutions
                      ne=0
                      DO nonp=1,NENP(np,0)
                        IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                      ENDDO !nonp
                      CALL ASSERT(ne.NE.0,
     '                  '>>Could not find an element in the region',
     '                  ERROR,*9999)
                      CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                  DIPOLE_DIR_NTIME,1,NDIPOLES,np,NP_INTERFACE,
     '                  nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                  XP(1,1,1,np),EXACT,XP(1,1,1,np),
     '                  XP(1,1,2,np),XP(1,1,3,np),ERROR,*9999)
                    ELSE
                      ERROR='>>Analytic solution type not implemented'
                      GOTO 9999
                    ENDIF
                  ELSE IF(nk.EQ.2.OR.nk.EQ.3) THEN
C cpb 28/10/96 Get tangent from function
C                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C                      TANGENT(nj,nk-1)=XP(nk,nv,nj,np)
C                    ENDDO
                    CALL GET_TNVECTOR(IBT,IDO,INP,NBJ,NENP,NKJE,np,
     '                NPF,NP_INTERFACE,NPNE,nr,NRE,NVJE,nx,
     '                NORMAL,SE,TANGENT,XA,XE,XP,ERROR,*9999)
                    IF(ANAL_CHOICE(nr).EQ.1) THEN !K(x-y)
                      GRADF(1)=ANAL_K
                      GRADF(2)=-ANAL_K
                      EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN !K(z)
                      GRADF(3)=ANAL_K
                      EXACT=GRADF(3)*TANGENT(3,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN !K(x^2-y^2)
                      GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                      GRADF(2)=-2.0d0*ANAL_K*XP(1,nv,2,np)
                      EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN !K(x^2-z^2)
                      GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                      GRADF(3)=-2.0d0*ANAL_K*XP(1,nv,3,np)
                      EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(3)*
     '                  TANGENT(3,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.5) THEN !K(y^2-z^2)
                      GRADF(2)=2.0d0*ANAL_K*XP(1,nv,2,np)
                      GRADF(3)=-2.0d0*ANAL_K*XP(1,nv,3,np)
                      EXACT=GRADF(2)*TANGENT(2,nk-1)+GRADF(3)*
     '                  TANGENT(3,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !K(x^2+y^2-2z^2)
                      GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                      GRADF(2)=2.0d0*ANAL_K*XP(1,nv,2,np)
                      GRADF(3)=-4.0d0*ANAL_K*XP(1,nv,3,np)
                      EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)+GRADF(3)*TANGENT(3,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.7) THEN !K(x^2-2y^2+z^2)
                      GRADF(1)=2.0d0*ANAL_K*XP(1,nv,1,np)
                      GRADF(2)=-4.0d0*ANAL_K*XP(1,nv,2,np)
                      GRADF(3)=2.0d0*ANAL_K*XP(1,nv,3,np)
                      EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)+GRADF(3)*TANGENT(3,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.8) THEN !K(-2x^2+y^2+z^2)
                      GRADF(1)=-4.0d0*ANAL_K*XP(1,nv,1,np)
                      GRADF(2)=2.0d0*ANAL_K*XP(1,nv,2,np)
                      GRADF(3)=2.0d0*ANAL_K*XP(1,nv,3,np)
                      EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)+GRADF(3)*TANGENT(3,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN !K(x^2+2xy-y^2)
                      GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,1,np)+
     '                  XP(1,nv,2,np))
                      GRADF(2)=2.0d0*ANAL_K*(XP(1,nv,1,np)-
     '                  XP(1,nv,2,np))
                      EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(2)*
     '                  TANGENT(2,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN !K(x^2+2xz-z^2)
                      GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,1,np)+
     '                  XP(1,nv,3,np))
                      GRADF(3)=2.0d0*ANAL_K*(XP(1,nv,1,np)-
     '                  XP(1,nv,3,np))
                      EXACT=GRADF(1)*TANGENT(1,nk-1)+GRADF(3)*
     '                  TANGENT(3,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN !K(y^2+2yz-z^2)
                      GRADF(1)=2.0d0*ANAL_K*(XP(1,nv,2,np)+
     '                  XP(1,nv,3,np))
                      GRADF(2)=2.0d0*ANAL_K*(XP(1,nv,2,np)-
     '                  XP(1,nv,3,np))
                      EXACT=GRADF(2)*TANGENT(2,nk-1)+GRADF(3)*
     '                  TANGENT(3,nk-1)
                    ELSE IF(ANAL_CHOICE(nr).GE.12.AND
     '                  .ANAL_CHOICE(nr).LE.15) THEN
C                     Dipole solutions
                      ne=0
                      DO nonp=1,NENP(np,0)
                        IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                      ENDDO !nonp
                      CALL ASSERT(ne.NE.0,
     '                  '>>Could not find an element in the region',
     '                  ERROR,*9999)
                      CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                  DIPOLE_DIR_NTIME,nk,NDIPOLES,np,NP_INTERFACE,
     '                  nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                  XP(1,1,1,np),EXACT,
     '                  XP(1,1,1,np),XP(1,1,2,np),XP(1,1,3,np),
     '                  ERROR,*9999)
                    ELSE
                      ERROR='>>Analytic solution type not implemented'
                      GOTO 9999
                    ENDIF
                  ELSE IF(nk.EQ.4) THEN
                    IF(ANAL_CHOICE(nr).GE.12.AND.
     '                ANAL_CHOICE(nr).LE.15) THEN
                      ne=0
                      DO nonp=1,NENP(np,0)
                        IF(NRE(NENP(np,nonp)).EQ.nr) ne=NENP(np,nonp)
                      ENDDO !nonp
                      CALL ASSERT(ne.NE.0,
     '                  '>>Could not find an element in the region',
     '                  ERROR,*9999)
                      CALL DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,
     '                  DIPOLE_DIR_NTIME,7,NDIPOLES,np,NP_INTERFACE,
     '                  nr,nx,CE(1,ne),DIPOLE_CEN,DIPOLE_DIR,
     '                  XP(1,1,1,np),EXACT,
     '                  XP(1,1,1,np),XP(1,1,2,np),XP(1,1,3,np),
     '                  ERROR,*9999)
                    ENDIF
                  ENDIF
                ENDIF !NJT
C cpb 30/10/96 Now taken care of in the get_tnvector function.
C                IF((NP_INTERFACE(np,0).GT.1).AND.(nr.NE.
C     '            NP_INTERFACE(np,1))) THEN
CC                 Interface node and not in the first region.
C                  FEMFEMCOUP=ITYP4(nr,nx).EQ.1.AND.
C     '              ITYP4(NP_INTERFACE(np,1),nx).EQ.1
C                  IF(.NOT.FEMFEMCOUP) THEN
C                    IF(NJT.EQ.2) THEN
C                      IF(nk.GE.2) THEN
CC                       Reverse sign of s1 derivatives at interface
C                        EXACT=-EXACT
C                      ENDIF
C                    ELSE IF(NJT.EQ.3) THEN
C                      IF(nk.GE.3) THEN
CC                       Reverse sign of s2 and s1s2 derivatives at
CC                       interface
C                        EXACT=-EXACT
C                      ENDIF
C                    ENDIF
C                  ENDIF
C                ENDIF
                YP(ny,7)=EXACT
              ENDDO !nk
            ENDDO !nv
          ENDDO !nonode (np)

        ELSE
          ERROR='>>Not implemented'
          GOTO 9999
        ENDIF
      ELSE IF(ITYP5(nr,nx).EQ.5.AND.NJT.EQ.3) THEN !wavefront propagation
C       These are just to test the diffusion term.
        FORMAT='('' Enter coefficients of x, y, z [0 0 0]: ''/'//
     '    '$,''    '',3D12.4)'
        IF(IOTYPE.EQ.3) THEN
          DO njj=1,3
            RDATA(njj)=ANAL_P1K(njj)
          ENDDO !njj
        ENDIF !IOTYPE=3
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO njj=1,3
            ANAL_P1K(njj)=RDATA(njj)
          ENDDO !njj
        ENDIF !IOTYPE!=3
        FORMAT='('' Enter coefficients of yz, zx, xy [0 0 0]: ''/'//
     '    '$,''    '',3D12.4)'
        IF(IOTYPE.EQ.3) THEN
          DO njj=1,3
            RDATA(njj)=ANAL_P2AK(njj)
          ENDDO !njj
        ENDIF !IOTYPE=3
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO njj=1,3
            ANAL_P2AK(njj)=RDATA(njj)
          ENDDO !njj
        ENDIF !IOTYPE!=3
        FORMAT=
     '    '('' Enter coefficients of y^2-z^2, z^2-x^2, x^2-y^2 '//
     '    '[0 0 0]: ''/$,''    '',3D12.4)'
        IF(IOTYPE.EQ.3) THEN
          DO njj=1,3
            RDATA(njj)=ANAL_P2BK(njj)
          ENDDO !njj
        ENDIF !IOTYPE=3
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO njj=1,3
            ANAL_P2BK(njj)=RDATA(njj)
          ENDDO !njj
        ENDIF !IOTYPE!=3
        FORMAT='('' Enter coefficient of xyz [0]: ''/'//
     '    '$,''    '',3D12.4)'
        IF(IOTYPE.EQ.3) THEN
          RDATA(1)=ANAL_P3K
        ENDIF !IOTYPE=3
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          ANAL_P3K=RDATA(1)
        ENDIF !IOTYPE!=3
        nh=NH_LOC(1,nx)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nv=1,NVHP(nh,np,1)
            DO nk=1,MAX(NKH(nh,np,1)-KTYP93(1,nr),1)

              EXACT=0.0d0
C KAT 22May99: Not very efficient but the code is short and easy to change.
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)

                nj=NJ_LOC(NJL_GEOM,njj,nr)
                njj2=1+MOD(njj,3)
                nj2=NJ_LOC(NJL_GEOM,njj2,nr)
                njj2=1+MOD(njj2,3)
                nj3=NJ_LOC(NJL_GEOM,njj2,nr)

                EXACT=EXACT+ANAL_P1K(njj)*XP(nk,nv,nj,np)

                DO nks1=1,NKST(nk)
                  nk1=NKS(nks1,nk) !subderivs
                  nk2=nk-nk1+1 !rest of derivs

                  EXACT=EXACT
     '              +ANAL_P2AK(njj)
     '                *XP(nk1,nv,nj2,np)*XP(nk2,nv,nj3,np)
     '              +ANAL_P2BK(njj)
     '                *(XP(nk1,nv,nj2,np)*XP(nk2,nv,nj2,np)
     '                  -XP(nk1,nv,nj3,np)*XP(nk2,nv,nj3,np))

                ENDDO !nk1

              ENDDO !nj

              DO nks1=1,NKST(nk)
                nk1=NKS(nks1,nk) !subderivs
                nk1r=nk-nk1+1 !rest of derivs

                DO nks2=1,NKST(nk1r)
                  nk2=NKS(nks2,nk1r) !subderivs
                  nk3=nk1r-nk2+1 !rest of derivs

                  EXACT=EXACT
     '              +ANAL_P3K*XP(nk1,nv,nj,np)
     '                *XP(nk2,nv,nj2,np)*XP(nk3,nv,nj3,np)

                ENDDO !nk2
              ENDDO !nk1

              ny=NYNP(nk,nv,nh,np,0,1,nr) !global var number
              YP(ny,1)=EXACT
            ENDDO !nk
          ENDDO !nv
        ENDDO !nonode (np)

      ELSE IF(ITYP5(nr,nx).EQ.2) THEN !time integration
C Time integration advection-diffusion analytic solutions
C DMAL 22-MAY-2002
        IF (ITYP2(nr,nx).EQ.3) THEN !advection-diffusion
c            FORMAT='('' Specify the analytic solution [1]:'''//
c     '        '/''   (1) 2D diffusion, rectangular plate,'''//
c     '        '/''       BCs: T(0,y,t)=T(x,0,t)=T(x,b,t)=0'''//
c     '        '/''       ICs: T(a,y,t)=0'''//
c     '        '/''   (2) 1D diffusion, rod of length, L,'''//
c     '        '/''       BCs: T(0,t)=T1, T(L,t)=T2,'''//
c     '        '/''       ICs: T(x,0)=sin(PI*x/L)'''//
c     '        '/$,''   '',I1)'
c              MAXOPT=2
              FORMAT='('' Specify the analytic solution [1]:'''//
     '        '/''   (1) 2D rectangular plate of size, axb with '//
     '           'T(0,y,t)=T(x,0,t)=T(x,b,t)=0 and T(a,y,t)='//
     '           'sin(pi*y/b)'''//
     '        '/''   (2) 1D rod of length, L with T(0,t)=T1, '//
     '          'T(L,t)=T2 and T(x,0)=A'''//
     '        '/$,''   '',I2)'
            MAXOPT=2
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        MAXOPT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) ANAL_CHOICE(nr)=IDATA(1)


C          IF(IOTYPE.EQ.3) RDATA(1)=ANAL_CHOICE(nr)
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,MAXOPT,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) ANAL_CHOICE(nr)=IDATA(1)

          IF(ANAL_CHOICE(nr).EQ.1) THEN
            FORMAT='(/$,'' Enter the width, a (xi=1) of the plate '//
     '               '[1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=ANAL_A
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_A=RDATA(1)
            FORMAT='(/$,'' Enter the height, b (xi=2) of the plate '//
     '               '[1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=ANAL_B
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_B=RDATA(1)
            FORMAT='(/$,'' Enter the diffusion constant, D '//
     '               '[1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=ANAL_D
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_D=RDATA(1)
            FORMAT='(/$,'' Enter the time '//
     '               '[1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=ANAL_TIME
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_TIME=RDATA(1)

C calculating analytic solution at nodes
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              ny1=NYNP(1,1,1,np,0,1,nr)
              ny2=NYNP(2,1,1,np,0,1,nr)
              ny3=NYNP(3,1,1,np,0,1,nr)
              ny4=NYNP(4,1,1,np,0,1,nr)
              YP(ny1,7)=0.0d0
              YP(ny2,7)=0.0d0
              YP(ny3,7)=0.0d0
              YP(ny4,7)=0.0d0
              CALL ANS001(YP(ny4,7),YP(ny2,7),YP(ny3,7),
     '          YP(ny1,7),ANAL_TIME,XP(1,1,1,np),XP(1,1,2,np),ERROR,
     '          *9999)
            ENDDO

          ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN
            FORMAT='(/$,'' Enter the length, L (xi=1) of the rod '//
     '               '[1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=ANAL_A
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_A=RDATA(1)
            FORMAT='(/$,'' Enter the temperature, T1 at x=0 '//
     '               '[1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=ANAL_B
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_B=RDATA(1)
            FORMAT='(/$,'' Enter the temperature, T2 at x=L '//
     '               '[1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=ANAL_C
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_C=RDATA(1)
            FORMAT='(/$,'' Enter the diffusion constant, D '//
     '               '[1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=ANAL_D
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_D=RDATA(1)
            FORMAT='(/$,'' Enter the time '//
     '               '[1.0]: '',D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=ANAL_TIME
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ANAL_TIME=RDATA(1)

C calculating analytic solution at nodes
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              ny1=NYNP(1,1,1,np,0,1,nr)
              ny2=NYNP(2,1,1,np,0,1,nr)
              YP(ny1,7)=0.0d0
              YP(ny2,7)=0.0d0
              CALL ANS002(YP(ny2,7),YP(ny1,7),ANAL_TIME,
     '          XP(1,1,1,np),ERROR,*9999)
            ENDDO

          ELSE
            ERROR='>>Not implemented'
            GOTO 9999
          ENDIF
        ELSE
          ERROR='>>Not implemented'
          GOTO 9999
        ENDIF
      ELSE
        ERROR='>>Not implemented'
        GOTO 9999
      ENDIF

      CALL EXITS('IPANA3')
      RETURN
 9999 CALL ERRORS('IPANA3',ERROR)
      CALL EXITS('IPANA3')
      RETURN 1
      END



