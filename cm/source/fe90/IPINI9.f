      SUBROUTINE IPINI9(IBT,NBH,NEELEM,NHP,NKH,NPLIST,NP_INTERFACE,
     '  NPNODE,NPNY,nr,NVHP,NW,nx,NYNP,NYNR,XP,YP,ACTIVATION,
     '  ALL_REGIONS,FIX,FIX_ZERO,GENER,NOFIX,ERROR,*)

C#### Subroutine: IPINI9
C###  Description:
C###    IPINI9 inputs initial conditions and boundary conditions for
C###    BE problems.

C**** New version which takes no account of mixed boundary conditions
C**** and includes cornering using the new nv variable.  For mixed bc
C**** see the old version of this subroutine.
C**** If node np lies at a corner then the user must specify
C**** (for elliptic pdes) either
C****  1.  All corner values of du/dn (2 for 2d, 3 for 3d)   OR
C****  2.  Value of dependent variable and all but one value of du/dn
C****  Gener is .true. means that the initial conditions are
C****  generated according to some test problem.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'anal00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh03.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHP(NPM),
     '  NKH(NHM,NPM,NCM),NPLIST(0:NPM),
     '  NP_INTERFACE(0:NPM,0:3),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),nr,
     '  NVHP(NHM,NPM,NCM),NW(NEM,3),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
CC AJPs - 191297 - rgb
      LOGICAL ACTIVATION,ALL_REGIONS,FIX(NYM,NIYFIXM),FIX_ZERO,GENER,
     '  NOFIX
CC AJPe
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,iy,loop,n,N1,nc,nh,nhx,nk,nonode,no_nynr,
     '  NOQUES,np,nrr,nv,ny,ny1,ny2,ny3,ny4,ny_first,NYTOT
      REAL*8 R
      LOGICAL CORNER,FILEIP,INLIST
      CHARACTER CHAR1*1

      CALL ENTERS('IPINI9',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(iotype.NE.3) THEN
C***    Need to only initialise those belonging to the
C***    current region.
        DO nc=1,NCT(nr,nx)
          DO no_nynr=1,NYNR(0,0,nc) ! loop over global variables
            ny=NYNR(no_nynr,0,nc) ! global variable number
            IF(NPNY(0,ny,0,nx).EQ.1) THEN
              np=NPNY(4,ny,0,nx)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            ENDIF
            DO iy=1,NIYM
              IF(iy.NE.7) THEN !iy=7 is analytic
                YP(ny,iy)=0.0d0
              ENDIF
            ENDDO !iy
            DO iy=1,NIYFIXM
              FIX(ny,iy)=.FALSE.
            ENDDO !iy
          ENDDO !no_nynr (ny)
        ENDDO !nc
      ENDIF

      IF(GENER.AND..NOT.ACTIVATION) THEN

        IF(CALL_ANAL_TRSF.AND.NJT.EQ.3.AND.ANAL_CHOICE(nr).EQ.7) THEN

          nrr=2

          nc=1
          np=NPNODE(NPNODE(0,nrr),nrr)
          ny1=NYNP(1,1,NH_LOC(1,nx),np,0,nc,nrr)
          YP(ny1,1)=ANAL_REF
          FIX(ny1,1)=.TRUE.

          nc=2
          DO nonode=NPNODE(0,1)+1,NPNODE(0,nrr)
            np=NPNODE(nonode,nrr)
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,nc)
                DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nrr),1)
                  ny1=NYNP(nk,nv,nh,np,0,nc,nrr)
                  YP(ny1,1)=0.0d0
                  FIX(ny1,1)=.TRUE.
                ENDDO !nk
              ENDDO !nv
            ENDDO !nhx
          ENDDO !nonode

        ELSE
C***      If we are dealing with more than one region then
C***      only set up the initial conditions for the nodes
C***      in the first and last regions that don't belong to
C***      an interface
C***      Need to check if hermite interpolation is used
          CALL EQTYPE(IBT,NBH,NEELEM,nr,NW,nx,ERROR,*9999)
          DO nc=1,NCT(nr,nx)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              DO nhx=1,NHP(np)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,nc)
                  DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                    ny1=NYNP(nk,nv,nh,np,0,nc,nr)
                    IF(nc.EQ.1) THEN
C cpb 22/1/97 Changing test
C                    IF((.NOT.ALL_REGIONS).OR.
C     '                ((nr.EQ.1).AND.(NP_INTERFACE(np,0)
C     '                .EQ.1)).OR.((nr.EQ.NRT).AND.(NP_INTERFACE(np,0)
C     '                .EQ.1))) THEN
C!                     Only dealing with one region, or dealing with them
C!                     all and np is in first reigon and no interface
C!                     (i.e. epicardium) or in the last region and no
C!                     interface (i.e outer torso).
                      IF((.NOT.ALL_REGIONS).OR.
     '                  (NP_INTERFACE(np,0).EQ.1)) THEN
C                       If we are on a 'free' surface or we want to
C                       specify a particular region(s)
                        IF(nk.EQ.1) THEN
C CPB 23/1/96 Changing dipole cases to actually use a dipole
C                        IF(NJT.EQ.3.AND.ANAL_CHOICE(nr).GE.5.AND.
C     '                    ANAL_CHOICE(nr).LE.8) THEN
C                          R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
C     '                      XP(1,nv,3,np)**2)
C                          IF(DABS(R-MESH3_RAD(NSPHERES))/
C     '                      MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL) THEN
C                            !Not on outer sphere
C                            YP(ny1,1)=YP(ny1,7)
C                            FIX(ny1,1)=.TRUE.
C                          ENDIF
                          IF(NJT.EQ.2) THEN
                            IF(ANAL_CHOICE(nr).EQ.5) THEN !single circle
                              IF(np.EQ.ANAL_FIXEDNODE) THEN
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !mulit circ
                              R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2)
                              IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                          MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL)
     '                          THEN
                                IF(np.EQ.ANAL_FIXEDNODE) THEN
C                                 First node on outer circle
                                  YP(ny1,1)=YP(ny1,7)
                                  FIX(ny1,1)=.TRUE.
                                ENDIF
                              ENDIF
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ELSE IF(NJT.EQ.3) THEN
                            IF(ANAL_CHOICE(nr).EQ.12.OR.
     '                        ANAL_CHOICE(nr).EQ.14) THEN !single sphere
                              IF(np.EQ.ANAL_FIXEDNODE) THEN
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE IF(ANAL_CHOICE(nr).EQ.13.OR.
     '                          ANAL_CHOICE(nr).EQ.15) THEN !multi sph
                              R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
     '                          XP(1,nv,3,np)**2)
                              IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                          MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL)
     '                          THEN
                                IF(np.EQ.ANAL_FIXEDNODE) THEN
C                               First node on outer sphere
                                  YP(ny1,1)=YP(ny1,7)
                                  FIX(ny1,1)=.TRUE.
                                ENDIF
                              ENDIF
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ENDIF
                        ELSE IF(HERMITE) THEN
                          IF(NKH(NH_LOC(1,nx),np,1).EQ.
     '                      NKH(NH_LOC(1,nx),np,2)) THEN
C CPB 23/1/96 Changing dipole cases to actually use a dipole
C!                         !Hermite interpolation used for normal
C!                         !derivative as well. Set the values of
C!                         !arclength derivatives.
C                          IF(NJT.EQ.3.AND.ANAL_CHOICE(nr).GE.5.AND.
C     '                      ANAL_CHOICE(nr).LE.8) THEN
C                            R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
C     '                        XP(1,nv,3,np)**2)
C                            IF(DABS(R-MESH3_RAD(NSPHERES))/
C '                            MESH3_RAD(NSPHERES).GT.SPHERE_RAD_TOL)
C     '                        THEN
C                              YP(ny1,1)=YP(ny1,7)
C                              FIX(ny1,1)=.TRUE.
C                            ENDIF
C                          ELSE IF(NJT.EQ.2.AND.ANAL_CHOICE(nr).EQ.6)
                            IF(NJT.EQ.3.AND.ANAL_CHOICE(nr).GE.12.AND.
     '                        ANAL_CHOICE(nr).LE.15) THEN
C                             Do Nothing
                            ELSE IF(NJT.EQ.2.AND.(ANAL_CHOICE(nr).GE.5
     '                          .AND.ANAL_CHOICE(nr).LE.6)) THEN
C                           Do Nothing
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDIF
                    ELSE
C                     Dipole solutions - set normal derivatives to zero
C                     on outer boundary
                      IF(NJT.EQ.2) THEN
                        IF(ANAL_CHOICE(nr).EQ.5) THEN !single circle
                          IF(FIX_ZERO) THEN
                            YP(ny1,1)=YP(ny1,7)
                            FIX(ny1,1)=.TRUE.
                          ELSE
                            IF(np.EQ.ANAL_FIXEDNODE.AND.nk.EQ.1) THEN
                              IF(POT_BC_TYPE.EQ.2) THEN !fix deriv,val
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ENDIF
                        ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !multi circles
                          R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2)
                          IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                      MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL) THEN
                            IF(FIX_ZERO) THEN
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ELSE
                              IF(np.EQ.ANAL_FIXEDNODE.AND.nk.EQ.1) THEN
                                IF(POT_BC_TYPE.EQ.2) THEN !fix deriv,val
                                  YP(ny1,1)=YP(ny1,7)
                                  FIX(ny1,1)=.TRUE.
                                ENDIF
                              ELSE
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDIF
                      ELSE IF(NJT.EQ.3) THEN
                        IF(ANAL_CHOICE(nr).EQ.12.OR.
     '                    ANAL_CHOICE(nr).EQ.14) THEN
                          IF(FIX_ZERO) THEN
                            YP(ny1,1)=YP(ny1,7)
                            FIX(ny1,1)=.TRUE.
                          ELSE
                            IF(np.EQ.ANAL_FIXEDNODE.AND.nk.EQ.1) THEN
                              IF(POT_BC_TYPE.EQ.2) THEN !fix deriv,val
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ENDIF
                        ELSE IF(ANAL_CHOICE(nr).EQ.13.OR.
     '                      ANAL_CHOICE(nr).EQ.15) THEN
                          R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
     '                      XP(1,nv,3,np)**2)
                          IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                      MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL) THEN
                            IF(FIX_ZERO) THEN
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ELSE
                              IF(np.EQ.ANAL_FIXEDNODE.AND.nk.EQ.1) THEN
                                IF(POT_BC_TYPE.EQ.2) THEN !fix deriv,val
                                  YP(ny1,1)=YP(ny1,7)
                                  FIX(ny1,1)=.TRUE.
                                ENDIF
                              ELSE
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDIF !nc
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nh
            ENDDO !nonode (np)
          ENDDO !nc

        ENDIF

      ELSEIF(GENER.AND.ACTIVATION) THEN
C***    Set up the YP(ny,1) from the analytic activation times
C***    stored in YP(ny,7).  Note that this really only makes sense
C***    to be done for the regions in which the activation times have
C***    been evaluated.  The following is set up only for nc=1 vars.
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            DO nv=1,NVHP(nh,np,1)
              DO nk=1,MAX(NKH(nh,np,1)-KTYP93(1,nr),1)
                ny1=NYNP(nk,nv,nh,np,0,1,nr)
                YP(ny1,1)=YP(ny1,7)
                IF(.NOT.NOFIX)FIX(ny1,1)=.TRUE.
              ENDDO !nk
            ENDDO !nv
          ENDDO !nhx
        ENDDO !nonode



      ELSE !not generating initial conditions
        DO loop=1,2
          nonode=0
          IF(loop.EQ.1) THEN
            FORMAT='(/'' Essential boundary conditions'//
     '       ' defined at nodes:'')'
            nc=1
          ELSE IF(loop.EQ.2) THEN
            FORMAT='(/'' Flux boundary conditions'//
     '       ' defined at nodes:'')'
            nc=2
          ENDIF
          CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '       LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
 6100     FORMAT='($,'' Enter node #s/name [EXIT]: '',I5)'
          CORNER=.FALSE.
          IF(IOTYPE.EQ.3) THEN
            nonode=nonode+1
            IF(nonode.LE.NPNODE(0,nr)) THEN
              np=NPNODE(nonode,nr)
              IDATA(1)=np
              NPLIST(0)=1
              NPLIST(1)=np
            ELSE
              IDATA(1)=0
            ENDIF
          ENDIF

 6500     CDATA(1)='NODES' !for use with group input
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).NE.0) THEN !not default exit
            IF(iotype.NE.3) THEN
              NPLIST(0)=IDATA(0)
              DO n=1,IDATA(0)
                NPLIST(n)=IDATA(n)
                np=IDATA(n)
                IF(.NOT.INLIST(np,NPNODE(1,nr),NPNODE(0,nr),N1)) THEN
                  WRITE(OP_STRING,*)'>>This node is not in the',
     '              ' current region'
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 6500
                ENDIF
              ENDDO !n
            ENDIF !iotype.NE.3
            IF(NVHP(1,np,nc).GT.1) THEN
              CORNER=.TRUE.
            ENDIF

C           Define bdry condition for first node in group
            np=NPLIST(1) !rest of group is filled at end of nh loop
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              WRITE(CHAR1,'(I1)') nhx
              FORMAT='('' Dependent variable number '//CHAR1(1:1)//
     '          ' :'')'
              CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                nv=1 !temporary
                NYTOT=NYNP(nk,nv,nh,np,0,nc,nr)
                RDEFLT(1)=1.1d6
                IF(nk.EQ.1) THEN
                  IF((.NOT.CORNER).OR.(loop.NE.2)) THEN
                    FORMAT='($,'' The value of the dependent variable'//
     '                ' is [no b.c.]: '',G12.5)'
                  ELSE ! Loop=2 and at a corner node
! find ne           CHAR2=CFROMI(NVHP(1,np,nc),'(I3)')
C ***               lowest number corner element
c                    FORMAT='($,'' The value of the dependent'//
c     '                ' variable on '//
c     '                'element '//CHAR2(1:3)//' is [no b.c.]: '',G12.5)'
                  ENDIF
                ELSE IF(nk.GT.1) THEN
                  WRITE(CHAR1,'(I1)') nk
                  IF((.NOT.CORNER).OR.(loop.NE.2)) THEN
                    FORMAT='($,'' The value of derivative number '//
     '                CHAR1(1:1)//' is [no b.c.]: '',G12.5,1X,'//
     '                'G12.5)'
                  ELSE !Request information on all values for du/dn
!find ne            CHAR2=CFROMI(NVHP(1,np,nc),'(I3)')
                    !lowest number corner element
c                    FORMAT='($,'' The value of derivative number '//
c     '                CHAR1(1:1)//' on element '//
c     '                CHAR2(1:3)//' is [no b.c.]: '',G12.5,1X,'//
c     '                'G12.5)'
                  ENDIF ! End of no corner loop
                ENDIF ! End of nk.gt.1 loop
                IF(IOTYPE.EQ.3) THEN
                  IF(loop.EQ.1) THEN
                    IF(nk.LE.MAX(NKH(nh,np,2)-KTYP93(2,nr),1)) THEN
                      ny2=NYNP(nk,nv,nh,np,0,2,nr)
                    ELSE !no ny2 value exists
                      ny2=0
                    ENDIF
                    IF((.NOT.CORNER)) THEN
                      IF(ny2.GT.0) THEN
                        IF(FIX(NYTOT,1).AND.(.NOT.FIX(ny2,1))) THEN
                          RDATA(1)=YP(NYTOT,1)
                        ELSE
                          RDATA(1)=RDEFLT(1)
                        ENDIF
                      ELSE
                        RDATA(1)=RDEFLT(1)
                      ENDIF
                    ELSE !Corner
                      IF(ny2.GT.0) THEN
                        IF(FIX(NYTOT,1).AND.(.NOT.FIX(ny2,1))) THEN
                          RDATA(1)=YP(NYTOT,1)
                        ELSE
                          RDATA(1)=RDEFLT(1)
                        ENDIF
                      ELSE
                        RDATA(1)=RDEFLT(1)
                      ENDIF
                    ENDIF
                  ELSE IF(loop.EQ.2) THEN
                    IF(nk.LE.MAX(NKH(nh,np,1)-KTYP93(1,nr),1)) THEN
                      ny1=NYNP(nk,nv,nh,np,0,1,nr)
                    ELSE !no ny1 value exists
                      ny1=0
                    ENDIF
                    IF(ny1.GT.0) THEN
                      IF(FIX(NYTOT,1).AND.(.NOT.FIX(ny1,1))) THEN
                        RDATA(1)=YP(NYTOT,1)
                      ELSE
                        RDATA(1)=RDEFLT(1)
                      ENDIF
                    ELSE
                      RDATA(1)=RDEFLT(1)
                    ENDIF
                  ENDIF
                ENDIF
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '            FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '            ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '            RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
                IF(RDATA(1).LT.1.0D6) THEN
                  IF(iotype.NE.3) THEN

C LKC 12-JUL-2000 Setup nofixing ....
C
C                    IF(loop.EQ.1) THEN
C                      YP(NYTOT,1)=RDATA(1)
C                      FIX(NYTOT,1)=.TRUE.
C                    ELSE IF(loop.EQ.2) THEN
C                      ny2=NYNP(nk,nv,nh,np,0,2,nr)
C                      YP(ny2,1)=RDATA(1)
C                      FIX(ny2,1)=.TRUE.
C                    ENDIF

                    IF(loop.EQ.1) THEN
                      YP(NYTOT,1)=RDATA(1)
                      IF(NOFIX) THEN
                        FIX(NYTOT,1)=.FALSE.
                      ELSE
                        FIX(NYTOT,1)=.TRUE.
                      ENDIF
                    ELSE IF(loop.EQ.2) THEN
                      ny2=NYNP(nk,nv,nh,np,0,2,nr)
                      YP(ny2,1)=RDATA(1)
                      IF(NOFIX) THEN
                        FIX(ny2,1)=.FALSE.
                      ELSE
                        FIX(ny2,1)=.TRUE.
                      ENDIF
                    ENDIF

C LKC 3-MAR-2000 adding assert
                    CALL ASSERT(NIYFIXM.GE.5,
     '                '>> Increase NIYFIXM to 5',ERROR,*9999)


                    IF(FILEIP) THEN
                      FIX(NYTOT,5)=.TRUE.
                    ELSE IF(.NOT.FILEIP) THEN
                      FIX(NYTOT,5)=.FALSE.
                    ENDIF
                  ENDIF
                ENDIF
                IF((loop.EQ.2).AND.CORNER) THEN
C ***             Enter other corner values of du/dn
                  CORNER=.FALSE.
!find ne          CHAR2=CFROMI(NVHP(1,np,nc),'(I3)')
                  IF(IOTYPE.EQ.3) THEN
                    NY3=NYNP(nk,2,nh,np,0,1,nr)
                    IF(FIX(NY3,1).AND.(.NOT.FIX(NYTOT,1))) THEN
                      RDATA(1)=YP(NY3,1) !Middle corner element
                    ELSE
                      RDATA(1)=RDEFLT(1)
                    ENDIF
                  ENDIF
c                  FORMAT='($,'' The value of the dependent variable '//
c     '              'on element '//CHAR2(1:3)//' is [no b.c.]: '','//
c     '              'G12.5)'
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '              FILEIP,FORMAT,1,
     '              ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     '              IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '              ERROR,*9999)
                  IF(RDATA(1).LT.1.0d6) THEN
                    IF(iotype.NE.3) THEN
                      NY3=NYNP(nk,2,nh,np,0,1,nr)
                      YP(NY3,1)=RDATA(1)
                      FIX(NY3,1)=.TRUE.
                      IF(FILEIP) THEN
                        FIX(NY3,5)=.TRUE.
                      ELSE IF(.NOT.FILEIP) THEN
                        FIX(NY3,5)=.FALSE.
                      ENDIF
                    ENDIF
                  ENDIF
                  IF(NVHP(nh,np,nc).GT.2) THEN !3d corner
!find   ne          CHAR2=CFROMI(NVHP(1,np,nc),'(I3)')
c                    FORMAT='($,'' The value of the dependent '//
c     '                'variable on element '//
c     '                CHAR2(1:3)//' is [no b.c.]: '',G12.5)'
                    IF(IOTYPE.EQ.3) THEN
                      NY4=NYNP(nk,3,nh,np,0,1,nr)
                      IF(FIX(NY4,1).AND.(.NOT.FIX(NYTOT,1))) THEN
                        RDATA(1)=YP(NY4,1)
                      ELSE
                        RDATA(1)=RDEFLT(1)
                      ENDIF
                    ENDIF
                    CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '                FILEIP,FORMAT,1,
     '                ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     '                IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '                ERROR,*9999)
                    IF(RDATA(1).LT.1.0D6) THEN
                      IF(iotype.NE.3) THEN
                        NY4=NYNP(nk,3,nh,np,0,1,nr)
                        YP(NY4,1)=RDATA(1)
                        FIX(NY4,1)=.TRUE.
                        IF(FILEIP) THEN
                          FIX(NY4,5)=.TRUE.
                        ELSE
                          FIX(NY4,5)=.FALSE.
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF ! End of 3D corner loop
                ENDIF !End of corner loop
              ENDDO ! End of nk loop
            ENDDO !nh

C           Apply bdry conditions to rest of group
            DO n=2,NPLIST(0)
              np=NPLIST(n)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              IF(DOP) THEN
C$OMP CRITICAL(IPINI9_1)
                WRITE(*,'('' Group data: n='',I3,'' np='',I5)') n,np
C$OMP END CRITICAL(IPINI9_1)
              ENDIF
              DO nhx=1,NHP(np)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,nc)
                  DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                    IF(DOP) WRITE(*,'('' nh='',I1,'' nv='',I1,'
     '                //''' nk='',I1)') nh,nv,nk
                    ny=NYNP(nk,nv,nh,np,0,nc,nr)
                    ny_first =NYNP(nk,nv,nh,NPLIST(1),0,nc,nr)
                    YP(ny,1) = YP(ny_first,1)
                    FIX(ny,1)=FIX(ny_first,1)
                    IF(DOP) WRITE(*,'('' FIX('',I4,'',1)='',L1)')
     '                ny,FIX(ny,1)
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nh
            ENDDO !n

            GO TO 6100 !for more nodes
          ENDIF !End of IDATA(1) loop

          IF(loop.EQ.1) THEN !essential bdry conditions
C CPB ??? 6/12/94
            DO no_nynr=1,NYNR(0,0,nc) ! loop over global variables
              ny=NYNR(no_nynr,0,nc) ! global variable number
              IF(NPNY(0,ny,0,nx).EQ.1) THEN
                np=NPNY(4,ny,0,nx)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              ENDIF
              YP(ny,4)=YP(ny,1) !to include any essential bcs in YP(4)
            ENDDO
          ENDIF
        ENDDO !End of loop loop
      ENDIF !gener

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(IPINI9_2)
        FORMAT='(/'' YP array:'')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        FORMAT='('' YP('',I4,'',1)= '',D12.6)'
        DO nc=1,NCT(nr,nx)
          DO no_nynr=1,NYNR(0,0,nc) ! loop over global variables
            ny=NYNR(no_nynr,0,nc) ! global variable number
            WRITE(OP_STRING,FORMAT) ny,YP(ny,1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !no_nynr
        ENDDO !nc
        FORMAT='(/'' FIX array:'')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        FORMAT='('' FIX('',I4,'',1)= '',L1)'
        DO nc=1,NCT(nr,nx)
          DO no_nynr=1,NYNR(0,0,nc) ! loop over global variables
            ny=NYNR(no_nynr,0,nc) ! global variable number
            WRITE(OP_STRING,FORMAT) ny,FIX(ny,1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !no_nynr
        ENDDO !nc
CC$OMP END CRITICAL(IPINI9_2)
      ENDIF

      CALL EXITS('IPINI9')
      RETURN
 9999 CALL ERRORS('IPINI9',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPINI9')
      RETURN 1
      END


