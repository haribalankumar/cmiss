      SUBROUTINE NETTFLUX(NFVC,NODENVC,NPNODE,nr,NVCB,NVCNODE,nx,
     '  NYNP,VC,XNFV,YP,ZNFV,ZNFVMSH,ERROR,*)

C#### Subroutine: NETTFLUX
C###  Description:
C###    <HTML>
C###    This subroutine calculates the nett flux of momentum into each
C###    cell by advection and diffusion. This includes contributions
C###    from advective boundaries and shear stress boundaries.
C###    Calculated by:
C###    <PRE>
C###
C###    NettFlux_(i) = (1/Re) *
C###                   Sum_j(Un_(j) - Un(i))*(A(i,j)/D(i,j) -
C###                   Sum_j<Un_(i,j)>*(RCFlux(i,j) - MeshFlux(i,j))
C###
C###
C###    Where Re       = Reynolds number
C###          Sum_j    = sum over neighbouring j cells
C###          Un       = Fluid velocity at time step n
C###          Area     = Area of Voronoi face
C###          D        = distance between two nodes
C###          RCFlux   = Rhie Chow Flux
C###          MeshFlux = Flux of the moving mesh
C###    </PRE>
C###    </HTML>
      IMPLICIT NONE
      INCLUDE 'b12.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'voro00.inc'
!     Parameter Variables
      INTEGER NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NVCB(-1:3,NVCBM),NVCNODE(2,NP_R_M),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 VC(0:NVCM),XNFV(-(NJM+1):NJM,NFVM),YP(NYM,NIYM),
     '  ZNFV(NFVM),ZNFVMSH(NFVM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nonode,np,nhx,nh,nfvl,nvc,cnonode,cnp,nfv,ny,cny,bnvc,
     '  INODE,i,nfvl2,cnonode2,cnp2,nfv2,cpny,pny,cnvc
      REAL*8 DIFFLUX(3),VELFLUX(3),VELC(3),TOTFLUX(3),PROJ,A,DDOT,
     '  NORMALV,QUICK,ULTRAQUICK,QUICKEST,ULTRAQUICKEST,CL,PE
      LOGICAL OPPNODE

      CALL ENTERS('NETTFLUX',*9999)

      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)

C       ..Initialise diffusive flux, flux due to velocity, and the
C       combined total flux
        DO nhx=1,nh_loc(0,nx)-1
          DIFFLUX(nhx)=0.d0
          VELFLUX(nhx)=0.d0
          TOTFLUX(nhx)=0.d0
        ENDDO

C       ..Loop over the adjoing nodes of Voronoi cell nvc
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)
          cnp=NPNODE(cnonode,nr)
          cnvc=NVCNODE(MAP,cnonode)
          nfv=NFVC(2,nfvl,nvc)

C         ..Internal connections
          IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN


C           ..Calculate the face normal component of velocity
C           =Rhie Chow face flux (minus the mesh flux if mesh moving)
            IF(.NOT.MESHFIXD) THEN
              NORMALV=ZNFV(nfv)-ZNFVMSH(nfv)
            ELSE
              NORMALV=ZNFV(nfv)
            ENDIF

C           ..Loop over the velocities
            DO nhx=1,nh_loc(0,nx)-1 ! Loop over velocities
              nh=nh_loc(nhx,nx)
              ny=NYNP(1,1,nh,np,0,1,nr)
              cny=NYNP(1,1,nh,cnp,0,1,nr)

C             ..Calculate the diffusive flux through the cell walls
              DIFFLUX(nhx)=INV_REYNOLD*(YP(cny,8)-YP(ny,8))
     '          *XNFV(IDIST,nfv)*XNFV(FAREA,nfv)

C             ..Calculate the advective velocity term

C             ..Upwind Differencing
              IF(KTYP61.EQ.0) THEN
                IF(NORMALV.GT.0.d0) THEN
                  VELFLUX(nhx)=YP(ny,8)
                ELSE
                  VELFLUX(nhx)=YP(cny,8)
                ENDIF

C               ..Central Differencing
              ELSEIF(KTYP61.EQ.1) THEN
                VELFLUX(nhx)=0.5d0*(YP(ny,8)+YP(cny,8))

C               ..Hybrid Differencing
              ELSEIF(KTYP61.EQ.2) THEN
                PE=DABS(0.5d0*REYNOLD*NORMALV/
     '            (XNFV(FAREA,nfv)*XNFV(IDIST,nfv)))
                IF(PE.GT.2.d0) THEN
                  IF(NORMALV.GT.0.d0) THEN
                    VELFLUX(nhx)=YP(ny,8)
                  ELSE
                    VELFLUX(nhx)=YP(cny,8)
                  ENDIF
C                 ..Diffusive cut off
                  DIFFLUX(nhx)=0.d0
                ELSE
                  VELFLUX(nhx)=0.5d0*(YP(ny,8)+YP(cny,8))
                ENDIF

C               ..Interpolated Donor Differencing
              ELSEIF(KTYP61.EQ.3) THEN
                CL=NORMALV*DT*XNFV(IDIST,nfv)/XNFV(FAREA,nfv)
                VELFLUX(nhx)=
     '            0.5d0*(YP(ny,8)+YP(cny,8)-CL*(YP(cny,8)-YP(ny,8)))

C               ..QUICK Differencing
              ELSEIF(KTYP61.EQ.4) THEN
                IF(NORMALV.GT.0.d0) THEN
                  VELFLUX(nhx)=QUICK(nfv,NFVC,nh,NPNODE,nr,nvc,
     '              NVCNODE,NYNP,YP(ny,8),YP(cny,8),VC,XNFV,YP)
                ELSE
                  VELFLUX(nhx)=QUICK(nfv,NFVC,nh,NPNODE,nr,cnvc,
     '              NVCNODE,NYNP,YP(cny,8),YP(ny,8),VC,XNFV,YP)
                ENDIF

C               ..ULTRA-QUICK Differencing
              ELSEIF(KTYP61.EQ.5) THEN
                IF(NORMALV.GT.0.d0) THEN
                  VELFLUX(nhx)=ULTRAQUICK(nfv,NFVC,nh,NPNODE,nr,nvc,
     '              NVCNODE,NYNP,YP(ny,8),YP(cny,8),VC,XNFV,YP)
                ELSE
                  VELFLUX(nhx)=ULTRAQUICK(nfv,NFVC,nh,NPNODE,nr,cnvc,
     '              NVCNODE,NYNP,YP(cny,8),YP(ny,8),VC,XNFV,YP)
                ENDIF

C               ..QUICK-EST differencing
              ELSEIF(KTYP61.EQ.6) THEN
                CL=NORMALV*DT*XNFV(IDIST,nfv)/XNFV(FAREA,nfv)
                IF(NORMALV.GT.0.d0) THEN
                  VELFLUX(nhx)=QUICKEST(nfv,NFVC,nh,NPNODE,nr,nvc,
     '              NVCNODE,NYNP,DABS(CL),YP(ny,8),YP(cny,8),VC,XNFV,YP)
                ELSE
                  VELFLUX(nhx)=QUICKEST(nfv,NFVC,nh,NPNODE,nr,cnvc,
     '              NVCNODE,NYNP,DABS(CL),YP(cny,8),YP(ny,8),VC,XNFV,YP)
                ENDIF

C               ..ULTRA-QUICK-EST differencing
              ELSEIF(KTYP61.EQ.7) THEN
                CL=NORMALV*DT*XNFV(IDIST,nfv)/XNFV(FAREA,nfv)
                IF(NORMALV.GT.0.d0) THEN
                  VELFLUX(nhx)=ULTRAQUICKEST(nfv,NFVC,nh,NPNODE,nr,nvc,
     '              NVCNODE,NYNP,DABS(CL),YP(ny,8),YP(cny,8),VC,XNFV,YP)
                ELSE
                  VELFLUX(nhx)=ULTRAQUICKEST(nfv,NFVC,nh,NPNODE,nr,cnvc,
     '              NVCNODE,NYNP,DABS(CL),YP(cny,8),YP(ny,8),VC,XNFV,YP)
                ENDIF

C               ..FULL-QUICK differencing
              ELSEIF(KTYP61.EQ.8) THEN
                IF(NORMALV.GT.0.d0) THEN
                  VELFLUX(nhx)=QUICK(nfv,NFVC,nh,NPNODE,nr,nvc,
     '              NVCNODE,NYNP,YP(ny,8),YP(cny,8),VC,XNFV,YP)
                ELSE
                  VELFLUX(nhx)=QUICK(nfv,NFVC,nh,NPNODE,nr,cnvc,
     '              NVCNODE,NYNP,YP(cny,8),YP(ny,8),VC,XNFV,YP)
                ENDIF

C               ..ULTRA-FULL-QUICK differencing
              ELSEIF(KTYP61.EQ.9) THEN
                IF(NORMALV.GT.0.d0) THEN
                  VELFLUX(nhx)=ULTRAQUICK(nfv,NFVC,nh,NPNODE,nr,nvc,
     '              NVCNODE,NYNP,YP(ny,8),YP(cny,8),VC,XNFV,YP)
                ELSE
                  VELFLUX(nhx)=ULTRAQUICK(nfv,NFVC,nh,NPNODE,nr,cnvc,
     '              NVCNODE,NYNP,YP(cny,8),YP(ny,8),VC,XNFV,YP)
                ENDIF
              ENDIF

C             ..Total flux=diffusive flux minus the advective flux times
C             the face fluxes
              TOTFLUX(nhx)=TOTFLUX(nhx)+DIFFLUX(nhx)
     '          -NORMALV*VELFLUX(nhx)
            ENDDO !nhx

C           ..Boundary connection - apply the BC's
          ELSE
            bnvc=NVCNODE(MAP,cnonode)

C           ..First deal with velocity BC's
            IF(NVCB(BCTYPE,bnvc).EQ.INLET) THEN

              DO nhx=1,nh_loc(0,nx)-1
                nh=nh_loc(nhx,nx)
                cny=NYNP(1,1,nh,cnp,0,1,nr)
                ny=NYNP(1,1,nh,np,0,1,nr)

C rgb if you want to see the same results as CJ Were's 2D code,
C comment out this line to skip the additional Laplacian calculation
                TOTFLUX(nhx)=TOTFLUX(nhx)+INV_REYNOLD*
     '            (YP(cny,3)-YP(ny,8))*XNFV(FAREA,nfv)*XNFV(IDIST,nfv)
                IF(MESHFIXD) THEN
                  TOTFLUX(nhx)=TOTFLUX(nhx)-YP(cny,3)*ZNFV(nfv)
                ELSE
                  TOTFLUX(nhx)=TOTFLUX(nhx)-YP(cny,3)*
     '              (ZNFV(nfv)-ZNFVMSH(nfv))
                ENDIF

              ENDDO

            ELSEIF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN

C             ..Outflow
              IF(ZNFV(nfv).GT.0.d0) THEN
                DO nhx=1,nh_loc(0,nx)-1
                  nh=nh_loc(nhx,nx)
                  ny=NYNP(1,1,nh,np,0,1,nr)
                  TOTFLUX(nhx)=TOTFLUX(nhx)-YP(ny,8)*ZNFV(nfv)
                ENDDO

C               ..Inflow
              ELSE
                PROJ=0.d0
                DO nhx=1,nh_loc(0,nx)-1
                  nh=nh_loc(nhx,nx)
                  ny=NYNP(1,1,nh,np,0,1,nr)
                  PROJ=PROJ+YP(ny,8)*XNFV(nhx,nfv)
                ENDDO
                DO nhx=1,nh_loc(0,nx)-1
                  TOTFLUX(nhx)=TOTFLUX(nhx)-PROJ*XNFV(nhx,nfv)*
     '              ZNFV(nfv)
                ENDDO
              ENDIF

            ELSEIF(NVCB(BCTYPE,bnvc).EQ.WALL.OR.
     '          NVCB(BCTYPE,bnvc).EQ.DRIVING) THEN

              DO nhx=1,nh_loc(0,nx)-1
                nh=nh_loc(nhx,nx)
                cny=NYNP(1,1,nh,cnp,0,1,nr)
                ny=NYNP(1,1,nh,np,0,1,nr)

                TOTFLUX(nhx)=TOTFLUX(nhx)+INV_REYNOLD*2.d0*
     '            (YP(cny,3)-YP(ny,8))*
     '            XNFV(FAREA,nfv)*XNFV(IDIST,nfv)
              ENDDO
            ENDIF

C           ..Now deal with pressure BC's
            IF(NVCB(BCTYPE,bnvc).EQ.INLET.OR.
     '        NVCB(BCTYPE,bnvc).EQ.WALL.OR.
     '        NVCB(BCTYPE,bnvc).EQ.FREESLIP.OR.
     '        NVCB(BCTYPE,bnvc).EQ.DRIVING) THEN

              pny=NYNP(1,1,nh_loc(0,nx),np,0,nr,nx)
              A=-0.5d0*YP(pny,8)*XNFV(FAREA,nfv)
              DO nhx=1,nh_loc(0,nx)-1
                TOTFLUX(nhx)=TOTFLUX(nhx)+A*XNFV(nhx,nfv)
              ENDDO

            ELSEIF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN

C             ..Ensure that we only work with neighbouring cell
              OPPNODE=.FALSE.
              DO i=1,NVCB(0,bnvc)
                INODE=NVCB(i,bnvc)
                IF(NPNODE(INODE,nr).EQ.np) THEN
                  OPPNODE=.TRUE.
                  GOTO 100
                ENDIF
              ENDDO
 100          CONTINUE
              IF(OPPNODE) THEN

                pny=NYNP(1,1,nh_loc(0,nx),np,0,nr,nx)
                A=-0.5d0*YP(pny,8)*XNFV(FAREA,nfv)
                DO nhx=1,nh_loc(0,nx)-1
                  TOTFLUX(nhx)=TOTFLUX(nhx)+A*XNFV(nhx,nfv)

C                 ..Use the same loop to initialise velc
                  VELC(nhx)=0.d0
                ENDDO

                DO nfvl2=1,NFVC(1,0,nvc)
                  cnonode2=NFVC(1,nfvl2,nvc)

                  IF(NVCNODE(TYPE,cnonode2).NE.BOUNDARY) THEN
                    cnp2=NPNODE(cnonode2,nr)
                    nfv2=NFVC(2,nfvl2,nvc)
                    cpny=NYNP(1,1,nh_loc(0,nx),cnp2,0,1,nr)
                    A=0.5d0*YP(cpny,8)*XNFV(FAREA,nfv2)
                    DO nhx=1,nh_loc(0,nx)-1
                      VELC(nhx)=VELC(nhx)-A*XNFV(nhx,nfv2)
                    ENDDO
                  ENDIF
                ENDDO

                A=0.5d0*YP(pny,8)*XNFV(FAREA,nfv)
                DO nhx=1,nh_loc(0,nx)-1
                  VELC(nhx)=VELC(nhx)-A*XNFV(nhx,nfv)
                ENDDO

                PROJ=DDOT(NJ_LOC(NJL_GEOM,0,nr),VELC,1,XNFV(1,nfv),1)
                DO nhx=1,nh_loc(0,nx)-1
                  TOTFLUX(nhx)=TOTFLUX(nhx)+PROJ*XNFV(nhx,nfv)
                ENDDO

              ENDIF !Oppnode
            ENDIF !BC type
          ENDIF !Boundary or internal connection
        ENDDO !nfvl
        DO nhx=1,nh_loc(0,nx)-1
          nh=nh_loc(nhx,nx)
          ny=NYNP(1,1,nh,np,0,1,nr)
          YP(ny,1)=TOTFLUX(nhx)
        ENDDO

      ENDDO

      CALL EXITS('NETTFLUX')
      RETURN
 9999 CALL ERRORS('NETTFLUX',ERROR)
      CALL EXITS('NETTFLUX')
      RETURN 1
      END


