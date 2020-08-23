      SUBROUTINE CORRVELO(NFVC,NODENVC,NPNODE,nr,NVCB,NVCBNODE_OUTLET,
     '  NVCNODE,nx,NYNP,N_OUTLET,P_OUTLET,VC,XNFV,YP,ERROR,*)

C#### Subroutine: CORRVELO
C###  Description:
C###    <HTML>
C###    Correct velocities using the cell centred pressure correction
C###    <PRE>
C###
C###    Un+1_(i) = U'_(i) -
C###               [dt/2V(i)]*
C###               Sum_j[P(j)*A(i,j).N(i,j)]
C###
C###    Where Un+1  = new velocity at time step n+1
C###          A     = Voronoi face area
C###          U'    = momentum marched non-solenoidal estimates of
C###                  velocity
C###          N     = unit outward normal from cell i
C###          dt    = time step
C###          V     = cell volume
C###          P     = pressure
C###          Sum_j = sum over neighbouring j cells
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NVCB(-1:3,NVCBM),NVCBNODE_OUTLET(NP_R_M),NVCNODE(2,NP_R_M),
     '  nx,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),N_OUTLET
      REAL*8 P_OUTLET(N_OUTLET),VC(0:NVCM),XNFV(-(NJM+1):NJM,NFVM),
     '  YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nvc,nonode,np,nhx,nfvl,cnonode,cnp,cny,ny,nfv,nh,
     '  bnvc,bnfv,cpny,nfvl2,nfv2,cnonode2,cnp2,pny,i,njj,j,nhx2,
     '  nh2,ny2
      REAL*8 VELCOR(3),DUDX(3,3),UC(3),GRADSQUN,DUNDT,DPDN,A,
     '  SUM,TEMP,AREA,UN,PROJ
      LOGICAL OPPNODE

      CALL ENTERS('CORRVELO',*9999)

C     ..Calculate internal fluxes..

      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)
        DO nhx=1,nh_loc(0,nx)-1
          VELCOR(nhx)=0.d0
        ENDDO

C       ..Loop over the adjoining nodes, but not the boundary ones,
C       computing the velocity correction VELCOR
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)
          IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
            cnp=NPNODE(cnonode,nr)
            nfv=NFVC(2,nfvl,nvc)
            cny=NYNP(1,1,nh_loc(0,nx),cnp,0,1,nr)
            A=0.5d0*YP(cny,1)*XNFV(FAREA,nfv)
            DO nhx=1,nh_loc(0,nx)-1
              VELCOR(nhx)=VELCOR(nhx)-A*XNFV(nhx,nfv)
            ENDDO
          ENDIF
        ENDDO

C       ..Apply the velocity correction to the non-solenoidal estimates
        DO nhx=1,nh_loc(0,nx)-1
          nh=nh_loc(nhx,nx)
          VELCOR(nhx)=VELCOR(nhx)/VC(nvc)
          ny=NYNP(1,1,nh,np,0,1,nr)
          YP(ny,1)=YP(ny,1)+DT*VELCOR(nhx)
        ENDDO

      ENDDO

C     ..Calculate boundary fluxes..

      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)

C       ..Pressure ny location of nonode
        pny=NYNP(1,1,nh_loc(0,nx),np,0,1,nr)

        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)

C         ..Make sure that this connection is a boundary connection
          IF(NVCNODE(TYPE,cnonode).EQ.BOUNDARY) THEN
            bnvc=NVCNODE(MAP,cnonode)
            bnfv=NFVC(2,nfvl,nvc)

C           ..Initialisation
            DO nhx=1,nh_loc(0,nx)-1
              VELCOR(nhx)=0.d0
            ENDDO

            IF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN

C             ..Check that only the exact opposite node is used
              OPPNODE=.FALSE.
              DO i=1,NVCB(0,bnvc)
                IF(NVCB(i,bnvc).EQ.nonode) THEN
                  OPPNODE=.TRUE.
                  GOTO 100
                ENDIF
              ENDDO
 100          CONTINUE

              IF(OPPNODE) THEN

                IF(DUCTFLOW) THEN
                  A=-0.5d0*P_OUTLET(NVCBNODE_OUTLET(NVCB(1,bnvc)))*
     '              XNFV(FAREA,bnfv)
                  DO nhx=1,nh_loc(0,nx)-1
                    VELCOR(nhx)=A*XNFV(nhx,bnfv)
                  ENDDO

                ELSE
                  A=-0.5d0*YP(pny,1)*XNFV(FAREA,bnfv)
                  DO nhx=1,nh_loc(0,nx)-1
                    VELCOR(nhx)=A*XNFV(nhx,bnfv)
                  ENDDO

                  DO nhx=1,nh_loc(0,nx)-1
                    UC(nhx)=0.d0
                  ENDDO

C                 ..Loop over adjoining internal nodes
                  DO nfvl2=1,NFVC(1,0,nvc)
                    cnonode2=NFVC(1,nfvl2,nvc)

                    IF(NVCNODE(TYPE,cnonode2).NE.BOUNDARY) THEN
                      cnp2=NPNODE(cnonode2,nr)
                      nfv2=NFVC(2,nfvl2,nvc)
                      cpny=NYNP(1,1,nh_loc(0,nx),cnp2,0,1,nr)
                      A=0.5d0*YP(cpny,1)*XNFV(FAREA,nfv2)
                      DO nhx=1,nh_loc(0,nx)-1
                        UC(nhx)=UC(nhx)-A*XNFV(nhx,nfv2)
                      ENDDO
                    ENDIF
                  ENDDO
                  A=0.5d0*YP(pny,1)*XNFV(FAREA,bnfv)
                  DO nhx=1,nh_loc(0,nx)-1
                    UC(nhx)=UC(nhx)-A*XNFV(nhx,bnfv)
                  ENDDO
                  PROJ=0.d0

                  DO nhx=1,nh_loc(0,nx)-1
                    PROJ=PROJ+UC(nhx)*XNFV(nhx,bnfv)
                  ENDDO

                  DO nhx=1,nh_loc(0,nx)-1
                    VELCOR(nhx)=VELCOR(nhx)+PROJ*XNFV(nhx,bnfv)
                  ENDDO

                ENDIF !ductflow
              ENDIF !oppnode

            ELSEIF(NVCB(BCTYPE,bnvc).EQ.WALL.OR.
     '          NVCB(BCTYPE,bnvc).EQ.FREESLIP.OR.
     '          NVCB(BCTYPE,bnvc).EQ.DRIVING) THEN

C             ..Initialise gradients
              DO i=1,nh_loc(0,nx)-1
                DO j=1,NJ_LOC(NJL_GEOM,0,nr)
                  DUDX(i,j)=0.d0
                ENDDO
              ENDDO
              GRADSQUN=0.d0 !Viscous term

C             ..Calculate normal component of velocity
              UN=0.d0
              DO nhx=1,nh_loc(0,nx)-1
                nh=nh_loc(nhx,nx)
                ny=NYNP(1,1,nh,np,0,1,nr)
                UN=UN+YP(ny,8)*XNFV(nhx,bnfv)
              ENDDO

C             ..Calculate velocity gradients
              DO nfvl2=1,NFVC(1,0,nvc)
                cnonode2=NFVC(1,nfvl2,nvc)

                IF(NVCNODE(TYPE,cnonode2).NE.BOUNDARY) THEN
                  cnp2=NPNODE(cnonode2,nr)
                  nfv2=NFVC(2,nfvl2,nvc)
                  AREA=0.5d0*XNFV(FAREA,bnfv)
                  TEMP=0.d0

                  DO nhx=1,nh_loc(0,nx)-1
                    nh=nh_loc(nhx,nx)
                    ny=NYNP(1,1,nh,cnp2,0,1,nr)
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      DUDX(nhx,njj)=DUDX(nhx,njj)+AREA*YP(ny,8)*
     '                  XNFV(njj,nfv2)
                    ENDDO
                    TEMP=TEMP+YP(ny,8)*XNFV(nhx,bnfv)
                  ENDDO
                  GRADSQUN=GRADSQUN+(TEMP-UN)*
     '              XNFV(FAREA,nfv2)*XNFV(IDIST,nfv2)

                ENDIF

              ENDDO

C             ..Divide velocity gradients through by the volume
              DO i=1,nh_loc(0,nx)-1
                DO j=1,NJ_LOC(NJL_GEOM,0,nr)
                  DUDX(i,j)=DUDX(i,j)/VC(nvc)
                ENDDO
              ENDDO
              GRADSQUN=GRADSQUN/VC(nvc)

C             ..Calculate tangential component of velocity
              DUNDT=0.d0
              DO nhx=1,nh_loc(0,nx)-1
                nh=nh_loc(nhx,nx)
                ny=NYNP(1,1,nh,np,0,1,nr)
                DUNDT=DUNDT+(YP(ny,8)-YP(ny,12))*XNFV(nhx,bnfv)
              ENDDO

C             ..Calculate normal derivative of pressure
C             dp/dn = (1/Re)*Laplacian(Un) - d(Un)/dt - U.Grad(Un)
              SUM=0.d0
              DO nhx=1,nh_loc(0,nx)-1
                nh=nh_loc(nhx,nx)
                ny=NYNP(1,1,nh,np,0,1,nr)
                TEMP=0.d0
                DO nhx2=1,nh_loc(0,nx)-1
                  nh2=nh_loc(nhx2,nx)
                  ny2=NYNP(1,1,nh2,np,0,1,nr)
                  TEMP=TEMP+YP(ny2,8)*DUDX(nhx,nhx2)
                ENDDO
                SUM=SUM+TEMP*XNFV(nhx,bnfv)
              ENDDO
              DPDN=INV_REYNOLD*GRADSQUN-SUM-DUNDT

C             ..Calculate velocity correction
              A=-0.5d0*(YP(pny,1)+DPDN/XNFV(IDIST,bnfv))*
     '          XNFV(FAREA,bnfv)
              DO nhx=1,nh_loc(0,nx)-1
                VELCOR(nhx)=A*XNFV(nhx,bnfv)
              ENDDO

C             ..Velocity fixed - inlet
            ELSEIF(NVCB(BCTYPE,bnvc).EQ.INLET) THEN
              A=0.5d0*YP(pny,1)*XNFV(FAREA,bnfv)
              DO nhx=1,nh_loc(0,nx)-1
                VELCOR(nhx)=-A*XNFV(nhx,bnfv)
              ENDDO

C             ..End of boundary condition "case" statement
            ENDIF

C           ..Divide velocities through by the volume
            DO nhx=1,nh_loc(0,nx)-1
              VELCOR(nhx)=VELCOR(nhx)/VC(nvc)
            ENDDO

C           ..Correct the velocities
            DO nhx=1,nh_loc(0,nx)-1
              nh=nh_loc(nhx,nx)
              ny=NYNP(1,1,nh,np,0,1,nr)
              YP(ny,1)=YP(ny,1)+DT*VELCOR(nhx)
            ENDDO
          ENDIF !Boundary
        ENDDO ! nfvl
      ENDDO

C     Store velocities for next iteration & calculate field change
      DV=0.d0
      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)
        DO nhx=1,nh_loc(0,nx)-1
          nh=nh_loc(nhx,nx)
          ny=NYNP(1,1,nh,np,0,1,nr)
          DV=DV+(YP(ny,1)-YP(ny,8))**2
          YP(ny,12)=YP(ny,8)
          YP(ny,8)=YP(ny,1)
        ENDDO
      ENDDO
      DV=DSQRT(DV)/DT

      CALL EXITS('CORRVELO')
      RETURN
 9999 CALL ERRORS('CORRVELO',ERROR)
      CALL EXITS('CORRVELO')
      RETURN 1
      END


