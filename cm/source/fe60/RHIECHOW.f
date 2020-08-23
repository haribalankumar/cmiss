      SUBROUTINE RHIECHOW(NFVC,NODENVC,NPNODE,nr,NVCB,
     '  NVCNODE,nx,NYNP,XNFV,YP,ZNFV,ZNFVMSH,ERROR,*)

C#### Subroutine: RHIECHOW
C###  Description:
C###    <HTML>
C###    Calculates the scalar field of solenoidal Rhie Chow fluxes
C###    <PRE>
C###
C###    RCFlux(i,j) = A(i,j) *
C###                  {0.5*[U'_(i) + U'_(j)].N_(i,j) -
C###                   dt*[1/D(i,j)]*[P(j)-P(i)]}
C###
C###    Where A  = Voronoi face area
C###          U' = momentum marched non-solenoidal estimates of
C###               velocity
C###          N  = unit outward normal from cell i
C###          dt = time step
C###          D  = internodal distance
C###          P  = pressure
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVCB(-1:3,NVCBM),NVCNODE(2,NP_R_M),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 XNFV(-(NJM+1):NJM,NFVM),YP(NYM,NIYM),ZNFV(NFVM),
     '  ZNFVMSH(NFVM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nvc,nfvl,cnonode,nfv,nonode,np,cnp,nh,nhx,ny,cny,
     '  bnp,bnvc,bny,i,nfvl2,cnonode2,nfv2
      REAL*8 DIV,LINCOMP
      LOGICAL DIRECT_OPP

      CALL ENTERS('RHIECHOW',*9999)

C     ..Calculate internal fluxes..

      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)

C       ..Loop over the faces of Voronoi cell nvc to get adjoining nodes
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)

C         ..Make sure that the adjoining node is not a boundary node
          IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
            cnp=NPNODE(cnonode,nr)
            nfv=NFVC(2,nfvl,nvc)
            LINCOMP=0.d0

C           ..Compute [U'_(i) + U'_(j)].N_(i,j)
            DO nhx=1,nh_loc(0,nx)-1
              nh=NH_LOC(nhx,nx)
              ny=NYNP(1,1,nh,np,0,1,nr)
              cny=NYNP(1,1,nh,cnp,0,1,nr)
              LINCOMP=LINCOMP+(YP(ny,1)+YP(cny,1))*XNFV(nhx,nfv)
            ENDDO

C           ..Multiply this by 0.5*A(i,j) to get the linear component
            LINCOMP=0.5d0*XNFV(FAREA,nfv)*LINCOMP

C           ..Add on the pressure term A(i,j)*dt*[1/D(i,j)]*[P(j)-P(i)]
            ny=NYNP(1,1,nh_loc(0,nx),np,0,1,nr)
            cny=NYNP(1,1,nh_loc(0,nx),cnp,0,1,nr)
            ZNFV(nfv)=LINCOMP-DT*(YP(cny,1)-YP(ny,1))*
     '        XNFV(FAREA,nfv)*XNFV(IDIST,nfv)
          ENDIF
        ENDDO
      ENDDO

C     ..Calculate boundary fluxes..

      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)

C       ..Loop over the faces of cell nvc to get adjoining nodes
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)

C         ..Make sure the adjoining node is a boundary node
          IF(NVCNODE(TYPE,cnonode).EQ.BOUNDARY) THEN
            bnp=NPNODE(cnonode,nr)
            bnvc=NVCNODE(MAP,cnonode)
            nfv=NFVC(2,nfvl,nvc)

C           ..Inlets, noting that dp/dn = 0
            IF(NVCB(BCTYPE,bnvc).EQ.INLET) THEN

C             ..Compute [U'_(i) + U'_(j)].N_(i,j)
              LINCOMP=0.d0
              DO nhx=1,nh_loc(0,nx)-1
                nh=NH_LOC(nhx,nx)
                ny=NYNP(1,1,nh,np,0,1,nr)
                bny=NYNP(1,1,nh,bnp,0,1,nr)
                LINCOMP=LINCOMP+(YP(ny,1)+YP(bny,1))*XNFV(nhx,nfv)
              ENDDO

C             ..Multiply this by 0.5*A(i,j) to get the linear part
              ZNFV(nfv)=0.5d0*XNFV(FAREA,nfv)*LINCOMP

C             ..Outlets
            ELSEIF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN

C             ..Ensure that only the directly adjoining node is used
              DIRECT_OPP=.FALSE.
              DO i=1,NVCB(0,bnvc)
                IF(nonode.EQ.NVCB(i,bnvc)) THEN
                  DIRECT_OPP=.TRUE.
                  GOTO 100
                ENDIF
              ENDDO
 100          CONTINUE
              IF(DIRECT_OPP) THEN

C               ..Calculate flux through internal face walls
                DIV=0.d0
                DO nfvl2=1,NFVC(1,0,nvc)
                  cnonode2=NFVC(1,nfvl2,nvc)

C                 ..Make sure the adjoining node is an internal node
                  IF(NVCNODE(TYPE,cnonode2).NE.BOUNDARY) THEN
                    nfv2=NFVC(2,nfvl2,nvc)
                    DIV=DIV+ZNFV(nfv2)
                  ENDIF
                ENDDO

C               ..The flux through the outlet cell face will be equal
C               & opposite to the flux through the internal cell faces
                ZNFV(nfv)=-DIV
              ELSE
                ZNFV(nfv)=0.d0
              ENDIF
            ELSEIF(.NOT.MESHFIXD) THEN
              ZNFV(nfv)=ZNFVMSH(nfv)
            ENDIF
          ENDIF
        ENDDO

      ENDDO

C     ..Output variables
      INFLUX=0.d0
      OUTFLUX=0.d0
      VSWEPT=0.d0
      DIVMAG=0.d0
      DO nvc=1,NVCT
        DIV=0.d0
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)
          nfv=NFVC(2,nfvl,nvc)
          IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
            DIV=DIV+ZNFV(nfv)
          ELSE
            bnvc=NVCNODE(MAP,cnonode)
            IF(NVCB(BCTYPE,bnvc).EQ.INLET) THEN
              DIV=DIV+ZNFV(nfv)
              INFLUX=INFLUX-ZNFV(nfv)
            ELSEIF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN
              DIV=DIV+ZNFV(nfv)
              OUTFLUX=OUTFLUX-ZNFV(nfv)
            ENDIF
            IF(.NOT.MESHFIXD) THEN
              DIV=DIV+ZNFVMSH(nfv)
              VSWEPT=VSWEPT+ZNFVMSH(nfv)
            ENDIF
          ENDIF
        ENDDO
        DIVMAG=DIVMAG+DIV**2
      ENDDO
      DIVMAG=DSQRT(DIVMAG)
      OUTFLUX=-1.d0*OUTFLUX*DT
      INFLUX=-1.d0*INFLUX*DT
      VSWEPT=VSWEPT*DT

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(RHIECHOW_1)
        WRITE(OP_STRING,'(/$,'' Face flux listing:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ##################'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nvc=1,NVCT
          WRITE(OP_STRING,'('' Cell:      '',I12,'' (np = '',
     '      I5,'')'')') nvc,NPNODE(NODENVC(nvc),nr)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Opp node:  '',60I12)')
     '      (NPNODE(NFVC(1,nfvl,nvc),nr),nfvl=1,NFVC(1,0,nvc))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Face no.:  '',60I12)')
     '      (NFVC(2,nfvl,nvc),nfvl=1,NFVC(1,0,nvc))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Face flux: '',60D12.4)')
     '      (ZNFV(NFVC(2,nfvl,nvc)),nfvl=1,NFVC(1,0,nvc))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          IF(.NOT.MESHFIXD) THEN
            WRITE(OP_STRING,'('' Mesh flux: '',60D12.4)')
     '        (ZNFVMSH(NFVC(2,nfvl,nvc)),nfvl=1,NFVC(1,0,nvc))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !nvc
CC$OMP END CRITICAL(RHIECHOW_1)
      ENDIF !dop

      CALL EXITS('RHIECHOW')
      RETURN
 9999 CALL ERRORS('RHIECHOW',ERROR)
      CALL EXITS('RHIECHOW')
      RETURN 1
      END


