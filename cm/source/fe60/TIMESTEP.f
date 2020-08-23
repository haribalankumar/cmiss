      SUBROUTINE TIMESTEP(NFVC,NODENVC,NPNODE,nr,NVCNODE,nx,NYNP,XNFV,
     '  YP,ERROR,*)

C#### Subroutine: TIMESTEP
C###  Description:
C###    <HTML>
C###    Calculates the time step, using a Voronoi mesh, represented
C###    here by the topology arrays NFVC and XNFV. This subroutine
C###    sweeps through the Voronoi cells, and then uses an inner loop
C###    to calculate the advective and diffusive stability limits for
C###    each face of each Voronoi cell, using explicit time stepping.
C###    <PRE>
C###    timestep < |D(i,j)/Un_(i).N_(i,j)|  (advective stability) and
C###    timestep < Re*4*D(i,j)^2            (viscous stability)
C###
C###    Where D  = distance between node i and j
C###          Un = fluid velocity at time step n
C###          N  = unit outward normal from node i
C###          Re = Reynolds number
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVCNODE(2,NP_R_M),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 XNFV(-(NJM+1):NJM,NFVM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nvc,nonode,np,nfvl,cnonode,nfv,nhx,nh,ny
      REAL*8 INVTMAX,INVT,VELNRM,BVEL

      CALL ENTERS('TIMESTEP',*9999)
      INVTMAX=1.d-6

C     ..This section checks the advective stability criterion

      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)

C       ..Loop over the nodes adjoining the Voronoi cell
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)
          nfv=NFVC(2,nfvl,nvc)

C         ..Internal nodes only
          IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN

C           .Compute the normal component of velocity = (U.n)
            VELNRM=0.d0
            DO nhx=1,NH_LOC(0,nx)-1
              nh=nh_loc(nhx,nx)
              ny=NYNP(1,1,nh,np,0,1,nr)
              VELNRM=VELNRM+YP(ny,8)*XNFV(nhx,nfv)
            ENDDO

C           ..The Courant condition udt/dx <= 1
            INVT=DABS(XNFV(IDIST,nfv)*VELNRM)
            IF(INVT.GT.INVTMAX)THEN
              INVTMAX=INVT
            ENDIF

C         ..Otherwise a boundary connection
          ELSE

C           ..Just Velocity for boundaries
            BVEL=0.d0
            DO nhx=1,NH_LOC(0,nx)-1
              nh=nh_loc(nhx,nx)
              ny=NYNP(1,1,nh,np,0,1,nr)
              BVEL=BVEL+DABS(YP(ny,3))
            ENDDO

C           ..The Courant condition udt/dx <= 1
            INVT=DABS(XNFV(IDIST,nfv)*BVEL)
            IF(INVT.GT.INVTMAX)THEN
              INVTMAX=INVT
            ENDIF

          ENDIF
        ENDDO
      ENDDO

      DO nvc=1,NVCT

C       ..Loop over adjoining nodes
        DO nfvl=1,NFVC(1,0,nvc)

C         ..The Peclet number restriction
          nfv=NFVC(2,nfvl,nvc)
          INVT=4.d0*INV_REYNOLD*XNFV(IDIST,nfv)**2
          IF(INVT.GT.INVTMAX)THEN
            INVTMAX=INVT
          ENDIF

        ENDDO
      ENDDO

C     ..Scale time step by factor of safety
      DT=COURANT/INVTMAX

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/$,'' Time step DT = '',D14.6)') DT
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$     call mp_unsetlock()
      ENDIF

      CALL EXITS('TIMESTEP')
      RETURN
 9999 CALL ERRORS('TIMESTEP',ERROR)
      CALL EXITS('TIMESTEP')
      RETURN 1
      END


