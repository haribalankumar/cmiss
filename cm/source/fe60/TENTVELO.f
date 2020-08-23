      SUBROUTINE TENTVELO(NODENVC,NPNODE,nr,nx,NYNP,OLDVC,VC,YP,ERROR,*)

C#### Subroutine: TENTVELO
C###  Description:
C###    <HTML>
C###    This subroutine calculates the predicted velocities, combining
C###    the old and new velocities to estimate new ones.
C###    Calculated by:
C###    <PRE>
C###
C###    U'_(i) = (Vn(i)*Un_(i) + dt*NettFlux_(i))/Vn+1(i)
C###
C###    Where  U'       = momentum marched non-solenoidal estimates of
C###                      velocity
C###           Vn       = velocity at time step n
C###           Vn+1     = Velocity at time step n+1
C###           dt       = time step
C###           NettFlux = Nett flux of momentum into cell
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b12.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER NPNODE(0:NP_R_M,0:NRM),nr,NODENVC(NVCM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 OLDVC(0:NVCM),VC(0:NVCM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nonode,np,nhx,nh,nvc,ny
      REAL*8 INV_VOL

      CALL ENTERS('TENTVELO',*9999)

      IF(MESHFIXD) THEN
        DO nvc=1,NVCT
          nonode=NODENVC(nvc)
          np=NPNODE(nonode,nr)
          INV_VOL=1.d0/VC(nvc)
          DO nhx=1,nh_loc(0,nx)-1
            nh=nh_loc(nhx,nx)
            ny=NYNP(1,1,nh,np,0,1,nr)
            YP(ny,1)=YP(ny,8)+DT*INV_VOL*YP(ny,1)
          ENDDO
        ENDDO
      ELSE !Mesh moving
        DO nvc=1,NVCT
          nonode=NODENVC(nvc)
          np=NPNODE(nonode,nr)
          INV_VOL=1.d0/VC(nvc)
          DO nhx=1,nh_loc(0,nx)-1
            nh=nh_loc(nhx,nx)
            ny=NYNP(1,1,nh,np,0,1,nr)
            YP(ny,1)=(OLDVC(nvc)*YP(ny,8)+DT*YP(ny,1))*INV_VOL
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('TENTVELO')
      RETURN
 9999 CALL ERRORS('TENTVELO',ERROR)
      CALL EXITS('TENTVELO')
      RETURN 1
      END


