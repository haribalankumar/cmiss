      SUBROUTINE MESH_FLOW(NBJ,NEELEM,NPNE,NVJE,nx,NXI,XP,ERROR,*)

C#### Subroutine: MESH_FLOW
C###  Description:
C###    MESH_FLOW calculates flow in 1D mesh.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'moti00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),
     '  NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  nx,NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,ne0,nn,noelem,np,np2,nv
      REAL*8 Re,rho,visc,radius,radius_trachea
      CHARACTER STRING*255

      CALL ENTERS('MESH_FLOW',*9999)

      ne_trachea=NEELEM(1)
      nb=NBJ(1,ne_trachea)
      radius_trachea=XP(1,1,nj_radius,NPNE(1,nb,ne_trachea))
      WRITE(STRING,'(''>>Must specify field for coefficient''
     &  '' in DEFINE MOTION;c command'')')
      CALL ASSERT(nj_coeff.GT.0,STRING,ERROR,*9999)
      visc=1.82d-5 !viscosity of air (g/mm/s)
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(nj_radius,ne)
        np2=NPNE(2,nb,ne)
        radius=0.d0
        DO nn=1,2
          np=NPNE(nn,nb,ne)
          nv=NVJE(nn,nb,nj_radius,ne)
 !temporary calculation for radius
          radius=radius+0.5d0*XP(1,nv,nj_radius,np)
        ENDDO !nn
        ne0=NXI(-1,1,ne)
        rho=PULMAT(3)
        Re=DABS(INLET_FLOW(nx)*XP(1,1,nj_flow,np2)/(PI*radius**2.d0))*
     &    2.d0*radius/(visc/rho)
 !coefficient is proportional to Re and roughness
 !ET tube is smooth, airways 'rough' proportional to trachea (ratio of radius)
        DO nn=1,2
          np=NPNE(nn,nb,ne)
          nv=NVJE(nn,nb,nj_coeff,ne)
c          XP(1,nv,nj_coeff,np)=Re/150.d0 !/(radius_trachea/radius)
          XP(1,nv,nj_coeff,np)=Re/50.d0 !/(radius_trachea/radius)
c          IF(FLOW.GT.0.d0)THEN
c            XP(1,nv,nj_coeff,np)=MIN(XP(1,nv,nj_coeff,np)*3.d0,40.d0
c     &        *radius/radius_trachea)
c          ELSE
c            XP(1,nv,nj_coeff,np)=MIN(XP(1,nv,nj_coeff,np)*3.d0,120.d0
c     &        *radius/radius_trachea)
c          ENDIF
c          XP(1,nv,nj_coeff,np)=MAX(XP(1,nv,nj_coeff,np),2.d0)

          IF(INLET_FLOW(nx).GE.0.d0)THEN
            XP(1,nv,nj_coeff,np)=MIN(XP(1,nv,nj_coeff,np)*WHT(5),60.d0)
          ELSE
            XP(1,nv,nj_coeff,np)=MIN(XP(1,nv,nj_coeff,np)*WHT(6),60.d0)
          ENDIF
c          XP(1,nv,nj_coeff,np)=MAX(XP(1,nv,nj_coeff,np),2.d0)

        ENDDO !nn
      ENDDO !noelem

      CALL EXITS('MESH_FLOW')
      RETURN

 9999 CALL ERRORS('MESH_FLOW',ERROR)
      CALL EXITS('MESH_FLOW')
      RETURN 1
      END



