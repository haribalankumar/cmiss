      SUBROUTINE POISEUILLE_R(nb,ne,NPNE,NVJE,CE,XP,ERROR,*)


C#### Subroutine: POISEUILLE_R
C###  Description:
C###    POISEUILLE_R calculates the resistance to flow assuming Poiseuille
C###    flow. Also uses an empirically derived model (Kiani:1991) to
C###    calculate the apparent viscosity in each segment based on
C###    vessel diameter and discharge hematocrit.

C*** Created by KSB, February 2004.
C***  This has been set up to obtain an approximate solution for flow
C***  in the pulmonary arteriole and venule networks. This resistance
C***  term is used if the vessel diameter is greater than 100 micrometers.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter list
      INTEGER nb,ne,NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 CE(NMM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nj,np1,np2,nv1,nv2
      REAL*8 Dm,D_star,length,u_app,u_cyto,u_plasma

      CALL ENTERS('POISEUILLE_R',*9999)
      
      u_plasma=1.2d0*1.d-3 !units (cP->Pa.s)
      Dm=0.00270d0 !(mm)=diam of smallest vessel RBC can pass through
        u_cyto=exp(0.48d0+2.35d0*CE(nm_Hd,ne))*1.0d-3 !cP->Pa.s (RBC)
        D_star=(2.03d0-2.0d0*CE(nm_Hd,ne))*1.0d-3 !um -> mm
        np1=NPNE(1,nb,ne)
        np2=NPNE(2,nb,ne)
        length=0.d0
        DO nj=1,NJT
          length=length+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2.d0
        ENDDO
        length=DSQRT(length) !length of element
        nv1=NVJE(1,nb,nj_radius,ne)
        nv2=NVJE(2,nb,nj_radius,ne)
        IF(CE(nm_Dh,ne).EQ.0.d0) THEN !element field not yet defined
          CE(nm_Dh,ne)=(XP(1,nv2,nj_radius,np2)+XP(1,nv1,nj_radius,np1))
     &      /2.d0 !this used in solution output procedure
        ENDIF
        u_app=u_plasma/((1.d0-((1.d0-u_plasma/u_cyto)*(1.d0-
     '    (2.d0*D_star/CE(nm_Dh,ne)))**4.d0))*(1.d0-(Dm/CE(nm_Dh,ne))))
        CE(nm_Rseg,ne)=(128.d0*u_app*length)/(PI*CE(nm_Dh,ne)**4.d0)
        !Poiseuille resistance term (Pa.s/mm**3)
        !currently neglecting overlap of elements at junctions
      
      CALL EXITS('POISEUILLE_R')
      RETURN
      
 9999 CALL ERRORS('POISEUILLE_R',ERROR)
      CALL EXITS('POISEUILLE_R')
      RETURN 1
      END

      
