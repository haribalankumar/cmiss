      SUBROUTINE CALC_R_SEGMENT(nb,ne,NENP,NPNE,CE,R_factor,XP,ERROR,*)

C#### Subroutine: CALC_R_SEGMENT
C###  Description:
C###   CALC_R_SEGMENT calculates resistance in each capillary segment.
C###   The apparent viscosity of blood in the pulmonary microcirculation
C###   is first calculated, taking into account blood hematocrit
C###   in each vessel.
C***  Resistance stored in CE(nm_Rseg,ne).

C***  Created by KSB, October 2001
C***  Details for calculation from Huang et. al, 2001.
C***  NB/ all u's represent mu's (i.e. viscosity terms).
C*** For the element which the resistance is being calculated, the
C*** resistance is proportional to the effective surface
C*** area which is the total area (length*perimeter) less
C*** the surface area of each overlapping triangular wedge.
C*** Effective S.A=perimeter*(total ne length-1/4*sum of overlap length)

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter List
      INTEGER nb,ne,NENP(NPM,0:NEPM),NPNE(NNM,NBFM,NEM)
      REAL*8 CE(NMM,NEM),R_factor,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ne2,nj,noelem,noelem2,np0,np1,np2,nn
      REAL*8 A(3),ANG_BETWEEN,ANG3P,B(3),BETA,C(3),Dm,D_star,LEN_W,
     &  LENGTH,LENGTH_SUM,MIN_LENGTH,perimeter,SA_JUNCTION,TOTAL_SA,
     &  u_app,u_cyto,u_plasma

      CALL ENTERS('CALC_R_SEGMENT',*9999)

      u_plasma=1.2d0*1.d-3 !units (cP->Pa.s)
      Dm=0.00270d0 !(mm)=diam of smallest vessel RBC can pass through
      u_cyto=exp(0.48d0+2.35d0*CE(nm_Hd,ne))*1.0d-3 !cP->Pa.s (RBC)
      D_star=(2.03d0-2.0d0*CE(nm_Hd,ne))*1.0d-3 !um -> mm
      u_app=u_plasma/((1.d0-((1.d0-u_plasma/u_cyto)*(1.d0-
     '  (2.d0*D_star/CE(nm_Dh,ne)))**4.d0))*(1.d0-(Dm/CE(nm_Dh,ne))
     &  **4.d0))
      CE(nm_Rseg,ne)=CE(nm_length,ne)*u_app*R_factor
C..  !units Pa.s/mm^3
      LENGTH_SUM=0.d0
      DO nn=1,NNT(nb) !each node of an element
        np0=NPNE(nn,nb,ne)
        np1=NPNE(3-nn,nb,ne) !other node of ne
        DO noelem2=1,NENP(np0,0)
          ne2=NENP(np0,noelem2)
          IF(ne.NE.ne2) THEN !only evaluate connecting elements
            np2=NPNE(1,nb,ne2)
            IF(np2.EQ.np0) np2=NPNE(2,nb,ne2)
            DO nj=1,3
              A(nj)=XP(1,1,nj,np1) !stores XP co-ords
              B(nj)=XP(1,1,nj,np0) !for ANG_BETWEEN calculation
              C(nj)=XP(1,1,nj,np2)
            ENDDO !nj
            ANG_BETWEEN=ANG3P(A,B,C) !angle btwn 2 elems (ne & ne2)
            BETA=PI-ANG_BETWEEN
            MIN_LENGTH=(0.5d0*CE(nm_a,ne2))/DCOS(BETA)
            IF(MIN_LENGTH.LE.(0.5d0*CE(nm_a,ne)).AND.
     '        BETA.NE.(0.5d0*PI).OR.BETA.EQ.0.d0) THEN
              LENGTH=0.d0
            ELSE IF(BETA.EQ.(0.5d0*PI)) THEN !elements orthogonal
              LENGTH=CE(nm_a,ne2)*0.5d0
            ELSE
              LEN_W=DSQRT((0.5d0*CE(nm_a,ne))**2.d0+
     '          (0.5d0*CE(nm_a,ne2))**2.d0-2.d0*(0.5d0*CE(nm_a,ne))
     '          *(0.5d0*CE(nm_a,ne2))*DCOS(BETA)) !cos rule
              LENGTH=LEN_W/DSIN(ANG_BETWEEN)*DSIN(PI*0.5d0- !sine
     '          (ASIN(DSIN(BETA)*0.5d0*CE(nm_a,ne2)/LEN_W))) !rule
C!! LENGTH=length of ne which is overlapping into ne2
            ENDIF
            LENGTH_SUM=LENGTH_SUM+LENGTH !sum of overlap for ne
          ENDIF !ne.NE.ne2
        ENDDO
      ENDDO !ns
C... Calculation of surface areas & recalculation of resistance for ne
      perimeter=2.d0*PI*DSQRT(0.5d0*(CE(nm_a,ne)**2.d0+
     '  CE(nm_b,ne)**2.d0))
      TOTAL_SA=CE(nm_length,ne)*perimeter
      SA_JUNCTION=0.25d0*LENGTH_SUM*perimeter
      IF(SA_JUNCTION.GE.TOTAL_SA) THEN
        CE(nm_Rseg,ne)=1.d0 !effectively no resistance in segment
      ELSE
        CE(nm_Rseg,ne)=(1.d0-(SA_JUNCTION/TOTAL_SA))*CE(nm_Rseg,ne)
      ENDIF
      DO noelem=1,NE_BLOCK(0) !segments blocked with WBCs
        ne2=NE_BLOCK(noelem)
        IF(ne2.EQ.ne) CE(nm_Rseg,ne2)=1.d10 !increase resistance to high level
      ENDDO !noelem
      
      CALL EXITS('CALC_R_SEGMENT')
      RETURN
 9999 CALL ERRORS('CALC_R_SEGMENT',ERROR)
      CALL EXITS('CALC_R_SEGMENT')
      RETURN 1
      END


