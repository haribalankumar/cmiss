      SUBROUTINE CALC_HEMODYNAMICS(ITYP3,nb,NEELEM,NENP,NPNE,NVJE,
     '  CE,XP,ERROR,*)

C#### Subroutine: CALC_HEMODYNAMICS
C###  Description:
C###    CALC_HEMODYNAMICS calculates hemodynamic properties to
C###    be used when solving for the blood flow in the
C###    pulmonary microcirculation. These calculations only dependent
C###    on vessel dimensions, therefore changing time.

C***  Created by KSB, October 2001.
C***
C***  Calculates: Approximation for Reynolds number*friction factor.
C***  CE(nm_a,ne) and CE(nm_b,ne) stores the major & minor axes
C***  of each segment (element).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter List
      INTEGER ITYP3,nb,NEELEM(0:NE_R_M),NENP(NPM,0:NEPM),
     '  NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 CE(NMM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER noelem,ne,nj,nv1,nv2
      REAL*8 Refd,R_factor

      CALL ENTERS('CALC_HEMODYNAMICS',*9999)

      IF(ITYP3.EQ.3) THEN !only capillaries
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
C         ny_flow=NYNE(1,1,1,1,ne) !assumes na=1
          IF(CE(nm_Dh,ne).NE.0.d0) THEN
            Refd=80.2d0-30.3d0*(CE(nm_b,ne)/CE(nm_a,ne))+3.45d0*
     '        (CE(nm_b,ne)/CE(nm_a,ne))**2.0d0+10.6d0*
     '        (CE(nm_b,ne)/CE(nm_a,ne))**3.0d0
            R_factor=Refd/(2.0d0*PI*CE(nm_a,ne)*CE(nm_b,ne)
     &        *CE(nm_Dh,ne)**2.0d0)
            CALL CALC_R_SEGMENT(nb,ne,NENP,NPNE,CE,R_factor,XP,ERROR,
     &        *9999)
          ELSE
            CE(nm_Rseg,ne)=1.d10 !vessel collapse therefore zero flow so setting high resistance
            WRITE(OP_STRING,'('' Diameter=0 for element '',I6)') ne
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO        
      ELSEIF(ITYP3.EQ.6) THEN !couple arteriole-cap-venule problem
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          IF(CE(nm_Dh,ne).EQ.0.d0.AND.CE(nm_a0,ne).EQ.0.d0) THEN
            nv1=NVJE(1,nb,nj_radius,ne)
            nv2=NVJE(2,nb,nj_radius,ne)
            IF(nv1.EQ.0.OR.nv2.EQ.0) THEN
              CALL ASSERT(.FALSE.,'>>Diameter versions not set up',
     &          ERROR,*9999)
            ENDIF
            IF(XP(1,nv1,nj_radius,NPNE(1,nb,ne)).EQ.0.d0.AND.
     &        XP(1,nv2,nj_radius,NPNE(2,nb,ne)).EQ.0.d0) THEN
              CALL ASSERT(.FALSE.,'>>Diameters not set up',ERROR,*9999)
            ELSE
              CE(nm_Dh,ne)=(XP(1,nv1,nj_radius,NPNE(1,nb,ne))+
     &          XP(1,nv2,nj_radius,NPNE(2,nb,ne))) !diameter
            ENDIF
          ENDIF
          nv1=NVJE(1,nb,nj_radius,ne)
          IF(XP(1,nv1,nj_radius,NPNE(1,nb,ne)).NE.0.d0) THEN
            IF(CE(nm_length,ne).EQ.0.d0) THEN
              DO nj=1,NJT
                CE(nm_length,ne)=CE(nm_length,ne)+
     &            (XP(1,1,nj,NPNE(2,nb,ne))-XP(1,1,nj,NPNE(1,nb,ne)))
     &            **2.d0
              ENDDO
              CE(nm_length,ne)=DSQRT(CE(nm_length,ne))
            ENDIF
            CALL POISEUILLE_R(nb,ne,NPNE,NVJE,CE,XP,ERROR,*9999)
          ELSE !alveolar-capillary vessel
            IF(CE(nm_a,ne).NE.0.d0.AND.CE(nm_b,ne).NE.0.d0) THEN
              Refd=80.2d0-30.3d0*(CE(nm_b,ne)/CE(nm_a,ne))+3.45d0*
     '          (CE(nm_b,ne)/CE(nm_a,ne))**2.0d0+10.6d0*
     '          (CE(nm_b,ne)/CE(nm_a,ne))**3.0d0
              R_factor=Refd/(2.0d0*PI*CE(nm_a,ne)*CE(nm_b,ne)
     &          *CE(nm_Dh,ne)**2.0d0)
              CALL CALC_R_SEGMENT(nb,ne,NENP,NPNE,CE,R_factor,XP,ERROR,
     &          *9999)
            ELSE !vessel collapsed
              CE(nm_Rseg,ne)=1.d10 !vessel collapse therefore zero flow so setting high resistance
              WRITE(OP_STRING,'('' Diameter=0 for element '',I6)') ne
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF !capillary vessel collapsed
          ENDIF
        ENDDO !noelem
      ENDIF !ITYP3      
            
      CALL EXITS('CALC_HEMODYNAMICS')
      RETURN
 9999 CALL ERRORS('CALC_HEMODYNAMICS',ERROR)
      CALL EXITS('CALC_HEMODYNAMICS')
      RETURN 1
      END



