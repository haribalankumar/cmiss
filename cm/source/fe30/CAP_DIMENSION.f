      SUBROUTINE CAP_DIMENSION(ITYP3,KOUNT,KOUNT_NE,nb,NEELEM,NENP,NPNE,
     &  NVJE,NYNP,CE,XP,YP,ERROR,*)

C#### Subroutine: CAP_DIMENSION
C###  Description:
C###    CAP_DIMENSION calculates the dimensional changes, resulting
C###    from changes in transmural pressure (Ptm=P_vessel-P_alveolar)
C###    and transpulmonary pressure (Ptp=lung volume=P_alveolar-
C###    P_intrapleural), in the pulmonary capillary network.

C***   Method from Huang et al., 2001, equations 1-13.
C***   Created by Kelly Burrowes, March 2002.

C*** Will keep all length in micrometers & pressures in cm H2O for
C*** calculation then convert back to mm & Pa.

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter list
      INTEGER ITYP3,KOUNT,KOUNT_NE(NEM),nb,NEELEM(0:NE_R_M),
     &  NENP(NPM,0:NEPM),NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM)
      REAL*8 CE(NMM,NEM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER count,i,IT,ne,ne2,ne_adj,nn,noelem,noelem2,np1,
     '  ny1,nv1,nv2,ny2
      REAL*8 A,a1,a2,B,C,c1,c2,c3,c4,d,D01,D02,dx,E,E_star,f,fmid,h,
     &  L,Lo,L_Lo,M,P_atm,P_vessel,Ptm,Ptp,perimeter,
     &  radius,rtbis,theta1,theta2,tol,Vo
      CHARACTER CHAR*5,STRING*255
      LOGICAL FEED,FOUND

      CALL ENTERS('CAP_DIMENSION',*9999)

C.....Setting default values for parameters to ensure initialized
      A=2.90d5 !um**3
      B=6.43d0 !cm H2O
      M=2.14d0
      
      count=0
      tol=1.0d-6 !zero convergence tolerance
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        FEED=.FALSE.
        IF(ITYP3.EQ.3) THEN
          DO i=1,INLETS(0)
            IF(ne.EQ.INLETS(i)) FEED=.TRUE.
          ENDDO
          DO I=1,OUTLETS(0)
            IF(ne.EQ.OUTLETS(i)) FEED=.TRUE.
          ENDDO
        ENDIF !only exclude feed vessels for pure capillary model
        IF(.NOT.FEED) THEN !only do capillary vessels not feed vessels
 !this applies a damping effect - if vessel is collapsed it cannot reopen straight away -
 !waits for 5 iterations before it is allowed to open.
          IF(CE(nm_Dh,ne).EQ.0.d0.AND.KOUNT_NE(ne).LT.5) THEN
            KOUNT_NE(ne)=KOUNT_NE(ne)+1 
          ELSE
            KOUNT_NE(ne)=0
            ny1=NYNP(1,1,1,NPNE(1,nb,ne),0,1) !nk=nv=nh=nc=1
            ny2=NYNP(1,1,1,NPNE(2,nb,ne),0,1)
            D01=0.d0
            IF(ITYP3.EQ.6) THEN
              nv1=NVJE(1,nb,nj_radius,ne)
              nv2=NVJE(2,nb,nj_radius,ne)
              D01=XP(1,nv1,nj_radius,NPNE(1,nb,ne))*2.d0 !diameter
              D02=XP(1,nv2,nj_radius,NPNE(2,nb,ne))*2.d0
            ELSE
              D02=0.d0
            ENDIF
C...  average pressure in capillary segment
            P_vessel=((YP(ny1,1)+YP(ny2,1))/2.d0)*0.0102d0 !Pa -> cm H20
C           IF(CE(nm_Dh,ne).LT.0.011d0.OR.ITYP3.EQ.3) THEN !alveolar capillary vessel
            IF(D01.EQ.0.d0) THEN
              IF(KOUNT.LT.3.AND.P_vessel.LT.0.d0) then
                P_vessel=1.d0 !make pressure +ve ->KSB 7/7/05 so can set Pin & Qin
              ENDIF
              P_atm=0.0d0 !atmospheric defined as 0 cm H2O
              Ptm=P_vessel-P_alveolar !should be 10-14 cm H2O
C... During inspiration, intrapleural pressure decreases, alveolar
C... pressure should be calculated from this & passed into this
C... subroutine.
              Ptp=P_alveolar-P_pleural !cm H2O
              
!NB! 21st JUNE for PPV simulations -> +ve alveolar pressure corresponds to inspiration!   
              IF(P_alveolar.LT.P_atm) THEN !inflation
                A=2.90d5 !um**3
                B=6.43d0 !cm H2O (--------> This done 15 August 2005 for PPV simulations !)
                M=2.14d0
              ELSEIF(P_alveolar.GE.P_atm) THEN !deflation and at rest
                A=3.07d5 !um**3
                B=3.42d0 !cm H2O
                M=1.03d0
              ENDIF
              
              Vo=20.7d3 !um**3 Inflation pressure=0 cm H2O (Mercer:1987a) rat
C             Vo=5.44d6 !alveolar volume from Weibel:1963, pg 68
              L_Lo=(1.d0+(A/Vo-1.d0)*EXP(-(B/Ptp)**M))**(1.d0/3.d0) !dimensionless
              Lo=0.07d0 !mm from Huang, 2001. test sensitivity
              L=L_Lo*Lo !mm
              E=(0.5d0*Ptp*L)/(L_Lo-1.d0) !cm H2O*mm, elastic coefficient
C... the following calculates d (distance btwn capillary segments)
C... d is approximated by average length of all ne's adjacent to ne
              d=0.d0
              ne_adj=0
              DO nn=1,NNT(nb)
                np1=NPNE(nn,nb,ne)
                DO noelem2=1,NENP(np1,0)
                  ne2=NENP(np1,noelem2)
                  IF(ne2.NE.ne) THEN
                    d=d+CE(nm_length,ne2) !sums lengths of adjacent caps
                    ne_adj=ne_adj+1 !counts # of adjacent elements
                  ENDIF
                ENDDO !noelem2
              ENDDO !nn
              d=(d/ne_adj)*2.d0 !mm
C... test senstitivity of variations in d, Vo and Lo
              IF(CE(nm_a0,ne).EQ.0.d0) THEN
C               WRITE(*,*),"Capillary diameter",CE(nm_Dh,ne)
                WRITE(CHAR,'(I5)') ne
                WRITE(STRING,'(''>>Capillary diam values=0:'//CHAR//
     &            ' '')')
                CALL ASSERT(.FALSE.,STRING,ERROR,*9999)            
              ENDIF
              a1=CE(nm_a0,ne)*L_Lo !(mm) length of septal part @ Ptp, Ptm=0
              E_star=E*((1.d0/a1)+(1.d0/(d-a1))) !cm H20
              c1=2.d0*a1*E_star !(cm H2O*mm) coefficients to determine
              c2=-CE(nm_C0,ne)*Ptm !mm the angle, theta, (see equation
              c3=-2.d0*CE(nm_C0,ne)*E_star !13 Huang et al, 2001).
C             kc=PULMAT(3)
              kc=CE(3,ne) !material properties now stored in CE
              c4=-(a1*E_star*Ptm)/kc !cm H2O*um
C... PULMAT(3)=kc= (stiffness of capillary wall (cm H2O)
C... find theta using bisection method:
              theta1=-2.d0*PI
              theta2=2.d0*PI
              f=c1*theta1+c2*DCOS(theta1)+c3*DSIN(theta1)+c4
              fmid=c1*theta2+c2*DCOS(theta2)+c3*DSIN(theta2)+c4
              IF(f.LT.0.d0) THEN
                rtbis=theta1
                dx=theta2-theta1
              ELSE
                rtbis=theta2
                dx=theta1-theta2
              ENDIF
              FOUND=.FALSE.
              IT=0 !# iterations
              DO WHILE(.NOT.FOUND.AND.IT.LT.150) !Bisection loop
                IT=IT+1
                dx=dx*0.5d0
                theta2=rtbis+dx
                fmid=c1*theta2+c2*DCOS(theta2)+c3*DSIN(theta2)+c4
                IF(fmid.LE.0.d0) rtbis=theta2
                IF(fmid.GT.0.d0-tol.AND.fmid.LT.0.d0+tol) FOUND=.TRUE.
              ENDDO !WHILE
              IF(.NOT.FOUND) THEN
                WRITE(OP_STRING,'('' Theta not found for ne #: '',I5)')
     &            ne
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ENDIF

!              WRITE(OP_STRING,'('' test cap_dimen '',2(D12.4,X))')
!     &          P_vessel,theta2
!              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              
              IF(theta2.LE.0.d0.OR..NOT.FOUND) THEN !Diameter of vessel=0
                count=count+1
                theta2=0.d0
                CE(nm_Dh,ne)=0.d0 !vessel collapsed
                CE(nm_b,ne)=0.d0 
                CE(nm_a,ne)=0.d0
              ELSE
                radius=CE(nm_C0,ne)/(2.d0*theta2-(Ptm/kc)) !mm
                C=(2.d0*radius*theta2) !new C at Ptp and Ptm (um->mm)
                a2=(2.d0*radius*DSIN(theta2)) !length septal portion (mm)
                h=(radius*(1.d0-DCOS(theta2))) !mm
C               IF(theta2.LT.0.d0.AND.h.LT.0.d0) h=-h
                perimeter=C+a2 !mm
C.. perimeter of capillary used to find elliptical dimensions a* and b*
C.. elliptical dimensions used for solution
                CE(nm_b,ne)=h/2.d0 !mm          
                CE(nm_a,ne)=DSQRT(perimeter**2.d0/(2.d0*PI**2.d0)-
     '            CE(nm_b,ne)**2.d0) !mm
C               CE(nm_Ac,ne)=PI*CE(nm_b,ne)*CE(nm_a,ne) !mm**2, area of ellise
C               CE(nm_Dh,ne)=4.d0*CE(nm_Ac,ne)/perimeter !hydraulic diameter,mm
                CE(nm_Dh,ne)=4.d0*PI*CE(nm_b,ne)*CE(nm_a,ne)/perimeter !hydraulic diameter,mm
              ENDIF !theta.NOT.FOUND
C             IF(KOUNT.EQ.1) THEN !write some stuff out...
C             IF(theta2.LT.(PI/2.d0).OR.theta2.GT.PI) THEN
C             WRITE(OP_STRING,'('' theta,E,E* '',3(D10.4,X))') theta2,E,
C             &        E_star
C             CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C             ENDIF
            ELSE
C larger vessel: use compliance relationship: Diam/D0=1+alpha
C             *delta(Pressure)
C             nv1=NVJE(1,nb,nj_radius,ne)
C             nv2=NVJE(2,nb,nj_radius,ne)
C             D01=XP(1,nv1,nj_radius,NPNE(1,nb,ne))*2.d0 !diameter
C             D02=XP(1,nv2,nj_radius,NPNE(2,nb,ne))*2.d0 2
C             P1=YP(ny1,1)*0.0102d0 !Pa -> cm H20
C             P2=YP(ny2,1)*0.0102d0 
C             XP(1,nv1,nj_radius,NPNE(1,nb,ne))=D01+CE(3,ne)*P1
C             XP(1,nv2,nj_radius,NPNE(2,nb,ne))=D02+CE(3,ne)*P2
C             CE(nm_Dh,ne)=(XP(1,nv1,nj_radius,NPNE(1,nb,ne))+XP(1,nv2,
C             &      nj_radius,NPNE(2,nb,ne)))/2.d0
              IF(D01.EQ.0.d0.AND.D02.EQ.0.d0) THEN
C               WRITE(*,*),"DIAMETER",CE(nm_Dh,ne)
                WRITE(CHAR,'(I5)') ne
                WRITE(STRING,'(''>>Diameters=0 for element:'//CHAR//
     &            ' '')')
                CALL ASSERT(.FALSE.,STRING,ERROR,*9999)
              ENDIF
              CE(nm_Dh,ne)=((D01+D02)/2.d0)*(1.d0+CE(3,ne)*
     &          (P_vessel-P_pleural))
C             write(*,*),"diam",ne,D01,D02,P_vessel,CE(nm_Dh,ne)
 !average diameter over element
 !currently leaving XP unchanged as D0
            ENDIF
          ENDIF
          
        ENDIF !NOT.FEED
      ENDDO !noelem
      
      CALL EXITS('CAP_DIMENSION')
      RETURN
 9999 CALL ERRORS('CAP_DIMENSION',ERROR)
      CALL EXITS('CAP_DIMENSION')
      RETURN 1
      END


