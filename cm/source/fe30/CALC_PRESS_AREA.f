      SUBROUTINE CALC_PRESS_AREA(FIRST_IT,NBJ,NEELEM,NORD,NPNE,NVJE,
     '  NVJP,NYNP,XAB,XP,YP,ERROR,*)

C#### Subroutine: CALC_PRESS_AREA
C###  Description:  Calculates new area based upon change in pressure
C#### From LAmbert et al
C###  Created by KLH May 2003
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'ptr00.cmn'

!     Parameter List
      INTEGER maxord,NBJ(NJM,NEM),NEELEM(0:NE_R_M),nb,NORD(5,NE_R_M),
     '  NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM)
      REAL*8 XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL FIRST_IT

!     Local Variables
      INTEGER gen,ne,nj,noelem,nn,np,np1,np2,nv,ny1,
     &  ny2,nps(2)
      REAL*8 alpha1(2),beta,beta2,Go,beta_vein,Go_vein,
     &  G_PLEURAL,HEIGHT(NJT),Pblood,Ppl,press1,
     '  press2,pl1,pl2,PLEURAL_DENSITY,Ptm,R0,maxpl,minpl,
     &  stretch,k_factor
      
      CALL ENTERS('CALC_PRESS_AREA',*9999)

      nb=NBJ(1,NEELEM(1))
      maxord=NORD(2,NEELEM(1))
      maxpl=0.d0
      minpl=7000.d0
      nps(1)=npne(1,nb,neelem(1))
      nps(2)=npne(1,nb,neelem(neelem(0)))

      IF (COMPLIANCE.NE.3.AND.COMPLIANCE.NE.6)THEN
        IF(FIRST_IT)then
          DO gen=1,24
            LAMB_P2(gen)=-nlam2(gen)*(1.d0-alpha0(gen))/dash0(gen)
            LAMB_P1(gen)=alpha0(gen)*nlam1(gen)/dash0(gen)
            alphastrt(gen)=1.d0-(1.d0-alpha0(gen))*(1.d0-30.d0/
     '        LAMB_P2(gen))**(-nlam2(gen))
          ENDDO
        ENDIF

        IF(COMPLIANCE.eq.1)THEN
          pl1=ptrans
          pl2=ptrans
        ENDIF
        IF(COMPLIANCE.EQ.5)THEN
          press1=ptrans
          press2=ptrans
        ENDIF
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          gen=nord(1,ne)
          IF(gen.gt.24)gen=24
          np1=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
          ny1=NYNP(1,1,1,np1,0,1)
          ny2=NYNP(1,1,1,np2,0,1)
          IF(COMPLIANCE.NE.5)THEN
            IF(FIRST_IT)then
              press1=pl1-0.d0
              press2=pl2-0.d0
            ELSE
              press1=pl1-YP(ny1,1)/98.07d0
              press2=pl2-YP(ny2,1)/98.07d0
            ENDIF
          ENDIF
          IF(press1.GE.0.0d0)THEN
            alpha1(1)=1.d0-(1.d0-alpha0(gen))*(1.d0-press1/
     '        LAMB_P2(gen))**(-nlam2(gen))
          ELSE
            alpha1(1)=alpha0(gen)*(1.d0-press1/
     '        LAMB_P1(gen))**(-nlam1(gen))
          ENDIF
          IF(press2.GE.0.0d0)THEN
            alpha1(2)=1.d0-(1.d0-alpha0(gen))*(1.d0-press2/
     '        LAMB_P2(gen))**(-nlam2(gen))
          ELSE
            alpha1(2)=alpha0(gen)*(1.d0-press2/
     '        LAMB_P1(gen))**(-nlam1(gen))
          ENDIF
          DO nn=1,2
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_radius,ne)
            XP(1,nv,nj_radius,np)=sqrt(XAB(nej_strtrad,ne)**2.0d0
     '        *alpha1(nn)/alphastrt(gen))
          ENDDO
        ENDDO
      ELSEIF(COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6) THEN !Blood vessel pressure-area relationship:
        IF(FIRST_IT)THEN ! AJS Sep 2007 store unstressed radius values (initial values)
          DO noelem=1,NEELEM(0)
            ne=NEELEM(noelem)
            nb=NBJ(nj_radius,ne)
            DO nn=1,NNT(nb)
              np=NPNE(nn,nb,ne)
              DO nv=1,NVJP(nj_radius,np) !Loops over each nodal radius field
                XP(1,nv,nj_radius_R0,np)=XP(1,nv,nj_radius,np)
              ENDDO
            ENDDO !nn
          ENDDO !noelem
        ENDIF
        Go=PULMAT(3)   !kPa 
        beta=PULMAT(4) !unitless 
        Go_vein=PULMAT(6)   !kPa 
        beta_vein=PULMAT(7) !unitless 

C... KSB: Using inlet node as gravitational reference height        
        CALL ASSERT(np_in.NE.0,
     &    '>>Inlet node not defined: use INLET_REF_NODE [np#]',ERROR,
     &    *9999)        
        PLEURAL_DENSITY=0.25d0*0.1d-05 !kg/mm**3
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          np1=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
          DO nn=1,2
            np=NPNE(nn,nb,ne)
            ny1=NYNP(1,1,1,np,0,1)
            IF(COMPLIANCE.EQ.3) THEN
C...  KSB: Adding gravity-dependent pleural pressure            
              DO nj=1,NJT
                HEIGHT(nj)=XP(1,1,nj,np)-XP(1,1,nj,np_in)
              ENDDO
              G_PLEURAL=0.d0    !gravitational force
              DO nj=1,NJT
                G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*grav_vect(nj)
     &            *gravfact*9810.d0*HEIGHT(nj) !kg
              ENDDO
!            ptrans is set up in INIT_2_11.f
              Ppl=(ptrans*98.07d0/1000.0d0)-G_PLEURAL !cmH2O->kPa
              IF(nj_pleural.NE.0) THEN
                XP(1,1,nj_pleural,np)=Ppl/0.09806d0 !Convert to cmH2O for viewing
              ELSE
                XP(1,1,7,np)=Ppl/0.09806d0 
              ENDIF
            ELSEIF(COMPLIANCE.EQ.6) THEN
              Ppl=XP(1,1,nj_pleural,np)*98.07d0/1000.0d0 !cmH2O->kPa
              IF(nj_stretch.EQ.0) THEN
                stretch=1.d0   !Set to 1 i.e. no effect if undefined
              ELSE
                 stretch=XP(1,1,nj_stretch,np)
              ENDIF
              beta2=PULMAT(5)
              !See Smith:2004, incl vessel stretch factor, lambda
              beta=beta2*stretch+PULMAT(4)
            ENDIF
            Pblood=YP(ny1,1)/1000.0d0 ! Pa->kPa
            Ptm=Pblood+Ppl     ! kPa 
            DO nv=1,NVJP(nj_radius,np) !Update each nodal radius version
              R0=XP(1,nv,nj_radius_R0,np)
              IF(COMPLIANCE.EQ.6) R0=R0*sqrt(1.0d0/stretch) 
                 IF(XAB(nej_cap,ne).EQ.1.d0)THEN ! its a capillary
                 ELSEIF(XAB(nej_cap,ne).EQ.2.d0)THEN ! its a vein
                   IF(Ptm.LT.Ptm_max.and.Go_vein.gt.0.d0)THEN
                    XP(1,nv,nj_radius,np)=R0*((Ptm/Go_vein)+1.d0)
     &              **(1.d0/beta_vein)
                   ELSEIF(Ptm.lt.0.or.Go_vein.eq.0.d0)THEN
          IF(Ptm.lt.0)write(*,*) 'Transmural pressure less than zero',ne
                      XP(1,nv,nj_radius,np)=R0
                   ELSE
                       XP(1,nv,nj_radius,np)=R0*((Ptm_max/Go_vein)+1.d0)
     &                    **(1.d0/beta_vein)
                    ENDIF
                 ELSE !its an artery or no distinction
C...ARC: giving a maximum distension
                   IF(Ptm.LT.Ptm_max.and.Go_vein.gt.0.d0)THEN
                    XP(1,nv,nj_radius,np)=R0*((Ptm/Go)+1.d0)
     &              **(1.d0/beta)
                   ELSEIF(Ptm.lt.0.or.Go.eq.0.d0)THEN
          IF(Ptm.lt.0)write(*,*) 'Transmural pressure less than zero',ne
                      XP(1,nv,nj_radius,np)=R0
                   ELSE
                       XP(1,nv,nj_radius,np)=R0*((Ptm_max/Go)+1.d0)
     &                    **(1.d0/beta)
               write(*,*) 'Transmural pressure less than zero',ne
                   ENDIF
               ENDIF
              IF(nj_hypoxia.NE.0) THEN
                IF(XP(1,nv,nj_radius,np).LE.0.25d0) THEN !If vessels < 500 um in diam then apply HPV
                 k_factor=XP(1,1,nj_hypoxia,np2) !Arterial constriction factor, NB k factor only defined in version 1
c                 if (k_factor.NE.1.d0) THEN
c                    write(*,*),"CHECK",ne,k_factor,XP(1,nv,nj_radius,np)
c     &                  ,nj_hypoxia
c                 endif
                 XP(1,nv,nj_radius,np)=XP(1,nv,nj_radius,np)*k_factor
               ENDIF
              ENDIF
            ENDDO !nv
          ENDDO                 !nn
       ENDDO                    !noelem
      ELSE
        WRITE(OP_STRING,'('' Compliance model not implemented '')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF


      CALL EXITS('CALC_PRESS_AREA')
      RETURN
 9999 CALL ERRORS('CALC_PRESS_AREA',ERROR)
      CALL EXITS('CALC_PRESS_AREA')
      RETURN 1
      END
      


