      SUBROUTINE CALC_RESIS_FLOW(KOUNT,NBJ,NEELEM,NENP,NORD,NPNE,NYNE,
     &  NVJE,
     &  XAB,XP,YP,ERROR,*)

C#### Subroutine: CALC_RESIS_FLOW
C###  Description:  Calculates airway resistance only at the moment
C#### Equation from Anafi and Wilson,  J.Appl.Physiol,91:1185-1192,2001
C###  Created by KLH April 2003
C     Correction addded for turbulent flow in bifurcation.  See Pedley
C     et al.
C###  Updated ARC 01-2011 to include blood vessel resistance        
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter List
      INTEGER KOUNT,NBJ(NJM,NEM),NEELEM(0:NE_R_M),NENP(NPM,0:NEPM),
     &  NORD(5,NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),YP(NYM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER nb,noelem,ne,next,nj,np1,np2,nv1,nv2,ny
      REAL*8 avrad,area,reynolds,length,elength,blength,max_reynolds,
     &  ne_reynolds,zeta,check
      LOGICAL end1
      
      CALL ENTERS('CALC_RESIS_FLOW',*9999)

      max_reynolds=0.d0 !initialise - this used to find maximum Re # in system
      nb=NBJ(1,NEELEM(1))
      DO noelem=1,NEELEM(0)
        end1=.FALSE.
        length=0.d0
        ne=NEELEM(noelem)
        IF((COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6).AND.
     &     COUPLE_VIA_LPM.EQ.'Y'.AND.XAB(nej_cap,ne).EQ.1.d0)THEN
            XAB(nej_resis,ne)=INITIAL_LPM !Perfusion problems initial cap resistance is an input 
        ELSE
          np1=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
          nv1=NVJE(1,nb,nj_radius,ne)
          nv2=NVJE(2,nb,nj_radius,ne)
          avrad=0.5d0*(XP(1,nv1,nj_radius,np1)+XP(1,nv2,nj_radius,np2))
          area=pi*avrad**2.d0
          ny=NYNE(1,1,1,1,ne)
          reynolds=DABS(YP(ny))*avrad*2.d0*PULMAT(1)/(area*PULMAT(2))
          IF(reynolds.GT.max_reynolds) THEN
           max_reynolds=reynolds
           ne_reynolds=ne
          ENDIF
          DO nj=1,NJT
            length=length+(XP(1,1,nj,np1)-XP(1,1,nj,np2))**2.d0
          ENDDO
          elength=DSQRT(length)
          IF(SCALE_TURBULENCE.EQ.0.d0)THEN
            zeta=1.d0 !Allowing turbulence correction factor to be switched off
          ELSE
            length=0.d0
c Calculation of branch length to end of current element
            DO WHILE(.NOT.end1)
              IF(NENP(np1,0).EQ.2)THEN
                next=NENP(np1,1)
                np1=NPNE(1,nb,next)
              ELSE
                end1=.TRUE.
              ENDIF
            ENDDO
            DO nj=1,NJT
              length=length+(XP(1,1,nj,np1)-XP(1,1,nj,np2))**2.d0
            ENDDO
            blength=DSQRT(length)
            zeta=SQRT(2.d0*avrad*reynolds/blength)*1.85d0/(4.d0
     &        *SQRT(2.d0))
            zeta=zeta*SCALE_TURBULENCE !scale turbulence correction factor
            IF(zeta.LT.1.d0)zeta=1.d0
           ENDIF !SCALE_TUEBULENCE.EQ.0 
           XAB(nej_resis,ne)=128.0d0*elength*PULMAT(2)/
     &     (pi*(2.d0*avrad)**4.d0)*zeta !Turbulent Resistance
         ENDIF 
      ENDDO
 
      IF(max_reynolds.GT.2000)THEN !Checking for possible turbulence         
        write(*,*) "CHANGE CODE Max reynolds # (ne#=)=",
     '    max_reynolds,ne_reynolds
        !Use WRITE(OP_STRING...) and CALL WRITES for this
      ENDIF
      CALL EXITS('CALC_RESIS_FLOW')
      RETURN
 9999 CALL ERRORS('CALC_RESIS_FLOW',ERROR)
      CALL EXITS('CALC_RESIS_FLOW')
      RETURN 1
      END

