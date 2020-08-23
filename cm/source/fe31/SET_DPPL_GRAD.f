      SUBROUTINE SET_DPPL_GRAD(Gdirn,NBJ,NEELEM,NPNE,NVJE,NXI,CE,
     &  Pmax,Pmin,Ppl_step,XP,ERROR,*)

C#### Subroutine: SET_DPPL_GRAD
C###  Description:
C###    SET_DPPL_GRAD sets a linear gradient for the change in
C###    pleural pressure over a breath. Currently this is set
C###    to be equal everywhere. Needs setting up to use the option
C###    for distributing dPpl.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung_nej00.cmn' 
      INCLUDE 'pulm00.cmn' 

!     Parameter List
      INTEGER Gdirn,NBJ(NJM,NEM),NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM),Pmax,Pmin,Ppl_step,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)

      INTEGER nb,ne,noelem,np1,np2,nv1,nv2
      REAL*8 max_z,min_z,range_z,Xi

      CALL ENTERS('SET_DPPL_GRAD',*9999)

c      Pmax = Ppl_step !change in pressure at base of lungs
c      Pmin = Ppl_step !change in pressure at apex of lungs
      
      max_z=-1.d6
      min_z=1.d6
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).EQ.0)THEN
          np2=NPNE(2,1,ne)
          max_z=MAX(max_z,XP(1,1,Gdirn,np2))
          min_z=MIN(min_z,XP(1,1,Gdirn,np2))
        ENDIF
      ENDDO
      range_z=DABS(max_z-min_z)
      IF(DABS(range_z).LE.1d-5) range_z=1.d0

      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        np1=NPNE(1,1,ne)
        np2=NPNE(2,1,ne)
        Xi=(XP(1,1,Gdirn,np2)-min_z)/range_z
c        CE(nm_dPl,ne)=((1.d0-Xi)*(Pmax-Pmin)+Pmin)*Ppl_step
        CE(nm_dPl,ne)=((1.d0-Xi)*(Pmax-Pmin)+Pmin)
      ENDDO

      CALL EXITS('SET_DPPL_GRAD')
      RETURN
 9999 CALL ERRORS('SET_DPPL_GRAD',ERROR)
      CALL EXITS('SET_DPPL_GRAD')
      RETURN 1
      END

