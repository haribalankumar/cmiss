      SUBROUTINE UPDATE_BBM(XAB,GAS,ERROR,*)

C#### Subroutine: UPDATE_BBM
C###  Description:
C###    UPDATE_BBM updates the BBA model parameters during breath
C###    simulations in a lung model.
C***  Created by Merryn Howatson Tawhai, April 1998

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'moti00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      REAL*8 XAB(NORM,NEM)
      CHARACTER ERROR*(*)
      LOGICAL GAS
!     Local Variables
      INTEGER Nhumi,nlpm,Ntemp

      CALL ENTERS('UPDATE_BBM',*9999)

      IF(GAS)THEN
        IF(NTB.GT.0)THEN
          DO nlpm=1,NTB !for each BBM
c            nh=NH_LOC(1,nx)
c            IF(INSP)THEN !start of inspiration
c              IF(LTYP_R(4).NE.3)THEN
c                ne=NXI(1,1,NINT(XAB(1,nlpm)))
c              ELSE
c                ne=NINT(XAB(1,nlpm))
c              ENDIF
c            ELSE IF((.NOT.EXPN).AND.(.NOT.INSP))THEN !start BH
c              CALL CALC_BBM(XAB(1,nlpm),reg_conc,ERROR,*9999)
c            ELSE IF(EXPN)THEN !start of expiration
c              IF(MIXING.EQ.1)THEN !perfect mixing in acini
                XAB(9,nlpm)=XAB(5,nlpm)
c              ELSE
c                ne=NINT(XAB(1,nlpm))
c                np=NPNE(2,nb,ne) !terminal node
c                ny=NYNP(1,1,nh,np)
c                XAB(12,nlpm)=YP(ny) !store conc at end of insp
c                XAB(3,nlpm)=XAB(2,nlpm)+XAB(7,nlpm)*Vt
c                XAB(4,nlpm)=XAB(3,nlpm)
c                CALL CALC_XAB(XAB(1,nlpm),reg_conc,ERROR,*9999)
c              ENDIF
c            ENDIF
          ENDDO !nlpm
        ENDIF
      ELSE !temperature and humidity
c        IF((.NOT.EXPN).AND.(.NOT.INSP))THEN !start BH
          Ntemp=0
          Nhumi=0
          DO nlpm=1,NTB
            IF(NTB.NE.5)THEN !BBMs are acini
              XAB(9,nlpm)=XAB(5,nlpm)
              XAB(10,nlpm)=XAB(6,nlpm)/XAB(3,nlpm)
              IF(XAB(9,nlpm).LE.0.98d0*CORE_TEMP) Ntemp=Ntemp+1
              IF(XAB(10,nlpm).LE.0.98d0) Nhumi=Nhumi+1
            ENDIF !NTB.NE.5
            XAB(9,nlpm)=CORE_TEMP
            XAB(10,nlpm)=100.d0
          ENDDO !nlpm
c        ELSE IF(EXPN)THEN !start of expiration
c          nh=NH_LOC(1,nx) !nh for temperature
c          DO nlpm=1,NTB
c            ne=NINT(BBM(1,nlpm))
c            np=NPNE(2,nb,ne) !terminal node
c            ny=NYNP(1,1,nh,np)
c            BBM(9,nlpm)=CORE_TEMP
c            BBM(10,nlpm)=CONC_EVAL(CORE_TEMP)
c            YP(ny)=BBM(9,nlpm)
c            YP(ny+1)=BBM(10,nlpm)
c          ENDDO !nlpm
c        ENDIF
      ENDIF !GAS

      CALL EXITS('UPDATE_BBM')
      RETURN

 9999 CALL ERRORS('UPDATE_BBM',ERROR)
      CALL EXITS('UPDATE_BBM')
      RETURN 1
      END


