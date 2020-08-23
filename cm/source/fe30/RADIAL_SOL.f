      SUBROUTINE RADIAL_SOL(NBJ,NENP,NPNE,NPNODE,nr,nx,NXI,NYNP,CE,CP,
     &  ERR,XP,YP,RET_ERROR,*)

C#### Subroutine: RADIAL_SOL
C###  Description:
C###    RADIAL_SOL

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'tol00.cmn'

      !Parameter list
      INTEGER NBJ(NJM,NEM),NENP(NPM,0:NEPM),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M),nr,nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  NYNP(NKM,NVM,NHM,NPM)
      REAL*8 CE(NMM,NEM),CP(NMM,NPM),ERR,XP(NKM,NVM,NJM,NPM),
     &  YP(NYM,NIYM)
      CHARACTER RET_ERROR*(*)
      !Local variables
      INTEGER nb,ne,ne0,nh,nj,nonode,np,np1,np2,ny1,ny2
      REAL*8 beta,Clumen,Cmax,flux,LENGTH,Lpvc,tissue_depth,Tlumen,
     &  Tlumen0,Twall,Vevap
      REAL*8 CONC_EVAL
      LOGICAL ET_TUBE
      CHARACTER ERROR*(ERRSTRLEN)
      LOGICAL ERROR_FLAG

      CALL ENTERS('RADIAL_SOL',*9999)

      ERR=0.d0
      ERROR_FLAG=.FALSE.
      
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        ne=NENP(np,1)
        nb=NBJ(1,ne)
        np1=NPNE(1,nb,ne)
        np2=NPNE(2,nb,ne)
        length=0.d0
        DO nj=1,NJT
          length=length+(XP(1,1,nj,np1)-XP(1,1,nj,np2))**2
        ENDDO !nj
        length=DSQRT(length)

c        IF(.NOT.ERROR_FLAG)THEN
c          ne=NENP(np,1)
          ny1=NYNP(1,1,NH_LOC(1,nx),np)
          ny2=NYNP(1,1,NH_LOC(2,nx),np)
          IF(PULMAT(5).NE.0.d0)THEN !heat transfer across wall
            ET_TUBE =.FALSE.
c            IF(ne.LT.BEGIN_ELEM(5))THEN
c              ET_TUBE=.TRUE.
c              Lpvc=2.d0
c            ELSE
              Lpvc=0.d0
c            ENDIF
c            IF(ne.LT.BEGIN_ELEM(4)) THEN !In humi circuit
c              Twall=AMBIENT
c            ELSE
c              ne0=NXI(-1,1,ne) !parent
c              LENGTH=0.d0
c              DO WHILE(ne0.NE.0)
c                LENGTH=LENGTH+CE(3,ne0) !add parent length
c                ne0=NXI(-1,1,ne0)
c              ENDDO !WHILE (ne0.NE.0)
c              IF(ne.LE.BEGIN_ELEM(5))THEN
c                Twall=MOUTH_TEMP
c     '            +(CORE_TEMP-MOUTH_TEMP)/350.d0*LENGTH
c              ELSE
c                IF(ne_trachea.EQ.1)THEN !assuming no et tube
c                  Twall=MOUTH_TEMP
c     '              +(CORE_TEMP-MOUTH_TEMP)/100.d0*LENGTH
c                ELSE
              Twall=CORE_TEMP
c                ENDIF
c              ENDIF
c              IF(Twall.GT.CORE_TEMP) Twall=CORE_TEMP
c            ENDIF !ne.LT.BEGIN_ELEM(2)
c            Clumen=YP(ny2,1) ! this gives convergence, don't change!
            Clumen=YP(ny2,3) !actually doesn't seem to make a diff
            Tlumen=YP(ny1,3) !previous solution
            Tlumen0=YP(ny1,1) !current solution
            tissue_depth=2.d0*MAX(0.25d0*XP(1,1,nj_radius,np),2.d0)
            beta = XP(1,1,nj_coeff,np)
            CALL CALC_WALLTEMP(np,nr,nx,beta,Clumen,PULMAT(6),flux,
     &        PULMAT(1),length,Lpvc,tissue_depth,Tlumen,Tlumen0,Twall,
     &        CP(1,np),Vevap,XP(1,1,1,np),ET_TUBE,ERROR,*9999)
c            Vevap_sum=Vevap_sum+Vevap
            IF(YP(ny1,2).NE.0.d0)THEN
              ERR=ERR+DABS(((YP(ny1,2)-CP(10,np))/YP(ny1,2))**2.d0)
            ENDIF
            IF(YP(ny2,2).NE.0.d0)THEN   
              ERR=ERR+DABS(((YP(ny2,2)-CP(9,np))/YP(ny2,2))**2.d0)
            ENDIF       
            YP(ny1,2)=CP(10,np) !wall temperature
            YP(ny2,2)=CP(9,np) !wall concentration
          ENDIF
c          GO TO 102
C         This statement is designed to be skipped if no error
C         occurs. However if a error occurs within a subroutine
C         the alternate return points to line 100 to set the flag
c 100      CONTINUE
c          ERROR_FLAG=.TRUE.
c          WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
c          CALL WRITES(IOER,OP_STRING,ERROR,*101)
c          WRITE(OP_STRING,'(/'' >>An error occurred - '
c     '      //'results may be unreliable!'')')
c          CALL WRITES(IODI,OP_STRING,ERROR,*101)
c 101      CONTINUE
c 102      CONTINUE
c        ENDIF !.NOT.ERROR_FLAG
      ENDDO !nonode


      nh=NH_LOC(1,nx)
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        ny1=NYNP(1,1,nh,np) !corresp to temperature
        ny2=NYNP(1,1,nh+1,np) !corresp to water vapour
        IF(YP(ny1,1).GT.TMAX) YP(ny1,1)=TMAX
        IF(YP(ny1,1).LT.TMIN) YP(ny1,1)=TMIN
        Cmax=CONC_EVAL(YP(ny1,1))
        IF(YP(ny2,1).GT.Cmax) YP(ny2,1)=Cmax
        YP(ny2,8)=YP(ny2,1)/Cmax*100.d0 !relative humidity
      ENDDO !nonode

      CALL EXITS('RADIAL_SOL')
      RETURN
 9999 CALL ERRORS('RADIAL_SOL',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('RADIAL_SOL')
      RETURN 1
      END



