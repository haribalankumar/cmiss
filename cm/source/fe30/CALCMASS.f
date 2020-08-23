      SUBROUTINE CALCMASS(NBJ,NEELEM,nh,NORD,NPNE,NVJE,NXI,NYNP,
     &  BBM,MASS,XAB,XP,YP,ERROR,*)

C#### Subroutine: CALCMASS
C###  Description:
C###    CALCMASS calculates the mass of inspired washin gas in a
C###    branching lung model.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      !Parameter list
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),nh,NORD(5,NE_R_M),
     &  NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 BBM(2,NEM),MASS,XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),YP(NYM)
      CHARACTER ERROR*(*)
      !local variables
      INTEGER nb,ne,ne1,ne2,ngen,nj,nlpm,noelem,np1,np2,nv1,nv2,ny1,ny2,
     &  GENELEM(17)
      REAL*8 AV_CE,cumulative,STOREMASS(16),STOREVOL(16),GenMass,length,
     &  radius,volume,GenVol,CumulVol

      DATA GENELEM/1,11,19,22,24,26,27,28,29,30,31,32,33,34,35,36,37/

      CALL ENTERS('CALCMASS',*9999)

c      DO ngen=1,16
c        STOREMASS(ngen)=0.d0
c        STOREVOL(ngen)=0.d0
c      ENDDO
      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        np1=NPNE(1,nb,ne)
        ny1=NYNP(1,1,nh,np1)
        np2=NPNE(2,nb,ne)
        ny2=NYNP(1,1,nh,np2)
        nv1=NVJE(1,nb,nj_radius,ne)
        nv2=NVJE(2,nb,nj_radius,ne)
        length=0.d0
        DO nj=1,NJT
          length=length+(XP(1,1,nj,np1)-XP(1,1,nj,np2))**2
        ENDDO !nj
        length=DSQRT(length)
        radius=0.5d0*(XP(1,nv1,nj_radius,np1)+XP(1,nv2,nj_radius,np2))
        volume=length*PI*radius**2
        
        AV_CE=(YP(ny1)+YP(ny2))*0.5d0 !average concentration

C... AJS 02/2011 Concentrations should never be negative. This can cause 
C... big mass errors. We need to add something to the solution methodology 
C... to identify where this is happening, and then force these nodes 
C... back to positive values. Will require additional iterations to find the 
C... steady state solution.
C... I have found this particularly problematic in small test geometries
C... with zero initial concentrations. If initial concentrations are >>0
C... then it is not a problem b/c if the solution steps past the initial value.
        WRITE(OP_STRING,'(''>> Concentration is negative for node '',
     &    I6, ''. Either modify initial concentrations, parameter '//
     &    'values, or modify solution methodology in code.'')') np1
        CALL ASSERT(YP(ny1).GE.0.d0,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''>> Concentration is negative for node '',
     &    I6, ''. Either modify initial concentrations, parameter '//
     &    'values, or modify solution methodology in code.'')') np2
        CALL ASSERT(YP(ny2).GE.0.d0,OP_STRING,ERROR,*9999)

        IF(ne.EQ.1)THEN
          XAB(6,ne)=volume*YP(ny2)
        ELSE
          XAB(6,ne)=volume*AV_CE
        ENDIF !ne.EQ.1
        XAB(7,ne)=volume
        IF(nord(1,ne).LT.5)THEN        
          STOREMASS(NORD(1,ne))=STOREMASS(NORD(1,ne))+XAB(6,ne)
          STOREVOL(NORD(1,ne))=STOREVOL(NORD(1,ne))+volume
        ENDIF
      ENDDO !noelem
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).EQ.0)THEN
          XAB(6,ne)=XAB(6,ne)+BBM(1,ne)*BBM(2,ne)
        ENDIF
      ENDDO
      DO noelem=NEELEM(0),1,-1
        ne1=NEELEM(noelem)
        ne2=NXI(-1,1,ne1)
        IF(ne2.NE.0)THEN !not top of trachea
          IF(NORD(5,ne1).EQ.1)THEN !start of a 'half' branch
            XAB(6,ne2)=XAB(6,ne2)+2.d0*XAB(6,ne1) !mass
            XAB(7,ne2)=XAB(7,ne2)+2.d0*XAB(7,ne1) !volume
          ELSE !within a tube branch
            XAB(6,ne2)=XAB(6,ne2)+XAB(6,ne1) !mass
            XAB(7,ne2)=XAB(7,ne2)+XAB(7,ne1) !volume
          ENDIF !ce
        ENDIF !ne2
      ENDDO !noelem
      MASS=XAB(6,NEELEM(1))

      IF(DOP)THEN
        IF(MASS.LT.1.d6)THEN !small values
          WRITE(OP_STRING,'('' Current mass='',E12.6)') MASS
        ELSE
          WRITE(OP_STRING,'('' Current mass='',E12.6)') MASS/1.d6
        ENDIF
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('CALCMASS')
      RETURN

 9999 CALL ERRORS('CALCMASS',ERROR)
      CALL EXITS('CALCMASS')
      RETURN 1
      END



