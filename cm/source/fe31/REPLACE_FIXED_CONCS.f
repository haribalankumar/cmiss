      SUBROUTINE REPLACE_FIXED_CONCS(NBJ,NEELEM,NPNE,NVJE,nx,
     '  NXI,NYNP,XP,YP,FIX,ERROR,*)

C#### Subroutine: REPLACE_FIXED_CONCS
C###  Description:

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,ne0,nn,noelem,np,np1,npnn(2),nv,nv2,nvnn(2),ny,ny1,
     &  kount
      REAL*8 distance,length,radius,time,total_time,volume,xi
      REAL*8 LENGTH_1D,RADIUS_1D

      CALL ENTERS('REPLACE_FIXED_CONCS',*9999)
      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        np=NPNE(2,nb,ne) ! the first element for the node
        ny=NYNP(1,1,1,np)
        IF(FIX(ny,1))THEN
           total_time=0.d0
c evaluate the concentration by back tracking
           length=LENGTH_1D(NBJ,ne,NPNE,NVJE,XP)
           radius=RADIUS_1D(NBJ,ne,NPNE,NVJE,XP)
           volume=length*PI*radius**2
           time = volume/ (XP(1,1,nj_flow,np)*INLET_FLOW(nx))
           IF(time.GE.DT)THEN !within this element
              ny1=NYNP(1,1,1,NPNE(1,nb,ne))
              distance=DT*XP(1,1,nj_flow,np)*INLET_FLOW(nx)/
     &          (PI*radius**2)
              xi = 1.d0 - DT/time !i.e. we want to end of element
              YP(ny) = (1.d0-xi)*YP(ny1)+xi*YP(ny)
           ELSE
              ne0=ne
              total_time=total_time+time
              DO WHILE(total_time.LE.DT)
                 ne0=NXI(-1,1,ne0)
                 IF(ne0.NE.0)THEN
                    np=NPNE(2,nb,ne0)
                    length=LENGTH_1D(NBJ,ne0,NPNE,NVJE,XP)
                    radius=RADIUS_1D(NBJ,ne0,NPNE,NVJE,XP)
                    volume=length*PI*radius**2
                    time = volume/ (XP(1,1,nj_flow,np)*INLET_FLOW(nx))
                    total_time = total_time + time
                 ELSE
                    total_time=DT*2.d0
                 ENDIF
              ENDDO
              IF(ne0.EQ.0)THEN
                 YP(ny)=YP(1)
              ELSE
                 ! within element ne0
                 ny1=NYNP(1,1,1,NPNE(1,nb,ne))
                 xi = (total_time-DT)/time
c              IF(ne.EQ.15.OR.ne.EQ.21)THEN
c                 write(*,*) ne,'at xi=',xi,' in elem',ne0
c                 write(*,*) 'conc1=',YP(ny1),', conc2=',YP(ny)
c                 write(*,*) 'final conc=',(1.d0-xi)*YP(ny1)+xi*YP(ny)
c              ENDIF
                 YP(ny) = (1.d0-xi)*YP(ny1)+xi*YP(ny)


c              IF(ne.EQ.11.OR.ne.EQ.12.OR.ne.EQ.15.OR.ne.EQ.21)THEN
c                 write(*,*) ne,'at xi=',xi,' in elem',ne0,YP(ny)
c                 write(*,*) 'time=',total_time,', dt=',dt
c                 write(*,*) 'vol=',volume
c                 write(*,*) 'radius=',radius
c                 write(*,*) 'length=',length
c              ENDIF
c               write(*,*) ne,'at xi=',xi,' in elem',ne0,YP(ny)
              ENDIF
           ENDIF
        ELSE
c           write(*,*) 'not fixed at',np
        ENDIF
      ENDDO !noelem

      CALL EXITS('REPLACE_FIXED_CONCS')
      RETURN

 9999 CALL ERRORS('REPLACE_FIXED_CONCS',ERROR)
      CALL EXITS('REPLACE_FIXED_CONCS')
      RETURN 1
      END



