      SUBROUTINE IPBAS6(NBJ,NBH,NEELEM,NHE,NHP,
     '  NPNODE,nr,NW,nx,ERROR,*)

C#### Subroutine: IPBAS6
C###  Description:
C###    IPBAS6 inputs basis functions for dependent variables for
C###    FE60 problems. This routine was copied from ipbas3.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM,0:NRM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NW(NEM,3),nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,INITIAL(1),ne,nh,nhx,nhx_MAX,noelem,NOQUES,
     '  nonode,np
      LOGICAL FILEIP
      CHARACTER CHAR1*2,CHAR2*2,CHAR3*8

      CALL ENTERS('IPBAS6',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        NW(ne,1)=1 ! Buuh !!? Dunno what this does
      ENDDO

      CALL ASSERT(NCM.GE.2,'>> Increase NCM >= 2',ERROR,*9999)

! Set NHE for the problem.
      nhx_MAX=0
      DO noelem=1,NEELEM(0,nr) !to set #dependent variables
        ne=NEELEM(noelem,nr)
        IF(ITYP5(nr,nx).EQ.2) THEN     !Time integration
          IF(ITYP2(nr,nx).EQ.5) THEN   !Navier-Stokes equations
            IF(ITYP3(nr,nx).EQ.3) THEN !General Navier Stokes
              NHE(ne)=NJ_LOC(NJL_GEOM,0,nr)+1 !Velocity + pressure
            ELSE IF(ITYP3(nr,nx).EQ.4) THEN !Stokes flow (no advection)
              NHE(ne)=NJ_LOC(NJL_GEOM,0,nr)+1 !Velocity + pressure
              ERROR='Stokes flow not implemented yet'
              GOTO 9999
            ELSE
              ERROR='Only Stokes and General-NS set up in ipbas6'
              GOTO 9999
            ENDIF
          ELSE
            ERROR='Only Navier-Stokes equations set up in ipbas6'
            GOTO 9999
          ENDIF
        ELSE
          ERROR='Only time integration equations set up in ipbas6'
          GOTO 9999
        ENDIF
        IF(NHE(ne).GT.nhx_MAX) nhx_MAX=NHE(ne)
      ENDDO !noelem

      CALL CALC_NH_LOC(nhx_MAX,nx,ERROR,*9999)

      DO nhx=1,NH_LOC(0,nx)
        nh=nh_loc(nhx,nx)
        WRITE(CHAR1,'(I2)') nhx
        WRITE(CHAR2,'(I2)') 1
        IF(nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN
          WRITE(CHAR3,'(A8)') 'pressure'
        ELSE
          WRITE(CHAR3,'(A8)') 'velocity'
        ENDIF
        FORMAT='($,'' The basis type number for dependent'//
     '    ' variable '//CHAR1(1:2)//' ('//CHAR3(1:8)//') is '//
     '    '['//CHAR2(1:2)//']: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,1,NEELEM(1,nr))
        INITIAL(1)=NBJ(1,NEELEM(1,nr))
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,INITIAL,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NBH(nh,1,ne)=IDATA(1) !variables
            NBH(nh,2,ne)=IDATA(1) !equations
          ENDDO !noelem
        ENDIF !iotype
      ENDDO !nhx

C     Set NHP for the problem.
      IF(IOTYPE.NE.3) THEN
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          NHP(np,nr)=NH_LOC(0,nx)
        ENDDO
      ENDIF !iotype.ne.3

      CALL EXITS('IPBAS6')
      RETURN
 9999 CALL ERRORS('IPBAS6',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPBAS6')
      RETURN 1
      END


