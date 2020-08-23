      SUBROUTINE CALCDIVU(NFVC,NODENVC,NPNODE,nr,NVCB,NVCNODE,nx,
     '  NYNP,GRR,XNFV,YP,ZNFVMSH,ERROR,*)

C#### Subroutine: CALCDIVU
C###  Description:
C###    <HTML>
C###    Forms the rhs of the pressure Poisson eqn
C###    <PRE>
C###
C###    Divergence(U') = rhs
C###
C###    Where U' = momentum marched non-solenoidal estimates of
C###               velocity
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NFVC(2,0:NFVCM,NVCM),NODENVC(NVCM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVCB(-1:3,NVCBM),NVCNODE(2,NP_R_M),
     '  nx,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 GRR(NOM),XNFV(-(NJM+1):NJM,NFVM),YP(NYM,NIYM),ZNFVMSH(NFVM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nvc,nonode,np,nfvl,cnonode,cnp,nfv,nh,nhx,cny,bnvc,ny
      REAL*8 DIV,VELNRM,DDOT


      CALL ENTERS('CALCDIVU',*9999)

      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)
        DIV=0.d0
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)
          cnp=NPNODE(cnonode,nr)
          nfv=NFVC(2,nfvl,nvc)

          IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
            VELNRM=0.d0

            DO nhx=1,nh_loc(0,nx)-1
              nh=nh_loc(nhx,nx)
              ny=NYNP(1,1,nh,np,0,1,nr)
              cny=NYNP(1,1,nh,cnp,0,1,nr)
              VELNRM=VELNRM+(YP(ny,1)+YP(cny,1))*XNFV(nhx,nfv)
            ENDDO

            DIV=DIV+0.5d0*XNFV(FAREA,nfv)*VELNRM

          ELSE !Boundary
            bnvc=NVCNODE(MAP,cnonode)

            IF(NVCB(BCTYPE,bnvc).EQ.INLET) THEN

              VELNRM=0.d0
              DO nhx=1,nh_loc(0,nx)-1
                nh=nh_loc(nhx,nx)
                ny=NYNP(1,1,nh,np,0,1,nr)
                cny=NYNP(1,1,nh,cnp,0,1,nr)
                VELNRM=VELNRM+(YP(ny,1)+YP(cny,3))*XNFV(nhx,nfv)
              ENDDO

              DIV=DIV+0.5d0*XNFV(FAREA,nfv)*VELNRM

            ELSEIF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN

              VELNRM=0.d0
              DO nhx=1,nh_loc(0,nx)-1
                nh=nh_loc(nhx,nx)
                ny=NYNP(1,1,nh,np,0,1,nr)
                VELNRM=VELNRM+YP(ny,1)*XNFV(nhx,nfv)
              ENDDO

              DIV=DIV+XNFV(FAREA,nfv)*VELNRM

            ENDIF

            IF(.NOT.MESHFIXD) THEN
              DIV=DIV+ZNFVMSH(nfv)
            ENDIF
          ENDIF
        ENDDO
        GRR(nvc)=DIV/DT
      ENDDO

C     ..Calculate whole field divergence magnitude
      DIVMAG=DT*DDOT(NVCT,GRR,1,GRR,1)
      DIVMAG=DSQRT(DIVMAG)

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(CALCDIVU_1)
        WRITE(OP_STRING,'(/$,'' RHS listing of PPE:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ###################'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nvc=1,NVCT
          WRITE(OP_STRING,'('' Voronoi cell: '',I7,'' RHS:'',D16.8)')
     '      nvc,GRR(nvc)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nvc
CC$OMP END CRITICAL(CALCDIVU_1)
      ENDIF !dop

      CALL EXITS('CALCDIVU')
      RETURN
 9999 CALL ERRORS('CALCDIVU',ERROR)
      CALL EXITS('CALCDIVU')
      RETURN 1
      END


