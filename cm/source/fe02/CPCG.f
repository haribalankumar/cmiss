      SUBROUTINE CPCG(ie,nb,NPNE,nr,nx,CE,CG,CGE,CP,PG,ERROR,*)

C#### Subroutine: CPCG
C###  Description:
C###    CPCG transfers element parameters CE(il,ne) or interpolates
C###    nodal parameters CP(il,np) with linear basis functions (nb=NMB),
C###    for il=1,ILT(ie,nr,nx), in element ne to return Gauss pt values
C###    CG(il,ng). ILP(il,ie,nr,nx) indicates whether parameter IL is
C###    constant(1),piecewise constant(2) or piecewise linear(3).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ie,nb,NPNE(NNM,NBFM),nr,nx
      REAL*8 CE(NMM),CG(NMM,NGM),CGE(NMM,NGM),
     '  CP(NMM,NPM),PG(NSM,NUM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER il,ILMAX,NBL,ng,nk,nm,nn,np,ns
      REAL*8 CPE(64),SUM

      CALL ENTERS('CPCG',*9999)
      DO nm=1,NMM
        DO ng=1,NGM
          CG(nm,ng)=0.0D0
        ENDDO
      ENDDO
      ILMAX=ILT(ie,nr,nx)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' nr='',I2,'' nx='',I2,'' ie='',I4,'
     '    //''' ILMAX='',I4)') nr,nx,ie,ILMAX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      DO il=1,ILMAX
        IF(ABS(ILP(il,ie,nr,nx)).EQ.1.OR.ABS(ILP(il,ie,nr,nx)).EQ.2)
     '    THEN
C       Constant spatially or piecewise constant
C       (Use ABS for KUNLE problem 2-7-91)
          DO ng=1,NGT(nb)
            CG(il,ng)=CE(il)
          ENDDO

        ELSE IF(ABS(ILP(il,ie,nr,nx)).EQ.3) THEN
C       Piecewise linear (nodal parameters)
          NBL=NMB(il,ie,nx)
          ns=0
          DO nn=1,NNT(NBL)
            np=NPNE(nn,NBL)
            DO nk=1,NKT(nn,NBL)
              ns=ns+1
              IF(NKT(nn,NBL).EQ.1) THEN
                CPE(ns)=CP(il,np)
              ELSE IF(NKT(nn,NBL).EQ.2) THEN
!!NOTE: assumes nodes are numbered sequentially - should have nk in CP
                CPE(ns)=CP(il,2*(np-1)+nk)
              ELSE
                ERROR='Only 2 derivs allowed for material param interp'
                GO TO 9999
              ENDIF
            ENDDO
          ENDDO !nn
          DO ng=1,NGT(NBL)
            SUM=0.0D0
            DO ns=1,NST(NBL)
              SUM=SUM+PG(ns,1,ng,NBL)*CPE(ns)
            ENDDO
            CG(il,ng)=SUM
          ENDDO !ng

        ELSE IF((ABS(ILP(il,ie,nr,nx)).EQ.4).OR.
     '      (ABS(ILP(il,ie,nr,nx))).EQ.5) THEN
C       Defined by Gauss points or Grid points
          DO ng=1,NGT(nb)
            CG(il,ng)=CGE(il,ng)
          ENDDO


        ELSE IF(ABS(ILP(il,ie,nr,nx)).EQ.5) THEN
C       Defined elsewhere
          DO ng=1,NGT(nb)
            CG(il,ng)=CE(il)
          ENDDO
        ENDIF
      ENDDO !il

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO ng=1,NGT(nb)
          WRITE(OP_STRING,'('' CG(1..,ng='',I3,''): '',6E11.3,'
     '      //'/(17X,6E11.3))') ng,(CG(il,ng),il=1,ILMAX)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF !dop

      CALL EXITS('CPCG')
      RETURN
 9999 CALL ERRORS('CPCG',ERROR)
      CALL EXITS('CPCG')
      RETURN 1
      END


