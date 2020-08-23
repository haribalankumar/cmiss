      SUBROUTINE CPXI(ie,IBT,IDO,INP,NPNE,nr,nx,CE,CP,CXI,XI,ERROR,*)

C#### Subroutine: CPXI
C###  Description:
C###    CPXI interpolates element parameters CE(il,ne) or nodal
C###    parameters CP(il,np) with linear basis functions (nb=NMB),
C###    for il=1,ILT(ie,nr,nx), to return Xi position params CXI(il).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ie,IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NPNE(NNM,NBFM),nr,nx
      REAL*8 CE(NMM),CP(NMM,NPM),CXI(NMM),XI(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER il,ILMAX,nb,nk,nn,np,ns
      REAL*8 CPE(64),PSI1

      CALL ENTERS('CPXI',*9999)

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
        CXI(il)=0.0d0 !initialise
        IF(ILP(il,1,nr,nx).EQ.4) THEN !defined by Gauss points
          CALL ASSERT(CALL_UPNODE_MAT,' >>> Gauss pt mat '
     '      //'param variation not implemented, '
     '      //'transfer to nodes',ERROR,*9999)
          NMB(il,ie,nx)=1
        ENDIF

        IF(ILP(il,1,nr,nx).EQ.1.OR.ILP(il,1,nr,nx).EQ.2) THEN
C         constant spatially or defined by elements
          CXI(il)=CE(il)

        ELSE IF((ILP(il,1,nr,nx).EQ.3)
     '    .OR.(ILP(il,1,nr,nx).EQ.4)) THEN !defined by nodes
          nb=NMB(il,ie,nx)
C         Put global node material params into element param array CPE
          ns=0
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb)
            DO nk=1,NKT(nn,nb)
              ns=ns+1
              IF(NKT(nn,nb).EQ.1) THEN
                CPE(ns)=CP(il,np)
              ELSE IF(NKT(nn,nb).EQ.2) THEN
C!!!NOTE: assumes nodes are numbered sequentially - should have nk in CP
                CPE(ns)=CP(il,2*(np-1)+nk)
              ELSE
                ERROR='Only 2 derivs allowed for material param interp'
                GO TO 9999
              ENDIF
            ENDDO !nk
          ENDDO !nn
C         Interpolate local element params CPE at Xi
          ns=0
          DO nn=1,NNT(nb)
            DO nk=1,NKT(nn,nb)
              ns=ns+1
              CXI(il)=CXI(il)+PSI1(IBT,IDO,INP,nb,1,nk,nn,XI)*CPE(ns)
            ENDDO !nk
          ENDDO !nn
        ENDIF
      ENDDO !il

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' CXI(1..): '',6D12.4,:/(11X,6D12.4))')
     '    (CXI(il),il=1,ILMAX)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF !dop

      CALL EXITS('CPXI')
      RETURN
 9999 CALL ERRORS('CPXI',ERROR)
      CALL EXITS('CPXI')
      RETURN 1
      END

      
