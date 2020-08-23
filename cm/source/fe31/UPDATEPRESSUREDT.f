      SUBROUTINE UPDATEPRESSUREDT(steps,NBJ,NEELEM,NPNE,NXI,CE,dPl,
     &  dt,time,ttime,T_interval,XAB,XP,DIAG_OP,INIT_ITER,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn' 
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),steps
      REAL*8 CE(NMM,NEM),dPl,dt,time,ttime,T_interval,XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL INIT_ITER
!     Local Variables
      INTEGER nb,ncount,ne,noelem,np,np2,NC1,NC2,NC3,NTOTAL
      REAL*8 dPla,dP_c,Palv,Q,Qinit,est,new_XAB
      LOGICAL DIAG_OP

      CALL ENTERS('UPDATEPRESSUREDT',*9999)

      ncount=0
      NC1=0
      NC2=0
      NC3=0
      NTOTAL=0
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        np=NPNE(1,nb,ne)       !want to store pressure at start node
        np2=NPNE(2,nb,ne)      !alveolar node
        Q=XP(2,1,nj_flow,np2)
        Qinit=XAB(2,ne)
        dPla=dPl*CE(nm_dPl,ne)
        IF(time.EQ.0.d0)THEN
           XAB(1,ne)=0.d0
        ELSE
          IF(NXI(1,0,ne).EQ.0)THEN !a terminal
            NTOTAL=NTOTAL+1

            IF(INIT_ITER)THEN
! linear estimate
               est=(XP(1,1,nj_pressure,np2)
     &              -XP(2,1,nj_pressure,np2))/dt 
! include previous two iterations
               new_XAB=0.75d0*XAB(3,ne)+0.25d0*
     &              (XAB(1,ne)+est)*0.5d0  

               XAB(3,ne)=XAB(1,ne)
               XAB(1,ne)=new_XAB

               XAB(5,ne)=XP(1,1,nj_pressure,np2)-XAB(4,ne)
               XAB(4,ne)=XP(1,1,nj_pressure,np2)

               dP_c = Q*dt/CE(nm_C,ne)

             IF(ne.EQ.8074.AND.DIAG_OP)THEN
               WRITE(OP_STRING,'(''element'',I6,'': Q ='',F8.3,'
     &           //''' mm^3/s, dpladt ='',F8.3,'
     &           //''' Pa/s, dpdt ='',F8.3,'' Pa/s, P_t = '',F8.3,'
     &           //''' Pa, dP_C ='',F8.3,'' Pa'')') ne,Q,
     &           dPla/dt,XAB(1,ne),XP(1,1,nj_pressure,np2),dP_c
               CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

            ENDIF
                  
         ENDIF                  !INIT_ITER
      ENDIF !terminal
      ENDIF !time
      ENDDO

      CALL EXITS('UPDATEPRESSUREDT')
      RETURN
 9999 CALL ERRORS('UPDATEPRESSUREDT',ERROR)
      CALL EXITS('UPDATEPRESSUREDT')
      RETURN 1
      END
