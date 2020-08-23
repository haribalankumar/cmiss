      SUBROUTINE OPCOUP(NEELEM,NP_INTERFACE,NPNODE,ERROR,*)

C#### Subroutine: OPCOUP
C###  Description:
C###    OPCOUP outputs coupling data.

      IMPLICIT NONE
      INCLUDE 'aero00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NP_INTERFACE(0:NPM,0:3),
     '  NPNODE(0:NP_R_M,0:NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ne,noelem,NO_INTERFACE,nonode,NO_TIP,np,nr,NR_INTERFACE

      CALL ENTERS('OPCOUP',*9999)

      IF(KTYP90.EQ.1) THEN      !Coupled saturated-unsaturated
        WRITE(OP_STRING,'('' Coupling is coupled saturated-'','
     '    //'''unsaturated'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ELSE IF(KTYP90.EQ.2) THEN      !Coupled Laplace
        WRITE(OP_STRING,'('' Coupling is coupled Laplace'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(KTYP91.EQ.1)THEN
          WRITE(OP_STRING,'('' Full equation generation at'//
     '      ' interface'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSEIF(KTYP91.EQ.2)THEN
          WRITE(OP_STRING,'('' Derivative equations generated'//
     '      ' at first interface region only'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        WRITE(OP_STRING,'('' Interface nodes :'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nr=1,NRT
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            IF(NP_INTERFACE(np,0).GT.1)THEN !np is an interface node
              WRITE(OP_STRING,'(''  Node '',I3,'
     '          //''' shared by regions :'',20I5)')
     '          np,(NP_INTERFACE(np,NR_INTERFACE),NR_INTERFACE=1,
     '          NP_INTERFACE(np,0))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !nonode (np)
        ENDDO !nr

      ELSE IF(KTYP90.EQ.3) THEN !Coupled tree growth problem
        WRITE(OP_STRING,'('' Coupling is coupled tree growth prob.'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Growing tip nodes'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NEELEM(0,1)
          ne=NEELEM(noelem,1)
          WRITE(OP_STRING,'('' Node'',I5,'': '',20I6)')
     '      (NP_GROWING_TIP(NO_TIP,ne),NO_TIP=1,NP_GROWING_TIP(0,ne))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''   Flows: '',10E12.3)')
     '      (FLOW_GROWING_TIP(NO_TIP,ne),NO_TIP=1,NP_GROWING_TIP(0,ne))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

      ELSE IF(KTYP90.EQ.4) THEN !Fluid interface stability problem
        WRITE(OP_STRING,'('' Coupling is fluid interface stability'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Fluid interface nodes'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nr=1,3
          WRITE(OP_STRING,'('' Nodes for region '',I1,'':'',20I5)')
     '      nr,(NP_INTERFACE(NO_INTERFACE,nr),NO_INTERFACE=1,
     '      NP_INTERFACE(0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

      ELSE IF(KTYP90.EQ.5) THEN !Coupled aerofoil flow & stress
        WRITE(OP_STRING,'('' Coupling is coupled aerofiol flow '','
     '    //'''& stress'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Node triplets on aerofoil:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''   Upper surface flow field nodes:'','
     '    //'20I5)')
     '    (NP_INTERFACE(NO_INTERFACE,1),
     '    NO_INTERFACE=1,NP_INTERFACE(0,1))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''   Lower surface flow field nodes:'','
     '    //'20I5)')
     '    (NP_INTERFACE(NO_INTERFACE,2),
     '    NO_INTERFACE=1,NP_INTERFACE(0,2))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''   Aerofoil stress field nodes   :'','
     '    //'20I5)')
     '    (NP_INTERFACE(NO_INTERFACE,3),
     '    NO_INTERFACE=1,NP_INTERFACE(0,3))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Basis# for aerofoil pressure is '',I2)')
     '    NB_AERO_PRESS
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('OPCOUP')
      RETURN
 9999 CALL ERRORS('OPCOUP',ERROR)
      CALL EXITS('OPCOUP')
      RETURN 1
      END


