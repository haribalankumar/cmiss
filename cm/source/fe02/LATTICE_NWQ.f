      SUBROUTINE LATTICE_NWQ(NLATNE,NLATNQ,NLATPNQ,NQNLAT,NQS,NQXI,
     '  NWQ,nq,ERROR,*)

C#### Subroutine: LATTICE_NWQ
C###  Description:
C###    LATTICE_NWQ sets up the NWQ array with reference values
C###    that show what boundary type a grid point nq lies on. 
C###    These values are stored in NWQ(1,nq,na) and can
C###    be interpreted as:
C###    LOW Xi 1: 1,7,8,11,12,19,20,21,22
C###    HIGH Xi 1: 2,9,10,13,14,23,24,25,26
C###    LOW Xi 2: 3,7,9,15,16,19,20,23,24
C###    HIGH Xi 2: 4,8,10,17,18,21,22,25,26
C###    LOW Xi 3: 5,11,13,15,17,19,21,23,25
C###    HIGH Xi 3: 6,12,14,16,18,20,22,24,26
C###    Two examples to illustrate this are:
C###    Value of 1 is only low Xi 1 so is an Xi 1 face.
C###    Value of 24 is high Xi 1 and Xi 3 and low Xi 2
C###    so is a corner of those three faces.      
C###  See-Also: NWQ
      
      IMPLICIT NONE
      
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
!     Parameter List
      INTEGER NLATNE(NEQM+1),NLATNQ(NEQM*NQEM),NLATPNQ(NQM),
     '  NWQ(8,0:NQM,NAM),nq,NQNLAT(NEQM*NQEM),NQS(NEQM),
     '  NQXI(0:NIM,NQSCM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER EXT_TYPE,bound,i,i1,i2,i3,ii,j,j1,j2,j3,jj,k,k1,k2,k3,
     '  kk,loc,loc1,loc2,nbest,nbound,ne,ne1,ne2,ne3,nlat,nqsc,nqsc1
      LOGICAL NOTFOUND,EXT_I,EXT_J,EXT_K,FIRST
      
      CALL ENTERS('LATTICE_NWQ',*9999)

C     Work with the principal point first
      loc=NLATPNQ(nq)
      
C     Find the index of the NLATNE entry that is less than or
C     equal to nlat. ie, which element we are in.

      bound=0
      FIRST=.TRUE.
      
      nbest=0
      DO WHILE(NLATNQ(loc).NE.0)
        IF(FIRST) THEN
          FIRST=.FALSE.
        ELSE
          loc=NLATNQ(loc)
        ENDIF
        nbound=0
        CALL FIND_LATT_NEIJK(i,j,k,ne,loc,NLATNE,NQXI,NQS,nqsc,
     '    ERROR,*9999)
        IF(i.EQ.1.OR.i.EQ.NQXI(1,nqsc)) THEN
          nbound=1
        ENDIF
        IF(j.EQ.1.OR.j.EQ.NQXI(2,nqsc)) THEN
          nbound=nbound+1
        ENDIF
        IF(k.EQ.1.OR.k.EQ.NQXI(3,nqsc)) THEN
          nbound=nbound+1
        ENDIF
        IF(nbound.GT.nbest) THEN
          nbest=nbound
          bound=loc
          IF(nbound.EQ.3) THEN
            GOTO 100
          ENDIF
        ENDIF
      ENDDO

      IF(bound.EQ.0) bound=loc
      
      CALL FIND_LATT_NEIJK(i,j,k,ne,bound,NLATNE,NQXI,NQS,nqsc,
     '  ERROR,*9999)
 100  EXT_I=.FALSE.
      EXT_J=.FALSE.
      EXT_K=.FALSE.      
C     Check the boundaries of the lattice cube and assign proper
C     boundary support.
      IF((i.EQ.1.OR.i.EQ.NQXI(1,nqsc)).OR.
     '  (j.EQ.1.OR.j.EQ.NQXI(2,nqsc)).OR.
     '  (k.EQ.1.OR.k.EQ.NQXI(3,nqsc))) THEN !at edge of lattice cube
        IF((i.EQ.1.OR.i.EQ.NQXI(1,nqsc))
     '    .AND.(NQXI(1,nqsc).NE.1))THEN !check this face
          DO jj=1,NQXI(2,nqsc)
            DO kk=1,NQXI(3,nqsc)
              nlat=NLATNE(ne)+((kk-1)*NQXI(2,nqsc)+jj-1)*NQXI(1,nqsc)
     '          +i-1
              IF(NLATNQ(NLATPNQ(NQNLAT(nlat))).EQ.0) THEN !external face
                EXT_I=.TRUE.
                GOTO 200
              ENDIF
            ENDDO
          ENDDO
C         Also need to check if nq has grid points collapsed onto it and
C         check those points to see if they associate with boundary
C         faces.
          loc1=NLATPNQ(nq)
          DO WHILE(NLATNQ(loc1).NE.0) !have sub points
            loc1=NLATNQ(loc1)
            CALL FIND_LATT_NEIJK(i1,j1,k1,ne1,loc1,NLATNE,NQXI,NQS,
     '        nqsc1,ERROR,*9999)
            IF((i1.EQ.1.OR.i1.EQ.NQXI(1,nqsc1))
     '        .AND.(NQXI(1,nqsc1).NE.1)) THEN !check this face
              DO jj=1,NQXI(2,nqsc1)
                DO kk=1,NQXI(3,nqsc1)
                  nlat=NLATNE(ne1)+((kk-1)*NQXI(2,nqsc1)+jj-1)*NQXI(1,
     '              nqsc1)+i1-1
                  IF(NLATNQ(NLATPNQ(NQNLAT(nlat))).EQ.0) THEN !external face
                    EXT_I=.TRUE.
                    GOTO 200
                  ENDIF
                  loc2=NLATPNQ(NQNLAT(nlat))
                  CALL FIND_LATT_NEIJK(i3,j3,k3,ne3,loc2,NLATNE,
     '              NQXI,NQS,nqsc1,ERROR,*9999)
                  NOTFOUND=.TRUE.
                  DO WHILE(NLATNQ(loc2).NE.0.AND.NOTFOUND)
                    loc2=NLATNQ(loc2)
                    CALL FIND_LATT_NEIJK(i2,j2,k2,ne2,loc2,NLATNE,
     '                NQXI,NQS,nqsc1,ERROR,*9999)
                    IF(ne3.EQ.ne2) THEN
                      EXT_I=.TRUE.
                    ELSE
                      EXT_I=.FALSE.
                      NOTFOUND=.FALSE.
                    ENDIF
                  ENDDO
                  IF(EXT_I) GOTO 200
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
 200    IF((j.EQ.1.OR.j.EQ.NQXI(2,nqsc))
     '    .AND.(NQXI(2,nqsc).NE.1)) THEN !check this face         
          DO ii=1,NQXI(1,nqsc)
            DO kk=1,NQXI(3,nqsc)
              nlat=NLATNE(ne)+((kk-1)*NQXI(2,nqsc)+j-1)*NQXI(1,nqsc)
     '          +ii-1
              IF(NLATNQ(NLATPNQ(NQNLAT(nlat))).EQ.0) THEN !external face
                EXT_J=.TRUE.
                GOTO 400
              ENDIF
            ENDDO
          ENDDO
          loc1=NLATPNQ(nq)
          DO WHILE(NLATNQ(loc1).NE.0) !have sub points
            loc1=NLATNQ(loc1)
            CALL FIND_LATT_NEIJK(i1,j1,k1,ne1,loc1,NLATNE,NQXI,NQS,
     '        nqsc1,ERROR,*9999)
            IF((j1.EQ.1.OR.j1.EQ.NQXI(2,nqsc1))
     '        .AND.(NQXI(2,nqsc1).NE.1)) THEN !check this face
              DO ii=1,NQXI(1,nqsc1)
                DO kk=1,NQXI(3,nqsc1)
                  nlat=NLATNE(ne1)+((kk-1)*NQXI(2,nqsc1)+j1-1)*NQXI(1,
     '              nqsc1)+ii-1
                  IF(NLATNQ(NLATPNQ(NQNLAT(nlat))).EQ.0) THEN !external face
                    EXT_J=.TRUE.
                    GOTO 400
                  ENDIF
                  loc2=NLATPNQ(NQNLAT(nlat))
                  CALL FIND_LATT_NEIJK(i3,j3,k3,ne3,loc2,NLATNE,
     '              NQXI,NQS,nqsc1,ERROR,*9999)
                  NOTFOUND=.TRUE.                   
                  DO WHILE(NLATNQ(loc2).NE.0.AND.NOTFOUND)
                    loc2=NLATNQ(loc2)
                    CALL FIND_LATT_NEIJK(i2,j2,k2,ne2,loc2,NLATNE,
     '                NQXI,NQS,nqsc1,ERROR,*9999)
                    IF(ne3.EQ.ne2) THEN
                      EXT_J=.TRUE.
                    ELSE
                      EXT_J=.FALSE.
                      NOTFOUND=.FALSE.                       
                    ENDIF
                  ENDDO
                  IF(EXT_J) GOTO 400
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
 400    IF((k.EQ.1.OR.k.EQ.NQXI(3,nqsc))
     '    .AND.(NQXI(3,nqsc).NE.1)) THEN !check this face
          DO ii=1,NQXI(1,nqsc)
            DO jj=1,NQXI(2,nqsc)
              nlat=NLATNE(ne)+((k-1)*NQXI(2,nqsc)+jj-1)*NQXI(1,nqsc)
     '          +ii-1
              IF(NLATNQ(NLATPNQ(NQNLAT(nlat))).EQ.0) THEN !external face
                EXT_K=.TRUE.
                GOTO 600
              ENDIF
            ENDDO
          ENDDO
          loc1=NLATPNQ(nq)
          DO WHILE(NLATNQ(loc1).NE.0) !have sub points
            loc1=NLATNQ(loc1)
            CALL FIND_LATT_NEIJK(i1,j1,k1,ne1,loc1,NLATNE,NQXI,NQS,
     '        nqsc1,ERROR,*9999)
            IF((k1.EQ.1.OR.k1.EQ.NQXI(3,nqsc1))
     '        .AND.(NQXI(3,nqsc1).NE.1)) THEN !check this face
              DO ii=1,NQXI(1,nqsc1)
                DO jj=1,NQXI(2,nqsc1)
                  nlat=NLATNE(ne1)+((k1-1)*NQXI(2,nqsc1)+jj-1)*NQXI(1,
     '              nqsc1)+ii-1
                  IF(NLATNQ(NLATPNQ(NQNLAT(nlat))).EQ.0) THEN !external face
                    EXT_K=.TRUE.
                    GOTO 600
                  ENDIF
                  loc2=NLATPNQ(NQNLAT(nlat))
                  CALL FIND_LATT_NEIJK(i3,j3,k3,ne3,loc2,NLATNE,
     '              NQXI,NQS,nqsc1,ERROR,*9999)
                  NOTFOUND=.TRUE.                    
                  DO WHILE(NLATNQ(loc2).NE.0.AND.NOTFOUND)
                    loc2=NLATNQ(loc2)
                    CALL FIND_LATT_NEIJK(i2,j2,k2,ne2,loc2,NLATNE,
     '                NQXI,NQS,nqsc1,ERROR,*9999)
                    IF(ne3.EQ.ne2) THEN
                      EXT_K=.TRUE.
                    ELSE
                      EXT_K=.FALSE.
                      NOTFOUND=.FALSE.                       
                    ENDIF
                  ENDDO
                  IF(EXT_K) GOTO 600
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDIF
 600  IF(EXT_I.OR.EXT_J.OR.EXT_K) THEN !external grid point
        EXT_TYPE = 0
        IF(EXT_I) EXT_TYPE=EXT_TYPE+1
        IF(EXT_J) EXT_TYPE=EXT_TYPE+1
        IF(EXT_K) EXT_TYPE=EXT_TYPE+1
        IF(EXT_TYPE.EQ.1) THEN !face boundary
          IF(EXT_I) THEN
            IF(i.EQ.1) THEN
              NWQ(1,nq,1)=1
            ELSE
              NWQ(1,nq,1)=2
            ENDIF
          ELSEIF(EXT_J) THEN
            IF(j.EQ.1) THEN
              NWQ(1,nq,1)=3
            ELSE
              NWQ(1,nq,1)=4
            ENDIF
          ELSE !K_BOUND
            IF(k.EQ.1) THEN
              NWQ(1,nq,1)=5
            ELSE
              NWQ(1,nq,1)=6
            ENDIF
          ENDIF
        ELSEIF(EXT_TYPE.EQ.2) THEN !edge boundary
          IF(EXT_I.AND.EXT_J) THEN !IJ_EDGE
            IF(i.EQ.1) THEN
              IF(j.EQ.1) THEN
                NWQ(1,nq,1)=7
              ELSE !j==NQXI(2,nqsc)
                NWQ(1,nq,1)=8                  
              ENDIF
            ELSE
              IF(j.EQ.1) THEN
                NWQ(1,nq,1)=9                  
              ELSE !j==NQXI(2,nqsc)
                NWQ(1,nq,1)=10                  
              ENDIF                
            ENDIF
          ELSEIF(EXT_I.AND.EXT_K) THEN !IK_EDGE
            IF(i.EQ.1) THEN
              IF(k.EQ.1) THEN
                NWQ(1,nq,1)=11
              ELSE !k==NQXI(3,nqsc)
                NWQ(1,nq,1)=12                 
              ENDIF
            ELSE
              IF(k.EQ.1) THEN
                NWQ(1,nq,1)=13                 
              ELSE !k==NQXI(3,nqsc)
                NWQ(1,nq,1)=14                  
              ENDIF                
            ENDIF 
          ELSE !EXT_J AND EXT_K JK_EDGE
            IF(j.EQ.1) THEN
              IF(k.EQ.1) THEN
                NWQ(1,nq,1)=15
              ELSE !k==NQXI(3,nqsc)
                NWQ(1,nq,1)=16                 
              ENDIF
            ELSE
              IF(k.EQ.1) THEN
                NWQ(1,nq,1)=17                  
              ELSE !k==NQXI(3,nqsc)
                NWQ(1,nq,1)=18                  
              ENDIF                
            ENDIF
          ENDIF
        ELSE !EXT_TYPE==3 corner boundary
          IF(i.EQ.1) THEN
            IF(j.EQ.1) THEN
              IF(k.EQ.1) THEN
                NWQ(1,nq,1)=19
              ELSE !k==NQXI(3,nqsc)
                NWQ(1,nq,1)=20                  
              ENDIF
            ELSE !j=NQXI(2,nqsc)
              IF(k.EQ.1) THEN
                NWQ(1,nq,1)=21                  
              ELSE
                NWQ(1,nq,1)=22                  
              ENDIF
            ENDIF
          ELSE !i==NQXI(1,nqsc)
            IF(j.EQ.1) THEN
              IF(k.EQ.1) THEN
                NWQ(1,nq,1)=23                  
              ELSE
                NWQ(1,nq,1)=24                  
              ENDIF
            ELSE !j==NQXI(2,nqsc)
              IF(k.EQ.1) THEN
                NWQ(1,nq,1)=25                  
              ELSE !k==NQXI(3,nqsc)
                NWQ(1,nq,1)=26                  
              ENDIF
            ENDIF              
          ENDIF         
        ENDIF
      ELSE !internal grid point
        NWQ(1,nq,1)=0
        NWQ(2,nq,1)=0
      ENDIF


      CALL EXITS('LATTICE_NWQ')
      RETURN
 9999 CALL ERRORS('LATTICE_NWQ',ERROR)
      CALL EXITS('LATTICE_NWQ')
      RETURN 1
      END


