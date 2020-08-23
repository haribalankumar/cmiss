      SUBROUTINE XGMG(IP,JAC,nb,nr,DXIX,GL,GU,RG,XG,ERROR,*)

C#### Subroutine: XGMG
C###  Description:
C###    XGMG evaluates the covariant (GL) & contravariant (GU) metric
C###    tensors wrt the Xi-coordinate system  and  the derivs  of
C###    the Xi-coords wrt the Xj-coords (DXIX) -if NIT=NJT only - at
C###    current Gauss pt.
C**** If IP=0 DXIX contains derivatives of Xi wrt X(ref)-coords.
C**** If IP=1 DXIX contains derivatives of Xi wrt Nu(fibre)-coords.
C**** If IP=-1 DXIX is not touched.
C**** The Jacobian RG for a length,area or volume integral is returned
C****   if JAC=1,2 or 3, respec.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IP,JAC,nb,nr
      REAL*8 DXIX(3,3),GL(3,3),GU(3,3),RG,XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,ni,NITB,njj,nu
      REAL*8 AA,D,DXXI(3,3),G,G1,G3,R,RC,RR,RRC,SLX,SMX

      CALL ENTERS('XGMG',*9999)

C     Calculate derivatives of X wrt Xi
      NITB=NIT(nb)
      DO ni=1,NITB
        nu=1+ni*(1+ni)/2
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          DXXI(njj,ni)=XG(njj,nu)
        ENDDO !njj
      ENDDO !ni

C     Initialise metric tensors
      DO mi=1,3
        DO ni=1,3
          GL(mi,ni)=0.0d0
          GU(mi,ni)=0.0d0
        ENDDO !ni
        GL(mi,mi)=1.0d0
        GU(mi,mi)=1.0d0
      ENDDO !mi

C     Calculate covariant metric tensor GL(i,j)
      IF(ITYP10(nr).EQ.2) THEN       !cyl polar
        R=XG(1,1)
        RR=R*R
      ELSE IF(ITYP10(nr).EQ.3) THEN  !sph polar
        R=XG(1,1)
        RR=R*R
        RC=R*DCOS(XG(3,1))
        RRC=RC*RC
      ELSE IF(ITYP10(nr).EQ.4) THEN  !prolate sph
        IF( DABS(XG(2,1)).LT.RDELTA) THEN
          RG=0.0d0
          WRITE(OP_STRING,'('' >>Warning: mu is zero in XGMG'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          GO TO 9998
        ENDIF
        AA=FOCUS*FOCUS
        SLX=DSINH(XG(1,1))
        SMX=DSIN(XG(2,1))
        G1=AA*(SLX*SLX+SMX*SMX)
        G3=AA* SLX*SLX*SMX*SMX
      ELSE IF(ITYP10(nr).EQ.5) THEN  !oblate sph
      ENDIF

      DO mi=1,NITB
        DO ni=1,NITB
          IF(ITYP10(nr).NE.4) GL(mi,ni)=DXXI(1,mi)*DXXI(1,ni)
          IF(NJ_LOC(NJL_GEOM,0,nr).GT.1) THEN
            IF(ITYP10(nr).EQ.1) THEN       !rect cart
              DO njj=2,NJ_LOC(NJL_GEOM,0,nr)
                GL(mi,ni)=GL(mi,ni)+DXXI(njj,mi)*DXXI(njj,ni)
              ENDDO !njj
            ELSE IF(ITYP10(nr).EQ.2) THEN  !cyl polar
              GL(mi,ni)=GL(mi,ni)+RR*DXXI(2,mi)*DXXI(2,ni)
              IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3)
     '          GL(mi,ni)=GL(mi,ni)+DXXI(3,mi)*DXXI(3,ni)
            ELSE IF(ITYP10(nr).EQ.3) THEN  !sph polar
              GL(mi,ni)=GL(mi,ni)+RRC*DXXI(2,mi)*DXXI(2,ni)
     '                             +RR *DXXI(3,mi)*DXXI(3,ni)
            ELSE IF(ITYP10(nr).EQ.4) THEN  !prolate sph
              GL(mi,ni)=G1*(DXXI(1,mi)*DXXI(1,ni)
     '                               +DXXI(2,mi)*DXXI(2,ni))
              IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) GL(mi,ni)=GL(mi,ni)+G3
     '                               *DXXI(3,mi)*DXXI(3,ni)
            ELSE IF(ITYP10(nr).EQ.5) THEN  !oblate sph
            ENDIF
          ENDIF
        ENDDO !ni
      ENDDO !mi

C new MPN 17-Apr-96: calc GU and DXIX with calls to INVERT and DXIDXM

C     Calculate contravariant metric tensor GU(i,j)
      CALL INVERT(NITB,GL,GU,G)
      IF(DABS(G).LT.RDELTA) THEN
        RG=0.0d0
        WRITE(OP_STRING,'('' >>Warning: zero G in XGMG. G='',D12.5,'
     '    //''' RDELTA='',D12.5)') G,RDELTA
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        GOTO 9998
      ENDIF

C     Calculate derivs DXIX(i,j) of Xi wrt X (IP=0) or Nu (IP=1)
      IF(IP.EQ.0) THEN !DXIX is based on reference coords, X
        IF(NITB.EQ.NJ_LOC(NJL_GEOM,0,nr)) CALL INVERT(NITB,DXXI,DXIX,D)
      ELSE IF(IP.GE.1) THEN !DXIX is based on material fibre coords, Nu
        CALL DXIDXM(NITB,nr,DXIX,RG,XG,'Fibre',ERROR,*9999)
      ENDIF

C old way of handling fibre/sheets
C old DATA M/1,2,3,1,2/
C
C      IF(NITB.EQ.1) THEN
C        G=GL(1,1)
C        IF(DABS(G).GT.RDELTA) THEN
C          GU(1,1)=1.0d0/G
C        ELSE
C          RG=0.0d0
C          WRITE(OP_STRING,'('' >>Warning: zero G in XGMG. G='',D12.5,'
C     '      //''' RDELTA='',D12.5)') G,RDELTA
C          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C          GOTO 9998
C        ENDIF
C        IF(IP.EQ.0) THEN !DXIX is based on reference coords
C          IF(DABS(DXXI(1,1)).GT.RDELTA) THEN
C            DXIX(1,1)=1.D0/DXXI(1,1)
C          ENDIF
C        ELSE IF(IP.GE.1) THEN !DXIX is based on material coords
C                              !.. ie arclength along 1D element
C          DXIX(1,1)=DSQRT(GL(1,1))
C        ENDIF
C      ELSE IF(NITB.EQ.2) THEN
C        G=GL(1,1)*GL(2,2)-GL(1,2)*GL(2,1)
C        IF(DABS(G).LT.RDELTA) THEN
C          RG=0.0d0
C          WRITE(OP_STRING,'('' >>Warning: zero G in XGMG. G='',D12.5,'
C     '      //''' RDELTA='',D12.5)') G,RDELTA
C          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C          GOTO 9998
C        ENDIF
C        GU(1,1)= GL(2,2)/G
C        GU(1,2)=-GL(1,2)/G
C        GU(2,1)=-GL(2,1)/G
C        GU(2,2)= GL(1,1)/G
C        GU(3,3)=1.0d0
C        IF(IP.EQ.0) THEN !DXIX is based on reference coords
C          IF(NITB.EQ.NJ_LOC(NJL_GEOM,0,nr)) THEN
C            D=DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1)
C            IF(DABS(D).GT.RDELTA) THEN
CC             calc DXIX(ni,nj) from DXXI(nj,ni)
C              DXIX(1,1)= DXXI(2,2)/D
C              DXIX(1,2)=-DXXI(1,2)/D
C              DXIX(2,1)=-DXXI(2,1)/D
C              DXIX(2,2)= DXXI(1,1)/D
C            ENDIF
C          ENDIF
C        ELSE IF(IP.GE.1) THEN !DXIX is based on material coords
C          IF(JTYP9.EQ.0) THEN
C            ETA1=0.0D0
C            ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))
C            C1=1.0d0
C            C2=DCOS(ETA2)
C          ELSE IF(JTYP9.GE.1) THEN
C            IF(JTYP12.EQ.1) THEN
C              ETA1=XG(NJ_LOC(NJL_FIBR,1,nr),1)
C              ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA1
C            ELSE IF(JTYP12.EQ.2) THEN
C              ETA2=-XG(NJ_LOC(NJL_FIBR,1,nr),1)
C              ETA1=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))+ETA2
C            ENDIF
C            C1=DCOS(ETA1)
C            C2=DCOS(ETA2)
C          ENDIF
C          RK1=DSQRT(GL(1,1))*C1
C          RK2=DSQRT(GL(2,2))*C2
C          RDSQ=GL(1,1)*GL(2,2)-GL(1,2)*GL(1,2)
C          RD=DSQRT(RDSQ)
C          RGU33=DSQRT(GU(3,3))
C          DXIX(1,1)=(GL(2,2)*RK1-GL(1,2)*RK2)/RDSQ
C          DXIX(2,1)=(GL(1,1)*RK2-GL(1,2)*RK1)/RDSQ
C          DXIX(3,1)= 0.0D0
C          DXIX(1,2)=-RK2/RD
C          DXIX(2,2)= RK1/RD
C          DXIX(3,2)= 0.0D0
C          DXIX(1,3)= 0.0D0
C          DXIX(2,3)= 0.0D0
C          DXIX(3,3)= RGU33
C        ENDIF
C      ELSE IF(NITB.EQ.3) THEN
C        G=DET(GL)         !Note DET should handle dble precision.
C        IF(DABS(G).LT.RDELTA) THEN
C          RG=0.0D0
C          WRITE(OP_STRING,'('' >>Warning: zero G in XGMG. G='',D12.5,'
C     '      //''' RDELTA='',D12.5)') G,RDELTA
C          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C          GOTO 9998
C        ENDIF
C        DO mi=1,3
C          DO ni=1,3
C            GU(mi,ni)=(GL(M(ni+1),M(mi+1))*GL(M(ni+2),M(mi+2))
C     '                -GL(M(ni+2),M(mi+1))*GL(M(ni+1),M(mi+2)))/G
C          ENDDO
C        ENDDO
C        IF(IP.EQ.0) THEN
C          IF(NITB.EQ.NJ_LOC(NJL_GEOM,0,nr)) THEN
C            D=DXXI(1,1)*
C     '        (DXXI(2,2)*DXXI(3,3)-DXXI(3,2)*DXXI(2,3))
C     '       +DXXI(1,2)*
C     '        (DXXI(2,3)*DXXI(3,1)-DXXI(3,3)*DXXI(2,1))
C     '       +DXXI(1,3)*
C     '        (DXXI(2,1)*DXXI(3,2)-DXXI(3,1)*DXXI(2,2))
C            DO ni=1,NITB
C              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
C                nj=NJ_LOC(NJL_GEOM,njj,nr)
C                DXIX(ni,nj)=(DXXI(M(nj+1),M(ni+1))*DXXI(M(nj+2),
C     '              M(ni+2))-DXXI(M(nj+2),M(ni+1))*DXXI(M(nj+1),
C     '              M(ni+2)))/D
C              ENDDO
C            ENDDO
C          ENDIF
C        ELSE IF(IP.GE.1) THEN
C          IF(JTYP9.EQ.0) THEN
C            ETA1=0.0D0
C            ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))
C            C1=1.0d0
C            C2=DCOS(ETA2)
C          ELSE IF(JTYP9.GE.1) THEN
C            IF(JTYP12.EQ.1) THEN
C              ETA1=XG(NJ_LOC(NJL_FIBR,1,nr),1)
C              ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA1
C            ELSE IF(JTYP12.EQ.2) THEN
C              ETA2=-XG(NJ_LOC(NJL_FIBR,1,nr),1)
C              ETA1=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))+ETA2
C            ENDIF
C            C1=DCOS(ETA1)
C            C2=DCOS(ETA2)
C          ENDIF
C          RK1=DSQRT(GL(1,1))*C1
C          RK2=DSQRT(GL(2,2))*C2
C          RDSQ=(GL(1,1)*GL(2,2)-GL(1,2)*GL(1,2))
C          RD=DSQRT(RDSQ)
C          RGU33=DSQRT(GU(3,3))
C          DXIX(1,1)=(GL(2,2)*RK1-GL(1,2)*RK2)/RDSQ
C          DXIX(2,1)=(GL(1,1)*RK2-GL(1,2)*RK1)/RDSQ
C          DXIX(3,1)= 0.0D0
C          DXIX(1,2)=-RK2/RD
C          DXIX(2,2)= RK1/RD
C          DXIX(3,2)= 0.0D0
C          DXIX(1,3)=RGU33*(GL(1,2)*GL(2,3)-GL(2,2)*GL(1,3))/RDSQ
C          DXIX(2,3)=RGU33*(GL(1,2)*GL(1,3)-GL(1,1)*GL(2,3))/RDSQ
C          DXIX(3,3)=RGU33
C        ENDIF
C      ENDIF

C     Calculate Jacobian RG
      IF(JAC.GT.0) THEN
        IF(JAC.EQ.1) RG=DSQRT(DABS(GL(1,1)))
        IF(JAC.EQ.2) RG=DSQRT(DABS(G*GU(3,3)))
        IF(JAC.EQ.3) RG=DSQRT(DABS(G))
      ENDIF

 9998 IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO mi=1,NITB
          WRITE(OP_STRING,'('' DXXI('',I1,'',ni): '',3(1X,D12.4))')
     '      mi,(DXXI(mi,ni),ni=1,NITB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !mi
        DO mi=1,3
          WRITE(OP_STRING,'('' GL('',I1,'',1..): '',3(1X,D12.4))')
     '      mi,(GL(mi,ni),ni=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !mi
        DO mi=1,3
          WRITE(OP_STRING,'('' GU('',I1,'',1..): '',3(1X,D12.4))')
     '      mi,(GU(mi,ni),ni=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !mi
        DO mi=1,NITB
          WRITE(OP_STRING,'('' DXIX('',I1,'',ni): '',3(1X,D12.4))')
     '      mi,(DXIX(mi,ni),ni=1,NITB)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !mi
        WRITE(OP_STRING,'('' RG: '',D12.4)') RG
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('XGMG')
      RETURN
 9999 CALL ERRORS('XGMG',ERROR)
      CALL EXITS('XGMG')
      RETURN 1
      END


