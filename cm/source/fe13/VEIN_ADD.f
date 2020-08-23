      SUBROUTINE VEIN_ADD(nb,NBJ,noelem_vein,nonode_vein,ne,NEELEM,NENP,
     &  ne_vein,ne_vein1,NKJ,NKJE,np1,np2,NPNE,NPNODE,np_vein,np_vein1,
     &  NRE,nr_vein,NVJE,NVJP,NXI,CE,SE,XP,HALF,ERROR,*)

C#### Subroutine: VEIN_ADD
C###  Description:
C###    VEIN_ADD adds new venous vessels in direction of host counterpart.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),noelem_vein,nonode_vein,ne,
     &  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),ne_vein,ne_vein1,
     &  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),np1,np2,NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),np_vein,np_vein1,NRE(NEM),nr_vein,
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM),SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL HALF
!     Local Variables
      INTEGER ne_v1,ne_v2,nj,np_v1,np_v2
      REAL*8 DOT_PROD,length,theta1,theta2,V0(NJT),V1(NJT),V2(NJT)
      
      CALL ENTERS('VEIN_ADD',*9999)
      
      np_vein=np_vein+1
      nonode_vein=nonode_vein+1
      CALL ASSERT(nonode_vein.LE.NP_R_M,'>>Increase NP_R_M',ERROR,*9999)
      CALL ASSERT(np_vein.LE.NPM,'>>Increase NPM',ERROR,*9999)              
      NPNODE(nonode_vein,nr_vein)=np_vein
      ne_v1=NXI(1,1,ne)
      ne_v2=NXI(1,2,ne)
      np_v1=NPNE(2,nb,ne_v1)
      np_v2=NPNE(2,nb,ne_v2)
      DO nj=1,NJT
        V0(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1)
        V1(nj)=XP(1,1,nj,np_v1)-XP(1,1,nj,np2)
        V2(nj)=XP(1,1,nj,np_v2)-XP(1,1,nj,np2)
      ENDDO
C... Now determine which vector to use???
      theta1=DCOS(DOT_PROD(V0,V1))
      theta2=DCOS(DOT_PROD(V0,V2))
      IF(.NOT.HALF) THEN
        IF(theta1.LT.theta2) THEN !use V1 vector and length1
          DO nj=1,NJT
            XP(1,1,nj,np_vein)=XP(1,1,nj,np_vein1)+V1(nj)
          ENDDO
        ELSE !use V2 vector and length2
          DO nj=1,NJT
            XP(1,1,nj,np_vein)=XP(1,1,nj,np_vein1)+V2(nj)
          ENDDO
        ENDIF
      ELSE
        IF(theta2.LT.theta1) THEN !use V1 vector and length1
          DO nj=1,NJT
            XP(1,1,nj,np_vein)=XP(1,1,nj,np_vein1)+V1(nj)
          ENDDO
        ELSE !use V2 vector and length2
          DO nj=1,NJT
            XP(1,1,nj,np_vein)=XP(1,1,nj,np_vein1)+V2(nj)
          ENDDO
        ENDIF
      ENDIF
      ne_vein=ne_vein+1
      noelem_vein=noelem_vein+1
      CALL ASSERT(noelem_vein.LE.NE_R_M,'>>Increase NE_R_M',ERROR,*9999)
      CALL ASSERT(ne_vein.LE.NEM,'>>Increase NEM',ERROR,*9999)
      NEELEM(noelem_vein,nr_vein)=ne_vein
      NPNE(1,nb,ne_vein)=np_vein1
      NPNE(2,nb,ne_vein)=np_vein
      CALL GN1DNEJ(nb,NBJ,ne_vein,NKJ,NKJE,NPNE(2,nb,ne_vein),
     &  NPNE(1,nb,ne_vein),nr_vein,NRE,NVJE,NVJP,SE,ERROR,*9999)
      NENP(np_vein1,0,nr_vein)=NENP(np_vein1,0,nr_vein)+1
      NENP(np_vein1,NENP(np_vein1,0,nr_vein),nr_vein)=
     '  ne_vein
      NXI(1,0,ne_vein1)=NXI(1,0,ne_vein1)+1
      length=0.d0
      DO nj=1,NJT 
        length=length+(XP(1,1,nj,NPNE(2,nb,ne_vein))-
     '    XP(1,1,nj,NPNE(1,nb,ne_vein)))**2.d0
      ENDDO
      length=DSQRT(length)
      CE(1,ne_vein)=length !stores segment length
      
      CALL EXITS('VEIN_ADD')
      RETURN
 9999 CALL ERRORS('VEIN_ADD',ERROR)
      CALL EXITS('VEIN_ADD')
      RETURN 1
      END  
