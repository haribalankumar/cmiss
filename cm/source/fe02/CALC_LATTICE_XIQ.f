      SUBROUTINE CALC_LATTICE_XIQ(NLATNE,NLATNQ,NLATPNQ,NQNLAT,NQS,
     &  NQSCNB,NQXI,XIQ,ERROR,*)

C#### Subroutine: CALC_LATTICE_XIQ
C###  Description:
C###    CALC_LATTICE_XIQ queries the lattice grid scheme
C###    mapping arrays to build the array XIQ.
C###  See-Also: XIQ    

      IMPLICIT NONE
      
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'       
!     Parameter List
      INTEGER NLATNE(NEQM+1),NLATNQ(NEQM*NQEM),NLATPNQ(NQM),
     &  NQNLAT(NEQM*NQEM),NQS(NEQM),NQSCNB(NQSCM),NQXI(0:NIM,NQSCM)
      REAL*8 XIQ(NIM,NQM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER EXTREME(3),bi,bj,bk,c,ei,ej,ek,i,j,k,nb,ne,ne1,ni,NITB,
     &  nlat,nq,nqsc,nqsc1,reptype,reptype2,si,sj,sk,U(3),L(3)

!     Functions
      INTEGER NLATUNIQUE
      
      CALL ENTERS('CALC_LATTICE_XIQ',*9999)

      DO nq=1,NQT
C       Loop over repeated grid points and from these points choose
C       the extreme point in the following order or precedence:
C       vertex/edge/face/internal
        nlat=NLATPNQ(nq)
C       Get element number that hosts nlat as well as nlat's
C       lattice based i,j and k coordinates.
        CALL FIND_LATT_NEIJK(i,j,k,ne,nlat,NLATNE,NQXI,NQS,nqsc,
     &    ERROR,*9999)
        reptype = 0
        IF(i.EQ.1.OR.i.EQ.NQXI(1,nqsc)) reptype = reptype+1
        IF(j.EQ.1.OR.j.EQ.NQXI(2,nqsc)) reptype = reptype+1
        IF(k.EQ.1.OR.k.EQ.NQXI(3,nqsc)) reptype = reptype+1
        EXTREME(1)=i
        EXTREME(2)=j
        EXTREME(3)=k
        DO WHILE(NLATNQ(nlat).NE.0)
          nlat=NLATNQ(nlat)
          CALL FIND_LATT_NEIJK(i,j,k,ne1,nlat,NLATNE,NQXI,NQS,nqsc1,
     &      ERROR,*9999)
          reptype2 = 0
          IF(i.EQ.1.OR.i.EQ.NQXI(1,nqsc1)) reptype2 = reptype2+1
          IF(j.EQ.1.OR.j.EQ.NQXI(2,nqsc1)) reptype2 = reptype2+1
          IF(k.EQ.1.OR.k.EQ.NQXI(3,nqsc1)) reptype2 = reptype2+1
          IF(reptype.LT.reptype2.AND.ne.EQ.ne1) THEN
            reptype=reptype2
            EXTREME(1)=i
            EXTREME(2)=j
            EXTREME(3)=k           
          ENDIF
        ENDDO !looking for extreme point
C       Search in the ijk directions from the extreme point location
C       to determine the number of unique grid points from the
C       extreme point.
        c=(NLATNE(ne)-1)+((EXTREME(3)-1)*NQXI(2,nqsc)+EXTREME(2)-1)*
     &    NQXI(1,nqsc)+EXTREME(1)
        bi=(NLATNE(ne)-1)+((EXTREME(3)-1)*NQXI(2,nqsc)+EXTREME(2)-1)*
     &    NQXI(1,nqsc)+1
        ei=(NLATNE(ne)-1)+((EXTREME(3)-1)*NQXI(2,nqsc)+EXTREME(2)-1)*
     &    NQXI(1,nqsc)+NQXI(1,nqsc)
        si=1
        bj=(NLATNE(ne)-1)+((EXTREME(3)-1)*NQXI(2,nqsc))*
     &    NQXI(1,nqsc)+EXTREME(1)
        ej=(NLATNE(ne)-1)+((EXTREME(3)-1)*NQXI(2,nqsc)+NQXI(2,nqsc)-1)*
     &    NQXI(1,nqsc)+EXTREME(1)
        sj=NQXI(1,nqsc)
        bk=(NLATNE(ne)-1)+(EXTREME(2)-1)*
     &    NQXI(1,nqsc)+EXTREME(1)
        ek=(NLATNE(ne)-1)+((NQXI(3,nqsc)-1)*NQXI(2,nqsc)+EXTREME(2)-1)*
     &    NQXI(1,nqsc)+EXTREME(1)
        sk=NQXI(1,nqsc)*NQXI(2,nqsc)
        L(1)=-1*NLATUNIQUE(-EXTREME(1),bi,si,c-si,NQNLAT)
        U(1)=NLATUNIQUE(EXTREME(1),c,si,ei-si,NQNLAT)
        L(2)=-1*NLATUNIQUE(-EXTREME(2),bj,sj,c-sj,NQNLAT)
        U(2)=NLATUNIQUE(EXTREME(2),c,sj,ej-sj,NQNLAT)
        L(3)=-1*NLATUNIQUE(-EXTREME(3),bk,sk,c-sk,NQNLAT)
        U(3)=NLATUNIQUE(EXTREME(3),c,sk,ek-sk,NQNLAT)
        nb=NQSCNB(nqsc)
        NITB=NIT(nb)
        DO ni=1,NITB
          IF(U(ni).EQ.L(ni)) THEN
            XIQ(ni,nq)=0.0d0
          ELSE
            XIQ(ni,nq)=1.0d0*(EXTREME(ni)-L(ni))/(U(ni)-L(ni))
          ENDIF
        ENDDO
      ENDDO
      
      CALL EXITS('CALC_LATTICE_XIQ')
      RETURN
 9999 CALL ERRORS('CALC_LATTICE_XIQ',ERROR)
      CALL EXITS('CALC_LATTICE_XIQ')
      RETURN 1      
      END


