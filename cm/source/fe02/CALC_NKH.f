      SUBROUTINE CALC_NKH(NBH,NEELEM,NHP,NKH,NPNE,NPNODE,nr,NW,nx,
     '  ERROR,*)

C#### Subroutine: CALC_NKH
C###  Description:
C###    CALC_NKH calculates the array NKH.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NHP(NPM),
     '  NKH(NHM,NPM,NCM,0:NRM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NW(NEM,3),nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nc,ne,nh,nhx,nn,noelem,nonode,np

      CALL ENTERS('CALC_NKH',*9999)

C     Initialise NKH for nh's in this nx in current region
      DO nc=1,NCT(nr,nx)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
C MPN 3Mar96 DO nhx=1,NH_LOC(0,nx)
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            NKH(nh,np,nc,nr)=0
          ENDDO !nh
        ENDDO !np
      ENDDO !nc

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        DO nc=1,2
          DO nhx=1,NH_LOC(0,nx)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh,nc,ne)
            IF(nb.EQ.0.AND.nc.GT.1) THEN
              nb=NBH(nh,1,ne)
              NBH(nh,nc,ne)=nb
              WRITE(OP_STRING,'('' >>WARNING: NBH(nh='',I1,'',nc='','
     '          //'I1,'',ne='',I4,'') has not been set. Defaulted '','
     '          //''' to nc=1 value.'')') nh,nc,ne
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ENDIF
            DO nn=1,NNT(nb)
              np=NPNE(nn,nb,ne)
              IF(ITYP2(nr,nx).EQ.1.AND.NW(ne,1).EQ.3) THEN
                !linear elastic beam requires derivs for all nh
                NKH(nh,np,nc,nr)=2
              ELSE IF(ITYP2(nr,nx).EQ.1.AND.NW(ne,1).EQ.6) THEN
                !lin. elast plate may require derivs for all nh
                IF(NJT.EQ.3.AND.NBH(3,nc,ne).GT.0) THEN
                  !3D & bending incl
                  IF(NKT(nn,NBH(3,nc,ne)).EQ.4) THEN
                    NKH(nh,np,nc,nr)=4  !transv. defl.n is bicubic
                  ELSE
                    IF(NKT(nn,nb).GT.NKH(nh,np,nc,nr))
     '                NKH(nh,np,nc,nr)=NKT(nn,nb)
                  ENDIF
                ELSE
                  IF(NKT(nn,nb).GT.NKH(nh,np,nc,nr))
     '              NKH(nh,np,nc,nr)=NKT(nn,nb)
                ENDIF
              ELSE
                IF(NKT(nn,nb).GT.NKH(nh,np,nc,nr))
     '            NKH(nh,np,nc,nr)=NKT(nn,nb)
              ENDIF
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' nh='',I1,'' np='',I3,'
     '            //''' nc='',I1,'' NKH(nh,np,nc,nr)='',I2)')
     '            nh,np,nc,NKH(nh,np,nc,nr)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ENDDO !nn
          ENDDO !nh
        ENDDO !nc
      ENDDO !noelem (ne)

C     Find max NKH for zero'th nr location
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        DO nc=1,2
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            IF(NKH(nh,np,nc,nr).GT.NKH(nh,np,nc,0))
     '        NKH(nh,np,nc,0)=NKH(nh,np,nc,nr)
          ENDDO !nh
        ENDDO !nc
      ENDDO !nonode (np)

      CALL EXITS('CALC_NKH')
      RETURN
 9999 CALL ERRORS('CALC_NKH',ERROR)
      CALL EXITS('CALC_NKH')
      RETURN 1
      END


