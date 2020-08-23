      SUBROUTINE MELGEB(ISC,ISR,LGEB,nb,nc,NDTOT,nh,NHST,NHE,NKHE,NKH,
     '  np,NPNE,NPNY,nr,nx,NYNE,NYNP,NZ_MAX,ERROR,*)

C#### Subroutine: MELGEB
C###  Description:
C###    MELGEB calculates the row numbers and column numbers ny's and
C###    their corresponding nz numbers for BEM element stiffness matrix.
C###    The routine calculates two variables, LGEB(-nkm:nhsm,nkm,2),
C###    and NHST. LGEB contains the global ny variable numbers
C###    (LGEB(x /= 0,x,1)), ny row number (LGEB(0,x,1)) and the nz
C###    numbers (LGEB(x,x,2)). NHST contains the number of nhs
C###    variables.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ISC(*),ISR(*),LGEB(-NKM:NHM*NSM,NKM,2),nb,nc,NDTOT,nh,NHE,
     '  NHST,NKHE(NKM,NNM),NKH(NHM,NPM),np,NPNE(NNM),
     '  NPNY(0:6,NYM,0:NRCM),nr,nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NZ_MAX
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER GETNYR,mk,NDMAX,nh2,nhs,nhx,nk,nkk,nn,npnn,nv,ny1,ny2,
     '  nz,type

      CALL ENTERS('MELGEB',*9999)

      nv=1 !temporary

      DO type=1,2
        DO nkk=1,NDTOT+1
          DO nhs=-NDTOT-1,NHE*NST(nb)
            LGEB(nhs,nkk,type)=0
          ENDDO !nhs
        ENDDO !nkk
      ENDDO !type

      ny1=NYNP(1,nv,nh,np,1,nc,nr)
      nhs=0
      DO nhx=1,NHE
        nh2=NH_LOC(nhx,nx)
        DO nn=1,NNT(nb)
          npnn=NPNE(nn)
          DO nk=1,NKT(nn,nb)
            mk=NKHE(nk,nn)
            nhs=nhs+1
            IF(mk.GT.0.AND.nk.LE.MAX(NKH(nh2,npnn)-KTYP93(nc,nr),1))
     '        THEN
              LGEB(nhs,1,1)=NYNP(mk,nv,nh2,npnn,0,nc,nr)
              ny2=NYNP(mk,nv,nh2,npnn,2,nc,nr)
              CALL SPARSE(ny1,ny2,NYT(1,nc,nx),nz,NZ_MAX,NZT(nc,nx),
     '          ISC,ISR,KTYP24,ERROR,*9999)
              LGEB(nhs,1,2)=nz
            ENDIF
          ENDDO !nk
        ENDDO !nn
      ENDDO !nhx
      NHST=nhs
      LGEB(0,1,1)=ny1
      IF(NDTOT.GT.0) THEN !Derivative
        DO nkk=1,NDTOT+1
          LGEB(-nkk,1,1)=NYNP(nkk,nv,nh,np,0,nc,nr)
          ny2=NYNP(nkk,nv,nh,np,2,nc,nr)
          IF(ny2.NE.0) THEN
            CALL SPARSE(ny1,ny2,NYT(1,nc,nx),nz,NZ_MAX,NZT(nc,nx),
     '        ISC,ISR,KTYP24,ERROR,*9999)
            LGEB(-nkk,1,2)=nz
          ELSE
            LGEB(-nkk,1,2)=0
          ENDIF
        ENDDO !nkk
        DO nkk=2,NDTOT+1
          ny1=NYNP(nkk,nv,nh,np,1,nc,nr)
          DO nhs=-NDTOT-1,NHST
            IF(nhs.EQ.0) THEN
              LGEB(0,nkk,1)=ny1
              LGEB(0,nkk,2)=0
            ELSE
              IF(LGEB(nhs,1,1).NE.0) THEN
                LGEB(nhs,nkk,1)=LGEB(nhs,1,1)
                ny2=GETNYR(nc,NPNY,nr,2,0,LGEB(nhs,1,1),NYNE,NYNP)
                IF(ny2.NE.0) THEN
                  CALL SPARSE(ny1,ny2,NYT(1,nc,nx),nz,NZ_MAX,NZT(nc,nx),
     '              ISC,ISR,KTYP24,ERROR,*9999)
                  LGEB(nhs,nkk,2)=nz
                ELSE
                  LGEB(nhs,nkk,2)=0
                ENDIF
              ELSE
                LGEB(nhs,nkk,1)=0
                LGEB(nhs,nkk,2)=0
              ENDIF
            ENDIF
          ENDDO !nhs
        ENDDO !nkk
      ENDIF

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MELGEB_1)
        WRITE(OP_STRING,'('' nb='',I2,'', nc='',I1,'', NDTOT='',I1,'
     '    //''', nh='',I2,'', NHST='',I2,'', np='',I5,'', nr='',I2)')
     '    nb,nc,NDTOT,nh,NHST,np,nr
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        IF(NDTOT.GT.0) THEN
          NDMAX=NDTOT+1
        ELSE
          NDMAX=1
        ENDIF
        DO type=1,2
          IF(type.EQ.1) THEN
            WRITE(OP_STRING,'(/'' type=1, ny values'')')
          ELSE
            WRITE(OP_STRING,'(/'' type=2, nz values'')')
          ENDIF
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO nkk=1,NDMAX
            WRITE(OP_STRING,'('' LGEB: '',15I10)') (LGEB(nhs,nkk,type),
     '        nhs=-NDMAX,NHST)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
CC$OMP END CRITICAL(MELGEB_1)
      ENDIF

      CALL EXITS('MELGEB')
      RETURN
 9999 CALL ERRORS('MELGEB',ERROR)
      CALL EXITS('MELGEB')
      RETURN 1
      END


