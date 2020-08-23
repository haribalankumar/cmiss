      SUBROUTINE MELGE(LGE,NBH,nc,ne,NHE,NHST,NPNE,nr,NVHE,nx,
     '  NYNE,NYNP,ERROR,*)

C#### Subroutine: MELGE
C###  Description:
C###    MELGE calculates the row numbers (LGE(*,1)) and column numbers
C###    (LGE(*,2)) in the global matrix nc for element variables nhs
C###    in region nr. It also returns the total number of element
C###    variables NHST(nrc).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER LGE(NHM*NSM,NRCM),NBH(NHM,NCM),nc,ne,NHE,NHST(2),
     '  NPNE(NNM,NBFM),nr,NVHE(NNM,NBFM,NHM),
     '  nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER na,nb,nh,nhs,nhx,nk,nn,np,nrc,nv

      CALL ENTERS('MELGE',*9999)

      DO nrc=1,2
        NHST(nrc)=0
        DO nhx=1,NHE
          nh=NH_LOC(nhx,nx)
C !!!     Use the LHS (nc=1) basis to determine the # of equations
          IF(nrc.EQ.1) THEN
            nb=NBH(nh,1)
          ELSE
            nb=NBH(nh,nc)
          ENDIF
          DO nn=1,NNT(nb)            !nodal variables
            np=NPNE(nn,nb)
            nv=NVHE(nn,nb,nh)
C CPB 12/9/95 Using NKT instead of NKH
C            IF(nrc.EQ.1) THEN
C              NK_TOT=MAX(NKH(nh,np,1)-KTYP93(1,nr),NKH(nh,np,2)-
C     '          KTYP93(2,nr),1)
C            ELSE
C              NK_TOT=MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
C            ENDIF
C            DO nk=1,NK_TOT
            DO nk=1,NKT(nn,nb)
              NHST(nrc)=NHST(nrc)+1
              LGE(NHST(nrc),nrc)=NYNP(nk,nv,nh,np,nrc,nc,nr)
            ENDDO !nk
          ENDDO !nn
          DO na=1,NAT(nb)            !auxillary variables
            NHST(nrc)=NHST(nrc)+1
            LGE(NHST(nrc),nrc)=NYNE(na,nh,nrc,nc,ne)
          ENDDO !na
        ENDDO !nh
      ENDDO !nrc

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MELGE_1)
        DO nrc=1,2
          WRITE(OP_STRING,'('' NHST('',I1,'')='',I4,'' LGE(nhs,'',I1,'
     '      //'''):'',/(20I4))') nrc,NHST(nrc),nrc,(LGE(nhs,nrc),nhs=1,
     '      NHST(nrc))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$OMP END CRITICAL(MELGE_1)
      ENDIF

      CALL EXITS('MELGE')
      RETURN
 9999 CALL ERRORS('MELGE',ERROR)
      CALL EXITS('MELGE')
      RETURN 1
      END


