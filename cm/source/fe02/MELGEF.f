      SUBROUTINE MELGEF(LGE,NBH,ne,NHST,njj,NPNE,nr,NVHE,nx,NYNE,NYNP,
     '  ERROR,*)

C#### Subroutine: MELGEF
C###  Description:
C###    MELGEF calculates the row numbers (LGE(*,1)) and column numbers
C###    (LGE(*,2)) in the matrix for fitting for element variables nhs
C###    and fit variable njj in region nr.  It also returns the total
C###    number of element variables NHST(nrc).
C###  See-Also: MELGE

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER LGE(NHM*NSM,NRCM),NBH(NHM),ne,NHST(NRCM),njj,
     '  NPNE(NNM,NBFM),nr,NVHE(NNM,NBFM,NHM),
     '  nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER na,nb,nh,nhj,nhs,nhx,nk,nn,np,nrc,nv

      CALL ENTERS('MELGEF',*9999)

      DO nrc=1,2
        NHST(nrc)=0
        DO nhj=1,NUM_FIT(njj)
          nhx=NLH_FIT(nhj,3,njj)
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh)
          DO nn=1,NNT(nb) !nodal variables
C PM 14Aug2002: Section relevant to face fitting moved to a new routine
C               MELGEF_FACE.
            np=NPNE(nn,nb)
            nv=NVHE(nn,nb,nh)

C KAT 13Dec99: NKH is node dependent.  Here we are working with elements
C           so the appropriate array is the basis function dependent NKT
CC AJP 17/8/99 Need to use NKH and KTYP93
C            NK_TOT=MAX(NKH(nh,np,1)-KTYP93(1,nr),1)
C            DO nk=1,NK_TOT
            DO nk=1,NKT(nn,nb)
              NHST(nrc)=NHST(nrc)+1
              LGE(NHST(nrc),nrc)=NYNP(nk,nv,nh,np,nrc,1,nr)
            ENDDO !nk
          ENDDO !nn
          DO na=1,NAT(nb) !auxillary variables
            NHST(nrc)=NHST(nrc)+1
            LGE(NHST(nrc),nrc)=NYNE(na,nh,nrc,1,ne)
          ENDDO !na
        ENDDO !nh
      ENDDO !nrc

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nrc=1,2
          WRITE(OP_STRING,'('' NHST('',I1,'')='',I4,'' LGE(nhs,'',I1,'
     '      //'''):'',/(20I4))') nrc,NHST(nrc),nrc,(LGE(nhs,nrc),nhs=1,
     '      NHST(nrc))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('MELGEF')
      RETURN
 9999 CALL ERRORS('MELGEF',ERROR)
      CALL EXITS('MELGEF')
      RETURN 1
      END



