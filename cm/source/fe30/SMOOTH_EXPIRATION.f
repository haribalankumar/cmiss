      SUBROUTINE SMOOTH_EXPIRATION(NBJ,NEELEM,nh,NORD,NPNE,NVJE,NXI,
     '  NYNP,XP,YP,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),nh,NORD(5,NE_R_M),
     &  NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG2,IEND2,ne,ne_next,ngen,ngen_next,NLIST(0:300),noelem,
     '  nogrno
      CHARACTER LABEL*30
      LOGICAL GROUP,SMOOTHED(NEM)
      
      CALL ENTERS('SMOOTH_EXPIRATION',*9999)
      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        SMOOTHED(ne)=.FALSE.
      ENDDO
C Set all ETUBE elements to have no smoothing      
      GROUP=.FALSE.
      DO nogrno=1,NTGRNO
        CALL CUPPER(LAGRNO(nogrno),LABEL)
        CALL STRING_TRIM(LABEL,IBEG2,IEND2)
        IF('ETUBE'.EQ.LABEL(IBEG2:IEND2)) GROUP=.TRUE.
      ENDDO
      IF(GROUP)THEN
        CDATA(1)='ELEMENTS'
        CALL PARSILG(NLIST,NEM,CDATA(1),'ETUBE',ERROR,*9999)
        DO noelem=1,NLIST(0)
          ne=NLIST(noelem)
          SMOOTHED(ne)=.TRUE.
        ENDDO
      ENDIF
C Set all STEM elements to have no smoothing      
      GROUP=.FALSE.
      DO nogrno=1,NTGRNO
        CALL CUPPER(LAGRNO(nogrno),LABEL)
        CALL STRING_TRIM(LABEL,IBEG2,IEND2)
        IF('STEM'.EQ.LABEL(IBEG2:IEND2)) GROUP=.TRUE.
      ENDDO
      IF(GROUP)THEN
        CDATA(1)='ELEMENTS'
        CALL PARSILG(NLIST,NEM,CDATA(1),'STEM',ERROR,*9999)
        DO noelem=1,NLIST(0)
          ne=NLIST(noelem)
          SMOOTHED(ne)=.TRUE.
        ENDDO
      ENDIF
      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).EQ.1)THEN !do only for discretised branch
          IF(.NOT.SMOOTHED(ne))THEN !hasn't been done as part of branch
            ngen=NORD(1,ne) !generation of ne
            ne_next=NXI(1,1,ne) !element # of adjacent element
            ngen_next=NORD(1,ne_next) !generation of adjacent element
            IF(ngen.EQ.ngen_next)THEN !same generation, therefore same branch
              NLIST(0)=1
              NLIST(1)=ne
              DO WHILE(ne_next.NE.0)
                NLIST(0)=NLIST(0)+1
                NLIST(NLIST(0))=ne_next
                SMOOTHED(ne_next)=.TRUE.
                IF(NXI(1,0,ne_next).EQ.1)THEN
                  ne_next=NXI(1,1,ne_next)
                  ngen_next=NORD(1,ne_next)
                  IF(ngen.NE.ngen_next)THEN
                    ne_next=0
                  ENDIF
                ELSE
                  ne_next=0
                ENDIF
              ENDDO !WHILE
              CALL SMOOTH_EXPIRATION_LS(NBJ,nh,NLIST,NPNE,NVJE,NYNP,XP,
     &          YP,ERROR,*9999)
            ENDIF !gen
          ENDIF !SMOOTHED
        ELSE !(NXI=0 or 2)
          SMOOTHED(ne)=.TRUE.
        ENDIF !NXI
      ENDDO !noelem(ne)
        
      CALL EXITS('SMOOTH_EXPIRATION')
      RETURN
 9999 CALL ERRORS('SMOOTH_EXPIRATION',ERROR)
      CALL EXITS('SMOOTH_EXPIRATION')
      RETURN 1
      END

      
