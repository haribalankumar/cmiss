      SUBROUTINE DEXI_CLOSEST(IBT,IDO,INP,LD,NBJ,
     '  NBH,ND0,ND1,NELIST,NHE,nj1,
     '  NKHE,NKJE,NPF,NPNE,NRE,NVHE,NVJE,NW,nx,NXI,
     '  START_ELEMENT,CURVCORRECT,SE,SQ,SQMAX,XA,XE,XI,
     '  XI_1,XI_2,XI_3,XID,XP,ZA,ZD,ZP,DEFORM,
     '  NEW,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY,
     '  ERROR,*)

C#### Subroutine: DEXI_CLOSEST
C###  Description:
C###    DEXI_CLOSEST calculates the closest point on any of the elements
C###    in NELIST for each data point

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),NBJ(NJM,NEM),NBH(NHM,NCM,NEM),ND0,ND1,
     '  NELIST(0:NEM),NHE(NEM,NXM),nj1,NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),START_ELEMENT
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  SQ(NDM),SQMAX,XA(NAM,NJM,NEM),XE(NSM,NJM),XI(3),XI_1,XI_2,XI_3,
     '  XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM),
     '  ZD(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DEFORM,NEW,SET_XI_1,SET_XI_2,SET_XI_3,SPECIFY
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IT,ITMAX,nb,nd,ne,neadj,nelast,neold,ni,NITB,nj,no_ne
      REAL*8 SQND,TEMP
      LOGICAL FOUND,ERROR_FLAG
!     Functions
      REAL*8 PXI

      PARAMETER(ITMAX=20)

      CALL ENTERS('DEXI_CLOSEST',*9999)

      ERROR_FLAG=.FALSE.
      

      IF(NEW) THEN
C       Make a guess at the best element
        DO ni=1,3 !compare the distances to the centres of each element
          XI(ni)=0.5D0
        ENDDO
        IF(SET_XI_3) THEN
          XI(3)=XI_3
        ENDIF
        DO nd=ND0,ND1
          LD(nd)=0
          SQ(nd)=0.0D0
        ENDDO ! nd

C$OMP   PARALLEL DO
C$OMP&  PRIVATE(no_ne,ne,NITB,XE,nd,SQND,nj,nb,TEMP)
        DO no_ne=1,NELIST(0)
          IF(.NOT.ERROR_FLAG) THEN
            ne=NELIST(no_ne)
            NITB=NIT(NBJ(NJ1,ne))
            IF(DEFORM)THEN
c           Note that deformed coords --> XE
              CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     '          NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),NRE(ne),
     '          NVHE(1,1,1,ne),NW(ne,1),nx,
     '          CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '          XE,ZP,ERROR,*100)
            ELSE
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),
     '          NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*100)
            ENDIF
            DO nd=ND0,ND1
              SQND=0.0d0
              DO nj=1,NJT
                nb=NBJ(nj,ne)
                TEMP=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '            XE(1,nj))-ZD(nj,nd)
                SQND=SQND+TEMP*TEMP
              ENDDO
              IF(LD(nd).EQ.0.OR.SQND.LT.SQ(nd)) THEN
C$OMP CRITICAL(DEXI_CLOSEST1)
                DO ni=1,NITB !copy new xi values
                  XID(ni,nd)=XI(ni)
                ENDDO
                LD(nd)=ne
                SQ(nd)=SQND
C$OMP END CRITICAL(DEXI_CLOSEST1)
              ENDIF
            ENDDO ! nd
            GO TO 102
C             This statement is designed to be skipped if no error
C             occurs. However if a error occurs within a subroutine
C             the alternate return points to line 100 to set the flag
 100        CONTINUE
C$OMP CRITICAL(DEXI_CLOSEST2)
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
            CALL WRITES(IOER,OP_STRING,ERROR,*101)
            WRITE(OP_STRING,'(/'' >>An error occurred - '
     '        //'results may be unreliable!'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101        CONTINUE
C$OMP END CRITICAL(DEXI_CLOSEST2)
 102        CONTINUE

          ENDIF !error_flag
        ENDDO !no_ne
C$OMP   END PARALLEL DO

      ENDIF !NEW


C$OMP   PARALLEL DO
C$OMP&  PRIVATE(nd,nelast,ne,NITB,XI,FOUND,IT,XE,SQND,neold,ni,neadj)
      DO nd=ND0,ND1 !new AAY 23 March 95 limits set with ADD keyword
        IF(.NOT.ERROR_FLAG) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Data point '',I5)') nd
            CALL WRITES(IODI,OP_STRING,ERROR,*200)
          ENDIF

C       Start with the guess element
          nelast=0

C LKC 5-JUN-2002 move inside specify loop
C            ne=START_ELEMENT
          IF(SPECIFY) THEN !Set any xi values
            ne=START_ELEMENT
C LKC 5-JUN-2002 Also initialise NITB
            NITB=NIT(NBJ(NJ1,ne))
            IF(SET_XI_1) THEN
              XI(1)=XI_1
            ENDIF
            IF(SET_XI_2) THEN
              XI(2)=XI_2
            ENDIF
            IF(SET_XI_3) THEN
              XI(3)=XI_3
            ENDIF
          ELSE !start from previous values
            ne=LD(nd)
            NITB=NIT(NBJ(NJ1,ne))
            DO ni=1,NITB
              XI(ni)=XID(ni,nd)
            ENDDO
          ENDIF
          FOUND=.FALSE.

          IF(ne.NE.0) THEN
            IT=0
            DO WHILE(.NOT.FOUND.AND.IT.LT.ITMAX)
              IT=IT+1
C           Search the elements around the guess element
              IF(DOP) THEN
                WRITE(OP_STRING,'(/'' Element '',I4)') ne
                CALL WRITES(IODI,OP_STRING,ERROR,*200)
              ENDIF
C GMH 30/10/95 Adding call to orthog projection procedure, accounting
C for deformed coordinates correctly (I think)
C GMH 30/10/95 Should NBJ be NBH for deformed?
              IF(DEFORM)THEN
c             Note that deformed coords --> XE
                CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     '            NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),NRE(ne),
     '            NVHE(1,1,1,ne),NW(ne,1),nx,
     '            CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '            XE,ZP,ERROR,*200)
              ELSE
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '            NPF(1,1),NPNE(1,1,ne),
     '            NRE(ne),NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*200)
              ENDIF

              FOUND=.TRUE. !find nearest point in element
              CALL PROJ_ORTHOG(IBT,IDO,INP,NBJ(1,ne),SQND,
     '          XE,XI,ZD(1,nd),FOUND,ERROR,*200)

              IF(.NOT.SPECIFY) THEN !try neighbouring element
                neold=nelast
                nelast=ne
                DO ni=1,NITB
                  IF(XI(ni).EQ.0.0d0) THEN
                    neadj=NXI(-ni,1,ne)
                    IF(neadj.GT.0) THEN
                      ne=neadj
                      XI(ni)=1.0d0
                    ENDIF
                  ELSE IF(XI(ni).EQ.1.0d0) THEN
                    neadj=NXI(ni,1,ne)
                    IF(neadj.GT.0) THEN
                      ne=neadj
                      XI(ni)=0.0d0
                    ENDIF
                  ENDIF !in bounds
                ENDDO
                FOUND=ne.EQ.neold.OR.ne.EQ.nelast
              ENDIF !not FOUND
              IF(FOUND) THEN
                DO ni=1,NITB !copy new xi values
                  XID(ni,nd)=XI(ni)
                ENDDO
                LD(nd)=ne
                SQ(nd)=SQND
              ENDIF !FOUND

            ENDDO ! until found
          ENDIF
          GO TO 202
C             This statement is designed to be skipped if no error
C             occurs. However if a error occurs within a subroutine
C             the alternate return points to line 200 to set the flag
 200      CONTINUE
C$OMP CRITICAL(DEXI_CLOSEST3)
          ERROR_FLAG=.TRUE.
          WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
          CALL WRITES(IOER,OP_STRING,ERROR,*201)
          WRITE(OP_STRING,'(/'' >>An error occurred - '
     '      //'results may be unreliable!'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*201)
 201      CONTINUE
C$OMP END CRITICAL(DEXI_CLOSEST3)
 202      CONTINUE

        ENDIF !error_flag

C MPN 7-NOV-2014: use SQMAX to omit projections with magnitude greater
C than MAX
        IF(SQ(nd).GT.SQMAX) THEN
          LD(nd)=0
          SQ(nd)=0.0d0
          DO ni=1,NITB
            XID(ni,nd)=0.0d0
          ENDDO
        ENDIF

      ENDDO ! nd
C$OMP END PARALLEL DO

      CALL EXITS('DEXI_CLOSEST')
      RETURN
 9999 CALL ERRORS('DEXI_CLOSEST',ERROR)
      CALL EXITS('DEXI_CLOSEST')
      RETURN 1
      END



