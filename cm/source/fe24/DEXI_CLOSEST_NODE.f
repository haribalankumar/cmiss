      SUBROUTINE DEXI_CLOSEST_NODE(IBT,IDO,INP,NBJ,
     '  NBH,NELIST,NEP,NEW,NHE,nj1,
     '  NKHE,NKJE,NPF,NPNE,NPNODE,NRE,NRLIST2,NVHE,NVJE,NW,nx,NXI,
     '  START_ELEMENT,CURVCORRECT,SE,XA,XE,XI,
     '  XI_1,XI_2,XI_3,XIP,XP,ZA,ZP,
     '  ERROR,*)

C#### Subroutine: DEXI_CLOSEST_NODE
C###  Description:
C###    DEXI_CLOSEST_NODE calculates the closest point on any of the elements
C###    in NELIST for each nodes

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NBH(NHM,NCM,NEM),
     '  NELIST(0:NEM),NEP(NPM),NHE(NEM,NXM),nj1,NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRE(NEM),NRLIST2(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),START_ELEMENT
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XI(3),XI_1,XI_2,XI_3,
     '  XIP(NIM,NPM),XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL NEW
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IT,ITMAX,ix,nb,ne,neadj,NEC(NPM,3),nelast,neold,netemp,
     '  ni,NITB,nj,nk,no_ne,nonr,np,nplist,nr2,nv
      REAL*8 SQ1(NPM),SQND,TEMP,XITEMP(3),ZD(NJM)
      LOGICAL DEFORM,FOUND,ERROR_FLAG,SET_XI_1,SET_XI_2,SET_XI_3,
     '  SPECIFY
!     Functions
      REAL*8 PXI

      PARAMETER(ITMAX=20)

      CALL ENTERS('DEXI_CLOSEST_NODE',*9999)

C EWR 12-FEB-2003
C    This procedure is taken from DEXI_CLOSEST which finds xi-pos of
C    data-points in elements. It is modified so that it finds xi-pos
C    for nodes using the same algorythm. The procedure most likely 
C    needs some cleaning up as it probably have a lot of unnecessary
C    stuff, but then again if somebody later wants to add some more
C    functionalities I have let i stand as it is.
      ERROR_FLAG=.FALSE.
!      NEW=.TRUE.
      SPECIFY=.FALSE.
      DEFORM=.FALSE. 
!      SET_XI_1=.FALSE.
!      SET_XI_2=.FALSE.
      SET_XI_3=.FALSE.
!     Versions not taken care of
      nv=1
      nk=1
      nj1=1

      DO nonr=1,NRLIST2(0)
        nr2=NRLIST2(nonr)
        DO nplist=1,NPNODE(0,nr2)
          np=NPNODE(nplist,nr2)
          NEC(np,1)=0
          NEC(np,2)=0
          NEC(np,3)=0
        ENDDO ! nplist
      ENDDO ! nonr

      IF(NEW) THEN
C       Make a guess at the best element
        DO ni=1,3 !compare the distances to the centres of each element
          XI(ni)=0.5D0
        ENDDO
        IF(SET_XI_3) THEN
          XI(3)=XI_3
        ENDIF
        DO nonr=1,NRLIST2(0)
          nr2=NRLIST2(nonr)
c          noelem=0
          DO nplist=1,NPNODE(0,nr2)
            np=NPNODE(nplist,nr2)
            NEP(np)=0
            SQ1(np)=0.0D0
          ENDDO ! np
        ENDDO ! nonr


C$OMP   PARALLEL DO
C$OMP&  PRIVATE(no_ne,ne,NITB,XE,np,SQND,nj,nb,TEMP)
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
            DO nonr=1,NRLIST2(0)
              nr2=NRLIST2(nonr)
              DO nplist=1,NPNODE(0,nr2)
                np=NPNODE(nplist,nr2)
                SQND=0.0d0
                DO nj=1,NJT
                  nb=NBJ(nj,ne)
                  TEMP=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '              XI,XE(1,nj))-XP(nk,nv,nj,np)
                  SQND=SQND+TEMP*TEMP
                ENDDO
                IF(NEP(np).EQ.0.OR.SQND.LT.SQ1(np)) THEN
C$OMP CRITICAL(DEXI_CLOSEST_NODE1)
                  DO ni=1,NITB !copy new xi values
                    XIP(ni,np)=XI(ni)
                  ENDDO
                  NEC(np,3)=NEC(np,2) !NEC is closest element candidates to node
                  NEC(np,2)=NEC(np,1)
                  NEC(np,1)=NEP(np)
                  NEP(np)=ne
                  SQ1(np)=SQND
C$OMP END CRITICAL(DEXI_CLOSEST_NODE1)
                ENDIF
              ENDDO ! np
            ENDDO ! nonr
            GO TO 102
C             This statement is designed to be skipped if no error
C             occurs. However if a error occurs within a subroutine
C             the alternate return points to line 100 to set the flag
 100        CONTINUE
C$OMP CRITICAL(DEXI_CLOSEST_NODE2)
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
            CALL WRITES(IOER,OP_STRING,ERROR,*101)
            WRITE(OP_STRING,'(/'' >>An error occurred - '
     '        //'results may be unreliable!'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101        CONTINUE
C$OMP END CRITICAL(DEXI_CLOSEST_NODE2)
 102        CONTINUE

          ENDIF !error_flag
        ENDDO !no_ne
C$OMP   END PARALLEL DO

      ENDIF !NEW


C$OMP   PARALLEL DO
C$OMP&  PRIVATE(np,nelast,ne,NITB,XI,XITEMP,FOUND,IT,XE,SQND,neold,
C$OMP&          ni,neadj,ix,netemp)
      DO nonr=1,NRLIST2(0) !new AAY 23 March 95 limits set with ADD keyword
        nr2=NRLIST2(nonr)
c        noelem=0
        DO nplist=1,NPNODE(0,nr2)
          np=NPNODE(nplist,nr2)
          IF(.NOT.ERROR_FLAG) THEN
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Node '',I5)') np
              CALL WRITES(IODI,OP_STRING,ERROR,*200)
            ENDIF

C       Start with the guess element
            nelast=0
            netemp=0
            neadj=0
            ix=1
            XITEMP(1)=0.5d0
            XITEMP(2)=0.5d0
            XITEMP(3)=0.5d0

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
              ne=NEP(np)
              NITB=NIT(NBJ(NJ1,ne))
              DO ni=1,NITB
                XI(ni)=XIP(ni,np)
              ENDDO
            ENDIF
            FOUND=.FALSE.

            IF(ne.NE.0) THEN
              IT=0
              DO WHILE(.NOT.FOUND.AND.IT.LT.ITMAX)
                IT=IT+1
C             Search the elements around the guess element
                IF(DOP) THEN
                  WRITE(OP_STRING,'(/'' Element '',I4)') ne
                  CALL WRITES(IODI,OP_STRING,ERROR,*200)
                ENDIF
C GMH 30/10/95 Adding call to orthog projection procedure, accounting
C for deformed coordinates correctly (I think)
C GMH 30/10/95 Should NBJ be NBH for deformed?
                IF(DEFORM)THEN
c               Note that deformed coords --> XE
                  CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     '              NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),NRE(ne),
     '              NVHE(1,1,1,ne),NW(ne,1),nx,
     '              CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '              XE,ZP,ERROR,*200)
                ELSE
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '              NPF(1,1),NPNE(1,1,ne),
     '              NRE(ne),NVJE(1,1,1,ne),
     '              SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*200)
                ENDIF

                FOUND=.TRUE. !find nearest point in element
                DO nj=1,NJT
                  ZD(nj)=XP(nk,nv,nj,np)
                ENDDO
                CALL PROJ_ORTHOG(IBT,IDO,INP,NBJ(1,ne),SQND,
     '            XE,XI,ZD(1),FOUND,ERROR,*200)

                IF(.NOT.SPECIFY) THEN !try neighbouring element
                  neold=nelast
                  nelast=ne
                  DO ni=1,NITB
                    IF(XI(ni).EQ.0.0d0) THEN
                      neadj=NXI(-ni,1,ne)
                      IF(neadj.GT.0) THEN
                        ne=neadj
                        XI(ni)=1.0d0
                      ELSE
                        neadj=-1
                      ENDIF
                    ELSE IF(XI(ni).EQ.1.0d0) THEN
                      neadj=NXI(ni,1,ne)
                      IF(neadj.GT.0) THEN
                        ne=neadj
                        XI(ni)=0.0d0
                      ELSE
                        neadj=-1
                      ENDIF
                    ENDIF !in bounds
                  ENDDO
                  IF(ne.EQ.neold.OR.neadj.LT.0) THEN
                    neadj=0
                    IF(netemp.EQ.0) THEN
                      netemp=ne
                      DO ni=1,NITB !save xi values for netemp
                        XITEMP(ni)=XI(ni)
                      ENDDO
                    ENDIF
                    IF(ix.LT.4) THEN
                      ne=NEC(np,ix)
                      ix=ix+1
                      DO ni=1,3
                        XI(ni)=0.5D0
                      ENDDO
                    ELSE IF(ix.GE.4) THEN
                      ne=netemp !restore saved ne as this is probably better
                      DO ni=1,NITB !restore xi values from netemp
                        XI(ni)=XITEMP(ni)
                      ENDDO
                    ENDIF
                    IF(ne.EQ.0) THEN
                      ix=4
                      ne=netemp
                      DO ni=1,NITB !restore xi values from netemp
                        XI(ni)=XITEMP(ni)
                      ENDDO
                    ENDIF
                  ENDIF
                  FOUND=ne.EQ.neold.OR.ne.EQ.nelast
                ENDIF !not FOUND
                IF(FOUND) THEN
                  DO ni=1,NITB !copy new xi values
                    XIP(ni,np)=XI(ni)
                  ENDDO
                  NEP(np)=ne
                  SQ1(np)=SQND
                ENDIF !FOUND

              ENDDO ! until found
            ENDIF
            GO TO 202
C             This statement is designed to be skipped if no error
C             occurs. However if a error occurs within a subroutine
C             the alternate return points to line 200 to set the flag
 200        CONTINUE
C$OMP CRITICAL(DEXI_CLOSEST_NODE3)
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
            CALL WRITES(IOER,OP_STRING,ERROR,*201)
            WRITE(OP_STRING,'(/'' >>An error occurred - '
     '        //'results may be unreliable!'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*201)
 201        CONTINUE
C$OMP END CRITICAL(DEXI_CLOSEST_NODE3)
 202        CONTINUE

          ENDIF !error_flag
        ENDDO ! np
      ENDDO ! nonr
C$OMP END PARALLEL DO

      CALL EXITS('DEXI_CLOSEST_NODE')
      RETURN
 9999 CALL ERRORS('DEXI_CLOSEST_NODE',ERROR)
      CALL EXITS('DEXI_CLOSEST_NODE')
      RETURN 1
      END




