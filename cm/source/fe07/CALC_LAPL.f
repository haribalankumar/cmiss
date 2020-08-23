      SUBROUTINE CALC_LAPL(IBT,NBH,NENP,NPLIST3,NPNE,nx,NXI,
     '  LAPL,LAPLSQR,XP,ERROR,*)

C#### Subroutine: CALC_LAPL
C###  Description:
C###    Calculates the Laplacian matrix (LAPL and LAPLSQR)
C###    for regularising inverse activation problems.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NPLIST3(0:NP_R_M),NPNE(NNM,NBFM,NEM),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 LAPL(NY_TRANSFER_M,NY_TRANSFER_M),
     '  LAPLSQR(NY_TRANSFER_M,NY_TRANSFER_M),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER nb,ncollapsed,ne,NE1,ne2,ni,nolist,nolist2,nolist3,
     '  np,np2,NPNXI(NPM,0:25),NPLIST_MAP(NPM)
      REAL*8 DIST,TOT_DIST,TOT_LAPL
      LOGICAL COLLAPSED

!     Functions
      REAL*8 SQ_DIST

      CALL ENTERS('CALC_LAPL',*9999)

      DO nolist=1,NPLIST3(0) !heart nodes (NOPTI)
        np=NPLIST3(nolist)
        NPLIST_MAP(np)=nolist !setup NPLIST_MAP -- reverse mapping

C***    Determine the number of collapsed elements surrounding the node
        ncollapsed=0
        DO ne1=1,NENP(np,0,trsf_nr_first)
          ne=NENP(np,ne1,trsf_nr_first)
          nb=NBH(NH_LOC(1,nx),1,ne)
          IF((IBT(1,1,nb).EQ.5).OR.(IBT(1,1,nb).EQ.6)) THEN
            ncollapsed=ncollapsed+1
          ENDIF
        ENDDO

C***    Considered a "collapsed" node if it is only surrounded
C***    by collapsed elements
        IF(ncollapsed.GE.NENP(np,0,trsf_nr_first)) THEN
          COLLAPSED=.TRUE.
        ELSE
          COLLAPSED=.FALSE.
        ENDIF

        IF(.NOT.COLLAPSED) THEN
          ne=-1
          DO ne1=1,NENP(np,0,trsf_nr_first) !elem surrounding np
            ne2=NENP(np,ne1,trsf_nr_first)
            nb=NBH(NH_LOC(1,nx),1,ne2)
            IF(NPNE(1,nb,ne2).EQ.np) THEN !first node of ne
              ne=ne2
            ENDIF
          ENDDO

C***      Compute nodes in xi1=+1 dirn
          ne1=ne
          nb=NBH(NH_LOC(1,nx),1,ne1)
          NPNXI(np,1)=NPNE(2,nb,ne1)

C***      Compute nodes in xi1=-1 dirn
          ne1=NXI(-1,1,ne)
          nb=NBH(NH_LOC(1,nx),1,ne1)
          NPNXI(np,2)=NPNE(1,nb,ne1)

C***      Compute nodes in xi2=+1 dirn
          ne1=ne
          nb=NBH(NH_LOC(1,nx),1,ne1)
          NPNXI(np,3)=NPNE(3,nb,ne1)

C***      Compute nodes in xi2=-1 dirn
          ne1=NXI(-2,1,ne)
          nb=NBH(NH_LOC(1,nx),1,ne1)
          NPNXI(np,4)=NPNE(1,nb,ne1)
          NPNXI(np,0)=4

        ELSE
          CALL ASSERT(NENP(np,0,trsf_nr_first).LE.25,
     '      '>> Increase size of NPNXI in code',ERROR,*9999)
          nb=NBH(NH_LOC(1,nx),1,NENP(np,1,trsf_nr_first))
          IF(IBT(1,1,nb).EQ.5) THEN !collapsed at bottom
            ni=3
          ELSEIF(IBT(1,1,nb).EQ.6) THEN !collapsed at top
            ni=1
          ELSE
            ERROR='Inconsistent collapsed element'
            GOTO 9999
          ENDIF !collapsed top/bottom
          DO ne1=1,NENP(np,0,trsf_nr_first)
            ne=NENP(np,ne1,trsf_nr_first)
            NPNXI(np,ne1)=NPNE(ni,nb,ne)
          ENDDO
          NPNXI(np,0)=NENP(np,0,trsf_nr_first)
        ENDIF !collapsed
      ENDDO ! heart nodes

C*** NPNXI is the list of nodes surrounding node np

      IF(DOP) THEN
        WRITE(OP_STRING(1),'('' NPNXI Array: '')')
        WRITE(OP_STRING(2),'('' ------------ '')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        DO nolist=1,NPLIST3(0) !heart nodes (NOPTI)
          np=NPLIST3(nolist)
          WRITE(OP_STRING,'('' Node: '',I5,''   Neighbours: '',8I5)')
     '      np,(NPNXI(np,np2),np2=1,NPNXI(np,0))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF !DOP

C*** Calculate the Lapl matrix
      DO np=1,NPLIST3(0) !initialise
        DO np2=1,NPLIST3(0)
          LAPL(np,np2)=0.D0
          LAPLSQR(np,np2)=0.D0
        ENDDO !np2
      ENDDO !np

C***  Fill in the matrix elements
C***  Using equations from oostendorp:1989a
C***  Gives same values as old code where adjacent points are evenly spaced
      DO nolist=1,NPLIST3(0)
        np=NPLIST3(nolist)
        TOT_DIST=0.0d0
        TOT_LAPL=0.0d0
        DO nolist2=1,NPNXI(np,0)
          np2=NPNXI(np,nolist2)
          nolist3=NPLIST_MAP(np2)
          DIST=SQRT(SQ_DIST(np,np2,XP))
          TOT_DIST=TOT_DIST+DIST
          LAPL(nolist,nolist3)=1.0d0/DIST
          TOT_LAPL=TOT_LAPL+1.0d0/DIST
        ENDDO !np2
        LAPL(nolist,nolist)=-TOT_LAPL*4.0d0/TOT_DIST
        DO nolist2=1,NPNXI(np,0)
          np2=NPNXI(np,nolist2)
          nolist3=NPLIST_MAP(np2)
          LAPL(nolist,nolist3)=LAPL(nolist,nolist3)*4.0d0/TOT_DIST
        ENDDO !np2
C***    Old code - to be deleted
C        IF(NPNXI(np,0).EQ.4) THEN !std config
C          delta1=0.D0 !horizontal distance
C          delta1=delta1+SQ_DIST(np,NPNXI(np,1),XP)
C          delta1=delta1+SQ_DIST(np,NPNXI(np,2),XP)
C
C          delta2=0.D0 !vertical distance
C          delta2=delta2+SQ_DIST(np,NPNXI(np,3),XP)
C          delta2=delta2+SQ_DIST(np,NPNXI(np,4),XP)
C
C          LAPL(nolist,nolist)=-4.D0/delta1 -4.D0/delta2
C          np2=NPLIST_MAP(NPNXI(np,1))
C          LAPL(nolist,np2)=2.D0/delta1
C          np2=NPLIST_MAP(NPNXI(np,2))
C          LAPL(nolist,np2)=2.D0/delta1
C          np2=NPLIST_MAP(NPNXI(np,3))
C          LAPL(nolist,np2)=2.D0/delta2
C          np2=NPLIST_MAP(NPNXI(np,4))
C          LAPL(nolist,np2)=2.D0/delta2
C        ELSE
CC*** Fill in the tops and bottom nodes latter ...
C          WRITE(OP_STRING(1),
C     '      '('' WARNING: Ignoring collapsed node number '',I6)') np
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        ENDIF ! std config
      ENDDO !np

C***  Compute LAPLSQR --  LAPL^T *LAPL
      CALL DGEMM('T','N', NPLIST3(0), NPLIST3(0), NPLIST3(0),
     '  1.D0, LAPL, NY_TRANSFER_M, LAPL, NY_TRANSFER_M, 1.D0,
     '  LAPLSQR, NY_TRANSFER_M)

      IF(DOP) THEN
        WRITE(OP_STRING(1),'('' '')')
        WRITE(OP_STRING(2),'('' LAPL Array: '')')
        WRITE(OP_STRING(3),'('' ------------ '')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO np=1,NPLIST3(0)
          IF(NPLIST3(0).GT.14) THEN
            WRITE(OP_STRING,'(''Row: '',I3)') np
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(10G12.4)')
     '        (LAPL(np,np2),np2=1,NPLIST3(0))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'(I3,'':'',14F8.3)') np,
     '        (LAPL(np,np2),np2=1,NPLIST3(0))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
        WRITE(OP_STRING(1),'('' '')')
        WRITE(OP_STRING(2),'('' LAPLSQR Array: '')')
        WRITE(OP_STRING(3),'('' -------------- '')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO np=1,NPLIST3(0)
          IF(NPLIST3(0).GT.14) THEN
            WRITE(OP_STRING,'(''Row:'', I3)') np
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(10G12.4)')
     '        (LAPLSQR(np,np2),np2=1,NPLIST3(0))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'(I3,'':'',14F8.3)') np,
     '        (LAPLSQR(np,np2),np2=1,NPLIST3(0))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF !DOP

      CALL EXITS('CALC_LAPL')
      RETURN
 9999 CALL ERRORS('CALC_LAPL',ERROR)
      CALL EXITS('CALC_LAPL')
      RETURN 1
      END


