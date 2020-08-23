      SUBROUTINE SURFACELAPLACIAN(IBT,LIST_RESID,
     '  NBH,NENP,NPLIST3,NPNE,nx,NXI,NYNP,
     '  RESID,RESJAC,TOT_LAPL,XP,YP,
     '  ERROR,*)

C#### Subroutine: SURFACELAPLACIAN
C###  Description:
C###    Calculates the surface laplcian for regularising inverse
C###    activation problems.

      IMPLICIT NONE
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
      INTEGER IBT(3,NIM,NBFM),LIST_RESID,
     '  NBH(NHM,NCM,NEM),NENP(NPM,0:NEPM,0:NRM),NPLIST3(0:NPM),
     '  NPNE(NNM,NBFM,NEM),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 RESID(*),RESJAC(NREM,*),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER icol,nb,nb1,ne,NE1,nh1,ni,nk1,
     '  nolist1,noopti,nores,np,NP1,NP2,NP3,
     '  nv1,ny1
      REAL*8 CENT_DIFF2,delta1,delta2,SUM2,SUM3,
     '  TOT_LAPL,yp0,yp1,yp2
      LOGICAL FOUND
!     Functions
      REAL*8 NODE_DIST,SQ_DIST
      LOGICAL INLIST

      CALL ENTERS('SURFACELAPLACIAN',*9999)

C*** assuming a linear field with 1 dependent variable
      nores=NT_RES !the last residual
      nk1=1
      nv1=1
      nh1=1
      noopti=1


C*** Calculating the neighbouring nodes
      DO nolist1=1,NPLIST3(0) ! heart nodes (NOPTI)
        np=NPLIST3(nolist1)
        ne1=1
        FOUND=.FALSE.
        DO WHILE(.NOT.FOUND
     '    .AND.ne1.LE.NENP(np,0,trsf_nr_first))

          ne=NENP(np,ne1,trsf_nr_first)
          nb=NBH(NH_LOC(1,nx),1,ne)
          IF(np.EQ.NPNE(1,nb,ne)) THEN !First node of ne
            FOUND=.TRUE.
            ne=NENP(np,ne1,trsf_nr_first)
          ENDIF
          ne1=ne1+1
        ENDDO

        IF(.NOT.FOUND) THEN
!              WRITE(OP_STRING,'(''ERROR: No element found for node '',
!     '          I5,'' - sector node'' )') np
!              CALL WRITES(IOER,OP_STRING,ERROR,*9999)

          IF(LIST_RESID.GE.2) THEN
            WRITE(OP_STRING,'(''Closing sector node '',I5)') np
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          RESJAC(nores,nolist1)=0.0d0 ! at top
        ELSE
          nb=NBH(NH_LOC(1,nx),1,ne)

C*** Exclude sectors (only exclude the bottom ring as the top
C*** ring of elements is still required as we are picking the node
C*** at the (0,0) position.
C*** ASSUMES a certain mesh layout.

          sum2=0.0d0 !sum of jacobian of residuals

          IF((IBT(1,1,nb).EQ.5).OR.(IBT(1,2,nb).EQ.5)) THEN
            IF(LIST_RESID.GE.3) THEN
              WRITE(OP_STRING,
     '          '(''Excluding: ne,np* '',2I5,'' Sector at xi=0'')')
     '          ne,np
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!                  WRITE(*,*) 'Excluding: ne,npstar',ne,np,
!     '              ' Sector at xi=0'
            ENDIF
            sum2=0.0d0 !at bottom
          ELSE
            IF(INLIST(NPNE(1,nb,ne),NPLIST3(1),NPLIST3(0),icol)) THEN

              IF(LIST_RESID.GE.3) THEN
                WRITE(OP_STRING(1),
     '            '(''Searching: ne,npstar '',2I5)') ne,np
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF


C*** The center node
              ny1=NYNP(nk1,nv1,nh1,np,0,1,trsf_nr_first)
              yp0=YP(ny1,1)

C*** For each xi direction
              DO ni=1,NIT(nb) !should be 2 (BEM surface)
                np1=NPNE(ni+1,nb,ne) !positive dirn
                ny1=NYNP(nk1,nv1,nh1,np1,0,1,trsf_nr_first)
                yp1=YP(ny1,1)

                ne1=NXI(-ni,1,ne) !negative dirn
                nb1=NBH(NH_LOC(1,nx),1,ne1)
                np2=NPNE(1,nb1,ne1)
                ny1=NYNP(nk1,nv1,nh1,np2,0,1,trsf_nr_first)
                yp2=YP(ny1,1)

                IF(LIST_RESID.GE.4) THEN
                  WRITE(OP_STRING(1),'(''  Cent2 nodes '',3I5)')
     '              np2,np,np1
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF

                delta1=NODE_DIST(nk1,nk1,np,np1,nv1,nv1,XP)
                delta2=NODE_DIST(nk1,nk1,np,np2,nv1,nv1,XP)

C*** Residual - Squared central differences
                TOT_LAPL=TOT_LAPL
     '            +CENT_DIFF2(yp1,yp2,yp0,delta1,delta2)

C*** Jacobian of Residual
                sum3=SQ_DIST(np,np1,XP) ! Surrounding nodes
                sum3=sum3+SQ_DIST(np,np2,XP)
                sum2=sum2+(-4)/sum3

C*** Nodes in positve direction
                IF(ni.EQ.2..AND.
     '            (IBT(1,1,nb).EQ.6).OR.(IBT(1,2,nb).EQ.6)) THEN
                  np3=NPNE(3,nb,ne) ! Top elements use ne=ne1
                  IF(LIST_RESID.GE.3) THEN
!                        WRITE(*,*) '  Sector element: Top Ring'
!                        WRITE(*,*) '  ** Ignoring Resjac **'

                    WRITE(OP_STRING(1),
     '                '(''    Top Ring Sector elem'')')
                    WRITE(OP_STRING(2),
     '                '(''    ** IGNORING Resjac'')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSE
                  ne1=NXI(ni,1,ne)
                  nb1=NBH(NH_LOC(1,nx),1,ne1)
                  IF(ni.EQ.2..AND.
     '              (IBT(1,1,nb1).EQ.6).OR.(IBT(1,2,nb1).EQ.6)) THEN
                    IF(LIST_RESID.GE.3) THEN
!                          WRITE(*,*) '  Sector element: ',
!     '                      'Not shuffling up - 2nd to top'
                      WRITE(OP_STRING(1),
     '                  '(''    Sector elem - 2nd to top ring '')')
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    np3=NPNE(3,nb1,ne1)
                  ELSE
                    ne1=NXI(ni,1,ne1) !shuffle across 2
                    nb1=NBH(NH_LOC(1,nx),1,ne1)
                    np3=NPNE(1,nb1,ne1)
                  ENDIF
                  IF(LIST_RESID.GE.4) THEN
                    WRITE(OP_STRING(1),'(''  Resjac nodes+ '',3I5)')
     '                np,np1,np3
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  sum3=SQ_DIST(np,np1,XP)
                  sum3=sum3+SQ_DIST(np1,np3,XP)
                  sum2=sum2+2/sum3
                ENDIF

*** Nodes in negative direction
                ne1=NXI(-ni,1,ne)
                nb1=NBH(NH_LOC(1,nx),1,ne1)
                IF((IBT(1,1,nb1).EQ.5).OR.(IBT(1,2,nb1).EQ.5)) THEN
                  IF(LIST_RESID.GE.3) THEN
!                        WRITE(*,*) '  Sector element: Not suffling down'
!                        WRITE(*,*) '  *** ignoring resjac ***'
                    WRITE(OP_STRING(1),
     '                '(''    Sector elem: Not shuffling down'')')
                    WRITE(OP_STRING(2),
     '                '(''    ** IGNORING Resjac'')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  ne1=ne
                  nb1=NBH(NH_LOC(1,nx),1,ne1)
                ELSE
                  ne1=NXI(-ni,1,ne1) !shuffle across 2
                  nb1=NBH(NH_LOC(1,nx),1,ne1)
                  np3=NPNE(1,nb1,ne1)

                  IF(LIST_RESID.GE.4) THEN
                    WRITE(OP_STRING(1),'(''  Resjac nodes- '',3I5)')
     '                np,np2,np3
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  sum3=SQ_DIST(np,np2,XP)
                  sum3=sum3+SQ_DIST(np2,np3,XP)
                  sum2=sum2+2/sum3
                ENDIF
              ENDDO !ni

            ELSE
              IF(LIST_RESID.GE.2) THEN
                WRITE(OP_STRING(1),
     '            '(''Excluding: ne,np* '',2I5,'' NOT inlist'')')
     '            ne,np
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF !INLIST
          ENDIF ! Exclude sectors
          noopti=noopti+1
!              RESJAC(nores,nolist1)=noopti !test filling pattern
          RESJAC(nores,nolist1)=sum2*ACTN_REG_PARAM_LAPLACE
        ENDIF !Found
      ENDDO !np
      RESID(nores)=TOT_LAPL*ACTN_REG_PARAM_LAPLACE

      CALL EXITS('SURFACELAPLACIAN')
      RETURN
 9999 CALL ERRORS('SURFACELAPLACIAN',ERROR)
      CALL EXITS('SURFACELAPLACIAN')
      RETURN 1
      END


