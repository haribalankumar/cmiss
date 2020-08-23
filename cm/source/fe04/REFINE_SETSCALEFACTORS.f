      SUBROUTINE REFINE_SETSCALEFACTORS(IBT,IDO,IDRN,ne,ne_new,NIT1,SE,
     '  XII,ERROR,*)

C#### Subroutine: REFINE_SETSCALEFACTORS
C###  Description:
C###    REFINE_SETSCALEFACTORS adjusts the scale factors in the
C###    element that has been refined and the new element created in
C###    the refine process.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),IDRN,ne,ne_new,
     '  NIT1
      REAL*8 SE(NSM,NBFM,NEM),XII
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,nb1,nk,nn,ns,ns2

      CALL ENTERS('REFINE_SETSCALEFACTORS',*9999)

C CS 9 Oct 98 new If refinement is done on a mesh containing hanging
      DO nb1=1,NBFT
        IF(NKT(0,nb1).GT.1) THEN
          IF((IBT(1,IDRN,nb1).EQ.2).OR.((IBT(1,IDRN,nb1).EQ.5.OR.
     '        IBT(1,IDRN,nb1).EQ.6).AND.IBT(2,IDRN,nb1).EQ.4)) THEN
            IF(NBI(nb1).EQ.1) THEN !unit scale factors
C             Change scale factors from unit.
              NBI(nb1)=2 !change to element scale factors
              WRITE(OP_STRING,
     '          '('' >>Warning: scale factors are no longer unit'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            IF(NBI(nb1).EQ.2.OR.JTYP2C.EQ.1) THEN
!            IF(NBI(nb1).EQ.2) THEN !element scale factors
C CS 16/7/2001 This is not correct for lower dimensions
C CS 30/1/2004 Not sure why I changed my mind in 16/7/2001 comment
C              reverting back to allowing all dimensions.              
C KAT 13Jan00: Trying for lower dimensions
C KAT 30/1/2004: Seems to work for simple 2d cases and seems better than
C                not setting scale factors at all.
C               IF(NIT(nb1).EQ.3.OR.IDRN.EQ.3) THEN
              ! Only set scale factors on bases that are of the same
              ! dimension as the element.  (There is no point and the old
              ! scale factors may not have been set.)
              IF(NIT(nb1).EQ.NIT1) THEN
C CS July 99  This should work for 2D and 3D elements with 4 or 8 nodes
C             Has only been thoughly tested for 3D tricubic hermite
                ! Set up scale factors on old boundary of new element,
                ! based on the old scale factors on that boundary. 
                ns=0
                DO nn=1,NNT(nb1)
C                  DO nk=1,NKT(nn,nb1)
                  DO nk=1,NKT(nn,nb1)
                    ns=ns+1
                    IF(nk.EQ.1) THEN
                      SE(ns,nb1,ne_new)=1.0d0
                    ELSE IF(IDO(nk,nn,IDRN,nb1).EQ.1) THEN !not deriv in refn dirn
                      IF(IDRN.EQ.1) THEN
                        IF(MOD(nn,2).EQ.0) THEN
                          SE(ns,nb1,ne_new)=SE(ns,nb1,ne)
                        ENDIF
                      ELSE IF(IDRN.EQ.2) THEN
                        IF((nn.EQ.3).OR.(nn.EQ.4)
     '                    .OR.(nn.EQ.7).OR.(nn.EQ.8)) THEN
                          SE(ns,nb1,ne_new)=SE(ns,nb1,ne)
                        ENDIF
                      ELSE
                        IF((nn.EQ.5).OR.(nn.EQ.6)
     '                    .OR.(nn.EQ.7).OR.(nn.EQ.8)) THEN
                          SE(ns,nb1,ne_new)=SE(ns,nb1,ne)
                        ENDIF
                      ENDIF
                    ELSE
                      IF(IDRN.EQ.1) THEN
                        IF(MOD(nn,2).EQ.0) THEN
                          SE(ns,nb1,ne_new)=SE(ns,nb1,ne)*(1.0d0-XII)
                        ENDIF
                      ELSE IF(IDRN.EQ.2) THEN
                        IF((nn.EQ.3).OR.(nn.EQ.4)
     '                    .OR.(nn.EQ.7).OR.(nn.EQ.8)) THEN
                          SE(ns,nb1,ne_new)=SE(ns,nb1,ne)*(1.0d0-XII)
                        ENDIF
                      ELSE
                        IF((nn.EQ.5).OR.(nn.EQ.6)
     '                    .OR.(nn.EQ.7).OR.(nn.EQ.8)) THEN
                          SE(ns,nb1,ne_new)=SE(ns,nb1,ne)*(1.0d0-XII)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
                ns=0
                ! Set up scale factors on new boundary of old element,
                ! based on the new scale factors on the equivalent boundary
                ! of the new element,
                ! and, on old boundary of old element,
                ! based on the old scale factors on that boundary. 
                DO nn=1,NNT(nb1) ! Adjust old element
                  DO nk=1,NKT(nn,nb1)
                    ns=ns+1
                    IF(IDO(nk,nn,IDRN,nb1).EQ.1) THEN !not deriv in refn dirn
                      IF(IDRN.EQ.1) THEN
                        IF(nn.EQ.2) THEN
                          ns2=nk
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ELSE IF(nn.EQ.4) THEN
                          ns2=nk
                          DO i=1,2
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ELSE IF(nn.EQ.6) THEN
                          ns2=nk
                          DO i=1,4
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ELSE IF(nn.EQ.8) THEN
                          ns2=nk
                          DO i=1,6
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ENDIF
                      ELSE IF(IDRN.EQ.2) THEN
                        IF(nn.EQ.3) THEN
                          ns2=nk
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ELSE IF(nn.EQ.4) THEN
                          ns2=NKT(1,nb1)+nk
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ELSE IF(nn.EQ.7) THEN
                          ns2=nk
                          DO i=1,4
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ELSE IF(nn.EQ.8) THEN
                          ns2=nk
                          DO i=1,5
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ENDIF
                      ELSE
                        IF(nn.EQ.5) THEN
                          ns2=nk
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ELSE IF(nn.EQ.6) THEN
                          ns2=NKT(1,nb1)+nk
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ELSE IF(nn.EQ.7) THEN
                          ns2=nk
                          DO i=1,2
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ELSE IF(nn.EQ.8) THEN
                          ns2=nk
                          DO i=1,3
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
                        ENDIF
                      ENDIF
                    ELSE
                      IF(IDRN.EQ.1) THEN
                        IF(nn.EQ.2) THEN
                          ns2=nk
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE IF(nn.EQ.4) THEN
                          ns2=nk
                          DO i=1,2
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE IF(nn.EQ.6) THEN
                          ns2=nk
                          DO i=1,4
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE IF(nn.EQ.8) THEN
                          ns2=nk
                          DO i=1,6
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE
                          SE(ns,nb1,ne)=XII*SE(ns,nb1,ne)
                        ENDIF
                      ELSE IF(IDRN.EQ.2) THEN
                        IF(nn.EQ.3) THEN
                          ns2=nk
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE IF(nn.EQ.4) THEN
                          ns2=NKT(1,nb1)+nk
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE IF(nn.EQ.7) THEN
                          ns2=nk
                          DO i=1,4
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE IF(nn.EQ.8) THEN
                          ns2=nk
                          DO i=1,5
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE
                          SE(ns,nb1,ne)=XII*SE(ns,nb1,ne)
                        ENDIF
                      ELSE
                        IF(nn.EQ.5) THEN
                          ns2=nk
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE IF(nn.EQ.6) THEN
                          ns2=NKT(1,nb1)+nk
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE IF(nn.EQ.7) THEN
                          ns2=nk
                          DO i=1,2
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE IF(nn.EQ.8) THEN
                          ns2=nk
                          DO i=1,3
                            ns2=ns2+NKT(i,nb1)
                          ENDDO
                          SE(ns,nb1,ne)=
     '                      XII*SE(ns2,nb1,ne_new)/(1.0d0-XII)
                        ELSE
                          SE(ns,nb1,ne)=XII*SE(ns,nb1,ne)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO

C CS JULY 99 I think the above code and new SE setting code in
C REFINE_SETNODE should supersede this.
C This is not accurate anyway.
C To be removed once new code has been put though its paces.
C              ELSE
C                ns=0
C                DO nn=1,NNT(nb1)
C                  DO nk=1,NKT(nn,nb1)
C                    ns=ns+1
C                    IF(IDO(nk,nn,IDRN,nb1).EQ.1) THEN !not deriv in refine dirn
C CS new 21.12.98 this only works for 2D BEM biqubic elements at the
C moment
C                      IF((ns.EQ.1).AND.(IDRN.EQ.1)) ns2=5
C                      IF((ns.EQ.3).AND.(IDRN.EQ.1)) ns2=7
C                      IF((ns.EQ.9).AND.(IDRN.EQ.1)) ns2=13
C                      IF((ns.EQ.11).AND.(IDRN.EQ.1)) ns2=15
C
C                      IF((ns.EQ.1).AND.(IDRN.EQ.2)) ns2=9
C                      IF((ns.EQ.2).AND.(IDRN.EQ.2)) ns2=10
C                      IF((ns.EQ.5).AND.(IDRN.EQ.2)) ns2=13
C                      IF((ns.EQ.6).AND.(IDRN.EQ.2)) ns2=14

C                      IF(IDRN.EQ.2) THEN
C                        IF(nn.LE.2) THEN
C                          SE(ns,nb1,ne_new)=SE(ns,nb1,ne)*XII
C     '                      +SE(ns2,nb1,ne)*(1.0d0-XII)
C                        ELSE
C                          SE(ns,nb1,ne_new)=SE(ns,nb1,ne)
C                        ENDIF
C                      ELSE
C                        IF((nn.EQ.1).OR.(nn.EQ.3)) THEN
C                          SE(ns,nb1,ne_new)=SE(ns,nb1,ne)*XII
C     '                      +SE(ns2,nb1,ne)*(1.0d0-XII)
C                        ELSE
C                          SE(ns,nb1,ne_new)=SE(ns,nb1,ne)
C                        ENDIF
C                      ENDIF
C                    ELSE
C                      SE(ns,nb1,ne_new)=SE(ns,nb1,ne)*XII
C                    ENDIF
C                  ENDDO
C                ENDDO
C                ns=0
C                DO nn=1,NNT(nb1)
C                  DO nk=1,NKT(nn,nb1)
C                    ns=ns+1
C                    IF(IDO(nk,nn,IDRN,nb1).EQ.1) THEN !not deriv in refine dirn
C                      IF((ns.EQ.5).AND.(IDRN.EQ.1)) ns2=1
C                      IF((ns.EQ.7).AND.(IDRN.EQ.1)) ns2=3
C                      IF((ns.EQ.13).AND.(IDRN.EQ.1)) ns2=9
C                      IF((ns.EQ.15).AND.(IDRN.EQ.1)) ns2=11

C                      IF((ns.EQ.9).AND.(IDRN.EQ.2)) ns2=1
C                      IF((ns.EQ.10).AND.(IDRN.EQ.2)) ns2=2
C                      IF((ns.EQ.13).AND.(IDRN.EQ.2)) ns2=5
C                      IF((ns.EQ.14).AND.(IDRN.EQ.2)) ns2=6

C                      IF(IDRN.EQ.2) THEN
C                        IF(nn.GT.2) THEN
C                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
C                        ENDIF
C                      ELSE
C                        IF((nn.EQ.2).or.(nn.eq.4)) THEN
C                          SE(ns,nb1,ne)=SE(ns2,nb1,ne_new)
C                        ENDIF
C                      ENDIF
C                    ELSE
C                      SE(ns,nb1,ne)=SE(ns,nb1,ne)*(1.0d0-XII)
C                    ENDIF
C                  ENDDO
C                ENDDO

C KAT 7Dec98: Otherwise scale factors are calculated in lincal (from refine)
C            ELSE IF(NKT(0,nb).LE.2) THEN
C              DO nn=1,NNT(nb1)
C                ns=(nn-1)*NKT(0,nb1)+2
C                SE(ns,nb1,ne)=SE(ns,nb1,ne)/2.0d0
C              ENDDO !ns
C            ELSE IF(NKT(0,nb1).EQ.4) THEN
C              IF(NIT(nb1).EQ.2) THEN
C                NTOT=1
C              ELSE IF(NIT(nb1).EQ.3) THEN
C                NTOT=2
C              ENDIF
C              DO n=1,NTOT
C                k=16*(n-1)
C                IF(IDRN.EQ.1) THEN
CC ***             Xi(1) derivs
C                  SE( 2+k,nb1,ne_NEW)=0.25D0*(SE( 2+k,nb1,ne)
C     '              +SE( 6+K,nb1,ne))
C                  SE( 6+k,nb1,ne_NEW)=0.50D0* SE( 6+k,nb1,ne)
C                  SE(10+k,nb1,ne_NEW)=0.25D0*(SE(10+k,nb1,ne)
C     '              +SE(14+k,nb1,ne))
C                  SE(14+k,nb1,ne_NEW)=0.50D0* SE(14+k,nb1,ne)
C                  SE( 6+k,nb1,ne)    =0.25D0*(SE( 2+k,nb1,ne)
C     '              +SE( 6+k,nb1,ne))
C                  SE( 2+k,nb1,ne)    =0.50D0* SE( 2+k,nb1,ne)
C                  SE(14+k,nb1,ne)    =0.25D0*(SE(10+k,nb1,ne)
C     '              +SE(14+k,nb1,ne))
C                  SE(10+k,nb1,ne)    =0.50D0* SE(10+k,nb1,ne)
CC ***             Xi(2) derivs
C                  SE( 3+k,nb1,ne_NEW)=0.50D0*(SE( 3+k,nb1,ne)
C     '              +SE( 7+k,nb1,ne))
C                  SE( 7+k,nb1,ne_NEW)=      SE( 7+k,nb1,ne)
C                  SE(11+k,nb1,ne_NEW)=0.50D0*(SE(11+k,nb1,ne)
C     '              +SE(15+k,nb1,ne))
C                  SE(15+k,nb1,ne_NEW)=      SE(15+k,nb1,ne)
C                  SE( 7+k,nb1,ne)    =0.50D0*(SE( 3+k,nb1,ne)
C     '              +SE( 7+k,nb1,ne))
C                  SE( 3+k,nb1,ne)    =      SE( 3+k,nb1,ne)
C                  SE(15+k,nb1,ne)    =0.50D0*(SE(11+k,nb1,ne)
C     '              +SE(15+k,nb1,ne))
C                  SE(11+k,nb1,ne)    =      SE(11+k,nb1,ne)
C                ELSE IF(IDRN.EQ.2) THEN
CC ***             Xi(1) derivs
C                  SE( 2+k,nb1,ne_NEW)=0.50D0*(SE( 2+k,nb1,ne)
C     '              +SE(10+k,nb1,ne))
C                  SE( 6+k,nb1,ne_NEW)=0.50D0*(SE( 6+k,nb1,ne)
C     '              +SE(14+k,nb1,ne))
C                  SE(10+k,nb1,ne_NEW)=      SE(10+k,nb1,ne)
C                  SE(14+k,nb1,ne_NEW)=      SE(14+k,nb1,ne)
C                  SE(10+k,nb1,ne)    =0.50D0*(SE( 2+k,nb1,ne)
C     '              +SE(10+k,nb1,ne))
C                  SE(14+k,nb1,ne)    =0.50D0*(SE( 6+k,nb1,ne)
C     '              +SE(14+k,nb1,ne))
C                  SE( 2+k,nb1,ne)    =      SE( 2+k,nb1,ne)
C                  SE( 6+k,nb1,ne)    =      SE( 6+k,nb1,ne)
CC ***             Xi(2) derivs
C                  SE( 3+k,nb1,ne_NEW)=0.25D0*(SE( 3+k,nb1,ne)
C     '              +SE(11+k,nb1,ne))
C                  SE( 7+k,nb1,ne_NEW)=0.25D0*(SE( 7+k,nb1,ne)
C     '              +SE(15+k,nb1,ne))
C                  SE(11+k,nb1,ne_NEW)=0.50D0* SE(11+k,nb1,ne)
C                  SE(15+k,nb1,ne_NEW)=0.50D0* SE(15+k,nb1,ne)
C                  SE(11+k,nb1,ne)    =0.25D0*(SE( 3+k,nb1,ne)
C     '              +SE(11+k,nb1,ne))
C                  SE(15+k,nb1,ne)    =0.25D0*(SE( 7+k,nb1,ne)
C     '              +SE(15+k,nb1,ne))
C                  SE( 3+k,nb1,ne)    =0.50D0* SE( 3+k,nb1,ne)
C                  SE( 7+k,nb1,ne)    =0.50D0* SE( 7+k,nb1,ne)
C                ELSE IF(IDRN.EQ.3) THEN
C                  ERROR='>>Not implemented'
C                  GOTO 9999
C                ENDIF
CC ***           Xi(1),Xi(2) derivs
C                SE( 4+k,nb1,ne_NEW)=SE( 2+k,nb1,ne_NEW)
C     '            *SE( 3+k,nb1,ne_NEW)
C                SE( 8+k,nb1,ne_NEW)=SE( 6+k,nb1,ne_NEW)
C     '            *SE( 7+k,nb1,ne_NEW)
C                SE(12+k,nb1,ne_NEW)=SE(10+k,nb1,ne_NEW)
C     '            *SE(11+k,nb1,ne_NEW)
C                SE(16+k,nb1,ne_NEW)=SE(14+k,nb1,ne_NEW)
C     '            *SE(15+k,nb1,ne_NEW)
C                SE( 4+k,nb1,ne)  =SE( 2+k,nb1,ne) *SE( 3+k,nb1,ne)
C                SE( 8+k,nb1,ne)  =SE( 6+k,nb1,ne) *SE( 7+k,nb1,ne)
C                SE(12+k,nb1,ne)  =SE(10+k,nb1,ne) *SE(11+k,nb1,ne)
C                SE(16+k,nb1,ne)  =SE(14+k,nb1,ne) *SE(15+k,nb1,ne)
C              ENDDO
C            ELSE
C              WRITE(OP_STRING,'('' >>Scale factors not '
C     '          //'calculated'')')
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

              ENDIF !Dimension
            ENDIF !calc scale factor
          ELSE !NBI,ibt
            ns=0
            DO nn=1,NNT(nb1)
              DO nk=1,NKT(nn,nb1)
                ns=ns+1
                SE(ns,nb1,ne_NEW)=SE(ns,nb1,ne)
              ENDDO
            ENDDO
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' nb1='',I2,'', SE(ns,nb1,ne)='','
     '        //'4(1X,D12.5),/:(23X,4(1X,D12.5)))') nb1,
     '        (SE(ns,nb1,ne),ns=1,NST(nb1))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !dop
        ELSE !NKT=1
          DO ns=1,NST(nb1)+NAT(nb1)
            SE(ns,nb1,ne_NEW)=SE(ns,nb1,ne)
          ENDDO
        ENDIF !NKT

      ENDDO !nb1
C      ENDIF !idrn

      CALL EXITS('REFINE_SETSCALEFACTORS')
      RETURN
 9999 CALL ERRORS('REFINE_SETSCALEFACTORS',ERROR)
      CALL EXITS('REFINE_SETSCALEFACTORS')
      RETURN 1
      END




