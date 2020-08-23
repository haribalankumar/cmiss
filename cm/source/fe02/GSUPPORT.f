      SUBROUTINE GSUPPORT(NLATNE,NLATNQ,NLATPNQ,NQGP,NQGP_PIVOT,NQNLAT,
     &  NQS,NQXI,NRLIST,NWQ,ERROR,*)

C#### Subroutine: GSUPPORT
C###  Description:
C###    GSUPPORT calculates the support for the lattice grid points.
C###    This support is generated using templates which match the
C###    support arrangement used in the collocation method. The
C###    support points are stored in ascending numerical order this
C###    meaning that no pivot information needs to be stored.      
C###    As the lattice method is element based, support points
C###    need to be implemented for a finite difference scheme
C###    by applying spatially related weightings. This is done
C###    by applying a Moore-Penrose Inverse.
C### See-Also CALC_LATTICE_WEIGHTS      
      
      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'    
!     Parameter List      
      INTEGER NLATNE(NEQM+1),NLATNQ(NEQM*NQEM),NLATPNQ(NQM),
     &  NQGP(0:NQGM,NQM),NQGP_PIVOT(NQGM,NQM),NQNLAT(NEQM*NQEM),
     &  NQS(NEQM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NWQ(8,0:NQM)
      CHARACTER ERROR*(*)
      
!     Local Variables
      INTEGER count,first,grid_offset,i,j,k,loc,mi,mj,mk,ne,ni,nj,nk,
     &  nlat,nq,nq1,nqsc,nr,nrr,p

      CALL ENTERS('GSUPPORT',*9999)

      DO nrr=1,NRLIST(0)
        nr=NRLIST(nrr)
        DO nq=NQR(1,nr),NQR(2,nr) !create support for each nq
C         Work with the principal point first.
          loc=NLATPNQ(nq)
          
C         Initialise some values.
          NQGP(0,nq)=0
          DO count=1,NQGM
            NQGP(count,nq)=0        
          ENDDO

          DO WHILE(loc.NE.0)
C           Find the index of the NLATNE entry that is less than or
C           equal to nlat. ie, which element we are in.
            CALL FIND_LATT_NEIJK(i,j,k,ne,loc,NLATNE,NQXI,NQS,nqsc,
     &        ERROR,*9999)
            ni=NQXI(1,nqsc)
            nj=NQXI(2,nqsc)
            nk=NQXI(3,nqsc)
            DO mk=-1,1
              DO mj=-1,1
                DO mi=-1,1
                  IF(i+mi.GE.1.AND.i+mi.LE.ni.AND.
     &               j+mj.GE.1.AND.j+mj.LE.nj.AND.
     &               k+mk.GE.1.AND.k+mk.LE.nk) THEN
                    grid_offset=(mk*nj+mj)*ni+mi
                    nlat=loc+grid_offset
                    count=1
                    DO WHILE(count.LE.NQGP(0,nq).AND.
     &                NQGP(count,nq).NE.NQNLAT(nlat))
                      count=count+1
                    ENDDO
                    IF(NQGP(count,nq).NE.NQNLAT(nlat)) THEN
                      NQGP(0,nq)=NQGP(0,nq)+1
                      CALL ASSERT(NQGP(0,nq).LE.NQGM,
     &                  '>>Increase NQGM',ERROR,*9999)
                      NQGP(NQGP(0,nq),nq)=NQNLAT(nlat)
                    ENDIF
                  ENDIF
                ENDDO !mi
              ENDDO !mj
            ENDDO !mk
            loc=NLATNQ(loc)
          ENDDO !MORE

C***      Add all supporting points of the supporting points if on the bdy
C***      NOT DONE AT THE MOMENT
          IF(NWQ(1,nq).NE.0.AND..FALSE.) THEN
            first=NQGP(0,nq)
            DO p=1,first
              nq1=NQGP(p,nq)
              IF(nq1.NE.nq) THEN
C               Work with the principal point first
                loc=NLATPNQ(nq1)

                DO WHILE(loc.NE.0)
C                 Find the index of the NLATNE entry that is less than or
C                 equal to nlat. ie, which element we are in.
                  CALL FIND_LATT_NEIJK(i,j,k,ne,loc,NLATNE,NQXI,NQS,
     &              nqsc,ERROR,*9999)
                  ni=NQXI(1,nqsc)
                  nj=NQXI(2,nqsc)
                  nk=NQXI(3,nqsc)
                  DO mk=-1,1
                    DO mj=-1,1
                      DO mi=-1,1
                        IF(i+mi.GE.1.AND.i+mi.LE.ni.AND.
     &                     j+mj.GE.1.AND.j+mj.LE.nj.AND.
     &                     k+mk.GE.1.AND.k+mk.LE.nk) THEN
                          grid_offset=(mk*nj+mj)*ni+mi
                          nlat=loc+grid_offset
                          count=1
                          DO WHILE(count.LE.NQGP(0,nq).AND.
     &                      NQGP(count,nq).NE.NQNLAT(nlat))
                            count=count+1
                          ENDDO
                          IF(NQGP(count,nq).NE.NQNLAT(nlat)) THEN
                            NQGP(0,nq)=NQGP(0,nq)+1
                            CALL ASSERT(NQGP(0,nq).LE.NQGM,
     &                        '>>Increase NQGM',ERROR,*9999)
                            NQGP(NQGP(0,nq),nq)=NQNLAT(nlat)
                          ENDIF
                        ENDIF
                      ENDDO !mi
                    ENDDO !mj
                  ENDDO !mk
                  loc=NLATNQ(loc)
                ENDDO
              ENDIF !not original nq   
            ENDDO !support
          ENDIF !external

C         Order the support. This should not be necessary but is being
C         done for compatibility with current CMISS grid support.
          CALL ISORTP(NQGP(0,nq),NQGP(1,nq),NQGP_PIVOT(1,nq))
        ENDDO !nq
      ENDDO !nr

      CALL EXITS('GSUPPORT')
      RETURN
 9999 CALL ERRORS('GSUPPORT',ERROR)
      CALL EXITS('GSUPPORT')
      RETURN 1      
      END
