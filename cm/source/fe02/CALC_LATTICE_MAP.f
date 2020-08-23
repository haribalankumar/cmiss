      SUBROUTINE CALC_LATTICE_MAP(NEAM,NEES,NENQ,NEVS,NLATNE,NLATNQ,
     &  NLATPNQ,NQNE,NQNLAT,NQS,NQXI,ERROR,*)

C#### Subroutine: CALC_LATTICE_MAP
C###  Description:
C###    CALC_LATTICE_MAP takes the local element sharing arrays
C###    NEES and NEVS and builds the global lattic mapping arrays:
C###    NQNLAT, NLATPNQ and NLATNQ.      

C#### Variable: NLATNQ
C###  Type: INTEGER
C###  Set_up: CALC_LATTICE_MAP       
C###  Description:
C###    NLATNQ is a linked list that returns the lattice grid points
C###    that map onto global grid point nq. The first entry in the
C###    list for grid point nq is NLATNQ(NLATPNQ(nq)). This entry
C###    will either be zero if this is the last lattice point that
C###    maps onto nq, or will have the next position in NLATNQ for
C###    this nq.      
C###  See-Also: NLATPNQ
      
C#### Variable: NLATPNQ(nq)
C###  Type: INTEGER
C###  Set_up: CALC_LATTICE_MAP
C###  Description:      
C###    NLATPNQ(nq) returns the principal lattice grid point nlat
C###    for a given global grid point nq. The principal lattice
C###    grid point is the one which was first used in the construction
C###    of the global grid point nq by CALC_LATTICE_MAP.      
C###  See-Also: NLATNQ
      
C#### Variable: NQNLAT(nlat)
C###  Type: INTEGER
C###  Set_up: CALC_LATTICE_MAP
C###  Description:
C###    NQNLAT(nlat) returns the global grid point nq
C###    for a given lattice grid point nlat.      

      IMPLICIT NONE
      
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'      
!     Parameter List
      INTEGER NEAM,NEES(8,NEAM),NENQ(0:8,NQM),NEVS(0:12,NEAM),
     &  NLATNE(NEQM+1),NLATNQ(NEQM*NQEM),NQNE(NEQM,NQEM),
     &  NQNLAT(NEQM*NQEM),NLATPNQ(NQM),NQS(NEQM),NQXI(0:NIM,NQSCM) 
      CHARACTER ERROR*(*)      

!     Local Variables
      INTEGER adjacent_element,ADJN(3),counter,GPI,GPJ,GPK,
     &  HIGHLIM(3),i,IJK(8,3),IJK_NEW(8,3),IJKMAX(3),IJKMIN(3),
     &  j,k,LIM(3),line,llat,loc,LOWLIM(3),LPI,MAXS(3),min_grid,MINS(3),
     &  ne,ne_adj,ni,ni_adj,nlat,nqsc,shar_ne
      
!     Functions
      INTEGER INT_LATT_JUMP
      
      CALL ENTERS('CALC_LATTICE_MAP',*9999)

C     Initialise some variables.

      ne=NEVS(1,1)
      nqsc=NQS(ne)

      DO ni=1,3
        MINS(ni)=1
        MAXS(ni)=NQXI(ni,nqsc)
        ADJN(ni)=NQXI(ni,nqsc)
      ENDDO
      
      CALL VERT_LATT_INDICES(IJK,NEVS(0,1),MAXS(1),MAXS(2),MAXS(3),
     &  ERROR,*9999)

C     Loop over all the lattice grid points in the current
C     element and determine which ones to use as global grid
C     points, considering adjacent elements and collapsing
C     information.

      DO k=1,NQXI(3,nqsc)
        DO j=1,NQXI(2,nqsc)
          DO i=1,NQXI(1,nqsc)
            !calculate the lattice point index
            llat=((k-1)*NQXI(2,nqsc)+j-1)*NQXI(1,nqsc)+i
            LPI=NLATNE(ne)-1+llat
            adjacent_element=ne
            !isolated element limits
            IJKMIN(1)=INT_LATT_JUMP(1,j,k,MINS,MAXS,IJK(1,1),1)
            IJKMIN(2)=INT_LATT_JUMP(i,1,k,MINS,MAXS,IJK(1,2),1)
            IJKMIN(3)=INT_LATT_JUMP(i,j,1,MINS,MAXS,IJK(1,3),1)
            IJKMAX(1)=INT_LATT_JUMP(NQXI(1,nqsc),j,k,MINS,MAXS,IJK(1,1),
     &        1)
            IJKMAX(2)=INT_LATT_JUMP(i,NQXI(2,nqsc),k,MINS,MAXS,IJK(1,2),
     &        1)
            IJKMAX(3)=INT_LATT_JUMP(i,j,NQXI(3,nqsc),MINS,MAXS,IJK(1,3),
     &        1)
            !loop over the adjacent elements
            DO ne_adj=2,NEVS(0,1)
              !determine shared element number
              shar_ne=INT_LATT_JUMP(i,j,k,IJKMIN,IJKMAX,NEES(1,ne_adj),
     &          0)
              shar_ne=ne+(NEVS(1,ne_adj)-ne)*
     &          INT(shar_ne/NEVS(1,ne_adj))*INT(NEVS(1,ne_adj)/shar_ne)
              adjacent_element=min(shar_ne,adjacent_element)
              !determine range limitations
              DO ni=1,3 !loop over the three xi directions for the
                        !adjacent element
                IF(NEVS(9+ni,ne_adj).NE.0.
     &            .AND.NEVS(9+ni,ne_adj).NE.ne) THEN
                  ADJN(ni)=min(ADJN(ni),
     &              NQXI(abs(NEVS(9+ni,ne_adj)),NQS(shar_ne)))
                ELSE
                  ADJN(ni)=min(ADJN(ni),NQXI(ni,nqsc))
                ENDIF
              ENDDO
            ENDDO !end adjacent loop
            DO ni=1,3
              IJKMAX(ni)=min(IJKMAX(ni),(IJKMIN(ni)+ADJN(ni)-1))
              ADJN(ni)=NQXI(ni,nqsc)
            ENDDO
            DO ne_adj=2,NEVS(0,1) !second loop over adjacent elements
              !determine shared element number
              shar_ne=INT_LATT_JUMP(i,j,k,IJKMIN,IJKMAX,NEES(1,ne_adj),
     &          0)
              shar_ne=ne+(NEVS(1,ne_adj)-ne)*
     &          INT(shar_ne/NEVS(1,ne_adj))*INT(NEVS(1,ne_adj)/shar_ne)
              adjacent_element=min(shar_ne,adjacent_element)
              !determine range limitations
              DO ni=1,3
                IF(NEVS(9+ni,ne_adj).NE.0.AND.(shar_ne.NE.ne)) THEN
                  ADJN(ni)=min(ADJN(ni),
     &              NQXI(abs(NEVS(9+ni,ne_adj)),NQS(shar_ne)))
                ELSE
                  ADJN(ni)=min(ADJN(ni),NQXI(ni,nqsc))
                ENDIF
              ENDDO
            ENDDO !end second adjacent loop
            ne_adj=adjacent_element
C           ne_adj is a global element number. NEVS needs to be
C           searched to find the column which holds ne_adj.
            DO counter=1,NEVS(0,1)
              IF(NEVS(1,counter).EQ.ne_adj) THEN
                line=counter
              ENDIF
            ENDDO
            IF(adjacent_element.LT.ne) THEN
              DO ni_adj=1,3
                IF(abs(NEVS(10,line)).EQ.ni_adj.OR.
     &            abs(NEVS(11,line)).EQ.ni_adj
     &            .OR.abs(NEVS(12,line)).EQ.ni_adj) THEN
                  !direction ni_adj in adjacent element is shared with
                  !direction ni in current element and range is set to
                  !be minimum of grid points in shared directions
                  !and grid points from other higher numbered elements
                  min_grid=NLATXIM
                  IF(abs(NEVS(10,line)).EQ.ni_adj) THEN
                    ni=1
                    min_grid=NQXI(ni,nqsc)
                  ENDIF
                  IF(abs(NEVS(11,line)).EQ.ni_adj) THEN
                    ni=2
                    min_grid=min(min_grid,NQXI(ni,nqsc))
                  ENDIF
                  IF(abs(NEVS(12,line)).EQ.ni_adj) THEN
                    ni=3
                    min_grid=min(min_grid,NQXI(ni,nqsc))
                  ENDIF
                  ADJN(ni_adj)=min(NQXI(ni_adj,NQS(adjacent_element)),
     &              min_grid,ADJN(ni))
                ELSE
                  !direction ni_adj in adjacent element is not shared
                  !with current element
                  ADJN(ni_adj)=NQXI(ni_adj,NQS(adjacent_element))
                ENDIF
              ENDDO
            ENDIF
            !determine the ranges for calculating the vertex lattice
            !indices
            IF(adjacent_element.EQ.ne) THEN
              DO counter=1,3
                IF(ADJN(counter).EQ.1) THEN
                  LIM(counter) = 1
                ELSE
                  LIM(counter) = NQXI(counter,nqsc)
                ENDIF
              ENDDO
            ELSE
              DO counter=1,3
                LIM(counter)=ADJN(counter)
              ENDDO
            ENDIF

            CALL VERT_LATT_INDICES(IJK_NEW,NEVS(0,line),LIM(1),LIM(2),
     &        LIM(3),ERROR,*9999)
            !calculate the ijk grid pt indices and xi location
            DO ni=1,3
              LOWLIM(ni)=1+INT((1+NQXI(ni,nqsc))/2)-INT((1+ADJN(ni))/2)
              HIGHLIM(ni)=LOWLIM(ni)+ADJN(ni)-1
            ENDDO
            GPI=INT_LATT_JUMP(i,j,k,LOWLIM,HIGHLIM,IJK_NEW(1,1),1)
            GPJ=INT_LATT_JUMP(i,j,k,LOWLIM,HIGHLIM,IJK_NEW(1,2),1)
            GPK=INT_LATT_JUMP(i,j,k,LOWLIM,HIGHLIM,IJK_NEW(1,3),1)
            nlat=(NLATNE(adjacent_element)-1)+((GPK-1)*
     &        NQXI(2,NQS(adjacent_element))+GPJ-1)*
     &        NQXI(1,NQS(adjacent_element))+GPI
C           Store the new global grid point or update lattice to
C           global grid point map.
            IF(NQNLAT(nlat).EQ.0) THEN !new global grid point
              NQT=NQT+1
              CALL ASSERT(NQT.LE.NQM, '>> Increase NQM',ERROR,*9999)
              NQNLAT(nlat)=NQT
              NQNLAT(LPI)=NQT
              NQNE(ne,llat)=NQT
              NLATPNQ(NQT)=nlat
              NLATNQ(nlat)=0
              NENQ(1,NQT)=ne
            ELSE !global grid point exists.
              NQNLAT(LPI)=NQNLAT(nlat)
              NQNE(ne,llat)=NQNLAT(nlat)
              loc=NLATPNQ(NQNLAT(LPI))
              DO WHILE (NLATNQ(loc).NE.0)
                loc=NLATNQ(loc)
              ENDDO
              NLATNQ(loc)=LPI
              NLATNQ(LPI)=0
            ENDIF
          ENDDO !end i loop
        ENDDO !end j loop
      ENDDO ! end k loop
      
      CALL EXITS('CALC_LATTICE_MAP')
      RETURN
 9999 CALL ERRORS('CALC_LATTICE_MAP',ERROR)
      CALL EXITS('CALC_LATTICE_MAP')
      RETURN 1      
      END

      
