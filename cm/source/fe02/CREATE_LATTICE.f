      SUBROUTINE CREATE_LATTICE(NBJ,NEELEM,NENP,NENQ,NLATNE,NLATNQ,
     &  NLATPNQ,NNB,NPNE,NQNE,NQNLAT,NQS,NQSCNB,NQXI,NRLIST,NWQ,XIQ,
     &  ERROR,*)

C#### Subroutine: CREATE_LATTICE
C###  Description:
C###    CREATE_LATTICE is the parent subroutine for the
C###    set up of the lattice grid scheme. Create lattice steps
C###    through regions and elements.      

C#### Variable: NEAM
C###  Type: INTEGER
C###  Set_up: CREATE_LATTICE
C###  Description:
C###    NEAM is the maximum number of elements adjacent to
C###    element ne. It is used to dimension the element
C###    sharing arrays NEVS and NEES for the lattice grid
C###    method. NEAM takes it value from the sum of adjacent
C###    elements for all the nodes of ne. This value is
C###    always at least as big as needed, but may be larger
C###    than needed and a more suitable value may be able
C###    be found.   

C#### Variable: NLATNE(ne)
C###  Type: INTEGER      
C###  Set_up: CREATE_LATTICE
C###  Description:
C###    NLATNE(ne) is the mapping between element ne
C###    and the first lattice grid point for that element.      
C###    A zero entry in NLATNE indicated that less than NEM
C###    elements have been defined, and there is no element
C###    for this element number.   

      IMPLICIT NONE
      
      INCLUDE 'cbdi02.cmn'      
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     &  NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),NLATNE(NEQM+1),
     &  NLATNQ(NEQM*NQEM),NLATPNQ(NQM),NNB(4,4,4,NBFM),
     &  NPNE(NNM,NBFM,NEM),NQNE(NEQM,NQEM),NQNLAT(NEQM*NQEM),NQS(NEQM),
     &  NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NWQ(8,0:NQM,NAM)
      REAL*8 XIQ(NIM,NQM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER i,nb,ne,NEAM,nlatnet,nn,noelem,nonrlist,nq,nqsc,nqsc_old,
     &  nr

      CALL ENTERS('CREATE_LATTICE',*9999)
      ! Initialise the mapping arrays

      DO i=1,NEQM*NQEM
        NLATNQ(i)=0
        NQNLAT(i)=0
      ENDDO
      
      DO nonrlist=1,NRLIST(0) !loop over regions
        nr=NRLIST(nonrlist)
        NQR(1,nr)=NQT+1
        !set up NLATNE
        DO ne=1,NEQM+1 !initialise NLATNE.
          NLATNE(ne)=0
        ENDDO
        ne=NEELEM(1,nr)
        NLATNE(ne)=1
        nqsc=NQS(ne)
        nlatnet=NQXI(1,nqsc)*NQXI(2,nqsc)*NQXI(3,nqsc)        
        DO noelem=2,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nqsc_old=nqsc
          nqsc=NQS(ne)
          NLATNE(ne)=NLATNE(NEELEM(noelem-1,nr))+
     &      NQXI(1,nqsc_old)*NQXI(2,nqsc_old)*NQXI(3,nqsc_old)
          nlatnet=nlatnet+NQXI(1,nqsc)*NQXI(2,nqsc)*NQXI(3,nqsc)
        ENDDO  
        NLATNE(ne+1)=nlatnet+1
        DO noelem=1,NEELEM(0,nr) !loop over elements
          ne=NEELEM(noelem,nr)
          !calculate NEAM for ne.
          NEAM=0
          nb=NBJ(1,ne)
          DO nn=1,NNT(nb)
            NEAM=NEAM+NENP(NPNE(nn,nb,ne),0,nr)
          ENDDO
          CALL SIZE_LATT_ARRAYS(NBJ,ne,NEAM,NENP,NENQ,NNB,NLATNE,NLATNQ,
     &      NLATPNQ,NPNE,NQNE,NQNLAT,NQS,NQXI,nr,ERROR,*9999)
        ENDDO !end loop elements
        NQR(2,nr)=NQT
      ENDDO !end loop regions

      DO nq=1,NQT
        CALL LATTICE_NWQ(NLATNE,NLATNQ,NLATPNQ,NQNLAT,NQS,NQXI,NWQ,nq,
     &    ERROR,*9999)
      ENDDO
      CALL CALC_LATTICE_XIQ(NLATNE,NLATNQ,NLATPNQ,NQNLAT,NQS,NQSCNB,
     &  NQXI,XIQ,ERROR,*9999)
      
      CALL EXITS('CREATE_LATTICE')
      RETURN
 9999 CALL ERRORS('CREATE_LATTICE',ERROR)
      CALL EXITS('CREATE_LATTICE')
      RETURN 1
      END
