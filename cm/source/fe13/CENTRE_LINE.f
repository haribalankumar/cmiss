      SUBROUTINE CENTRE_LINE(AXIS_XI_DIVISIONS,AXIS_XI,CIRCUM_XI,IBT,      
     '  IDO,INP,NBJ,NDDATA,NDP,NEELEM,NKJE,NPF,NPNE,nr,NVJE,NXI,SE,
     '  START_ELEMENT,WD,WG,XA,XE,XIG,XP,ZD,ERROR,*)
      
C#### Subroutine: CENTRE_LINE
C###  Description:
C###   Calculates the centre line for a tubular mesh
C###  Requires NXI to be set up (use NENXI if not done)
C*** Created by Peter Bier, Mar 2003.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'call00.cmn'

!     Parameter list
      INTEGER AXIS_XI,AXIS_XI_DIVISIONS,CIRCUM_XI,IBT(3,NIM,NBFM),      
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),NBJ(NJM,NEM),
     '  NDDATA(0:NDM,0:NRM),NDP(NDM),NEELEM(0:NE_R_M,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),nr,
     '  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM),START_ELEMENT      
      REAL*8 SE(NSM,NBFM,NEM),WD(NJM,NDM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER AXIS_XI_INDEX,CIRCUM_ELEM(0:NEM),com,ne,nb,nbgauss,nd,ng,      
     '  ni,nj      
      REAL*8 AXIS_XI_INC, ARC_COM(3),COM_SUM(3),X(3),XI(3)
      LOGICAL AXIS_DONE,FOUND_GAUSS_BASIS
!     Function
      REAL*8 PXI
      CALL ENTERS('CENTRE_LINE',*9999)
      
C     Find a 1D linear lagrange basis fn to use for gaussian quadrature
      FOUND_GAUSS_BASIS=.FALSE.
      DO nb=1,NBT
        IF( (NIT(nb).EQ.1) .AND. (IBT(1,1,nb).EQ.1) .AND.
     '    IBT(2,1,nb).EQ.1) THEN
          nbgauss=nb
          FOUND_GAUSS_BASIS=.TRUE.
        ENDIF !1D linear lagrange
      ENDDO !nb

      CALL ASSERT( FOUND_GAUSS_BASIS,'>>need to define 1D ' //
     '  'linear lagrange basis function', ERROR,*9999)
      
C     existing data will be clobbered
      NDT = 0
      NDDATA(0,0) = NDT
      NDDATA(0,nr) = 0
      nd = NDT
      
C     each xi value corresponds to a com point
      AXIS_XI_INC = 1.0d0/DBLE(AXIS_XI_DIVISIONS)
      
      DO nj=1,3
        XI(nj) = 0.d0
      ENDDO
      
C     Beginning with the start element, increment around the circumference
C     building up a list of circumferential points, then increment the
C     xi coord and repeat for the entire tube axis
      
      com = 0 ! could use nd as data increment variable
C     com is used in case later want to avoid clobbering existing data
      
      AXIS_XI_INDEX = 0
      
      ne=NEELEM(START_ELEMENT,nr)
      CALL GET_CIRCUM_ELEM(CIRCUM_XI,ne,NXI,CIRCUM_ELEM,ERROR,*9999)
      
      AXIS_DONE=.FALSE.
      
      DO WHILE( .NOT.AXIS_DONE ) ! not processed entire length of axis

C       Set variables for a pass around the circumference
        DO nj=1,3
          COM_SUM(nj)=0.d0
        ENDDO
        
        ni = 1
                
        DO WHILE(ni .LE. CIRCUM_ELEM(0)) ! more elements in circumference
C         For each elem calculate com
          ne = CIRCUM_ELEM(ni) ! 
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,          
     '      NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)                           
          nb = NBJ(1,ne)
          
          DO nj=1,3
            ARC_COM(nj) = 0.d0 !
          ENDDO
          
C         Use gaussian quadrature to determine c.o.m. for line
C         This will give an exact value for cubic elements
          
          DO ng=1,NGT(nbgauss) ! for both gauss points
            XI(CIRCUM_XI) = XIG(1,ng,nbgauss)
            XI(AXIS_XI) = DBLE(AXIS_XI_INDEX) * AXIS_XI_INC
            DO nj=1,3
              X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          XI,XE(1,nj)) !global coordinates of surface point
              
C             Do not need to divide by arc length to calculate com as
C             arc length is 1
              ARC_COM(nj) = ARC_COM(nj) + WG(ng,nbgauss) * X(nj)
            ENDDO                
          ENDDO ! for both gauss points
          
C         Add contribution from arc com to total com for circumference
          DO nj=1,3
            COM_SUM(nj) = COM_SUM(nj) + ARC_COM(nj) 
          ENDDO
          
          ni = ni + 1 !increment to the next element
        ENDDO
        
C       Add new com data point
        com = com + 1
        nd = nd + 1
        DO nj=1,3            
          ZD(nj,nd) = COM_SUM(nj) / DBLE(CIRCUM_ELEM(0))
          WD(nj,nd) = 1.0d0
        ENDDO !nj
        
        NDDATA(com,nr) = nd
        NDDATA(0,nr) = NDDATA(0,nr) + 1
        NDDATA(0,0) = NDDATA(0,0) + 1
        NDT = NDDATA(0,0)
        NDP(nd)=com

        AXIS_XI_INDEX = AXIS_XI_INDEX + 1
        
C       Only go to the next element AFTER having done XI=1
        IF( AXIS_XI_INDEX .GT. AXIS_XI_DIVISIONS) THEN ! increment to next element            
          AXIS_XI_INDEX = 0
          ne = NXI(AXIS_XI,1,ne) ! fetch next element in axis_xi dir
          IF( ne .EQ. 0 ) THEN ! no more elements
            AXIS_DONE = .TRUE.
          ELSE
C           more elements so calculate next circumference
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '        nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,
     '        ERROR,*9999)               
            CALL GET_CIRCUM_ELEM(CIRCUM_XI,ne,NXI,
     '        CIRCUM_ELEM,ERROR,*9999)      
            AXIS_XI_INDEX = 1
          ENDIF ! no more elements
        ENDIF ! move to next element

      ENDDO                     ! not processed entire axis

      CALL EXITS('CENTRE_LINE')
      RETURN
 9999 CALL ERRORS('CENTRE_LINE',ERROR)
      CALL EXITS('CENTRE_LINE')
      RETURN 1
      END


