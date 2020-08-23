      SUBROUTINE MESH_DEFORM(nb,NBJ,NEELEM,NELIST,NORD,NPNE,nr,NVJE,nx,
     '  NXI,NYNP,BBM,CE,XAB,XP,YP,FIX,GAS,ERROR,*)

C#### Subroutine: MESH_DEFORM
C###  Description:
C###    MESH_DEFORM calculates volume deformations in 1D mesh.
C###    AJS update: Calculates new BBM concentration using mass balance
C###    (mass inspired/expired and mass exchanged with blood).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M),NELIST(0:NEM),
     &  NORD(5,NE_R_M),NPNE(NNM,NBFM,NEM),nr,NVJE(NNM,NBFM,NJM,NEM),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM),YP(NYM)
      LOGICAL FIX(NYM,NIYFIXM),GAS
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb_radius,ne,ne0,nn,noelem,np1,npnn(2),nv,nv2,nvnn(2),ny1,
     &  kount,np2,nj
      REAL*8 DV,length,radius,volume_old,total_old,dist,BBM_mass,
     &  mass_exch,mass_in,new_mass,new_volume,length0
      REAL*8 LENGTH_1D,RADIUS_1D

      CALL ENTERS('MESH_DEFORM',*9999)

      DV=INLET_FLOW(nx)*DT
      total_old=current_volume(nx)
      NELIST(0)=0
      kount=0

      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        length=LENGTH_1D(NBJ,ne,NPNE,NVJE,XP)
c        radius=RADIUS_1D(NBJ,ne,NPNE,NVJE,XP)
        nb=NBJ(nj_radius,ne)
        np2=NPNE(2,nb,ne)
        radius=XP(1,1,nj_radius,np2)
        volume_old=length*PI*radius**2
        XAB(9,ne)=volume_old+XAB(1,ne)*DV
        IF(NXI(1,0,ne).EQ.0)THEN
          nb=NBJ(nj_flow,ne)
          np1=NPNE(2,nb,ne) !terminal node
          dist=XP(1,1,nj_flow,np1)*DV/(PI*radius**2)
        ENDIF
       
c        length0=0.d0
c        np1=NPNE(1,1,ne)
c        np2=NPNE(2,1,ne)
c        DO nj=1,3
c          length0=length0+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2
c        ENDDO
c        length0=DSQRT(length0)
c        np1=NPNE(2,1,ne)
c        write(*,*) length,length0
c        pause
        
        XAB(2,ne)=volume_old
        IF(DABS(XAB(1,ne)).GT.ZERO_TOL)THEN
          NELIST(0)=NELIST(0)+1
c          CALL ASSERT(NELIST(0).LE.NEM,'>> Increase NEM ',ERROR,*9999)
          NELIST(NELIST(0))=ne
          BBM(1,ne) = 0.d0
          BBM(2,ne) = 0.d0
        ENDIF

        nb=NBJ(1,ne)
        nb_radius=NBJ(nj_radius,ne)
        DO nn=1,NNT(nb)
          npnn(nn)=NPNE(nn,nb,ne) !node # at start of current branch
          nvnn(nn)=NVJE(nn,nb_radius,nj_radius,ne) !versions for geometry
        ENDDO !nn

        DO nn=1,2 !adjust the radius
          XP(1,nvnn(nn),nj_radius,npnn(nn))=DSQRT(XAB(9,ne)/(PI*length))
        ENDDO !nn
      ENDDO !noelem

      O2_EXCHANGED=0.0d0 !change in mass due to gas exchange
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(nj_flow,ne)
        nv=NVJE(2,nb,nj_flow,ne)
        IF(NTB.GT.0.AND.NXI(1,0,ne).EQ.0)THEN
          np1=NPNE(2,nb,ne) !terminal node
          IF(XP(1,1,nj_alveoli,np1).EQ.1.d0)THEN !can have LPM attached
C... This condition is used to cope with having a combination of LPMs
C... and multi-branching acini. Without this condition there would
C... effectively be LPMs on the ends of the terminal alveolar ducts.
C... Not infallible: LPM can be attached if NO alveoli.
            kount=kount+1 !increment LPM counter
            ny1=NYNP(1,1,1,np1)
C AJS - mass lost to blood
            BBM_mass=BBM(2,ne)*BBM(1,ne)
            new_volume=BBM(1,ne)+XP(1,nv,nj_flow,np1)*DV !volume
            IF(ITYP7(nr,nx).GT.1)THEN
              nv2=NVJE(2,nb,nj_pob,ne)
              mass_exch=XP(1,nv2,nj_source,np1)*DT
              O2_EXCHANGED=O2_EXCHANGED+mass_exch
            ELSE
              mass_exch=0.0d0
            ENDIF
C.........Mass into / out of BBM
            IF(DV.GE.0.d0) mass_in=YP(ny1)*XP(1,nv,nj_flow,np1)*DV !use T.B concentration
            IF(DV.LT.0.d0) mass_in=BBM(2,ne)*XP(1,nv,nj_flow,np1)*DV !use BBM concentration (flow -ve)
            new_mass=BBM_mass+mass_in+mass_exch
c..........Update BBM volume 
            BBM(1,ne)=new_volume !BBM(1,ne)+XP(1,nv,nj_flow,np1)*DV !volume
            XAB(9,ne)=XAB(9,ne)+BBM(1,ne) !the volume
C.........If inspiration set the BBM concentration
            IF(DV.GT.0.d0)THEN
              IF(new_volume.GT.0.0d0) BBM(2,ne)=new_mass/new_volume !(BBM_mass+mass_in+mass_exch)/new_volume
C.........If expiration set the terminal node concentration (for gas exchange also update BBM conc)
            ELSE
              IF(ITYP7(nr,nx).GT.1.AND.new_volume.GT.0.0d0) 
     &          BBM(2,ne)=new_mass/new_volume !(BBM_mass+mass_in+mass_exch)/new_volume
              IF(FIX(ny1,1)) YP(ny1)=BBM(2,ne) ! ??? this command is repeated 
              YP(ny1)=BBM(2,ne)
            ENDIF

! C.........If inspiration set the BBM concentration
!             IF(DV.GE.0.d0)THEN
!               !conc = (conc(0)*vol(0)+concTB*volflowTB)/(vol(0)+volflowTB)
! !                 BBM(2,ne)=(BBM(2,ne)*BBM(1,ne)+YP(ny1)*XP(1,nv,nj_flow,
! !      &            np1)*DV)/(BBM(1,ne)+XP(1,nv,nj_flow,np1)*DV)
!               mass_in=YP(ny1)*XP(1,nv,nj_flow,np1)*DV
!               BBM(2,ne)=(BBM_mass+mass_in+mass_exch)/new_volume
! C.........If expiration set the terminal node concentration
!             ELSE
!               mass_in=BBM(2,ne)*XP(1,nv,nj_flow,np1)*DV ! negative b/c expired
! c              IF(GAS_EXCHANGE_TERM)
! c     &          BBM(2,ne)=(BBM_mass+mass_in+mass_exch)/new_volume
!               YP(ny1)=BBM(2,ne)
!             ENDIF
!             BBM(1,ne)=BBM(1,ne)+XP(1,nv,nj_flow,np1)*DV !volume
!             XAB(9,ne)=XAB(9,ne)+BBM(1,ne) !the volume
            IF(kount.EQ.1)THEN
!               WRITE(OP_STRING,'('' Mass: current= '',D10.3,
!      &          '' in/out='',D10.3,'' exch='',D10.3,
!      &          '' Conc: term= '',D10.3,'' BBM='',D10.3)')
!      &          BBM_mass,mass_in,mass_exch,YP(ny1),BBM(2,ne)
!               CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDIF
      ENDDO

C Update the volumes below and path length from stem
      DO noelem=NEELEM(0),1,-1
        ne=NEELEM(noelem)
        ne0=NXI(-1,1,ne) !parent
        IF(ne0.NE.0)THEN
            
          IF(NORD(5,ne).EQ.1)THEN !symmetry, add twice as much
            XAB(9,ne0)=XAB(9,ne0)+2.d0*XAB(9,ne) !sum volumes below
c            XAB(1,ne0)=XAB(1,ne0)+2.d0*XAB(1,ne) !sum volumes below
            XAB(2,ne0)=XAB(2,ne0)+2.d0*XAB(2,ne) !sum volumes below
          ELSE
            XAB(9,ne0)=XAB(9,ne0)+XAB(9,ne)
c            XAB(1,ne0)=XAB(1,ne0)+XAB(1,ne)
            XAB(2,ne0)=XAB(2,ne0)+XAB(2,ne)
          ENDIF
        ENDIF !ne0.NE.0
      ENDDO !noelem

      CURRENT_VOLUME(nx)=XAB(9,NEELEM(1))
      DV_TOTAL(nx)=DV_TOTAL(nx)+DV

      CALL EXITS('MESH_DEFORM')
      RETURN

 9999 CALL ERRORS('MESH_DEFORM',ERROR)
      CALL EXITS('MESH_DEFORM')
      RETURN 1
      END



