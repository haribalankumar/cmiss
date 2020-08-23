      SUBROUTINE UPMATE(IBT,IDO,INP,IRCQS_SPATIAL,NBJ,
     &  NEELEM,NELIST,NENQ,NKJ,NPLIST,NPNODE,NRLIST,NKH,NKJE,NPF,NPNE,
     &  NQET,NQNE,NQS,NVHP,NVJE,NVJP,NXI,NXLIST,NYNP,CE,CELL_RCQS_VALUE,
     &  CP,CQ,RCQS_SPATIAL,PAOPTI,SE,XA,XE,XP,XQ,YP,STRING,ERROR,*)

C#### Subroutine: UPMATE
C###  Description:
C###    UPMATE updates the material parameters of the bidomain problem
C###    from the monodomain parameters

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  IRCQS_SPATIAL(0:NQRSVM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NENQ(0:8,NQM),NKJ(NJM,NPM),NPLIST(0:NPM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),NQET(NQSCM),
     '  NQNE(NEQM,NQEM),NQS(NEQM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NXLIST(0:NXM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CELL_RCQS_VALUE(NQRM,NQVM),
     '  CP(NMM,NPM,NXM),CQ(NMM,NQM,NXM),
     '  RCQS_SPATIAL(NQRSVM,NQM),PAOPTI(NOPM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,opt_var,param,PART2,
     '  nb,ne,N3CO,
     '  nj,nj1,nm,noelem,nonode,no_nrlist,
     '  np,nq,nqele,nqq,nqrs,nqrsv,nr,NTLIST,
     '  nx,nxc,nxc2,nx_s,nx_t,
     '  XI_COUNT,variant
      REAL*8 border,CENTRE(3),DENSITY,ELLIPSE(3),FRACTION,MATERIAL(2),
     '  radius1,radius2,Ro,RFROMC,SEMIAXIS(3),theta,Xi,XI1(3),PXI
      CHARACTER OPERATION*16,TYPE*16,UPDATE*16
      LOGICAL ALL_REGIONS,CBBREV,CELL,SPATIALLY_VARYING
C     INTEGER ILT_TOT1,ILT_TOT2,nh1,nh2,nhc1,nhc2,niy2,njc2,nktot1,
C    '  nm1,nm2,nmc1,nmc2,nx1,nx2,nxc1,nx_opt,ny2,nj2
C     REAL*8 CONST,SCALE
C     CHARACTER ATPOINTS*16,WITH*16

      CALL ENTERS('UPMATE',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

        UPDATE='HELPMATERIAL'
        CALL UPFG(UPDATE,%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     '    %VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     '    %VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     &    STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update material bidomain
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <from_class #[1]>
C###    Specify the class of the monodomain problem
C###  Parameter:      <to_class #[2]>
C###    Specify the class of the bidomain problem
C###  Description: Updates the material parameters of the
C###    bidomain problem from the monodomain parameters
C###

        OP_STRING(1)=STRING(1:IEND)//' bidomain'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<from_class #[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<to_class #[2]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C PM 26-JUL-01 : command added
C#### Command: FEM update material field
C###  Parameter:   <material_param #[1]>
C###    Specify the material parameter number
C###  Parameter:   <no_field #[1]>
C###    Specify the field number in ipfiel

        OP_STRING(1)=STRING(1:IEND)//' field'
        OP_STRING(2)=BLANK(1:IEND)//'<material_param #[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<no_field #[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(5)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update material infarct
C###  Parameter:   material_param #[1]
C###    Specify the material parameter number
C###  Parameter:   centre xc,yc
C###    Specify the ellipse centre
C###  Parameter:   semiaxes xd,yd
C###    Specify the ellipse semiaxes
C###  Parameter:   values p1,p2
C###    Specify the material parameter values
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Description: Updates grid material parameters for an elliptical
C###    infarct with centre (xc,yc), major and minor semiaxes xd,yd
C###    and material params varying linearly from p1 at the centre to
C###    p2 at the periphery.

        OP_STRING(1)=STRING(1:IEND)//' infarct'
        OP_STRING(2)=BLANK(1:IEND)//'material_param #[1]'
        OP_STRING(3)=BLANK(1:IEND)//'centre xc,yc'
        OP_STRING(4)=BLANK(1:IEND)//'semiaxes xd,yd'
        OP_STRING(5)=BLANK(1:IEND)//'values p1,p2'
        OP_STRING(6)=BLANK(1:IEND)//'<border FRACTION[0.1]>'
        OP_STRING(7)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C LKC 11-AUG-2001 new section
C
C#### Command: FEM update material nodes
C###  Parameter:      <material_param #[1]>
C###    Specify the material parameter number
C###  Parameter:      <from_class #[2]>
C###    Specify the class of the material to update from
C###  Parameter:      <to_class #[1]>
C###    Specify the class of the material to update to
C###  Description:
C###    Update the material parameters from one problem class to
C###    another.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
        OP_STRING(1)=STRING(1:IEND)//' nodes'
        OP_STRING(2)=BLANK(1:IEND)//'material_param #[1]'
        OP_STRING(3)=BLANK(1:IEND)//'<from_class #[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<to_class #[2]>'
        OP_STRING(5)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM update material density
C###  Parameter:      <(value)[1]>
C###    Specify the density.
C###  Parameter:      <element (#s/all)[all]>
C###    Specify the elements in which to update the density.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
        OP_STRING(1)=STRING(1:IEND)//' density'
        OP_STRING(2)=BLANK(1:IEND)//'<(value)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPMATE',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'SUBSTITUTE',4,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='SUBSTITUTE'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'ADD',2,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='ADD'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'SUBTRACT',4,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='SUBTRACT'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'MULTIPLY',2,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='MULTIPLY'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'DIVIDE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='DIVIDE'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'BIDOMAIN',1,noco+1,NTCO,N3CO)) THEN
          TYPE='BIDOMAIN'
C PM 26-JUL_01: Added field option
        ELSEIF(CBBREV(CO,'FIELD',1,noco+1,NTCO,N3CO)) THEN
          TYPE='FIELD'
        ELSEIF(CBBREV(CO,'INFARCT',1,noco+1,NTCO,N3CO)) THEN
          TYPE='INFARCT'
        ELSEIF(CBBREV(CO,'NODES',1,noco+1,NTCO,N3CO)) THEN
          TYPE='NODES'
        ELSEIF(CBBREV(CO,'OPTIMISE',2,noco+1,NTCO,n3co)) THEN
          TYPE='OPTIMISE'
        ELSEIF(CBBREV(CO,'DENSITY',2,noco+1,NTCO,n3co)) THEN
          TYPE='DENSITY'
          DENSITY=RFROMC(CO(N3CO+1))
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(TYPE(1:9).EQ.'OPERATION') THEN
          ! DO NOTHING: see subroutine UPFG
        ELSEIF(TYPE(1:8).EQ.'BIDOMAIN') THEN
          IF(CBBREV(CO,'FROM_CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_s,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_s.NE.0,'Invalid source class',ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL ASSERT(ITYP5(nr,nx_s).EQ.2.AND.ITYP19(nr,nx_s)
     '        .EQ.1.AND.ITYP2(nr,nx_s).EQ.9,
     '        'Source equation must be activation model',ERROR,*9999)
          ENDDO

          IF(CBBREV(CO,'TO_CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=2
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_t,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_t.NE.0,'Invalid target class',ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL ASSERT(ITYP5(nr,nx_t).EQ.1.AND.ITYP2(nr,nx_t).EQ.5,
     '        'Target equation must be div(k.grad(u))=f',ERROR,*9999)
          ENDDO

          CALL ASSERT(nx_s.NE.nx_t,
     '      'Source and target classes must be different',ERROR,*9999)

C PM 26-JUL-01 field option
        ELSE IF(TYPE(1:5).EQ.'FIELD') THEN

          IF(CBBREV(CO,'CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF

          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.NE.0,'Invalid class',ERROR,*9999)

          IF(CBBREV(CO,'MATERIAL_PARAM',1,noco+1,NTCO,N3CO)) THEN
            nm=IFROMC(CO(N3CO+1))
          ELSE
            nm=1
          ENDIF
          IF(CBBREV(CO,'NO_FIELD',1,noco+1,NTCO,N3CO)) THEN
            nj=IFROMC(CO(N3CO+1))
          ELSE
            nj=1
          ENDIF

        ELSE IF(TYPE(1:7).EQ.'INFARCT') THEN
          IF(CBBREV(CO,'CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF

          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.NE.0,'Invalid class',ERROR,*9999)

          IF(CBBREV(CO,'MATERIAL_PARAM',1,noco+1,NTCO,N3CO)) THEN
            nm=IFROMC(CO(N3CO+1))
          ELSE
            ERROR='no material parameter specified'
            GOTO 9999
          ENDIF
          IF(CBBREV(CO,'CENTRE',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),2,NTLIST,CENTRE,ERROR,*9999)
          ELSE
            ERROR='no ellipse centre specified'
            GOTO 9999
          ENDIF
          IF(CBBREV(CO,'SEMIAXES',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),2,NTLIST,SEMIAXIS,ERROR,*9999)
          ELSE
            ERROR='no ellipse semiaxes specified'
            GOTO 9999
          ENDIF
          IF(CBBREV(CO,'VALUES',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),2,NTLIST,MATERIAL,ERROR,*9999)
          ELSE
            ERROR='no material values specified'
            GOTO 9999
          ENDIF
          IF(CBBREV(CO,'BORDER',1,noco+1,NTCO,N3CO)) THEN
            FRACTION=RFROMC(CO(N3CO+1))
          ELSE
            FRACTION=0.1d0
          ENDIF

C LKC 11-AUG-2001 new section for nodes
        ELSEIF(TYPE(1:7).EQ.'NODES') THEN

          IF(CBBREV(CO,'FROM_CLASS',2,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=2
          ENDIF
          IF(CBBREV(CO,'TO_CLASS',2,noco+1,NTCO,N3CO)) THEN
            nxc2=IFROMC(CO(N3CO+1))
          ELSE
            nxc2=1
          ENDIF

          IF(CBBREV(CO,'MATERIAL_PARAM',3,noco+1,NTCO,N3CO)) THEN
            nm=IFROMC(CO(N3CO+1))
          ELSE
            nm=1
          ENDIF

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              CP(nm,np,nxc2)=CP(nm,np,nxc)
            ENDDO
          ENDDO

        ELSE IF(TYPE(1:8).EQ.'OPTIMISE') THEN
C          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
C          CALL ASSERT(nx.GT.0,
C     '      '>>No nx defined for this optimisation class',
C     '      ERROR,*9999)
           nx=1

         ELSE IF(TYPE(1:7).EQ.'DENSITY') THEN
           CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '       ERROR,*9999) !gives an element list of parents
           CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
           nxc=NXLIST(1)
           DO noelem=1,NELIST(0)
             ne=NELIST(noelem)
             CE(IL_DENSITY,ne,nxc)=DENSITY
           ENDDO
        ELSE
          ERROR='Code needs updating'
          GOTO 9999
        ENDIF !type


        IF(TYPE(1:9).EQ.'OPERATION')THEN

          UPDATE='MATERIAL'
          CALL UPFG(UPDATE,CP,%VAL(0),%VAL(0),NKH,NKJ,NPLIST,
     '      NPNODE,NRLIST,NVHP,NVJP,NYNP,OPERATION,PART2,%VAL(0),XP,
     &      %VAL(0),YP,%VAL(0),STRING,ERROR,*9999)

C PM 26-JUL-01: Field option
        ELSEIF(TYPE(1:5).EQ.'FIELD') THEN
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL NX_LOC(NX_INQUIRE,1,nx_s,NX_SOLVE,ERROR,*9999)
            DO nq=NQR(1,nr),NQR(2,nr)
              ne=NENQ(1,nq)
              nb=NBJ(1,ne)
              nj1=NJ_LOC(NJL_FIEL,nj,nr)
!             nj2=NJ_LOC(NJL_FIEL,nj+1,nr)
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,
     '          *9999)
              XI1(1)=0.0d0
              nqele=NQNE(ne,1)
              XI_COUNT=1
              DO WHILE ((nqele.NE.nq).AND.
     '          (nqele.NE.NQNE(ne,NQET(NQS(ne)))))
                XI_COUNT=XI_COUNT+1
                nqele=NQNE(ne,XI_COUNT)
              ENDDO
              XI1(1)=(1.0d0/(NQET(NQS(ne))-1.0d0))*(XI_COUNT-1.0d0)
              nb=NBJ(nj1,ne)
              CQ(nm,nq,nx)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '          INP(1,1,nb),nb,1,XI1(1),XE(1,nj1))
C... (KSB 2003)
C... For Navier stokes blood flow solution is elastic tubes at bifurcations:-
c              IF(ITYP2(nr,nx).EQ.5.AND.ITYP3(nr,nx).EQ.1.AND.
c     '          (NXI(-1,0,ne).GT.1)) THEN!.OR.NXI(1,0,ne).GT.1)) THEN
c                IF(ne.NE.NEELEM(1,nr)) THEN !don't do for 1st element
cC... the element is adjacent to a bifurcation thus a discontinuity in
cC... radius and the trace exists between elements.
c                  radius1=XP(1,1,nj1,NPNE(1,nb,ne))
c                  radius2=XP(1,1,nj1,NPNE(NNT(nb),nb,ne))
c                  Ro=MIN(radius1,radius2)                
cC                 TRACE_COMP1=1.0d0-XP(1,1,nj2,NPNE(1,nb,ne))
cC                 TRACE_COMP2=1.0d0-XP(1,1,nj2,NPNE(NNT(nb),nb,ne))
cC                 TRACE_COMP=MAX(TRACE_COMP1,TRACE_COMP2)
c                  CQ(nm,nq,nx)=Ro !unstressed radius
c                ENDIF
c              ENDIF
            ENDDO
          ENDDO
          
        ELSE IF(TYPE(1:8).EQ.'BIDOMAIN') THEN
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)

C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
            IF(ITYP4(nr,nx_t).EQ.4.OR.ITYP4(nr,nx_t).EQ.6.OR.
     '         ITYP4(nr,nx_t).EQ.7) THEN !target class is grid
              DO nq=1,NQT
                CQ(2,nq,nx_t)=CQ(3,nq,nx_s)+CQ(6,nq,nx_s)
                CQ(3,nq,nx_t)=CQ(4,nq,nx_s)+CQ(7,nq,nx_s)
                CQ(4,nq,nx_t)=CQ(5,nq,nx_s)+CQ(8,nq,nx_s)
              ENDDO !nq

            ELSE !target class is finite element mesh
              IF(ILP(3,1,nr,nx_s).EQ.1.OR.ILP(3,1,nr,nx_s).EQ.2)
     '          THEN !elem
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  CE(2,ne,nx_t)=CE(3,ne,nx_s)+CE(6,ne,nx_s)
                  CE(3,ne,nx_t)=CE(4,ne,nx_s)+CE(7,ne,nx_s)
                  CE(4,ne,nx_t)=CE(5,ne,nx_s)+CE(8,ne,nx_s)
                ENDDO !noelem
              ELSE IF(ILP(3,1,nr,nx_s).EQ.3) THEN !node
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  CP(2,np,nx_t)=CP(3,np,nx_s)+CP(6,np,nx_s)
                  CP(3,np,nx_t)=CP(4,np,nx_s)+CP(7,np,nx_s)
                  CP(4,np,nx_t)=CP(5,np,nx_s)+CP(8,np,nx_s)
                ENDDO !nonode
              ENDIF !ILP
            ENDIF !ityp4

          ENDDO !no_nrlist

        ELSE IF(TYPE(1:7).EQ.'INFARCT') THEN
          IF(DOP) THEN
            WRITE(OP_STRING,'(''   CENTRE:'',2E12.3)')
     '        CENTRE(1),CENTRE(2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' SEMIAXIS:'',2E12.3)')
     '        SEMIAXIS(1),SEMIAXIS(2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' MATERIAL:'',2E12.3)')
     '        MATERIAL(1),MATERIAL(2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' FRACTION:'', E12.3)') FRACTION
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !dop
C         Choose border to be fraction of largest semiaxis
          IF(SEMIAXIS(1).GE.SEMIAXIS(2)) THEN
            border=FRACTION*SEMIAXIS(1)
          ELSE IF(SEMIAXIS(2).GT.SEMIAXIS(1)) THEN
            border=FRACTION*SEMIAXIS(2)
          ENDIF
          DO nq=1,NQT
            theta=ATAN2((XQ(2,nq)-CENTRE(2))/SEMIAXIS(2),
     '                  (XQ(1,nq)-CENTRE(1))/SEMIAXIS(1))
C            WRITE(*,'('' theta='',E12.3)') theta
            ELLIPSE(1)=SEMIAXIS(1)*DCOS(theta) !x-coord on ellipse
            ELLIPSE(2)=SEMIAXIS(2)*DSIN(theta) !y-coord on ellipse
C            WRITE(*,'(''  ELLIPSE:'',2E12.3)') ELLIPSE(1),ELLIPSE(2)
            radius1=DSQRT(ELLIPSE(1)**2+ELLIPSE(2)**2) !is outer infarct border
            radius2=DSQRT((XQ(1,nq)-CENTRE(1))**2
     '                   +(XQ(2,nq)-CENTRE(2))**2)
            IF(radius2.LT.radius1) THEN
              WRITE(OP_STRING,'('' nq='',I7,'
     '          //''' is inside outer border of infarct'')') nq
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              IF(radius2.LT.radius1-border) THEN !inside inner border
                CQ(nm,nq,nx)=MATERIAL(1)
              ELSE !within infarct border zone
                WRITE(OP_STRING,'('' & inside infarct border zone'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                Xi=(radius2-radius1+border)/border
                CQ(nm,nq,nx)=(1.d0-Xi)*MATERIAL(1)+Xi*MATERIAL(2)
              ENDIF !infarct border zone
            ENDIF !inside ellipse
          ENDDO !nq

        ELSE IF(TYPE(1:8).EQ.'OPTIMISE') THEN
          nr=1 !tmp
          IF(CBBREV(CO,'CELL',3,noco+1,NTCO,N3CO)) THEN
            CELL=.TRUE.
          ELSE
            CELL=.FALSE.
          ENDIF

          IF(.NOT.CELL) THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,1)
              IF(KTYP27.EQ.1) THEN
                CE(1,ne,nx)=PAOPTI(1)
                WRITE(OP_STRING,'('' Updated opt '
     '            //'var 1 to '',F12.4)') PAOPTI(1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ELSE IF(KTYP27.EQ.2) THEN
                CE(2,ne,nx)=PAOPTI(1)
                WRITE(OP_STRING,'('' Updated opt '
     '            //'var 1 to '',F12.4)') PAOPTI(1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ELSE IF(KTYP27.EQ.3) THEN
                CE(1,ne,nx)=PAOPTI(1)
                CE(2,ne,nx)=PAOPTI(2)
                WRITE(OP_STRING,'('' Updated opt '
     '            //'var 1 to '',F12.4)') PAOPTI(1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' Updated opt '
     '            //'var 2 to '',F12.4)') PAOPTI(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ELSE IF(KTYP27.EQ.4) THEN
                CE(3,ne,nx)=PAOPTI(1)
                WRITE(OP_STRING,'('' Updated opt '
     '            //'var 1 to '',F12.4)') PAOPTI(1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ELSE IF(KTYP27.EQ.5) THEN
                CE(1,ne,nx)=PAOPTI(1)
                CE(2,ne,nx)=PAOPTI(2)
                CE(3,ne,nx)=PAOPTI(3)
                WRITE(OP_STRING,'('' Updated opt '
     '            //'var 1 to '',F12.4)') PAOPTI(1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' Updated opt '
     '            //'var 2 to '',F12.4)') PAOPTI(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' Updated opt '
     '            //'var 3 to '',F12.4)') PAOPTI(3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !neelem

          ELSE

            IF(CBBREV(CO,'VARIANT',3,noco+1,NTCO,N3CO)) THEN
              variant=IFROMC(CO(N3CO+1))
            ELSE
              variant=1
            ENDIF

            IF(CBBREV(CO,'PARAMETER',3,noco+1,NTCO,N3CO)) THEN
              param=IFROMC(CO(N3CO+1))
            ELSE
              param=1
            ENDIF

            IF(CBBREV(CO,'OPTIMISATION_VARIABLE',16,
     '        noco+1,NTCO,N3CO)) THEN
              opt_var=IFROMC(CO(N3CO+1))
            ELSE
              opt_var=1
            ENDIF

            !set grid point values
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              DO nqq=1,NQET(NQS(ne))
                nq=NQNE(ne,nqq) ! global grid point number
                SPATIALLY_VARYING=.FALSE.
                DO nqrsv=1,IRCQS_SPATIAL(0)
                  IF(IRCQS_SPATIAL(nqrsv).EQ.param) THEN
                    SPATIALLY_VARYING=.TRUE.
                    nqrs=nqrsv
                  ENDIF
                ENDDO
                IF(SPATIALLY_VARYING) THEN
                  RCQS_SPATIAL(nqrs,nq)=PAOPTI(opt_var)
                ELSE
                  CELL_RCQS_VALUE(param,variant)=PAOPTI(opt_var)
                ENDIF
              ENDDO ! nqq
            ENDDO !noelem
          ENDIF

        ENDIF
      ENDIF

      CALL EXITS('UPMATE')
      RETURN
 9999 CALL ERRORS('UPMATE',ERROR)
      CALL EXITS('UPMATE')
      RETURN 1
      END


