      SUBROUTINE IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,il,ILPIN,
     '  IRCQS_SPATIAL,NBJ,NEELEM,NELIST,NGLIST,NMBIN,NPLIST,NPNODE,NQET,
     '  NQNE,NQS,NQXI,nr,nx,CE,CELL_RCQS_VALUE,CGE,CIN,CP,RCQS_SPATIAL,
     '  G_RDEFLT,XIG,G_FORMAT,SETCONSTIT,ERROR,*)

C#### Subroutine: IPMAT5_SINGLE
C###  Description:
C###    IPMAT5_SINGLE inputs a single material parameter with index il

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp50.cmn'
!     Parameter List
      INTEGER GRNGLIST(0:NEGM),ICQS_SPATIAL(NQISVM,NQM),
     '  il,ILPIN(NMM),IRCQS_SPATIAL(0:NQRSVM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NGLIST(0:NGM),
     '  NMBIN(NMM),NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),
     '  NQNE(NEQM,NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),nr,nx
      REAL*8 CE(NMM,NEM),CELL_RCQS_VALUE(NQRM,NQVM),
     '  CGE(NMM,NGM,NEM),CIN(NMM,0:NGM,NNEPM),
     '  CP(NMM,NPM),G_RDEFLT(IORMX),
     '  RCQS_SPATIAL(NQRSVM,NQM),XIG(NIM,NGM,NBM)
      CHARACTER*(*) G_FORMAT,ERROR
      LOGICAL SETCONSTIT
!     Local Variables
      INTEGER IBEG2,ICHAR,IEND2,IJF,IJS,IKF,IKS,INFO,
     '  index,n,nb,N1,ne,ne1,ng,ni,ni1,ni2,ni3,noelem,
     '  nong,nonode,NOQUES,np,np1,nq,nqrs,nqq,nqrsv,num_points,
     '  param_num,variant
      REAL*8 DS,DSH,XI(3)
      CHARACTER CHAR2*100
      LOGICAL ALLSET,FILEIP,INLIST,SPATIALLY_VARYING

      CALL ENTERS('IPMAT5_SINGLE',*9999)

      ICHAR=0 !temporary (needs deleting later)

      FILEIP=.FALSE.
      NOQUES=0

C CS 21/8/00 new assert
      CALL ASSERT(USE_NONLIN.GT.0,'>>Set USE_NONLIN to 1 '
     '  //'in parameters file',ERROR,*9999)



C     G_FORMAT for parameter question is set in IPMAT5
      IF(IOTYPE.EQ.3) IDATA(1)=ILPIN(il)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  G_FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,5,
     '  LDATA,LDEFLT,RDATA,G_RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
C CS new 17/5/99 should work now
C        CALL ASSERT(IDATA(1).LE.3,'>>Material params '
C     '    //'cannot vary by Gauss points',ERROR,*9999)
        ILPIN(il)=IDATA(1)
        IF(SETCONSTIT) ILP(il,1,nr,nx)=ILPIN(il)
      ENDIF
C AWCL & MPN 26/07/2010: check ippara file
      IF(ILPIN(il).EQ.4.OR.ILPIN(il).EQ.5) THEN
        CALL ASSERT(USE_GAUSS_PT_MATERIALS.GT.0,
     '    '>>Set USE_GAUSS_PT_MATERIALS to 1 in parameters file',
     '    ERROR,*9999)
      ENDIF

      IF(ILPIN(il).EQ.1) THEN !constant spatially
C       G_RDEFLT for parameter value is set in IPMAT5
C        CHAR2=CFROMR(G_RDEFLT(1),'(D12.3)')
        WRITE(CHAR2,'(D12.3)') G_RDEFLT(1)
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        FORMAT='($,'' The value is ['//CHAR2(IBEG2:IEND2)
     '    //']: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=CIN(il,0,NEELEM(1,nr))
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,G_RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            CIN(il,0,ne)=RDATA(1)
            IF(SETCONSTIT) CE(il,ne)=CIN(il,0,ne)
          ENDDO !noelem (ne)
        ENDIF

      ELSE IF(ILPIN(il).EQ.2) THEN !defined by elements
C       Prompt for element parameters allowing for element group input
        IF(IOTYPE.EQ.3) THEN
          noelem=0 !init element loop for writing params
        ELSE IF(IOTYPE.EQ.3) THEN
C         init CIN so can check that each element has been set later
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            CIN(il,0,ne)=-RMAX
          ENDDO !noelem (ne)
        ENDIF
 6800   FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
        IF(IOTYPE.EQ.3) THEN
          noelem=noelem+1
          IF(noelem.LE.NEELEM(0,nr)) THEN
            ne=NEELEM(noelem,nr)
            IDATA(1)=ne
          ELSE
            IDATA(1)=0 !to terminate element loop
          ENDIF
        ENDIF
 6900   CDATA(1)='ELEMENTS' !for use with group input
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NET(nr),
     '    LDATA,LDEFLT,RDATA,G_RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IDATA(1).NE.0) THEN !not default exit
          NELIST(0)=IDATA(0)
          DO n=1,IDATA(0)
            NELIST(n)=IDATA(n)
            ne=IDATA(n)
            IF(.NOT.INLIST(ne,NEELEM(1,nr),
     '        NEELEM(0,nr),N1)) THEN
              WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '          //'in the current region'')') ne
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              GOTO 6900
            ENDIF
          ENDDO !n
C         Define parameter for first element in group
          ne=NELIST(1) !rest of group filled at end
          FORMAT='($,'' The value is [0]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=CIN(il,0,ne)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '      IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            CIN(il,0,ne)=RDATA(1)
            IF(SETCONSTIT) CE(il,ne)=CIN(il,0,ne)
C           Set parameters for rest of elements in group
            DO n=2,NELIST(0)
              ne1=NELIST(n)
              CIN(il,0,ne1)=CIN(il,0,ne)
              IF(SETCONSTIT) CE(il,ne1)=CIN(il,0,ne1)
            ENDDO !n
          ENDIF
          GO TO 6800 !for more elements
        ENDIF !idata(1).NE.0

        IF(IOTYPE.NE.3) THEN
C         check that CIN has been set for each element
          ALLSET=.TRUE.
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(CIN(il,0,ne).EQ.-RMAX) THEN
              ALLSET=.FALSE.
              CIN(il,0,ne)=0.0d0
            ENDIF
          ENDDO !noelem (ne)
          IF(.NOT.ALLSET) THEN
            WRITE(OP_STRING,'('' >>WARNING: Parameter has been '
     '        //'assigned zero for elements not specified'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(ILPIN(il).EQ.3) THEN !defined by nodes
        FORMAT='($,'' Enter basis type number [1]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NMBIN(il)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBT,
     '    LDATA,LDEFLT,RDATA,G_RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NMBIN(il)=IDATA(1)
          IF(SETCONSTIT) NMB(il,nr,nx)=NMBIN(il)
        ENDIF

C       Prompt for nodal parameters allowing for node group input
        IF(IOTYPE.EQ.3) THEN
          nonode=0 !init node loop for writing params
        ELSE IF(IOTYPE.NE.3) THEN
C         init CIN so can check that each node has been set later
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            CIN(il,0,np)=-RMAX
          ENDDO !nonode (np)
        ENDIF
 7000   FORMAT='($,'' Enter node #s/name [EXIT]: '',I5)'
        IF(IOTYPE.EQ.3) THEN
          nonode=nonode+1
          IF(nonode.LE.NPNODE(0,nr)) THEN
            np=NPNODE(nonode,nr)
            IDATA(1)=np
          ELSE
            IDATA(1)=0 !to terminate node loop
          ENDIF
          IDATA(0)=1 !write out one node at a time
        ENDIF !iotype=3
 7100   CDATA(1)='NODES' !for use with group input
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NPT(nr),
     '    LDATA,LDEFLT,RDATA,G_RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IDATA(1).NE.0) THEN !not default exit
          NPLIST(0)=IDATA(0)
          DO n=1,IDATA(0)
            NPLIST(n)=IDATA(n)
            np=IDATA(n)
            IF(.NOT.INLIST(np,NPNODE(1,nr),
     '        NPNODE(0,nr),N1)) THEN
              WRITE(OP_STRING,'('' Node '',I5,'' does not '
     '          //'belong to the current region'')') np
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              GOTO 7100
            ENDIF
          ENDDO !n
C         Define parameter for first node in group
          np=NPLIST(1) !rest of group filled at end
          FORMAT='($,'' The value is [0]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=CIN(il,0,np)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            CIN(il,0,np)=RDATA(1)
            IF(SETCONSTIT) CP(il,np)=CIN(il,0,np)
C           Set parameters for rest of nodes in group
            DO n=2,NPLIST(0)
              np1=NPLIST(n)
              CIN(il,0,np1)=CIN(il,0,np)
              IF(SETCONSTIT) CP(il,np1)=CIN(il,0,np1)
            ENDDO !n
          ENDIF
          GO TO 7000 !for more nodes
        ENDIF !idata(1).NE.0

        IF(IOTYPE.NE.3) THEN
C         check that CIN has been set for each node
          ALLSET=.TRUE.
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            IF(CIN(il,0,np).EQ.-RMAX) THEN
              ALLSET=.FALSE.
              CIN(il,0,np)=0.0d0
            ENDIF
          ENDDO !nonode (np)
          IF(.NOT.ALLSET) THEN
            WRITE(OP_STRING,'('' >>WARNING: Parameter has been '
     '      //'assigned zero for nodes not specified'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF((ILPIN(il).EQ.4).OR.(ILPIN(il).EQ.5)) THEN !defined by Gauss or grid points
        IF(IOTYPE.EQ.3) THEN
C         Init CIN so can check that each element has been set later
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NBJ(1,ne)
            DO ng=1,NGT(nb)
              CIN(il,ng,ne)=-RMAX
            ENDDO !ng
          ENDDO !noelem (ne)
        ENDIF


        IF(ILPIN(il).EQ.4) THEN ! defined by Gauss points
          FORMAT='('' Select Gauss points by [1]: '''//
     '      '/''   (1) Named Gauss point groups'''//
     '      '/''   (2) Elements and local Gauss point numbers'''//
     '      '/$,''    '',I1)'
          IF (IOTYPE.EQ.3) THEN
            IDEFLT(1)=2
          ELSE
            IDEFLT(1)=1
          ENDIF
          IDATA(1)=IDEFLT(1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '      LDATA,LDEFLT,RDATA,G_RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP5E(nr)=IDATA(1)

          IF(KTYP5E(nr).EQ.1) THEN
            IF(IOTYPE.NE.3) THEN
              IDATA(1)=0
            ENDIF
 9800       FORMAT='($,'' Enter Gauss point group name [EXIT]: '',I5)'
            CDATA(1)='GAUSS' !for use with group input
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '        IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,G_RDEFLT,RMIN,
     '        RMAX,INFO,ERROR,*9999)
            IF(IDATA(1).NE.0) THEN !not default exit
              GRNGLIST(0)=IDATA(0)
              index=1
              DO noelem=1,GRNGLIST(0)
                GRNGLIST(index)=IDATA(index)
                IF(.NOT.INLIST(GRNGLIST(index),NEELEM(1,nr),
     '            NEELEM(0,nr),N1)) THEN
                  WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '              //'in the current region'')') ne
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  GOTO 9800
                ENDIF
                GRNGLIST(index+1)=IDATA(index+1)
                DO num_points=1,IDATA(index+1)
                  IDATA(index+1+num_points)=IDATA(index+1+num_points)
                ENDDO !num_points
                index=index+2+num_points
              ENDDO !noelem

              WRITE(CHAR2,'(D12.3)') G_RDEFLT(1)
              CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
              FORMAT='($,'' The value is '
     '          //'['//CHAR2(IBEG2:IEND2)
     '          //']: '',D11.4)'
              IF(IOTYPE.EQ.3) RDATA(1)=CIN(il,ng,ne)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '          FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '          IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,G_RDEFLT,
     '          -RMAX,RMAX,INFO,ERROR,*9999)
              index=1
              DO noelem=1,GRNGLIST(0)
                ne=GRNGLIST(index)
                DO num_points=1,GRNGLIST(index+1)
                  ng=IDATA(index+1+num_points)
C               Set up CIN
                  IF(IOTYPE.NE.3) THEN
                    CIN(il,ng,ne)=RDATA(1)
                    IF(SETCONSTIT) CGE(il,ng,ne)=CIN(il,ng,ne)
                  ENDIF
                ENDDO !num_points
                index=index+2+num_points
              ENDDO !noelem
              GOTO 9800
            ENDIF

          ELSE
C         Prompt for element parameters allowing for element group input
            IF(IOTYPE.EQ.3) THEN
              noelem=0 !init element loop for writing params
            ENDIF
 7800       FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
            IF(IOTYPE.EQ.3) THEN
              noelem=noelem+1
              IF(noelem.LE.NEELEM(0,nr)) THEN
                ne=NEELEM(noelem,nr)
                IDATA(0)=1
                IDATA(1)=ne
              ELSE
                IDATA(0)=0 !to terminate element loop
                IDATA(1)=0 !to terminate element loop
              ENDIF
            ENDIF
 7900       CDATA(1)='ELEMENTS' !for use with group input
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,
     '        NET(nr),LDATA,LDEFLT,RDATA,G_RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
            IF(IDATA(1).NE.0) THEN !not default exit
              NELIST(0)=IDATA(0)
              DO n=1,IDATA(0)
                NELIST(n)=IDATA(n)
                ne=IDATA(n)
                IF(.NOT.INLIST(ne,NEELEM(1,nr),
     '            NEELEM(0,nr),N1)) THEN
                  WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '              //'in the current region'')') ne
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  GOTO 7900
                ENDIF
              ENDDO !n

C           Prompt for Gauss pt params allowing for Gauss pt group input
              IF(IOTYPE.EQ.3) THEN
                nong=0
              ENDIF
 8800         FORMAT='($,'' Gauss Point #s[EXIT]: '',I5)'
              IF(IOTYPE.EQ.3) THEN
                nong=nong+1
                IF(nong.LE.NGT(nb)) THEN
                  ng=nong
                  IDATA(0)=1
                  IDATA(1)=ng
                ELSE
                  IDATA(1)=0 !to terminate gauss point loop
                  IDATA(0)=0 !to terminate gauss point loop
                ENDIF
              ENDIF
 8900         CDATA(1)='GAUSS_POINTS' !for use with group input
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NGM,
     '          LDATA,LDEFLT,RDATA,G_RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IDATA(1).NE.0) THEN !not default exit
                NGLIST(0)=IDATA(0)
                DO n=1,IDATA(0)
                  ng=IDATA(n)
                  NGLIST(n)=ng
                  nb=NBJ(1,NELIST(1))
                  IF(ng.GT.NGT(nb)) THEN
                    WRITE(OP_STRING,'('' Gauss Point '',I5,'' is not '
     '                //'in the current basis'')') ng
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    GOTO 8900
                  ENDIF
                ENDDO

                WRITE(CHAR2,'(D12.3)') G_RDEFLT(1)
                CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
                FORMAT='($,'' The value is '
     '            //'['//CHAR2(IBEG2:IEND2)
     '            //']: '',D11.4)'
                IF(IOTYPE.EQ.3) THEN
                  CIN(il,ng,ne)=CGE(il,ng,ne)
                  RDATA(1)=CGE(il,ng,ne)
                ENDIF
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '            FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '            IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,G_RDEFLT,
     '            -RMAX,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  DO n=1,NGLIST(0)
                    ng=NGLIST(n)
C                   Set up CIN
                    DO noelem=1,NELIST(0)
                      ne=NELIST(noelem)
                      IF(IOTYPE.NE.3) THEN
                        CIN(il,ng,ne)=RDATA(1)
                        IF(SETCONSTIT) CGE(il,ng,ne)=CIN(il,ng,ne)
                      ENDIF
                    ENDDO !noelem (ne)
                  ENDDO
                ENDIF
                GOTO 8800 ! Get more Gauss points
              ELSE
                GOTO 7800
              ENDIF
            ENDIF !idata(1).NE.0
          ENDIF

        ELSE ! define by Grid points
          FORMAT='(/$,'' Enter the grid points '
     '      //'parameter number [1]: '',I3)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=1
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '      LDATA,LDEFLT,RDATA,G_RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            param_num=IDATA(1)
C This could be done this using gengrid map but this seems unnecessarily
C complicated and inefficient. For now just use the closest grid point
!            IF(.NOT.CALL_GMAP) THEN
!              CALL GEN_GRID_MAP(IBT,IDO,INP,NQSCNB,NQXI,PGNQE,XIG,
!     '          ERROR,*9999)
!              CALL_GMAP=.TRUE.
!            ENDIF

            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              DO ng=1,NGT(nb)
C             Set up CIN
                CIN(il,ng,ne)=0.0d0
                DO nqq=1,NQET(NQS(ne))
C This is largely pulled from GEN_GRID_MAP()
! Find the grid point which is closest to the current gauss point
! in xi space
                  DSH=RMAX
                  IJF=1
                  IKF=1
                  IJS=1
                  IKS=1
                  IF(NIT(nb).GE.2) THEN
                    IJF=NQXI(2,NQS(ne))-1
                    IJS=2
                  ENDIF
                  IF(NIT(nb).EQ.3) THEN
                    IKF=NQXI(3,NQS(ne))-1
                    IKS=2
                  ENDIF
! Loop over the internal grid points, if we restrict ourselves to
! the internal points in each scheme then we can always make a
! local quadratic element without going over global element boundaries.
                  DO ni3=IKS,IKF
                    DO ni2=IJS,IJF
                      DO ni1=2,NQXI(1,NQS(ne))-1
                        XI(1)=DBLE(ni1-1)/DBLE(NQXI(1,NQS(ne))-1)
                        IF(IJS.NE.1) XI(2)=DBLE(ni2-1)/
     '                    DBLE(NQXI(2,NQS(ne))-1)
                        IF(IKS.NE.1) XI(3)=DBLE(ni3-1)/
     '                    DBLE(NQXI(3,NQS(ne))-1)
                        DS=0.0d0
                        DO ni=1,NIT(nb)
                          DS=DS+DABS(XIG(ni,ng,nb)-XI(ni))**2
                        ENDDO !ni
                        DS=DSQRT(DS)
                        IF(DS.LT.DSH) THEN
                          DSH=DS
                          nq=ni1+((ni2-1)*NQXI(1,NQS(ne)))+
     '                      ((ni3-1)*NQXI(1,NQS(ne))*
     '                      NQXI(2,NQS(ne))) !store local grid point number
                        ENDIF
                      ENDDO !ni1
                    ENDDO !ni2
                  ENDDO !ni3
                ENDDO ! nqq

                CALL ASSERT(DSH.LT.RMAX,'>>No grid points found',
     '            ERROR,*9999)
                IF(DOP) THEN
                  WRITE(OP_STRING,'(''Grid: '',I6,'' Gauss: '',I6)')
     '              nq,ng
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF

                nqq=NQNE(ne,nq) ! global grid point number
                variant=ICQS_SPATIAL(1,nqq)
                SPATIALLY_VARYING=.FALSE.
                DO nqrsv=1,IRCQS_SPATIAL(0)
                  IF(IRCQS_SPATIAL(nqrsv).EQ.param_num) THEN
                    SPATIALLY_VARYING=.TRUE.
                    nqrs=nqrsv
                  ENDIF
                ENDDO
                IF(SPATIALLY_VARYING) THEN
                  CIN(il,ng,ne)=RCQS_SPATIAL(nqrs,nqq)
                ELSE
                  CIN(il,ng,ne)=CELL_RCQS_VALUE(param_num,variant)
                ENDIF

                IF(SETCONSTIT) CGE(il,ng,ne)=CIN(il,ng,ne)
              ENDDO !ng
            ENDDO !noelem
          ENDIF !IOTYPE.NE.3
        ENDIF

        IF(IOTYPE.NE.3) THEN
C         Check that CIN has been set for each element and Gauss pt
          ALLSET=.TRUE.
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NBJ(1,ne)
            DO ng=1,NGT(nb)
              IF(CIN(il,ng,ne).EQ.-RMAX) THEN
                ALLSET=.FALSE.
                CIN(il,ng,ne)=0.0d0
              ENDIF
            ENDDO !ng
          ENDDO !noelem (ne)
          IF(.NOT.ALLSET) THEN
            WRITE(OP_STRING,'('' >>WARNING: Parameter has been '
     '        //'assigned zero for elements not specified'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('IPMAT5_SINGLE')
      RETURN
 9999 CALL ERRORS('IPMAT5_SINGLE',ERROR)
      CALL EXITS('IPMAT5_SINGLE')
      RETURN 1
      END


