      SUBROUTINE IPACTI(NBH,NBJ,NEELEM,NELIST,NGLIST,nr,nx,FEXT,ERROR,*)

C#### Subroutine: IPACTI
C###  Description:
C###    IPACTI inputs active muscle contraction parameters.
C###    Set up in subroutine IPACTI.
      
      IMPLICIT NONE
      INCLUDE 'acti00.cmn'
      INCLUDE 'acti01.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ptr00.cmn'

!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NGLIST(0:NGM),nr,nx
      REAL*8 FEXT(NIFEXTM,NGM,NEM)
      CHARACTER ERROR*(*)
      
!     Local Variables
      INTEGER  IBEG1,IBEG2,IBEGA,IBEG_T1,ICHAR,IEND1,IEND2,IENDA,
     '  IEND_T1,INFO,n,N1,nb,ne,ng,NO_ACTI,noelem,
     '  nong,NOQUES
      CHARACTER ACTIVE_STR*30,TITLE1*200,CHAR2*100
      LOGICAL ALLSET,FILEIP,INLIST

      CALL ENTERS('IPACTI',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      nb=0
      nong=0
      KTYP59S(nr)=1             ! DEFAULT IS THAT THE ACTIVE STRESS 
                                ! COMPONENT IS DEFINED AS A CAUCHY 
                                ! (TRUE) TYPE OF STRESS

      IF(ITYP1(nr,nx).EQ.5) THEN !finite deformation problem
        CALL ASSERT(KTYP53(nr).EQ.3,
     '    ' >>Stresses must be referred to body coords with active '
     '    //'fibre stress in Define material',ERROR,*9999)
      ENDIF

C EWR 04Mar2004: Adding regional variation of activation
      TITLE1='/''  (1) Constant spatially'''//
     '  '/''  (2) Piecewise constant (defined by elements)'''//
     '  '/''  (3) Defined by Gauss points'''
      CALL STRING_TRIM(TITLE1,IBEG_T1,IEND_T1)

C OR 15-08-06 
C     Added in additional functionality such that one can choose
C     from a  cellml function definition and add that particular
C     value to a position within the 2nd PKST.
C     '  //'/''   (3) Unused'''
C news MPN 23May2000: adding HMT (several changes below)
C news MPN/VYW 6Aug2014: changed option 4 for TCa input
      FORMAT='('' Specify type of contraction model [2]: '''
     '  //'/''   (1) SS tension-length-Ca relation (set Cai)'''
     '  //'/''   (2) Get active stress from cell/grid model '
     '  //' (e.g. Hunter/McCulloch/ter Keurs model)'''
     '  //'/''   (3) Enhance 2nd PKST by user-defined functions'''
     '  //'/''   (4) SS TCa-length relation (set TCa)'''
     '  //'/$,''    '',I1)'
      IDEFLT(1)=2
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP59(nr)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP59(nr)=IDATA(1)

      IF(KTYP59(nr).EQ.1.OR.
     '  KTYP59(nr).EQ.4) THEN      !SS tension-length-Ca relation

      IF(KTYP59(nr).EQ.1) THEN
         RDEFLT(1)=1.0d2
        FORMAT='($,'' Enter max isometric tension at ext.ratio=1 '
     '    //'(Tref) [100kPa]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=Tref
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Tref=RDATA(1)
      ENDIF

        RDEFLT(1)=1.45d0
        FORMAT='($,'' Enter non-dimensional slope parameter (beta)'
     '    //' [1.45]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=T0_beta
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) T0_beta=RDATA(1)

      IF(KTYP59(nr).EQ.1) THEN
        RDEFLT(1)=0.5d0
        FORMAT='($,'' Enter c50 for [Ca]i saturation curve (0<c<1)'
     '    //' [0.5]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=Ca_c50
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Ca_c50=RDATA(1)

        RDEFLT(1)=3.0d0
        FORMAT='($,'' Enter Hill coeff. for [Ca]i saturation curve'
     '    //' (h) [3.0]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=Ca_h
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Ca_h=RDATA(1)

        RDEFLT(1)=1.0d0
        FORMAT='($,'' Enter max [Ca]i (Ca_max) [1]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=Ca_max
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Ca_max=RDATA(1)

        ACTIVE_STR='calcium level [Ca]i'
      ELSE
        ACTIVE_STR='tension level TCa'
      ENDIF
      CALL STRING_TRIM(ACTIVE_STR,IBEGA,IENDA)

C EWR 04Mar2004 Regional variation of activation
        IDEFLT(1)=1 ! Default: Constant spatially
        FORMAT='('' Specify whether the initial '
     '    //ACTIVE_STR(IBEGA:IENDA)//' is [1]:'''
     '    //TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP5J(nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP5J(nr)=IDATA(1)
        IF(KTYP5J(nr).EQ.1) THEN !constant spatially
          RDEFLT(1)=0.0d0
          FORMAT='($,'' Enter initial '
     '      //ACTIVE_STR(IBEGA:IENDA)//' [0]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=FEXT(4,1,1) !1st Gauss pt in 1st elem
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     '      IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO ng=1,NGT(NBH(NH_LOC(1,nx),1,ne))
                FEXT(4,ng,ne)=RDATA(1)
              ENDDO !ng
            ENDDO !noelem (ne)
          ENDIF
        ELSE IF(KTYP5J(nr).EQ.2) THEN !defined by elements
C       Prompt for element parameters allowing for element group input
          IF(IOTYPE.EQ.3) THEN
            noelem=0 !init element loop for writing params
          ELSE IF(IOTYPE.NE.3) THEN
C           init FEXT so can check that each element has been set later
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO ng=1,NGT(NBH(NH_LOC(1,nx),1,ne))
                FEXT(4,ng,ne)=-RMAX
              ENDDO !ng
            ENDDO !noelem (ne)
          ENDIF
 6800     FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
          IF(IOTYPE.EQ.3) THEN
            noelem=noelem+1
            IDATA(0)=1 ! write out one element at a time
            IF(noelem.LE.NEELEM(0,nr)) THEN
              ne=NEELEM(noelem,nr)
              IDATA(1)=ne
            ELSE
              IDATA(1)=0 !to terminate element loop
            ENDIF
          ENDIF
 6900     CDATA(1)='ELEMENTS' !for use with group input
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NET(nr),
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).NE.0) THEN !not default exit
            NELIST(0)=IDATA(0)
            DO n=1,IDATA(0)
              NELIST(n)=IDATA(n)
              ne=IDATA(n)
              IF(.NOT.INLIST(ne,NEELEM(1,nr),
     '          NEELEM(0,nr),N1)) THEN
                WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '            //'in the current region'')') ne
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                GOTO 6900
              ENDIF
            ENDDO !n
C           Define parameter for first element in group
            FORMAT='($,'' The value is [0]: '',D11.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=FEXT(4,1,NELIST(1)) ! changed ne to NELIST(1)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '        IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO n=1,NELIST(0)
                ne=NELIST(n)
                DO ng=1,NGT(NBH(NH_LOC(1,nx),1,ne))
                  FEXT(4,ng,ne)=RDATA(1)
                ENDDO !ng
              ENDDO !n
            ENDIF
            GO TO 6800 !for more elements
          ENDIF !idata(1).NE.0

          IF(IOTYPE.NE.3) THEN
C         check that FEXT has been set for each element
            ALLSET=.TRUE.
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO ng=1,NGT(NBH(NH_LOC(1,nx),1,ne))
                IF(FEXT(4,ng,ne).EQ.-RMAX) THEN
                  ALLSET=.FALSE.
                  FEXT(4,ng,ne)=0.0d0
                ENDIF
              ENDDO !ng
            ENDDO !noelem (ne)
            IF(.NOT.ALLSET) THEN
              WRITE(OP_STRING,'('' >>WARNING: Parameter has been '
     '          //'assigned zero for elements not specified'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF

C        ELSE IF(ILPIN(il).EQ.3) THEN !defined by nodes - Not implemented
        ELSE IF(KTYP5J(nr).EQ.3) THEN ! defined by Gauss points
C        Prompt for element parameters allowing for element group input
          IF(IOTYPE.EQ.3) THEN
            noelem=0 !init element loop for writing params
          ENDIF
 7800     FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
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
 7900     CDATA(1)='ELEMENTS' !for use with group input
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,
     '      NET(nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '      ERROR,*9999)
          IF(IDATA(1).NE.0) THEN !not default exit
            NELIST(0)=IDATA(0)
            DO n=1,IDATA(0)
              NELIST(n)=IDATA(n)
              ne=IDATA(n)
              IF(.NOT.INLIST(ne,NEELEM(1,nr),
     '          NEELEM(0,nr),N1)) THEN
                WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '            //'in the current region'')') ne
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                GOTO 7900
              ENDIF
            ENDDO !n

C         Prompt for Gauss pt params allowing for Gauss pt group input
            IF(IOTYPE.EQ.3) THEN
              nong=0
              nb=NBJ(1,ne) ! Is this right? need to define nb
            ENDIF
 8800       FORMAT='($,'' Gauss Point #s[EXIT]: '',I5)'
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
 8900       CDATA(1)='GAUSS_POINTS' !for use with group input
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NGM,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IDATA(1).NE.0) THEN !not default exit
              NGLIST(0)=IDATA(0)
              DO n=1,IDATA(0)
                ng=IDATA(n)
                NGLIST(n)=ng
                nb=NBJ(1,NELIST(1))
                IF(ng.GT.NGT(nb)) THEN
                  WRITE(OP_STRING,'('' Gauss Point '',I5,'' is not '
     '              //'in the current basis'')') ng
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  GOTO 8900
                ENDIF
              ENDDO

              WRITE(CHAR2,'(D12.3)') RDEFLT(1)
              CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
              FORMAT='($,'' The value is '
     '          //'['//CHAR2(IBEG2:IEND2)
     '          //']: '',D11.4)'
              IF(IOTYPE.EQ.3) RDATA(1)=FEXT(4,ng,ne)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '          FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '          IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '          -RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                DO n=1,NGLIST(0)
                  ng=NGLIST(n)
                  DO noelem=1,NELIST(0)
                    ne=NELIST(noelem)
                    IF(IOTYPE.NE.3) FEXT(4,ng,ne)=RDATA(1)
                  ENDDO !noelem (ne)
                ENDDO
              ENDIF
              GOTO 8800 ! Get more Gauss points
            ELSE
              GOTO 7800
            ENDIF
          ENDIF !idata(1).NE.0

C        ELSE ! define by Grid points - not implemented
        ENDIF
C EWR 04Mar2004 Endo of regional variation of activation section

      ELSE IF(KTYP59(nr).EQ.2) THEN  !Dynamic HMT (from grid problem)

C       No active material properties are defined here as they are
C       set up in IPCELL for HMT.

      ELSE IF(KTYP59(nr).EQ.3) THEN  !Adding userspecified values to
                                ! the 2nd PKST. The user defined values
                                ! are defined within cellml.
C OR 15-08-06
C 
        FORMAT='(/'' Specify type of stress being added to the'
     &       //' 2nd PKST [1]: '''
     &       //'/''   (1) Cauchy-type stress'''
     &       //'/''   (2) 2nd Piola-Kirchhoff-type stress'''
     &       // '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP59S(nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '       ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     &       LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP59S(nr)=IDATA(1)

  
        FORMAT='(/$,'' Specify tensor component the stress should be' //
     &       ' added to [1,1]: '',I1)'
        
        IDEFLT(1)=1
        IDEFLT(2)=1
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=KTYP59S1(nr)
          IDATA(2)=KTYP59S2(nr)
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
     '       ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     &       LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          KTYP59S1(nr)=IDATA(1)
          KTYP59S2(nr)=IDATA(2)
        ENDIF
        
        CALL ASSERT( ((KTYP59S1(nr).LE.3).AND.(KTYP59S1(nr).GE.1)).AND
     &       .((KTYP59S2(nr).LE.3).AND.(KTYP59S2(nr).GE.1))
     &       ,'>>Not a valid stress tensor component',ERROR,*9999)

        CDATA(1)=' '
        IF (KTYP59S(nr).EQ.1) THEN
          FORMAT='(/$,'' Enter CellML fumction name for' //
     &         ' Cauchy stress component: '',I5)'
        ELSE
          FORMAT='(/$,'' Enter CellML fumction name for' //
     &         ' 2nd Piola Kirchhoff stress component: '',I5)'        
        ENDIF
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &       ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4, LDATA
     &       ,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        CALL STRING_TRIM(CDATA(1),IBEG1,IEND1)
        IF(IOTYPE.NE.3) ASCOMP(nr)=CDATA(1)

        CALL DETERMASC(nr,%VAL(CELL_ICQS_NAMES_PTR)
     &       ,%VAL(CELL_RCQS_NAMES_PTR),%VAL(CELL_YQS_NAMES_PTR),ASCOMP
     &       ,ASC_ARRAYNAME,ASC_CELLVARINDEX,ERROR,*9999)
     
C
C OR 15-08-06
C     NEEDS MAYBE MORE GENERALITY! E.G ASSOCIATE A PRTICULAR STRESS
C     COMPONENT WITH GAUSS/GRID POINTS (GROUPS). 

C old MPN 6Aug2014: the hair cell option is not used anywhere
C      ELSE IF(KTYP59(nr).EQ.4) THEN  !Outer hair cell
C        RDEFLT(1)=1.0D0
C        FORMAT='($,'' Enter the resonant frequency [1 /ms]: '',D11.4)'
C        IF(IOTYPE.EQ.3) RDATA(1)=ACTIVE_FREQUENCY
C        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) ACTIVE_FREQUENCY=2.0D0*PI*RDATA(1) !Conversion to
C                    !  rads CPB 5/4/92

C        RDEFLT(1)=0.0D0
C        FORMAT='($,'' Enter the damping factor [0]: '',D11.4)'
C        IF(IOTYPE.EQ.3) RDATA(1)=ACTIVE_DAMPING
c        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
c     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) ACTIVE_DAMPING=RDATA(1)

C        RDEFLT(1)=0.0D0
C        FORMAT='($,'' Enter the phase shift [0 rad]: '',D11.4)'
C        IF(IOTYPE.EQ.3) RDATA(1)=ACTIVE_PHASE
C        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) ACTIVE_PHASE=RDATA(1)

C        RDEFLT(1)=1.0D0
c        RDEFLT(2)=1.0D0
C        FORMAT='($,'' Enter x & y dir.n amplitudes [1 kPa/m]: '','
C     '    //'2D11.4)'
C        IF(IOTYPE.EQ.3) THEN
c          RDATA(1)=ACTIVE_AMPLITUDE_1
C          RDATA(2)=ACTIVE_AMPLITUDE_2
c        ENDIF
C        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) THEN
C          ACTIVE_AMPLITUDE_1=RDATA(1)
C          ACTIVE_AMPLITUDE_2=RDATA(2)
C        ENDIF

C        NO_ACTI=0
C 100      FORMAT='($,'' Enter an element number for active stiffness'
C     '      //'[exit]: '',I5)'
C          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
c     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NET(1),
C     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c          IF(IDATA(1).NE.0) THEN
C            NO_ACTI=NO_ACTI+1
C            CALL ASSERT(NO_ACTI.LE.10,'>>Max number allowed is 10',
C     '        ERROR,*9999)
C            NE_ACTI(NO_ACTI)=IDATA(1)
C            GO TO 100
C          ENDIF
C        NE_ACTI(0)=NO_ACTI

      ENDIF !KTYP59

      CALL EXITS('IPACTI')
      RETURN
 9999 CALL ERRORS('IPACTI',ERROR)
      CALL EXITS('IPACTI')
      RETURN 1
      END


