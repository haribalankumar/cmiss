      SUBROUTINE GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne_start,NKJ,NKJE,
     '  noelem,nonode,np,NP_INTERFACE,np_start,NPNE,NPNODE,nr,NRE,NVJE,
     &  NVJP,NXI,SE,MAKE,ERROR,*)

C#### Subroutine: GN1DNE
C###  Description:
C###    GN1DNE sets up NEELEM and NPNODE for a new mesh 'branch'.  Then
C###    GN1DNEJ is called to set up region, basis function, and scale
C###    factor arrays for the new element.


      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'

!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),ne,NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),ne_start,NKJ(NJM,NPM),
     '  NKJE(NKM,NNM,NJM,NEM),noelem,nonode,np,NP_INTERFACE(0:NPM,0:3),
     &  np_start,NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM)
      LOGICAL MAKE
      CHARACTER ERROR*(*)
!     Local variables
      CHARACTER STRING*255

      CALL ENTERS('GN1DNE',*9999)

      IF(MAKE)THEN
        ne=ne+1
        noelem=noelem+1
        CALL ASSERT(noelem.LE.NE_R_M,'>>Increase NE_R_M',ERROR,*9999)
        CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
        NEELEM(noelem,nr)=ne        
        NPNE(1,nb,ne)=np_start
        NENP(np_start,0,nr)=NENP(np_start,0,nr)+1
        WRITE(STRING,'(''>>Increase NEPM to at least  '',I6)')
     &    NENP(np_start,0,nr)
        CALL ASSERT(NENP(np_start,0,nr).LE.NEPM,STRING,ERROR,*9999)
        NENP(np_start,NENP(np_start,0,nr),nr)=ne
        IF(DOP)THEN
          WRITE(OP_STRING,'('' NPNE(1) '',I6,'' NENP(0,np_start)'//
     '      ''',I6,'' NENP(n,np_start) '',I6)') np_start,NENP(np_start,
     '      0,nr),NENP(np_start,NENP(np_start,0,nr),nr)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        np=np+1
        nonode=nonode+1
        CALL ASSERT(nonode.LE.NP_R_M,'>>Increase NP_R_M',ERROR,*9999)
        CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
        NPNODE(nonode,nr)=np !create a new node
        NENP(np,0,nr)=0 !initialise
        NPNE(2,nb,ne)=np !end node of new element
        NENP(np,0,nr)=0 !initialise
        NENP(np,0,nr)=NENP(np,0,nr)+1
        NENP(np,NENP(np,0,nr),nr)=ne
        NP_INTERFACE(np,0)=1
        NP_INTERFACE(np,1)=nr
        IF(DOP)THEN
          WRITE(OP_STRING,'('' NPNE(2) '',I6,'' NENP(0,np)'//
     '      ''',I6)') np,ne
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        NXI(1,0,ne)=0 !initialise number of proximal branches
        IF(ne_start.NE.0)THEN
c          NXI(-1,0,ne)=NXI(-1,0,ne)+1
          NXI(-1,0,ne)=1
          NXI(-1,NXI(-1,0,ne),ne)=ne_start
          NXI(1,0,ne_start)=NXI(1,0,ne_start)+1
          IF(DOP)THEN
            WRITE(OP_STRING,'('' ne_start '',I6,'' number'//
     '        ''',I6)') ne_start,NXI(1,0,ne_start)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          NXI(1,NXI(1,0,ne_start),ne_start)=ne
        ENDIF
      ENDIF
      IF(DOP)THEN
        WRITE(OP_STRING,'('' NXI(-1) '',I6,'' NXI(1)'//
     '    ''',I6)') ne_start,ne
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,np,np_start,nr,NRE,NVJE,NVJP,
     '  SE,ERROR,*9999)

      CALL EXITS('GN1DNE')
      RETURN
 9999 CALL ERRORS('GN1DNE',ERROR)
      CALL EXITS('GN1DNE')
      RETURN 1
      END


