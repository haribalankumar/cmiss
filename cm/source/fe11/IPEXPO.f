      SUBROUTINE IPEXPO(ERROR,*)

C#### Subroutine: IPEXPO
C###  Description:
C###    IPEXPO does input for export parameters.

C**** CPB 28/6/95  Adding export signal parameters
C**** AJP 7-4-94   Currently sets up options for local refining to output
C                  information in MAP3D format.
C**** LKC 4-SEP-98 Adding signal output to Cmgui .exgobj file

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cmgui01.cmn'
      INCLUDE 'emap00.cmn'
      INCLUDE 'expo00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,ICHAR,IEND,INFO,nolist,NOQUES,nr,nrow,ncol
      CHARACTER CHAR1*1
      LOGICAL FILEIP

      CALL ENTERS('IPEXPO',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

CC AJPs 191297
      FORMAT='(/'' Specify the export type [1]: '''
     '  //'/''   (1) Signal'''
     '  //'/''   (2) ZCROSSING'''
     '  //'/$,''    '',I1)'
CC AJPe
      IF(IOTYPE.EQ.3) IDATA(1)=DEFEXPO_TYPE
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) DEFEXPO_TYPE=IDATA(1)

      IF(DEFEXPO_TYPE.EQ.1) THEN !Signal

        FORMAT='(/'' Specify the signal export type [1]: '''
     '    //'/''   (1) UNEMAP'''
     '    //'/''   (2) CMGUI'''
     '    //'/''   (3) Data file'''
     '    //'/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=SIGEXPO_TYPE
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SIGEXPO_TYPE=IDATA(1)


C LKC 13-SEP-1999
C*** Add the ability to set the frequency
C***  The overcomes the problem of mapping between solving
C***  in time steps and viewing in real time
        FORMAT='(/$,'' Do you wish to set the'//
     '    ' frequency to map between time steps & real time [N]? '',A)'

        IF(IOTYPE.EQ.3) THEN
          IF(.NOT.SET_FREQUENCY) THEN
            ADATA(1)='N'
          ELSE
            ADATA(1)='Y'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          SET_FREQUENCY=.TRUE.
        ELSE
          SET_FREQUENCY=.FALSE.
        ENDIF

        IF(SET_FREQUENCY) THEN
C LKC/AJP 13-SEPT-1999 adding freq to map between real time & time steps
          FORMAT='($,'' Enter frequency to map between'//
     '      ' time and tstep [1000.0]: '', D12.4)'
          RDEFLT(1)=1000.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=CM_FREQUENCY
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) CM_FREQUENCY=RDATA(1)
        ENDIF


        IF(SIGEXPO_TYPE.EQ.1) THEN !UNEMAP
          FORMAT='(/'' Specify the rig type [1]: '''
     '      //'/''   (1) Sock'''
     '      //'/''   (2) Patch'''
     '      //'/''   (3) Torso'''
     '      //'/''   (4) Mixed'''
     '      //'/''   (5) Unused'''
     '      //'/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) THEN
            IF(EMAP_RIGTYPE(0).EQ.EMAP_SOCK) THEN
              IDATA(1)=1
            ELSE IF(EMAP_RIGTYPE(0).EQ.EMAP_PATCH) THEN
              IDATA(1)=2
            ELSE IF(EMAP_RIGTYPE(0).EQ.EMAP_TORSO) THEN
              IDATA(1)=3
            ELSE IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED) THEN
              IDATA(1)=4
            ENDIF
          ENDIF
          IF(IOTYPE.EQ.3) IDATA(1)=EMAP_RIGTYPE(0)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(IDATA(1).EQ.1) THEN
              EMAP_RIGTYPE(0)=EMAP_SOCK
            ELSE IF(IDATA(1).EQ.2) THEN
              EMAP_RIGTYPE(0)=EMAP_PATCH
            ELSE IF(IDATA(1).EQ.3) THEN
              EMAP_RIGTYPE(0)=EMAP_TORSO
            ELSE IF(IDATA(1).EQ.4) THEN
              EMAP_RIGTYPE(0)=EMAP_MIXED
            ENDIF
          ENDIF

          IF((IDATA(1).EQ.1).OR.(IDATA(1).EQ.3).OR.(IDATA(1).EQ.4)) THEN
            CALL ASSERT(NJT.EQ.3,'>>Use PATCH rig type for NJT.LT.3',
     '        ERROR,*9999)
          ENDIF
          IF(IDATA(1).EQ.2) THEN
            CALL ASSERT(NJT.EQ.2,
     &        '>>Use PATCH rig type only for 2d problem',ERROR,*9999)
          ENDIF

          
          IF(EMAP_RIGTYPE(0).EQ.EMAP_SOCK) THEN
            FORMAT='(/$,'' Enter the focus for the sock [1.0]: '','
     '        //'D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=EMAP_SOCKFOCUS(1)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) EMAP_SOCKFOCUS(1)=RDATA(1)
          ENDIF

          FORMAT='($,'' Enter the rig name [CMISS]: '',A)'
          CDEFLT(1)='CMISS'
          IF(IOTYPE.EQ.3) CDATA(1)=EMAP_RIGNAME
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) EMAP_RIGNAME=CDATA(1)(1:30)

          FORMAT='(/$,'' Enter the number UNEMAP regions [1]: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=EMAP_NUMREGIONS
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,10,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) EMAP_NUMREGIONS=IDATA(1)


C!!!  LKC 3-SEP-1999 adding a region list.
C!!! This will need change when data has region dependence

          FORMAT='($,'' Enter the regions: '',10I4,/(1X,I1))'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) THEN
            DO nolist=1,EMAP_NUMREGIONS
              IDATA(nolist)=EXPORT_NRLIST(nolist)
            ENDDO
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      EMAP_NUMREGIONS,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,1,10,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '      ERROR,*9999)

          IF(IOTYPE.NE.3) THEN
            DO nolist=1,EMAP_NUMREGIONS
              EXPORT_NRLIST(nolist)=IDATA(nolist)
            ENDDO
          ENDIF

C*** There should only be 1 region for ZCROSSING exporting
          IF(DEFEXPO_TYPE.NE.1) THEN
            CALL ASSERT(EMAP_NUMREGIONS.EQ.1,
     '        '>> Should only export 1 region for ZCROSSING'
     '        ,ERROR,*9999)
          ENDIF



          DO nr=1,EMAP_NUMREGIONS

            WRITE(CHAR1,'(I1)') nr
            IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED) THEN !Mixed rig type
              FORMAT='(/'' Specify the rig type for region '//CHAR1
     '          //' [1]: '''
     '          //'/''   (1) Sock'''
     '          //'/''   (2) Patch'''
     '          //'/''   (3) Torso'''
     '          //'/''   (4) Unused'''
     '          //'/$,''    '',I1)'
              IF(IOTYPE.EQ.3) THEN
                IF(EMAP_RIGTYPE(nr).EQ.EMAP_SOCK) THEN
                  IDATA(1)=1
                ELSE IF(EMAP_RIGTYPE(nr).EQ.EMAP_PATCH) THEN
                  IDATA(1)=2
                ELSE IF(EMAP_RIGTYPE(nr).EQ.EMAP_TORSO) THEN
                  IDATA(1)=3
                ENDIF
              ENDIF
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                IF(IDATA(1).EQ.1) THEN
                  EMAP_RIGTYPE(nr)=EMAP_SOCK
                ELSE IF(IDATA(1).EQ.2) THEN
                  EMAP_RIGTYPE(nr)=EMAP_PATCH
                ELSE IF(IDATA(1).EQ.3) THEN
                  EMAP_RIGTYPE(nr)=EMAP_TORSO
                ENDIF
              ENDIF
              IF(EMAP_RIGTYPE(nr).EQ.EMAP_SOCK) THEN
                FORMAT='(/$,'' Enter the focus for the sock [1.0]: '','
     '            //'D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=EMAP_SOCKFOCUS(nr)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,
     '            INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) EMAP_SOCKFOCUS(nr)=RDATA(1)
              ENDIF
            ELSE
              EMAP_RIGTYPE(nr)=EMAP_RIGTYPE(0)
            ENDIF
            FORMAT='(/$,'' Enter the region name [CMISSregion'//CHAR1
     '        //']: '',A)'
            CDEFLT(1)='CMISSregion'//CHAR1
            IF(IOTYPE.EQ.3) CDATA(1)=EMAP_REGIONNAME(nr)
            CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) EMAP_REGIONNAME(nr)=CDATA(1)(1:30)
            FORMAT='(/$,'' Enter the start electrode number for '
     '        //'region '//CHAR1//' [1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=EMAP_STARTELEC(nr)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) EMAP_STARTELEC(nr)=IDATA(1)
            FORMAT='(/$,'' Enter the stop electrode number for '
     '        //'region '//CHAR1//' [1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=EMAP_STOPELEC(nr)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) EMAP_STOPELEC(nr)=IDATA(1)
          ENDDO !nr
C LKC 4-SEP-98 No longer used
C
C        ELSE IF(SIGEXPO_TYPE.EQ.2) THEN !Map3d
C          FORMAT='(/$,'' Enter the number of local subdivisions for '
C     '      //'each element [1]: '',I1)'
C          IF(IOTYPE.EQ.3) IDATA(1)=NT_SUB_DIV
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,10,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) NT_SUB_DIV=IDATA(1)


C*** Export to CMGUI
        ELSEIF(SIGEXPO_TYPE.EQ.2) THEN ! export to CMGUI

          FORMAT='(/$,'' Enter graphical object name [trace]: '')'
          CDEFLT(1)='trace'
          IF(IOTYPE.EQ.3) CDATA(1)=NAME_EXGOBJ(1)
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            CALL STRING_TRIM(CDATA(1),IBEG,IEND)
            NAME_EXGOBJ(1)=CDATA(1)(IBEG:IEND)
          ENDIF

          FORMAT='(/$,'' Enter pointer name [trace_arrow]: '')'
          CDEFLT(1)='trace_arrow'
          IF(IOTYPE.EQ.3) CDATA(1)=NAME_EXGOBJ(2)
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            CALL STRING_TRIM(CDATA(1),IBEG,IEND)
            NAME_EXGOBJ(2)=CDATA(1)(IBEG:IEND)
          ENDIF

          FORMAT='(/$,'' Enter the electrode number to export [1]: '')'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=ELEC_EXGOBJ
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ELEC_EXGOBJ=IDATA(1)

          FORMAT='(/$,'' Enter signal material name '
     '      //' [default]: '')'
          CDEFLT(1)='default'
          IF(IOTYPE.EQ.3) CDATA(1)=MATERIAL_EXGOBJ(1)
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            CALL STRING_TRIM(CDATA(1),IBEG,IEND)
            MATERIAL_EXGOBJ(1)=CDATA(1)(IBEG:IEND)
          ENDIF

          FORMAT='(/$,'' Enter arrow pointer material name '
     '      //' [default]: '')'
          CDEFLT(1)='default'
          IF(IOTYPE.EQ.3) CDATA(1)=MATERIAL_EXGOBJ(2)
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            CALL STRING_TRIM(CDATA(1),IBEG,IEND)
            MATERIAL_EXGOBJ(2)=CDATA(1)(IBEG:IEND)
          ENDIF

C Initialise the transfer matrix
          DO nrow=1,4
            DO ncol=1,4
              IF(nrow.EQ.ncol) THEN
                TRANSFORM_EXGOBJ(nrow,ncol)=1.D0
              ELSE
                TRANSFORM_EXGOBJ(nrow,ncol)=0.D0
              ENDIF
            ENDDO
          ENDDO

          FORMAT='(/$,'' Enter time scaling factor [0.15]: '')'
          RDEFLT(1)=0.15D0
          IF(IOTYPE.EQ.3) RDATA(1)=TRANSFORM_EXGOBJ(1,1)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TRANSFORM_EXGOBJ(1,1)=RDATA(1)

          FORMAT='(/$,'' Enter potential scaling factor [25]: '')'
          RDEFLT(1)=25.D0
          IF(IOTYPE.EQ.3) RDATA(1)=TRANSFORM_EXGOBJ(3,3)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) TRANSFORM_EXGOBJ(3,3)=RDATA(1)

          FORMAT='(/$,'' Enter a displacement for object'
     '      //' [0.0 0.0 0.0]: '',3E12.5)'
          DO nrow=1,3
            RDEFLT(nrow)=0.D0
          ENDDO
          IF(IOTYPE.EQ.3) THEN
            DO nrow=1,3
              RDATA(nrow)=TRANSFORM_EXGOBJ(4,nrow)
            ENDDO
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nrow=1,3
              TRANSFORM_EXGOBJ(4,nrow)=RDATA(nrow)
            ENDDO
          ENDIF

C LKC 7-MAY-1999 new export type
C*** Export to a datafile
        ELSEIF(SIGEXPO_TYPE.EQ.3) THEN ! export to datafile

C!!! LKC 7-MAY-1999 Awaiting data point region dependence
          DO nr=1,1
            WRITE(CHAR1,'(I1)') nr

            IF(IOTYPE.NE.3) EMAP_REGIONNAME(nr)=CDATA(1)(1:30)
            FORMAT='(/$,'' Enter the start electrode number for '
     '        //'region '//CHAR1//' [1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=EMAP_STARTELEC(nr)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) EMAP_STARTELEC(nr)=IDATA(1)
            FORMAT='(/$,'' Enter the stop electrode number for '
     '        //'region '//CHAR1//' [1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=EMAP_STOPELEC(nr)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) EMAP_STOPELEC(nr)=IDATA(1)
          ENDDO
        ELSE
          ERROR='Unknown export signal type'
          GOTO 9999
        ENDIF


C*** Export a zero-crossing signal
      ELSEIF(DEFEXPO_TYPE.EQ.2) THEN !Zero Crossing

        FORMAT='(/'' Specify the signal export type [1]: '''
     '    //'/''   (1) UNEMAP'''
     '    //'/''   (2) Zcrossing ranges to a Datafile'''
     '    //'/''   (3) Datafile'''
     '    //'/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=SIGEXPO_TYPE
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SIGEXPO_TYPE=IDATA(1)

C LKC 2-AUG-2000 option 2 now implemented
CC LKC 9-MAY-1999 added assert and added export type 3
C        CALL ASSERT(SIGEXPO_TYPE.NE.2,
C     '    '>> Zero-crossing export type 2 not implemented',
C     '    ERROR,*9999)


C LKC 1-MAR-1999 adding freq to map between real time & time steps
        FORMAT='(/$,'' Enter frequency to map between'//
     '    ' time and tstep [1000.0]: '', D12.4)'
        RDEFLT(1)=1000.0d0
        IF(IOTYPE.EQ.3) RDATA(1)=CM_FREQUENCY
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) CM_FREQUENCY=RDATA(1)


C*** Export to unemap
        IF(SIGEXPO_TYPE.EQ.1) THEN !UNEMAP

          FORMAT='(/'' Specify the rig type [1]: '''
     '      //'/''   (1) Sock'''
     '      //'/''   (2) Patch'''
     '      //'/''   (3) Torso'''
     '      //'/''   (4) Mixed'''
     '      //'/''   (5) Unused'''
     '      //'/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) THEN
            IF(EMAP_RIGTYPE(0).EQ.EMAP_SOCK) THEN
              IDATA(1)=1
            ELSE IF(EMAP_RIGTYPE(0).EQ.EMAP_PATCH) THEN
              IDATA(1)=2
            ELSE IF(EMAP_RIGTYPE(0).EQ.EMAP_TORSO) THEN
              IDATA(1)=3
            ELSE IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED) THEN
              IDATA(1)=4
            ENDIF
          ENDIF
          IF(IOTYPE.EQ.3) IDATA(1)=EMAP_RIGTYPE(0)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(IDATA(1).EQ.1) THEN
              EMAP_RIGTYPE(0)=EMAP_SOCK
            ELSE IF(IDATA(1).EQ.2) THEN
              EMAP_RIGTYPE(0)=EMAP_PATCH
            ELSE IF(IDATA(1).EQ.3) THEN
              EMAP_RIGTYPE(0)=EMAP_TORSO
            ELSE IF(IDATA(1).EQ.4) THEN
              EMAP_RIGTYPE(0)=EMAP_MIXED
            ENDIF
          ENDIF
C LKC 25-FEB-2003 SOCK, TORSO and MIXED are all 3d rig types          
C LKC 26-AUG-97 add assert statement for 2d problems
C         IF((IDATA(1).EQ.3).OR.(IDATA(1).EQ.4)) THEN
          
          IF(EMAP_RIGTYPE(0).NE.EMAP_PATCH) THEN
            CALL ASSERT(NJT.EQ.3,'>>Use PATCH rig type for 2D',
     '        ERROR,*9999)
          ELSE !patch
            CALL ASSERT(NJT.EQ.2,'>>Can only use PATCH rig type for 2D',
     '        ERROR,*9999)            
          ENDIF

          IF(EMAP_RIGTYPE(0).EQ.EMAP_SOCK) THEN
            FORMAT='(/$,'' Enter the focus for the sock [1.0]: '','
     '        //'D12.4)'
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) RDATA(1)=EMAP_SOCKFOCUS(1)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) EMAP_SOCKFOCUS(1)=RDATA(1)
          ENDIF
          FORMAT='($,'' Enter the rig name [CMISS]: '',A)'
          CDEFLT(1)='CMISS'
          IF(IOTYPE.EQ.3) CDATA(1)=EMAP_RIGNAME
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) EMAP_RIGNAME=CDATA(1)(1:30)

          EMAP_NUMREGIONS=1
          EXPORT_NRLIST(1)=1

          DO nr=1,EMAP_NUMREGIONS
            WRITE(CHAR1,'(I1)') nr
            IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED) THEN !Mixed rig type
              FORMAT='(/'' Specify the rig type for region '//CHAR1
     '          //' [1]: '''
     '          //'/''   (1) Sock'''
     '          //'/''   (2) Patch'''
     '          //'/''   (3) Torso'''
     '          //'/''   (4) Unused'''
     '          //'/$,''    '',I1)'
              IF(IOTYPE.EQ.3) THEN
                IF(EMAP_RIGTYPE(nr).EQ.EMAP_SOCK) THEN
                  IDATA(1)=1
                ELSE IF(EMAP_RIGTYPE(nr).EQ.EMAP_PATCH) THEN
                  IDATA(1)=2
                ELSE IF(EMAP_RIGTYPE(nr).EQ.EMAP_TORSO) THEN
                  IDATA(1)=3
                ENDIF
              ENDIF
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                IF(IDATA(1).EQ.1) THEN
                  EMAP_RIGTYPE(nr)=EMAP_SOCK
                ELSE IF(IDATA(1).EQ.2) THEN
                  EMAP_RIGTYPE(nr)=EMAP_PATCH
                ELSE IF(IDATA(1).EQ.3) THEN
                  EMAP_RIGTYPE(nr)=EMAP_TORSO
                ENDIF
              ENDIF
              IF(EMAP_RIGTYPE(nr).EQ.EMAP_SOCK) THEN
                FORMAT='(/$,'' Enter the focus for the sock [1.0]: '','
     '            //'D12.4)'
                RDEFLT(1)=1.0d0
                IF(IOTYPE.EQ.3) RDATA(1)=EMAP_SOCKFOCUS(nr)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,
     '            INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) EMAP_SOCKFOCUS(nr)=RDATA(1)
              ENDIF
            ELSE
              EMAP_RIGTYPE(nr)=EMAP_RIGTYPE(0)
            ENDIF
            FORMAT='(/$,'' Enter the region name [CMISSregion'//CHAR1
     '        //']: '',A)'
            CDEFLT(1)='CMISSregion'//CHAR1
            IF(IOTYPE.EQ.3) CDATA(1)=EMAP_REGIONNAME(nr)
            CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) EMAP_REGIONNAME(nr)=CDATA(1)(1:30)
          ENDDO !nr


C*** Export to a datafile (but not zcrossing ranges)
        ELSEIF(SIGEXPO_TYPE.EQ.2) THEN !Datafile ranges

C Do nothing

        ELSEIF(SIGEXPO_TYPE.EQ.3) THEN !Datafile

C LKC 7-MAY-1999 new export type
C!!! LKC 7-MAY-1999 Awaiting data point region dependence
          DO nr=1,1
            WRITE(CHAR1,'(I1)') nr

            IF(IOTYPE.NE.3) EMAP_REGIONNAME(nr)=CDATA(1)(1:30)
            FORMAT='(/$,'' Enter the start electrode number for '
     '        //'region '//CHAR1//' [1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=EMAP_STARTELEC(nr)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) EMAP_STARTELEC(nr)=IDATA(1)
            FORMAT='(/$,'' Enter the stop electrode number for '
     '        //'region '//CHAR1//' [1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=EMAP_STOPELEC(nr)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) EMAP_STOPELEC(nr)=IDATA(1)
          ENDDO

        ELSE
          ERROR='>>Zcrossing export type not implemented yet'
          GOTO 9999
        ENDIF
      ENDIF

      CALL_EXPO=.TRUE.

      CALL EXITS('IPEXPO')
      RETURN
 9999 CALL ERRORS('IPEXPO',ERROR)
      CALL EXITS('IPEXPO')
      RETURN 1
      END


