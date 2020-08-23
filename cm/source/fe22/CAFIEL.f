      SUBROUTINE CAFIEL(ISEG,ISFIEL,NEELEM,NRLIST,STRING,ERROR,*)

C#### Subroutine: CAFIEL
C###  Description:
C###    CAFIEL cancels field segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ISEG(*),ISFIEL(NWM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NRLIST(0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),ne,nj,njj,njj1,njj2,noelem,noiw,
     '  no_nrlist,nr,NTIW
      LOGICAL ABBREV,ALL_REGIONS,SEGME

      CALL ENTERS('CAFIEL',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel field<;s>
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to cancel the
C###    field on.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Description:
C###    Cancel fields. If the s qualified is present the field segments
C###    are canceled on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //';s'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAFIEL',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        CALL CHECKQ(' S',noco,1,CO,COQU,STRING,*1)
        IF(ABBREV(COQU(noco,1),'S',1)) THEN
          SEGME=.TRUE.
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        ELSE
          SEGME=.FALSE.
        ENDIF

        IF(SEGME) THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO nr=1,NRT
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IF(ISFIEL(iw,ne).GT.0) THEN
                    CALL DELETE_SEGMENT(ISFIEL(iw,ne),ISEG,iw,ERROR,
     '                *9999)
                  ENDIF
                ENDDO
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO

        ELSE IF(.NOT.SEGME) THEN
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
              nj=NJ_LOC(NJL_FIEL,njj,nr)
              NJ_TYPE(nj,1)=0
              NJ_TYPE(nj,2)=0
              NJ_LOC(NJL_FIEL,njj,nr)=0
            ENDDO
            NJ_LOC(NJL_FIEL,0,nr)=0
            DO njj1=1,3
              DO njj2=1,NJ_LOC(njj1,0,nr)
                IF(NJ_LOC(njj1,njj2,nr).GT.NJ_LOC(0,0,nr))
     '            NJ_LOC(0,0,nr)=NJ_LOC(njj1,njj2,nr)
              ENDDO !njj2
            ENDDO !njj1
          ENDDO
          NJ_LOC(0,0,0)=0
          NJ_LOC(NJL_FIEL,0,0)=0
          DO nr=1,NRT
            IF(NJ_LOC(0,0,nr).GT.NJ_LOC(0,0,0))
     '        NJ_LOC(0,0,0)=NJ_LOC(0,0,nr)
            IF(NJ_LOC(NJL_FIEL,0,nr).GT.NJ_LOC(NJL_FIEL,0,0))
     '        NJ_LOC(NJL_FIEL,0,0)=NJ_LOC(NJL_FIEL,0,nr)
          ENDDO !nr
          CALL_DATA_FIELD=.FALSE.
        ENDIF
      ENDIF

      CALL EXITS('CAFIEL')
      RETURN
 9999 CALL ERRORS('CAFIEL',ERROR)
      CALL EXITS('CAFIEL')
      RETURN 1
      END


