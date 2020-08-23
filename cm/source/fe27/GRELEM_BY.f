      SUBROUTINE GRELEM_BY(GROUP_TYPE,NELIST,NORD,ne_parent,NXI,
     &  STRING2,AS,ERROR,*)

C#### Subroutine: GRELEM_BY
C###  Description:
C###    GRELEM_BY groups elements by generation, order, lobe, or terminal.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER GROUP_TYPE,NELIST(0:NEM),NELIST2(0:NEM),NORD(5,NE_R_M),
     &  ne_parent,NXI(-NIM:NIM,0:NEIM,0:NEM)
      LOGICAL AS
      CHARACTER STRING2*(255),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,M,N,ne,ne_count,negen,NE_OLD(NE_R_M),
     &  NE_TEMP(NE_R_M),ne0,ngen,NGEN_MAX,noelem,NT_BNS,NUM_NODES
      CHARACTER CHAR2*4,STRINGT*(255)

      CALL ENTERS('GRELEM_BY',*9999)

      NELIST2(0)=0
      IF(GROUP_TYPE.EQ.1)THEN !BY GENERATION
        ne=NELIST(NELIST(0))
        NGEN_MAX=NORD(1,ne)!assuming last element=highest generation
        DO ngen=1,NGEN_MAX
          NELIST2(0)=0
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            negen=NORD(1,ne)
            IF(negen.EQ.ngen)THEN
              NELIST2(0)=NELIST2(0)+1
              NELIST2(NELIST2(0))=ne
            ENDIF
          ENDDO !noelem
          CALL STRING_TRIM(STRING2,IBEG,IEND)
          IF(AS)THEN
            STRINGT=STRING2(IBEG:IEND)
            CALL APPENDI(IEND,ngen,STRINGT)
          ELSE
            WRITE(CHAR2,'(I4)') ngen
            STRINGT='gen'//CHAR2(IBEG:IEND)
          ENDIF
          CALL GRELEM_SUB(NELIST2,STRINGT,.TRUE.,ERROR,*9999)
        ENDDO !ngen

      ELSE IF(GROUP_TYPE.EQ.2)THEN !BY HORSFIELD ORDER
        ne=NELIST(1)
        NGEN_MAX=NORD(2,ne)!assuming first element=highest order
        DO ngen=1,NGEN_MAX
          NELIST2(0)=0
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            negen=NORD(2,ne)
            IF(negen.EQ.ngen)THEN
              NELIST2(0)=NELIST2(0)+1
              NELIST2(NELIST2(0))=ne
            ENDIF
          ENDDO !noelem
          CALL STRING_TRIM(STRING2,IBEG,IEND)
          IF(AS)THEN
            STRINGT=STRING2(IBEG:IEND)
            CALL APPENDI(IEND,ngen,STRINGT)
          ELSE
            WRITE(CHAR2,'(I4)') ngen
            STRINGT='ord'//CHAR2(IBEG:IEND)
          ENDIF
          CALL GRELEM_SUB(NELIST2,STRINGT,.TRUE.,ERROR,*9999)
        ENDDO !ngen

      ELSE IF(GROUP_TYPE.EQ.4)THEN !BY TERMINAL BRANCHES
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          IF(NXI(1,0,ne).EQ.0)THEN
            NELIST2(0)=NELIST2(0)+1
            NELIST2(NELIST2(0))=ne
          ENDIF !NXI
        ENDDO !noelem
        IF(.NOT.AS) STRING2='terminal'
        CALL GRELEM_SUB(NELIST2,STRING2,.TRUE.,ERROR,*9999)

      ELSE IF(GROUP_TYPE.EQ.5)THEN !BY PARENT LIST
        NT_BNS=1
        NE_OLD(1)=ne_parent
        ne_count=1
        NELIST2(ne_count)=ne_parent
        DO WHILE(NT_BNS.NE.0)
          NUM_NODES=NT_BNS
          NT_BNS=0
          DO M=1,NUM_NODES
            ne0=NE_OLD(M) !parent global element #
            DO N=1,NXI(1,0,ne0) !for each daughter branch
              NT_BNS=NT_BNS+1
              NE_TEMP(NT_BNS)=NXI(1,N,ne0)
            ENDDO !N
          ENDDO !M
          DO N=1,NT_BNS
            NE_OLD(N)=NE_TEMP(N) !updates list of prev gen elem#s
            ne_count=ne_count+1
            NELIST2(ne_count)=NE_TEMP(N)
          ENDDO !N
        ENDDO !WHILE
        NELIST2(0)=ne_count
        CALL STRING_TRIM(STRING2,IBEG,IEND)
        IF(AS)THEN
          STRINGT=STRING2(IBEG:IEND)
          CALL APPENDI(IEND,ne_parent,STRINGT)
        ELSE
          STRING2='parent'
          CALL STRING_TRIM(STRING2,IBEG,IEND)
          STRINGT=STRING2(IBEG:IEND)
          CALL APPENDI(IEND,ne_parent,STRINGT)
        ENDIF
        CALL GRELEM_SUB(NELIST2,STRINGT,.TRUE.,ERROR,*9999)

      ENDIF !GROUP_TYPE

      CALL EXITS('GRELEM_BY')
      RETURN
 9999 CALL ERRORS('GRELEM_BY',ERROR)
      CALL EXITS('GRELEM_BY')
      RETURN 1
      END


