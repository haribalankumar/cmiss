      SUBROUTINE IPMESH1(NBJ,NEELEM,NENP,NKJE,NKJ,NP_INTERFACE,NPNE,
     &  NPNODE,nr,NRE,NVJE,NVJP,SE,XP,CALCU,ERROR,*)

C#### Subroutine: IPMESH1
C###  Description:
C###    IPMESH1 defines mesh parameters for regular mesh.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh01.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NP_INTERFACE(0:NPM,0:3),nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJP(NJM,NPM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      LOGICAL CALCU
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IBEG2,IBEG3,IEND2,IEND3,ICHAR,IFROMC,INFO,M,m1,m2,m3,
     '  N,N12,N3CO,nb,nb1,nc,ne,ne1,NE2,NE3,ne_count,ne_start,ni,nj,nk,
     &  nn,nnp,ns,np,noelem,nonode,NOQUES,np_count,np_start
      REAL*8 RDEFAULT(3,8),RFROMC,SIDE(3,3),SIDE_COARSE(3,3),
     &  SIDE_FINE(3,3),SIDE_LENGTH(3),SIDE_TOT(3,3),SIDE_X(3)
      CHARACTER CDEFAULT1(2)*1,CDEFAULT2(4)*3,CDEFAULT3(8)*5,
     '  CHAR1*1,CHAR2*2,CHAR3*5
      LOGICAL CBBREV,FILEIP,SAMETYPE,SAMENUMXI
      DATA RDEFAULT/0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0,
     '              1.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0,
     '              1.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0,
     '              1.0d0, 1.0d0, 1.0d0/
      DATA CDEFAULT1/'0','1'/,
     '     CDEFAULT2/'0,0','1,0','0,1','1,1'/,
     '     CDEFAULT3/'0,0,0','1,0,0','0,1,0','1,1,0',
     '               '0,0,1','1,0,1','0,1,1','1,1,1'/

      CALL ENTERS('IPMESH1',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      ne_start=NET(0)
      ne_count=0
      np_start=NPT(0)
      np_count=np
      
      IF(.NOT.CALCU)THEN
        FORMAT='('' Enter mesh type [1]:'''//
     '    '/''   (1) Even spacing'''//
     '    '/''   (2) Spacing specified by blocks'''//
     '    '/''   (3) Spacing specified by position'''//
     '    '/''   (4) Unused'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=MESH1_TYPE
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MESH1_TYPE=IDATA(1)
        
        FORMAT='($,'' Enter basis function # for mesh [1]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=MESH1_NB
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBFM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MESH1_NB=IDATA(1)
        
        nb=MESH1_NB
        DO ni=1,3
          DO i=0,3
            MESH1_S(ni,i)=0 !initialize #elements in Xi directions
          ENDDO
        ENDDO
        
C Define corner node positions
        WRITE(CHAR1,'(I1)') NJT
        DO N=1,2**NIT(nb) !loops over corner nodes
          WRITE(CHAR2,'(I1)') N
          DO nj=1,NJT
            RDEFLT(nj)=RDEFAULT(nj,N)
          ENDDO
          IF(NJT.EQ.1) THEN
            CHAR3=CDEFAULT1(N)
          ELSE IF(NJT.EQ.2) THEN
            CHAR3=CDEFAULT2(N)
          ELSE IF(NJT.EQ.3) THEN
            CHAR3=CDEFAULT3(N)
          ENDIF
          CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
          FORMAT='($,'' Enter the '//CHAR1//' coords of corner node '
     '      //CHAR2(1:1)//' ['//CHAR3(IBEG3:IEND3)//']: '',3D12.4)'
          IF(IOTYPE.EQ.3) THEN
            DO nj=1,NJT
              RDATA(nj)=MESH1_COORD(N,nj)
            ENDDO
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nj=1,NJT
              MESH1_COORD(N,nj)=RDATA(nj)
            ENDDO
          ENDIF
        ENDDO !n
        
        DO ni=1,NIT(nb)
          WRITE(CHAR1,'(I1)') ni
          IF(MESH1_TYPE.EQ.1) THEN !even spacing
            FORMAT='($,'' Enter #elements in s('//CHAR1//') [1]: '',I4)'
            IF(IOTYPE.EQ.3) IDATA(1)=MESH1_S(ni,1)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) MESH1_S(ni,1)=IDATA(1)
            MESH1_S(ni,0)=MESH1_S(ni,1)
            
          ELSE IF(MESH1_TYPE.EQ.2) THEN !Spacing specified by blocks
            FORMAT='($,'' Enter #s of elements in 3 blocks in s('
     '        //CHAR1//') [1,1,1]: '',I4)'
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=MESH1_S(ni,1)
              IDATA(2)=MESH1_S(ni,2)
              IDATA(3)=MESH1_S(ni,3)
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              MESH1_S(ni,1)=IDATA(1)
              MESH1_S(ni,2)=IDATA(2)
              MESH1_S(ni,3)=IDATA(3)
            ENDIF
            MESH1_S(ni,0)=MESH1_S(ni,1)+MESH1_S(ni,2)+MESH1_S(ni,3)
            
            IDEFLT(1)=3
            FORMAT=
     &        '($,'' Enter (integer) ratio of coarse to fine spacing'//
     &        ' [1]: '',I4)'
            IF(IOTYPE.EQ.3) IDATA(1)=MESH1_R(ni)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) MESH1_R(ni)=IDATA(1)
            
          ELSE IF(MESH1_TYPE.EQ.3) THEN !Spacing specified by position
            FORMAT='($,'' Enter #elements in s('//CHAR1//') [1]: '',I4)'
            IF(IOTYPE.EQ.3) IDATA(1)=MESH1_S(ni,1)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,99,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) MESH1_S(ni,1)=IDATA(1)
            MESH1_S(ni,0)=MESH1_S(ni,1)
            
            MESH1_X(ni,0)=0.0d0 !to set reference point
            MESH1_X(ni,MESH1_S(ni,0))=1.0d0 !final point
            IF(MESH1_S(ni,0)-1.GT.0) THEN
              WRITE(CHAR2,'(I2)') MESH1_S(ni,0)-1
              CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
              FORMAT='($,'' Enter the '//CHAR2(IBEG2:IEND2)
     '          //' relative positions of interior nodes along s('
     &          //CHAR1//') [0]: '',D12.4)'
              IF(IOTYPE.EQ.3) THEN
                DO M=1,MESH1_S(ni,0)-1
                  RDATA(M)=MESH1_X(ni,M)
                ENDDO
              ENDIF
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,MESH1_S(ni,0)-1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &          IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,
     &          1.0d0,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                DO M=1,MESH1_S(ni,0)-1
                  MESH1_X(ni,M)=RDATA(M)
                ENDDO
                MESH1_X(ni,0)=0.0d0 !to set reference point
                MESH1_X(ni,MESH1_S(ni,0))=1.0d0 !final point
              ENDIF
            ENDIF
          ENDIF !mesh1_type
        ENDDO !ni

      ELSE IF(CALCU)THEN
        MESH1_TYPE=0
        IF(CBBREV(CO,'EVEN_SPACING',4,noco+1,NTCO,N3CO))THEN
          MESH1_TYPE=1 !for mesh with even spacing
          IF(CBBREV(CO,'BASIS',3,noco+1,NTCO,N3CO))THEN
            MESH1_NB=IFROMC(CO(N3CO+1))
          ELSE
            MESH1_NB=1
          ENDIF
          nb=MESH1_NB
          DO ni=1,3
            DO i=0,3
              MESH1_S(ni,i)=0 !initialize #elements in Xi directions
            ENDDO !i
          ENDDO !ni
          IF(CBBREV(CO,'START_COORDINATES',4,noco+1,NTCO,N3CO))THEN
            DO nj=1,NJT
              MESH1_COORD(1,nj)=RFROMC(CO(N3CO+nj))
            ENDDO !nj
          ELSE
            DO nj=1,NJT
              MESH1_COORD(1,nj)=0.d0
            ENDDO !nj
          ENDIF
          DO N=1,2**NJT
            DO nj=1,NJT
              MESH1_COORD(N,nj)=MESH1_COORD(1,nj)
            ENDDO !nj
          ENDDO !N
          IF(CBBREV(CO,'END_COORDINATES',3,noco+1,NTCO,N3CO))THEN
            DO nj=1,NJT
              MESH1_COORD(2**NJT,nj)=RFROMC(CO(N3CO+nj))
            ENDDO !nj
          ELSE
            DO nj=1,NJT
              MESH1_COORD(2**NJT,nj)=0.d0
            ENDDO !nj
          ENDIF
C... To tidy up!          
          MESH1_COORD(2,1)=MESH1_COORD(2**NJT,1)
          MESH1_COORD(3,2)=MESH1_COORD(2**NJT,2)
          MESH1_COORD(4,1)=MESH1_COORD(2**NJT,1)
          MESH1_COORD(4,2)=MESH1_COORD(2**NJT,2)
          MESH1_COORD(5,3)=MESH1_COORD(2**NJT,3)
          MESH1_COORD(6,1)=MESH1_COORD(2**NJT,1)
          MESH1_COORD(6,3)=MESH1_COORD(2**NJT,3)
          MESH1_COORD(7,2)=MESH1_COORD(2**NJT,2)
          MESH1_COORD(7,3)=MESH1_COORD(2**NJT,3)

          IF(CBBREV(CO,'SPACING',3,noco+1,NTCO,N3CO))THEN
            DO ni=1,NIT(nb)
              MESH1_S(ni,1)=IFROMC(CO(N3CO+ni))
              MESH1_S(ni,0)=MESH1_S(ni,1)
            ENDDO !ni
          ENDIF
          
        ELSE
          CALL ASSERT(MESH1_TYPE.GT.0,'Not implemented, use file',
     &      ERROR,*9999)
        ENDIF
        
      ENDIF !CALCU
      
!     Calculate lengths of elements along side ni
      DO ni=1,NIT(nb)
        nn=2+ni*(ni-1)/2 !node position# (2,3 or 5) for end of side ni

        DO nj=1,NJT
          !Compute x(nj) distance along a side
          SIDE_X(nj)=MESH1_COORD(nn,nj)-MESH1_COORD(1,nj)

          !Compute fine element length and coarse element length
          IF(MESH1_TYPE.EQ.1) THEN      !even spacing
            SIDE_COARSE(ni,nj)=SIDE_X(nj)/MESH1_S(ni,0)
            SIDE_FINE(ni,nj)  =SIDE_COARSE(ni,nj)
          ELSE IF(MESH1_TYPE.EQ.2) THEN !spacing specified by blocks
            SIDE_FINE(ni,nj)=SIDE_X(nj)
     '        /((MESH1_S(ni,1)+MESH1_S(ni,3))*MESH1_R(ni)+MESH1_S(ni,2))
            SIDE_COARSE(ni,nj)=MESH1_R(ni)*SIDE_FINE(ni,nj)
          ELSE IF(MESH1_TYPE.EQ.3) THEN !spacing specified by position
          ENDIF
        ENDDO !nj

        !Compute total length of side ni
        SIDE_LENGTH(ni)=0.0D0
        DO nj=1,NJT
          SIDE_LENGTH(ni)=SIDE_LENGTH(ni)+SIDE_X(nj)**2
        ENDDO
        SIDE_LENGTH(ni)=DSQRT(SIDE_LENGTH(ni))
        IF(DOP) THEN
          WRITE(OP_STRING,'('' SIDE_LENGTH('',I1,''):'',D12.4)')
     '      ni,SIDE_LENGTH(ni)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO !ni

      NE1=MESH1_S(1,0) !is total #elements in s(1) direction
      NE2=MESH1_S(2,0) !is total #elements in s(2) direction
      NE3=MESH1_S(3,0) !is total #elements in s(3) direction

!Loop over the nodes and elements
      DO nj=1,NJT
        SIDE_TOT(3,nj)=0.0D0
      ENDDO
!  Step in Xi(3) direction
      DO m3=1,NE3+1
        DO nj=1,NJT
          IF(MESH1_TYPE.LE.2) THEN !spacing even or specified by blocks
            IF(m3.LE.MESH1_S(3,1)) THEN                     !1st block
              SIDE(3,nj)=SIDE_COARSE(3,nj)
            ELSE IF(m3.LE.MESH1_S(3,1)+MESH1_S(3,2)) THEN   !2nd block
              SIDE(3,nj)=SIDE_FINE(3,nj)
            ELSE                                            !3rd block
              SIDE(3,nj)=SIDE_COARSE(3,nj)
            ENDIF
          ELSE IF(MESH1_TYPE.EQ.3) THEN !spacing specified by position
            SIDE(3,nj)=(MESH1_X(3,m3)-MESH1_X(3,m3-1))
     '                *(MESH1_COORD(5,nj)-MESH1_COORD(1,nj))
            IF(DOP) THEN
              WRITE(OP_STRING,'('' SIDE(3,nj='',I1,''):'',D12.4)')
     '          nj,SIDE(3,nj)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO !nj

        DO nj=1,NJT
          SIDE_TOT(2,nj)=0.0D0
        ENDDO
!  Step in Xi(2) direction
        DO m2=1,NE2+1
          DO nj=1,NJT
            IF(MESH1_TYPE.LE.2) THEN !spacing even or by blocks
              IF(m2.LE.MESH1_S(2,1)) THEN                   !1st block
                SIDE(2,nj)=SIDE_COARSE(2,nj)
              ELSE IF(m2.LE.MESH1_S(2,1)+MESH1_S(2,2)) THEN !2nd block
                SIDE(2,nj)=SIDE_FINE(2,nj)
              ELSE                                          !3rd block
                SIDE(2,nj)=SIDE_COARSE(2,nj)
              ENDIF
            ELSE IF(MESH1_TYPE.EQ.3) THEN !spacing by position
              SIDE(2,nj)=(MESH1_X(2,m2)-MESH1_X(2,m2-1))
     '                  *(MESH1_COORD(3,nj)-MESH1_COORD(1,nj))
              IF(DOP) THEN
                WRITE(OP_STRING,'('' SIDE(2,nj='',I1,''):'',D12.4)')
     '            nj,SIDE(2,nj)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
          ENDDO !nj

          DO nj=1,NJT
            SIDE_TOT(1,nj)=0.0D0
          ENDDO
!  Step in Xi(1) direction
          DO m1=1,NE1+1
            DO nj=1,NJT
              IF(MESH1_TYPE.LE.2) THEN !spacing even or by blocks
                IF(m1.LE.MESH1_S(1,1)) THEN                 !1st block
                  SIDE(1,nj)=SIDE_COARSE(1,nj)
                ELSE IF(m1.LE.MESH1_S(1,1)+MESH1_S(1,2))THEN!2nd block
                  SIDE(1,nj)=SIDE_FINE(1,nj)
                ELSE                                        !3rd block
                  SIDE(1,nj)=SIDE_COARSE(1,nj)
                ENDIF
              ELSE IF(MESH1_TYPE.EQ.3) THEN !spacing by position
                SIDE(1,nj)=(MESH1_X(1,m1)-MESH1_X(1,m1-1))
     '                    *(MESH1_COORD(2,nj)-MESH1_COORD(1,nj))
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' SIDE(1,nj='',I1,'
     '              //'''):'',D12.4)') nj,SIDE(1,nj)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO !nj

C MHT 28/6/05 np=NPT(0)+..., so that increases from current
            np=NPT(0)+m1+(m2-1)*(NE1+1)+(m3-1)*(NE1+1)*(NE2+1)
            np_count=np_count+1
            CALL ASSERT(NPM.GE.np,'>>Increase NPM',ERROR,*9999)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)

            DO nj=1,NJT
              XP(1,1,nj,np)=MESH1_COORD(1,nj)
     '          +SIDE_TOT(1,nj)+SIDE_TOT(2,nj)+SIDE_TOT(3,nj)
            ENDDO

            IF(DOP) THEN
              WRITE(OP_STRING,'('' m3='',I3,'' m2='',I3,'
     '          //''' m1='',I3)') m3,m2,m1
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' np='',I5,'' xp:'',3E12.3)')
     '          np,(XP(1,1,nj,np),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(NIT(nb).LT.3.OR.m3.LE.NE3) THEN
              IF(NIT(nb).LT.2.OR.m2.LE.NE2) THEN
                IF(m1.LE.NE1) THEN
C MHT 28/6/05 ne=NET(0)+... so that increase from current
                  ne=NET(0)+m1+(m2-1)*NE1+(m3-1)*NE1*NE2
                  ne_count=ne_count+1
C LKC 7-NOV-97 added assert
                  CALL ASSERT(NEM.GE.ne,'>>Increase NEM',ERROR,*9999)
                  IF(DOP) write(*,'('' ne='',i4)') ne
                  NPNE(1,nb,ne)=np             !1st node of element
                  NPNE(2,nb,ne)=np+1           !2nd node of element
                  IF(NIT(nb).GT.1) THEN
                    NPNE(3,nb,ne)=np+NE1+1     !3rd node of element
                    NPNE(4,nb,ne)=np+NE1+2     !4th node of element
                  ENDIF
                  IF(NIT(nb).GT.2) THEN
                    N12=(NE1+1)*(NE2+1)
                    NPNE(5,nb,ne)=np+N12       !5th node of element
                    NPNE(6,nb,ne)=np+N12+1     !6th node of element
                    NPNE(7,nb,ne)=np+N12+NE1+1 !7th node of element
                    NPNE(8,nb,ne)=np+N12+NE1+2 !8th node of element
                  ENDIF
                ENDIF
              ENDIF
            ENDIF !nit

            IF(m1.LE.NE1) THEN
              DO nj=1,NJT
                SIDE_TOT(1,nj)=SIDE_TOT(1,nj)+SIDE(1,nj)
              ENDDO
            ENDIF
          ENDDO !m1

          IF(m2.LE.NE2) THEN
            DO nj=1,NJT
              SIDE_TOT(2,nj)=SIDE_TOT(2,nj)+SIDE(2,nj)
            ENDDO
          ENDIF
        ENDDO !m2

        IF(m3.LE.NE3) THEN
          DO nj=1,NJT
            SIDE_TOT(3,nj)=SIDE_TOT(3,nj)+SIDE(3,nj)
          ENDDO
        ENDIF
      ENDDO !m3

      NET(nr)=ne      !highest element# in region nr
      NET(0) =ne      !highest element# in all regions
      NEELEM(0,nr)=ne_count !number elements in region nr
c      NEELEM(0,nr)=ne !number elements in region nr
      NEELEM(0,0) =NEELEM(0,0)+ne_count !number elements in all regions
c      NEELEM(0,0) =ne !number elements in all regions
      NPT(nr)=np      !highest node# in region nr
      NPT(0) =np      !highest node# in all regions
      NPNODE(0,nr)=np_count !number nodes in region nr
c      NPNODE(0,nr)=np !number nodes in region nr
      NPNODE(0,0) =NPNODE(0,0)+np_count !number nodes in all regions
c      NPNODE(0,0) =np !number nodes in all regions

      DO noelem=1,NEELEM(0,nr)
        ne=ne_start+noelem
c        ne=noelem
        NEELEM(noelem,nr)=ne
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          NBJ(nj,ne)=nb
          DO nn=1,NNT(nb)
            DO nk=1,NKT(nn,nb)
              NKJE(nk,nn,nj,ne)=nk
            ENDDO !nk
          ENDDO !nn
        ENDDO
        NRE(ne)=nr
      ENDDO

      DO nonode=1,NPNODE(0,nr)
        np=np_start+nonode
c        np=nonode
        NPNODE(nonode,nr)=np
        DO nj=1,NJT
          NKJ(nj,np)=NKT(0,nb)
          DO nc=1,NCM
            !one version per nj per node
            NVJP(nj,np)=1
          ENDDO !nc
        ENDDO !nj
      ENDDO !nonode (np)

      DO nb1=1,NBFT
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO ns=1,NST(nb1)+NAT(nb1)
            SE(ns,nb1,ne)=1.0D0
          ENDDO !ns
          DO nn=1,NNT(nb1)
C KAT 23Feb01: now handled by NKJE above
C            DO nk=1,NKT(nn,nb1)
C              NKE(nk,nn,nb1,ne)=NK
C            ENDDO !nk
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              NVJE(nn,nb1,nj,ne)=1 !version one of nn,nj in elem ne
            ENDDO !nj
          ENDDO !nn
C Update alternate bases
          SAMETYPE=.FALSE.
          SAMENUMXI=.FALSE.
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            IF(NBC(nb1).EQ.NBC(NBJ(nj,ne)).OR.NBC(nb1).EQ.7)
     '        SAMETYPE=.TRUE. !Same basis type or extended basis
            IF(NIT(nb1).EQ.NIT(NBJ(nj,ne))) SAMENUMXI=.TRUE.
          ENDDO
          IF(NNT(nb1).GT.0.AND.SAMETYPE.AND.SAMENUMXI.AND.
     '      nb1.NE.nb) THEN
            DO nnp=1,8
              NPNE(nnp,nb1,ne)=NPNE(nnp,nb,ne)
            ENDDO
          ENDIF
        ENDDO !noelem (ne)
      ENDDO !nb1

      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

      CALL EXITS('IPMESH1')
      RETURN
 9999 CALL ERRORS('IPMESH1',ERROR)
      CALL EXITS('IPMESH1')
      RETURN 1
      END


