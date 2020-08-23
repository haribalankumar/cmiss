      SUBROUTINE OPEG50(COORDS,IPOINTTYP,NBH,ng,nr,nx,
     '  DZDX,EG,PHI,PST,R,RI1,RI2,RI3,RM,U,XG,XI,ZG,FULL,
     '  HIGH_PRECISION,ERROR,*)

C#### Subroutine: OPEG50
C###  Description:
C###    OPEG50 outputs strain fields at current solution wrt 'Reference'
C###    or 'Fibre' coordinates as specified in COORDS.

C**** RI1,RI2,RI3 are principal invariants of AZL
C**** DZDX are components of the deformation gradient tensor
C**** XG   are undeformed theta coords and derivs wrt Xi
C**** ZG   are deformed theta coords and derivs wrt undeformed COORDS
C**** EG   are physical cpts of Green's strain
C**** EXR  are extension ratios wrt COORDS
C**** PHI  are the Euler angles wrt COORDS of the principal extensions
C**** PST  are principal strains
C**** R    is the orthogonal rotation tensor
C**** RM   is the modal matrix whose cols are the eigenvectors
C**** U    is the right stretch tensor

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'stra00.cmn'
!     Parameter List
      INTEGER IPOINTTYP,NBH(NHM),ng,nr,nx
      REAL*8 DZDX(3,3),EG(3,3),PHI(3),PST(3),R(3,3),RI1,RI2,RI3,
     '  RM(3,3),U(3,3),XG(NJM,NUM),XI(3),ZG(NHM,NUM)
      CHARACTER COORDS*(*),ERROR*(*)
      LOGICAL FULL,HIGH_PRECISION
!     Local Variables
      INTEGER i0,i1,IBEG,IEND,IOSTAT,LENI,LENR,LENS,mi,
     '  nb,nhx,ni,NITB,njj1,njj2
      REAL*8 EXR,RI
      PARAMETER (LENI=4,LENR=12,LENS=1)
      CHARACTER CHAR2*1,
     '  CI*(LENI),CR*(LENI),FMTI*(LENI+3),FMTR*(LENI+3),POINTLABEL*3

      CALL ENTERS('OPEG50',*9999)
      nb=NBH(NH_LOC(1,nx))
      NITB=NIT(nb)
      IF(TABLE) THEN
        WRITE(CI,'(I2)') LENI
        WRITE(CR,'(I2)') LENR
        FMTI='(I'//CI//')'
C Chnaging to fixed format as no cfromr call now
C        FMTR='(*'//CR//')'
        FMTR='(F10.4)'
        FORMAT=' '
        i0=1
        i1=LENI
        WRITE(FORMAT(i0:i1),FMTI) ng
        DO ni=1,NITB
          DO mi=1,ni
            i0=i1+LENS+1
            i1=i0+LENR-1
C            FORMAT(i0:i1)=CFROMR(EG(mi,ni),FMTR)
            WRITE(FORMAT(i0:i1),FMTR) EG(mi,ni)
          ENDDO
        ENDDO
        DO ni=1,NITB
          i0=i1+LENS+1
          i1=i0+LENR-1
C          FORMAT(i0:i1)=CFROMR(PST(ni),FMTR)
          WRITE(FORMAT(i0:i1),FMTR) PST(ni)
        ENDDO
        DO ni=1,NITB
          i0=i1+LENS+1
          i1=i0+LENR-1
C          FORMAT(i0:i1)=CFROMR(PHI(ni),FMTR)
          WRITE(FORMAT(i0:i1),FMTR) PHI(ni)
        ENDDO
        WRITE(OP_STRING,FMT='('''//FORMAT//''')',IOSTAT=IOSTAT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(IOSTAT.NE.0) THEN
          ERROR=' Write error'
          GOTO 9999
        ENDIF

      ELSE
C MPN 25May2000: various changes below to tidy up format of output
        IF(IPOINTTYP.EQ.1) THEN
          POINTLABEL='ng='
        ELSE IF(IPOINTTYP.EQ.2) THEN
          POINTLABEL='nd='
        ELSE
          POINTLABEL='pt.'
        ENDIF
        FORMAT='(/'' '//POINTLABEL//''',I9,''  X(j,0):'',6D12.4)'
        WRITE(OP_STRING,FMT=FORMAT) ng,((XG(NJ_LOC(njj1,njj2,nr),1),
     '    njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,2)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( '' ======'', '' Z(j,0):'',6D12.4)'
        WRITE(OP_STRING,FMT=FORMAT) (ZG(NH_LOC(nhx,nx),1),
     '    nhx=1,NH_LOC(0,nx))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( ''       '', ''  Xi(i):'',3D12.4)'
        WRITE(OP_STRING,FMT=FORMAT) (XI(ni),ni=1,NITB)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        CALL STRING_TRIM(COORDS,IBEG,IEND)
        FORMAT='(7X,''Physical Greens strains and extension ratios'
     '    //' wrt '//COORDS(IBEG:IEND)//' coords, and principal'
     '    //' invariants:'')'
        WRITE(OP_STRING,FMT=FORMAT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO mi=1,NITB
          EXR=DSQRT(2.0d0*EG(mi,mi)+1.0d0)
          IF(mi.EQ.1) THEN
C TVK 16Nov1999: option to raise output precision
            IF(HIGH_PRECISION) THEN
              FORMAT='(7X,''EG(1,j):'', D18.10,37X,''EXR(1):'','
     '        //'D18.10,''  RI1:'',D18.10)'
            ELSE
              FORMAT='(7X,''EG(1,j):'', D12.4, 25X,''EXR(1):'','
     '        //'D12.4, ''  RI1:'',D12.4)'
            ENDIF
            RI=RI1
          ELSE IF(mi.EQ.2) THEN
C TVK 16Nov1999: option to raise output precision
            IF(HIGH_PRECISION) THEN
              FORMAT='(7X,''   2    '',2D18.10,19X,''    2  '','
     '        //'D18.10,''  RI2:'',D18.10)'
            ELSE
              FORMAT='(7X,''   2    '',2D12.4, 13X,''    2  '','
     '        //'D12.4, ''  RI2:'',D12.4)'
            ENDIF
            RI=RI2
          ELSE IF(mi.EQ.3) THEN
C TVK 16Nov1999: option to raise output precision
            IF(HIGH_PRECISION) THEN
              FORMAT='(7X,''   3    '',3D18.10, 1X,''    3  '','
     '        //'D18.10,''  RI3:'',D18.10)'
            ELSE
              FORMAT='(7X,''   3    '',3D12.4,  1X,''    3  '','
     '        //'D12.4, ''  RI3:'',D12.4)'
            ENDIF
            RI=RI3
          ENDIF
          WRITE(OP_STRING,FMT=FORMAT) (EG(mi,ni),NI=1,mi),EXR,RI
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

        IF(FULL) THEN

        CALL STRING_TRIM(COORDS,IBEG,IEND)
        FORMAT='(7X,''Deformation gradient, rotation & right stretch'//
     '    ' tensors wrt '//COORDS(IBEG:IEND)//' coords:'')'
        WRITE(OP_STRING,FMT=FORMAT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(CHAR2,'(I1)') NITB
        DO mi=1,NITB
          IF(mi.EQ.1) THEN
            FORMAT='(7X,''F(1,i): '','//CHAR2//'D12.4,'' R(1,i):'','
     '        //CHAR2//'D12.4,'' U(1,i):'',D12.4)'
          ELSE IF(MI.EQ.2) THEN
            FORMAT='(7X,''  2     '','//CHAR2//'D12.4,''   2    '','
     '        //CHAR2//'D12.4,''   2    '',2D12.4)'
          ELSE IF(mi.EQ.3) THEN
            FORMAT='(7X,''  3     '','//CHAR2//'D12.4,''   3    '','
     '        //CHAR2//'D12.4,''   3    '',3D12.4)'
          ENDIF
          WRITE(OP_STRING,FMT=FORMAT) (DZDX(mi,ni),ni=1,NITB),
     '      (R(mi,ni),ni=1,NITB),(U(mi,ni),ni=1,mi)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        CALL STRING_TRIM(COORDS,IBEG,IEND)
        FORMAT='(7X,''Principal strains, Euler angles (deg)'//
     '    ' and eigenvectors wrt '//COORDS(IBEG:IEND)//' coords:'')'
        WRITE(OP_STRING,FMT=FORMAT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO mi=1,NITB
          IF(mi.EQ.1) THEN
            FORMAT='(7X,''PST(1): '',D12.4,''  PHI(1):'','//
     '        'D12.4,''   RM(1,i):'',3D12.4)'
          ELSE IF(mi.EQ.2) THEN
            FORMAT='(7X,''    2   '',D12.4,''      2  '','//
     '        'D12.4,''      2    '',3D12.4)'
          ELSE IF(mi.EQ.3) THEN
            FORMAT='(7X,''    3   '',D12.4,''      3  '','//
     '        'D12.4,''      3    '',3D12.4)'
          ENDIF
          WRITE(OP_STRING,FMT=FORMAT) PST(MI),PHI(MI),
     '      (RM(mi,ni),ni=1,NITB)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

        ENDIF

      ENDIF

      CALL EXITS('OPEG50')
      RETURN
 9999 CALL ERRORS('OPEG50',ERROR)
      CALL EXITS('OPEG50')
      RETURN 1
      END


