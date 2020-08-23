      SUBROUTINE APINVE(ISIZE_TBH,NHP,NHQ,NKH,NPLIST3,NPLIST4,NPNY,
     '  NQNY,NVHP,NXLIST,NYNP,NYNR,NYQNR,S,STRING,T_BH,T_BH_INV,
     '  WK1_INV,WK4_INV,YP,YQ,YQS,ERROR,*)

C#### Subroutine: APINVE
C###  Description:
C###    APINVE applies the inverse of the transfer matrix T_BH_INV to
C###    a set of nodal values stored in YP or a history file.
C###    If the inverse matrix has not been
C###    explicitly evaluated, but rather the user has specified the
C###    inverse solution is to be constructed by solving the
C###    underdetermined matrix equations, then that solution is
C###    performed here.
C**** List of nodes to which T_bh_inv is to be applied are stored in
C**** NPLIST4, results are for nodes in NPLIST3.  Define transfer
C**** should have been called to set these up.
C**** If no history file is specified then T_BH_INV is applied to
C**** values stored in YP - new calculated
C**** values are also stored in YP (old values at this location
C**** are stored in YP(ny,7,nx) - the analytic location).
C**** If a history file is specified then T_BH_INV is applied to
C**** the solution at each time instance and the solution written
C**** out to another history file.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER ISIZE_TBH(2),
     '  NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPLIST3(0:NPM),NPLIST4(0:NPM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NQNY(2,NYQM,0:NRCM,NXM),NVHP(NHM,NPM,NCM,0:NRM),NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  T_BH_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK1_INV(NY_TRANSFER_M),WK4_INV(NY_TRANSFER_M),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM)
      CHARACTER STRING*(MXCH),ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,IBEG1,icol,IEND,IEND1,IFAIL,IRANK,irow,isize,M,N,
     '  na,N1LIST,N3CO,nh1,nh2,nhx1,nhx2,NIQLIST(0:1),NIQSLIST(0:1),
     '  NIYLIST(0:16),nk1,nk2,nolist1,nolist2,nonr,no_nynr,notime,
     '  np,np1,np2,nr_inner,NRLIST_LOC(0:9),NRLIST2_LOC(0:9),nr_outer,
     '  NUMTIMEDATA,NUMTIMEDATA1,nv1,nv2,nx,nxc,ny1,ny2,
     '  TRSF_NRLIST2(0:9)
      REAL*8 S(NY_TRANSFER_M),SUM,TIME,TOL,YPMIN(16),
     '  YPMAX(16)
      CHARACTER ERROR1*255,FILEFORMAT*6,INFILE*(MXCH),OUTFILE*(MXCH)
      LOGICAL CBBREV,ENDFILE,HISTORY,INLIST,YPDATA,YQDATA,YQSDATA

      CALL ENTERS('APINVE',*9999)
      CALL STRING_TRIM(FILE00,IBEG1,IEND1)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM apply inverse
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Apply the inverse transfer matrix T_BH_INV to a set of nodal
C###    values.

C LKC 25-OCT-97 Made redundent by apnois (apply noise)
C C###  Parameter:      <noise NOISE_LEVEL#[0.0]>
C        OP_STRING(1)=STRING(1:IEND)//' <noise NOISE_LEVEL#[0.0]>'
C        OP_STRING(2)=BLANK(1:15)//'<class #[1]>'

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM apply inverse history
C###  Parameter:      <infile FILENAME[$current]>
C###    The input file of the nodal values
C###  Parameter:      <outfile FILENAME[$current_new]>
C###    The output file of the updated nodal values
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Apply the inverse transfer matrix T_BH_INV to a set of nodal
C###    values stored in a history file and write out the result
C###    to another history file.

C LKC 25-OCT-97 Made redundent by apnois (apply noise)
C C###  Parameter:      <noise NOISE_LEVEL#[0.0]>

        OP_STRING(1)=STRING(1:IEND)//' history <infile FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//' <outfile FILENAME['
     '    //FILE00(IBEG1:IEND1)//'_new]>'
        OP_STRING(3)=BLANK(1:15)//' <(ascii/binary)[ascii]>'
C LKC 25-OCT-97 Made redundent by apnois (apply noise)
C        OP_STRING(4)=BLANK(1:15)//'<noise NOISE_LEVEL#[0.0]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','APINVE',ERROR,*9999)
      ELSE
        HISTORY=.FALSE.
C*** Setup INPUT and OUTPUT files
        IF(CBBREV(CO,'HISTORY',2,noco+1,NTCO,N3CO)) THEN
          HISTORY=.TRUE.
          IF(CBBREV(CO,'INFILE',2,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            INFILE=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            INFILE=FILE00(IBEG1:IEND1)
          ENDIF

          IF(CBBREV(CO,'OUTFILE',2,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            OUTFILE=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            OUTFILE=FILE00(IBEG1:IEND1)//'_new'
          ENDIF

          CALL STRING_TRIM(INFILE,IBEG,IEND)
          CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)
          IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
     '      INFILE(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
            ERROR='>>Input and output file names must be different'
            GOTO 9999
          ENDIF

C LKC 17-JAN-97 Code repeated 3 times
C          CALL STRING_TRIM(INFILE,IBEG,IEND)
C          CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)
C          IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
C     '      INFILE(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
C            ERROR='>>Input and output file names must be different'
C            GOTO 9999
C          ENDIF
C          CALL STRING_TRIM(INFILE,IBEG,IEND)
C          CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)
C          IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
C     '      INFILE(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
C            ERROR='>>Input and output file names must be different'
C            GOTO 9999
C          ENDIF

          CALL STRING_TRIM(INFILE,IBEG,IEND)
          CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)
          IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
     '      INFILE(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
            ERROR='>>Input and output file names must be different'
            GOTO 9999
          ENDIF
C*** File Type
          IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF
        ENDIF
C LKC 25-OCT-97 Made redundent by apnois (apply noise)
C
CC!!! The variable NOISE_LEVEL has been put into the 'set but not
CC!!! used' .omit file to prevent it showing in the parsed fortran
CC!!! check output. When this piece of code has been developed remove
CC!!! NOISE_LEVEL from the .omit file  CS 12-FEB-1997
C        IF(CBBREV(CO,'NOISE',1,noco+1,NTCO,N3CO)) THEN
C          NOISE_LEVEL=RFROMC(CO(N3CO+1))
C          WRITE(OP_STRING,'('' >>Noise not yet implemented'')')
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C        ELSE
C          NOISE_LEVEL=0.0d0
C        ENDIF

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

C LKC for AJP NPLIST1 -> NPLIST3 and NPLIST -> NPLIST4
        CALL ASSERT(NPLIST3(0).GT.0,'>>Call define transfer first',
     '    ERROR,*9999)
        CALL ASSERT(NPLIST4(0).GT.0,'>>No body surface nodes ??',
     '    ERROR,*9999)
        CALL ASSERT(EVALUATE_INVERSE,'>>Evaluate inverse first',
     '    ERROR,*9999)
C LKC for AJP 191297 new asserts
        CALL ASSERT(ISIZE_TBH(1).GT.1,'>>Zero size inverse matrix?',
     '    ERROR,*9999)
        CALL ASSERT(ISIZE_TBH(2).GT.1,'>>Zero size inverse matrix?',
     '    ERROR,*9999)
C LKC for AJP 191297
C        nr_inner=TRSF_NR_INNER
C        nr_outer=TRSF_NR_OUTER
CC AJPs
        nr_inner=TRSF_NR_FIRST
        nr_outer=TRSF_NR_SECOND
CC AJPe

        IF(.NOT.HISTORY) THEN
          CALL ASSERT(NIYM.GE.7,'>> NIYM must be at least 7',ERROR,
     '      *9999)
C***      Move existing solution (if any) to the analytic solution
C***      location
          DO no_nynr=1,NYNR(0,0,1,nr_inner,nx)
            ny1=NYNR(no_nynr,0,1,nr_inner,nx)
            IF(NPNY(0,ny1,0,nx).EQ.1) THEN
              np=NPNY(4,ny1,0,nx)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            ENDIF
C AJPn 1/4/98
            IF(INLIST(np,NPLIST3(1),NPLIST3(0),N1LIST)) THEN
              !np is on the inner surface
              YP(ny1,7,nx)=YP(ny1,1,nx)
              YP(ny1,1,nx)=0.0d0
            ENDIF !inlist
C            YP(ny1,7,nx)=YP(ny1,1,nx)
C AJPn 1/4/98
          ENDDO !no_nynr

          IF(ICALC_TRANSFER.EQ.1) THEN
C***        !Transfer matrix has been inverted explicitly
C           !Apply T_bh_inv to heart surface nodes.
C           !Firstly check that T_BH_INV is not all zero (e.g. has
C           !not been read in from an incorrect file)
            IF(DABS(T_BH_INV(1,1)).LE.RDELTA) THEN
              sum=0.0d0

C              DO i=1,NY_TRANSFER_M
C                sum=sum+T_BH_INV(1,i)*T_BH_INV(1,i)
C              ENDDO
CC AJPs 191297
              DO irow=1,ISIZE_TBH(2)
                DO icol=1,ISIZE_TBH(1)
                  sum=sum+T_BH_INV(irow,icol)*T_BH_INV(irow,icol)
                ENDDO
              ENDDO
CC AJPe
              CALL ASSERT(SUM.GE.RDELTA,
     '          '>>T_BH_INV matrix is zero?',ERROR,*9999)
            ENDIF
            irow=0
            DO nolist1=1,NPLIST3(0) !list of heart nodes
              np1=NPLIST3(nolist1)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np1,.FALSE.,ERROR,*9999)
              DO nhx1=1,NHP(np1,nr_inner,nx)
                nh1=NH_LOC(nhx1,nx)
                DO nv1=1,NVHP(nh1,np1,1,nr_inner)
                  DO nk1=1,
     '              MAX(NKH(nh1,np1,1,nr_inner)-KTYP93(1,nr_inner),1)
! AJP 18/8/99
!                    ny1=NYNP(nk1,nv1,nh1,np1,1,1,nr_inner)
                    ny1=NYNP(nk1,nv1,nh1,np1,0,1,nr_inner)
C                   heart row number
                    irow=irow+1
                    icol=0
                    SUM=0.0d0
                    DO nolist2=1,NPLIST4(0) !list of torso nodes
                      np2=NPLIST4(nolist2)
                      DO nhx2=1,NHP(np2,nr_outer,nx)
                        nh2=NH_LOC(nhx2,nx)
                        DO nv2=1,NVHP(nh2,np2,1,nr_outer)
                          DO nk2=1,MAX(NKH(nh2,np2,1,nr_outer)-
     '                      KTYP93(1,nr_outer),1)
! AJP 18/8/99
!                            ny2=NYNP(nk2,nv2,nh2,np2,2,1,nr_outer)
                            ny2=NYNP(nk2,nv2,nh2,np2,0,1,nr_outer)
C                           body column number
                            icol=icol+1
                            SUM=SUM+T_BH_INV(irow,icol)*YP(ny2,1,nx)
                          ENDDO !nk2
                        ENDDO !nv2
                      ENDDO !nh2
                    ENDDO !np2
                    YP(ny1,1,nx)=SUM
                  ENDDO !nk1
                ENDDO !nv1
              ENDDO !nh1
            ENDDO !np1
          ELSE !ICALC_TRANSFER
C           Need to solve system of equations T_bh_inv*phi_h=phi_b
            IFAIL=0
            IF(IREGULARISE.EQ.1) THEN !No regularisation
C           !Store T_BH in T_BH_INV

CC AJPs 191297
              DO irow=1,ISIZE_TBH(1)
                DO icol=1,ISIZE_TBH(2)
c              DO irow=1,ISIZE_GK_BB
c                DO icol=1,ISIZE_GK_HH
                  T_BH_INV(irow,icol)=T_BH(irow,icol)
                ENDDO
              ENDDO
              M=ISIZE_TBH(1)
              N=ISIZE_TBH(2)
c              M=ISIZE_GK_BB
c              N=ISIZE_GK_HH
CC AJPe

              TOL=1.0D-5
            ELSE IF(IREGULARISE.EQ.2) THEN !SVD
C           !Store T_BH in T_BH_INV
C              DO irow=1,ISIZE_GK_BB
C                DO icol=1,ISIZE_GK_HH
C                  T_BH_INV(irow,icol)=T_BH(irow,icol)
C                ENDDO
C              ENDDO
C              M=ISIZE_GK_BB
C              N=ISIZE_GK_HH

CC AJPs
              DO irow=1,ISIZE_TBH(1)
                DO icol=1,ISIZE_TBH(2)
c              DO irow=1,ISIZE_GK_BB
c                DO icol=1,ISIZE_GK_HH
                  T_BH_INV(irow,icol)=T_BH(irow,icol)
                ENDDO
              ENDDO
              M=ISIZE_TBH(1)
              N=ISIZE_TBH(2)
c              M=ISIZE_GK_BB
c              N=ISIZE_GK_HH
CC AJPe

              TOL=RMIN_SVD
            ELSE IF(IREGULARISE.EQ.3) THEN !Twomey
              ERROR='>>Not implemented yet'
              GOTO 9999
            ELSE IF(IREGULARISE.GE.4.AND.IREGULARISE.LE.6) THEN !Tikhonov
              IF(IREGULARISE.EQ.4) THEN !zero-order Tikhonov
C             !Store [T_bh,TIKH_VALUE*Identity] in T_bh_inv

CC AJPs 191297
c                CALL ASSERT((ISIZE_GK_BB+ISIZE_GK_HH).LE.NY_TRANSFER_M,
c     '            '>>T_BH_INV too small',ERROR,*9999)
                CALL ASSERT((ISIZE_TBH(1)+ISIZE_TBH(2)).LE.
     '            NY_TRANSFER_M,
     '            '>>T_BH_INV too small',ERROR,*9999)
CC AJPe

CC AJPs
              DO irow=1,ISIZE_TBH(1)
                DO icol=1,ISIZE_TBH(2)
c              DO irow=1,ISIZE_GK_BB
c                DO icol=1,ISIZE_GK_HH
                    T_BH_INV(irow,icol)=T_BH(irow,icol)
                  ENDDO
                ENDDO
c                isize=ISIZE_GK_BB
                isize=ISIZE_TBH(1)
c                DO irow=1,ISIZE_GK_HH
                DO irow=1,ISIZE_TBH(2)
                  T_BH_INV(irow+isize,irow)=TIKH_VALUE
                ENDDO
c                M=ISIZE_GK_BB+ISIZE_GK_HH
c                N=ISIZE_GK_HH
                M=ISIZE_TBH(1)+ISIZE_TBH(2)
                N=ISIZE_TBH(2)
CC AJPe
                TOL=1.0D-5
              ELSE
                ERROR='>>That option not implemented yet'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>Not implemented yet'
              GOTO 9999
            ENDIF
C           Store torso surface solutions in WK4_INV
            irow=0
            DO nolist2=1,NPLIST4(0) !list of torso nodes
              np2=NPLIST4(nolist2)
              DO nhx2=1,NHP(np2,nr_outer,nx)
                nh2=NH_LOC(nhx2,nx)
                DO nv2=1,NVHP(nh2,np2,1,nr_outer)
                  DO nk2=1,
     '              MAX(NKH(nh2,np2,1,nr_outer)-KTYP93(1,nr_outer),1)
! AJP 18/8/99
!                    ny2=NYNP(nk2,nv2,nh2,np2,2,1,nr_outer)
                    ny2=NYNP(nk2,nv2,nh2,np2,0,1,nr_outer)
C                 body column number of GK
                    irow=irow+1
                    WK4_INV(irow)=YP(ny2,1,nx)
                  ENDDO !nk2
                ENDDO !nv2
              ENDDO !nh2
            ENDDO !np2
C           Least squares solution
C           CALL F04JGF(M,N,T_BH_INV,NY_TRANSFER_M,WK4_INV,TOL,SVD,SIGMA,
C     '       IRANK,WK1_INV,4*NY_TRANSFER_M,IFAIL)
            CALL DGELSS(M,N,1,T_BH_INV,NY_TRANSFER_M,WK4_INV,
     '        NY_TRANSFER_M,S,TOL,IRANK,WK1_INV,4*NY_TRANSFER_M,IFAIL)

            IF(IFAIL.NE.0)THEN
              WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in inverse solution'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ERROR='>>IFAIL<>0 in DGELSS'
              GOTO 9999
            ENDIF
            IF(DOP) THEN
              DO i=1,N
                WRITE(OP_STRING,'('' SOLUTION '',D12.4)') WK4_INV(i)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
            WRITE(OP_STRING,'('' Rank of matrix is '',I5)')IRANK
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C***        Put solution into YP
            irow=0
            DO nolist1=1,NPLIST3(0) !list of heart nodes
              np1=NPLIST3(nolist1)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np1,.FALSE.,ERROR,*9999)
              DO nhx1=1,NHP(np1,nr_inner,nx)
                nh1=NH_LOC(nhx1,nx)
                DO nv1=1,NVHP(nh1,np1,1,nr_inner)
                  DO nk1=1,
     '              MAX(NKH(nh1,np1,1,nr_inner)-KTYP93(1,nr_inner),1)
                    ny1=NYNP(nk1,nv1,nh1,np1,1,1,nr_inner)
C                   heart row number of GK
                    irow=irow+1
                    YP(ny1,1,nx)=WK4_INV(irow)
                  ENDDO !nk1
                ENDDO !nv1
              ENDDO !nh1
            ENDDO !np1
          ENDIF !ICALC_TRANSFER
        ELSE !History file to be transformed.
          CALL ASSERT(ICALC_TRANSFER.EQ.1,
     '      '>>Can only apply an explicit inverse to a history file',
     '      ERROR,*9999)
C***      !Firstly check that T_BH_INV is not all zero (e.g. has
C***      !not been read in from an incorrect file)
          IF(DABS(T_BH_INV(1,1)).LE.RDELTA) THEN
            sum=0.0d0

C            DO i=1,NY_TRANSFER_M
C              sum=sum+T_BH_INV(1,i)*T_BH_INV(1,i)
C            ENDDO

CC AJPs 191297
              DO irow=1,ISIZE_TBH(2)
                DO icol=1,ISIZE_TBH(1)
                  sum=sum+T_BH_INV(irow,icol)*T_BH_INV(irow,icol)
                ENDDO
              ENDDO
CC AJPe

            CALL ASSERT(SUM.GE.RDELTA,
     '        '>>T_BH_INV matrix is zero?',ERROR,*9999)
          ENDIF

C LKC 24-APR-1998 Initialise NIYLIST
          NIYLIST(0)=0
          NIQLIST(0)=0
          NIQSLIST(0)=0
          na=1
C CPB 15/10/98 Default to 1-1 region mapping
          TRSF_NRLIST2(0)=TRSF_NRLIST(0)
          DO nonr=1,TRSF_NRLIST(0)
            TRSF_NRLIST2(nonr)=TRSF_NRLIST(nonr)
          ENDDO

C***      Open the INPUT history file
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '      NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,
     '      INFILE,'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

C***      Check on the number of regions in the history file and
C***      decide what to do if different.
          IF(NRLIST_HIST(0,IOFILE1).LT.TRSF_NRLIST(0)) THEN
C***        Less regions than current setup in the history file.  This
C           can be obtained by e.g. fitting measured signals on the
C           outside of a homogeneous torso and then applying a
C           multi-region inverse matrix to it.
            IF(NRLIST_HIST(0,IOFILE1).EQ.1) THEN
C***          In this case assume that the data should be read into the
C             last location
              NRLIST_LOC(0)=1
! AJP 17/8             NRLIST_LOC(1)=NRLIST_HIST(0,IOFILE1)
              NRLIST_LOC(1)=NRLIST_HIST(1,IOFILE1)
            ELSE
              WRITE(OP_STRING,
     '          '('' >> WARNING: Do not know where to store data'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ELSEIF(NRLIST_HIST(0,IOFILE1).GT.TRSF_NRLIST(0)) THEN
C***        More regions were used to generate the data than the
C           current setup has.
            IF(NRLIST_HIST(0,IOFILE1).EQ.TRSF_NRLIST(0)+1) THEN
C***          Will occur when e.g. the data has been generated using
C             a moving dipole inside the heart - the interior of the
C             heart being the extra region not needed in the inverse
              CALL ASSERT(NRLIST_HIST(0,IOFILE1).LE.9,
     '          '>>NRLIST_LOC not big enough',ERROR,*9999)
              NRLIST_LOC(0)=TRSF_NRLIST(0)
              DO nolist1=1,NRLIST_LOC(0)
                NRLIST_LOC(nolist1)=TRSF_NRLIST(nolist1)
              ENDDO
              DO nolist1=NRLIST_LOC(0)+1,NRLIST_HIST(0,IOFILE1)
C               !initialising rest of NRLIST_LOC - needed in IOHIST.
                NRLIST_LOC(nolist1)=0
              ENDDO
            ELSE
              WRITE(OP_STRING,
     '          '('' >> WARNING: Do not know where to store data'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
C***        Problem setup and history file have the same number of
C           regions.  Could check that the number of nys etc are the
C           same, but for now just assume it is okay.
            NRLIST_LOC(0)=TRSF_NRLIST(0)
            DO nolist1=1,NRLIST_LOC(0)
              NRLIST_LOC(nolist1)=TRSF_NRLIST(nolist1)
            ENDDO
          ENDIF
          NIYLIST(0)=1
          NIYLIST(1)=1

C LKC 31-OCT-97 new check for ascii
          IF(FILEFORMAT.EQ.'ASCII') THEN !need to find NUMTIMEDATA
            NUMTIMEDATA=0
            ENDFILE=.FALSE.
            DO WHILE(.NOT.ENDFILE)
              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,
     '          TRSF_NRLIST2,NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),
     '          NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,YPMIN,
     '          YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,INFILE,'TIME_DATA',
     '          ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
              IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
            ENDDO
C*** RESET the input file
            CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '        NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '        TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '        FILEFORMAT,INFILE,'RESET',ENDFILE,.TRUE.,YPDATA,YQDATA,
     '        YQSDATA,ERROR,*9999)
          ENDIF
          CALL ASSERT(NUMTIMEDATA.GT.0,
     '      '>>No signal, NUMTIMEDATA.LE.0',ERROR,*9999)

C***      Open the OUTPUT history file
          CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '      NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',FILEFORMAT,
     '      OUTFILE,'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

C 30-OCT-97 LKC This section only works for ASCII
C 31-OCT-97 Rewritten to handle binary files
C          ENDFILE=.FALSE.
C          DO WHILE(.NOT.ENDFILE)

C CPB 15/10/98 Default to 1-1 region mapping
          NRLIST2_LOC(0)=NRLIST_LOC(0)
          DO nonr=1,NRLIST_LOC(0)
            NRLIST2_LOC(nonr)=NRLIST_LOC(nonr)
          ENDDO


          DO notime=1,NUMTIMEDATA !loop over all the time steps
C***        Read YP for each time series
            CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST_LOC,NRLIST2_LOC,
     '        NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '        TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '        FILEFORMAT,INFILE,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,
     '        YQDATA,YQSDATA,ERROR,*9999)

            irow=0
            DO nolist1=1,NPLIST3(0) !list of heart nodes
              np1=NPLIST3(nolist1)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np1,.FALSE.,ERROR,*9999)
              DO nhx1=1,NHP(np1,nr_inner,nx)
                nh1=NH_LOC(nhx1,nx)
                DO nv1=1,NVHP(nh1,np1,1,nr_inner)
                  DO nk1=1,
     '              MAX(NKH(nh1,np1,1,nr_inner)-KTYP93(1,nr_inner),1)
! AJP 18/8/99
!                    ny1=NYNP(nk1,nv1,nh1,np1,1,1,nr_inner)
                    ny1=NYNP(nk1,nv1,nh1,np1,0,1,nr_inner)
C                   heart row number
                    irow=irow+1
                    icol=0
                    SUM=0.0d0
                    DO nolist2=1,NPLIST4(0) !list of torso nodes
                      np2=NPLIST4(nolist2)
                      DO nhx2=1,NHP(np2,nr_outer,nx)
                        nh2=NH_LOC(nhx2,nx)
                        DO nv2=1,NVHP(nh2,np2,1,nr_outer)
                          DO nk2=1,MAX(NKH(nh2,np2,1,nr_outer)-
     '                      KTYP93(1,nr_outer),1)
! AJP 18/8/99
!                            ny2=NYNP(nk2,nv2,nh2,np2,2,1,nr_outer)
                            ny2=NYNP(nk2,nv2,nh2,np2,0,1,nr_outer)
C                           body column number
                            icol=icol+1
                            SUM=SUM+T_BH_INV(irow,icol)*YP(ny2,1,nx)
                          ENDDO !nk2
                        ENDDO !nv2
                      ENDDO !nh2
                    ENDDO !np2
                    YP(ny1,1,nx)=SUM
                  ENDDO !nk1
                ENDDO !nv1
              ENDDO !nh1
            ENDDO !np1
C***        Write out new solution

C 12-APR-98 Pass in Variables not constants
            YPDATA=.TRUE.
            YQDATA=.FALSE.

            CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '        NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '        TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',
     '        FILEFORMAT,OUTFILE,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,
     '        YQDATA,YQSDATA,ERROR,*9999)
          ENDDO !nottime
C          ENDDO !Hit end of file

C*** Close both history files
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST_LOC,NRLIST2_LOC,
     '      NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '      INFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
          CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '      NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '      OUTFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

        ENDIF !History file transformed.
      ENDIF

      APPLY_INVERSE=.TRUE.
      CALL EXITS('APINVE')
      RETURN
 9999 CALL ERRORS('APINVE',ERROR)
C*** Close both history files correctly
      CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '  NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST_LOC,NRLIST2_LOC,
     '  NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '  YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '  INFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR1,*9998)
 9998 CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '  NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '  NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '  YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '  OUTFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR1,*9997)
 9997 CALL EXITS('APINVE')
      RETURN 1
      END


