      SUBROUTINE IPMAPPING(NBJ,NEELEM,NKJ,
     '  NPNODE,NPNY,nr,NVJP,nx,NYNP,NYNY,
     '  CYNY,XP,ZP,ERROR,*)

C#### Subroutine: IPMAPPING
C###  Description:
C###    IPMAPPNG inputs ny->ny mapping

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'parameters.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKJ(NJM,NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),nr,NVJP(NJM,NPM),
     '  nx,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNY(0:NYYM,NYM,NRM,NXM)
      REAL*8 CYNY(0:NYYM,NYM,NRM,NXM),XP(NKM,NVM,NJM,NPM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,MAP_NODES,nj,nj2,nk,nonode,NOQUES,np,
     '  nu,NUK(8),nv,nv2,NV_MAX,ny,ny_target,
     '  nk_t,nv_t,nj_t,np_t,nk_m,nv_m,nj_m,np_m
      CHARACTER CHAR1*5,CHAR2*1,CHAR3*3
      LOGICAL FILEIP

      DATA NUK/1,2,4,6,7,9,10,11/

      CALL ENTERS('IPMAPPING',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='($,'' Define node position mapping [N]? '',A)'
      ADEFLT(1)='N'
      ADATA(1)='N'
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(ADATA(1).EQ.'Y') THEN
        ITYP21(nr)=3
      ELSE
        ITYP21(nr)=2
      ENDIF

      IDEFLT(1)=1
      WRITE(CHAR1,'(I5)') IDEFLT(1)
      IF(IOTYPE.EQ.3) THEN
        IDATA(1)=1
      ENDIF
      FORMAT='($,'' The number of nodes with '
     '  //'special mappings is ['//CHAR1(1:5)//']: '',I5)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NPM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        MAP_NODES=IDATA(1)
      ENDIF

      IF(MAP_NODES.GT.0) THEN
        NOQUES=0
        DO nonode=1,MAP_NODES
          IDEFLT(1)=1
          WRITE(CHAR1,'(I5)') IDEFLT(1)
          IF(IOTYPE.EQ.3) THEN
            np=NPNODE(nonode,nr)
            IDATA(1)=NP
          ENDIF
          FORMAT='(/$,'' Node number ['//CHAR1(1:5)//']: '',I5)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            np=IDATA(1)
          ENDIF

          DO nj=1,NJT
            WRITE(CHAR1,'(I1)') nj
            FORMAT='('' For the Xj('//CHAR1(1:1)//') '
     '        //'coordinate: '')'
              CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '          ERROR,*9999)

            NV_MAX=NVJP(nj,np)
            DO nv=1,NV_MAX
              IF(NV_MAX.GT.1) THEN !ask for diff nj versions
                WRITE(CHAR1,'(I2)') nv
                FORMAT='('' For version number'//CHAR1(1:2)//':'')'
                CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '            0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
              ENDIF

              WRITE(CHAR1,'(I1)') nj

              DO nk=1,NKJ(nj,np)
                nu=NUK(nk)
                IF((nu.EQ.1).AND.(ITYP21(nr).EQ.3)) THEN
                  FORMAT='($,'' Is the nodal position '
     '              //'mapped out [N]? '',A)'
                ELSE IF(nu.EQ.2.OR.nu.EQ.4.OR.nu.EQ.7) THEN
                  IF(nu.EQ.2) THEN
                    CHAR2='1'
                  ELSE IF(nu.EQ.4) THEN
                    CHAR2='2'
                  ELSE IF(nu.EQ.7) THEN
                    CHAR2='3'
                  ENDIF
                  FORMAT='($,'' Is the derivative wrt direction '//
     '              CHAR2(1:1)//' is mapped out [N]? '',A)'
                ELSE IF(nu.EQ.6.OR.nu.EQ.9.OR.nu.EQ.10) THEN
                  IF(nu.EQ.6) THEN
                    CHAR2='1'
                    CHAR3='2'
                  ELSE IF(nu.EQ.9) THEN
                    CHAR2='1'
                    CHAR3='3'
                  ELSE IF(nu.EQ.10) THEN
                    CHAR2='2'
                    CHAR3='3'
                  ENDIF
                  FORMAT='($,'' Is the derivative wrt directions '//
     '              CHAR2(1:1)//' & '//CHAR3(1:1)//
     '              ' is mapped out [N]? '',A)'
                ELSE IF(nu.EQ.11) THEN
                  FORMAT='($,'' Is the derivative wrt directions '//
     '              '1, 2 & 3 is mapped out [N]? '',A)'
                ENDIF

                IF((nu.NE.1).OR.(ITYP21(nr).EQ.3)) THEN
                  IF(IOTYPE.EQ.3) ADATA(1)='N'
                  CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,
     '            NOQUES,FILEIP,FORMAT,1,
     '            ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '            LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                  IF(ADATA(1).EQ.'Y') THEN
                    ny=NYNP(nk,nv,nj,np,0,1,nr)
                       IF(DOP) THEN
                      WRITE(OP_STRING,
     &                  '(''INFO:'
     &                  //' ny for node '',I4,'' vers '',I4,
     &                  '' deriv '',I4,'' dir '',I4,'' ny '',I5)')
     &                  np,nv,nk,nj,ny
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                    ! warn if there is an attempt to map to a non-existant dof
                    IF(ny.LE.0) THEN
                      WRITE(OP_STRING,
     &                  '(''>>WARNING: ignoring non-existant source'
     &                  //' dof for node '',I4,'' vers '',I4,
     &                  '' deriv '',I4,'' dir '',I4,'' - equation'
     &                  //' must be defined before mapping'')')
     &                  np,nv,nk,nj
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
!                     CALL ASSERT(ny.GT.0,
!      '                '>>Attempt to map undefined mesh dof (ny=0)',
!      '                ERROR,*9999)
                    FORMAT='($,'' Enter node, version, direction, '
     '                //'derivative numbers to '
     '                //'map to [1,1,1,1]: '',4I5)'
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,
     '                NOQUES,FILEIP,FORMAT,1,
     '                ADATA,ADEFLT,CDATA,CDEFLT,
     '                ICHAR,IDATA,IDEFLT,0,NPM,
     '                LDATA,LDEFLT,RDATA,RDEFLT,
     '                RMIN,RMAX,INFO,ERROR,*9999)

                    ny_target=NYNP(IDATA(4),IDATA(2),
     '                  IDATA(3),IDATA(1),0,1,nr)

                    IF(ny.GT.0) THEN
                      CALL ASSERT(NYNY(0,ny,nr,nx)+1.LE.NYYM,
     '                   '>>Increase NYYM',ERROR,*9999)
                      ! 0th element of NYNY is the number of array elements
                      NYNY(0,ny,nr,nx)=NYNY(0,ny,nr,nx)+1 ! increment array size
                      NYNY(NYNY(0,ny,nr,nx),ny,nr,nx)=ny_target
                    ENDIF

                    FORMAT='($,'' Enter the mapping '
     '                //'coefficient [1]: '',G12.5)'
                    RDATA(1)=1.0d0
                    RDEFLT(1)=1.0d0
                    CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,
     '                CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     '                IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '                -RMAX,RMAX,INFO,ERROR,*9999)

                    IF(ny.GT.0) THEN
                      CYNY(0,ny,nr,nx)=0.0d0
                      CYNY(NYNY(0,ny,nr,nx),ny,nr,nx)=RDATA(1)

C*** One should set up the initial node file with the correct initial values
C*** for the mapped out dofs, however that is extremely tedious so set them
C*** here for the simple one to one case.
C*** this is not at all general as the there could non unit coefficients
C*** but will suffice for the envisioned use in the near future.

                      nk_t=NPNY(1,ny_target,0,1)
                      nv_t=NPNY(2,ny_target,0,1)
                      nj_t=NPNY(3,ny_target,0,1)
                      np_t=NPNY(4,ny_target,0,1)
                      nk_m=NPNY(1,ny,0,1)
                      nv_m=NPNY(2,ny,0,1)
                      nj_m=NPNY(3,ny,0,1)
                      np_m=NPNY(4,ny,0,1)

                      XP(nk_m,nv_m,nj_m,np_m)=
     '                  XP(nk_t,nv_t,nj_t,np_t)*RDATA(1)

                      ZP(nk_m,nv_m,nj_m,np_m,1)=
     '                  XP(nk_t,nv_t,nj_t,np_t)*RDATA(1)

                      DO nv2=1,NV_MAX
                         DO nj2=1,NJT
                           ZP(1,nv2,nj2,np_m,1)=
     '                       ZP(1,nv_t,nj2,np_t,1)
                         ENDDO
                      ENDDO ! nv2
                    ENDIF ! ny.GT.0

                  ENDIF ! ADATA(1).EQ.'Y'
                ENDIF !(nu.NE.1).OR.(ITYP21(nr).EQ.3)
              ENDDO !nk
            ENDDO !nv
          ENDDO !nj
        ENDDO !End of nonode loop
      ENDIF

      IF(ITYP21(nr).EQ.2) THEN !Set the node coords the same as version 1
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nj=1,NJT
            NV_MAX=NVJP(nj,np)
            DO nv2=2,NV_MAX
               XP(1,nv2,nj,np)=XP(1,1,nj,np)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('IPMAPPING')
      RETURN
 9999 CALL ERRORS('IPMAPPING',ERROR)
      CALL EXITS('IPMAPPING')
      RETURN 1
      END


