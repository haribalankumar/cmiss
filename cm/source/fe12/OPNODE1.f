      SUBROUTINE OPNODE1(NHP,NKH,NKJ,np,NP_INTERFACE,nr,NVJP,NWP,nx,
     '  XP,ZP,EXTEND,INTERFACE,ERROR,*)

C#### Subroutine: OPNODE1
C###  Description:
C###    OPNODE1 outputs nodal coordinates for node np.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NHP(NPM,0:NRM,NXM),NKH(NHM,NPM),
     '  NKJ(NJM),np,NP_INTERFACE(0:NPM,0:3),
     '  nr,NVJP(NJM,NPM),NWP(NPM,2),nx
      REAL*8 XP(NKM,NVM,NJM),ZP(NKM,NVM,NHM)
      CHARACTER ERROR*(*)
      LOGICAL EXTEND,INTERFACE

!     Local Variables
      INTEGER NJ_LOC_MX_LOCAL
      PARAMETER (NJ_LOC_MX_LOCAL=3)
      INTEGER i,IB4,IDIGITS,IE4,MAXNJLEN,nh,nhx,nj,njj1,njj2,nk,
     '  NKJMAX,nkjPERLINE,nkjTHISLINE,nkjTOWRITE,
     '  NKHMAX,nkhPERLINE,nkhTHISLINE,nkhTOWRITE,
     '  NUM_NH_NK,nv
      CHARACTER CHAR1*1,CHAR2*1,CHAR3*2,CHAR4*16,
     '  FIBRETYPE(NJ_LOC_MX_LOCAL)*6,GEOMTYPE(NJ_LOC_MX_LOCAL,5)*6,
     '  LEAD_STRING*30,NJTYPE(3,NJ_LOC_MX)*9,
     '  NODE_STRING*9,NUMSPACES*2,NV_STRING*20

      DATA GEOMTYPE /'x     ','y     ','z     ',
     '               'r     ','theta ','z     ',
     '               'r     ','theta ','phi   ',
     '               'lambda','mu    ','theta ',
     '               'lambda','mu    ','theta '/
      DATA FIBRETYPE/'fibre ','imbric','sheet '/

      CALL ENTERS('OPNODE1',*9999)

      DO njj2=1,MIN(NJ_LOC_MX_LOCAL,NJ_LOC_MX)
        NJTYPE(1,njj2)=GEOMTYPE(njj2,ITYP10(nr))
        NJTYPE(2,njj2)=FIBRETYPE(njj2)
      ENDDO
      DO njj2=NJ_LOC_MX_LOCAL+1,NJ_LOC_MX
        NJTYPE(1,njj2)='      '
        NJTYPE(2,njj2)='      '
      ENDDO
      DO njj2=1,NJ_LOC_MX
        WRITE(CHAR4,'(I2)') njj2
        CALL STRING_TRIM(CHAR4,IB4,IE4)
        NJTYPE(3,njj2)='field'//CHAR4(IB4:IE4)
      ENDDO
      MAXNJLEN=5+IDIGITS(NJ_LOC(NJL_FIEL,0,nr))

      IF(EXTEND) THEN
        NKJMAX=0
        DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj2,nr)
          IF(NKJ(nj).GT.NKJMAX) NKJMAX=NKJ(nj)
        ENDDO !njj2
C new MPN 17Apr97: better formatting for tricubic nj's
        nkjPERLINE=MIN(NKJMAX,4)
        NKHMAX=0
        DO nhx=1,NHP(np,nr,nx)
          nh=NH_LOC(nhx,nx)
          NUM_NH_NK=MAX(NKH(nh,np)-KTYP93(1,nr),1)
          IF(NUM_NH_NK.GT.NKHMAX) NKHMAX=NUM_NH_NK
        ENDDO !nhx
        nkhPERLINE=MIN(NKHMAX,4)

        CHAR3=' '
        DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj2,nr)
          NV_STRING='['//CHAR3(1:2)//']'
          NUM_NH_NK=0
          IF(njj2.LE.NHP(np,nr,nx)) THEN
            nh=NH_LOC(njj2,nx)
            NUM_NH_NK=MAX(NKH(nh,np)-KTYP93(1,nr),1)
          ENDIF
          IF(njj2.EQ.1) THEN
            WRITE(NODE_STRING,'(A4,I5)') 'Node',np
          ELSE
            NODE_STRING='         '
          ENDIF
          DO nv=1,NVJP(nj,np)
            IF(NVJP(nj,np).GT.1) THEN
              WRITE(CHAR3,'(I2)') nv
              NV_STRING='['//CHAR3(1:2)//']'
            ELSE
              NV_STRING='    '
            ENDIF
C new MPN 17Apr97: better formatting for tricubic nj's
            nkjTOWRITE=NKJ(nj)
            nkhTOWRITE=NUM_NH_NK
            DO WHILE((nkjTOWRITE.GT.0).OR.(nkhTOWRITE.GT.0))
              nkjTHISLINE=MIN(nkjTOWRITE,nkjPERLINE)
              nkhTHISLINE=MIN(nkhTOWRITE,nkhPERLINE)
              IF(nkjTOWRITE.EQ.NKJ(nj)) THEN !first line
                WRITE(LEAD_STRING,'(2X,A9,'' nj='',I2,'' ('','
     '            //'A6,'')'',A4)') NODE_STRING,nj,
     '            NJTYPE(NJ_TYPE(nj,1),NJ_TYPE(nj,2))(:MAXNJLEN),
     '            NV_STRING
              ELSE
                WRITE(LEAD_STRING,'(30X)')
              ENDIF
              WRITE(CHAR1,'(I1)') nkjTHISLINE
              WRITE(NUMSPACES,'(I2)') 12*(nkjPERLINE-nkjTHISLINE)+1
              WRITE(CHAR2,'(I1)') nkhTHISLINE
C C MPN Nov2004: the following 'fix' doesn't work for 2D/3D cubic Hermite 
C C              mechanics reactions.
C RGB 29/4/98   For fluid stuff-want to solve for velocity and pressure
C               Assume only one derivative for now
             IF(njj2.EQ.1
     &          .AND.NKJ(nj).EQ.1 ! not implemented for > 1
     &          .AND.ITYP5(nr,nx).EQ.2 ! time integration
     &          .AND.ITYP2(nr,nx).EQ.5 ! Navier Stokes
     &          .AND.ITYP4(nr,nx).EQ.5 ! Finite Volume. Any others?
     &          .AND.NHP(np,nr,nx).EQ.(NJT+1)) THEN
C MJS Fixed for linux, format could be set to 0D12.4
C which is illegal
                 IF(CHAR1.EQ.'0'.AND.CHAR2.EQ.'0') THEN
                    FORMAT='(A30,1X,D12.4,'//NUMSPACES//'X,D12.4,
     '                   D12.4)'
                    WRITE(OP_STRING,FORMAT) LEAD_STRING,
     '                   ZP(1,nv,NHP(np,nr,nx))
                 ELSE IF(CHAR1.EQ.'0') THEN
                    FORMAT='(A30,1X,'//NUMSPACES//'X,'//CHAR2//'D12.4,
     '                   D12.4)'
                    WRITE(OP_STRING,FORMAT) LEAD_STRING,
     '                   ZP(1,nv,nh),ZP(1,nv,NHP(np,nr,nx))
                 ELSE IF(CHAR2.EQ.'0') THEN
                    FORMAT='(A30,1X,'
     '                   //CHAR1//'D12.4,'//NUMSPACES//'X,
     '                   D12.4)'
                    WRITE(OP_STRING,FORMAT) LEAD_STRING,
     '                   XP(1,nv,nj),ZP(1,nv,NHP(np,nr,nx))
                 ELSE
                    FORMAT='(A30,1X,'
     '                   //CHAR1//'D12.4,'//NUMSPACES//'X,'
     '                   //CHAR2//'D12.4,D12.4)'
                    WRITE(OP_STRING,FORMAT) LEAD_STRING,
     '                   XP(1,nv,nj),ZP(1,nv,nh),ZP(1,nv,NHP(np,nr,nx))
                 ENDIF

                 CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                 nkjTOWRITE=0
                 nkhTOWRITE=0
               ELSE ! not fluid-stuff
                 IF(CHAR1.EQ.'0'.AND.CHAR2.EQ.'0') THEN
                    FORMAT='(A30)'
                    WRITE(OP_STRING,FORMAT) LEAD_STRING
                 ELSE IF(CHAR1.EQ.'0') THEN
                    FORMAT='(A30,1X,'//NUMSPACES//'X,'//CHAR2//'D12.4)'
                    WRITE(OP_STRING,FORMAT) LEAD_STRING,
     '                   (ZP(nk+NUM_NH_NK-nkhTOWRITE
     '                   ,nv,nh),nk=1,nkhTHISLINE)
                 ELSE IF(CHAR2.EQ.'0') THEN
                    FORMAT='(A30,1X,'//CHAR1//'D12.4)'
                    WRITE(OP_STRING,FORMAT) LEAD_STRING,
     '                   (XP(nk+NKJ(nj)-nkjTOWRITE,nv,nj)
     '                   ,nk=1,nkjTHISLINE)
                 ELSE
                    FORMAT='(A30,1X,'
     '                   //CHAR1//'D12.4,'//NUMSPACES//'X,'
     '                   //CHAR2//'D12.4)'
                    WRITE(OP_STRING,FORMAT) LEAD_STRING,
     '                   (XP(nk+NKJ(nj)-nkjTOWRITE,nv,nj)
     '                   ,nk=1,nkjTHISLINE),(ZP(nk+NUM_NH_NK-nkhTOWRITE
     '                   ,nv,nh),nk=1,nkhTHISLINE)
                 ENDIF
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                nkjTOWRITE=nkjTOWRITE-nkjTHISLINE
                nkhTOWRITE=nkhTOWRITE-nkhTHISLINE
              ENDIF ! fluid-stuff/not
            ENDDO !while
C old MPN 17Apr97: bad formatting for tricubic nj's
C            WRITE(CHAR1,'(I1)') MIN(NKJ(nj),4)
C            WRITE(CHAR2,'(I2)') 12*(NKJMAX-NKJ(nj))
C            FORMAT='(2X,A9,'' nj='',I2,'' ('
C     '        //NJTYPE(NJ_TYPE(nj,1),NJ_TYPE(nj,2))
C     '        //')'//NV_STRING(1:4)//' '','
C     '        //CHAR1//'D12.4,'//CHAR2//'X,1X,4D12.4,(/31X,'
C     '        //CHAR1//'D12.4,'//CHAR2//'X,1X,4D12.4))'
C            WRITE(OP_STRING,FORMAT) NODE_STRING,nj,
C     '        (XP(nk,nv,nj),nk=1,NKJ(nj)),
C     '        (ZP(nk,nv,nh),nk=1,NUM_NH_NK)
C            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO  !nv
        ENDDO  !nj

        DO njj1=2,3                   !Loop over fibres and fields
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            DO nv=1,NVJP(nj,np)
              IF(NVJP(nj,np).GT.1) THEN
                WRITE(CHAR3,'(I2)') nv
                NV_STRING='['//CHAR3(1:2)//']'
              ELSE
                NV_STRING='    '
              ENDIF
              FORMAT='(12X,''nj='',I2,'' ('
     '          //NJTYPE(NJ_TYPE(nj,1),NJ_TYPE(nj,2))(:MAXNJLEN)
     '          //')'//NV_STRING(1:4)//' '',4D12.4,'
     '          //'(/31X,4D12.4))'
              WRITE(OP_STRING,FORMAT) nj,(XP(nk,nv,nj),nk=1,NKJ(nj))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !nv
          ENDDO !njj2
        ENDDO !njj1

      ELSE IF(INTERFACE) THEN
        IF(NP_INTERFACE(np,1).EQ.0) THEN
          WRITE(OP_STRING,'(1X,''Node'',I5,'' belongs to no'
     '      //' interfaces'')') np
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'(1X,''Node'',I5,'' belongs to '',I2,'
     '      //''' region(s): '',9I3)') np,NP_INTERFACE(np,0),
     '      (NP_INTERFACE(np,i),i=1,NP_INTERFACE(np,0))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE
        DO njj1=1,3      !Loop over geometry, fibres and fields
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            DO nv=1,NVJP(nj,np)
              IF(NVJP(nj,np).GT.1) THEN
                WRITE(CHAR3,'(I2)') nv
                NV_STRING='['//CHAR3(1:2)//']'
              ELSE
                NV_STRING='    '
              ENDIF
              IF(nj.EQ.1.AND.nv.EQ.1) THEN
                FORMAT='(2X,''Node'',I5,'' nj= 1 ('
     '            //NJTYPE(NJ_TYPE(1,1),NJ_TYPE(1,2))(:MAXNJLEN)
     '            //')'//NV_STRING(1:4)//' '',4D12.4,'
     '            //'(/31X,4D12.4))'

C!!! LKC 17-NOV-2005 Shouldn't NKJ(1) really be NKJ(nj) below?!
                WRITE(OP_STRING,FORMAT) np,
     '            (XP(nk,nv,1),nk=1,NKJ(1))
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                
              ELSE
                IF((ITYP10(nr).EQ.2.AND.nj.EQ.2).OR.
     '             (ITYP10(nr).GT.2.AND.nj.LE.NJT).OR.
     '            (nj.EQ.NJ_LOC(NJL_FIBR,1,nr)).OR.
     '            (nj.EQ.NJ_LOC(NJL_FIBR,2,nr)).OR.
     '            (nj.EQ.NJ_LOC(NJL_FIBR,3,nr))) THEN
                  FORMAT='(12X,''nj='',I2,'' ('
     '              //NJTYPE(NJ_TYPE(nj,1),NJ_TYPE(nj,2))(:MAXNJLEN)
     '              //')'//NV_STRING(1:4)//' '',D12.4,'
     '              //'1X,''('',F7.1,''deg)'',4D12.4,'
     '              //'(/56X,4D12.4))'
                  WRITE(OP_STRING,FORMAT)
     '              nj,XP(1,nv,nj),XP(1,nv,nj)*180.0d0/PI,
     '              (XP(nk,nv,nj),nk=2,NKJ(nj))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ELSE
                  FORMAT='(12X,''nj='',I2,'' ('
     '              //NJTYPE(NJ_TYPE(nj,1),NJ_TYPE(nj,2))(:MAXNJLEN)
     '              //')'//NV_STRING(1:4)//' '',4D12.4,'
     '              //'(/31X,4D12.4))'
                  WRITE(OP_STRING,FORMAT) nj,
     '              (XP(nk,nv,nj),nk=1,NKJ(nj))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO  !nv
          ENDDO  !njj2
        ENDDO  !njj1
      ENDIF

      IF(JTYP2C.EQ.1) THEN ! hanging nodes
        IF(NWP(np,1).GT.0) THEN
          WRITE(OP_STRING,'(''       Hangs in element '',I4)') NWP(np,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('OPNODE1')
      RETURN
 9999 CALL ERRORS('OPNODE1',ERROR)
      CALL EXITS('OPNODE1')
      RETURN 1
      END


