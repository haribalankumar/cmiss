      SUBROUTINE OPSEGM(ISAXES,ISBASE,ISCONO,ISCONT,ISDANO,ISDAPR,
     '  ISDATA,ISDATR,ISEG,ISELNO,ISFACE,ISFANO,ISFIBR,ISGAUS,
     '  ISGRAD,ISGRID,ISHIST,ISINCR,ISLINE,ISLINO,ISMAP,
     '  ISMATE,ISNONO,ISREAC,ISRULE,ISSECT,ISSTRA,ISSTRE,
     '  ISSTRM,ISSURF,ISVELO,NEELEM,NPNODE,
     '  CSEG,FULL,ERROR,*)

C#### Subroutine: OPSEGM
C###  Description:
C###    OPSEGM outputs segments.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISAXES(NWM),ISBASE(99),ISCONO(NHM,NEM),
     '  ISCONT(NHM,NEM,NGRSEGM),ISDANO(NWM,NEM),ISDAPR(NWM,NEM),
     '  ISDATA(NWM,NGRSEGM),ISDATR(NWM,NEM),ISEG(*),ISELNO(NWM,NEM),
     '  ISFACE(NWM,NFM),ISFANO(NWM,NFM),ISFIBR(NWM,NEM,NGRSEGM),
     '  ISGAUS(NWM),ISGRAD(NEM,NGRSEGM),
     '  ISGRID(NWM),ISHIST(0:NPM),ISINCR(NWM),
     '  ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),ISMAP(NGRSEGM),
     '  ISMATE(NWM,NEM),ISNONO(NWM,NPM),ISREAC(NWM),ISRULE(NWM),
     '  ISSECT(NGRSEGM),ISSTRA(NEM,NGRSEGM),
     '  ISSTRE(NEM,NGRSEGM),ISSTRM(NEM,NGRSEGM),
     '  ISSURF(NWM,NGRSEGM),ISVELO(NEM,NGRSEGM),
     '  NEELEM(0:NE_R_M,0:NRM),NPNODE(0:NP_R_M,0:NRM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL FULL
!     Local Variables
      INTEGER IW,nb,ne,nf,nh,nocont,nodata,noelem,nofibr,
     '  NOGRAD,noline,NOMAP,nonode,NOSECT,nosg,NOSTRE,NOSTRA,NOSTRM,
     '  NOVELO,np,nr

      CALL ENTERS('OPSEGM',*9999)
      DO nosg=1,NTSG
        WRITE(OP_STRING,'('' CSEG('',I4,'')= '',A,''  ISEG('',I4,'
     '    //''')='',I2)') nosg,CSEG(nosg)(1:60),nosg,ISEG(nosg)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO

      IF(FULL) THEN
        WRITE(OP_STRING,'('' ISAXES(iw): '',6I4)')
     '    (ISAXES(IW),IW=1,6)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISGAUS(iw): '',6I4)')
     '    (ISGAUS(IW),IW=1,6)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISINCR(iw): '',6I4)')
     '    (ISINCR(IW),IW=1,6)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISLINO(iw): '',6I4)')
     '    (ISLINO(IW),IW=1,6)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISREAC(iw): '',6I4)')
     '    (ISREAC(IW),IW=1,6)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISRULE(iw): '',6I4)')
     '    (ISRULE(IW),IW=1,6)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISBASE(nb): '',9I4)')
     '    (ISBASE(nb),nb=1,NBFT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nodata=1,NTDATA
          WRITE(OP_STRING,'('' ISDATA(iw,nodata='',I2,''): '',3I4)')
     '        nodata,(ISDATA(IW,nodata),IW=1,6)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        DO nr=1,NRT
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            WRITE(OP_STRING,'('' ISELNO(iw,ne='',I3,''): '',6I4)')
     '        ne,(ISELNO(IW,ne),IW=1,6)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ISDANO(iw,ne='',I3,''): '',6I4)')
     '        ne,(ISDANO(IW,ne),IW=1,6)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ISDAPR(iw,ne='',I3,''): '',6I4)')
     '        ne,(ISDAPR(IW,ne),IW=1,6)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ISDATR(iw,ne='',I3,''): '',6I4)')
     '        ne,(ISDATR(IW,ne),IW=1,6)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ISMATE(iw,ne='',I3,''): '',6I4)')
     '        ne,(ISMATE(IW,ne),IW=1,6)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
        DO nf=1,NFT
          WRITE(OP_STRING,'('' ISFACE(iw,nf='',I3,''): '',6I4)')
     '      nf,(ISFACE(IW,nf),IW=1,6)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' ISFANO(iw,nf='',I3,''): '',6I4)')
     '      nf,(ISFANO(IW,nf),IW=1,6)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        DO nr=1,NRT
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            WRITE(OP_STRING,'('' ISNONO(iw,np='',I3,''): '',6I4)')
     '        np,(ISNONO(IW,np),IW=1,6)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
        WRITE(OP_STRING,'('' ISGRID(iw=1,6): '',6I4)')
     '    (ISGRID(IW),IW=1,6)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISHIST(0)= '',I4)') ISHIST(0)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nr=1,NRT
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            WRITE(OP_STRING,'('' ISHIST('',I3,'')= '',I4)')ISHIST(np)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
        DO noline=1,NTLINE
          WRITE(OP_STRING,'('' ISLINE(iw,noline='',I2,''): '',6I4)')
     '      noline,(ISLINE(IW,noline),IW=1,6)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        DO nr=1,NRT
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            WRITE(OP_STRING,'('' ISSURF(iw,ne='',I3,''): '',6I4)')
     '        ne,(ISSURF(IW,ne),IW=1,6)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ISCONO(nh,ne='',I3,''): '',5I4)')
     '        ne,(ISCONO(nh,ne),nh=1,NHM)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nocont=1,NTCONT
              WRITE(OP_STRING,'('' ISCONT(nh,ne='',I3,''nocont='','
     '          //'I1): '',5I4)')
     '          ne,nocont,(ISCONT(nh,ne,nocont),nh=1,NHM)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
            WRITE(OP_STRING,'('' ISGRAD(ne='',I3,'',nograd): '','
     '        //'10I4)') ne,(ISGRAD(ne,NOGRAD),NOGRAD=1,NTGRAD)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ISSTRE(ne='',I3,'',nostre): '','
     '        //'10I4)') ne,(ISSTRE(ne,NOSTRE),NOSTRE=1,NTSTRE)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ISSTRA(ne='',I3,'',nostra): '','
     '        //'10I4)') ne,(ISSTRA(ne,NOSTRA),NOSTRA=1,NTSTRA)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ISSTRM(ne='',I3,'',nostrm): '','
     '        //'10I4)') ne,(ISSTRM(ne,NOSTRM),NOSTRM=1,NTSTRM)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ISVELO(ne='',I3,'',novelo): '','
     '        //'10I4)') ne,(ISVELO(ne,NOVELO),NOVELO=1,NTVELO)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nofibr=1,NTFIBR
              WRITE(OP_STRING,'('' ISFIBR(iw,ne='',I3,'
     '          //''',nofibr='',I2,''): '',6I4)')
     '          ne,nofibr,(ISFIBR(IW,ne,nofibr),IW=1,6)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
C Array misdimensioned !  AJP 15-1-96
C            DO ng=1,NGT(1)
C              WRITE(OP_STRING,'('' ISISOC(iw,ng='',I3,'',ne='',I3,'
C     '          //'''): '',6I4)') ng,ne,(ISISOC(IW,ng,ne),IW=1,6)
C              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C            ENDDO
          ENDDO
        ENDDO
        WRITE(OP_STRING,'('' ISMAP(nomap): '',20I4)')
     '    (ISMAP(NOMAP),NOMAP=1,NTMAP)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ISSECT(nosect): '',20I4)')
     '    (ISSECT(NOSECT),NOSECT=1,NTSECT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('OPSEGM')
      RETURN
 9999 CALL ERRORS('OPSEGM',ERROR)
      CALL EXITS('OPSEGM')
      RETURN 1
      END


