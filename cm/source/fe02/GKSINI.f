      SUBROUTINE GKSINI(ISALIG,ISAXES,ISBASE,ISCLOC,ISCONO,ISCONT,
     '  ISCROS,ISDANO,ISDAPR,ISDATA,ISDATR,ISDIPA,ISDIPO,ISEG,
     '  ISELNO,ISERR,ISFACE,ISFANO,ISFIBR,ISFIEL,ISGAUS,ISGRAD,ISGRID,
     '  ISHIST,ISINCR,ISISOC,ISLEAD,ISLINE,ISLINO,ISL2BE,
     '  ISL3BE,ISMAP,ISMATE,ISNONO,ISN2BE,ISN3BE,ISOBJE,ISPLIN,
     '  ISPLOT,ISPLOTXY,ISPMAR,ISPROF,ISREAC,ISRESI,ISRULE,ISSCAL,
     '  ISSECT,ISSHEE,ISSTRA,ISSTRE,ISSTRM,ISSURF,
     '  ISVELO,CSEG,ERROR,*)

C#### Subroutine: GKSINI
C###  Description:
C###    GKSINI clears workstations.

      IMPLICIT NONE
      INCLUDE 'colo00.cmn'
      INCLUDE 'geom00.cmn'
C      INCLUDE 'cmiss$reference:grou00.cmn'
      INCLUDE 'lead00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER ISALIG(NWM),ISAXES(NWM),ISBASE(99),ISCLOC(NWM),
     '  ISCONO(NHM,NEM),ISCONT(NHM,NEM,NGRSEGM),ISCROS(NWM,NGRSEGM),
     '  ISDANO(NWM,NEM),ISDAPR(NWM,NEM),ISDATA(NWM,NGRSEGM),
     '  ISDATR(NWM,NEM),ISDIPO(NWM,NDIPOLEM,NGRSEGM),
     '  ISDIPA(NWM,NDIPOLEM,NGRSEGM),ISEG(*),
     '  ISELNO(NWM,NEM),ISERR(NWM,NEM),ISFACE(NWM,NFM),ISFANO(NWM,NFM),
     '  ISFIBR(NWM,NEM,NGRSEGM),ISFIEL(NWM,NEM),ISGAUS(NWM,NGM,NEM),
     '  ISGRID(NWM),ISGRAD(NEM,NGRSEGM),
     '  ISHIST(0:NPM),ISINCR(NWM),ISISOC(NWM),
     '  ISLEAD(0:NUMLEADMX),ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),
     '  ISL2BE(NLM),ISL3BE(NLM),ISMAP(NGRSEGM),ISMATE(NWM,NEM),
     '  ISNONO(NWM,NPM),ISN2BE(NLM),ISN3BE(NLM),
     '  ISOBJE(NWM,NGRSEGM,NGRSEGM),
     '  ISPLIN(NWM,NGRSEGM),ISPLOT(NHM,0:NEM,NGRSEGM),ISPLOTXY(2),
     '  ISPMAR(NWM),ISPROF(2),
     '  ISREAC(NWM),ISRESI(NWM),ISRULE(NWM),ISSCAL(NWM,NGRSEGM),
     '  ISSECT(NGRSEGM),ISSHEE(NWM,NEM,NGRSEGM),
     '  ISSTRA(NEM,NGRSEGM),
     '  ISSTRE(NEM,NGRSEGM),ISSTRM(NEM,NGRSEGM),
     '  ISSURF(NWM,NGRSEGM),
     '  ISVELO(NEM,NGRSEGM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER i,iw,ne,nf,ng,nh,ndipole,nl,nlead,nosg,np,nseg,nseg2

      CALL ENTERS('GKSINI',*9999)

C 21/2/97 archived section : cpb 7/2/96 Initialise segment arrays
C cpb 2/4/96 Adding this for the moment
      COLOUR_WS=.TRUE.

      DO iw=1,NWM
        ISALIG(iw)=0
        ISAXES(iw)=0
        ISCLOC(iw)=0
        IF(USE_GRID.NE.0) ISGRID(iw)=0
        ISINCR(iw)=0
        ISISOC(iw)=0
        ISLINO(iw)=0
        ISPMAR(iw)=0
        ISREAC(iw)=0
        ISRESI(iw)=0
        ISRULE(iw)=0
        DO ne=1,NEM
          ISDANO(iw,ne)=0
          ISDAPR(iw,ne)=0
          ISDATR(iw,ne)=0
          ISELNO(iw,ne)=0
          ISERR(iw,ne)=0
          ISFIEL(iw,ne)=0
          ISMATE(iw,ne)=0
          DO nseg=1,NGRSEGM
            ISFIBR(iw,ne,nseg)=0
            ISSHEE(iw,ne,nseg)=0
          ENDDO !nseg
          DO ng=1,NGM
            ISGAUS(iw,ng,ne)=0
          ENDDO !ng
        ENDDO !ne
        DO nf=1,NFM
          ISFACE(iw,nf)=0
          ISFANO(iw,nf)=0
        ENDDO !nf
        DO np=1,NPM
          ISNONO(iw,np)=0
        ENDDO !np
        DO nseg=1,NGRSEGM
          ISCROS(iw,nseg)=0
          IF(USE_DATA.NE.0) ISDATA(iw,nseg)=0
          ISLINE(iw,nseg)=0
          ISLINE(iw,nseg+NGRSEGM)=0 !used for deformed lines
          ISPLIN(iw,nseg)=0
          ISSCAL(iw,nseg)=0
          ISSURF(iw,nseg)=0
          DO ndipole=1,NDIPOLEM*USE_DIPOLE
            ISDIPO(iw,ndipole,nseg)=0
            ISDIPA(iw,ndipole,nseg)=0
          ENDDO !ndipole
          DO nseg2=1,NGRSEGM
            ISOBJE(iw,nseg,nseg)=0
          ENDDO !ne
        ENDDO !nseg
      ENDDO !iw
      ISPROF(1)=0
      ISPROF(2)=0
      DO i=1,99
        ISBASE(i)=0
      ENDDO !i
      DO nlead=0,NUMLEADMX
        ISLEAD(nlead)=0
      ENDDO !nlead
      DO nl=1,NLM
        ISL2BE(nl)=0
        ISL3BE(nl)=0
        ISN2BE(nl)=0
        ISN3BE(nl)=0
      ENDDO !nl
      DO ne=1,NEM
        DO nh=1,NHM
          ISCONO(nh,ne)=0
          DO nseg=1,NGRSEGM
            ISCONT(nh,ne,nseg)=0
          ENDDO !nseg
        ENDDO !nh
      ENDDO !ne
      DO nseg=1,NGRSEGM
        ISMAP(nseg)=0
        ISSECT(nseg)=0
        DO ne=1,NEM
          ISGRAD(ne,nseg)=0
          ISSTRA(ne,nseg)=0
          ISSTRE(ne,nseg)=0
          ISSTRM(ne,nseg)=0
          ISVELO(ne,nseg)=0
        ENDDO !ne
        DO ne=0,NEM
          DO nh=1,NHM
            ISPLOT(nh,ne,nseg)=0
          ENDDO !nh
        ENDDO !ne
      ENDDO !nseg
      DO np=0,NPM
        ISHIST(np)=0
      ENDDO !np
      DO i=1,2
        ISPLOTXY(i)=0
      ENDDO !i
      DO nosg=1,NTSG
        ISEG(nosg)=0
        CSEG(nosg)=' '
      ENDDO

      IMAP=0
      NTSG=0
      NTCHOR=0
      NTCONT=0
      NTCROS=0
      NTDATA=0
      NTFIBR=0
      NTGRAD=0
      NTHIST=0
      NTISOC=0
      NTLINE=0
      NTMAP =0
      NTOBJE=0
      NTPLIN=0
      NTPLOT=0
      NTSCAL=0
      NTSECT=0
      NTSHEE=0
      NTSTRA=0
      NTSTRE=0
      NTSTRM=0
      NTSURF=0
      NTVELO=0

      CALL EXITS('GKSINI')
      RETURN
 9999 CALL ERRORS('GKSINI',ERROR)
      CALL EXITS('GKSINI')
      RETURN 1
      END


