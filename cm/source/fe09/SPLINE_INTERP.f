      REAL*8 FUNCTION SPLINE_INTERP(DIRN,START_NODE,NBJ,nn,nelec,
     '  layer,P,ZD,
     '  ne,XID,ni,NPNE,SNODE_NUM,SELEM_NUM)

C#### Function: SPLINE_INTERP
C###  Type: REAL*8
C###  Description: Interpolates between 2 knots given
C###    2nd derivative information about the knots
C###


C*** APRIL 1998 NEED TO TIDY UP
C NOTE: slice = the horiz/vert spline being constructed
C     : nelec = the number of electrodes in a slice -> nknots ?
C     : ne    = LD(nd) the element the dp is in
C     : row   = where the node we want is located
C     : SNODE_NUM = eg 66 for slice
C     : SELEM_NUM   eg 46 for slice

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'

      INTEGER NKNOTSM
      PARAMETER (NKNOTSM=40)       !change to NP_NRM ??

! Parameters
      INTEGER NBJ(NJM,NEM),ne,nelec,ni,nn,NPNE(NNM,NBFM,NEM),
     '  SNODE_NUM,SELEM_NUM,START_NODE
      REAL*8  P(NKNOTSM),XID(NIM,NDM),ZD(NJM,NDM)
      CHARACTER*5 DIRN
! Local Variables
      INTEGER col,offset,row,layer,nn_bottom,tmp1,tmp2
      REAL*8 DIST,xi,xi_dash
! Function
      REAL*8 DATA_DIST


      xi=XID(ni,nn-START_NODE+1)
      xi_dash=1-xi

C*** offset for each ring
      IF(DIRN(1:4).EQ.'VERT') THEN
        offset=NDT                         !data stored in NDT++
C        COL=MOD(ne-SELEM_NUM+1,nelec-1)   !this for plane mesh
        COL=MOD(ne-SELEM_NUM+1,nelec)      !this for closed mesh
        IF(COL.EQ.0) THEN !end column
          COL=nelec-1
        ENDIF

        ROW=INT((ne-SELEM_NUM+1-COL)/NELEC) +1 ! ROW < NSLICE !!!

        tmp2=ROW+offset
        tmp1=tmp2+1

        DIST=DATA_DIST(tmp1,tmp2,ZD)
        SPLINE_INTERP=xi*ZD(NJT+1,tmp1) + xi_dash*ZD(NJT+1,tmp2) +
     '    (
     '    DIST* ((xi**3-xi)*P(COL+1)
     '    +(xi_dash**3-xi_dash)*P(COL))
     '    ) /6.D0

      ELSE ! DIRN.EQ.'HORIZ'
        tmp1=NPNE(1,NBJ(ni,ne),ne)-SNODE_NUM+1 !the node
        nn_bottom=MOD(tmp1,nelec)

        IF(nn_bottom.EQ.0) THEN             !an end node
          ROW=INT(tmp1/nelec)                !finds which row
          nn_bottom=INT(tmp1/ROW)           !the bottom node

          tmp2=layer*nelec
          tmp1=(layer-1)*nelec+1

        ELSE
          tmp2=nn_bottom+(layer-1)*nelec ! +START_NODE-1
          tmp1=tmp2+1
        ENDIF

        DIST=DATA_DIST(tmp1,tmp2,ZD)
        SPLINE_INTERP=xi*ZD(NJT+1,tmp1) + xi_dash*ZD(NJT+1,tmp2) +
     '    (
     '    DIST*
     '    ((xi**3-xi)*P(tmp1)+(xi_dash**3-xi_dash)*P(tmp2))
     '    ) /6.D0

      ENDIF !DIRN

C The original formula
C     '   DIST* ((xi**3-xi)*P(ne-SELEM_NUM+1+1)
C     '   +(xi_dash**3-xi_dash)*P(ne-SELEM_NUM+1))

      RETURN
      END


