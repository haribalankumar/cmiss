!  voro00.inc
!  This include file is used to index arrays, identify the type of
!  Voronoi node, and identify boundary condition types, and also
!  index the Voronoi face centroids in XNFV

      INTEGER BOUNDARY,INTBOUN,INTERNAL,FAREA,IDIST,TYPE,MAP,BCTYPE,
     '  WALL,INLET,OUTLET,FREESLIP,DRIVING

      PARAMETER(BOUNDARY = 1, INTBOUN = 2,INTERNAL = 3,   FAREA =-1,
     '             IDIST = 0,    TYPE = 1,     MAP = 2,  BCTYPE =-1,
     '              WALL = 0,   INLET = 1,  OUTLET = 2,FREESLIP = 3,
     '           DRIVING = 4)
