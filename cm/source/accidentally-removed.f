C This documentation was accidentally removed at cvs tag cm-split and
C hasn't been replaced yet.

C#### Module: FE00
C###  Description:
C###    Routines for read/write input & output files.
C###  Routine: MATR       (block data)
C###  Routine: ISBINFILEOPEN returns whether a binary file is open
C###  Routine: ISENDBINFILE returns whether a binary file is eof
C###  Routine: BINCLOSEFILE Closes a binary file
C###  Routine: BINOPENFILE Opens a binary file
C###  Routine: BINREADFILE Read data from a binary file
C###  Routine: BINRESETNUMTAG Resets the # of tags in a binary file
C###  Routine: BINSETFILE Sets the position of a binary file
C###  Routine: BINSKIPCMHEADER Skips cmiss binary file header bytes
C###  Routine: BINSKIPFILE Skip bytes from a binary file
C###  Routine: BINWRITEFILE Write data to a binary file
C###  Routine: CHECK_IP_CHAR Check for * or ! in ip files
C###  Routine: DISPLAY_PROMPT displays a prompt with specified format
C###  Routine: GINOUT   handle general input/output
C###  Routine: GINPUT   read general input
C###  Routine: IODATA   handle i/o of field data
C###  Routine: IOELEM   handle i/o of element data
C###  Routine: IOGEOM   handle i/o of f.e. metafile data (IOD files)
C###  Routine: IOHIST   handle i/o of history files
C###  Routine: IO_MATRIX handle i/o of a matrix, inc. sparsity
C###  Routine: IONODE   handle i/o of node data
C###  Routine: IO_VECTOR handle i/o of a vector, inc. sparsity
C###  Routine: IOMATR   handle i/o of matrix files
C###  Routine: IOPOTN   ARCHIVED
C###  Routine: IOSIGN   handle i/o of signal files
C###  Routine: IPGEOM   ARCHIVED
C###  Routine: IOXI     handle i/o of Xi data
C###  Routine: IOXI_NODE handle i/o of Xi data for nodes in a host mesh
C###  Routine: OPEN_FILE opens direct i/o file and read/writes heading
C###  Routine: OPEN_BIN_FILE opens binary file and read/writes heading
C###  Routine: OPEN_SEQ_FILE opens sequen i/o file and r/ws heading
C###  Routine: READ_BIN_TAG_HEADER reads a tag header from a bin file
C###  Routine: READ_LINE read character string from output unit
C###  Routine: RECEIVE_DATA_ASSEMBLE5  recieves results from ZEES proc
C###  Routine: SEND_ARRAYS_ASSEMBLE5  pass global arrays for ZEES proc
C###  Routine: SEND_COMMON_ASSEMBLE5  pass common blocks for ZEES proc
C###  Routine: SEND_DATA_ASSEMBLE5    pass element data  to ZEES proc
C###  Routine: WRITES   write string
C###  Routine: WRITE_BIN_TAG_HEADER reads a tag header from a bin file
C###  Routine: WRITE_CHAR writes a Fortran character to output unit
C###  Routine: WRITE_INT writes an integer to output unit
C###  Routine: WRITE_LINE writes Fortran character and newline to output unit
C###  Routine: WRITE_STRING writes specified # of characters to output unit
C###  Routine: WRITE_STRING_WRAPPER calls WRITE_STRING with same args
C###  Routine: WRITE_LONG writes a long vector to IUNIT
C###  Routine: WRITE_LONG_IDX  writes a long indexed vector to IUNIT
C###  Routine: WRITE_POLYLINE  ARCHIVED
C###  Routine: WRITE_SOL_MATRIX writes solution matrices to a file.


C#### Module: FE01
C###  Description:
C###    General purpose routines

C###  Routine: BDTR       (block data) data for tracing routines
C###  Routine: BDUS       (block data) data for user defined names
C###  Routine: BLK00      (block data) data for common block B00
C###  Routine: BLKIO      (block data) data for i/o common blocks
C###  Routine: ABBREV     (fn) true if 1st argument is abbreviation of 2nd
C###  Routine: ABBRV      (fn) true if 1st argument is abbreviation of 2nd
C###  Routine: AFFIRM     (fn) written to cmiss_achive.f
C###  Routine: BINSEARCH  (fn) written to cmiss_achive.f
C###  Routine: CBBREV     (fn) true if 1st arg contains abbreviation of 2nd
C###  Routine: CELEM      (fn) true if single character arg belongs to set
C###  Routine: CFROMI     (fn) character string from integer variable
C###  Routine: CFROML     (fn) written to cmiss_achive.f
C###  Routine: CFROMR     (fn) character string from REAL*8 variable
C###  Routine: CI1        (fn) written to cmiss_achive.f
C###  Routine: CLOCAT     (fn) returns position of string within string
C###  Routine: CLOWER     (fn) converts string to lower case
C###  Routine: CUPPER     converts string to upper case
C###  Routine: DISTANCE   (fn) computes distance between two points
C###  Routine: EXIST      (fn) written to cmiss_achive.f
C###  Routine: FSOLV3     fast 3x3 system solver
C###  Routine: IDIGITS    (fn) returns number of digits in an integer
C###  Routine: IFROMC     (fn) integer from character
C###  Routine: ILISTLOC   (fn) returns a pointer to an integer in a list
C###  Routine: ILISTMBR   (fn) returns an integer from a list
C###  Routine: ILOG10     (fn) returns truncated log base 10 of REAL*8 arg
C###  Routine: INLIST     (fn) true if 1st arg belongs to list in 2nd arc
C###  Routine: IVALID     (fn) true if arg valid char. represent.n of integer
C###  Routine: LEN_TRIM   (fn) length of string minus trailing blanks
C###  Routine: LEFT_MULT  (fn) left multiplicity of the ith knot for
C###  OSLO Routine: LFROMC     (fn) logical variable from character
C###  string Routine: LOWCAS     (fn) true if argument lowercase
C###  Routine: MADEOF     (fn) true if 1st arg made up of characters in 2nd
C###  Routine: MAXD       (fn) written to cmiss_achive.f
C###  Routine: MAXIM      (fn) written to cmiss_achive.f
C###  Routine: MAXIM2     (fn) written to cmiss_achive.f
C###  Routine: MINIM      (fn) written to cmiss_achive.f
C###  Routine: MINIM2     (fn) written to cmiss_achive.f
C###  Routine: NEGATE     (fn) written to cmiss_achive.f
C###  Routine: NKNU       (fn) written to cmiss_achive.f
C###  Routine: QUAD       (fn) written to cmiss_achive.f
C###  Routine: RIGHT_MULT (fn) right multiplicity of the ith knot for OSLO
C###  Routine: RVALID     (fn) true if arg valid char repres.n of REAL*8 #
C###  Routine: TRIANGLE_AREA (fn) returns area of triangle      
C###  Routine: UPCAS      (fn) true if argument is uppercase
C###  Routine: ADD        add two real-valued strings
C###  Routine: ADD_SPARSE_ELEMENT  Adds a element to the matrix
C###  Routine: ADDSTRTOBUFF  Adds a command string to the buffer
C###  Routine: ALLOCATE_MEMORY Frees (if necessary) and allocates mem
C###  Routine: APPENDC    write a string onto the end of a fortran character
C###  Routine: APPENDCA   write a string onto the end of a character array
C###  Routine: APPENDI    write an integer onto the end of a fortran character
C###  Routine: ASSERT     return error if logical expression false
C###  Routine: ASSIGN     put character string in list of user variables
C###  Routine: CALCUL     calculate arithmetic expression
C###  Routine: CALC_SPARSE calculate the sparsity indicies
C###  Routine: CHECKF     checks a file name is defined
C###  Routine: CHECKQ     checks a command qualifier is member of list
C###  Routine: CMDESTROY  Destroys global CM variables
C###  Routine: CMINITIALISE Initialises global CM variables
C###  Routine: CMS        written to cmiss_achive.f
C###  Routine: CPU_TIMER  calls a 'c' routine to get CPU times
C###  Routine: CRASH      access violation routine
C###  Routine: CROSS      return the vector cross product of two vectors
C###  Routine: CROSS_SEC_AREA_FROM_ELEM calculates elem cross sec contribution
C###  Routine: CROSS_SECTION_AREA calculates cross sec area for surface mesh
C###  Routine: D_CROSS    cross product of two vectors with deriv
C###  Routine: D_NORMALISE normalises a vector and finds its deriv
C###  Routine: D_ROT_COORDSYS orthog rotation of orientations with deriv
C###  Routine: DCPU_TIMER  calls a 'c' routine to get delta CPU times
C###  Routine: DREAL_TIMER calls a 'c' routine to get delta wallclock times
C###  Routine: DEASSI     deassigns a previously defined user variable
C###  Routine: DIMCHK     checks array dimensions
C###  Routine: DIVIDE     divides two real-valued strings
C###  Routine: DO_COMMAND executes an operating system command
C###  Routine: EIGEN1     find e.values and e.vectors for symm matrix
C###  Routine: ELEM_PLANE_INTERSECT find intersection points for fixed xi
C###  Routine: ELEM_PLANE_INTERSECT_XIAXIS find intersect on elem boundary
C###  Routine: ENTERS        traces entry to a subprogram
C###  Routine: ERRORIN    stores and outputs call stack after error
C###  Routine: ERRORS     handles output of error message and call stack
C###  Routine: EVALUE     finds eigenvalues of 2*2 or 3*3 symm matrix
C###  Routine: EVECTR     returns normalised eigenvector of symm matrix
C###  Routine: EXITS      traces exit from a subprogram
C###  Routine: EXPON      raises real-valued string to a string exponent
C###  Routine: FIND       locates first occurrence of a string in a file
C###  Routine: FLAG_ERROR indicate an error, displaying MESSAGE to the user
C###  Routine: FREELU     written to cmiss_achive.f
C###  Routine: FREE_MEMORY frees dynamically allocated memory.
C###  Routine: FREETIMER  written to cmiss_achive.f
C###  Routine: GET_DATE_TIME uses rtl call to return date and time
C###  Routine: GET_SYSTEM    uses rtl call to return system parameters
C###  Routine: GETCHARFROMPTR *REMOVED* Gets a char from pointer+offset
C###  Routine: GETDPFROMPTR *REMOVED* Gets a real*8 from pointer+offset
C###  Routine: GETINTFROMPTR *REMOVED* Gets a integer from pointer+offset
C###  Routine: GETLOGFROMPTR *REMOVED* Gets a logical from pointer+offset
C###  Routine: GET_MEMORY dynamically allocates memory
C###  Routine: GETNEXTCOMFILELINE gets next line from command file
C###  Routine: GETPREVCOMFILELINE gets prev line from command file
C###  Routine: GETRAN     gets transformation matrix from a file
C###  Routine: GETSTR1    rets str typed by user from keybd mapping
C###  Routine: GETTIMER   written to cmiss_achive.f
C###  Routine: GRDXXY     draw grid lines of const x in x,y coords
C###  Routine: GRDYXY     draw grid lines of const y in x,y coords
C###  Routine: GRDXXZ     draw grid lines of const x in x,z coords
C###  Routine: GRDZXZ     draw grid lines of const z in x,z coords
C###  Routine: GRDYYZ     draw grid lines of const y in y,z coords
C###  Routine: GRDRRT     draw grid lines of const r in r,theta coords
C###  Routine: GRDTRT     draw grid lines of const theta in r,theta
C###  Routine: GRDLLM     draw grid lines of const lambda in lambda,mu
C###  Routine: GRDLLT     draw grid lines of con lambda in lambda,theta
C###  Routine: GRDMLM     draw grid lines of con mu in lambda,mu coords
C###  Routine: IHEAPSORT   Numerical Recipes heapsort
C###  Routine: ILIST_COPY  copies a list of integers to another array
C###  Routine: INIT_SPARSE_MATRIX  Initialises the matrix
C###  Routine: INSERT     *REMOVED* replace every occurr. of one string
C###  Routine: INVERT     invert a 2*2 or 3*3 matrix
C###  Routine: ILISTRMDUP ISORT's list and removes duplicates
C###  Routine: ISHELLSORT   Numerical Recipes shell sort
C###  Routine: ISORT      sort integer values into non-decr. sequence
C###  Routine: ISORTP      sort integers into asc. seq. & stores pivots
C###  Routine: LIST_COMMANDS list CMISS commands when use help
C###  Routine: MAINCMLOOP    Main CMISS Loop.
C###  Routine: MAP4       written to cmiss_achive.f
C###  Routine: MAP10      set world coords for history display
C###  Routine: MAQ_LOC    manages auxiliary grid parameter indicies
C###  Routine: MULT       written to cmiss_achive.f
C###  Routine: MULTLY     mult two real-valued strings
C###  Routine: NORMALISE normalises a vector
C###  Routine: NIQ_LOC     return niq# for a particular grid property
C###  Routine: NX_LOC     return nx#  for a particular prob type
C###  Routine: ORTHOG_VEC return a vector orthog to 1 or 2 i/p vectors
C###  Routine: OSLO       Oslo algorithm for refining B-splines
C###  Routine: PARSE      parse string to extract commands & qualifiers
C###  Routine: PARSE_CLASS       parse class
C###  Routine: PARSE_DATA        parse data group
C###  Routine: PARSE_ELEMENTS    parse list of elements
C###  Routine: PARSE_FACES       parse list of faces
C###  Routine: PARSE_GRID        parse list of grid points
C###  Routine: PARSE_LINES       parse list of lines
C###  Routine: PARSE_NODES       parse list of nodes
C###  Routine: PARSE_QUALIFIERS  parse list of qualifiers
C###  Routine: PARSE_REGIONS     parse list of regions
C###  Routine: PARSE_WINDOWS     parse list of windows
C###  Routine: PARSIA     written to cmiss_achive.f
C###  Routine: PARSIL     put item list in string into integer vector
C###  Routine: PARSILG    put list (incl groups) into integer vector
C###  Routine: PARSIN     put item in string into integer variable
C###  Routine: PARSLO     written to cmiss_achive.f
C###  Routine: PARSRA     put nested item list into REAL*8 array
C###  Routine: PARSRE     put item in string into REAL*8 variable
C###  Routine: PARSRL     put list of items in string into REAL*8 vector
C###  Routine: PARSSL     put list of items in string into string vector
C###  Routine: PARSTR     parse string
C###  Routine: PERIODICTASK  handles periodic tasks
C###  Routine: PLANE_LINE_INTERSECT returns intersect of line and plane
c###  Routine: POLAR      eval polar decomposition of a matrix
C###  Routine: POLY       written to cmiss_achive.f
C###  Routine: POST       written to cmiss_achive.f
C###  Routine: PROJAB     calculates the projection of vector B onto A
C###  Routine: PROMPT_REGION_ALL  prompt for region      
C###  Routine: PUTRAN     written to cmiss_achive.f
C###  Routine: QUAD_LINE_INTERSECT calculate intersection of line with quad
C###  Routine: QUIT       perform a graceful end to the plotting session
C###  Routine: RALPHA2D   calcs angle and axis of rotation from R
C###  Routine: RALPHA3D   calcs angle and axis of rotation from R
C###  Routine: READC      reads command file by calling READCOM
C###  Routine: READCOM    reads command file
C###  Routine: REAL_TIMER  calls a 'c' routine to get wallclock times
C###  Routine: RESET      reset transform.n matrix to identity matrix
C###  Routine: RESET_HANDLER reset old error handling routines
C###  Routine: ROOTSC     *** ARCHIVED ***
C###  Routine: ROT_COORDSYS returns orthog rotation of coord system
C###  Routine: ROT_GENERAL returns rotation of vector about arbitrary axis
C###  Routine: ROTATION returns orthogonal rotation of input vector
C###  Routine: RSORT      sorts REAL*8 values into non-decr sequence
C###  Routine: RSORT_REV  sorts REAL*8 values into decreasing sequence
C###  Routine: SCALE0     written to cmiss_achive.
C###  Routine: SCALE1     update transformation matrix by scale factor
C###  Routine: SCALE2     update transformation matrix by scale factor
C###  Routine: SCALE3     update transformation matrix by scale factor
C###  Routine: SET_HANDLER   sets new error handling routines
C###  Routine: SETCMISSCOMMANDPARAMS sets command parameters
C###  Routine: SETEXAMPLEDIR Sets the example directory.
C###  Routine: SET_NUM_THREADS set the number of threads to use
C###  Routine: SETPR      set the prompt string
C###  Routine: SETSTR     reset the prompt string if become corrupted
C###  Routine: SETTRANSMAT sets the transfomation matrix
C###  Routine: SHIFT      update transformation matrix by shift factor
C###  Routine: SHIFT1     update transformation matrix by shift factor
C###  Routine: SHIFT2     update transformation matrix by shift factor
C###  Routine: SHIFT3     update transformation matrix by shift factor
C###  Routine: SHOWCO     print parsed commands and qualifiers
C###  Routine: SLEEPER       sleeps for an integer period of seconds
C###  Routine: SPARSE     accesses sparse array information
C###  Routine: SPLIT      splits string into commands and qualifier
C###  Routine: STACK      written to cmiss_achive.
C###  Routine: STAND      handles standard commands
C###  Routine: STORE_COMMAND adds a command to the command buffer.
C###  Routine: SUBSTR     breaks string into substrings
C###  Routine: SUBTR      string subtraction
C###  Routine: TRACE      turns trace facility on or off
C###  Routine: TRAINV     inverts transformation matrix
C###  Routine: TRAN       transform a 2nd order 3*3 tensor
C###  Routine: TRANSVECT  geometrically transform a vector
C###  Routine: TRIANGLE_LINE_INTERSECT find intersection of line with triangle
C###  Routine: TRIM       find 1st & last non-blank char.s in a string
C###  Routine: TWIST      update transformation matrix by twist factor
C###  Routine: TWIST1     update transformation matrix by twist factor
C###  Routine: TWIST2     update transformation matrix by twist factor
C###  Routine: TWIST3     update transformation matrix by twist factor
C###  Routine: UNKNOW     gives error if command ambiguous or unknown
C###  Routine: USER       substitutes user defined names
C###  Routine: XZ         transform curvilinear coords to cartesian
C###  Routine: XZ_DERIV   transform curvilin derivs to cartesian derivs
C###  Routine: WS_LIST    return list of workstations for seg. creation
C###  Routine: ZERO_SPARSE_MATRIX  *** ARCHIVED  ***
C###  Routine: ZOOM       archived 18-Feb-99
C###  Routine: ZX         transform cartesian coords to curvilinear
C###  Routine: ZZ         apply transform.n matrix to cartesian coords




C#### Variable: IWKDEF(0:noiw)
C###  Type: INTEGER
C###  Set_up: WS_LIST
C###  Description:
C###    IWKDEF(0) is the number of defined workstations.
C###    IWKDEF(noiw),noiw=1,IWKDEF(0) is the list of defined
C###    workstations.

C#### Variable: IWKG(iw)
C###  Type: INTEGER
C###  Set_up: WS_LIST
C###  Description:
C###    IWKG(iw) is 0 for nongraphics window (eg menu), and 1 for
C###    graphics output window.

C#### Variable: IWKS(iw)
C###  Type: INTEGER
C###  Set_up: WS_LIST
C###  Description:
C###    IWKS(iw) is 0 when workstation iw is not defined; 1 when
C###    workstation iw is defined but not active; and 2 when workstation
C###    iw is defined and active.

C#### Variable: IWKT(iw)
C###  Type: INTEGER
C###  Set_up: WS_LIST
C###  Description:
C###    IWKT(iw) is 1 for GX or GKS workstation; 2 for PHIGS
C###    workstation; and 3 for Frame-grabber display screen.

C#### Variable: MAQ_LIST(maq)
C###  Type: INTEGER
C###  Set_up: MAQ_LOC
C###  Description:
C###    MAQ_LIST(0) is the number of maq's defined.
C###    MAQ_LIST(maqq=1..MAQ_LIST(0)) is the actual maq number at
C###    position maqq.

C#### Variable: MAQ_LOCKS(maq)
C###  Type: INTEGER
C###  Set_up: MAQ_LOC
C###  Description:
C###    MAQ_LOCKS(maq) is 1,2 if maq is/is not locked.

C#### Variable: MAQ_SUBTYPE(maq)
C###  Type: INTEGER
C###  Set_up: MAQ_LOC
C###  Description:
C###    MAQ_SUBTYPE(maq) defines the subtype of maq requested, a list
C###    is given in the description for MAQ_LOC

C#### Variable: MAQ_TYPE(maq)
C###  Type: INTEGER
C###  Set_up: MAQ_LOC
C###  Description:
C###    MAQ_TYPE(maq) defines the type of maq requested, a list is
C###    given in the description for MAQ_LOC

C#### Variable: NIQ_LIST(niq)
C###  Type: INTEGER
C###  Set_up: NIQ_LOC
C###  Description:
C###    NIQ_LIST(0) is the number of niq'sdefined.
C###    NIQ_LIST(niqq=1..NIQ_LIST(0)) is the actual niq number at
C###    position niqq.

C#### Variable: NIQ_LOCKS(niq)
C###  Type: INTEGER
C###  Set_up: NIQ_LOC
C###  Description:
C###    NIQ_LOCKS(niq) is 1,2 if niq is/is not locked.

C#### Variable: NIQ_SUBTYPE(niq)
C###  Type: INTEGER
C###  Set_up: NIQ_LOC
C###  Description:
C###    NIQ_SUBTYPE(niq) defines the subtype of niq requested, a list
C###    is given in the description for NIQ_LOC

C#### Variable: NIQ_TYPE(niq)
C###  Type: INTEGER
C###  Set_up: NIQ_LOC
C###  Description:
C###    NIQ_TYPE(niq) defines the type of niq requested, a list is
C###    given in the description for NIQ_LOC

C#### Variable: NRLIST(0:nr)
C###  Type: INTEGER
C###  Set_up: PARSE_REGIONS
C###  Description:
C###    NRLIST(0:nr) is a temporary region list.

C#### Variable: NXLIST(0:nx)
C###  Type: INTEGER
C###  Set_up: PARSE_CLASS
C###  Description:
C###    NXLIST(0:nx) is a temporary region list. The number of entries
C###    in the list are given in NXLIST(0). The list is created by
C###    parsing input from the command line. The default is one class
C###    assigned to be class 1. This is different from NX_LIST.
C###  See-Also: NX_LIST.

C#### Variable: NX_LIST(0:nx)
C###  Type: INTEGER
C###  Set_up: NX_LOC
C###  Description:
C###    NX_LIST(0) is the number of nx problem types defined.
C###    NX_LIST(nxx=1..NX_LIST(0)) is the actual nx number at
C###    position nxx.

C#### Variable: NX_CLASS(nx)
C###  Type: INTEGER
C###  Set_up: NX_LOC
C###  Description:
C###    NX_CLASS(nx) is the nxc class number for problem type nx.

C#### Variable: NX_LOCKS(nx)
C###  Type: INTEGER
C###  Set_up: NX_LOC
C###  Description:
C###    NX_LOCKS(nx) is 1,2 if nx is/is not locked.

C#### Variable: NX_TYPE(nx)
C###  Type: INTEGER
C###  Set_up: NX_LOC
C###  Description:
C###    NX_TYPE(nx) is 1,2,3 for fit/optimisation/solution for problem
C###    type nx.

C#### Variable: XGRC(nj,nu)
C###  Type: REAL*8
C###  Set_up: XZ_DERIV
C###  Description:
C###    Rectagular cartesian version of XG(nj,nu)
C###  See-Also: XG

C#### Module: FE02
C###  Description:
C###    Basic finite element routines.

C###  Routine: CAL_NUM_SITE (fn) calculates number of daughter vessels
C###  Routine: EXISTN     (fn) true if world coords lie inside node segment
C###  Routine: GETNYR     (fn) returns the corresponding rhs ny for a lhs ny
C###  Routine: ISATCOLLPASE (fn) returns if local node is collapsed
C###  Routine: ISSECTOR (fn) *REMOVED* returns if element is a sector element
C###  Routine: LIADJ      (fn) find lines adjacent to current line
C###  Routine: NLATUNIQUE (fn) find number of distinct lattice gp in xi
C###  Routine: ALLOCATE_LATTICE allocates memory for lattice grids
C###  Routine: ANGLE_CHECK *REMOVED* checks angle coronary vessels
C###  Routine: ANGSCA   calc global scale factors from angle change
C###  Routine: ARCDER   calc arc length derivatives wrt Xi
C###  Routine: ARCLEN   calc arc lengths DL of global line segments
C###  Routine: ARCSCA   calc arc lengths DL of global line segments
C###  Routine: AREA     calc area DF of face segments
C###  Routine: BASE_XI  calc base vectors in Xi coords at given Xi
C###  Routine: BIL_XI_FNDR determine if pt is inside bilinear element
C###  Routine: CALC_dydNu calc def derivs of vessel coords wrt nu coords
C###  Routine: CALC_ELEM_SHAR_VERT Calculates lattice element sharing arrays
C###  Routine: CALC_FACE_BASIS_DEP Calculates dep var face bases.
C###  Routine: CALC_FACE_INFORMATION_DEP Calculates face information.
C###  Routine: CALC_FACE_INFORMATION_IND Calculates face information.
C###  Routine: CALC_GLOBAL_NODES Updates global list of nodes from reg
C###  Routine: CALC_LATTICE_MAP Calculate the mappings for lattice grids
C###  Routine: CALC_LATTICE_WEIGHTS form the vector of grid point weights
C###  Routine: CALC_LATTICE_XIQ Calc the XIQ array for lattice grids
C###  Routine: CALC_NENP calc the list of elements surrounding a
C###  np Routine: CALC_NENP_1D calc the elements surrounding a np for 1D
C###  ne Routine: CALC_NH_LOC sets up the NH_LOC array.
C###  Routine: CALC_NHP calc the NHP array.
C###  Routine: CALC_NP_XI calculate the xi pos and elem num of a node
C###  Routine: CALC_NVNEP *** ARCHIVED ***
C###  Routine: CALC_NKH calc the NKH array
C###  Routine: CALC_NKH_FIT calc NKH for fitting problems
C###  Routine: CALC_NP_XI calc xi positions of nodes
C###  Routine: CALC_NUNK calc the NUNK mapping array
C###  Routine: CALC_NY_GRID_DEP calc the NYNQ & NQNY mapping arrays
C###  Routine: CALC_NY_MAPS_DEP calc the NYNP & NPNY mapping arrays
C###  Routine: CALC_NY_MAPS_DEP2 calc the NYNP & NPNY mapping arrays
C###  Routine: CALC_NY_MAPS_IND calc the NYNP & NPNY mapping arrays
C###  Routine: CALC_STRIPE_PTS *** ARCHIVED ***
C###  Routine: CALC_SPARSE_FIT calc sparsity patterns for fitting
C###  Routine: CALC_SPARSE_GKK calc sparsity pattern for GKK
C###  Routine: CALC_SPARSE_GKK_1DTREE calc sparsity pattern for tree
C###  Routine: CALC_SPARSE_SOLVE calc sparsity patterns for solving
C###  Routine: CALC_SPARSE_SOLVE_1DTREE calc sparsity patterns for tree
C###  Routine: CALC_VERSIONS_DEP calc versions for dep vars
C###  Routine: CalculateNLQ set NLQ=1 when local gradient exceeds threshold
C###  Routine: ConstructNAQ construct NAQ from the fine grid NXQ matrix
C###  Routine: ConstructNLQ construct NLQ interpolating connections
C###  Routine: ConstructNXQ construct NXQ connectivity array
C###  Routine: COORD    transform coords from one system to another
C###  Routine: COPY_AND_INIT_INT copies data from one array to another
C###  Routine: CPCG     transfer material params to Gauss point array
C###  Routine: CPCGF    transfer material params to face Gauss point array
C###  Routine: CPCP2    *** ARCHIVED ***
C###  Routine: CPXI     interpolate material params at Xi position
C###  Routine: CREATE_LATTICE sets up the lattice grid scheme
C###  Routine: D_FIBRE_REF_VECS dirn of fibre ref axes with deriv
C###  Routine: D_MAT_VEC dirn cosines of undef matl vects with deriv
C###  Routine: D_MAT_VEC_ROTATE rotates fib ref vecs and deriv into matl vectors
C###  Routine: DEFMGRADRC eval defm gradient tensor wrt rc coords
C###  Routine: DELAUNAY_FACE computes whether two adjoined simplices are delaunay
C###  Routine: DERIV_INFO adjust local derivative numbers
C###  Routine: DLSE     evaluate element scale factors from line info
C###  Routine: DLZJDX   eval covar derivs of deformed theta coords
C###  Routine: DSTATS   calc data point book-keeping variables
C###  Routine: DSTATS_FACE calc data points info based on global face nos.
C###  Routine: DXIDL    calc arc lengths of global line segments
C###  Routine: DXIDXM   calc derivs of Xi wrt undef Nu/Wall coords
C###  Routine: DXIDZM   calc derivs of Xi wrt def nu/wall crds & inverse
C###  Routine: DXRCDX   calc derivs of r.c. crds wrt reference crds
C###  Routine: DXRCDXI  calc 1st derivs of r.c. crds wrt Xi crds
C###  Routine: D2XRCDXI calc 1st and 2nd derivs of r.c. crds wrt Xi crds
C###  Routine: EXIT_FACES *REMOVED* calc faces a coronary element passes though
C###  Routine: EXT_FACEFIT_LIST determines external faces in a volume fit.
C###  Routine: FACCAL   calc face parameters (incl area)
C###  Routine: FACE_INT_LGE  calcs multipliers for face integral terms
C###  Routine: FACE_INT_PREP find info for a face integral
C###  Routine: FACSEG   calc global face topology
C###  Routine: FEMINI   initialise finite element variables & arrays
C###  Routine: FIBRE_REF_VECS calc dirn of fibre ref axes at ng or xi
C###  Routine: FIBRE_REF_VECS_DEF calc def fibre ref axes at ng or xi
C###  Routine: FIND_NEW_ELEM *REMOVED* calculates adjacent ne for coronary mesh
C###  Routine: FIND_FACE_BASIS calculates the face basis number.
C###  Routine: FIND_LATT_NEIJK find the element that holds nlat.
C###  Routine: GET_ELEMENT_VERT finds the vertices of an element.      
C###  Routine: GET_TNVECTOR calc the tangent/normal vector at a node pt
C###  Routine: GET_TNVECTOR2 calc tang./norm. vect. at a given xi pt
C###  Routine: GETEQVNONY finds an ny with the same no
C###  Routine: GKSINI   initialize graphical arrays
C###  Routine: GLOBALC  calc various mapping arrays for coupled probs
C###  Routine: GLOBALF  calc various mapping arrays for fitting
C###  Routine: GLOBALH  calc various mappings arrays from elem arrays
C###  Routine: GLOBALJ  calc various global arrays from elem arrays
C###  Routine: GLOBALJ_1D *REMOVED* calc various arrays from 1D elem arrays
C###  Routine: GLOBALO  calc various global arrays from elem arrays
C###  Routine: GLOBAL_LUNG  calc mappings arrays for lung model
C###  Routine: GSUPPORT generates the support for a lattice based nq
C###  Routine: INTERFACE find regions that a node belongs to
C###  Routine: INIT_NJ_LOC initialise NJ_LOC for a region
C###  Routine: LINCAL   calc global line parameters
C###  Routine: LINCAL_1D   calc global line parameters for 1D mesh only
C###  Routine: LINSCA   calc global line lengths and scale factors
C###  Routine: LINSEG   define global line parameters
C###  Routine: LATTICE_NWQ Calculate NWQ for lattice grid points
C###  Routine: MAT_VEC calc dirn cosines of undef matl vects at ng or xi
C###  Routine: MAT_VEC_DEF calc deformed material vectors at ng or xi
C###  Routine: MAT_VEC_NG calc undeformed material vectors at ng
C###  Routine: MAT_VEC_ROTATE rotates fib ref vecs into matl vectors
C###  Routine: MAT_VEC_XI calc undeformed material vectors at xi
C###  Routine: MELGE    calc #element variables for solution
C###  Routine: MELGEF   calc #element variables for fitting
C###  Routine: MELGEF_FACE   determines variables for face fitting
C###  Routine: MELGEG   calc #element variables for grid
C###  Routine: NENXI    find elements surrounding element ne
C###  Routine: NENXI_1D find elements surrounding 1D element ne
C###  Routine: NEW_POINTS *REMOVED* calc new xi position for coronary mesh
C###  Routine: NODE_CHANGE  CMGUI link - modify the node value/structure
C###  Routine: NODE_CHANGE_DYNAM  CMGUI link - modify the node value/structure
C###  Routine: NODE_CREATE  CMGUI link - create the node
C###  Routine: NODE_CREATE_DYNAM  CMGUI link - create the node
C###  Routine: NODE_DESTROY CMGUI link - destroy the node
C###  Routine: NODE_DESTROY_DYNAM CMGUI link - destroy the node
C###  Routine: OBJ_BEZIER  define object from GKS bezier curve
C###  Routine: OBJ_BOX     define object from GKS box
C###  Routine: OBJ_CIRCLE  define object from GKS circle
C###  Routine: OBJ_LINE    define object from GKS line
C###  Routine: OBJ_ELLIPSE define object from GKS ellipse
C###  Routine: OBJ_STROKE  define object from GKS stroke
C###  Routine: OPESTFMAT output element stiff matrix and/or vectors
C###  Routine: OPSTFMAT  output global stiffness matrix and/or vectors
C###  Routine: PSCOORD1 prolate sph. coords from map
C###  Routine: PSCOORD2 Xi coords from prolate sph. coords on map
C###  Routine: CM_RANGE find maximum range of field or dep. variable
C###  Routine: SEDL     line scale factors from element scale factors
C###  Routine: SIMPS    *** ARCHIVED ***
C###  Routine: SIZE_LATT_ARRAYS sets array sizes for lattice grids
C###  Routine: SOLV2    solve a 2*2 system of linear equations
C###  Routine: SPLINE   Numerical Recipes Version of Spline
C###  Routine: TOFFEL   eval cpts of Christoffel symbol of 2nd kind
C###  Routine: VERT_LATT_INDICES find the lattice indices for element vertices
C###  Routine: VOLUME   calculated the volume enclosed by an element gp
C###  Routine: WALL_VEC calc undef (cardiac) wall axes at ng or xi
C###  Routine: WALL_VEC_DEF calc deformed wall axes at ng or xi
C###  Routine: XCOORD   eval Xi-coords & elem no. of a r.c. point
C###  Routine: XECP     *REMOVED* transform elem params into control pt array
C###  Routine: XECURV   calc curvature at Xi
C###  Routine: XEXG     Gauss point array from elem node array
C###  Routine: XEXW     interp independ. var. array XE at Xi position
C###  Routine: XFXG     Gauss point array from face node array
C###  Routine: XGDXI    Copies derivs from XG to DXJDXI & inverts
C###  Routine: XGMG     metric tensor arrays from Gauss pt array
C###  Routine: XGYG     Gauss array of metric tensors from geom. and material
C###  Routine: XGYGF    face Gauss array of metric tensors from geom. and mat.
C###  Routine: XIXL     return mapped world coords at given Xi coords
C###  Routine: XMG      eval metric tensor arrays at arbitrary Xi pt
C###  Routine: XPXE     transfer global parameters to element params
C###  Routine: XPXF     transfer global parameters to face params
C###  Routine: XPXL     create world coord polyline array from XP
C###  Routine: XQXE     transfer global nq params to element params
C###  Routine: XQXE_REF transfer global ref. nq params to element params
C###  Routine: XYCOORD  find Xi coords of a pt in a bilinear element
C###  Routine: YPZP     transfer global vector to global node params
C###  Routine: YQZQE    transfer grid dep. vars to local elem params
C###  Routine: ZDZDL    put data parameters into element arrays
C###  Routine: ZDZDFACE puts data parameters into face arrays in face fitting
C###  Routine: ZEDS     *REMOVED* converts ZE from Xi to dS derivatives
C###  Routine: ZEEX50   strain cpts wrt fibre/ref coords
C###  Routine: ZEEX51   strain cpts wrt fibre/ref coords (short)
C###  Routine: ZEZG     eval Gauss point array from element array
C###  Routine: ZEZW     interp depend. var. array ZE at Xi position
C###  Routine: ZGMG     eval cpts of metric tensor in deformed state
C###  Routine: ZPXL     create world coord polyline array from ZP
C###  Routine: ZPYP     transfer global node params to global vector
C###  Routine: ZPZE     transfer global node params to element array
C###  Routine: ZPZF     transfer global node params to face array

C#### Variable: ANAL_CHOICE(nr)
C###  Type: INTEGER
C###  Set-up: IPANA3,IPANA9
C###  Description:
C###    <HTML><P>
C###    ANAL_CHOICE(nr) is the type of analytic function used in region
C###    nr. It's value dependends on the type of equation and NJT.</P>
C###    <P>For Laplaces equation:
C###    <PRE>
C###    NJT = 2; ANAL_CHOICE(nr) =
C###      1 : K(x-y)
C###      2 : K(x^2-y^2)
C###      3 : K(x^2+2xy-y^2)
C###      4 : K(a.r^n.cos(n.t)+b.r^m.sin(m.t)+c.r.cos(t)+d.r.sin(t)+e)
C###      5 : Centre dipole (single circle)
C###      6 : Centre dipole (multiple circles)
C###      7 : Eccentric dipole (single circle)
C###      8 : Eccentric dipole (multiple circles)
C###      9 : Anisotropic annulus
C###     10 : Anisotropic plate (f(z)=e^z)
C###     11 : Anisotropic plate (f(z)=a.sin(z)+b.cos(z))
C###    NJT=3; ANAL_CHOICE(nr) =
C###      1 : K(x-y)
C###      2 : K(z)
C###      3 : K(x^2-y^2)
C###      4 : K(x^2-z^2)
C###      5 : K(y^2-z^2)
C###      6 : K(x^2+y^2-2z^2)
C###      7 : K(x^2-2y^2+z^2)
C###      8 : K(-2x^2+y^2+z^2)
C###      9 : K(x^2+2xy-y^2)
C###     10 : K(x^2+2xz-z^2)
C###     11 : K(y^2+2yz-z^2)
C###     12 : Centre dipole (single sphere)
C###     13 : Centre dipole (multiple spheres)
C###     14 : Eccentric dipole (single sphere)
C###     15 : Eccentric dipole (multiple spheres)
C###    </PRE>
C###    </HTML>

C#### Variable: CONVERG_TOL
C###  Type: REAL*8
C###  Set_up: FEMINI
C###  Description:
C###    <HTML> <PRE>
C###    CONVERG_TOL=DLAMCH('EPS')*5.0d0
C###    where DLAMCH('EPS') is the relative machine precision
C###    CONVERG_TOL is for testing for convergence. Convergence tests
C###           should be of the form
C###
C###                  ABS(X   -  X )
C###                       i+1    i
C###           IF (  -------------------  <  CONVERG_TOL  ) THEN
C###                      1+ABS(X )
C###                             i
C###
C###           or for norms
C###
C###                       ||r||
C###                   ------------  < CONVERG_TOL
C###                   SQRT(n)+||b||
C###
C###    </PRE> </HTML>
C###  See-Also: LOOSE_TOL,ZERO_TOL

C#### Variable: ISC_GK(ny)
C###  Type: INTEGER
C###  Set_up: CALC_SPARSE_SOLVE,CALC_SPARSE_SOLVE_1DTREE
C###  Description:
C###    ISC_GK stores column information for sparsity pattern of GK
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: ISC_GKK(ny,nx)
C###  Type: INTEGER
C###  Set_up: CALC_SPARSE_GKK
C###  Description:
C###    ISC_GKK stores column information for sparsity pattern of GKK
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: ISC_GQ(ny)
C###  Type: INTEGER
C###  Set_up: CALC_SPARSE_SOLVE
C###  Description:
C###    ISC_GQ stores column information for sparsity pattern of GQ
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: ISR_GK(ny)
C###  Type: INTEGER
C###  Set_up: CALC_SPARSE_SOLVE,CALC_SPARSE_SOLVE_1DTREE
C###  Description:
C###    ISR_GK stores row information for sparsity pattern of GK
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: ISR_GKK(no,nx)
C###  Type: INTEGER
C###  Set_up: CALC_SPARSE_GKK
C###  Description:
C###    ISR_GKK stores row information for sparsity pattern of GKK
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: ISR_GQ(ny)
C###  Type: INTEGER
C###  Set_up: CALC_SPARSE_SOLVE
C###  Description:
C###    ISR_GQ stores row information for sparsity pattern of GQ
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: JTYP1
C###  Type: INTEGER
C###  Set_up: FEMINI
C###  Description:
C###    JTYP1 is 1,2 for elements defined by user/chosen from menu.

C#### Variable: JTYP3
C###  Type: INTEGER
C###  Set_up: FEMINI
C###  Description:
C###    JTYP3 is 1..5 for coordinates rectangular cartesian/cylindrical
C###    polar/spherical polar/prolate spheroidal/oblate spheroidal.

C#### Variable: JTYP5
C###  Type: INTEGER
C###  Set_up: FEMINI
C###  Description:
C###    JTYP5 is 1,2 for basis functions in Lagrange or
C###    Hermite/monomial form.

C#### Variable: JTYP7
C###  Type: INTEGER
C###  Set_up: FEMINI
C###  Description:
C###    JTYP7 is 1..5 for dependent variables rectangular cartesian/
C###    cylindrical polar/spherical polar/prolate spheroidal/oblate
C###    spheroidal.

C#### Variable: JTYP8
C###  Type: INTEGER
C###  Set_up: FEMINI
C###  Description:
C###    JTYP8 is 0,1,2 for no output/output/family output of basis
C###    functions.

C#### Variable: LDR(0:nd)
C###  Type: INTEGER
C###  Set_up: GLOBALO
C###  Description:
C###    LDR(0:nd) is the list of data points in the reduced system.
C###    ie data points that have been projected onto an element.

C#### Variable: LGE(nhs,nrc)
C###  Type: INTEGER
C###  Set_up: MELGE,MELGEF
C###  Description:
C###    LGE(nhs,nrc) is location of the element variables nhs
C###    (=1,NHST(nrc)) in global system. nrc=1 gives the row location,
C###    nrc=2 gives the column location.

C#### Variable: LGKE(-nk:nhs,nkk,2)
C###  Type: INTEGER
C###  Set_up: MELGEB
C###  Description:
C###    LGKE contains the ny and nz numbers for a BE stiffness matrix
C###    for the GK array. The first index specifies the variable in
C###    BE stiffness matrix. Negative first indicies specify the
C###    variables of the singular node e.g. -1 specifies the acutal
C###    singular node variable, -2 specifies the first derivative of
C###    the singular node etc. Positive first indicies specify the
C###    local variables (nhs) of the element currently being
C###    integrated as for the FEM case. The second index specifies
C###    the row for the columns e.g. a value of 1 specifies the row
C###    of the current singular node, a value of 2 specifies the row
C###    of the first derivative of the current singular node.
C###    The last index is for seperating nys and nzs.
C###    If it is 1 then the array will contain the variable ny
C###    numbers (and the first index is zero the row ny number).
C###    If it is 2 then the array will contain the nz numbers
C###    corresponding to those ny numbers. For example LGKE(-1,1,2)
C###    will contain the nz number of the element for the current
C###    singular row node and current singular node column.

C#### Variable: LGQE(-nk:nhs,nkk,2)
C###  Type: INTEGER
C###  Set_up: MELGEB
C###  Description:
C###    LGQE contains the ny and nz numbers for a BE stiffness matrix
C###    for the GQ array. The first index specifies the variable in
C###    BE stiffness matrix. Negative first indicies specify the
C###    variables of the singular node e.g. -1 specifies the acutal
C###    singular node variable, -2 specifies the first derivative of
C###    the singular node etc. Positive first indicies specify the
C###    local variables (nhs) of the element currently being
C###    integrated as for the FEM case. The second index specifies
C###    the row for the columns e.g. a value of 1 specifies the row
C###    of the current singular node, a value of 2 specifies the row
C###    of the first derivative of the current singular node.
C###    The last index is for seperating nys and nzs.
C###    If it is 1 then the array will contain the variable ny
C###    numbers (and the first index is zero the row ny number).
C###    If it is 2 then the array will contain the nz numbers
C###    corresponding to those ny numbers. For example LGQE(-1,1,2)
C###    will contain the nz number of the element for the current
C###    singular row node and current singular node column.

C#### Variable: LOOSE_TOL
C###  Type: REAL*8
C###  Set_up: FEMINI
C###  Description:
C###    <HTML> <PRE>
C###    LOOSE_TOL=DSQRT(DLAMCH('EPS'))
C###    where DLAMCH('EPS') is the relative machine precision
C###    LOOSE_TOL is to be used in the same manner as CONVERG_TOL
C###    when a looser criterion is desired
C###    </PRE> </HTML>
C###  See-Also: CONVERG_TOL,ZERO_TOL


C#### Variable: NBJF(nj,nf)
C###  Type: INTEGER
C###  Set_up: FACSEG
C###  Description:
C###    NBJF is the basis number for geometric variable nj on
C###    face nf.

C#### Variable: NEL(0:nelm,nl)
C###  Type: INTEGER
C###  Set_up: LINSEG
C###  Description:
C###    NEL(0,nl) is the number of elements adjoining line nl.
C###    NEL(1..,nl) are the element numbers of elements adjoining
C###    line nl.

C#### Variable: NENP(np,0:nep,0:nr)
C###  Type: INTEGER
C###  Set_up: CALC_NENP
C###  Description:
C###    NENP(np,0:nep,nr) is the reverse mapping of NPNE for region nr.
C###    NENP(np,0,nr) is the number of elements surrounding the global
C###    node np for region nr. NENP(np,1..,nr) is the list of elements
C###    surrounding global node np for region nr. For nr=0 the total
C###    list of elements is given.

C#### Variable: NFE(nb)
C###  Type: INTEGER
C###  Set_up: FACSEG
C###  Description:
C###    NFE(nb) is the number of faces for element basis type nb.

C#### Variable: NFF(nfe,ne)
C###  Type: INTEGER
C###  Set_up: FACSEG
C###  Description:
C###    NFF(nfe,ne) are global face numbers corresponding to side nfe
C###    of element ne.

C#### Variable: NFFACE(0:noface,nr)
C###  Type: INTEGER
C###  Set_up: FACSEG
C###  Description:
C###    NFFACE(0,nr) is the number of faces in region nr.
C###    NFFACE(noface,nr), noface=1..NFFACE(0,nr) are the face numbers
C###    in region nr.  (NFT is the total number of faces in all regions.)
C###  See-Also: NEELEM,NFFACE,NLLINE,NPNODE

C#### Variable: NFT
C###  Type: INTEGER
C###  Set_up: FEMINI
C###  Description:
C###    NFT is the total number of global element faces.

C#### Variable: NH_LOC(0:nhx,0:nx)
C###  Type: INTEGER
C###  Set_up: CALC_NH_LOC
C###  Description:
C###    NH_LOC(0:nhx,0:nx) is the nh number for dependent variable nhx
C###    of problem type nx.

C#### Variable: NH_TYPE(nh,2)
C###  Type: INTEGER
C###  Set_up: CALC_NH_LOC
C###  Description:
C###    NH_TYPE(nh,1..2) is the nhx/nx number for nh.

C#### Variable: NHP(np,nr,nx)
C###  Type: INTEGER
C###  Set_up: CALC_NHP
C###  Description:
C###    NHP(np,nr,nx) is the number of dependent variables defined at
C###    node np for problem type nx.

C#### Variable: NKEF(0:4,nn,nfe,nb)
C###  Type: INTEGER
C###  Set_up: FACSEG
C###  Description:
C###    NKEF(i,nn,nfe,nb) stores the element derivatives in a face for
C###    face node nn of face nfe and basis nb.  NKEF(0,nn,nfe,nb) is the
C###    number of element derivatives for face node nn of face nfe and
C###    basis nb. NKEF(1..,nn,nfe,nb) is the list of element derivatives
C###    for face node nn of face nfe and basis nb.

C#### Variable: NKH(nh,np,nc,nr)
C###  Type: INTEGER
C###  Set_up: CALC_NKH
C###  Description:
C###    NKH(nh,np,nc,nr) is the number of nodal derivs of dependent
C###    variable nh at node np.  nc=2 is its normal derivative (on the
C###    lowest corner element if np is at a corner), nc=3 is the normal
C###    derivative on the next highest corner element, if appropriate,
C###    and nc=4 is the normal derivative on the highest corner element
C###    if appropriate.

C#### Variable: NKJ(nj,np)
C###  Type: INTEGER
C###  Set_up: GLOBALJ
C###  Description:
C###    NKJ(nj,np) is number of derivatives for geometric variable nj
C###    at np.

C#### Variable: NLE(nb)
C###  Type: INTEGER
C###  Set_up: LINSEG
C###  Description:
C###    NLE(nb) is the number of element line segments for basis nb.

C#### Variable: NLF(naf,nf)
C###  Type: INTEGER
C###  Set_up: FACSEG
C###  Description:
C###    NLF(naf,nf) are the global line numbers corresponding to local
C###    arc NAF of face nf.

C#### Variable: NLL(nae,ne)
C###  Type: INTEGER
C###  Set_up: LINSEG
C###  Description:
C###    NLL(nae,ne) are the global line numbers corresponding to
C###    local arc NAE of ne.
C###  See-Also: NEELEM,NFFACE,NLLINE,NPNODE

C#### Variable: NLLINE(0:nl,0:nr)
C###  Type: INTEGER
C###  Set_up: LINSEG
C###  Description:
C###    NLLINE(0,nr) is the number of lines in region nr.
C###    NLLINE(l,nr), l=1..NLLINE(0,nr) are global nl line numbers
C###    in region nr.

C#### Variable: NLNO(no,nx)
C###  Type: INTEGER
C###  Set_up: GLOBALO
C###  Description:
C###    NLNO(no,nx) is the line number nl for optimisation variable no.

C#### Variable: NLT
C###  Type: INTEGER
C###  SET_UP: FEMINI
C###  Description:
C###    NLT is the number of global line segments.

C#### Variable: NMNO(1:2,noopti)
C###  Type: INTEGER
C###  Set_up: IPOPTI,GLOBALO
C###  Description:
C###    <HTML> <PRE>
C###    NMNO(1,0) is the total number of different material parameters
C###    (of a given material law) in the fit (i.e. the number of
C###    different nm indices of CE/CP that are used). Setup in IPOPTI.
C###    NMNO(1,noopti),noopti=1,NTOPTI are the material parameters
C###    (of a given material law) in the fit (i.e. the nm numbers in
C###    the CE/CP array) given optimisation variable number, noopti.
C###    NMNO(2,noopti),noopti=1,NTOPTI are the ne/np/ng numbers to which
C###    the optimisation variable noopti (a material param) is attached.
C###    e.g. If one is fitting nodally based material parameters,
C###    using a material law containing 4 material parameters, and only
C###    the second and third values are being fitted (at each node) then
C###    NMNO(1,0)  = 2 (fitting 2 of the 4 material law params).
C###    NMNO(1,no) = 2,3,2,3,2,3,2,3,2,3 ... (nm numbers in CP)
C###    NMNO(2,no) = 1,1,2,2,3,3,4,4,5,5 ... (node numbers)
C###    </PRE> </HTML>
C###  See-Also: NONM

C#### Variable: NONM(nm,nindex)
C###  Type: INTEGER
C###  Set_up: GLOBALO
C###  Description:
C###    NONM(nm,nindex) returns the no number of the material parameter
C###    nm in the fit. nindex ranges over elements, noelem
C###    (spatially/piecewise constant), or nodes, nonode (p'wise lin).
C###  See-Also: NMNO.

C#### Variable: NNF(0:17,nfe,nb)
C###  Type: INTEGER
C###  Set_up: FACSEG
C###  Description:
C###    NNF(0,nfe,nb) is the number of element nodes in face nfe of
C###    element with basis nb.  NNF(1,nfe,nb) is Xi direction normal to
C###    face.  NNF(2..17,nfe,nb) are the element node numbers in face.

C#### Variable: NNL(0:4,nae,nb)
C###  Type: INTEGER
C###  Set_up: LINSEG
C###  Description:
C###    NNL(0,nae,nb) is the number of nodes along arc nae of
C###    element with basis nb. NNL(1..4,nae,nb) are the
C###    element/face node numbers along local arc NAE.

C#### Variable:NONL(nl,nx)
C###  Type: INTEGER
C###  Set_up: GLOBALO
C###  Description:
C###    NONL(nl,nx) is the optimisation variable no for line
C###    variable nl.

C#### Variable: NPF(9,nf)
C###  Type: INTEGER
C###  Set_up: FACSEG
C###  Description:
C###    <HTML> <PRE>
C###    NPF(1,nf)    is the 1st Xi-direction of face segment nf
C###        2        is the basis function type for 1st Xi-direction
C###                   (1,2,3 or 4)
C###        3        is the 2nd Xi-direction
C###        4        is the basis function type for 2nd Xi-direction
C###        5        is the number of elements adjoining face
C###                   (1,2 for external,internal)
C###        6        is the first  element number joined to face
C###        7        is the second element number joined to face
C###        8        is the local face number of the first element
C###                   joined to face
C###        9        is the local face number of the second element
C###                   joined to face
C###
C###  For Simplex elements NPF is:
C###
C###    NPF(1,nf)    is 1
C###        2        is the basis function type of the face (Simplex)
C###        3        is 1
C###        4        is the basis function type of the face (Simplex)
C###        5..9     as above
C###
C###    </PRE> </HTML>


C#### Variable: NPL(5,0:3,nl)
C###  Type: INTEGER
C###  Set_up: LINSEG
C###  Description:
C###    <HTML> <PRE>
C###    The first index is 'i'
C###    The second index is nj, ie 0,nj=1,3

C###    nj =0  i=1 - is the Xi-direction of line segment nl
C###           i=2 - is the line number adjacent in the -XI direction
C###           i=3 - is the line number adjacent in the +XI direction
C###           i=4 - If Jtyp2B=1 NPL(4,0,nl) has special meaning (see bottom of ARCLEN) (line mapping)
C###           i=5 - spare
C###    nj =1+ i=1 - are the types of basis function (1,2,3 for linear, quadratic, cubic Lagrange, 4 for cubic Hermite or 6,7 for a quadratic Hermite (node 1/2))
C###           i=2..5 - are the global nodes along line nl in the direction of Xi
C###                    If line has any nj with Hermite basis:
C###                    i=2,3 are the global nodes along line nl in the direction of Xi
C###                    i=4,5 are nk #s (at each end)of 1st derivs wrt Xi
C###    </PRE> </HTML>


C#### Variable: NPNE(nn,nbf,ne)
C###  Type: INTEGER
C###  Set_up: IPELEM,IPMESH
C###  Description:
C###    NPNE(nn,nbf,ne) are global node numbers for basis nbf of
C###    element ne.
C###  See-Also: NNT

C#### Variable: NPNF(nn,nbf)
C###  Type: INTEGER
C###  Set_up: CALC_FACE_INFORMATION_DEP,CALC_FACE_INFORMATION_IND
C###  Description:
C###   NPNF(nn,nbf) are global node numbers for basis nbf.

C#### Variable: NQNY(2,nyq,0:nrc,nx)
C###  Type: INTEGER
C###  Set_up: CALC_NY_GRID_DEP
C###  Description:
C###   NQNY is the reverse mapping to NYNQ (for global nyq values).
C###   NQNY(1,nyq,0:nrc,nx)=nh and NQNY(2,ny,0:nrc,nx)=nq

C#### Variable: NYNQ(nh,nq,0:nrc,nx)
C###  Type: INTEGER
C###  Set_up: CALC_NY_GRID_DEP
C###  Description:
C###   NYNQ for problem nx is the mapping between
C###   the mesh grid points and mesh dofs ny

C#### Variable: NP_INTERFACE(0:np,0:i)
C###  Type: INTEGER
C###  Set_up: INTERFACE
C###  Description:
C###    <HTML> <PRE>
C###    NP_INTERFACE(np,0) = number of regions sharing node np (.ge.1)
C###    NP_INTERFACE(np,1) = nr1 if node is defined in region nr1
C###    NP_INTERFACE(np,2) = nr2 if node is defined in region nr2
C###    NP_INTERFACE(np,3) = nr3 if node is defined in region nr3
C###                         (nr1 < nr2 < nr3).
C###    OR
C###    NP_INTERFACE(0,1..3) is the number of node triplets on aerofoil.
C###    NP_INTERFACE(np,1..3) are the node numbers on the aerofoil.
C###    </PRE> </HTML>

C#### Variable: NTCNTR
C###  Type: INTEGER
C###  Set_up: GLOBALO
C###  Description:
C###    NTCNTR is the number of non-linear constraints.

C#### Variable: NTOPTI
C###  Type: INTEGER
C###  Set_up: GLOBALO
C###  Description:
C###    NTOPTI is the number of optimising variables.

C#### Variable: NVHE(nn,nbf,nh,ne)
C###  Type: INTEGER
C###  Set_up: CALC_VERSIONS_DEP
C###  Description:
C###    NVHE(nn,nbf,nh,ne) is the global version number used for
C###    dependent variable nh of element vertex nn, family basis nbf in
C###    element ne. For example: different versions for multiple thetas
C###    or corners for BEM problems.

C#### Variable: NVHF(nn,nbf,nh)
C###  Type: INTEGER
C###  Set_up: CALC_FACE_INFORMATION_DEP
C###  Description:
C###    NVHF(nn,nbf,nh) is the global version number used for
C###    dependent variable nh of face vertex nn, family basis nbf.

CC AJPs 191297
C#### Variable: NVHP(nh,np,nc,0:nr)
C###  Type: INTEGER
C###  Set_up: CALC_VERSIONS_DEP,IPFIT
C###  Description:
C###    NVHP(nh,np,nc,nr) is the number of versions of dependent
C###    variable nh at node np.
C###    For example: different versions for multiple
C###    thetas coordinates at the apex or corners for BEM problems.
CC AJPe

C#### Variable: NYQNR(0:ny,0:nrc,nc,0:nr,nx)
C###  Type: INTEGER
C###  Set_up: CALC_NY_GRID_DEP
C###  Description:
C###    NYQNR(0:ny,0:nrc,nc,0:nr,nx) is the list of grid
C###    nys for a region.
C###    NYQNR(0:ny,0,nc,0:nr,nx) gives the global list of nys for a
C###    region.  NYQNR(0:ny,1,nc,0:nr,nx) gives the list of local
C###    [=global in all cases] rows for region nr and matrix nc.
C###    NYQNR(0:ny,0,nc,0:nr,nx) gives the list of local column numbers
C###    for region nr and matrix nc. NYQNR(0,nrc,nc,nr,nx) is the number
C###    in the list and NYQNR(1..,nrc,nc,nr,nx) is the list of nys.
C###    NYQNR(0:ny,nrc,nc,0:nr,nx) gives the list of local row and
C###    column numbers for the matrix nc in the region nr.  nr=0 gives
C###    the mappings for the entire coupled problem.
C###    NOTE: This is not set up until after the solution mapping
C###    arrays have been set up.

C#### Variable: NYQT(nrc,nc,nx)
C###  Type: INTEGER
C###  Set_up: CALC_NY_GRID_DEP
C###  Description:
C###    NYQT(nrc,nc,nx) is the total number of grid ny values for
C###    rows (nrc=1) and columns (nrc=2) for matrix nc and problem
C###    type nx.

C#### Variable: NZT(nc,nx)
C###  Type: INTEGER
C###  Set_up:
C###  Description:
C###    NZT(nc,nx) is the number of entries in packed 1D global arrays
C###    for nc = 1..4 and problem nx.

C#### Variable: A_VECTOR
C###  Type: REAL*8
C###  Set_up: MAT_VEC_NG
C###  Description:
C###    A_VECTOR is fibre angle vector (in sheet).

C#### Variable: B_VECTOR
C###  Type: REAL*8
C###  Set_up: MAT_VEC_NG
C###  Description:
C###    B_VECTOR is sheet angle vector (in sheet orthogonal to fibres).

C#### Variable: C_VECTOR
C###  Type: REAL*8
C###  Set_up: MAT_VEC_NG
C###  Description:
C###    C_VECTOR is normal to sheet.

C#### Variable: DF(nf)
C###  Type: REAL*8
C###  Set_up: AREA
C###  Description:
C###    DF(nf) is the area of face segment nf.

C#### Variable: DL(i,nl)
C###  Type: REAL*8
C###  Set_up: LINCAL
C###  Description:
C###    DL(1..2,nl) are the scale factors used for
C###    derivatives in line segment nl.
C###    DL(3,nl) is arc-length of line segment nl.

C#### Variable: GZ
C###  Type: REAL*8
C###  Set_up: ZGMG
C###  Description:
C###    GZ is the determinant of GZL.

C#### Variable: GZL
C###  Type: REAL*8
C###  Set_up: ZGMG
C###  Description:
C###    GZL is the covariant component of the metric tensor.

C#### Variable: GZU
C###  Type: REAL*8
C###  Set_up: ZGMG
C###  Description:
C###    GZU is the contravariant component of the metric tensor.

C#### Variable: PAOPTY(noopti)
C###  Type: REAL*8
C###  Set_up: GLOBALO
C###  Description:
C###    PAOPTY(noopti) is 1..3 for nodal/line/material parameters to be
C###    optimised. For activation optimisations PAOPTY is 1..3 for
C###    mesh activation times/transmembrane jump/wave function width
C###    parameters to be optimised.

C#### Variable: RG(ng)
C###  Type: REAL*8
C###  Set_up: XGMG
C###  Description:
C###    RG(ng) is the Jacobian of element coordinate transformation at
C###    Gauss point ng.

C#### Variable: XE(ns,nj)
C###  Type: REAL*8
C###  Set_up: XPXE
C###  Description:
C###    XE(ns,nj) is the Xj geometric position or derivative for
C###    element dofs ns, coordinate ns.

C#### Variable: XL(nj,nodx)
C###  Type: REAL*8
C###  Set_up: ZPXL
C###  Description:
C###    XL(nj,nodx),nodx=1,ntdx is the world coordinate polyline array.

C#### Variable: XW(nj,nu)
C###  Type: REAL*8
C###  Set_up: XEXW
C###  Description:
C###    XW(nj,nu) are element dependent variables interpolated
C###    at an Xi location.

C#### Variable: ZA(na,nh,nc.ne)
C###  Type: REAL*8
C###  Set_up: YPZP
C###  Description:
C###    ZA(na,nh,nc,ne) is element based dependent variable.

C#### Variable: ZA1(na,nh,nc.ne)
C###  Type: REAL*8
C###  Set_up: YPZP
C###  Description:
C###    ZA1 is a copy of ZA
C###  See-Also: ZA

C#### Variable: ZE(ns,nhx)
C###  Type: REAL*8
C###  Set_up: ZPZE
C###  Description:
C###    ZE(ns,nhx) are element dependent variables.

C#### Variable: ZE1(ns,nhx)
C###  Type: REAL*8
C###  Set_up: ZPZE
C###  Description:
C###    ZE1 is a copy of ZE
C###  See-Also: ZE

C#### Variable: ZERO_TOL
C###  Type: REAL*8
C###  Set_up: FEMINI
C###  Description:
C###    <HTML> <PRE>
C###    ZERO_TOL=DLAMCH('EPS')*5.0d0
C###    where DLAMCH('EPS') is the relative machine precision
C###    ZERO_TOL is for testing real values
C###         eg. IF(DABS(X).GT.ZERO_TOL) THEN
C###    </PRE> </HTML>
C###  See-Also: CONVERG_TOL,LOOSE_TOL

C#### Variable: ZG(nhx,nu)
C###  Type: INTEGER
C###  Set_up: ZEZG
C###  Description:
C###    <HTML> <PRE>
C###    ZG(nhx,1) is value of dependent variable nhx at current Gauss
C###    point.
C###    ZG(nhx,2) are 1st derivatives wrt Xi material coords if JP=0
C###           4                       or nu or Theta    "   "  JP=1
C###           7
C###    </PRE> </HTML>

C#### Variable: ZP(nk,nv,nh,np,nc)
C###  Type: REAL*8
C###  Set_up: YPZP
C###  Description:
C###    ZP(nk,nv,nh,np,nc) are the dependent variables at global node
C###    np.

C#### Variable: ZP1(nk,nv,nh,np,nc)
C###  Type: REAL*8
C###  Set_up: YPZP
C###  Description:
C###    ZP1 is a copy of ZP
C###  See-Also: ZP

C#### Variable: ZW(nhx,nu)
C###  Type: REAL*8
C###  Set_up: ZEZW
C###  Description:
C###    ZW(nhx,nu) are element dependent variables interpolated
C###    at an Xi location


C#### Module:  FE03
C###  Description:
C###    Data fitting routines.

C NOTE: Real*8 variables not correctly set in nonlinear part
C Also check that WK1 is correctly dimensioned for nonlinear solution.

C###  Routine: NYPJK    (fn) *** ARCHIVED ***
C###  Routine: NYF      (fn) *** ARCHIVED ***
C###  Routine: NYPHK    (fn) *** ARCHIVED ***
C###  Routine: CLOS11   find closest Xi pt for 1D rect.cart elements
C###  Routine: CLOS12   written to cmiss_archive.f
C###  Routine: CLOS13   written to cmiss_archive.f
C###  Routine: CLOS14   written to cmiss_archive.f
C###  Routine: CLOS21   find closest Xi pt for 2D rect.cart elements
C###  Routine: CLOS22   written to cmiss_archive.f
C###  Routine: CLOS23   written to cmiss_archive.f
C###  Routine: CLOS24   find closest Xi pt for 2D prolate elements
C###  Routine: CLOS31   find closest Xi pt for 3D rect.cart elements
C###  Routine: CLOS3D   *REMOVED* find closest Xi pt for 3D curvilinear els.
C###  Routine: CONT_PLANE_X written to cmiss_archive.f
C###  Routine: FACORR    Fibre angle correction for sloping surface
C###  Routine: FITFLD   fit field variables to data points
C###  Routine: FIT_PATCH fit field vars to data points in a patch
C###  Routine: FITDEF   fits field variables with F smoothing
C###  Routine: FITMAT   fit material parameters to residuals
C###  Routine: FITFOU   *** ARCHIVED ***
C###  Routine: FITGAU   *** ARCHIVED ***
C###  Routine: FITGEO   *** ARCHIVED ***
C###  Routine: FITSIG   fits signal files
C###  Routine: FITSIG_SPLINE  fits signals using splines
C###  Routine: IOHESS   *** ARCHIVED ***
C###  Routine: IPDATA   *** ARCHIVED ***
C###  Routine: LDLDR    Calculate element data points for a data nr
C###  Routine: LNLXI    *** ARCHIVED ***
C###  Routine: NEWXID   *** ARCHIVED ***
C###  Routine: NL_NO    calculate optimisation dofs no from lines
C###  Routine: NWXID    calculate Xi coords of data point projection
C###  Routine: NY_NO    *** ARCHIVED ***
C###  Routine: PROJ21   find closest Xi proj pt for 2D r.c. elements
C###  Routine: PROJ_ORTHOG find closest point in element to given posn
C###  Routine: PROJ_FUNCT   *** ARCHIVED ***
C###  Routine: PROJ_HESS    *** ARCHIVED ***
C###  Routine: PROJ_MONIT   *** ARCHIVED ***
C###  Routine: SOBOLEV  evaluates Sobolev smoothing value and Jacobian
C###  Routine: XPES     *** ARCHIVED ***
C###  Routine: YGER     eval element load vector for least sqs fit
C###  Routine: YGES     eval element stiff matrix for least sqs fit
C###  Routine: ZDER     eval element load vector for least sqs fit
C###  Routine: ZDERF    *** ARCHIVED ***
C###  Routine: ZDES     eval element stiff matrix for least sqs fit
C###  Routine: ZDESF    *** ARCHIVED ***

C#### Input of data for parameter optimisation.
C###  NDT            is  total number of data points
C###  ZD(nj,nd)       "  rect. cart. coords of data point nd
C###  WD(nj,nd)       "  weighting factor for        "
C###  XID(ni,nd)      "  Xi-coordinate of            "
C###  SQ(nd)          "  square of dist from mesh to "
C###  LD(nd)          "  line or face l no. assoc. with data point nd
C###  NXI(ni,1,ne)    "  element number of first adjacent to element ne
C###  NDLT(ne)        "  no. data points within element ne
C###  NDDL(nl,nde)    "  global data pt no. of line data pt nde
C###  ZDL(nh,nde)     "  rect. cart. coords of       "
C###  WDL(nh,nde)     "  weighting factor for        "
C###  XIDL(ni,nde)    "  Xi-coordinate of            "
C###  NOT(nrc,nc,nr,nx)  "  #solution d.o.f.s
C###  NONY(  0,ny,nrc,nr,nx)  # solution d.o.f.s coupled to mesh var ny
C###  NONY(noy,ny,nrc,nr,nx) solution d.o.f. of var noy for mesh var ny
C###  CONY(noy,ny,nrc,nr,nx) coefficient for     "         "          "
C###  NYNO(  0,no,nrc,nr,nx)  # mesh d.o.f.s coupled to solution var no
C###  NYNO(nyo,no,nrc,nr,nx) mesh d.o.f. of var nyo for solution var no
C###  CYNO(nyo,no,nrc,nr,nx) coefficient for     "         "          "
C###  XO(no)          "  optimisation variable
C###  SS              "  total sum of squared distances SQ(nd)
C###  LN(0)           "  number of elements in fitting
C###  LN(l)           "  element number for l=1..LN(0)
C###  NPO(0)          "  number of global nodes in data fitting
C###  NPO(n)          "  global node number corr to n=1..NPO(0)

C#### Variable: NDDATA(nodata,nr)
C###  Type: INTEGER
C###  Set_up: IPDATA,DEDATA,IOSIGN
C###  Description:
C###    <HTML><UL>
C###    <LI>
C###      NDDATA(0,nr) is the total number of data points in region nr.
C###    <LI>
C###      NDDATA(nodata,nr), nodata=1..NDDATA(0,nr) are the data point
C###      numbers.
C###    <LI>
C###      NDDATA(0,0) is the total number of data points defined
C###      in all regions.
C###    </UL></HTML>
C###  See-Also: ZD,WD,NDT

C#### Module: FE04
C###  Description:
C###    Routines for mesh subdivision

C###  Routine: DIVH1 Subdivide Lagrange/Hermite tensor product element
C###  Routine: DIVS1 Subdivide special Hermite simplex elements
C###  Routine: DIVS2 Subdivide sector elements
C###  Routine: CALC_NNIP Calculates the local node numbers
C###  Routine: REFINE_FINDNODE Determines if a node exists at an xi location
C###  Routine: REFINE_SETNODE set nodal properties for new node
C###  Routine: REFINE_SETSCALEFACTORS adjust scale factors


C cpb 4/6/97 Rewritting DIVH1 - old DIVH1 *** ARCHIVED ***

C#### Module: FE05
C###  Description:
C###    Routines for basis functions.

C###  Routine: BASIS1   params for Lagrange/Hermite tensor prod basis
C###  Routine: BASIS2   params for Simplex/Serendipity/Sector basis
C###  Routine: BASIS3   params for B-spline tensor product basis
C###  Routine: BASIS4   params for Fourier basis
C###  Routine: BASIS5   params for Boundary Element Lagrange schemes
C###  Routine: BASIS6   params for Boundary Element Simplex schemes
C###  Routine: BASIS7   params for extended Lagrange basis
C###  Routine: BASIS8   params for auxilliary basis functions:lagrange
C###  Routine: CALC_CONTRIB_COEFF calc nodal contribution coefficients
C###  Routine: CALCPGG  *REMOVED* calc the integrated tensor prod of basis fns
C###  Routine: GAUSS1   Gauss coords & wgts for tensor prod elems
C###  Routine: GAUSS2   Gauss coords & wgts for simplex elements
C###  Routine: GAUSS2_HERMITE Hermite simplex elements
C###  Routine: GAUSS3   *** ARCHIVED ***
C###  Routine: GAUSS4   Gauss coords & wgts for tensor prod elems
C###  Routine: GAUSS5   set up BE tensor prod basis function info
C###  Routine: GAUSS6   set up BE hermite simplex basis fn info
C###  Routine: GAUSS7   Gauss coords & wgts for extended lagrange
C###  Routine: GAUSS8   Gauss coords & wgts for auxiliary elements
C###  Routine: GAUSS9   Gauss coords & wgts for fourier elements
C###  Routine: GAUSS10  define basis fn for BE problems
C###  Routine: GAUSS11  define basis fn for adaptive Telles scheme
C###  Routine: GAUSS12  Gauss coords & wgts for simplex/sector BEM elems
C###  Routine: GAUSSLEG Gauss coords & weights,Gauss-Legendre quadrature
C###  Routine: GAUSSLOG calc logarithmic gauss points and weights
C###  Routine: GAUSSPWB *** ARCHIVED ***

C#### Variable: NABTYP(nbf)
C###  Type: INTEGER
C###  Set_up: BASIS8
C###  Description:
C###    NABTYP is 1..3 for Legendre/Fourier/pressure auxiliary basis
C###    type for basis nbf.

C#### Variable: NAN(ni,na,nbf)
C###  Type: INTEGER
C###  Set_up: BASIS8
C###  Description:
C###    NAN is the polynomial degree in Xi(ni) direction for
C###    auxiliary variable na of basis nbf.

C#### Variable: NAT(nbf)
C###  Type: INTEGER
C###  Set_up:
C###  Description:
C###    NAT is the number of auxiliary or spline basis functions.

C#### Variable: NBCD(nbf)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NBCD(nbf) indicates whether or not cross derivatives for basis nbf
C###    are set to zero.  If NBCD(nbf) is 1, nbf is a Fergusson
C###    basis (Hermite lacking cross derivatives).  Otherwise NBCD(nbf) is 0.

C#### Variable: NBI(nbf)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    <HTML><PRE>
C###    NBI(nbf) is the index for the scaling type for basis nbf.
C###    It is
C###      = 1 for unit scale factors
C###      = 2 for element scale factors read in
C###      = 3 for global scale factors read in
C###      = 4 for scale factors calculated from angle change
C###      = 5 for scale factors calculated from arc length
C###      = 6 for scale factors calc'd from arithmetic mean of arc length
C###      = 7 for scale factors calc'd from harmonic mean of arc length
C###    </PRE></HTML>

C#### Variable: NBT
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NBT is the number of basis function types.

C#### Variable: NBFT
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NBFT is the number of basis function families.

C#### Variable: NGT(nb)
C###  Type: INTEGER
C###  Set_up: BASIS1
C###  Description:
C###    NGT is the number of Gauss points for basis nb.

C#### Variable: NIT(nbf)
C###  Type: INTEGER
C###  Set_up:
C###  Description:
C###    NIT(nbf) is the number of local Xi-coordinates for basis nbf.

C#### Variable: NKT(0:nn,nbf)
C###  Type: INTEGER
C###  Set_up: BASIS*
C###  Description:
C###    NKT(0,nbf) is maximum number of nodal derivatives or polynomial
C###    terms for splines for basis nbf. NKT(nn,nbf) is the number of
C###    nodal derivatives at node nn for basis nbf.

C#### Variable: NNT(nbf)
C###  Type: INTEGER
C###  Set_up: BASIS*
C###  Description:
C###    NNT(nbf) is the number of element nodes for basis nbf.

C#### Variable: NST(nbf)
C###  Type: INTEGER
C###  Set_up: BASIS*
C###  Description:
C###    NST(nbf) is the number of element parameters for basis
C###    nbf (= sum over nn of NKT(nn,nbf)).

C#### Variable: NUT(nbf)
C###  Type: INTEGER
C###  Set_up: BASIS1,2,3,4,5,6,7,8
C###  Description:
C###    NUT(nbf) is the number of Xi-coord derivatives for basis nbf.


C#### Module: FE06
C###  Description:
C###    Routines used in plotting.

C###  Routine: CFUNC    return field value and derivatives
C###  Routine: CONTOR   plot contours of a field
C###  Routine: CONTOR1  fast version for bilinear field
C###  Routine: CONTOR2  fast version for bicubic  field
C###  Routine: CROOTS   calc Xi-coords where contour crosses elem bdry
C###  Routine: CROOTS1  fast version for bilinear field
C###  Routine: CROOTS2  fast version for bicubic  field
C###  Routine: CROSSSECTION draws sagittal crosssection of heart
C###  Routine: FIBRE1   draw measured fibre field on curvilin plane
C###  Routine: FIBRE2   draw fitted fibre axis
C###  Routine: FIBRE3   draw fitted fibre sheet
C###  Routine: FIBRE4   draw fitted fibre helix
C###  Routine: GRADIENT draw gradient vectors of scalar field
C###  Routine: LDAT     draw lines from data pts to elem surface
C###  Routine: SHEET1   draw fitted sheet angles
C###  Routine: SHEET2   draw measured sheet angles (uncorrected)
C###  Routine: SHEET3   *** ARCHIVED ***
C###  Routine: STRAIN2  draw principal strains
C###  Routine: STRESS1  draw principal stresses


C#### Module: FE07
C###  Description:
C###    Routines to solve fem and bem equations

C###  Routine: A_COEFF     (fn) return coefficient of Uxx
C###  Routine: ALFA        (fn) return slope of 1st characteristic
C###  Routine: BETA        (fn) return slope of 2nd characteristic
C###  Routine: B_COEFF     (fn) return coefficient of Uxt
C###  Routine: C_COEFF     (fn) return coefficient of Utt
C###  Routine: F_COEFF     (fn) return 'F' in A Uxx + B Uxt + C Utt = F
C###  Routine: CALC_LAPL calculates Laplacian matrix
C###  Routine: CALC_SAMPLE_FROM_TIME calcs sample # from times
C###  Routine: CALC_TIME_FROM_SAMPLE calcs time from sample #s
C###  Routine: CALC_FORWARD_ACTN calculates the bsp from an actn field
C###  Routine: ASSEMBLE1 Assem matrices for Static linear FEM probs
C###  Routine: ASSEMBLE2 Assem matrices for Static linear BEM probs
C###  Routine: ASSEMBLE3 Assem matrices for Time dependent FEM probs
C###  Routine: ASSEMBLE4 Assem matrices for Time dependent BEM probs
C###  Routine: ASSEMBLE5 Assem matrices for Nonlinear FEM problems
C###  Routine: ASSEMBLE5_FACE Assem matrices for FEM face integrals
C###  Routine: ASSEMBLE6 Assem matrices for Modal analysis probs
C###  Routine: ASSEMBLE7 Assem matrices for Fourier transform probs
C###  Routine: ASSEMBLE8 *** ARCHIVED ***
C###  Routine: ASSEMBLE9 Assem matrices for Finite Difference problems
C###  Routine: ASSEMBLE10 Assem matricies for implicit ionic problems
C###  Routine: ASSEMBLE10_FE Assem matricies for implicit ionic problems
C###  Routine: ASSEMBLE11 Assem matrices for pulmonary transport
C###  Routine: BRENT     finds 1D minimum (Numerical Recipes)
C###  Routine: CHARAC    return sol.n at R from parameters at P & Q
C###  Routine: CFUNC2    *REMOVED* define contraint functs (called by MINLSSQP)
C###  Routine: CFUNC5    define contraint functs (called by MINOS)
C###  Routine: CONFUN    return contraints
C###  Routine: D_CONFUN  return derivs of contraints wrt opt params
C###  Routine: D_RESFUN  return derivs of resids wrt opt params
C###  Routine: D_ZPRP    calculates derivs of global residual vector
C###  Routine: FIN_SOL   final output and file closure for pulmonary
C###  Routine: FUNCT1    *** ARCHIVED ***
C###  Routine: FUNCT2    define fn to minimize (called by MINLSSQP)
C###  Routine: FUNCT3    *** ARCHIVED ***
C###  Routine: FUNCT4    define fn to minimize (called by MINLSSQP)
C###  Routine: FUNCT5    define fn to minimize (called by MINOS)
C###  Routine: GENASSEM  Assemble global arrays (GK, GM, etc.)
C###  Routine: GENSOL    call solution routines
C###  Routine: GEOMIN    *** ARCHIVED ***
C###  Routine: LSFUNC    used in nonlinear solution search
C###  Routine: MARCH1    solve parabolic  eqns by time integration
C###  Routine: MARCH2    solve hyperbolic eqns by time integration
C###  Routine: MARCH3    to cmiss_gridarchive.f
C###  Routine: MARCH4    explicit/inplicit finite difference problems
C###  Routine: MARCH5    to cmiss_gridarchive.f
C###  Routine: MARCH6    solve quasi-static problems
C###  Routine: MARCH7    to cmiss_gridarchive.f
C###  Routine: MARCH8    general grid/ionic current solver
C###  Routine: MARCH8_COUPLED general bidomain solver
C###  Routine: MARCH11   solve pulmonary (air/blood) transport
C###  Routine: MARCH19   -> cmiss_archive1
C###  Routine: MARCH60   Fluid dynamics explicit marcher
C###  Routine: MATMOD    matmod subroutine for MINOS
C###  Routine: MG_COLLOCATE collocation with multigrid acceler.n
C###  Routine: MG_INTERPOL  multigrid interpolat.n (prolongation)
C###  Routine: MG_RELAX     multigrid relaxation
C###  Routine: MG_RESIDUAL  multigrid residual
C###  Routine: MG_RESTRICT  multigrid restriction
C###  Routine: MG_SOLVE     solve multigrid coarse grid eqtns
C###  Routine: MINOSOPTI optimisation routine using minos (buffer)
C###  Routine: MINOSOPTIA optimisation routine using minos
C###  Routine: MNBRAK    bracket a minimum in 1D (Numerical Recipes)
C###  Routine: MODAL     *** ARCHIVED ***
C###  Routine: MONIT     *** ARCHIVED ***
C###  Routine: NAGMIN    minimisation routine using MINLSSQP (buffer)
C###  Routine: NAGMINA   minimisation routine using MINLSSQP
C###  Routine: NONLIN    control soln of nonlin,quasi-static problems
C###  Routine: NONLIN_CONT control soln of nonlin contact-mech problems
C###  Routine: OBJFUN    return objective function
C###  Routine: POST_SOL  post-solution update for pulmonary transport
C###  Routine: PRE_SOL   pre-solution update for pulmonary transport
C###  Routine: RESFUN    return residuals
C###  Routine: RESFUN_ACTN      return residuals for inverse problem
C###  Routine: RESFUN_ACTN_CC   return residuals for activation problem with CC
C###  Routine: RESFUN_DIPOLE    return residuals for dipole problem
C###  Routine: SOLVE1    solve symm, pos. def. system of eqtns
C###  Routine: SOLVE2    solves unsymm BEM equations
C###  Routine: SOLVE3    solve unsymm sparse linear equations
C###  Routine: SOLVE5    solve nonlinear FEM problems
C###  Routine: SOLVE6    solve modal analysis problems
C###  Routine: SOLVE8    solve static Laplace by collocation
C###  Routine: SOLVE9    solve coupled FEM-BEM problems
C###  Routine: SOLVE10   solve coupled Finite Difference problems
C###  Routine: SOLVE11   solve pulmonary transport problems
C###  Routine: STARTCON  set i.c.s in method of characteristics
C###  Routine: SURFACELAPLACIAN  calculates the surface laplcian residual
C###  Routine: ZEES      calculate element tangent stiffness matrix
C###  Routine: ZFFS      calculate face tangent stiffness matrix
C###  Routine: ZPGKK     *** renamed ASSEMBLE5 ***
C###  Routine: ZPOP      print solutions and residuals
C###  Routine: ZPRP      calculate global residual vector
C###  Routine: ZPRP_FACE calculate face residual vector component

C#### Variable: BC_POINTS(a,j,nq_list)
C###  Type: REAL*8
C###  Set_up: MARCH4
C###  Description:
C###    nq_list is an index for the list of bifurcation pts
C###    BC_POINTS(1,1,0) is total #bifurcations & terminal pts
C###    BC_POINTS(1,j,nq_list) is 1st segment of bifurcation
C###    BC_POINTS(2,j,nq_list) is 2nd segment of bifurcation
C###    BC_POINTS(3,j,nq_list) is 3rd segment of bifurcation
C###    j is 1 for pts adjacent to bifurcation
C###    j is 2 for pts at one grid point away from bifurcation
C###    j is 3 for location of half-step point information


C#### Variable: D_RP(ny,ny)
C###  Type: REAL*8
C###  Set_up: D_ZPRP
C###  Description:
C###    D_RP(ny,ny) are the derivatives of the global residuals
C###    with respect to material/geometric parameters.

C#### Variable: EIGVAL(nt,2)
C###  Type: REAL*8
C###  Set_up: SOLVE6
C###  Description:
C###    EIGVAL(nt,1) is the real part of eigenvalue nt.
C###    EIGVAL(nt,2) is the complex part of eigenvalue nt.

C#### Variable: EIGVEC(no,nt,2)
C###  Type: REAL*8
C###  Set_up: SOLVE6
C###  Description:
C###    EIGVEC(no=1..,nt,1) are the real eigenvector components for
C###    eigenvalue nt.
C###    EIGVEC(no=1..,nt,2) are the complex eigenvector components for
C###    eigenvalue nt.

C#### Variable: ER(nhs)
C###  Type: REAL*8
C###  Set_up: GENASSEM
C###  Description:
C###    ER(nhs) are the element RHS terms.

C#### Variable: GD(nz)
C###  Type: REAL*8
C###  Set_up: GENASSEM,XPGD
C###  Description:
C###    GD(nz)is the global damping matrix (first order time
C###    derivatives) or domain integrals or source terms for dipoles.

C#### Variable: GK(nz)
C###  Type: REAL*8
C###  Set_up: GENASSEM
C###  Description:
C###    GK(nz) is the global stiffness matrix (zeroth order time
C###    derivatives).

C#### Variable: GKK(nz,nx)
C###  Type: REAL*8
C###  Set_up: GENSOL
C###  Description:
C###    GKK(nz,nx) is the constraint reduced global stiffness matrix
C###    for problem type nx.

C#### Variable: GM(nz)
C###  Type: REAL*8
C###  Set_up: GENASSEM
C###  Description:
C###    GM(nz) is the global mass matrix (second order time
C###    derivatives).

C#### Variable: GMM(nz)
C###  Type: REAL*8
C###  Set_up: MODAL
C###  Description:
C###    GMM(nz) is the constraint reduced global mass matrix.

C#### Variable: GQ(nz)
C###  Type: REAL*8
C###  Set_up: GENASSEM
C###  Description:
C###    GQ(nz) is the right hand side matrix.

C#### Variable: GR(ny)
C###  Type: REAL*8
C###  Set_up: GENASSEMB
C###  Description:
C###    GR(ny) is the component of the unreduced RHS vector at row ny.

C#### Variable: GRR(no)
C###  Type: REAL*8
C###  Set_up: GENSOL
C###  Description:
C###    GRR(no) is the constraint reduced RHS vector.

C#### Variable: ISC_GD(ny)
C###  Type: INTEGER
C###  Set_up: ASSEMBLE3,ASSEMBLE6
C###  Description:
C###    ISC_GD stores column information for sparsity pattern of GD
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: ISC_GM(ny)
C###  Type: INTEGER
C###  Set_up: ASSEMBLE3,ASSEMBLE6
C###  Description:
C###    ISC_GM stores column information for sparsity pattern of GM
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: ISC_GMM(no)
C###  Type: INTEGER
C###  Set_up:
C###  Description:
C###    ISC_GMM stores column information for sparsity pattern of GMM
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: ISR_GD(ny)
C###  Type: INTEGER
C###  Set_up: ASSEMBLE3,ASSEMBLE6
C###  Description:
C###    ISR_GD stores row information for sparsity pattern of GD
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: ISR_GM(ny)
C###  Type: INTEGER
C###  Set_up: ASSEMBLE3,ASSEMBLE6
C###  Description:
C###    ISR_GM stores row information for sparsity pattern of GM
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: ISR_GMM(no)
C###  Type: INTEGER
C###  Set_up: ASSEMBLE3,ASSEMBLE6
C###  Description:
C###    ISR_GMM stores row information for sparsity pattern of GMM
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: LAPL(nytr,nytr)
C###  Type: REAL*8
C###  Set_up: CALC_LAPL
C###  Description:
C###    LAPL is the surface Laplacian (second derivative) used for
C###    regularisation of activation inverse problems.
C###  See-Also: LAPLSQR

C#### Variable: LAPLSQR(nytr,nytr)
C###  Type: REAL*8
C###  Set_up: CALC_LAPL
C###  Description:
C###    LAPLSQR is the square of LAPL (L^t.L).
C###  See-Also: LAPL

C#### Variable: RDF(ns,idoxf,nhx)
C###  Type: REAL*8
C###  Set_up: ZPRP_FACE
C###  Description:
C###    RDF(ns,idoxf,nhx) is the contribution from a face to the
C###    residual corresponding to face dof ns, cross-face derivative
C###    idoxf, and problem dependent variable number nhx.

C#### Variable: RQ(-1:1,-1:1,-1:1)
C###  Type: REAL*8
C###  Set_up: MG_RESTRICT
C###  Description:
C###    RQ are weights for restriction operator.

C#### Variable: WEIGHT(nitb)
C###  Type: REAL*8
C###  Set_up: MG_RESTRICT
C###  Description:
C###    WEIGHT is factor multiplying RQ to ensure unity sum.

C#### Variable: WK1(no,4)
C###  Type: REAL*8
C###  Set_up: GENSOL
C###  Description:
C###    WK1(no,4) is a work array.

C#### Variable: WK2(no,3)
C###  Type: REAL*8
C###  Set_up: GENSOL
C###  Description:
C###    WK2(no,3) is a work array.

C#### Variable: WK3(no)
C###  Type: REAL*8
C###  Set_up: GENSOL
C###  Description:
C###    WK3(no) is a work array.

C#### Variable: WK4(no,7)
C###  Type: REAL*8
C###  Set_up: GENSOL
C###  Description:
C###    WK4(no,7) is a work array.

C#### Variable: XO(no,nx)
C###  Type: REAL*8
C###  Set_up: GENSOL
C###  Description:
C###    XO(no,nx) is the solution vector.

C#### Variable: YP(ny,niy,nx)
C###  Type: REAL*8
C###  Set_up: GENSOL
C###  Description:
C###    <HTML>
C###    YP(ny,niy,nx) is the permenant store of the solution variables.
C###    The niy index determines the type of variable stored. The first
C###    7 variables are fixed for all problems. The niy variables from
C###    8 onwards depend on the problem. The variable types are:
C###    <PRE>
C###    YP(ny,1,nx) is current solution (solution at time T+DT for time
C###                   dependent problems)
C###          2     is incremental solution
C###          3     is initial solution/conditions
C###          4     is residual vector
C###          5     is temporary storage 1
C###          6     is temporary storage 2
C###          7     is analytic solution (if any)
C###    </PRE>
C###    For time dependent problems:
C###    <PRE>
C###    YP(ny,8,nx) is previous solution (solution at time T)
C###          9     is mean predicted displacement
C###         10     is estimate of error
C###         11     is velocity (at time T+DT)
C###         12     is previous velocity (at time T)
C###         13     is mean predicted velocity
C###         14     is acceleration (at time T+DT)
C###         15     is previous acceleration (at time T)
C###         16     is mean predicted acceleration
C###    </PRE>
C###    For nonlinear elasticity problems:
C###    <PRE>
C###    YP(ny,5,nx)  is used to store solution increments to
C###                 add to YP(ny,1);
C###    YP(ny,10,nx) is the current reference solution for cavity
C###                 elements (set up in IPINI5/UPSOLU).
C###    YP(ny,11,nx) is the current reference solution for active
C###                 contraction (set up with the "fem update solution"
C###                 command)
C###    </PRE>
C###    For nonlinear elasticity problems involving contact:
C###    <PRE>
C###    YP(ny,6,nx)  is used to store current estimate of the contact forces.
C###    </PRE>
C###    </HTML>

C#### Variable: ZE1(ns,nh)
C###  Type: REAL*8
C###  Set_up: D_RESFUN
C###  Description:


C#### Module: FE08
C###  Description:
C###    Numerical routines
C###  Routine: ADAMS           Adams-Moulton solver
C###  Routine: AM_DE           Used by Adams-Moulton solver
C###  Routine: AM_INTERPOLATE  Used by Adams-Moulton solver
C###  Routine: AM_STEP         Used by Adams-Moulton solver
C###  Routine: CHECK_FIRST_A   (fn) logical check of factorisation
C###  Routine: EIGENPROBLEM    Solves an eigenproblem
C###  Routine: ESOLVE          Evalues/vectors of a matrix (buffer)
C###  Routine: GETEIGENWORK    *REMOVED*
C###  Routine: IMPROVED_EULER  Improved Euler integrator
C###  Routine: INTEGRATOR integrates ionic currents from t to t+dt
C###  Routine: MP_PSEUDO_INVERSE Moore-Penrose pseude inverse.
C###  Routine: RUNGE_KUTTA     4th order Runge-Kutta integrator

C###  Routine: LSODA_WRAPPER   Wrapper for the lsoda ode integrator

C###  Routine: DSAUPD          ARPACK eigensolver
C###  Routine: DSEUPD          ARPACK eigensolver
C###  Routine: CMAV            ARPACK eigensolver
C###  Routine: DSTATSA         ARPACK eigensolver
C###  Routine: DSAUP2          ARPACK eigensolver
C###  Routine: DSESRT          ARPACK eigensolver
C###  Routine: DSORTR          ARPACK eigensolver
C###  Routine: DSCONV          ARPACK eigensolver
C###  Routine: DGETV0          ARPACK eigensolver
C###  Routine: DSAITR          ARPACK eigensolver
C###  Routine: DSEIGT          ARPACK eigensolver
C###  Routine: DSGETS          ARPACK eigensolver
C###  Routine: DSAPPS          ARPACK eigensolver
C###  Routine: DSTQRB          ARPACK eigensolver

C###  Routine: CURFIT          FITPACK fits a 2-D spline
C###  Routine: FPBACK          FITPACK auxililary routine used by CURFIT
C###  Routine: FPBSPL          FITPACK auxililary routine used by SPLDER
C###  Routine: FPCHEC          FITPACK auxililary routine used by CURFIT
C###  Routine: FPCURF          FITPACK auxililary routine used by CURFIT
C###  Routine: FPDISC          FITPACK auxililary routine used by CURFIT
C###  Routine: FPGIVS          FITPACK auxililary routine used by CURFIT
C###  Routine: FPKNOT          FITPACK auxililary routine used by CURFIT
C###  Routine: FPRATI          FITPACK auxililary routine used by CURFIT
C###  Routine: FPROTA          FITPACK auxililary routine used by CURFIT
C###  Routine: SPLDER          FITPACK evaluates the derivates of a B-spline
C###  Routine: SPLEV           FITPACK auxililary routine used by CURFIT

C###  Routine: BVALUE          PPPACK evaluates the derivates od a B-spline
C###  Routine: INTERV          PPPACK auxililary routine used by BVALUE

*
*     ================================================================
*
C#### Module: FE09
C###  Description:
C###    Non-machine specific real functions.

C###  Routine: EDP      elliptic integral of the second kind - E(m)
C###  Routine: I0       modified bessel function
C###  Routine: I1       modified bessel function
C###  Routine: KDP      elliptic integral of the first kind - K(m)
C###  Routine: K0       modified bessel function
C###  Routine: K1       modified bessel function
C###  Routine: M_PLGNDR_SIN   calculates m*plm/sin
C###  Routine: INT_LATT_JUMP (fn) interpolates the 3 lattice jump functions
C###  Routine: LATT_JUMP1 (fn) lattice jump function 1
C###  Routine: LATT_JUMP2 (fn) lattice jump function 2      
C###  Routine: LATT_JUMP3 (fn) lattice jump function 3   
C###  Routine: PAF      (fn) eval 1D Fourier basis functions
C###  Routine: PAL0     (fn) eval 1D Legendre aux basis fn  (constant)
C###  Routine: PAL1     (fn) eval 1D Legendre aux basis fn  (linear)
C###  Routine: PAL2     (fn) eval 1D Legendre aux basis fns (order>=2)
C###  Routine: PAP4     (fn) eval 1D Pressure aux quartic (hat) bas fn
C###  Routine: PF1      (fn) eval Fourier series basis function
C###  Routine: PFXI     (fn) interpolate nodal array XE at XI
C###  Routine: PGX      (fn) eval 1st deriv of basis fn wrt Xj
C###  Routine: PH2      (fn) eval 1D quadratic Hermite basis function
C###  Routine: PH3      (fn) eval 1D cubic Hermite basis function
C###  Routine: PL1      (fn) eval 1D linear Lagrange basis function
C###  Routine: PL2      (fn) eval 1D quadratic Lagrange basis function
C###  Routine: PL2S1    (fn) eval special 1D quadratic basis fn type 1
C###  Routine: PL2S3    (fn) eval special 1D quadratic basis fn type 3
C###  Routine: PL3      (fn) eval 1D cubic Lagrange basis funtion
C###  Routine: PLGNDR   associated Legendre polynomial
C###  Routine: PLXI     (fn) interpolates nodal values XP in a line segment
C###  Routine: PLXIZ    (fn) interpolates nodal values ZP in a line segment
C###  Routine: PP1      (fn) eval polynomial
C###  Routine: PS1      (fn) *** ARCHIVED ***
C###  Routine: PS2      (fn) *** ARCHIVED ***
C###  Routine: PS3      (fn) *** ARCHIVED ***
C###  Routine: PSE      (fn) *** ARCHIVED ***
C###  Routine: PSI      (fn) eval basis function at Xi coords
C###  Routine: PSI1     (fn) eval tensor prod. Lagrange & Hermite basis fns
C###  Routine: PSI2     (fn) eval Simplex basis functions
C###  Routine: PSI2_HERMITE (fn) eval Hermite Simplex basis functions
C###  Routine: PSI2_XI  (fn) eval simplex basis functions, Xi coords.
C###  Routine: PSI5     (fn) eval sector basis functions
C###  Routine: PSI8     (fn) eval auxiliary lagrange basis functions
C###  Routine: PSIM     (fn) eval contribution to Simplex basis function
C###  Routine: PXI      (fn) interpolates nodal array XE
C###  Routine: D2ZX     (fn) 2nd deriv of cart. coord wrt curvilinear
C###  Routine: DACOS_MOD(fn) returns corrected DACOS for out of range input
C###  Routine: DATAN_MOD  (fn) returns a corrected ATAN for zero inputs
C###  Routine: DET        (fn) determinant of 3*3 matrix
C###  Routine: DGREEN   normal der of Green's function
C###  Routine: DGREENA  normal der of axisymmetric Green's function
C###  Routine: DXZ      (fn) gradient of curvilinear coords wrt cartesian
C###  Routine: DZX      (fn) gradient of cart. coords wrt curvilinear
C###  Routine: DZXX     (fn) 2nd deriv of cart. coords wrt curvilinear
C###  Routine: GET_TV_VALUE_AT_TIME returns a T.V.'s value at a given time
C###  Routine: GREEN    Green's fns real BE problems (cart).
C###  Routine: GREENA   Green's fns for axisymmetric BE probs.
C###  Routine: HYPGREEN Green's function for hypersingular eqtn
C###  Routine: RFROMC   (fn) REAL*8 variable from character string
C###  Routine: SCALAR   (fn) returns scalar product of two vectors
C###  Routine: TIMER    (fn) *** ARCHIVED ***
C###  Routine: VTIME    (fn) returns the current time in seconds
C###  Routine: NORM_RAND_NUM (fn) normally distributed random number
C###  Routine: CM_RANDOM_NUMBER (fn) returns random number
C###  Routine: RANDSEEDS (fn) produces a random number from 3 seeds
C###  Routine: RAN0    (fn) returns a random num. with uniform distrib.
C###  Routine: RAN1    (fn) returns a random num. with uniform distrib.
C###  Routine: RAN10   (fn) returns a random num. of gaussian distrib.
C###  Routine: ERRORC  (fn) complementary error function
C###  Routine: DATA_DIST     (fn) compute dist between 2 data points
C###  Routine: NODE_DIST     (fn) compute dist between 2 node points
C###  Routine: SPLINE_INTERP (fn) interpolates a 1D spline
C###  Routine: CENT_DIFF2 (fn) 2nd order central diff at a node
C###  Routine: QUICK (fn) advective scheme
C###  Routine: ULTRAQUICK (fn) advection scheme
C###  Routine: QUICKEST (fn) advection scheme
C###  Routine: ULTRAQUICKEST (fn) advection scheme
C###  Routine: DMED (fn) double precision median
C###  Routine: DTSP (fn) *REMOVED* double precision triple scalar product
C###  Routine: ANALY_PHI_M (fn) Analytic transmembrane potential


C#### Module: FE10_LINUX
C###  Description:
C###    System specific routines (LINUX)

C###  Routine: CLOSEF        closes file
CC###  Routine: FIND_FILE     finds files in current directory
C###  Routine: OPENF         opens file
CC###  Routine: POST_FILE     posts file to printer
CC###  Routine: PURGE_FILE    purges file versions
C###  Routine: SETUPCMISS   setup initial CMISS constants

C#### Module: FE10_LINUX
C###  Description:
C###    System specific routines (LINUX)

C###  Routine: CLOSEF        closes file
CC###  Routine: FIND_FILE     finds files in current directory
C###  Routine: OPENF         opens file
CC###  Routine: POST_FILE     posts file to printer
CC###  Routine: PURGE_FILE    purges file versions
C###  Routine: SETUPCMISS   setup initial CMISS constants

C#### Module: FE10_SGI
C###  Description:
C###    System specific routines (SGI)

C###  Routine: CLOSEF        closes file
CC###  Routine: FIND_FILE     finds files in current directory
C###  Routine: OPENF         opens file
CC###  Routine: POST_FILE     posts file to printer
CC###  Routine: PURGE_FILE    purges file versions
C###  Routine: SETUPCMISS   setup initial CMISS constants

C#### Module: FE11
C###  Description:
C###    Routines for input of model parameters.

C###  Routine: IPACTI   input active muscle contraction properties
C###  Routine: IPACTV   input analytic activation time parameters
C###  Routine: IPAERO   input aerofoil parameters
C###  Routine: IPANAL   input analytic solution parameters
C###  Routine: IPBASE   input basis functions
C###  Routine: IPBOUN   input block boundary conditions
C###  Routine: IPCELL_READ   input cell parameters
C###  Routine: IPCELL_WRITE  output cell parameters
C###  Routine: IPCELL   input cell parameters
C###  Routine: IPCELL_CELLML   input cell parameters for a cellML model
C###  Routine: IPCELL_PROMPT   input cell parameters for simple models
C###  Routine: IPCOOR   input coordinate system
C###  Routine: IPCORN   input corner nodes (Bdry elements)
C###  Routine: IPCOUP   input coupling requirements between regions
C###  Routine: IPCUST   input mesh customisation parameters
C###  Routine: IPEIGE   input eigenvalue analysis parameters
C###  Routine: IPELEM   input element topology
C###  Routine: IPEQUA   input equation types
C###  Routine: IPEXPO   input export options
C###  Routine: IPFIBR   input fibre direction field
C###  Routine: IPFIEL   input additional field variables
C###  Routine: IPFIT    input opt.n params for geometry or field fit
C###  Routine: IPGAUS   input Gauss point array YG
C###  Routine: IPGRID   input grid parameters
C###  Routine: IPGRLV   input adaptive grid levels for residual calcs
C###  Routine: IPGROW   input growth parameters
C###  Routine: IPIMPO   input import options
C###  Routine: IPINIT   input initial and boundary conditions
C###  Routine: IPINIT_ELAS input init and bound conds for elasticity problems
C###  Routine: IPINVE   input inverse transfer matrix parameters
C###  Routine: IPITER   input iteration parameters
C###  Routine: IPLEAD   input electrocardiographic leads.
C###  Routine: IPLINE   input line scaling factors
C###  Routine: IPMATE   input material parameters
C###  Routine: IPMATE_COUP   input material parameters for coupling
C###  Routine: IPMOTI   input motion parameters
C###  Routine: IPMOTI_LUNG   input motion parameters for pulmonary
C###  Routine: IPNODE   input global node coordinates
C###  Routine: IPNOIS   input noise parameters
C###  Routine: IPNONL   input nonlinear analysis parameters
C###  Routine: IPNORM   input normal reversal parameters
C###  Routine: IPOPTI   input optimisation parameters
C###  Routine: IPPARA   input array dimension parameters
C###  Routine: IPPLIN   input polyline parameters
C###  Routine: IPREFE   input reference node/electrode location
C###  Routine: IPREFI   input refinement parameters
C###  Routine: IPREGI   input total number of regions
C###  Routine: IPQUAS   input quasi-static analysis parameters
C###  Routine: IPSHEE   *** ARCHIVED ***
C###  Routine: IPSING   input singularity pos.n (Bdry elements)
C###  Routine: IPSOLT   input time integration parameters
C###  Routine: IPSOLU   input solution parameters
C###  Routine: IPSOLV   input solution parameters
C###  Routine: IPSOUR   input source parameters
C###  Routine: IPTIME   input time variables
C###  Routine: IPTRAN   *REMOVED* input transformation parameters for Phigs
C###  Routine: IPTRSF   input transfer matrix parameters
C###  Routine: IPWIND   input window characteristics



C*** Integer

C#### Variable: CCR(MAX_COUP)
C###  Type: REAL
C###  Set_up: IPMATE_COUP
C###  Description:
C###    CCR(MAX_COUP) records real constants for coupling properties.
C###    NCCR is the total number of constants in CCR.
C###  See-Also: QUERY_COUP, UPCG_COUP coup00.cmn

C#### Variable: CELL_ICQS_VALUE(nmqi,nqv)
C###  Type: INTEGER
C###  Set_up: IPCELL
C###  Description:
C###    CELL_ICQS_VALUE(nmqi,nqv) stores the values of integer
C###    parameters for cellular models for each variant (nqv), as read
C###    from the ipcell file.
C###  See-Also: IPMAT3_CELL

C#### Variable: CELL_ICQS_SPATIAL(nmqi,nqv)
C###  Type: INTEGER
C###  Set_up: IPCELL
C###  Description:
C###    CELL_ICQS_SPATIAL(nmqi,nqv) stores a switch for each integer
C###    parameter, for each variant, which is used to determine whether
C###    the parameter is spatially varying or not. If the variable is
C###    spatially varying then its entry in this array will be
C###    initialised to -1 in IPCELL and then set to how the spatial
C###    variation is defined in IPMAT3_CELL1.
C###  See-Also: IPMAT3_CELL

C#### Variable: CELL_RCQS_VALUE(nmqr,nqv)
C###  Type: REAL*8
C###  Set_up: IPCELL
C###  Description:
C###    CELL_RCQS_VALUE(nmqr,nqv) stores the values of real
C###    parameters for cellular models for each variant (nqv), as read
C###    from the ipcell file.
C###  See-Also: IPMAT3_CELL

C#### Variable: CELL_RCQS_SPATIAL(nmqr,nqv)
C###  Type: INTEGER
C###  Set_up: IPCELL
C###  Description:
C###    CELL_RCQS_SPATIAL(nmqr,nqv) stores a switch for each real
C###    parameter, for each variant, which is used to determine whether
C###    the parameter is spatially varying or not. If the variable is
C###    spatially varying then its entry in this array will be
C###    initialised to -1 in IPCELL and then set to how the spatial
C###    variation is defined in IPMAT3_CELL1.
C###  See-Also: IPMAT3_CELL

C#### Variable: CELL_YQS_VALUE(niqs,nqv)
C###  Type: REAL*8
C###  Set_up: IPCELL
C###  Description:
C###    CELL_YQS_VALUE(niqs,nqv) stores the values of state and derived
C###    variables for cellular models for each variant (nqv), as read
C###    from the ipcell file.
C###  See-Also: IPMAT3_CELL

C#### Variable: CELL_YQS_SPATIAL(niqs,nqv)
C###  Type: INTEGER
C###  Set_up: IPCELL
C###  Description:
C###    CELL_YQS_SPATIAL(niqs,nqv) stores a switch for each state and
C###    derived variable, for each variant, which is used to determine
C###    whether the variable is spatially varying or not. If the
C###    variable is spatially varying then its entry in this array will
C###    be initialised to -1 in IPCELL and then set to how the spatial
C###    variation is defined in IPMAT3_CELL1.
C###  See-Also: IPMAT3_CELL

C#### Variable: COUP_CL(0:MAX_COUP,1:3)
C###  Type: INTEGER
C###  Set_up: IPCOUP
C###  Description:
C###    COUP_CL(no_coup,1:3) contains the coupling details for source terms
C###    in an advection-diffusion model of solutes. The coupling is
C###    through cellular processes like channels, pump and reactive
C###    pathways.
C###    COUP_CL(0,1) is the total number of couplings
C###    COUP_CL(no_coup,1) the class nx, of the coupling
C###    COUP_CL(no_coup,2) the type of coupling (channel, pump,
C###      reactive pathway)
C###    COUP_CL(no_coup,3) the sub-type of coupling eg. soduim channel,
C###      potassium channel, IKNA pump, etc..., each are indexed by a
C###      unique number.
C###  See-Also: IPMATE_COUP, coup00.cmn

C#### Variable: COUPLINKS(0:MAX_COUP,1:9)
C###  Type: INTEGER
C###  Set_up: IPMATE_COUP
C###  Description:
C###    COUPLINKS(no_coup,1:9) contains details on the information required
C###    by the coupling models specified in COUP_CL. The coupling is
C###    through cellular processes like channels, pump and reactive
C###    pathways. NCA is the total number of coupling details in COUPLINKS.
C###    COUPLINKS(no_coup,1) the class nx, of the coupling
C###    COUPLINKS(no_coup,2) question number for domain
C###    COUPLINKS(no_coup,3) question number for tissue properties
C###    COUPLINKS(no_coup,4) question number for solutes
C###    COUPLINKS(no_coup,5) question number for protein model properties
C###    COUPLINKS(no_coup,6) question number for other question
C###    COUPLINKS(no_coup,7) question number for directions
C###    COUPLINKS(no_coup,8) type of answer (constant, field in XP,
C###       dependent variable in YP)
C###    COUPLINKS(no_coup,8) value for COUPLINKS(no_coup,8). For a constant,
C###      the position of the constant in CRR is recorded. For a field the
C###      number of the field in XP, for a dependent variable the number of
C###      the dependent variable.
C###  See-Also: QUERY_COUP, UPCG_COUP coup00.cmn

C#### Variable: INV_APPROACH
C###  Type: INTEGER
C###  Set_up: IPINVE
C###  Description:
C###    Which approach taken to solve an inverse problem. ie.
C###    potential imaging or activation imaging.

C#### Variable: NPLIST3(np)
C###  Type: INTEGER
C###  Set_up: IPTRSF,APTRSF
C###  Description:
C###    NPLIST3(np) is the node list of the first surface nodes used
C###    in the construction of T_BH

C#### Variable: NPLIST4(np)
C###  Type: INTEGER
C###  Set_up: IPTRSF,APTRSF
C###  Description:
C###    NPLIST4(np) is the node list of the second surface nodes used
C###    in the construction of T_BH

C#### Variable: NPLIST5(np)
C###  Type: INTEGER
C###  Set_up: IPTRSF,APTRSF
C###  Description:
C###    NPLIST5(np) is the node list of the outer surface nodes used
C###    in the construction of T_BH

C#### Variable: CPLST(0:NCONMX,2)
C###  Type: INTEGER
C###  Set_up: IPCOUP
C###  Description:
C###    CPLST(0,1) contains the number of coupled grid points.
C###    CPLST(nq,1) contains a grid point number to be coupled.
C###    CPLST(nq,2) contains the number of the nearest grid point in the
C###    region to which CPLST(nq,1) is being coupled.

C#### Variable: ACTN_MIN(i)
C###  Type: REAL*8
C###  Set_up: IPOPTI
C###  Description:
C###    Change limits in optimisation for activation problems. The
C###    i indicies are 1|2|3|4 for min activation time change|global
C###    activation limit|transmembrane jump change|window width change.

C#### Variable: ACTN_MAX(i)
C###  Type: REAL*8
C###  Set_up: IPOPTI
C###  Description:
C###    Change limits in optimisation for activation problems. The
C###    i indicies are 1|2|3|4 for max activation time change|global
C###    activation limit|transmembrane jump change|window width change.

C#### Variable: NPMIN(nj)
C###  Type: REAL*8
C###  Set_up: IPOPTI
C###  Description:
C###    Change limits in optimisation

C#### Variable: NPMAX(nj)
C###  Type: REAL*8
C###  Set_up: IPOPTI
C###  Description:
C###    Change limits in optimisation

C#### Variable: CM_FREQUENCY
C###  Type: INTEGER
C###  Set_up: IPEXPO
C###  Description:
C###    Mapping between integer time steps and real time for
C###    inverse activation.

C#### Variable: DIPOLE_CEN(i,0:t,n,nr,nx)
C###  Type: REAL*8
C###  Set_up: IPSOUR
C###  Description:
C###    DIPOLE_CEN(i,t,n,nr) is the ith component of the
C###    centre of dipole n or region nr for time node t (integer).
C###    If the dipole is not time dependent then its centre is stored
C###    in position t=0. Note: the times of for the time points of
C###    a moving dipole are stored in the i=4 position.
C###  See-Also: NDIPOLES, DIPOLE_CEN_NTIME, DIPOLE_DIR, DIPOLE_DIR_NTIME

C#### Variable: DIPOLE_CEN_NTIME(n,nr)
C###  Type: INTEGER
C###  Set_up: IPSOUR
C###  Description:
C###    DIPOLE_CEN_NTIME(n,nr)=0 if the centre of dipole n
C###    in region nr is not time dependent, otherwise
C###    DIPOLE_CEN_TIME(n,nr) stores the number of time positions in
C###    the path.
C###  See-Also: NDIPOLES, DIPOLE_CEN, DIPOLE_DIR, DIPOLE_DIR_NTIME

C###  Variable: DIPOLE_DIR(i,0:t,n,nr,nx)
C###  Type: REAL*8
C###  Set_up: IPSOUR
C###  Description:
C###    DIPOLE_DIR is the ith component of the dipole
C###    vector of dipole n of regoin nr for time node t (integer).
C###    If the dipole is not time dependent then its direction is stored
C###    in position t=0. Note: the times of for the time points of
C###    a moving dipole are stored in the i=4 position.
C###  See-Also: NDIPOLES, DIPOLE_DIR_NTIME, DIPOLE_CEN, DIPOLE_CEN_NTIME

C#### Variable: DIPOLE_DIR_NTIME(n,nr)
C###  Type: INTEGER
C###  Set_up: IPSOUR
C###  Description:
C###    DIPOLE_DIR_NTIME(n,nr)=1 if the direction of dipole n in region
C###    nr is not time dependent, otherwise DIPOLE_DIR_NTIME(n,nr)
C###    stores the number of time positions in the path.
C###  See-Also: NDIPOLES, DIPOLE_CEN, DIPOLE_CEN_NTIME, DIPOLE_DIR_NTIME

C#### Variable: IBT(3,ni,nbf)
C###  Type: INTEGER
C###  Set_up: IPBASE,BASIS1,BASIS2,BASIS5,BASIS6
C###  Description:
C###    <HTML> <PRE>
C###    IBT is the index for family basis type nbf. IBT(1,ni,nbf) is the
C###    type of interpolation in the ni direction, IBT(2,ni,nbf) is the
C###    order of the interpolation.
C###    Interpolation types  are:
C###      IBT(1,ni,nbf) = 1  Lagrange
C###                    = 2  Hermite
C###                    = 3  Simplex
C###                    = 4  Seredipity
C###                    = 5  Sector (collapsed at xi=0)
C###                    = 6  Sector (collapsed at xi=1)
C###                    = 7  Transition
C###                    = 8  Singular
C###                    = 9  Fourier
C###
C###    Indicies for IBT(2,ni,nbf) are:
C###      1..3 for linear, quadr. or cubic Lagrange for IBT(1,ni,nbf)=1
C###      1..3 for cubic, quadratic1 or quadratic2 Hermite    "       2
C###                                                          "       3
C###                                                          "       4
C###      1..4 for linear, quadratic, cubic Lagrange, Hermite "     5,6
C###                                                          "       7
C###                                                          "       8
C###      N for number of Fourier terms (including constant). "       9
C###
C###    IBT(3,ni,nbf) is used for sector bases (IBT(1,ni,nbf) = 5 or 6)
C###      to store the direction perpendicular to collapsed face.
C###      If IBT(3,ni,nbf)=0, the element is not collapsed in direction ni.
C###    IBT(3,1,nbf) is used simplex bases (IBT(1,ni,nbf) = 3)
C###      To store the parametric coordinate type. IBT(3,1,nbf)=1 for
C###      area coordinates IBT(3,1,nbf)=2 for Xi coordinates.
C###    </PRE> </HTML>

C#### Variable: IWRIT5
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    Controls the output levels from the optimisation.

C#### Variable: MAX_MAJOR_ITER
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    Controls the maximum number of major iterations for E04UPF
C###    optimiser. Currently only setup for activtion problems.
C###    The value is -1 if the default value is to be used.
C###  See-Also: MAX_MINOR_ITER

C#### Variable: MAX_MINOR_ITER
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    Controls the maximum number of minor iterations for E04UPF
C###    optimiser. Currently only setup for activtion problems.
C###    The value is -1 if the default value is to be used.
C###  See-Also: MAX_MAJOR_ITER

C#### Variable: NENQ
C###  Type: INTEGER
C###  Set_up: IPGRID,CALC_LATTICE_MAP
C###  Description:
C###    NENQ(0:8,NQM) gives the list of elements which contain the grid
C###    point nq. For the lattice grid scheme NENQ(0,nq) stores the
C###    element which contains the principal lattice point for nq.      

C#### Variable: NCT
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    NCT(nr,nx) stores the number of variable types in each region.

C#### Variable: NLATXIM
C###  Type: INTEGER
C###  Set_up: IPGRID
C###  Description:
C###    NLATXIM is the maximum number of lattice grid points
C###    in any element, in any xi direction. The value of
C###    NLATXIM is derived from the maximum value of NQXI
C###    for all grid schemes.      
           
C#### Variable: NLQ
C###  Type: INTEGER
C###  Set_up: IPGRLV
C###  Description:
C###    NLQ(nq) is either the adaptive grid pt level na for calculating
C###    a residual at nq (0<na<9) or na//na+1 when nq at level na is
C###    interpolated from level na+1.

C#### Variable: NQET
C###  Type: INTEGER
C###  Set_up: IPGRID
C###  Description:
C###    NQET(NQSCM) stores the number of local grid points for each
C###    grid point scheme.

C#### Variable: NQLIST
C###  Type: INTEGER
C###  Set_up: PARSE_GRID
C###  Description:
C###    NQLIST is a temporay list of grid points.

C#### Variable: NQNE
C###  Type: INTEGER
C###  Set_up: IPGRID
C###  Description:
C###    NQNE(NEQM,NQEM) gives the global grid point number for an
C###    element ne and a local grid point number nq.

C#### Variable: NQS
C###  Type: INTEGER
C###  Set_up: IPGRID
C###  Description:
C###    NQS(NEQM) gives the grid grid scheme applied in element ne.

C#### Variable: NQSCNB
C###  Type: INTEGER
C###  Set_up: IPGRID
C###  Description:
C###    NQSCNB(NQSCM) gives the finite element basis function number
C###    to generate global grid point positions for each grid scheme.

C#### Variable: NQXI
C###  Type: INTEGER
C###  Set_up: IPGRID
C###  Description:
C###    NQXI(0:NIM,NQSCM), the 0th entry for a given scheme is the
C###    number of xi coordinates defined for the scheme. The other
C###    entries are the number of grid points defined in each xi 
C###    direction. If the number of grid points in an xi direction
C###    is not defined it is set to 1.      

C#### Variable: N_SITES
C###  Type: INTEGER
C###  Set_up: IPACTV
C###  Description:
C###    N_SITES is the number of initial activation sites

C#### Variable: ACTVN_SITE(3,10)
C###  Type: REAL
C###  Set_up: IPACTV
C###  Description:
C###    ACTVN_SITE(nj,n_site) are the coordinates of the initial
C###    activation site n_site.

C#### Variable: ICALC_TRANSFER
C###  Type: INTEGER
C###  Set_up: IPINVE
C###  Description:
C###    <HTML>
C###    ICALC_TRANSFER specifies the type of transfer to use for
C###    an imaging technique.
C###    <P> If Potential Imaging: <UL>
C###      <LI>  1 if the transfer matrix is to be inverted explicitly
C###      <LI>  2 if it is to be constructed from the matrix equations
C###    </UL>
C###    <P> If Activation Imaging: <UL>
C###      <LI>  1 if zero-crossing technique
C###    </UL>
C###    <HTML>

C#### Variable: INDEX_PLIN
C###  Type: INTEGER
C###  Set_up: IPPLIN
C###  Description:
C###    INDEX_PLIN is the index for a particular polyine.

C#### Variable: INDEX_PLIN_TYPE(index_plin)
C###  Type: INTEGER
C###  Set_up: IPPLIN
C###  Description:
C###    INDEX_PLIN_TYPE(index_plin) indicates the type (piecewise
C###    linear/Bezier).

C#### Variable: INP(nn,ni,nbf)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    INP gives the index for element node nn in each Xi-direction.
C###    Thus: INP(nn,1..,nb) = 1,2,2 indicates that node nn is the
C###    first node in the Xi(1) direction and second in the Xi(2) and
C###    Xi(3) directions; INP(nn,1..,nb) = 1,1,2 for a cubic Serendipity
C###    element indicates thet node nn is the first node in the Xi(1)
C###    and Xi(2) directions and second in the Xi(3) direction.
C###    INP(nn,1..,nb) = 1,1,3 for a quadratic tetrahedron is the apex.

C#### Variable: ACTN_IREGULARISE
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    ACTN_IREGULARISE is the type of additional constraints applied
C###    for activation inverse optimisations.

C#### Variable: IREGULARISE
C###  Type: INTEGER
C###  Set_up: IPINVE
C###  Description:
C###    IREGULARISE is the type of regularisation to be used for
C###    potential inverses.

C#### Variable: ITYP9(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPNONL
C###  Description:
C###    ITYP9(nr,nx) is 1..4 for full Newton/modified Newton
C###    /BFGS inverse/Sequential quadratic programming for nonlinear
C###    problems,or 1,2,3 for direct solve/Gauss-Seidel iterations
C###    /Gauss-Seidel with multigrid accel.n for linear problems.

C#### Variable: ITYP10(nr)
C###  Type: INTEGER
C###  Set_up: IPCOOR
C###  Description:
C###    ITYP10(nr) is 1..5 for coordinates rectangular cartesian/
C###    cylindrical polar/spherical polar/prolate spheroidal/oblate
C###    spheroidal.

C#### Variable: ITYP11(nr)
C###  Type: INTEGER
C###  Set_up: IPCOOR
C###  Description:
C###    ITYP11(nr) is 1..5 for dependent variable coordinate system.

C#### Variable: ITYP15(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    <HTML>
C###    Non-zero ITYP15(nr,nx) indicates a non-galerkin FEM.
C###    Possible values are:
C###    <LI>0 - standard Galerkin,
C###    <LI>1 - Petrov-Galerkin weighting functions based on material
C###    derivatives,
C###    <LI>1 - Petrov-Galerkin weighting functions based on local
C###    element xi derivatives,
C###    <LI>3 - Galerkin with derivative discontinuity stabilizing terms.
C###    </HTML>

C#### Variable: ITYP16(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP16(nr,nx) is 1..4 for fully explicit/fully implicit/
C###    Crank Nicholson/Lax Wendroff finite difference solution
C###    methods.

C#### Variable: ITYP19(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP19(nr,nx) stores the type of cellular model to be used,
C###    electrical, mechanical, metabolism, coupled, etc...

C#### Variable: ITYP20(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    ITYP20(nr,nx) indicates the contact method employed.
C###    Currently the possible options are:
C###    1. Frictionless Penalty Method Contact.

C#### Variable: ITYP21(nr)
C###  Type: INTEGER
C###  Set_up: IPMAPPING
C###  Description:
C###    ITYP21(nr,nx) indicates method of equivalent dof mapping
C###    1. Search for equivalent nodal positions and derivatives
C###    2. Search for equivalent nodal positions and use explicit derivative mapping
C###    3. Use explicit nodal position and derivative mapping

C#### Variable: IWRIT1(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    IWRIT1(nr,nx) controls output printing frequency for time
C###    stepping / nonlinear solutions.

C#### Variable: IWRIT2(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPNONL
C###  Description:
C###    IWRIT2(nr,nx) is 1 / 2 for equilibrium solution information only /
C###    intermediate solution information also.

C#### Variable: IWRIT3(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPNONL
C###  Description:
C###    IWRIT3(nr,nx) is 0 / 1 / 2 for solver progress only / solution
C###    vectors also / residual vectors also.

C#### Variable: IWRIT4(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    IWRIT4(nr,nx) controls output printing for the linear solution.
C###    The values are 0/1/2/3/4/5 for No output/Timing information/
C###    Solution information/Global solution matrices/Global
C###    stiffness matrices/Element matrices and information.

C#### Variable: JTYP2A
C###  Type: INTEGER
C###  Set_up: IPCOOR
C###  Description:
C###    JTYP2A  is 1 for non-standard degrees of freedom mapping
C###    (coincident versions of nodes) and 0 otherwise.

C#### Variable: JTYP2B
C###  Type: INTEGER
C###  Set_up: IPCOOR
C###  Description:
C###    JTYP2B  is 1 for all non-standard line mappings and 0 otherwise.

C#### Variable: JTYP2C
C###  Type: INTEGER
C###  Set_up: IPCOOR
C###  Description:
C###    JTYP2C  is 1 for non-standard degrees of freedom mapping
C###    (hanging nodes) and 0 otherwise.

C#### Variable: JTYP4
C###  Type: INTEGER
C###  Set_up: IPCOOR
C###  Description:
C###    JTYP4  is 1..3 for geometry unsymmetric/cylindrical symmetric/
C###    spherical symmetric.  (Note: for cylindrical symmetric, if
C###    ityp10=1 then radius is nj=2 coordinate; if 2 or 3 then radius
C###    is nj=1 coordinate; if 4 then r=focus*sinhX1*sinX2).

C#### Variable: JTYP6
C###  Type: INTEGER
C###  Set_up: IPCOOR
C###  Description:
C###    JTYP6  is 1,2  for global coordinate system constant/specified
C###    by elements.

C#### Variable: JTYP9
C###  Type: INTEGER
C###  Set_up: IPFIBR
C###  Description:
C###
C###    This variable has been replaced by NJ_LOC(NJL_FIBR,0,nr).
C###
C###    JTYP9 = 1 for fibre directions (lying in Xi(1),Xi(2) plane).
C###    Note: The fibre coord Nu(1) lies in the Xi(1)-Xi(2) plane, Nu(2)
C###    is orthogonal to the fibre coord and lies in the Xi(1)-Xi(2)
C###    plane and the remaining Nu coord is orthogonal to this plane.
C###    The fibre angle eta is eta(1), the angle between the fibre
C###    coordinate and the Xi(1) coordinate.  The Nu coordinates are
C###    stress coordinates and are orthonormal with metric a(i,j)=
C###    Kronecker delta. JTYP9 = 3 for sheet directions (in plane
C###    normal to Xi(1),Xi(2)).

C#### Variable: JTYP10
C###  Type: INTEGER
C###  Set_up: IPCOOR
C###  Description:
C###    JTYP10 is 1,2,3 for type of 'radial' interpolation
C###    (ITYP10(1)>1 only).

C#### Variable: JTYP11
C###  Type: INTEGER
C###  Set_up: IPFIEL
C###  Description:
C###
C###    This variable has been replaced by NJ_LOC(NJL_FIEL,0,nr).
C###
C###    JTYP11 is the number of additional field variables.

C#### Variable: JTYP12
C###  Type: INTEGER
C###  Set_up: IPFIBR
C###  Description:
C###    JTYP12 is 1,2 for fibres defined wrt Xi(1) or Xi(2) coordinates.

C#### Variable: JTYP13
C###  Type: INTEGER
C###  Set_up: IPFIBR
C###  Description:
C###    JTYP13 is 1,2 for fibre angles defined in degrees or radians.


C#### Variable: KTYP4
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    KTYP4 determines whether or not the solution matrices are
C###    written to a file. The possible values for KTYP4 are 0/1/2 for
C###    for No output/Solution matrix/+ RHS vector in ascii form.  If
C###    the matrices are to be written in binary form then the codes
C###    are the same except negative.

C#### Variable: KTYP7
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP7 is equation parameters constant wrt time/defined in
C###    subroutine USER/read from file IPC at each time step.

C#### Variable: KTYP8
C###  Type: INTEGER
C###  Set_up: DEFIT
C###  Description:
C###    <HTML>
C###    KTYP8 varies from 1..9.
C###    <UL>
C###      <LI> (1) Geometric fitting Note: ITYP6(nr) is 1,2 for
C###                 linear/non-linear
C###      <LI> (2) fibre-sheet fitting problem
C###      <LI> (3) field fitting
C###      <LI> (4) signal fitting
C###      <LI> (5) motion fitting with Fourier basis
C###      <LI> (6) fitting with optimisation
C###      <LI> (7) material fitting
C###      <LI> (8) patch fitting
C###      <LI> (9) spline fitting for signal
C###    </HTML>
C###    </UL>

C#### Variable: KTYP10
C###  Type: INTEGER
C###  Set_up: IPNONL
C###  Description:
C###    KTYP10 is 1,2 for solution with no search/linear search.

C#### Variable: KTYP1A
C###  Type: INTEGER
C###  Set_up: IPNONL
C###  Description:
C###    KTYP1A is 1,2 for series/parallel element stiffness matrix
C###    calculations.

C#### Variable: KTYP1D
C###  Type: INTEGER
C###  Set_up: IPNONL
C###  Description:
C###    <HTML>
C###    KTYP1D determines the method for calculation of residual
C###    derivatives for nonlinear problems:
C###    <PRE>
C###    KTYP1D=1 - algebraically.
C###    KTYP1D=2 - one-sided finite differences.
C###    KTYP1D=3 - central finite differences.
C###    </PRE>
C###    </HTML>

C#### Variable: KTYP14
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP14 is 0,>0 if material parameter is incremented.

C#### Variable: KTYP15
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP15 is equation type parameter.

C#### Variable: KTYP16
C###  Type: INTEGER
C###  Set_up: IPEIGE
C###  Description:
C###    KTYP16 is 1/2 for lowest/highest eigenvalues.

C#### Variable: KTYP17
C###  Type: INTEGER
C###  Set_up: IPEIGE
C###  Description:
C###    KTYP17 is number of eigenpairs required.

C#### Variable: KTYP20
C###  Type: INTEGER
C###  Set_up: IPITER
C###  Description:
C###    KTYP20 is type of iteration procedure (1:solution iteration).

C#### Variable: KTYP21
C###  Type: INTEGER
C###  Set_up: IPITER
C###  Description:
C###    KTYP21 is the specific iteration prodedure (1:forward problem
C###    calculation).

C#### Variable: KTYP22
C###  Type: INTEGER
C###  Set_up: IPSOLT
C###  Description:
C###    KTYP22 is 1..3 for linear/quadratic/cubic time-integration
C###    algorithm.

C#### Variable: KTYP23
C###  Type: INTEGER
C###  Set_up: IPSOLT
C###  Description:
C###    KTYP23 is 1..3 for fixed time step/automatic stepping/read
C###    from file.

C#### Variable: KTYP24
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP24 is the sparsity of the global matrices. If KTYP24 is 0
C###    the global matrices are not stored sparsly. If KTYP24 is 1 they
C###    are stored in a compressed-row sparsity format. If KTYP24 is
C###    2 they are stored in a row-column sparsity format. If KTYP24
C###    is 3 they are stored in a compressed-column sparsity format.
C###  See-Also: SPARSITY STRUCTURES

C#### Variable: KTYP25
C###  Type: INTEGER
C###  Set_up: IPSOLT
C###  Description:
C###    KTYP25 is type of driving function in Fourier analysis.

C#### Variable: KTYP1B
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    <HTML>
C###    KTYP1B determines how derivatives are optimised.  Possible
C###    values are:
C###    <LI>1:  Optimise as components with non-linear constraint on
C###        derivative magnitude.
C###    <LI>2:  Optimise in angular coordinates - no constraints.
C###    <LI>3:  Optimise as components, -1 <= x <= 1 for each component.
C###    </HTML>

C#### Variable: KTYP26
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    KTYP26 is 1,2, for optimising material parameters/geometric
C###    parameters/microstructure parameters.
C!!! CS 29/10/2000 Note true any more
CC###    Note: KTYP26 = 3 KTYP27 = 1 is used for stripe
CC###    intersection calculation.

C#### Variable: KTYP27
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    KTYP27 is 1..13 for type of minimisation objective function.
C!!! CS 29/10/2000 Note true any more
CC###    Note: KTYP26=3 KTYP27=1 is used for stripe intersection
CC###    calculation.


C#### Variable: KTYP27B
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    <HTML>
C###    KTYP27B type of field used in the minimisation objective function.
C###    Currently only set for KTYP27=13 (dipole inverses)
C###    <UL>
C###      <LI> 1 = Magnetic field
C###      <LI> 2 = Potential field
C###      <LI> 3 = Both magnetic and potential field
C###    </UL>
C###    </HTML>


C#### Variable: KTYP28
C###  Type: INTEGER
C###  Set_up: IPOPTI,IPFIT
C###  Description:
C###    KTYP28 is number of sets of measurements in fit.

C#### Variable: KTYP29
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    KTYP29 is the optimisation package used: 1=NPSOL, 2=MINOS.

C#### Variable: KTYP29B
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    KTYP29B determines the types of residuals for NPSOL.  At present
C###    this is only used for data fitting: 1 = components of data
C###    projections, 2 = Euclidean norm of data projection magnitudes, 3
C###    = Squares of data projection magnitudes.

C#### Variable: KTYP29C
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    KTYP29C determines the type/locations of residuals
C###    used in micro-structure optimisation problems.

C#### Variable: KTYP31
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP31 is 1,2 for activation model implemented forwards/
C###    backwards.

C#### Variable: KTYP32
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP32 is 1,2 for monodomain/bidomain.

C#### Variable: KTYP33
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP33 is 1,2 for cubic/quintic/heptic poly if ITYP3(nr,nx)=1;
C###    1,2 for standard/Rogers FHN if ITYP3(nr,nx)=2; or 1,2 for
C###    standard/Calif. VCD if ITYP3(nr,nx)=3; or 1,2,3 for BR/EJ/DR
C###    sodium kinetics if ITYP3(nr,nx)=4; or 1,2 for standard/Princeton
C###    JRW models ITYP3(nr,nx)=5.

C#### Variable: KTYP34
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP34 is 1,2 for normal/ischemic if ITYP3(nr,nx)=3 and
C###    KTYP33=2.

C#### Variable: KTYP36
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    KTYP36 is 1 if using Dynamic Tracking of Active Region (DTAR)
C###    for cardiac activation problems.

C#### Variable: KTYP37
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    KTYP37 is 1 if using Euler integration of the ionic currents;
C###    2 if using Improved Euler; 3 if using 4th order Runge-Kutta;
C###    4 if using Adams-Moulton (2nd Order) with variable time
C###    stepping (obselete); 5 if using Adams-Moulton (variable Order)
C###    with variable time stepping; 6 if using LSODA.

C#### Variable: KTYP3B
C###  Type: INTEGER
C###  Set_up: DEGRID
C###  Description:
C###    KTYP3B is 1 for regular grid schemes, 2 for a Gauss point
C###    grid scheme

C#### Variable: KTYP43
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP43 is 0..3 for thermal strain options.

C#### Variable: KTYP58(nr)
C###  Type: INTEGER
C###  Set_up: IPMOTI,IPMAT5
C###  Description:
C###    KTYP58(nr) is 1,2 for conventional/isochoric element, or
C###    1,2 for geometric coordinates/displacements of dependent
C###    variable, or 1..4 for motion type as Geometric coordinates/
C###    displacements/flow/lung gas flow.

C#### Variable: KTYP59(nr)
C###  Type: INTEGER
C###  Set_up: IPACTI,IPMAT5
C###  Description:
C###    KTYP59(nr) is for elastance/Hill-type/fading-memory
C###    formulation.

C#### Variable: KTYP5G(nr)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP5G(nr) indicates if the body in region nr uses contact mechanics.
C###    0. No Contact Mechanics (Default).
C###    1. Contact Mechanics used.

C#### Variable: KTYP5H(nr)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP5H(nr) indicates if the body in region nr is a host mesh.
C###    0. No (Default).
C###    1. Host mesh.

C#### Variable: KTYP5I(nr)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP5H(nr) indicates if the body in region nr is inertia dependent.
C###    0. No Inertia (Default).
C###    1. Inertia.

C#### Variable: KTYP5J(nr)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP5J(nr) indicates if the contact problem starts with active constraints.
C###    0. Initially inactive (Default).
C###    1. Initially active.

C#### Variable: KTYP5K(nr)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP5K(nr) indicates if the contact problem involves a psuedo-viscous term.
C###    0. No psuedo-viscosity (Default).
C###    1. Psuedo-viscosity.

C#### Variable: KTYP61
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    KTYP61 enumerates the scheme for the advection terms

C#### Variable: KTYP62
C###  Type: INTEGER(2)
C###  Set_up: IPSOLV
C###  Description:
C###    KTYP62(1), KTYP62(2) is how often, how much
C###    the Voronoi mesh is to be smoothed

C#### Variable: KTYP71
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    KTYP71 is 1 if pressure loads read from file PRESS.VSAERO
C###    (ID=14) after flow solution by VSAERO.

C#### Variable: KTYP90
C###  Type: INTEGER
C###  Set_up: IPCOUP
C###  Description:
C###    KTYP90 is 1/2/3 for saturated-unsaturated/heart-body/? coupling.

C#### Variable: KTYP91
C###  Type: INTEGER
C###  Description:
C###    KTYP91 is the extent of derivative equation generation for
C###    interface nodes (choice only available with hermite
C###    interpolation).

C#### Variable: KTYP92
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    KTYP92 indicates what combination of derivative equations are
C###    to be used when solving BEM problems using hypersingular
C###    equations. It can be 1/2/3 for BIE + TBIE/BIE + NBIE/NBIE + TBIE
C###    in 2D or BIE + s1, s2 and s1s2 TBIE/BIE + s1, s2 TBIE and NBIE/
C###    s1, s2, s1s2 TBIE + NBIE in 3d. (HYP is .TRUE. in this case
C###    and Hermite interpolation has been used for the dependent
C###    variable).

C#### Variable: KTYP94
C###  Type: INTEGER
C###  Set_up: IPTRSF
C###  Description:
C###    KTYP94 is 1/2 for direct/algebraic evaluation of the transfer
C###    matrix T_BH

CC AJPs 191297
C#### Variable: KTYP95
C###  Type: INTEGER
C###  Set_up: IPTRSF
C###  Description:
C###    KTYP95 is 1 for single layer (epi potls to body surface)
C###    KTYP95 is 2 for doulbe layer (trnasmembrane to body surface)
C###    KTYP95 is 3 for double layer (transmembrane to epi)
CC AJPe

C#### Variable: LEADELECS(0:nleadelec,nlead)
C###  Type: INTEGER
C###  Set_up: IPLEAD
C###  Description:
C###    LEADELECS(0:nleadelec,nlead) are the electrodes that are
C###    involved in lead nlead. LEADELECS(0,nlead) is the number of
C###    electrodes involed in lead nlead and LEADELECS(1..,nlead) are
C###    the electrode numbers.

C#### Variable: MAX_ITERATIVE_ITERS(nx)
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    MAX_ITERATIVE_ITERS(nx) is the max #iterations to use
C###    for the iterative solvers for problem class nx.

C#### Variable: MOTION_TYPE
C###  Type: INTEGER
C###  Set_up: IPMOTI
C###  Description:
C###    MOTION_TYPE is 1,2 for Fourier coefficients/Spreadsheet column.

C#### Variable: NBASEF(nbf,0:nbc)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NBASEF gives information about the children for family basis
C###    nbf. NBASEF(nbf,0) is the number of children in family basis
C###    nbf. NBASEF(nbf,1) is the global parent basis number of family
C###    basis nbf. NBASEF(nbf,2..) are global basis numbers of the
C###    child nbc in family basis nbf.

C#### Variable: NBC(nbf)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NBC is basis function type choice for family basis nbf.

C#### Variable: NBH(nh,nc,ne)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NBH identifies the family basis number nbf for dependent
C###    variable nh and derivative nc.

C#### Variable: NBH(nh,nc,ne)
C###  Type: INTEGER
C###  Set_up: IPBAS3,IPBAS4,IPBAS5,IPBAS9
C###  Description:
C###    NBH identifies the family basis number nbf for dependent
C###    variable nh and derivative nc.

C#### Variable: NBHF(nh,nc,nf)
C###  Type: INTEGER
C###  Set_up: CALC_FACE_BASIS_DEP
C###  Description:
C###    NBHF(nh,nc,nf) identifies the family basis number nbf for
C###    dependent variable nh and derivative nc for the global face
C###    nf.

C#### Variable: NBJ(nj,ne)
C###  Type: INTEGER
C###  Set_up: IPELEM
C###  Description:
C###    NBJ is the basis function type for geometric variable nj in
C###    element ne.

C#### Variable: NB_MOTION
C###  Type: INTEGER
C###  Set_up: IPMOTI
C###  Description:
C###    NB_MOTION is Fourier basis number.

C#### Variable: NCA
C###  Type: INTEGER
C###  Set_up: IPMATE_COUP
C###  Description:
C###    NA is the total number of constants in COUPLINK.
C###  See-Also: QUERY_COUP, UPCG_COUP coup00.cmn

C#### Variable: NCCR
C###  Type: INTEGER
C###  Set_up: IPMATE_COUP
C###  Description:
C###    NCCR is the total number of constants in CCR.
C###  See-Also: QUERY_COUP, UPCG_COUP coup00.cmn

C#### Variable: NCONMX
C###  Type: INTEGER
C###  Set_up: IPCOUP
C###  Description:
C###    NCONMX is a local parameter depicting the maximum number
C###    of grid points which may be coupled between regions.

C#### Variable: NDET(nbf,0:nn)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NDET is the number of BE subelements the (psi1,psi2)
C###    element is split into.

C#### Variable: NDIPOLES(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPSOUR
C###  Description:
C###    NDIPOLES(nr) is the number of dipole sources in region nr.
C###  See-Also: DIPOLE_CEN, DIPOLE_DIR

C#### Variable: NEELEM(0:ne,0:nr)
C###  Type: INTEGER
C###  Set_up: IPELEM
C###  Description:
C###    NEELEM(0,nr) is the total number of elements in region nr.
C###    NEELEM(noelem,nr), noelem=1..NEELEM(0,nr) are the element
C###    numbers.  NEELEM(0,0) is the total number of elements defined
C###    in all regions.
C###  See-Also: NEELEM,NFFACE,NLLINE,NPNODE


C#### Variable: NET(0:nr)
C###  Type: INTEGER
C###  Set_up: IPELEM
C###  Description:
C###    NET(nr) is the highest element number in region nr.
C###    NET(0)  is the highest element number in any region.

C#### Variable: NE_AERO(no_aero,i)
C###  Type: INTEGER
C###  Set_up: IPAERO
C###  Description:
C###    NE_AERO(0,1) is the number of elements adjacent to upper
C###    surface of aerofoil, NE_AERO(0,2) is the number of elements
C###    adjacent to lower surface of aerofoil, NE_AERO(no_aero,1) are
C###    element numbers adjacent to upper surface, NE_AERO(no_aero,2)
C###    are element numbers adjacent to lower surface.

C#### Variable: NE_WAKE(no_wake,i)
C###  Type: INTEGER
C###  Set_up: IPAERO
C###  Description:
C###    NE_WAKE(0,1) is the number of elements adjacent to upper
C###    surface of wake, NE_WAKE(0,2) is the number of elements
C###    adjacent to lower surface of wake, NE_WAKE(no_wake,1) are
C###    element numbers adjacent to upper surface, NE_WAKE(no_wake,2)
C###    are element numbers adjacent to lower surface.

C#### Variable: NFBASE(2,nb)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NFBASE gives information about the basis family that
C###    global basis number nb belongs to. NFBASE(1,nb) is the family
C###    basis number that global basis nb belongs to. NFBASE(2,nb) is
C###    local child number that corresponds to global basis nb.

C#### Variable: NGAP(ni,nb)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NGAP is number of Gauss points in Xi-coordinate direction ni
C###    for global basis nb.

C#### Variable: NG_FIT(nhj,njj)
C###  Type: INTEGER
C###  Set_up: IPFIT
C###  Description:
C###    NG_FIT(nhj,njj) is the Gauss variables to be fitted in Gauss
C###    point fitting for variable nhj of fit variable njj.

C#### Variable: NHE(ne,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    NHE(ne,nx) is the number of dependent variables defined in
C###    element ne for problem type nx.

C#### Variable: NHQ(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    NHQ(nr,nx) is the number of dependent variables defined for
C###    set of grid points in region nr delonging to equation class
C###    nx. Note: this array has no nq (grid point) index as currently
C###    there are no problems solved in cmiss which have different
C###    numbers of dependant variables associated with each grid point
C###    for a particular problem.

CC#### Variable: NH_FIT(2,nhx)
CC###  Type: INTEGER
CC###  Set_up: IPFIT
CC###  Description:
CC###    NH_FIT(1,nhx) stores the fit variable number njj associated with
CC###    the nhx.  NH_FIT(2,nhx) stores the variable number nhj
CC###    associated with the nhx.

C#### Variable: NJ_FIT(nhj,njj)
C###  Type: INTEGER
C###  Set_up: IPFIT
C###  Description:
C###    NJ_FIT(nhj,njj) is the nj variable that is being fitted for
C###    variable nhj of fit variable njj.

C#### Variable: NKB(ido1,ido2,ido3,nn,nbf)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NKB(ido1,ido2,ido3,nn,nbf) is the local derivative number nk for
C###    derivative orders ido1-1,ido2-1,ido3-1 in each xi direction at
C###    local node nn in basis nbf.  (Inverse is IDO.)

C#### Variable: NLH_FIT(nhj,3,njj)
C###  Type: INTEGER
C###  Set_up: IPFIT
C###  Description:
C###    NLH_FIT(nhj,1,njj) is the nj number in which to store the fit
C###    for variable nhj of fit variable njj.
C###    NLH_FIT(nhj,2,njj) is the field/geometric variable number for
C###    variable nhj of fit variable njj.
C###    NLH_FIT(nhj,3,njj) is the `actual nhx' variable number
C###    for variable nhj of fit variable njj.

C#### Variable: NL_AERO(no_aero,i)
C###  Type: INTEGER
C###  Set_up: IPAERO
C###  Description:
C###    NL_AERO(0,1) is #global lines on upper surface of aerofoil.
C###    NL_AERO(0,2) is #global lines on lower surface of aerofoil.
C###    NL_AERO(no_aero,1) are flow field line #s on upper surface.
C###    NL_AERO(no_aero,2) are flow field line #s on lower surface.

C#### Variable: NL_EXIT(no_exit,i)
C###  Type: INTEGER
C###  Set_up: IPAERO
C###  Description:
C###    NL_EXIT(0,1) is #global lines on downstream outflow face.
C###    NL_EXIT(no_exit,1) are line numbers on downstream outflow face.

C#### Variable: NL_WAKE(no_wake,i)
C###  Type: INTEGER
C###  Set_up: IPAERO
C###  Description:
C###    NL_WAKE(0,1) is #global lines on upper surface of wake.
C###    NL_WAKE(0,2) is #global lines on lower surface of wake.
C###    NL_WAKE(no_wake,1) are line #s on upper surface.
C###    NL_WAKE(no_wake,2) are line #s on lower surface.

C#### Variable: NMGT
C###  Type: INTEGER
C###  Set_up: IPGRID,IPSOLV
C###  Description:
C###    NMGT is the number of multigrid levels.

C#### Variable: NNB(inp1,inp2,inp3,nbf)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NNB(inp1,inp2,inp3,nbf) is the local node number nn in basis nbf
C###    corresponding to the node postion indices inp1,inp2,inp3 in
C###    each xi direction.  (Inverse is INP.)

C#### Variable: NPNODE(0:nonode,0:nr)
C###  Type: INTEGER
C###  Set_up: IPNODE
C###  Description:
C###    NPNODE(0,nr) is the total number of nodes in region nr.
C###    NPNODE(nonode,nr), nonode=1..NPNODE(0,nr) are the node numbers.
C###    NPNODE(0,0) is the total number of nodes defined in all regions.
C###  See-Also: NEELEM,NFFACE,NLLINE,NPNODE

C#### Variable: NP1OPT(no)
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    NP1OPT(no) is used for node numbers (np) of optimisation
C###    variables no.

C#### Variable: NP2OPT(no)
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    NP2OPT(no) is used for node numbers (np) of optimisation
C###    variables no.

C#### Variable: NP3OPT(no)
C###  Type: INTEGER
C###  Set_up: IPOPTI
C###  Description:
C###    NP3OPT(no) is used for node numbers (np) of optimisation
C###    variables no.

C#### Variable: NP_aero_LE
C###  Type: INTEGER
C###  Set_up: IPAERO
C###  Description:
C###    NP_aero_LE is node number of leading edge of aerofoil.

C#### Variable:NP_aero_TE1
C###  Type:INTEGER
C###  Set_up: IPAERO
C###  Description:
C###    NP_aero_TE1 is node number of trailing edge(upper surface node).

C#### Variable:NP_aero_TE2
C###  Type: INTEGER
C###  Set_up: IPAERO
C###  Description:
C###    NP_aero_TE2 is node number of trailing edge(lower surface node).

C#### Variable: NP_ENTRY(0:1..)
C###  Type: INTEGER
C###  Set_up: IPAERO
C###  Description:
C###    NP_ENTRY(0) is  the number of nodes on entry flow face,
C###    NP_ENTRY(1..) are the nodes on entry flow face.

C#### Variable: NP_MOTION
C###  Type: INTEGER
C###  Set_up: IPMOTI
C###  Description:
C###    NP_MOTION is node number for applying motion.

C#### Variable: NP_WAKE(1..,1)
C###  Type: INTEGER
C###  Set_up: IPAERO
C###  Description:
C###    NP_WAKE(1..,1) are nodes on wake surface.

CC#### Variable: NQGE(ng,ne,nb)
CC###  Type: INTEGER
CC###  Set_up: IPGRID
CC###  Description:
CC###    NQGE(ng,ne,nb) is:
CC###    for the extended basis: the global grid pt number nq for
CC###            local gauss point ng in element ne.
CC###    for all other bases   : the closest grid point in xi-space
CC###            to gauss point ng in element ne.
CC###    Note: should look to formulate a reverse mapping at some
CC###          stage.

C#### Variable: NQR(0:2,99)
C###  Type: INTEGER
C###  Set_up: IPGRID,CREATE_LATTICE
C###  Description:
C###    NQR(1,nr) contains the lowest grid point number in region nr
C###    NQR(2,nr) contains the highest grid point number in region nr

C#### Variable: NRE(ne)
C###  Type: INTEGER
C###  Set_up: IPELEM,IPMESH
C###  Description:
C###    NRE(ne) is the region nr number for element ne.

C#### Variable: NSB(nk,nn,nbf)
C###  Type: INTEGER
C###  Set_up: IPBASE
C###  Description:
C###    NSB(nk,nn,nbf) is the ns number (element dofs loop variable)
C###    corresponding to local derivative nk at local node nn in basis
C###    nbf.

C#### Variable: NTACTV
C###  Type: INTEGER
C###  Set_up: IPACTI
C###  Description:
C###    NTACTV is the number of dynamic terms in the active material
C###    response function.

C#### Variable: NTIME_INTERP(ntv)
C###  Type: INTEGER
C###  Set_up: IPTIME
C###  Description:
C###    <HTML> <PRE>
C###    For a given time variable, NTIME_INTERP will return an
C###    integer value which corresponds to the type of interpolation
C###    used with the time variable.
C###      1 = Constant
C###      2 = Linear Lagrange
C###      3 = Quadratic Lagrange
C###      4 = Cubic Lagrange
C###      5 = Cubic Hermite
C###      </PRE> </HTML>

C#### Variable: NTIME_POINTS(ntv)
C###  Type: INTEGER
C###  Set_up: IPTIME
C###  Description:
C###    For a given time variable, NTIME_POINTS will return the
C###    number of time points which have been set for that time
C###    variable.

C#### Variable: NT_PLIN
C###  Type: INTEGER
C###  Set_up: IPPLIN
C###  Description:
C###    NT_PLIN is the total number of polylines.

C#### Variable: NT_PLIN_SECTIONS(index_plin)
C###  Type: INTEGER
C###  Set_up: IPPLIN
C###  Description:
C###    NT_PLIN_SECTIONS(index_plin) is the number of sections to the
C###    polyline.

C#### Variable: NT_SUB_DIV
C###  Type: INTEGER
C###  Set_up: IPEXPO
C###  Description:
C###    NT_SUB_DIV is the number of local subdivisions per element for
C###    exporting to MAP3D.

C#### Variable: NUMLEADS
C###  Type: INTEGER
C###  Set_up: IPLEAD
C###  Description:
C###    NUMLEADS stores the number of electrocardiographic leads.

C#### Variable: NUM_FIT(0:njj)
C###  Type: INTEGER
C###  Set_up: IPFIT
C###  Description:
C###    NUM_FIT(0) is the total number of fit variables defined.
C###    NUM_FIT(njj) is the number of variables (nhj) for each fit
C###      variable njj.

C#### Variable: NUM_GMRES_ORTHOG(nx)
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    NUM_GMRES_ORTHOG(nx) is #stored orthogonalisations for the
C###    generalised minimum residual solver for problem class nx.

C#### Variable: NUM_INITS
C###  Type: INTEGER
C###  Set_up: IPCELL
C###  Description:
C###    NUM_INITS is the number of initial conditions being passed
C###    from CELL to CMISS through the IPCELL routine. There should
C###    be one initial value for each solution variable.

C#### Variable: NUM_PARAMS
C###  Type: INTEGER
C###  Set_up: IPCELL
C###  Description:
C###    NUM_PARAMS is the number of electrical parameters being
C###    passed from CELL to CMISS through the IPCELL routine

C#### Variable: NVJE(nn,nbf,nj,ne)
C###  Type: INTEGER
C###  Set_up: IPELEM
C###  Description:
C###    NVJE(nn,nbf,nj,ne) is the global version number used for
C###    coordinate nj of element vertex nn, family basis nbf
C###    in element ne.  For example: different versions for multiple
C###    thetas or corners for BEM problems.

C#### Variable: NVJF(nn,nbf,nj)
C###  Type: INTEGER
C###  Set_up: CALC_FACE_INFORMATION_IND
C###  Description:
C###    NVJF(nn,nbf,nj) is the global version number used for
C###    coordinate nj of face vertex nn, family basis nbf.

C#### Variable: NVJP(nj,np)
C###  Type: INTEGER
C###  Set_up: IPNODE
C###  Description:
C###    NVJP(nj,np) is the number of versions of a geometric variable
C###    nj at node np.  For example: different versions for multiple
C###    thetas or corners for BEM problems.

C#### Variable: NAQ(nq,na)
C###  Type: INTEGER
C###  Set_up: CONSTRUCTNAQ
C###  Description:
C###    <HTML> <PRE>
C###    NAQ(nq,na+1) is used for interpolating fine grid na from
C###    coarse grid na+1.
C###    NAQ(nq,na) is 0 if nq belongs to grid level na
C###    NAQ(nq,na+1) is
C###       1 if nq is between 2 grid na points in Xi(1) direction
C###       2 if nq is between 2 grid na points in Xi(2) direction
C###       3 if nq is between 2 grid na points in Xi(3) direction
C###       4 if nq is between 4 grid na points in Xi1,2 plane
C###       5 if nq is between 4 grid na points in Xi2,3 plane
C###       6 if nq is between 4 grid na points in Xi3,1 plane
C###       7 if nq is between 8 grid na points in 1,2,3 space
C###      -1 if nq does not belong to grids na or na+1
C###    </PRE> </HTML>

C#### Variable: NW(ne,3)
C###  Type: INTEGER
C###  Set_up: IPEQUA
C###  Description:
C###    <HTML>
C###    NW(ne,1) is the type number (= 1,2..12) for element ne.
C###    NW(ne,2) (BEM) is the type of interpolation used in element ne.
C###    The variable NW(ne,2) has been set according to the following
C###    list.  This variable determines what sort of basis function
C###    has been defined in element ne and hence what sort of special
C###    integration scheme is used when the Green's function
C###    singularity lies in or near the element ne. NW(ne,3) is
C###    used only for boundary element problems where NW(ne,3)=1 if
C###    the normal direction for element ne has been reversed
C###    (=0 otherwise).
C###    <PRE>
C###    NW(ne,2) is:
C###       1 for 1D or 2D constant dep variable interpolation
C###       2 for 1D linear dependent variable interpolation
C###       3 for 1D quadratic dependent variable interpolation
C###       4 for 1D cubic (Lagrange) dep variable interpolation
C###       5 for 1D cubic Hermite dep variable interpolation
C###       6 for 2D bilinear interpolation
C###       7 for 2D biquadratic interpolation
C###       8 for 2D bicubic (Lagrange) interpolation
C###       9 for 2D bicubic Hermite interpolation
C###      10 for 2D linear-quadratic interpolation
C###      11 for 2D linear-cubic Hermite interpolation
C###      12 for 2D quadratic-cubic Hermite interpolation
C###      13 for 2D cubic Hermite-linear interpolation
C###      14 for 2D cubic Hermite simplex (apex at node 1)
C###      15 for 2D cubic Hermite simplex (apex at node 3)
C###      16 for 2D 3 noded sectors
C###    </PRE> <HTML>

C#### Variable: NWQ(1:6,0:nq,na)
C###  Type: INTEGER
C###  Set_up: IPGRID,LATTICE_NWQ
C###  Description:
C###    <HTML> <PRE>
C###    NWQ(1,nq,na) is:
C###       0 if nq is internal
C###       mq1 if nq is on external boundary (for calc noflux bc)
C###       (mq1= 1st adjacent point in flux bc)
C###       If lattice grid method is being used NWQ(1,nq,na) stores
C###       the type of boundary nq lies upon:
C###       1,7,8,11,12,19,20,21,22 low i boundary
C###       2,9,10,13,14,23,24,25,26 high i boundary
C###       3,7,9,15,16,19,20,23,24 low j boundary
C###       4,8,10,17,18,21,22,25,26 high j boundary
C###       5,11,13,15,17,19,21,23,25 low k boundary
C###       6,12,14,16,18,20,22,24,26 high k boundary
C###    NWQ(2,nq,na) is:
C###       0 if nq is internal
C###       mq2 if nq is on external boundary (for calc noflux bc)
C###       (mq2= 2nd adjacent point in flux bc)
C###    NWQ(4,nq,na) is:
C###      -1 if nq was previously active and is now inactivated
C###       0 if nq is not active
C###       1 if nq is active
C###     >=2 if nq is surrounded by active points
C###    NWQ(5,nq,na) is:
C###       1 if nq is fixed value (for multigrid)
C###       2 if nq is fixed flux
C###    NWQ(5,nq,na) is:
C###       active list of points (for DTAR, bidomain)
C###       (5,0,na) is number of points in active list
C###    </PRE> </HTML>

C#### Variable: NXQ(-ni:ni,0:i,0:nq,na)
C###  Type: INTEGER
C###  Set_up: IPGRID
C###  Description:
C###    NXQ(-3..0..3,0:i,nq,na) are neighbouring points of
C###    grid point nq for multigrid level na. eg NXQ(-1,1,nq,na)
C###    is grid# in -ve Xi1 dir.n. NXQ(-1,0,nq,na)=#neighbouring pts.
C###    eg If NXQ(-1,0,nq,na)=2 then NXQ(-1,1,nq,na) & NXQ(-1,2,nq,na)
C###    store the 2 neighbouring points. NOTE: If update grid connectivity is
C###    called, then breaks in the intracellular domain are represented
C###    by giving the relevant grid numbers a negative sign.

C#### Variable: PRECON_CODE(nx)
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    PRECON_CODE(nx) determines the preconditioning applied for the
C###    iterative solver for problem class nx.
C###    The values are 0/1 for None/Diagonal preconditioning.

C#### Variable: SOLVEROPTION(nx)
C###  Type: INTEGER
C###  Set_up: IPSOLV
C###  Description:
C###    SOLVEROPTION(nx) controls what solution procedure is to be used.
C###    The values are 1/2/3/4/5 for LU Factorisation/SVD Factorisation
C###    /Biconjugate gradient/Generalised minimum Residual/
C###    Least Squares.

C#### Variable: SPARSEGKK(nx)
C###  Type: INTEGER
C###  Set_up: IPSOLU
C###  Description:
C###    SPARSEGKK(nx) is 0 if the global soln matrix GKK is stored as
C###    a fully population matrix, 1 if it is stored as a compressed
C###    row sparse matrix, 2 if it is stored as a row-column sparse
C###    matrix, 3 if it is stored in a compressed column sparse matrix,
C###    and 4 if it is stored as a sorted row-column sparse matrix.
C###  See-Also: SPARSITY STRUCTURES

CC AJPs 191297
C#### Variable: TRSF_NR_FIRST
C###  Type: INTEGER
C###  Set_up: IPTRSF
C###  Description:
C###    TRSF_NR_FIRST is the region number of the first/inner surface
C###    used in the transfer matrix construction.

C#### Variable: TRSF_NRLIST
C###  Type: INTEGER
C###  Set_up: IPTRSF
C###  Description:
C###    TRSF_NRLIST(0) is the total number of regions used in the
C###    transfer matrix construction. These regions numbers are stored
C###    in TRSF_NRLIST(1..trsf_nrlist(0)).

C#### Variable: TRSF_NR_SECOND
C###  Type: INTEGER
C###  Set_up: IPTRSF
C###  Description:
C###    TRSF_NR_SECOND is the region number of the second surface used
C###    inthe transfer matrix construction.

C#### Variable: TRSF_NR_OUTER
C###  Type: INTEGER
C###  Set_up: IPTRSF
C###  Description:
C###    TRSF_NR_OUTER is the region number of the outer surface used in
C###    the transfer matrix construction. NB.  Unless a double layer
C###    transfer matrix from transmembrane to epicaridal potentials is
C###    requested, TRSF_NR_OUTER and TRSF_NR_SECOND will be the same.



C*** Real *8


C#### Variable: ACOEFF(nactv)
C###  Type: REAL*8
C###  Set_up: IPACTI
C###  Description:
C###    ACOEFF(nactv),nactv=1,NTACTV are coefficients for linear
C###    dynamic terms.

C#### Variable: ALFA(nactv)
C###  Type: REAL*8
C###  Set_up: IPACTI
C###  Description:
C###    ALFA(nactv),nactv=1,NTACTV are time constants for linear
C###    dynamic terms.

C#### Variable: BINARYQUASIFILE
C###  Type: REAL*8
C###  Set_up: IPQUAS
C###  Description:
C###    BINARYQUASIFILE determines whether the quasi-static output
C###    file is stored as a binary file (=1) or an ascii file (=0).

C#### Variable: CE(nm,ne,nx)
C###  Type: REAL*8
C###  Set_up: IPMATE
C###  Description:
C###    CE(nm,ne,nx) are the values of the piecewise constant material
C###    parameters nm in element ne for problem type nx.

C#### Variable: CG(nm,ng)
C###  Type: REAL*8
C###  Set_up: IPMATE
C###  Description:
C###    CG(nm,ng) are the values of the material parameters nm at Gauss
C###    points ng.

C#### Variable: CGE(nm,ng,ne,nx)
C###  Type: REAL*8
C###  Set_up: IPMATE
C###  Description:
C###    CGE(nm,ng,ne,nx) are the values of the material parameters nm defined at Gauss
C###    points ng in element ne for problem nx.

C#### Variable: CP(nm,np,nx)
C###  Type: REAL*8
C###  Set_up: IPMATE
C###  Description:
C###    CP(nm,np,nx) are the values of the piecewise linear material
C###    parameters nm at node np for problem type nx.

C#### Variable: CPDST(NCONMX)
C###  Type: REAL*8
C###  Set_up: IPCOUP
C###  Description:
C###    CPDST(nq) contains the distance between a grid point
C###    from one region and the closest grid point in another region.

C#### Variable: DEL_T
C###  Type: REAL*8
C###  Set_up: IPACTI
C###  Description:
C###    DEL_T is time step used (taken from load stepping loop in FE07).

C#### Variable: DET(nbf,0:nn,ng,1..6)
C###  Type: REAL*8
C###  Set_up: IPBASE
C###  Description:
C###    DET is the Jacobian of the transformation from the ith split
C###    part of the BEM element (2d integrals) to (psi1',psi2')
C###    coordinates when the singularity is located at local node nn.

C#### Variable: dNUdXQ(ni,nj,nq)
C###  Type: REAL*8
C###  Set_up: UPGRID
C###  Description:
C###    dNUdXQ(ni,nj,nq) are derivatives of nu coordinatess wrt X
C###    coordinates at nq.

C#### Variable: dXdXiQ(nj,ni,nq)
C###  Type: REAL*8
C###  Set_up: UPGRID
C###  Description:
C###    dXdXiQ(nj,ni,nq) are derivatives of X coordinates wrt Xi at nq
C###    derived from local quadratic basis functions.

C#### Variable: dXdXiQ2(nj,ni,nq)
C###  Type: REAL*8
C###  Set_up: UPGRID
C###  Description:
C###    dXdXiQ(nj,ni,nq) are derivatives of X coordinates wrt Xi at nq.
C###    derived from global FE basis functions.

C#### Variable: FEXT(NIFEXTM,ng,ne)
C###  Type: REAL*8
C###  Set_up: IPACTI
C###  Description:
C###    <HTML> <PRE>
C###    FEXT(1,ng,ne) is current muscle fibre extension ratio
C###      "  2    "   "  previous   "     "       "       "
C###      "  3    "   "  muscle fibre ext. ratio at time of activation
C###      "  4    "   "  [Ca]i
C###      "  5    "   "  previous hereditary integral for 1st time constant
C###      "  6    "   "     "         "         "     "   2nd  "      "
C###      "  7    "   "     "         "         "     "   3rd  "      "
C###    </PRE> <HTML>

C#### Variable: FLOW_COEFFS(no_coeffs)
C###  Type: REAL*8
C###  Set_up: IPMOTI
C###  Description:
C###    FLOW_COEFFS(no_coeffs),no_coeffs=1,nt_coeffs are Fourier
C###    coefficients for flow.

C#### Variable: GCHQ(nk,nq)
C###  Type: REAL*8
C###  Set_up: UPGRID
C###  Description:
C###    GCHQ(nk,nq) are components of GUQ(i,j) *
C###    Christoffel_symbol(i,j,k).

C#### Variable: GUQ(ni,nj,nq)
C###  Type: REAL*8
C###  Set_up: UPGRID
C###  Description:
C###    GUQ(ni1,ni2,nq) are contravariant components of metric tensor
C###    at nq.

C#### Variable: ITERATIVE_TOL(nx)
C###  Type: REAL*8
C###  Set_up: IPSOLV
C###  Description:
C###    ITERATIVE_TOL(nx) is the absolute tolerance to use for the
C###    iterative solvers for problem class nx.

C#### Variable: LEADCOUP(0:nleadelec,nlead)
C###  Type: REAL*8
C###  Set_up: IPLEAD
C###  Description:
C###    LEADCOUP(0:nleadelec,nlead) is the coefficient on contribution
C###    that lead electrode nleadelec makes to lead nlead.
C###    LEADCOUP(0,nlead,nr) is an additive constant that can be
C###    applied to each lead (if necessary).

C#### Variable: OUTPUT_SOLTIMES
C###  Type: REAL*8
C###  Set_up: IPQUAS
C###  Description:
C###    OUPUT_SOLTIMES determines whether or not timing information is
C###    output for quasi-static solutions. Note: this is addition to
C###    IWRIT4.

C#### Variable: PAOPTI(noopti)
C###  Type: REAL*8
C###  Set_up: IPOPTI
C###  Description:
C###    PAOPTI(noopti) are initial values.

C#### Variable: PBOPTI
C###  Type: REAL*8
C###  Set_up: IPOPTI
C###  Description:
C###    PBOPTI(noopti) is used later to hold previous value of PAOPTI.

C#### Variable: PF(2,ne)
C###  Type: REAL*8
C###  Set_up: IPINIT
C###  Description:
C###    PF(1..2,ne) is the total pressure on element face.

C#### Variable: PG(ns,nu,ng,nb)
C###  Type: REAL*8
C###  Set_up: IPBASE
C###  Description:
C###    PG are the global basis function values evaluated
C###    at the Gauss points.

C#### Variable: PGNQE(ng,nq,nqsc)
C###  Type: REAL*8
C###  Set_up: GEN_GRID_MAP
C###  Description:
C###    PGNQE is the grid to gauss point mapping array. It is designed
C###    so for any transfer of data only a dot product is required.
C###    For a gauss point ng, this array gives the weighting of each
C###    local grid point in the current scheme.

C#### Variable: PLIN_DATA(nj,no_point,index_plin)
C###  Type: REAL*8
C###  Set_up: IPPLIN
C###  Description:
C###    PLIN_DATA(nj,no_point,index_plin),nj=1,njt are the coordinates
C###    of each point of the polyline.

C#### Variable: PMAX(noopti)
C###  Type: REAL*8
C###  Set_up: IPOPTI
C###  Description:
C###    PMAX(noopti) is maximum parameter value allowed.

C#### Variable: PMIN(noopti)
C###  Type: REAL*8
C###  Set_up: IPOPTI
C###  Description:
C###    PMIN(noopti) in minimum parameter value allowed.

C#### Variable: QUASIPROB
C###  Type: REAL*8
C###  Set_up: IPQUAS
C###  Description:
C###    QUASIPROB is the type of quasi-static problem. The values are
C###    1/2 for changing source/changing rhs.

C#### Variable: QUASITIMESTEP
C###  Type: REAL*8
C###  Set_up: IPQUAS
C###  Description:
C###    QUASITIMESTEP controls the type of time stepping for
C###    quasi-static problems. The values are 1/2 for fixed/automatic
C###    time stepping.

C#### Variable: QUASIT0
C###  Type: REAL*8
C###  Set_up: IPQUAS
C###  Description:
C###    QUASIT0 is the initial time for quasi-static problems.

C#### Variable: QUASIT1
C###  Type: REAL*8
C###  Set_up: IPQUAS
C###  Description:
C###    QUASIT1 is the final time for quasi-static problems.

C#### Variable: QUASITINCR
C###  Type: REAL*8
C###  Set_up: IPQUAS
C###  Description:
C###    QUASITINCR is the time step for quasi-static problems (initial
C###    if automatic stepping).

C#### Variable: REG_PARAM_LAPLACE
C###  Type: REAL*8
C###  Set_up: IPINVE
C###  Description:
C###    Regularisation parameter for surface laplacians.

C#### Variable: RMIN_SVD
C###  Type: REAL*8
C###  Set_up: IPINVE
C###  Description:
C###    RMIN_SVD is the cutoff level for the singular values if SVD is
C###    used to regularise.

C#### Variable: SNLPA
C###  Type: REAL*8
C###  Set_up: IPACTI
C###  Description:
C###    SNLPA  is static nonlinearity parameter "a".

C#### Variable: TIKH_VALUE
C###  Type: REAL*8
C###  Set_up: IPINVE
C###  Description:
C###    TIKH_VALUE is the Tikhonov regularisation parameter.

C#### Variable: TIME_VALUES(i,ntp,ntv)
C###  Type: REAL*8
C###  Set_up: IPTIME
C###  Description:
C###    <HTML> <PRE>
C###    This variable stores time and value information for all
C###    time variables. Given a time point ntp and a time variable
C###    ntv :
C###
C###    TIME_VALUES(2,0,ntv) = value before first time point
C###    TIME_VALUES(2,NTIMEPOINTSM+1,ntv = value after last time point
C###
C###    TIME_VALUES(1,ntp,ntv) = the time where point ntp is set
C###    TIME_VALUES(2,ntp,ntv) = the value of the variable at that time
C###    </PRE> </HTML>

C#### Variable: TV_SLO
C###  Type: REAL*8
C###  Set_up: IPACTI
C###  Description:
C###    TV_SLO is slope of force/velocity in stretching (before yield).

C#### Variable: WG(ng,nb)
C###  Type: REAL*8
C###  Set_up: IPBASE
C###  Description:
C###    WG are the global basis function Gauss point weights.

C#### Variable: WU(0:NUM+1,ne)
C###  Type: REAL*8
C###  Set_up: IPFIT
C###  Description:
C###    WU(nu,ne) are smoothing function weights for derivative nu on
C###    element ne.  WU(0,ne) is global weight scale factor for the
C###    element.  WU(NUM+1,ne) is the Sobolev value for that element.

C#### Variable: XA(na,nj,ne)
C###  Type: REAL*8
C###  Set_up: IPELEM
C###  Description:
C###    XA(na,nj,ne) are element based reference variable parameters.

C#### Variable: XIG(ni,ng,nb)
C###  Type: REAL*8
C###  Set_up: IPBASE
C###  Description:
C###    XIG is Xi position of Gauss point ng in global basis nb.

C#### Variable: XP(nk,nv,nj,np)
C###  Type: REAL*8
C###  Set_up: IPFIBR,IPMESH1,IPMESH2,IPMESH3,IPMESH5,IPNODE,IPSHEE
C###  Description:
C###    XP(nk,nv,nj,np) is the Xj position (nk=1) or derivative (nk>1)
C###    for version nv of coordinate nj at node np.

C#### Variable: YG(niyg,ng,ne)
C###  Type: REAL*8
C###  Set_up: IPGAUS
C###  Description:
C###    <html>
C###    <p>YG(niyg,ng,ne) is the Gauss point array.</p>
C###    <p>For nonlinear elasticity problems with cell-based active
C###       contraction:<br>
C###       <code>YG(1,ng,ne)</code> is expected to contain the active
C###         tension values computed from the cellular model, usually
C###         set with the command:
C###         <code>fem update gauss gridvars $T_GRIDARRAY $T_GRIDIDX YG 1</code><br>
C###       <code>YG(?,ng,ne)</code> is used as workspace to hold the last converged ..... command <code>fem update ....</code> XXXXXXXXXXXXX<br>
C###       <br>
C###       Need to be very careful when it comes to <code>fem update gauss stress/strain</code>
C###       not to overwrite this stuff.....
C###    </p>
C###    </html>

C#### Variable: YIELDR
C###  Type: REAL*8
C###  Set_up: IPACTI
C###  Description:
C###    YIELDR is ratio of yield tension to isometric tension.


C*** CHARACTERS


C#### Variable: LEADTITLE(nlead)
C###  Type: CHARACTER
C###  Set_up: IPLEAD
C###  Description:
C###    LEADTITLE is the title for lead nlead.

C#### Variable: SIGFNAME
C###  Type: CHARACTER
C###  Set_up: IPLEAD
C###  Description:
C###    SIGFNAME is the filename for the electrode signal file.

C#### Variable: TIME_VARIABLE_NAMES(ntv)
C###  Type: CHARACTER
C###  Set_up: IPTIME
C###  Description:
C###    For a given time variable, this array stores the unique
C###    name given to the time variable.


C*** Logicals



C#### Variable: ADAPINT
C###  Type: LOGICAL
C###  Set_up: IPSOLV
C###  Description:
C###    ADAPINT is .TRUE. if adaptive integration is to be
C###    used in BE routines.

C#### Variable: CALC_GLOBAL(nr,nx)
C###  Type: LOGICAL
C###  Set_up: IPEQUA
C###  Description:
C###    CALC_GLOBAL(nr,nx) controls whether or not the global matrices
C###    GK, GQ etc. are assembled or whether the problem is assembled
C###    directly into the solution matrices GKK etc.

C#### Variable: ETYP(ie)
C###  Type: LOGICAL
C###  Set_up: IPEQUA
C###  Description:
C###    ETYP(ie),ie=1,12 is .true. if element type ie (see TITLE2) used.

C#### Variable: FIX(ny,niy,nx)
C###  Type: LOGICAL
C###  Set_up: IPBOUN,IPFIT,IPMOTI,IPINI3,4,5,&9,IPMAT3,4,&5
C###  Description:
C###

C#### Variable: HERMITE
C###  Type: LOGICAL
C###  Set_up: IPSOLV
C###  Description:
C###    HERMITE is .TRUE. if hermite interpolation is used
C###    in a BE solution.

C#### Variable: LUMP
C###  Type: LOGICAL
C###  Set_up: IPSOLV
C###  Description:
C###    LUMP is .TRUE. if mass lumping is used.

C#### Variable: PRESSCORR
C###  Type: LOGICAL
C###  Set_up: IPSOLV
C###  Description:
C###    PROMPT is .TRUE. if there is to be a pressure correction
C###    equation

C#### Variable: PROMPT
C###  Type: LOGICAL
C###  Set_up: IPSOLV
C###  Description:
C###    PROMPT is .TRUE. if computations await prompt after
C###    IWRIT1(nr,nx) steps.

C#### Variable: RHIECHOW
C###  Type: LOGICAL
C###  Set_up: IPSOLV
C###  Description:
C###    RHIECHOW is .TRUE. if Rhie Chow interpolation is to be
C###    used


C AJPs 2nd July 1998
C#### Variable: TORSO_LENGTHS(0:2,99)
C###  Set_up: IPCUST
C###  Description:
C###    Used for customising torso.
C###    TORSO_LENGTHS(1,n) is height of landmark n on generic torso
C###    TORSO_LENGTHS(2,n) is height of landmark n on current torso

C#### Variable: SCLTYPE
C###  Set_up: IPCUST
C###  Description:
C###    Type of scaling used in torso customisation
C###    SCLTYPE=1 is simple scaling
C###    SCLTYPE=2 is varibale length scaling

C#### Variable: NTL
C###  Set_up: IPCUST
C###  Description:
C###    Number of length measurements.
C AJPe 2nd July 1998

C#### Module: FE12
C###  Description:
C###    Routines for output of model parameters.

C###  Routine: OPACTI  output active muscle model parameters
C###  Routine: OPAERO   output aerofoil parameters
C###  Routine: OPANAL   output analytic formula parameters
C###  Routine: OPBASE   output basis functions
C###  Routine: OPBASE1  output basis function nb
C###  Routine: OPCELL   output cell parameters
C###  Routine: OPCELL_PROMPT   output cell parameters
C###  Routine: OPCONS   *REMOVED* output constant
C###  Routine: OPCOOR   output coordinate data
C###  Routine: OPCORN   output corner node data (for BE problems)
C###  Routine: OPCOUP   output coupling data
C###  Routine: OPCUST   output customisation parameters
C###  Routine: OPDATA   output data
C###  Routine: OPEIGE   output eigenvalues
C###  Routine: OPELEM   output element topology for specified elements
C###  Routine: OPELEM_INTERFACE  output element interfaces
C###  Routine: OPELEMD  output element topology for def element ne
C###  Routine: OPELEM1  output element topology for undef element ne
C###  Routine: OPELEMG  output element groups
C###  Routine: OPEQUA   output equation parameters
C###  Routine: OPEXPO   output export parameters
C###  Routine: OPFACE   output face data
C###  Routine: OPFACE1  output one face
C###  Routine: OPFACEG  output face groups
C###  Routine: OPFIT    output fit
C###  Routine: OPFUNC   output objective function
C###  Routine: OPGAUSG  output Gauss point groups
C###  Routine: OPGAUSX  output Gauss point array XG
C###  Routine: OPGAUSY  *REMOVED* output Gauss point array YG
C###  Routine: OPGRID   output finite difference grid points
C###  Routine: OPGRIDG  output grid point groups
C###  Routine: OPGROW   output growth law parameters
C###  Routine: OPHEAD   output heading
C###  Routine: OPHIST   output history
C###  Routine: OPIMPO   output import parameters
C###  Routine: OPINCR   output increment vectors
C###  Routine: OPINIT   output initial & boundary conditions
C###  Routine: OPINVE   output inverse transfer matrix parameters
C###  Routine: OPITER   output iteration parameters
C###  Routine: OPLEAD   output electrocardiographic leads
C###  Routine: OPLINE   output line segments
C###  Routine: OPLINEG  output line groups
C###  Routine: OPMAP    output mapping arrays
C###  Routine: OPMATE   output materials
C###  Routine: OPMATR   output matrix
C###  Routine: OPMODA   output modal values
C###  Routine: OPMOTI   output motion parameters
C###  Routine: OPNODE   output nodal coordinates for all nodes
C###  Routine: OPNODE1  output nodal coordinates for node np
C###  Routine: OPNOIS   output noise parameters
C###  Routine: OPNODEG  output node groups
C###  Routine: OPNORM   output normal reversals
C###  Routine: OPOBJE   output objects
C###  Routine: OPOPTI   output optimisation parameters
C###  Routine: OPOUTP   output output
C###  Routine: OPPARA   output parameters
C###  Routine: OPPCAP   output pulmonary capillary model results      
C###  Routine: OPPLIN   output polyline parameters
C###  Routine: OPPLING  output polyline groups
C###  Routine: OPREGI   output regions data
C###  Routine: OPREGP   output regularisation parameters
C###  Routine: OPREFE   output reference node(s)/location(s)
C###  Routine: OPSEGM   output segments
C###  Routine: OPSING   output singularity location (BE problems)
C###  Routine: OPSIGN   output signal
C###  Routine: OPSOLV   output solution parameters
C###  Routine: OPSOUR   output source parameters
C###  Routine: OPSTRA   output strains
C###  Routine: OPSTRE   output stresses
C###  Routine: OPST80   *** ARCHIVED ***
C###  Routine: OPTCAP   output cell transit times through pulm capillaries
C###  Routine: OPTIME   output time variables
C###  Routine: OPTRSF   output transfer matrix parameters
C###  Routine: OPVARI   output YP variables
C###  Routine: OPVOLU   output volumes for regs bdd by selected elements.
C###  Routine: OPVORO   output voronoi cell information
C###  Routine: OPXI     output xi coordinates

C#### Variable: VOL(NBFM)
C###  Type: REAL*8
C###  Set_up: OPELEM1,OPELEMD
C###  Description:
C###    VOL(NBFM) stores the elements length/area/volume

C#### Variable: VOLT(NBFM)
C###  Type: REAL*8
C###  Set_up: OPELEM
C###  Description:
C###    VOLT(NBFM) stores the total length/area/volume

C#### Variable: VOLTC(NBFM)
C###  Type: INTEGER
C###  Set_up: OPELEM1,OPELEMD
C###  Description:
C###    Stores the number of element for each basis function.


C#### Module: FE13
C###  Description:
C###    Mesh Input and Output Routines

C###  Routine: BLK13    sets up morphometric values for meshes
C###  Routine: ALVANGLE Calculates alveolar angles, normal to opening      
C###  Routine: ALVPARA Calculates alveolar parameters for capillary mesh
C###  Routine: AREAFIELD generates data points with cross sec area values
C###  Routine: BRANCH_ORD calculates branch order for lung trees  
C###  Routine: CAP_IO   creates inlet and oulet vessles for cap network
C###  Routine: CAP_NE   defines parameters/arrays for capillary
C###  Routine: CAP_PROJ projects capillary mesh onto alveolar surface
C###  Routine: CAP_PROJ2 projects nodes back out to sphere
C###  Routine: CENTRE_LINE calculates centre line for a tube mesh
C###  Routine: IPMESH   input specialized mesh
C###  Routine: IPMESH1  input rectangular or cuboid mesh
C###  Routine: IPMESH2  input lung mesh parameters
C###  Routine: GENCIRC  generates pulmonary venous and arterial trees.
C###  Routine: GENPCAP  generates pulmonary capillary mesh
C###  Routine: GET_CIRCUM_ELEM calculates a list of circumferential elements
C###  Routine: GN1DNE   GeNerates a 1D element (NE) for a mesh
C###  Routine: GN1DNEJ  GeNerates 1D element (NE) arrays
C###  Routine: GNARTRY  GeNerates pulmonary arterial mesh from airway mesh
C###  Routine: GNART_SUPER GeNerates pulmonary arterial supernumerary branches
C###  Routine: GNBDMESH bifurcating distributive mesh in a host region
C###  Routine: GNMESH1  GeNerate conducting MESH based on IPMESH2
C###  Routine: GNMESH2  GeNerate respiratory MESH based on IPMESH2
C###  Routine: GNNECE   GeNerate NE-based CE array for mesh
C###  Routine: GNNEZD   GeNerate NE-hosted points, store in ZD
C###  Routine: GNSYM_R  GeNerate SYMmetric Respiratory airway tree
C###  Routine: GNVEIN   GeNerates pulmonary venous mesh from airway mesh
C###  Routine: IPMESH3  input eccentric spheres mesh
C###  Routine: IPMESH5  input cylindrical mesh
C###  Routine: IPMESH6  input dynamic memory wrapper for IPMESH6_SUB
C###  Routine: IPMESH6_SUB input coronary mesh
C###  Routine: IPMESH7  input
C###  Routine: IPMESH8_DYNAM  input Purkinje fibre mesh
C###  Routine: IPMESH8  dynamic memory wrapper for IPMESH8_DYNAM
C###  Routine: IPMESH10 input dynamic memory wrapper for IPMESH10_SUB
C###  Routine: IPMESH10_SUB input Coronary surface mesh
C###  Routine: IPMESH11 input voronoi mesh
C###  Routine: MESH_ANGLE calculates angle between parent and daughter
C###  Routine: MESH_ANGLE_CHECK checks branch angle limit
C###  Routine: MESH_BRANCH creates mesh branch from parent towards cofm
C###  Routine: MESH_COFM calculates cofm of a set of data points
C###  Routine: MESH_PLANE_ANGLE  calcs angle between parent and daughter
C###  Routine: MESH_REDUCE reduces a mesh to a specified order
C###  Routine: MESH_REPOINT reassigns data points to closest branch
C###  Routine: MESH_SPLIT divides data points in tow, using splitting plane
C###  Routine: MESH_TOSET makes sure that branch in contact with data set
C###  Routine: OPMESH   output specialized mesh parameters
C###  Routine: OPMESH1  output rectangular or cuboid mesh
C###  Routine: OPMESH2  output lung mesh
C###  Routine: OPMESH3  output eccentric spheres mesh
C###  Routine: OPMESH5  output open cylindrical mesh
C###  Routine: OPMESH6  output coronary mesh
C###  Routine: OPMESH8  output Purkinje fibre mesh
C###  Routine: OPMESH9  output regular mesh with incision
C###  Routine: REMESH   REfines a generated MESH
C###  Routine: REMESH2  reordering routine for REMESH

C#### Variable: AMAP(area)
C###  Type: INTEGER
C###  Set_up: IPMESH8_DYNAM
C###  Description:
C###    This array is dynamically allocated in IPMESH8.
C###    For each branch, AMAP contains the numbers of the two areas
C###    which the branch has split. Indexing with one number provides
C###    the other.

C#### Variable: AREANUM(nrc,NUMPLANES)
C###  Type: INTEGER
C###  Set_up: IPMESH8_DYNAM
C###  Description:
C###    This array is dynamically allocated in IPMESH8.
C###    For each random point, nrc, AREANUM contains the area
C###    number which the point currently lies in. If multiple
C###    branch planes are possible, it contains all the possible
C###    areas the point could be in for the next iteration.
      
C#### Variable: GENM
C###  Type: INTEGER
C###  Set_up: lung00.cmn
C###  Description:
C###  This is the maximum number of airway/pulmonary blood vessel
C###  generations.

C#### Variable: NPNE_ALV(0:NP_NE,NE_R_M)
C###  Type: INTEGER
C###  Set_up: CAP_PROJ
C###  Decription:
C###  This array is used for pulmonary capillary network generation. It
C###  stores which capillary nodes are projected onto which alveolar
C###  elements.
      
C#### Variable: NPQ
C###  Type: INTEGER
C###  Set_up: IPMESH7
C###  Description:
C###    NPQ(NQM) gives an np to nq mapping. It is used to make a fine
C###    resolution FEM mesh with grid points as element nodes.

C#### Variable: NRANDOM
C###  Type: INTEGER
C###  Set_up: IPMESH8_DYNAM
C###  Description:
C###    This paramter is used to specify the number of random points
C###    which are used to generate a tree. It does not relate to the
C###    number of nodes used.

C#### Variable: NUMPLANES
C###  Type: INTEGER
C###  Set_up: IPMESH8_DYNAM
C###  Description:
C###    This parameter stores the number of possible ways to perform
C###    the area/volume splitting at each iteration. In 2-d it will
C###    be set to 1 and in 3-d in will be 3 (x,y,z orthogonal planes)

C#### Variable: PTAR(area)
C###  Type: INTEGER
C###  Set_up: IPMESH8_DYNAM
C###  Description:
C###    This array is dunamically allocated in IPMESH8.
C###    For a given area, PTAR will return the number of the branch
C###    point which created the area by previously splitting a
C###    larger area.

C#### Variable: TCX(0:NP_R_M,2,0:BMAX)
C###  Type: INTEGER
C###  Set_up: IPMESH8_DYNAM
C###  Description:
C###    This array is dunamically allocated in IPMESH8.
C###    This array is the connectivity array to which branch information
C###    is written. The connectivity is nodally based.
C###    The first index is the node number, note that it is the
C###    number within the current region.
C###
C###    TCX(np,1,0) gives the parent node to np
C###    TCX(np,2,0) gives the number of chidren from np
C###    TCX(np,2,1) to TCX(np,2,TCX(np,2,0)) are the node numbers
C###    of the children (if any).

C#### Variable: A_SITE
C###  Type:  REAL*8
C###  Set_up:  IPMESH6
C###  Description:
C###    A_SITE stores the vector in the plane of fibres and sheets
C###    perpendicular to B (see B_SITE) at each free site in
C###    coronary network

C#### Variable: ALPHA_SITE
C###  Type:  REAL*8
C###  Set_up:  IPMESH6
C###  Description:
C###    ALPHA_SITE stores the normailized proportion of A added to
C###    B to achive the desired branch angle at each free site in
C###    coronary network

C#### Variable: B_SITE
C###  Type:  REAL*8
C###  Set_up:  IPMESH6
C###  Description:
C###    B_SITE stores the vector in the plane of fibres and sheets
C###    which is the projection of the parent vector onto that plane
C###    at each free site in the coronary network

C#### Variable: BBM(12,NORM)
C###  Type:  REAL*8
C###  Set_up:  IPMESH2
C###  Description:
C###    BBM stores information about Black Box Acinar models for full
C###    lung models.

C#### Variable: BRANCH_ARRAY
C###  Type:  INTEGER
C###  Set_up:  IPMESH6
C###  Description:
C###    Stores the two nodes at either end of a branch

C#### Variable: FREE_SITES
C###  Type:  INTEGER
C###  Set_up:  IPMESH6
C###  Description:
C###    FREE_SITES(point,1..4) stores node/branch/order
C###    /number of branches at each free sites in a coronary
C###    network


C#### Variable: JTYP14
C###  Type: INTEGER
C###  Set_up: IPMESH
C###  Description:
C###    <HTML>
C###    JTYP14 is mesh type for specialized meshes (1 = regular fractal
C###    tree; 2 = stochastic fractal tree) as follows:
C###    Fractal tree branch parameters:
C###    <UL>
C###    <LI>Ratio_Angle,Ratio_Length,Ratio_Diameter
C###    <LI>B_angle_y(no_gen) : mean angle branch makes with y-axis
C###    <LI>B_angle_xy(no_gen): mean angle branch makes with xy-plane
C###    <LI>B_angle_SD(no_gen): standard deviation of angles<BR>
C###    Note: angles are generated with Normal distribution using mean
C###    and standard deviation obtained by dividing previous generation
C###    values by Ratio_Angle.<BR>
C###    <LI>B_length(no_gen)
C###    <LI>B_diameter(no_gen)
C###    <LI>B_volume(no_gen)
C###    <LI>B_flow(no_gen)   ?? <-- not needed
C###    <LI>NW(ne,1)=generation number
C###    </UL>
C###    </HTML>

C#### Variable: MESH_TYPE
C###  Type: INTEGER
C###  Set_up: IPMESH1
C###  Description:
C###    MESH1_TYPE is 1/2/3 for even spacing/spacing specified by
C###    blocks/spacing specified by position.

C#### Variable: MESH1_S(ni,1:3)
C###  Type: INTEGER
C###  Set_up: IPMESH1
C###  Description:
C###    MESH1_S(ni,0) is total number of elements in s(ni) direction.
C###    MESH1_S(ni,1) is  number of elements in s(ni) direction for
C###    even spacing.  MESH1_S(ni,1..3) are number of elements in 3
C###    blocks in s(ni) direction for uneven mesh spacing.

C#### Variable: MESH3_NB
C###  Type: INTEGER
C###  Set_up: IPMESH3
C###  Description:
C###    <HTML>
C###    MESH3_NB(i,j) are the basis numbers used in the mesh. When j=1
C###    the bases are BE bases. When j=2 the bases are FE bases.
C###    <UL>
C###    <LI>i=1 gives the main basis function number
C###    <LI>i=2 gives the bottom sector (apex node 1) basis number
C###    <LI>i=3 gives the top sector (apex node 3) basis number
C###    </UL>
C###    </HTML>

C#### Variable: NEP(np)
C###  Type:  INTEGER
C###  Set_up:  IPMESH6,DEXI
C###  Description:
C###     NEP(np) is the global element number of a host mesh gobal
C###     node np co-ordinates can be interpolated within

C#### Variable: NPE(0:np_list,ne,nr)
C###  Type:  INTEGER
C###  Set_up:  IPMESH6,DEXI
C###  Description:
C###     NPE(0:np_list,ne,nr)  is the number of nodes in region
C###     nr contained within the host element ne.
C###     NPE(np_list,ne,nr) is the list of node in region nr contained
C###     within the host element ne.

C#### Variable: NSPHERES
C###  Type: INTEGER
C###  Set_up: IPMESH3
C###  Description:
C###    NSPHERES is the number of spheres (or circles).

C#### Variable: POINTS
C###  Type: REAL*8
C###  Set_up:  IPMESH6
C###  Description:
C###    Points is the xi co-ordinates of a points in a local coronary
C###    network

C#### Variable: SEG_STEP
C###  Type: REAL*8
C###  Set_up:  IPMESH6
C###    SEG_STEP is the current step size for a segment between two
C###    branch points

C#### Variable: SEG_SUM
C###  Type:  REAL*8
C###  Set_up:  IPMESH6
C###  Description:
C###    SEG_SUM is the acumulative total distance of a segments
C###    between two branch points

C#### Variable: COFM(nj,NUMPLANES)
C###  Type: REAL*8
C###  Set_up: IPMESH8_DYNAM
C###  Description:
C###    This array is dynamically allocated in IPMESH8.
C###    For a given area, COFM temporarily stores the position of
C###    the centre of mass for the area. If multiple planes are used,
C###    it stores all mass centres for the possible areas.

C#### Variable: MESH1_R(ni)
C###  Type: REAL*8
C###  Set_up: IPMESH1
C###  Description:
C###    MESH1_R(ni) is ratio of coarse to fine spacing in s(ni)
C###    direction. MESH1_COORD(n,nj) are coordinates of corner point n.
C###    MESH1_X(ni,1..) are positions relative to 1st along s(ni).

C#### Variable: MESH3_RAD(nosphere)
C###  Type: REAL*8
C###  Set_up: IPMESH3
C###  Description:
C###    MESH3_RAD(nosphere) is the radius of sphere nosphere (ordered
C###    from smallest to largest).

C#### Variable: MESH3_S(nosphere,1:2)
C###  Type: REAL*8
C###  Set_up: IPMESH3
C###  Description:
C###    MESH3_S(nosphere,1) is the number of elements around sphere
C###    nosphere (ie the number of elements in the theta direction).
C###    MESH3_S(nosphere,2) is the number of elements in the azimural
C###    direction of sphere nosphere.

C#### Variable: POINTCOORD(np,nj)
C###  Type: REAL*8
C###  Set_up: IPMESH8_DYNAM
C###  Description:
C###    This array is dynamically allocated in IPMESH8.
C###    This array provides temporary storage for nodal coordinates.
C###    It is necessary to have the coordinates stored when using
C###    element size smoothing.

C#### Variable: RANDOM_COORD(nrc,nj)
C###  Type: REAL*8
C###  Set_up: IPMESH8_DYNAM
C###  Description:
C###    This array is dynamically allocated in IPMESH8.
C###    RANDOM_COORD stores the coordinates of the random points
C###    generated using the fortran random number generator. The
C###    points may be transformed to suit the problem in which case
C###    only the transformed coordinates are stored.

C#### Variable: SCALE(no)
C###  Type: REAL*8
C###  Set_up: IPMESH2
C###  Description:
C###

C#### Variable: SIDE(ni,nj)
C###  Type: REAL*8
C###  Set_up: IPMESH1
C###  Description:
C###    SIDE(ni,nj) is the side length of the current element.

C#### Variable: SIDE_COARSE(ni,nj)
C###  Type: REAL*8
C###  Set_up: IPMESH1
C###  Description:
C###    SIDE_COARSE(ni,nj) is length of side for coarse elements.

C#### Variable: SIDE_FINE(ni,nj)
C###  Type: REAL*8
C###  Set_up: IPMESH1
C###  Description:
C###    SIDE_FINE(ni,nj) is length of side for fine elements.

C#### Variable: SIDE_LENGTH(ni)
C###  Type: REAL*8
C###  Set_up: IPMESH1
C###  Description:
C###    SIDE_LENGTH(ni) is total length of side ni.

C#### Variable: SIDE_TOT(ni,nj)
C###  Type: REAL*8
C###  Set_up: IPMESH1
C###  Description:
C###    SIDE_TOT(ni,nj) keeps track of the lengths along the sides.

C#### Variable: SIDE_X(nj)
C###  Type: REAL*8
C###  Set_up: IPMESH1
C###  Description:
C###    SIDE_X(nj) is x(nj) distance along a side.

C#### Variable: SIGMA(nosphere)
C###  Type: REAL*8
C###  Set_up: IPMESH3
C###  Description:
C###    SIGMA(nosphere) contains the conductivity inside sphere
C###    nosphere.

C#### Variable: XIP(ni,np)
C###  Type:  REAL*8
C###  Set_up:  IPMESH6, DEXI
C###  Description:
C###     XIP(ni,np) is the local xi positions for global node np
C###     within a host mesh
C###  See-Also: XID


C     rgb initialising NPE_PTR
C     over entire program
C#### Module: FE14
C###  Description:
C###    Routines which call graphics buffer routines.

C###  Routine: INDEX_FILL_AREA  *REMOVED* (fn) bundle index for fill-areas
C###  Routine: INDEX_POLYLINE   (fn) bundle index for polylines
C###  Routine: INDEX_POLYMARKER (fn) bundle index for polymarkers
C###  Routine: INDEX_TEXT       (fn) bundle index for text
C###  Routine: ACWK             activates workstation
C###  Routine: BEZIER           creates Bezier curve
C###  Routine: BEZIER_POINTS    calculates Bezier curve points
C###  Routine: CIRCLE           draws circle
C###  Routine: CLOSE_SEGMENT    closes graphics segment
C###  Routine: CLOSE_WS         closes graphics workstation
C###  Routine: CLWS             closes all workstations
C###  Routine: DAWK             deactivates a workstation
C###  Routine: DBOX             draws box
C###  Routine: DELETE_SEGMENT   deletes graphic segment
C###  Routine: DETECT           change segment detectability
C###  Routine: DOCUM            displays documentation window
C###  Routine: FILL_AREA        draws fill-area
C###  Routine: FILL_LINE        draws line, same scale as FILL_AREA
C###  Routine: LINE3D           *REMOVED* draws line on 3D viewport
C###  Routine: LOCATOR          calls GKS locator
C###  Routine: OPEN_SEGMENT     opens  graphics segment
C###  Routine: PICK             use pick input device
C###  Routine: POLYLINE         draws polyline
C###  Routine: POLYLINE_DYNAM   draws polyline (dynam allocate reals)
C###  Routine: POLYMARKER       draws polymarker
C###  Routine: POLYMARKER_DYNAM draws polymarker (dynam allocate reals)
C###  Routine: PRINT_SCREEN_TO_FILE saves current to screen to a file
C###  Routine: QUIT_GRAPHICS    close any workstat and gks or phigs
C###  Routine: REFRESH_GRAPHICS buffered call to GXWAIT
C###  Routine: SET_COLOUR_LUT   set colour lookup table
C###  Routine: SET_COLOUR_LUT_RANGE *REMOVED* set col lut for a given range
C###  Routine: SETUP            perf setup operations for workstations
C###  Routine: TEXT             draws text
C###  Routine: VECTOR           draws a vector
C###  Routine: VISIB            change segment visibility


C#### Module: FE15
C###  Description:
C###    University of Washington routines.
C#### Module: FE19
C###  Description:
C###    Cellular models (DFN,N98,JRW,LR-II,HMT,DM,ATR)

C ***   Used to solve membrane and
C ***   mechanics models, as well as coupled electro-mechanics models,
C ***   for single cells.
C ***   Note: for Fading memory model see Bergel and Hunter, The
C ***   Mechanics of the Heart 1979
C ***
C ***   Distribution-Moment model
C ***   *************************
C ***   see Zahalak,G.I. and Ma, S.-P. Muscle activation and
C ***   contraction:
C ***   Constitutive relations based directly on cross-bridge kinetics.
C ***   J.Biomech.Eng Vol 112, pp 52-62, 1990.
C ***
C ***   Ma, S.-P. and Zahalak,G.I. A distribution-moment model of
C ***   energetics in skeletal muscle.
C ***   J.Biomech. Vol 24, #1 pp 21-35,1991.

C###  Routine: ADAMS             -> cmiss_archive2.f
C###  Routine: ADAMS_MOULTON     -> cmiss_archive2.f
C###  Routine: ATR               computes human atrial eqns
C###  Routine: ATR_CHANGE        set the stimulus current
C###  Routine: ATR_CURRENTS      compute ionic currents for ATR model
C###  Routine: ATR_EQUILIB_POTS  computes equilib. pot's for ATR model
C###  Routine: BETA_DM       (fn) compute beta
C###  Routine: BR                rhs routine for Beeler-Reuter cell
C###  Routine: BR_INIT_GRID      setup routine for BR cell model
C###  Routine: BR_ION        to cmiss_gridarchive.f
C###  Routine: BR_RATES      to cmiss_gridarchive.f
! DPN 19/02/98 CAREV unused ???
C###  Routine: CAREV             find ICa reversal potential
C###  Routine: CH_INA            sodium channel model
C###  Routine: CH_ICL            chloride channel model
C###  Routine: CONCUR            transfer conc.s & currents to arrays
C###  Routine: CUBIC             rhs routine for cubic cell model
C###  Routine: CUBIC_INIT_GRID   setup routine for cubic model
C###  Routine: COUPLED_SYSTEM    the main routine for coupled systems
C###  Routine: DEFINE_ATR        initialises variables for ATR mem model
C###  Routine: DEFINE_JRW        initialises variables for JRW mem model
C###  Routine: DEFINE_LR         initialises variables for Luo-Rudy
C###  Routine: DEFINE_NOBLE98    -> cmiss_archive2.f
C###  Routine: DEOXS1            define Oxsoft Heart control modes
C###  Routine: DEOXS2            define Oxsoft Heart model modes
C###  Routine: DEOXS3            define Oxsoft Heart parameters
C###  Routine: DESOL             driver routine for integrating odes
C###  Routine: DFN_CHANGE        check time to change parameter values
C###  Routine: DFN_CONCENTRATIONS  compute ion concentrations for HEART
C###  Routine: DFN_CURRENTS      compute ionic currents for HEART
C###  Routine: DFN_CURRENTS_PLOT compute max (ss) currents for plotting
C###  Routine: DFN_RATES         compute the rate coeffs for HEART
C###  Routine: DFN_RATES_SUB     calc's alpha & beta rate constants
C###  Routine: DFN_RATES_PLOT    calc's alpha/beta for plotting
C###  Routine: DIFRANCESCO_NOBLE solve diFrancesco-Noble equations
C###  Routine: DISTRIBUTION_MOMENT set-up the DM model
C###  Routine: DM_INIT       (block data) distribution moment variables
C###  Routine: DN                compute diFrancesco-Noble eqns
C###  Routine: DN_ION        to cmiss_gridarchive.f
C###  Routine: FADING_MEMORY fading memory model
C###  Routine: FCN_DM        compute RHS for DM model
C###  Routine: FHN           rhs routine for FitzFugh-Nagumo cell model
C###  Routine: FHN_INIT_GRID setup routine for FHN model
C###  Routine: FHN_ION       to cmiss_gridarchive.f
C###  Routine: FMINIT        initialise fading memory parameters
C###  Routine: FM_INTEGRAND1 (fn) integrand for fm model
C###  Routine: FM_INTEGRAND2 (fn) integrand for integral of velocity
C###  Routine: FM_SOLVE      solves the fading memory model
C###  Routine: FN_CA         (fn) Ca twitch function
C###  Routine: FN_G          (fn) G=T-Load, where T is calc.d from L
C###  Routine: FN_Q          (fn) updated Q for current timestep dt
C###  Routine: FN_TN         (fn) troponin kinetics
C###  Routine: FN_TM         (fn) tropomyosin kinetics
C###  Routine: FN_TO         (fn) isometric tension
C###  Routine: FN_TO_SS      (fn) steady state isometric tension
C###  Routine: FN_ZSS        (fn) steady state z
C###  Routine: GATES             compute time constants and act.n vars
C###  Routine: GG            (fn) calculation of G(k,z,p,q)
C###  Routine: HH                compute RHS of Hodgkin-Huxley eqns
C###  Routine: HH_CELL           compute RHS of Hodgkin-Huxley eqns
C###  Routine: HH_INIT_CELL      initialise HH model
C###  Routine: HMT_CELL          compute RHS of HMT eqns
C###  Routine: HMT_SPM_CELL      compute RHS of HMT eqns coupled to the Colorado calcium model
C###  Routine: HODGKIN_HUXLEY    solve Hodgkin-Huxley eqns
C###  Routine: INQUIRE_CELL_VARIABLE find index of cell var in YQS etc
C###  Routine: INTERP            interpolate y variables
C###  Routine: JRW               compute Jafri-Rice-Winslow eqns
C###  Routine: JRW_CHANGE        switches the stimulus current on/off
C###  Routine: JRW_CURRENTS      evaluate currents for JRW model
C###  Routine: JRW_CURRENTS_PLOT calc SS currents for plotting
C###  Routine: JRW_EQUILIB_POTS  evaluate equilibrium potentials
C###  Routine: JRW_RATES         compute rate coefficients
C###  Routine: JRW_RATES_PLOT    calc alpha/beta for plotting
C###  Routine: JRWP              compute JRW (Princeton) eqns
C###  Routine: LIOXSPARA         list OXSOFT parameters
C###  Routine: LR                -> cmiss_archive2.f
C###  Routine: LR_CURRENTS       -> cmiss_archive2.f
C###  Routine: LR_CURRENTS_PLOT  compute currents for plotting
C###  Routine: LR_ION            to cmiss_gridarchive.f
C###  Routine: LR_RATES          to cmiss_gridarchive.f
C###  Routine: LR_CELL           compute RHS of LR model
C###  Routine: LR_CURRENTS_CELL  compute LR currents
C###  Routine: LR_INIT_CELL      initialise LR model
C###  Routine: LR_RATES_CELL     compute LR rates
! **  DPN 20/03/98 - resolve conflict with name LR_RATES
C###  Routine: L_R_RATES         -> cmiss_archive2.f
C###  Routine: LR_RATES_PLOT     calc's alpha/beta for plotting
C###  Routine: LUO_RUDY          solve Luo-Rudy eqns
C###  Routine: NOBLE98_CELL      compute RHS of Noble '98 equations
C###  Routine: NOBLE98_HMT_CELL  compute RHS of Noble '98 - HMT
C###  Routine: NOBLE98           -> cmiss_archive2.f
C###  Routine: NOBLE98_CHANGE    -> cmiss_archive2.f
C###  Routine: NOBLE98_CURRENTS  -> cmiss_archive2.f
C###  Routine: NOBLE98_RATES     -> cmiss_archive2.f
C###  Routine: OXSINI            initialize Oxsoft heart parameters
C###  Routine: OXSOLVE           solve set of odes from Oxsoft heart
C###  Routine: OXSPARAMS         define Oxsoft Heart parameters
C###  Routine: OXSPREP           set default params for preparations
C###  Routine: PHI1          (fn) compute phi1
C###  Routine: PHI2          (fn) compute phi2
C###  Routine: TP1   (fn) passive tension-length relation - fibre axis
C###  Routine: TP2   (fn) passive tension-length relation - sheet axis
C###  Routine: TP3   (fn) passive tension-length relation - sheet normal
C###  Routine: UPCG_COUP         updates the CG material array with coupling properties
C###  Routine: VCD               rhs routine for van Capelle-Durrer cell
C###  Routine: VCD_INIT_GRID     setup routine for VCD cell model
C###  Routine: VCD_ION           to cmiss_gridarchive.f

C???  Routine: OUTPUT        output


CC#### Comment: Cell Modelling Units
CC###  Description:
CC###    <html>
CC###    <p>
CC###    The following table describes the units which <b>must</b> be
CC###    used in <b>all</b> cellular modelling in CMISS!!
CC###    </p>
CC###    <TABLE BORDER="1">
CC###        <TR>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>Quantity</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>Symbol or equation</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>Currently</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>Proposal</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>Consistency</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>Typical value</B></P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD></TD>
CC###            <TD></TD>
CC###            <TD></TD>
CC###            <TD></TD>
CC###            <TD></TD>
CC###            <TD></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Length</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                l</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                cm or mm</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>mm</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Area</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                A</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                cm^2 or mm^2</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>mm^2</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Volume</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                V</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                cm^3 or mm^3</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>mm^3</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>=uL</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                1e-6mm^3</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Time</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                t</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                ms</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>ms</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                1ms</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Voltage</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                V</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                mV</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>mV</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                1mV</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Current</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD>uA</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>uA</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Conductivity</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>mS</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                =uA/mV</P>
CC###                </CENTER></TD>
CC###            <TD>.</TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Membrane conductance</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                g</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                uS/cm^2</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>mS.mm^-2</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                1e-4mS/mm^2</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Tissue conductivity</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                \sigma</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                mS.mm^-1</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>mS.mm^-1</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                0.1mS/mm</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Current density</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                I</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                uA/cm^2</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>uA.mm^-2</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                =mV.mS/mm^2</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                0.01uA/mm^2</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Volume current</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>uA.mm^-3</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Charge</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                q</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>nC</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                =uA.ms</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Capacitance</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                uF</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>uF</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                =mS.ms</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Specific capacitance</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                Cm</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                nF/mm^2 or F/cm^2</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>uF.mm^-2</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                0.01uF/mm^2 (=1uF/cm^2)</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>#ions</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>nmol</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Concentration</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                [x]</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>nmol.mm^-3</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>=mM</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                1mM=1nmol/mm^3</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Conc rate (for valence 1)</TD>
CC###            <TD>d[x]/dt=I.(Area/Vol)/F</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>mM.ms^-1</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Mass</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                m</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>g</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Force</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                F</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                kN or N</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>mN</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                =g.mm/ms</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Stress (pressure)</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                T or p</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                kPa</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>mN.mm^-2</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>=kPa</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Energy</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>pJ</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                =nC.mV</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Power</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>nW</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                =pJ/ms</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###        <TR>
CC###            <TD>Temperature</TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                <B>K</B></P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                .</P>
CC###                </CENTER></TD>
CC###            <TD><CENTER><P ALIGN="CENTER">
CC###                293K</P>
CC###                </CENTER></TD>
CC###        </TR>
CC###    </TABLE>
CC###    <p><b>Note:</b></p>
CC###    <p>Gas constant R = 8.314472e3 pJ.nmol^-1.K^-1 <BR>
CC###    Faraday's constant F = 9.6485341e4 nC.nmol^-1<BR>
CC###    R/F = 8.61734219e-2 mV.K^-1</p>
CC###    <p><b>Useful conversions:</b></p>
CC###    <p>
CC###    uA = mV . mS <br>
CC###    nC . ms^-1 = uA <br>
CC###    mM = mmol . L^-1 = nmol . mm^-3 <br>
CC###    uF = mS . ms
CC###    </p>
CC###    </html>

C#### Comment: Cell Modelling Units
C###  Description:
C###    <html>
C###    <p>
C###    The following tables describe the units which <b>should</b> be
C###    used in <b>all</b> cellular modelling in CMISS.
C###    </p>
C###    <TABLE BORDER="1">
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                <B>Base Unit</B></P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                <B>SI</B></P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                <B>CMISS</B></P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                <B>Multiplier</B></P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD></TD>
C###            <TD></TD>
C###            <TD></TD>
C###            <TD></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                Length</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m - meter</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                1x10^3</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                Mass</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                kg - kilogram</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                ng</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                1x10^12</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                Time</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                s - seconds</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                ms</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                1x10^3</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                Electric current</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                A - ampere</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                uA</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                1x10^6</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                Temperature</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                K - kelvin</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                K</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                1x10^0</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                Amount of substance</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mol - mole</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                nmol</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                1x10^9</P>
C###                </CENTER></TD>
C###        </TR>
C###    </TABLE>
C###    <BR></BR><BR></BR>
C###    <TABLE BORDER="1">
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                <B>Derived Unit</B></P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                <B>SI</B></P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                <B>SI Base Units</B></P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                <B>CMISS</B></P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                <B>CMISS Base Units</B></P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD></TD>
C###            <TD></TD>
C###            <TD></TD>
C###            <TD></TD>
C###            <TD></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                area</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^2</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                volume</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^3</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^3</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^3</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^3</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                voltage</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                V</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^2.kg.s^-3.A^-1</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mV</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^2.ng.ms^-3.uA^-1</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                conductivity</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                S</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^-2.kg^-1.s^3.A^2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mS</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^-2.ng^-1.ms^3.uA^2</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                membrane conductance</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                S.m^-2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^-4.kg^-1.s^3.A^2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mS.mm^-2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^-4.ng^-1.ms^3.uA^2</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                tissue conductivity</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                S.m^-1</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^-3.kg^-1.s^3.A^2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mS.mm^-1</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^-3.ng^-1.ms^3.uA^2</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                current density</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                A.m^-2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^-2.A</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                uA.mm^-2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^-2.uA</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                volume current</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                A.m^-3</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^-3.A</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                uA.mm^-3</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^-3.uA</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                charge</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                C</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                s.A</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                nC</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                ms.uA</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                capacitance</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                F</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^-2.kg^-1.s^4.A^2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                uF</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^-2.ng^-1.ms^4.uA^2</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                specific capacitance</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                F.m^-2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^-4.kg^-1.s^4.A^2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                uF.mm^-2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^-4.ng^-1.ms^4.uA^2</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                concentration</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                [x]</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mol.m^-3</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mM</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                nmol.mm^-3</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                conc. rate (for valence 1)</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                d[x]/dt</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mol.m^-3.s^-1</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mM.ms^-1</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                nmol.mm^-3.ms^-1</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                force</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                N</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m.kg.s^-2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                nN</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm.ng.ms^-2</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                stress (pressure)</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                Pa</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^-1.kg.s^-2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mPa</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^-1.ng.ms^-2</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                energy</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                J</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^2.kg.s^-2</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                pJ</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^2.ng.ms^-2</P>
C###                </CENTER></TD>
C###        </TR>
C###        <TR>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                power</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                W</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                m^2.kg.s^-3</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                nW</P>
C###                </CENTER></TD>
C###            <TD><CENTER><P ALIGN="CENTER">
C###                mm^2.ng.ms^-3</P>
C###                </CENTER></TD>
C###        </TR>
C###    </TABLE>
C###    <p><b>Note:</b></p>
C###    <p>Gas constant R = 8.314472e3 pJ.nmol^-1.K^-1 <BR>
C###    Faraday's constant F = 9.6485341e4 nC.nmol^-1<BR>
C###    R/F = 8.61734219e-2 mV.K^-1</p>
C###    <p><b>Useful conversions:</b></p>
C###    <p>
C###    uA = mV . mS <br>
C###    nC . ms^-1 = uA <br>
C###    mM = mmol . L^-1 = nmol . mm^-3 <br>
C###    uF = mS . ms
C###    </p>
C###    </html>


C#### Module: FE20
C###  Description:
C###    Main fem environment for setting up array sizes.
C###  Routine: COM_DOC3  handles documenting commands
C###  Routine: FEM       allocates finite element arrays
C###  Routine: FEM_DYNAM checks command string in FEM environment
C###  Routine: SYNTAX    executes valid commands

C#### Variable: IPIVOT(no)
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###

C#### Variable: ISEG(nosg)
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    0,1,2 if segment not yet created/created but not visible/
C###    created and visible.

C#### Variable: IWK1(5*no)
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###

C#### Variable: IWK2(8*no)
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###

C#### Variable: IWK3(5*no)
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    Integer work array, duplicate of IWK1 used in bidomain.

C#### Variable: na
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    na is the auxiliary parameters loop variable.

C#### Variable: nae
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nae is the element line segment loop variable.
      
C#### Variable: nb
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nb is the basis function loop variable, which defines the
C###    interpolation used over an element.

C#### Variable: nc
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    <pre>
C###    nc is the variable type loop variable, and could be thought
C###    of as the number of dependent variable types in the
C###    governing equation.
C###
C###    nc=1 represents your dependent variable and is sufficient for
C###    many problems with governing equations of the form:
C###    GK*u=f
C###
C###    If you are using boundary elements then the governing equations
C###    contains a matrix multiplied by du/dn (which we will call q):
C###    GK*u=GQ*q
C###    There are now 2 nc's, u and q. In the current implementation
C###    nc=2 refers to q
C###
C###    There is also the facility for nc=3 for equations involving
C###    time derivatives (du/dt = d):
C###    GD*d+GK*u=GQ*q
C###    Here there are now 3 nc's, d,u,q. If there was only d and u
C###    with no q then you have only nc=1,2 again.
C###
C###    The second order time derivative (d2u/dt2 = m) completes the
C###    set where now there are 4 nc's, m,d,u,q in the equation:
C###    GM*m+GD*d+GK*u=GQ*q
C###
C###    Care must be taken not to confuse these derivatives with
C###    basis function derivatives that are a property of the
C###    interpolation functions and not the governing equation.
C###    </pre>

C#### Variable: nd
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nd is the data points loop variable.

C#### Variable: ne
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    ne is the element loop variable.

C#### Variable: ne_adj
C###  Type: INTEGER      
C###  Set_up: FEM
C###  Description:
C###    ne_adj is the adjacent element loop variable. This is a
C###    variable local to each element.
      
C#### Variable: nf
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nf is the global face segments loop variable.

C#### Variable: ng
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    ng is the Gauss point loop variable

C#### Variable: nh
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nh is the dependent variable loop variable.

C#### Variable: ni
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    ni is the local Xi-coordinates loop variable.

C#### Variable: nj
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nj holds an independent variable identifier.  This variable may
C###    be a global reference coordinate, a fibre variable or a
C###    (miscellaneous) `field' variable.

C#### Variable: NJT
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    NJT is the number of global coordinates.

C#### Variable: nl
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nl is the global line segments loop variable.

C#### Variable: nlat
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:  
C###    nlat is the lattice grid point loop variable.
      
C#### Variable: nm
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nm is the material parameters loop variable.

C#### Variable: nn
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nn is the element node loop variable.

C#### Variable: no
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    no is the dofs loop variable.

C#### Variable: np
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    np is the global node loop variable.

C#### Variable: nq
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nq is the global grid point loop variable.

C#### Variable: nr
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nr is the region loop variable.

C#### Variable: ns
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    ns is the element dofs per var loop variable.  It identifies a
C###    parameter for the interpolation over an element.

C#### Variable: nt
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nt is the eigenvalue loop variable.

C#### Variable: nts
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nts is the time sample loop variable.

C#### Variable: nv
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nv identifies a "version" of a node.

C#### Variable: nw
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nw is the workstations or special element loop variable.

C#### Variable: NWP(np,2)
C###  Type: INTEGER
C###  Set_up: HANGING_NODE_DETECT
C###  Description:
C###    <HTML> <PRE>
C###    NWP is set up as described in HANGING_NODE_DETECT for problems
C###    containing hanging nodes.
C###    </PRE> </HTML>

C#### Variable: nx
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nx is the problem type loop variable.

C#### Variable: ny
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    <HTML> <PRE>
C###    ny is the mesh dof loop variable.
C###    The following convention is used for mapping arrays:
C###    ny1v = ny# for global variable number
C###    ny2v = ny# for global flux/reaction number
C###    ny1r = ny# for row of LHS (nc=1) matrix
C###    ny1c = ny# for col of LHS (nc=1) matrix
C###    ny2r = ny# for row of RHS (nc=2) matrix
C###    ny2c = ny# for col of RHS (nc=2) matrix
C###    </PRE> </HTML>

C#### Variable: nz
C###  Type: INTEGER
C###  Set_up: FEM
C###  Description:
C###    nz is the coefficients in the global stiffness matrix loop variable.


C#### Module: FE21
C###  Description:
C###    Routines called by miscellaneous commands.
C###  Routine: CHMESH_ARCLENGTH (fn) J.Crocombes change mesh fns
C###  Routine: CHMESH_DFDS   (fn) J.Crocombes change mesh fns
C###  Routine: CHMESH_DGDS   (fn) J.Crocombes change mesh fns
C###  Routine: CHMESH_FDERVI (fn) J.Crocombes change mesh fns
C###  Routine: CHMESH_FFUNC  (fn) J.Crocombes change mesh fns
C###  Routine: CHMESH_GDERIV (fn) J.Crocombes change mesh fns
C###  Routine: CHMESH_GFUNC  (fn) J.Crocombes change mesh fns
C###  Routine: ADDFIEL   "adds" two fields to produce a third field
C###  Routine: ADDSIG    "adds" two signals to produce a third signal
C###  Routine: APINVE    applies inverted transfer matrix
C###  Routine: APNOIS    applies noise to a signal file
C###  Routine: APREFE    applies a new reference to a set of signals
C###  Routine: APTRSF    applies transfer matrix
C###  Routine: CHBASE    changes basis functions
C###  Routine: CHCOOR    changes coordinate system
C###  Routine: CHDATA    changes data
C###  Routine: CHFOCU    changes focus
C###  Routine: CHLINE    changes lines
C###  Routine: CHMESH    changes mesh
C###  Routine: CHMESH_ALTER alters the mesh values
C###  Routine: CHMESH_CALC_ARCLENGTHS calculates cirumference of mesh
C###  Routine: CHMESH_FIND_ZETA calculates zeta
C###  Routine: CHMESH_NORMALISE normAlises mesh derivatives
C###  Routine: CHNODE    changes node
C###  Routine: CHNODS    changes nodes
C###  Routine: CHORD     defines individual chord
C###  Routine: CHPLIN    changes polylines
C###  Routine: CHREGP    changes epicardial inverse regularisation parameters
C###  Routine: CLOSE_FILES closes files
C###  Routine: COMBSIG   combines signal files
C###  Routine: COMPDAT   compares data files
C###  Routine: COMPGAUSSVAR comps GaussPt vals (YG) to (interp) YQS vals
C###  Routine: COMPGRIDVAR comps grid point values (YQS)
C###  Routine: COMPSIG   compares signal files
C###  Routine: COELEM    copies elements
C###  Routine: CONODE    copies nodes
C###  Routine: CONVSIG   convert signals
C###  Routine: CRHERM    creates Hermite lines with Bezier ctrl pts
C###  Routine: DATA3D    *REMOVED* 3D transformation for data points
C###  Routine: DATXY     calc width in x and y for given z in data
C###  Routine: FIT       fits with linear or nonlinear least-squares
C###  Routine: HANGING_NODE_DETECT detects hanging nodes and flags them in NWP
C###  Routine: ITERATE   handles iteration calls to fem (buffer)
C###  Routine: ITERATE_DYNAM  handles iteration calls to fem
C###  Routine: MAKESPLINE creates a 1D Spline
C###  Routine: MESHXY    calculates mesh width in x and y for given z
C###  Routine: MOMENTS   calculates 2nd moments of 2-D node & data set
C###  Routine: OPEN_FILES opens files
C###  Routine: PRSCRN    prints 2D graphics screen to a file
C###  Routine: READF     reads files
C###  Routine: REFINE    refines mesh
C###  Routine: REGRAP    refresh graphics
C###  Routine: SHAPE_CHORDS shapes chords
C###  Routine: SHAPEC    shapes chord
C###  Routine: SOLVE     solves a boundary value problem
C###  Routine: STEP      steps through time dependent nodal values
C###  Routine: TRACK     tracks current lines through the domain
C###  Routine: WRITEF    writes files

C#### Variable: ISL2BE(nl)
C###  Type: INTEGER
C###  Set_up: CRHERM
C###  Description:
C###    ISL2BE(nl) is segment number of Bezier tangent line 1 on nl.

C#### Variable: ISL3BE(nl)
C###  Type: INTEGER
C###  Set_up: CRHERM
C###  Description:
C###    ISL3BE(nl) is segment number of Bezier tangent line 2 on nl.

C#### Variable: ISN2BE(nl)
C###  Type: INTEGER
C###  Set_up: CRHERM
C###  Description:
C###    ISN2BE(nl) is segment number of Bezier control point 1 on nl.

C#### Variable: ISN3BE(nl)
C###  Type: INTEGER
C###  Set_up: CRHERM
C###  Description:
C###    ISN3BE(nl) is segment number of Bezier control point 2 on nl.

C#### Variable: NIYLIST(niy)
C###  Type: INTEGER
C###  Set_up: OPEN_FILES
C###  Description:
C###    NIYLIST is used to store the niy values which need to
C###    be written out to a history file using IOHIST.
C###    NIYLIST(0) contains the number of niy's stored in the list.

C#### Variable: NLCHOR(0:10,nr)
C###  Type: INTEGER
C###  Set_up: CHORD
C###  Description:
C###    NLCHOR(0:10,nr) was used in the sail work, and may be
C###    redundant now.

C#### Variable: NLINE1(nol1)
C###  Type: INTEGER
C###  Set_up: CHNODS
C###  Description:
C###    NLINE1(nol1) stores line numbers of lines which have
C###    NPL(3,1,nl)=np.

C#### Variable: NLINE2(nol2)
C###  Type: INTEGER
C###  Set_up: CHNODS
C###  Description:
C###    NLINE2(nol2) stores line numbers of lines which have
C###    NPL(2,1,nl)=np.

C#### Variable: NPLIST1(np)
C###  Type: INTEGER
C###  Set_up: .
C###  Description:
C###    NPLIST1(np) is a temporary node list
C###    is to be applied.

C#### Variable: NPLIST2(np)
C###  Type: INTEGER
C###  Set_up: .
C###  Description:
C###    NPLIST2(np) is a temporary node list

C#### Variable: DXI(2,nj)
C###  Type: REAL*8
C###  Set_up: CRHERM
C###  Description:
C###    DXI(1,nj) is derivative wrt Xi at 1st node.
C###    DXI(2,nj) is derivative wrt Xi at 2nd node.

C#### Variable: PROPQ
C###    <HTML>
C###  Set_up: UPGRID,CALC_FV_GRID_COEF,CACL_FV_GRID_SECFLUX 
C###  Description:
C###    <HTML>
C###    1. For grid based collocation solutions - grid point propagation 
C###       information:
C###    <PRE>PROPQ is conductivity tensors for activation at grid point nq.
C###    </PRE><PRE>
C###    PROPQ(i,j,1,1,nq,nx)   is C(I)ij at nq;
C###    PROPQ(i,j,1,2,nq,nx)   is C(E)ij at nq;
C###    PROPQ(i,j,k+1,1,nq,nx) is C(I)ij,k at nq; and
C###    PROPQ(i,j,k+1,2,nq,nx) is C(E)ij,k at nq,
C###      where C(I) and C(E) are intracellular and extracellular
C###      conductivities, respectively.
C###    </PRE>
C###    <PRE>2. For grid based finite volume solutions - grid point properties
C###            information:
C###    </PRE><PRE>
C###    U   : U geometric volume face vector (i.e. vectors in the face)
C###    V   : V geometric volume face vector (i.e. vectors in the face)
C###    INU : product of intracellular conductivity tensor with face normal
C###          and U face vector and face area
C###    INV : product of intracellular conductivity tensor with face normal
C###          and V face vector and face area
C###    ENU : product of extracellular conductivity tensor with face normal
C###          and U face vector and face area
C###    ENV : product of extracellular conductivity tensor with face normal
C###          and V face vector and face area
C###    DV  : derivative approximations of transmembrane potential
C###    DE  : derivative approximations of extracellular potential
C###    LP  : finite volume dimensions in +dir
C###    LN  : finite volume dimensions in -dir
C###    SF*1: secondary flux contribution to first equation rhs 
C###          contributions are stored for the three possible positive directions
C###          (SFP*1) and the three possible negative directions (SFN*1).
C###    SF*2: secondary flux contribution to second equation rhs
C###          contributions are stored for the three possible positive directions
C###          (SFP*1) and the three possible negative directions (SFN*1).
C###    </PRE><PRE>
C###    PROPQ(i,j,k,l):
C###    </PRE><PRE>
C###    l=1
C###    k=1           k=2          k=3    k=4
C###    +1            +2           +3             - face
C###    x   y   z     x   y   z    x y z
C###    U,  U,  U     U,  U,  U    U,U,U  -,-,-
C###    V,  V,  V     V,  V,  V    V,V,V  -,-,-
C###    INU,INU,INU   INV,INV,INV  -,-,-  -,-,-
C###    +1  +2  +3    +1  +2  +3                  - face
C###    </PRE><PRE>
C###    l=2
C###    k=1           k=2          k=3      k=4
C###    x   y   z     +1  +2  +3                  - direction/face
C###    DV, DV, DV    LP, LP, LP   SFP11,SFP21,SFP31  SFN11,SFN21,SFN31
C###    DE, DE, DE    LN, LN, LN   SFP12,SFP22,SFP32  SFN12,SFN22,SFN32
C###    ENU,ENU,ENU   ENV,ENV,ENV  -,  -,-  -,-,-
C###    +1  +2  +3    +1  +2  +3                  - face
C###    </PRE>
C###    </HTML>

C#### Variable: RE1(nh,ns)
C###  Type: REAL*8
C###  Set_up: SOLVE
C###  Description:
C###    RE1(nh,ns) is reaction from element variables.

C#### Variable: RE2(nh,ns)
C###  Type: REAL*8
C###  Set_up: SOLVE
C###  Description:
C###    RE2(nh,ns) is reaction from element variables.

C#### Variable: TRANS(3,4)
C###  Type: REAL*8
C###  Set_up: CHDATA
C###  Description:
C###    Transformation matrix


C#### Variable: XB(2,nj,nl)
C###  Type: REAL*8
C###  Set_up: CRHERM
C###  Description:
C###    XB(1,nj,nl) is first Bezier slope control point for line nl.
C###    XB(2,nj,nl) is second Bezier slope control point for line nl.

C#### Variable: XBEZ(i)
C###  Type: REAL*8
C###  Set_up: CRHERM
C###  Description:
C###    XBEZ(i), i=1 and 4, are x coordinates of nodes.
C###    XBEZ(i), i=2 and 3, are x coordinates of slope control points.

C#### Variable: YBEZ(i)
C###  Type: REAL*8
C###  Set_up: CRHERM
C###  Description:
C###    YBEZ(i), i=1 and 4, are y coordinates of nodes.
C###    YBEZ(i), i=2 and 3, are y coordinates of slope control points.


C#### Module: FE22
C###  Description:
C###    Routines called by 'cancel','hide' & 'show' commands

C###  Routine: CAALIG   cancel alignment
C###  Routine: CAAXES   cancel axes
C###  Routine: CABASE   cancel bases
C###  Routine: CACLOC   cancel clock
C###  Routine: CACONT   cancel contours
C###  Routine: CACROS   cancel cross-section
C###  Routine: CADATA   cancel data
C###  Routine: CAELEM   cancel elements
C###  Routine: CAFACE   cancel faces
C###  Routine: CAFIBR   cancel fibres
C###  Routine: CAFIEL   cancel field
C###  Routine: CAGAUS   cancel Gauss points
C###  Routine: CAGRAD   cancel gradient
C###  Routine: CAGRID   cancel grid
C###  Routine: CAHIST   cancel history plots
C###  Routine: CAINCR   cancel increment vectors
C###  Routine: CAISOC   cancel isochrones
C###  Routine: CALINE   cancel lines
C###  Routine: CAMAP    cancel map
C###  Routine: CAMATE   cancel materials
C###  Routine: CANODE   cancel nodes
C###  Routine: CAOBJE   cancel object
C###  Routine: CAPLIN   cancel polyline
C###  Routine: CAPLOT   cancel plot
C###  Routine: CAPMAR   *REMOVED* cancel polymarker
C###  Routine: CAPROF   cancel profile
C###  Routine: CAREAC   cancel reaction vectors
C###  Routine: CARESI   cancel residual vectors
C###  Routine: CARULE   cancel rule
C###  Routine: CASECT   cancel section plots
C###  Routine: CASHEE   cancel sheets
C###  Routine: CASTRA   cancel strain plots
C###  Routine: CASTRE   cancel stress plots
C###  Routine: CASURF   cancel surface plots
C###  Routine: CAVISI   *REMOVED* calculate visibility of segments
C###  Routine: HIALIG   hide alignment
C###  Routine: HIAXES   hide axes
C###  Routine: HICLOC   hide clock
C###  Routine: HICONT   hide contours
C###  Routine: HICROS   hide cross-section
C###  Routine: HIDATA   hide data
C###  Routine: HIDIPO   hide dipole
C###  Routine: HIELEM   hide elements
C###  Routine: HIFACE   hide faces
C###  Routine: HIFIBR   hide fibres
C###  Routine: HIFIEL   hide field
C###  Routine: HIGAUS   hide Gauss points
C###  Routine: HIGRAD   hide gradient
C###  Routine: HIGRID   hide grid
C###  Routine: HIHIST   hide history plots
C###  Routine: HIINCR   hide increment vectors
C###  Routine: HIISOC   hide isochrones
C###  Routine: HILINE   hide lines
C###  Routine: HIMAP    hide map
C###  Routine: HIMATE   hide materials
C###  Routine: HINODE   hide nodes
C###  Routine: HIOBJE   hide objects
C###  Routine: HIPLIN   hide polyline
C###  Routine: HIPMAR   *REMOVED* hide polymarker
C###  Routine: HIPROF   hide profile
C###  Routine: HIREAC   hide reaction vectors
C###  Routine: HIRESI   hide residual vectors
C###  Routine: HIRULE   hide rule
C###  Routine: HISCAL   hide scale
C###  Routine: HISECT   hide section plots
C###  Routine: HISHEE   hide sheets
C###  Routine: HISTRA   hide strain plots
C###  Routine: HISTRE   hide stress plots
C###  Routine: HISTRM   hide streamlines
C###  Routine: HISURF   hide surfaces
C###  Routine: HIVELO   hide velocity vectors
C###  Routine: PICONT   *** ARCHIVED  ***
C###  Routine: PIDATA   *** ARCHIVED ***
C###  Routine: PIELEC   *** ARCHIVED ***
C###  Routine: PIELEM   *** ARCHIVED ***
C###  Routine: PIFACE   *** ARCHIVED ***
C###  Routine: PIFIBR   *** ARCHIVED ***
C###  Routine: PILINE   *** ARCHIVED ***
C###  Routine: PIMATE   *** ARCHIVED ***
C###  Routine: PINODE   *** ARCHIVED ***
C###  Routine: SHALIG   show alignment
C###  Routine: SHAXES   show axes
C###  Routine: SHCLOC   show clock
C###  Routine: SHCONT   show contours
C###  Routine: SHCROS   show cross-section
C###  Routine: SHDATA   show data
C###  Routine: SHDIPO   show dipole
C###  Routine: SHELEM   show elements
C###  Routine: SHFACE   show faces
C###  Routine: SHFIBR   show fibres
C###  Routine: SHFIEL   show field
C###  Routine: SHGAUS   show Gauss points
C###  Routine: SHGRAD   show gradient
C###  Routine: SHGRID   show grid
C###  Routine: SHHIST   show history plots
C###  Routine: SHINCR   show increment vectors
C###  Routine: SHISOC   show isochrones
C###  Routine: SHLINE   show lines
C###  Routine: SHMAP    show map
C###  Routine: SHMATE   show materials
C###  Routine: SHNODE   show nodes
C###  Routine: SHOBJE   show objects
C###  Routine: SHPLIN   show polyline
C###  Routine: SHPMAR   *** ARCHIVED ***
C###  Routine: SHPROF   show profile
C###  Routine: SHREAC   show reaction vectors
C###  Routine: SHRESI   show residual vectors
C###  Routine: SHRULE   show rule
C###  Routine: SHSCAL   show scale
C###  Routine: SHSECT   show section plots
C###  Routine: SHSHEE   show sheets
C###  Routine: SHSTRA   show strain plots
C###  Routine: SHSTRE   show stress plots
C###  Routine: SHSTRM   show streamlines
C###  Routine: SHSURF   show surfaces
C###  Routine: SHVELO   show velocity vectors

C#### Module: FE23
C###  Description:
C###    Update Routines

C###  Routine: UPAERO    update aerofoil
C###  Routine: UPCOUP    update coupling
C###  Routine: UPDATA    update data
C###  Routine: UPDELA    update delaunay
C###  Routine: UPELEM    update elements
C###  Routine: UPFG      processes "fem upd solu/geom/fiel/mate" commands
C###  Routine: UPFG_OPERATE    calls the operations for the update command
C###  Routine: UPFGADD   operation FG=F+G
C###  Routine: UPFGDIVIDE      operation FG=F/G
C###  Routine: UPFGGAUSS transfers XG/YG/CG to/from F/FG/G
C###  Routine: UPFGMULTIPLY    operation FG=F*G
C###  Routine: UPFGNODES transfers XP/YP/CP to/from F/FG/G
C###  Routine: UPFGSUBSTITUTE  operation FG=G
C###  Routine: UPFGSUBTRACT    operation FG=F-G
C###  Routine: UPFIEL    update field
C###  Routine: UPFIEL_EIK update field for eikonal equation
C###  Routine: UPFLUX    update boundary condition fluxes
C###  Routine: UPGAUS    update Gauss point array
C###  Routine: UPGAUS_EIK updates Gauss point array for eikonal equation
C###  Routine: UPGEOM    update geometry
C###  Routine: UPGREL    update element groups
C###  Routine: UPGRID    update material parameters
C###  Routine: UPGRID_CONN update grid point connectivity for cleavage planes
C###  Routine: UPGRNO    update node groups
C###  Routine: UPGROW    update growth
C###  Routine: UPINIT    update initial conditions
C###  Routine: UPINVE    update inverse solution
C#### Routine: CALC_LNORM        (fn) computes solution norm
C#### Routine: CALC_RNORM        (fn) computes residual norm
C#### Routine: CALC_INVERSE_SOLN computes new inverse soln
C###  Routine: UPLINE    update scale_factors so they are the same for common lines
C###  Routine: SE2INPUT  enters scalefactors into SE2 for two nodes sharing line
C###  Routine: SE2CHANGE matches and changes the scale-factors in SE2
C###  Routine: SE2OUTPUT updates scale-factors with values from SE2 into SE
C###  Routine: UPMATE    update material
C###  Routine: UPMESH    update mesh
C###  Routine: UPNODE    update nodal coordinates
C###  Routine: UPOPTI    update fe vars from optimisation params
C###  Routine: UPORDR    dynamic mem. wrapper for UPORDR_DYNAM
C###  Routine: UPORDR_DYNAM update ordering of 1D mesh
C###  Routine: UPPHI     update PHI matrix
C###  Routine: UPPRES    update auxiliary pressure parameters
C###  Routine: UPRESI    update residual vector
C###  Routine: UPSCAL    update scale factors from line lengths
C###  Routine: UPSIGN    update a signal file
C###  Routine: UPSLAVE   updates the slave nodal values and their derivatives
C###  Routine: UPSOBE    update Sobolev weights and scaling factors
C###  Routine: UPSOLU    updates solution (YP) array
C###  Routine: UPSOUR    update GD from current dipole sources
C###  Routine: UPVERT    update vertex position, for Voronoi vertices
C###  Routine: UPVIEW    update view
C###  Routine: UPXI      update xi coordinates

C#### Variable: XQ(nj,nq)
C###  Type: REAL*8
C###  Set_up: UPGRID
C###  Description:
C###    XQ(nj,nq) is rectangular cartesian coordinates of grid point nq.

C#### Variable: Z_CONT_LIST(nd,i,j)
C###  Type: INTEGER
C###  Set_up: UPDATA
C###  Description:
C###    <HTML>
C###    Z_CONT_LIST contains contact point information (j=1,4) for either
C###    a slave surface (i=1) or master surface (i=2) and associated with
C###    the data point nd.
C###    <PRE>
C###    j=1     element #
C###    j=2     local face #
C###    j=3     contact indicator (1=Contact,0=No Contact)
C###    j=4     contact type (1=Frictionless,2=Tied,3=Friction)
C###    </PRE>
C###    </HTML>

C#### Variable: Z_CONT(nd,i,j)
C###  Type: REAL*8
C###  Set_up: UPDATA
C###  Description:
C###    <HTML>
C###    Z_CONT contains contact point information (j=1..22) for either a
C###    slave surface (i=1) or master surface (i=2) and associated with
C###    the data point nd.
C###    <PRE>
C###    j=1..3      XI locations [Xi1,Xi2,Xi3]
C###    j=4..6      Target body unit normal [x,y,z]
C###    j=7..9      Target body unit tangent(1) [x,y,z]
C###    j=10..12    Target body unit tangent(2)[x,y,z]
C###    j=13        signed normal gap
C###    j=14        signed tangential(1) gap
C###    j=15        signed tangential(2) gap
C###    j=16        jacobian at Xi location
C###    j=17        Current estimate of normal contact force F_hat
C###    j=18        New estimate of F_hat at load step convergence
C###    j=19        Current estimate of tangential contact force F_hatT1
C###    j=20        New estimate of F_hatT1 at load step convergence 
C###    j=21        Current estimate of tangential contact force F_hatT2
C###    j=22        New estimate of F_hatT2 at load step convergence 
C###    </PRE>
C###    </HTML>

C#### Module: FE24
C###  Description:
C###    Routines called by the 'define' command.

C###  Routine: DEACTI   define active muscle fibre properties
C###  Routine: DEANAL   define analytic solution
C###  Routine: DEBASE   define basis functions
C###  Routine: DEBOUN   define boundary conditions
C###  Routine: DECELL   define cell parameters
C###  Routine: DECHOR   define chords on sail
C###  Routine: DECLOC   define clock
C###  Routine: DECONS   *REMOVED* define constant
C###  Routine: DECOOR   define coordinates
C###  Routine: DECORN   define corners (for BE problems)
C###  Routine: DECOUP   define coupling between dep. vars in regions
C###  Routine: DECUST   define customisation parameters
C###  Routine: DEDATA   define data
C###  Routine: DEELEM   define elements
C###  Routine: DEEQUA   define equations
C###  Routine: DEEXPO   define export parameters
C###  Routine: DEFACE   define faces
C###  Routine: DEFIBR   define fibres
C###  Routine: DEFIDU   define fiducial markers
C###  Routine: DEFIEL   define field
C###  Routine: DEFILE   define i/p file for another program
C###  Routine: DEFIT    define fitting parameters
C###  Routine: DEGAUS   define Gauss points
C###  Routine: DEGRID   define finite difference grid
C###  Routine: DEGROW   define growth
C###  Routine: DEHEAD   define heading
C###  Routine: DEIMPO   define import parameters
C###  Routine: DEINCR   define increment
C###  Routine: DEINIT   define initial & boundary conditions
C###  Routine: DEINVE   define inverse matrix regularisation options
C###  Routine: DEISOC   define isochrones
C###  Routine: DEITER   define iteration
C###  Routine: DELEAD   define electrocardiographic leads
C###  Routine: DELINE   define lines
C###  Routine: DEMAP    define map window
C###  Routine: DEMATE   define materials
C###  Routine: DEMATR   *REMOVED* define matrix
C###  Routine: DEMESH   define mesh
C###  Routine: DEMOTI   define motion parameters
C###  Routine: DENODE   define node
C###  Routine: DENODS   define nodes
C###  Routine: DENOIS   define noise paramters
C###  Routine: DENORM   define normal reversals
C###  Routine: DEOPTI   define optimisation
C###  Routine: DEPARA   define array dimensions
C###  Routine: DEPOTE   *** archived 9-NOV-98 ***
C###  Routine: DEREFE   define reference
C###  Routine: DEREFI   define refinement
C###  Routine: DEREGI   define regions
C###  Routine: DESAIL   *REMOVED* define sail
C###  Routine: DESING   define location of singularities (BE)
C###  Routine: DESOLV   define solution
C###  Routine: DESOUR   define sources
C###  Routine: DESTRA   define strains
C###  Routine: DESURF   define calc of surface and tag subdivisions
C###  Routine: DETEXT   *** ARCHIVE ***
C###  Routine: DETIME   define time variables
C###  Routine: DETRSF   define transfer matrix
C###  Routine: DEVECT   *** ARCHIVE ***
C###  Routine: DEWIND   define window
C###  Routine: DEXI     define Xi coordinates
C###  Routine: DEXI_CLOSEST define Xi coordinates closest
C###  Routine: DEXI_CLOSEST_FACE define Xi coordinates closest_face
C###  Routine: DEXI_CLOSEST_NODE define Xi coordinates for nodes
C###  Routine: DEXI_POINT define Xi coordinates at a gobal coord
C###  Routine: DEXI_ORTHOG define Xi coordinates orthogonal
C###  Routine: DEXI_1D  define Xi coordinates 1D
C###  Routine: DEXI_NONLIN define Xi coordinates nonlinear
C###  Routine: DEXI_LINEAR define Xi coordinates linear
C###  Routine: DEXI_EXISTING define Xi previously calculated


C#### Variable: ASSEMBLE_GLOBAL(nr,nx)
C###  Type: INTEGER
C###  Set_up: DEEQUA
C###  Description:
C###    ASSEMBLE_GLOBAL indicates whether the global stiffness matrices
C###    have been assembled or not. It is reset to .FALSE. if deequa or
C###    demate is called.

C#### Variable: FD(nd)
C###  Type: INTEGER
C###  Set_up: DEXI_CLOSEST_FACE
C###  Description:
C###  FD(nd) is the local face number associated with data point nd.

C#### Variable: KTYP6
C###  Type: INTEGER
C###  Set_up: DEEQUA
C###  Description:
C###    KTYP6 is the number of boundary integral domains.

C#### Variable: MXI(2,ne)
C###  Type: INTEGER
C###  Set_up: DEMAP
C###  Description:
C###    MXI(2,ne) are the bottom left coords for Xi map projection

C#### Variable: NLLIST(0:nl)
C###  Type: INTEGER
C###  Set_up: DEELEM
C###  Description:
C###    NLLIST(0:nl) is a temporary array for list of lines.

C#### Variable: NPINTER(3,0:20)
C###  Type: INTEGER
C###  Set_up: DENODS
C###  Description:
C###

C#### Variable: NUNK(nk,nj,np)
C###  Type: INTEGER
C###  Set_up: CALC_NUNK
C###  Description:
C###    NUNK returns a nu description of a derivative nk for a
C###    particular nj and np
C###  See-Also: nu

C#### Variable: EDD(nd)
C###  Type: REAL*8
C###  Set_up: DEXI
C###  Description:
C###    EDD(nd) is the error at data point nd.

C#### Variable: SE(ns,nbf,ne)
C###  Type: REAL*8
C###  Set_up: DLSE
C###  Description:
C###    SE(ns,nbf,ne) is the sacling factor for element dof ns of basis
C###    nbf on element ne.

C#### Variable: SF(ns,nbf)
C###  Type: REAL*8
C###  Set_up: CALC_FACE_INFORMATION_DEP,CALC_FACE_INFORMATION_IND
C###  Description:
C###    SF(ns,nbf) is the sacling factor for face dof ns of basis nbf.

C#### Variable: SQ(nd)
C###  Type: REAL*8
C###  Set_up: DEXI
C###  Description:
C###    SQ(nd) is the square of the distance from mesh to data point nd.

C#### Variable: ZC(nj,ne)
C###  Type: REAL*8
C###  Set_up: DEELEM
C###  Description:
C###    ZC(nj,ne) are the coordinates of the center of element ne

C#### Variable: ZD2(nj,nd)
C###  Type: REAL*8
C###  Set_up: DENODS
C###  Description:
C###    ZD2(nj,nd) is a temporary duplicate of ZD(nj,nd) used in the
C###    calculation of XP in DENODS


C#### Module: FE25
C###  Description:
C###    Routines for creating graphical segments.

C###  Routine: SGALIG   create segment for alignment
C###  Routine: SGAXES   create segment for axes
C###  Routine: SGBASE   create segment for basis functions
C###  Routine: SGCLOC   create segment for clock
C###  Routine: SGCONT   create segment for contours
C###  Routine: SGCROS   create segment for cross-section
C###  Routine: SGDATA   create segment for data
C###  Routine: SGDIPO   create segment for dipole
C###  Routine: SGELEM   create segment for element numbers
C###  Routine: SGELEM   create segment for element errors
C###  Routine: SGFACE   create segment for face numbers
C###  Routine: SGFIBR   create segment for fibres
C###  Routine: SGFIEL   create segment for field
C###  Routine: SGFIPR   create segment for fibre profiles
C###  Routine: SGGAUS   create segment for Gauss points
C###  Routine: SGGRAD   create segment for gradient vectors
C###  Routine: SGGRID   create segment for grid point
C###  Routine: SGGRID_NORM create segment for grid normals
C###  Routine: SGHIST   create segment for history plots
C###  Routine: SGINCR   create segment for increment vectors
C###  Routine: SGISOC   create segment for isochrones
C###  Routine: SGLEAD   create segment for electrocardiographic leads
C###  Routine: SGLINE   create segment for lines
C###  Routine: SGMAP    create segment for map
C###  Routine: SGMATE   create segment for materials
C###  Routine: SGNODE   create segment for node numbers
C###  Routine: SGOBJE   create segment for objects
C###  Routine: SGPLIN   create segment for polyline
C###  Routine: SGPLOT   create segment for plot
C###  Routine: SGPLOTXY create segment for plotxy
C###  Routine: SGPLOTXYLABEL create axes-scales for SGPLOTXY
C###  Routine: SGPMAR   create segment for polymarker
C###  Routine: SGPROF   create segment for profile plots
C###  Routine: SGREAC   create segment for reaction vectors
C###  Routine: SGRESI   create segment for residual vectors
C###  Routine: SGRULE   create segment for ruled grid
C###  Routine: SGSCAL   create segment for scale
C###  Routine: SGSECT   create segment for section plots
C###  Routine: SGSHEE   create segment for sheets
C###  Routine: SGSIGN   create segment for signal
C###  Routine: SGSTRA   create segment for strain plots
C###  Routine: SGSTRE   create segment for stress plots
C###  Routine: SGSTRM   create segment for streamline
C###  Routine: SGSURF   create segment for surfaces
C###  Routine: SGTEXT   create segment for text
C###  Routine: SGTRAC   create segment for trace plots

C#### Variable: ISALIG(iw)
C###  Type: INTEGER
C###  Set_up: SGALIG
C###  Description:
C###    ISALIG(iw) is the segment number of alignment grid.

C#### Variable: ISAXES(iw)
C###  Type: INTEGER
C###  Set_up: SGAXES
C###  Description:
C###    ISAXES(iw) is the segment number of axes.

C#### Variable: ISBASE(nb)
C###  Type: INTEGER
C###  Set_up: SGBASE
C###  Description:
C###    ISBASE(nb) is segment number of basis function type nb.

C#### Variable: ISCLOC(iw)
C###  Type: INTEGER
C###  Set_up: SGCLOC
C###  Description:
C###    ISCLOC(iw) is the segment number of clock.

C#### Variable: ISCONO(nh,ne)
C###  Type: INTEGER
C###  Set_up: SGCONT
C###  Description:
C###    ISCONO(nh,ne) is the segment number of contour numbers in
C###    element ne.

C#### Variable: ISCONT(nh,ne,nocont)
C###  Type: INTEGER
C###  Set_up: SGCONT
C###  Description:
C###    ISCONT(nh,ne,nocont) is the segment number of contour nocont
C###    of variable nh in ne.

C#### Variable: ISCROS(iw,nocros)
C###  Type: INTEGER
C###  Set_up: SGCROS
C###  Description:
C###    ISCROS(iw,nocros) is the segment number of axes.

C#### Variable: ISDATA(iw)
C###  Type: INTEGER
C###  Set_up: SGDATA
C###  Description:
C###    ISDATA(iw) is the segment number of data points at segment
C###    nodata.

C#### Variable: ISDATR(iw,ne)
C###  Type: INTEGER
C###  Set_up: SGDATA
C###  Description:
C###    ISDATR(iw,ne) is the segment number of data point trace in
C###    element ne.

C#### Variable: ISDANO(iw,ne)
C###  Type: INTEGER
C###  Set_up: SGDATA
C###  Description:
C###    ISDANO(iw,ne) is the segment number of data point numbers.

C#### Variable: ISDAPR(iw,ne)
C###  Type: INTEGER
C###  Set_up: SGDATA
C###  Description:
C###    ISDAPR(iw,ne) is the segment number of data point projections.

C#### Variable: ISDIPO(iw,n,nr)
C###  Type: INTEGER
C###  Set_up: SGDIPO
C###  Description:
C###    ISDIPO(iw,n,nr) is the dipole segment on window iw for dipole
C###    n in region nr.

C#### Variable: ISDIPA(iw,n,nr)
C###  Type: INTEGER
C###  Set_up: SGDIPO
C###  Description:
C###    ISDIPA(iw,n,nr) is the dipole path segment on window iw for
C###    dipole n in region nr.

C#### Variable: ISELNO(iw,ne)
C###  Type: INTEGER
C###  Set_up: SGELEM
C###  Description:
C###    ISELNO(iw,ne) is the segment number of element numbers.

C#### Variable: ISERR(iw,ne)
C###  Type: INTEGER
C###  Set_up: SGERR
C###  Description:
C###    ISERR(iw,ne) is the segment number of element errors.

C#### Variable: ISFACE(iw,nf)
C###  Type: INTEGER
C###  Set_up: SGFACE
C###  Description:
C###    ISFACE(iw,nf) is the segment number of face nf.

C#### Variable: ISFANO(iw,nf)
C###  Type: INTEGER
C###  Set_up: SGFACE
C###  Description:
C###    ISFANO(iw,nf) is the segment number of face numbers.

C#### Variable: ISFIBR(iw,ne,nofibr)
C###  Type: INTEGER
C###  Set_up: SGFIBR
C###  Description:
C###    ISFIBR(iw,ne,nofibr) are segment numbers of fibres in element
C###    ne at set nofibr.

C#### Variable: ISFIEL(iw,ne)
C###  Type: INTEGER
C###  Set_up: SGFIEL
C###  Description:
C###    ISFIEL(iw,ne) is the segment number of field iw in element ne.

C#### Variable: ISFIPR
C###  Type: INTEGER
C###  Set_up: SGFIPR
C###  Description:
C###    ISFIPR is the segment number of fibre angle profile.

C#### Variable: ISGAUS(iw,ng,ne)
C###  Type: INTEGER
C###  Set_up: SGGAUS
C###  Description:
C###    ISGAUS(iw,ng,ne) are the segment numbers of Gauss points.

C#### Variable: ISGRAD(ne,nr)
C###  Type: INTEGER
C###  Set_up: SGGRAD
C###  Description:
C###    ISGRAD(ne,nr) are the segment numbers of element gradient.

C#### Variable: ISGRID(iw)
C###  Type: INTEGER
C###  Set_up: SGGRID
C###  Description:
C###    ISGRID(iw) are the segment numbers of grid.

C#### Variable: ISHIST(np)
C###  Type: INTEGER
C###  Set_up: SGHIST
C###  Description:
C###    ISHIST(np) is the segment number of time history at node np.
C###    ISHIST(0) is the segment number of time history axes and labels.

C#### Variable: ISINCR(iw)
C###  Type: INTEGER
C###  Set_up: SGINCR
C###  Description:
C###    ISINCR(iw) is the segment number of increments.

C#### Variable: ISISOC(iw,noisoc)
C###  Type: INTEGER
C###  Set_up: SGISOC
C###  Description:
C###    ISISOC(iw,noisoc) are the segment numbers of isochrones.

C#### Variable: ISLead
C###  Type: INTEGER
C###  Set_up: SGLEAD
C###  Description:
C###    ISLead is segment containing axes and labels.

C#### Variable: ISLINE(iw,noline)
C###  Type: INTEGER
C###  Set_up: SGLINE
C###  Description:
C###    ISLINE(iw,noline) is segment number of lines at set noline.

C#### Variable: ISLINO(iw)
C###  Type: INTEGER
C###  Set_up: SGLINE
C###  Description:
C###    ISLINO(iw) is the segment number of line numbers.

C#### Variable: ISMAP(nomap)
C###  Type: INTEGER
C###  Set_up: SGMAP
C###  Description:
C###    ISMAP(nomap) is segment number of map.

C#### Variable: ISMATE(iw,ne)
C###  Type: INTEGER
C###  Set_up: SGMATE
C###  Description:
C###    ISMATE(iw,ne) are segment numbers of material in element ne.

C#### Variable: ISNONO(iw,np)
C###  Type: INTEGER
C###  Set_up: SGNODE
C###  Description:
C###    ISNONO(iw,np) are segment numbers of node numbers.

C#### Variable: ISOBJE(iw,noobje,nopart)
C###  Type: INTEGER
C###  Set_up: SGOBJE
C###  Description:
C###    ISOBJE(iw,noobje,nopart) are segment numbers of graphical
C###    objects.

C#### Variable: ISPLIN(iw,nr)
C###  Type: INTEGER
C###  Set_up: SGPLIN
C###  Description:
C###    ISPLIN(iw,nr) are segment numbers of polylines.

C#### Variable: ISPLOT(nh,ne,nr)
C###  Type: INTEGER
C###  Set_up: SGPLOT
C###  Description:

C#### Variable: ISPMAR(iw)
C###  Type: INTEGER
C###  Set_up: SGPMAR
C###  Description:
C###    ISPMAR(iw) is segment number of polymarker.

C#### Variable: ISPROF
C###  Type: INTEGER
C###  Set_up: SGPROF
C###  Description:

C#### Variable: ISREAC(iw)
C###  Type: INTEGER
C###  Set_up: SGREAC
C###  Description:
C###    ISREAC(iw) are segment numbers of reactions.

C#### Variable: ISRESI(iw)
C###  Type: INTEGER
C###  Set_up: SGRESI
C###  Description:
C###    ISRESI(iw) are segment numbers of residuals.

C#### Variable: ISRULE(iw)
C###  Type: INTEGER
C###  Set_up: SGRULE
C###  Description:
C###    ISRULE(iw) are segment numbers of rulers.

C#### Variable: ISSCAL(iw,nr)
C###  Type: INTEGER
C###  Set_up: SGSCAL
C###  Description:
C###    ISSCAL(iw,nr) is segment number of scale.

C#### Variable: ISSECT(nosect)
C###  Type: INTEGER
C###  Set_up: SGSECT
C###  Description:
C###    ISSECT(nosect) is segment number of section.

C#### Variable: ISSHEE(nw,ne,nr)
C###  Type: INTEGER
C###  Set_up: SGSHEE
C###  Description:
C###    ISSHEE(iw,ne,nr) are segment numbers of sheets.

C#### Variable: ISSTRA(ne)
C###  Type: INTEGER
C###  Set_up: SGSTRA
C###  Description:
C###    ISSTRA(ne) are segment numbers of principal strains.

C#### Variable: ISSTRE(ne,nostre)
C###  Type: INTEGER
C###  Set_up: SGSTRE
C###  Description:
C###    ISSTRE(ne,nostre) are segment numbers of principal stresses at
C###    set nostre.

C#### Variable: ISSTRM(ne,nostrm)
C###  Type: INTEGER
C###  Set_up: SGSTRM
C###  Description:
C###    ISSTRM(ne,nostrm) are segment numbers of streamlines.

C#### Variable: ISSURF(ne)
C###  Type: INTEGER
C###  Set_up: SGSURF
C###  Description:
C###    ISSURF(ne) are segment numbers of surface grid in element ne.

C#### Variable: LD(nd)
C###  Type: INTEGER
C###  Set_up: SGDATA, DEXI
C###  Description:
C###    LD(nd) is the element number associated with data point nd.
C###  See-Also: XID, XID

C#### Variable: LN(0:nl)
C###  Type: INTEGER
C###  Set_up: SGDATA, DEXI
C###  Description:
C###    LN(nl) is the element number for l=1..LN(0) in fitting problem.
C###    LN(0) is the number of elements in fitting.

C#### Variable: NDLT(ne)
C###  Type: INTEGER
C###  Set_up: SGDATA,DSTATS
C###  Description:
C###    NDLT(ne) is the number data points within line or face nl.

C#### Variable: NDDL(ne,nde)
C###  Type: INTEGER
C###  Set_up: SGDATA,DSTATS
C###  Description:
C###    NDDL(ne,nde) is the global data point number of element data
C###    point nde in ne.

C#### Variable: NDP(np)
C###  Type: INTEGER
C###  Set_up: SGDATA,IODATA,IOSIGN
C###  Description:
C###    NDP(np) is the reference number for data point nd.
C###  See-Also: NDDATA,ZD

C#### Variable: NDT
C###  Type: INTEGER
C###  Set_up: SGDATA,IODATA,IOSIGN
C###  Description:
C###    NDT is  total number of data points.

C#### Variable: NELIST(nolist)
C###  Type: INTEGER
C###  Set_up: SGDATA
C###  Description:
C###    NELIST(nolist) nolist=1,NELIST(0) is list of element numbers
C###    where data numbers or projections or parameters are defined.

C#### Variable: WD(nj,nd)
C###  Type: REAL*8
C###  Set_up: SGDATA
C###  Description:
C###    WD(nj,nd) is the weighting factor for data point nd.

C#### Variable: WDL(nj,nde)
C###  Type: REAL*8
C###  Set_up: SGDATA
C###  Description:
C###    WDL(nj,nde) is the weighting factor for element data point nde.

C#### Variable: XID(ni,nd)
C###  Type: REAL*8
C###  Set_up: SGDATA, DEXI
C###  Description:
C###    XID(ni,nd) is the Xi-coordinate of data point nd.
C###  See-Also: XIP, LD

C#### Variable: XIDL(ni,nde)
C###  Type: REAL*8
C###  Set_up: SGDATA
C###  Description:
C###    XIDL(ni,nde) is the Xi-coordinate of element data point nde.

C#### Variable: ZD(nj,nd)
C###  Type: REAL*8
C###  Set_up: SGDATA,IODATA,IOSIGN
C###  Description:
C###    ZD(nj,nd) is rectangular cartesian coordinates of data point nd.
C###  See-Also: NDDATA,NDP

C#### Variable: ZDD(nj,nd)
C###  Type: REAL*8
C###  Set_up: SGDATA
C###  Description:
C###    ZDD(nj,nd) is temporary array - must be dimensioned ZDD(3,*).
C###    USE_DATA and USE_GRAPHICS must be set to 1.

C#### Variable: ZDL(nh,nde)
C###  Type: REAL*8
C###  Set_up: SGDATA
C###  Description:
C###    ZDL(nh,nde) are the rectandular cartesian coordinates of
C###    element data point nde.


C#### Module: FE26
C###  Description:
C###    Routines called by the 'display' & 'import'/'export' commands.

C###  Routine: BLKCMGUI  (block data) initialization for cmgui link
C###  Routine: CMGUI_LINK_DESTROY destroy cmgui link
C###  Routine: CMGUI_LINK_FIELD_INFO determines field information
C###  Routine: CMGUI_LINK_DATA_FIELD_INFO determines data field info
C###  Routine: CMGUI_LINK_GET_DATA gets information about changed nodes
C###  Routine: CMGUI_LINK_INITIALISE initialise cmgui link
C###  Routine: CMGUI_LINK_UPDATE update cmgui data link
C###  Routine: CMGUI_LINK_UPDATE_DATACHG changed nodes
C###  Routine: CMGUI_LINK_UPDATE_DATADEL deleted nodes
C###  Routine: CMGUI_LINK_UPDATE_REGIONCHG changed nodes
C###  Routine: CMGUI_LINK_UPDATE_REGIONDEL deleted nodes
C###  Routine: CMGUI_LINK_UPDATE_NODECHG changed nodes
C###  Routine: CMGUI_LINK_UPDATE_NODEDEL deleted nodes
C###  Routine: DATA_CHANGE  CMGUI link - modify the data value/struct
C###  Routine: DATA_CHANGE_DYNAM  As above
C###  Routine: DATA_CREATE  CMGUI link - create the data
C###  Routine: DATA_CREATE_DYNAM  As above
C###  Routine: DATA_DESTROY CMGUI link - destroy the data
C###  Routine: DATA_DESTROY_DYNAM As above
C###  Routine: DIBASE   *** ARCHIVED ***
C###  Routine: DIFIBR   display fibre angle profile plots
C###  Routine: DIHIST   display history of nodal time variation
C###  Routine: DILEAD   display electrocardiographic leads
C###  Routine: DIPROF   display profiles of stress or strain
C###  Routine: DISECT   display section
C###  Routine: DITRAC   *** ARCHIVED ***
C###  Routine: CALCULATE_BASIS_STRING calculates the basis string
C###  Routine: EXPORT_FIELD_HEADING export field info
C###  Routine: EXPORT_FIELD_HEADING_PROP export field info for property
C###  Routine: EXELEM   export element info from FE data base
C###  Routine: EXFIELDML export to FieldML
C###  Routine: EXGEOM   export geometry info from FE data base ***ARCHIVED***
C###  Routine: EXGEOM_NE export geometry info for element NE ***ARCHIVED***
C###  Routine: EXGRID   export grid connectivity in exelem format
C###  Routine: EXNODE   export node info from FE data base
C###  Routine: EXPOIN   export point info from FE data base
C###  Routine: EXPROP   export properties from CE to element field
C###  Routine: EXSIGN   export signal info from FE data base
C###  Routine: EXSOUR   export source info from FE data base
C###  Routine: EXVORO   export Voronoi mesh from simplex elements
C###  Routine: EXTEXT   export texture info from grid pt data base
C###  Routine: IMGRID   import grid directly from data file
C###  Routine: IMGRID_TRUEGRID   import grid from TrueGrid file
C###  Routine: IMSIGN   import signal info into FE data base
C###  Routine: REGION_CHANGE  CMGUI link - modify the region val/struct
C###  Routine: REGION_CHANGE_DYNAM  As above
C###  Routine: REGION_CREATE  CMGUI link - create the region
C###  Routine: REGION_CREATE_DYNAM  As above
C###  Routine: REGION_DESTROY CMGUI link - destroy the region
C###  Routine: REGION_DESTROY_DYNAM As above

C#### Variable: NFLIST(0:nf)
C###  Type: INTEGER
C###  Set_up: EXELEM
C###  Description:
C###    NFLIST(0:nf) is a temporary array for the list of faces.
C###    However, on entry to DEFIT, NFLIST is sometimes expected to
C###    contain a default list of faces in the fit.

C#### Variable: ZCROSSING_INDEX(nytr)
C###  Type: INTEGER
C###  Set_up: EXSIGN
C###  Description:
C###    Stores an index of the ZCROSSING nodes for outputing and
C###    ordered list.

C#### Variable: ZCROSSING_RANGE(nytr,3)
C###  Type: REAL*8
C###  Set_up: EXSIGN
C###  Description:
C###    <HTML>
C###    Stores the <UL>
C###      <LI> maximum - ZCROSSING_RANGE(nytr,1)
C###      <LI> minimum - ZCROSSING_RANGE(nytr,2) and
C###      <LI> ranges - ZCROSSING_RANGE(nytr,3)
C###    </UL> of the ZCROSSING array.
C###    </HTML>

C#### Module: FE27
C###  Description:
C###    Routines called by 'group' & 'list' commands

C###  Routine: GROUPS   (block data) initialization of group arrays
C###  Routine: DUMESH   duplicate mesh
C###  Routine: RNMESH   renumber mesh
C###  Routine: GRELEM   group elements from command line call
C###  Routine: GRELEM_SUB group elements directly (no command line call)
C###  Routine: GRFACE   group faces
C###  Routine: GRGAUS   group Gauss points
C###  Routine: GRGRID   group grid points
C###  Routine: GRLINE   group lines
C###  Routine: GRNODE   group nodes
C###  Routine: GRNODE_SUB groups nodes from command line call
C###  Routine: GRPLIN   group polylines
C###  Routine: LIACTI   list active muscle model parameters
C###  Routine: LIANAL   list analytic formula parameters
C###  Routine: LIASSI   list user assigned strings
C###  Routine: LIBASE   list basis functions
C###  Routine: LICELL   list cell
C###  Routine: LICOLO   list workstation index colours
C###  Routine: LICONS   *REMOVED* list constant
C###  Routine: LICOOR   list coordinates
C###  Routine: LICORN   list corner nodes (for BE problems)
C###  Routine: LICOUP   list coupling
C###  Routine: LICUST   list customisation coefficents
C###  Routine: LIDATA   list data
C###  Routine: LIDELA   list delaunay
C###  Routine: LIEIGE   list eigenvalues
C###  Routine: LIELEM   list elements
C###  Routine: LIEQUA   list equations
C###  Routine: LIEXPO   list export parameters
C###  Routine: LIFACE   list faces
C###  Routine: LIFIT    list fit
C###  Routine: LIFUNC   list objective function
C###  Routine: LIGAUS   list Gauss point values of arrays
C###  Routine: LIGRID   list grid points
C###  Routine: LIGROW   list growth law parameters
C###  Routine: LIHEAD   list heading
C###  Routine: LIHIST   list history
C###  Routine: LIIMPO   list import parameters
C###  Routine: LIINCR   list increment vectors
C###  Routine: LIINIT   list initial & boundary conditions
C###  Routine: LIINVE   list inverse matrix regularisation parameters
C###  Routine: LIITER   list iteration parameters
C###  Routine: LILEAD   list electrocardiographic leads
C###  Routine: LILINE   list lines
C###  Routine: LIMAP    list mapping arrays
C###  Routine: LIMATE   list materials
C###  Routine: LIMATR   list matrix
C###  Routine: LIMESH   list mesh
C###  Routine: LIMODE   list modes (eigenvalues and eigenvectors)
C###  Routine: LIMOTI   list motion parameters
C###  Routine: LINODE   list nodes
C###  Routine: LINOIS   list noise
C###  Routine: LINORM   list normals for BEM
C###  Routine: LIOBJE   list objects
C###  Routine: LIOPTI   list optimisation parameters
C###  Routine: LIOUTP   list output
C###  Routine: LIPARA   list parameters
C###  Routine: LIPLIN   list polyline parameters
C###  Routine: LIREGI   list region data
C###  Routine: LIREGP   list regularisation parameters
C###  Routine: LIREFE   list reference location(s)
C###  Routine: LISAIL   list sail shape parameters
C###  Routine: LISEGM   list segments
C###  Routine: LISING   list physical singularity location (BE)
C###  Routine: LISIGN   list signal
C###  Routine: LISOLV   list solution parameters
C###  Routine: LISOUR   list sources
C###  Routine: LISTRA   list strains
C###  Routine: LISTRE   list stresses
C###  Routine: LITIME   list time variable information
C###  Routine: LITRSF   list transfer matrix parameters
C###  Routine: LIVARI   list variables
C###  Routine: LIVOLU   list volumes enclosed by selected elements.
C###  Routine: LIVORO   list voronoi
C###  Routine: LIWIND   list window dimensions
C###  Routine: LIXI     list xi coordinates

C#### Variable: GRNGLIST(0:negm)
C###  Type: INTEGER
C###  Set_up: GRGAUS
C###  Description:
C###    GRNGLIST(0:negm) is a temporary array for the list of
C###    Gauss points in elements.
C###  See-Also: GRGAUS

C#### Variable: IPOINTTYP
C###  Type: INTEGER
C###  Set_up: LISTRA
C###  Description:
C###    <HTML>
C###    IPOINTTYP describes the type of points at which OPSTRA outputs
C###    strain information.
C###    <PRE>
C###    IPOINTTYP = 1 for Gauss point strains,
C###    IPOINTTYP = 2 for data point strains,
C###    IPOINTTYP = 0 for strains at point along a xi-coordinate line.
C###    </PRE> </HTML>

C#### Variable: NGLIST(0:ng)
C###  Type: INTEGER
C###  Set_up: LIGAUS,LISTRA,LISTRE
C###  Description:
C###    NGLIST(0:ng) is a temporary array for the list of Gauss points.

C#### Module: FE28
C###  Description:
C###    Routines called by the 'draw' command.
C###  Routine: DRALIG   draw alignment
C###  Routine: DRAXES   draw axes
C###  Routine: DRCLOC   draw clock
C###  Routine: DRCONT   draw contours
C###  Routine: DRCROS   draw cross-section
C###  Routine: DRDATA   draw data
C###  Routine: DRDIPO   draw dipole
C###  Routine: DRELEM   draw elements
C###  Routine: DRFIBR   draw fibres
C###  Routine: DRFIEL   draw field
C###  Routine: DRGAUS   draw Gauss points
C###  Routine: DRGRAD   draw gradient vector
C###  Routine: DRGRID   draw finite difference grid
C###  Routine: DRINCR   draw increment
C###  Routine: DRISOC   draw isochrones
C###  Routine: DRLCURVE draw L-curve
C###  Routine: DRLINE   draw lines
C###  Routine: DRMATE   draw materials
C###  Routine: DRNODS   draw nodes
C###  Routine: DROBJE   draw graphical object
C###  Routine: DRPLIN   draw polyline
C###  Routine: DRPLOT   draw 3D PHIGS plot
C###  Routine: DRPMAR   draw polymarker
C###  Routine: DRREAC   draw reactions
C###  Routine: DRREGP   draw reg-parameters
C###  Routine: DRRESI   draw residuals
C###  Routine: DRRULE   draw ruled lines
C###  Routine: DRSCAL   draw scale
C###  Routine: DRSHEE   draw sheet angles
C###  Routine: DRSTRA   draw strains
C###  Routine: DRSTRE   draw stresses
C###  Routine: DRSTRM   draw streamlines
C###  Routine: DRSURF   draw surface
C###  Routine: DRVELO   draw velocities
C###  Routine: DRXI     draw Xi coordinates

C#### Variable: NTCOVA(ne)
C###  Type: INTEGER
C###  Set_up: DRCONT
C###  Description:
C###    NTCOVA(ne) is number of contours in element ne

C#### Variable: COVA(ne,nocova)
C###  Type: REAL*8
C###  Set_up: DRCONT
C###  Description: COVA(ne,nocova) holds contour values
C###    nocova=1,NTCOVA(ne) for element ne.


C#### Module: FE29
C###  Description:
C###    Evaluate Routines

C###  Routine: EVAERO    evaluates aerofoil lift & drag & and wake
C###  Routine: EVCONT    evaluates contraints
C###  Routine: EVCORO    evaluates coronary flow
C###  Routine: EVCOUP    evaluates coupling model
C###  Routine: EVELEC    evaluates electrodes
C###  Routine: EVERROR   evaluates strain energy error norm
C###  Routine: EVEVENT   evaluates events from signals
C###  Routine: EVFIBR    evaluates fibre angle at a specified point
C###  Routine: EVFIEL    evaluates geom/field/dep vars at Xi points
C###  Routine: EVINTE    evaluates integral of a field
C###  Routine: EVINVE    evaluates inverse of transfer matrix
C###  Routine: EVINVE_DYNAM evalutes epicardial inverse PHI_H
C###  Routine: EVLAPL    evaluates the surface Laplacian (for inverse problems)
C###  Routine: EVNOIS    evaluates the noise levels in signals
C###  Routine: EVOBJE    evaluates objective function
C###  Routine: EVPHI     evaluates phi array and rank/svd of it
C###  Routine: EVREAC    evaluates reactions
C###  Routine: EVRESI    evaluates residuals
C###  Routine: EVSENS    evaluates sensitivity of optimisation params
C###  Routine: EVSENS_CALC  does the calculations of the sensitivity 
C###  Routine: EVSIGN    evaluates a signal file
C###  Routine: EVSOLU    evaluates domain solutions at pts in BEM
C###  Routine: EVTIME    evaluates time variables
C###  Routine: EVTIME_FCN evaluates time variables
C###  Routine: EVTRSF    buffer routine for EVTRSF_DYNAM
C###  Routine: EVTRSF_DYNAM evaluates transfer matrix
C###  Routine: EVZEROXING evaluates zerocrossing function (inverse soln)

C JMB 22-FEB-2000 Routines for inverse stuff

C###  Routine: CRESO              Evaluates CRESO function
C###  Routine: CRESOFUN      (fn) Auxillary routine for CRESO
C###  Routine: FMIN          (fn) Find global minimum
C###  Routine: FOURIERCOEFFS      Evaluates Fourier coefficients
C###  Routine: FZERO         (fn) Finds zero
C###  Routine: GCV                Evaluates GCV
C###  Routine: GCVFUN        (fn) Auxillary routine for GCV
C###  Routine: GSVALUES           Evaluates the generalised singular values
C###  Routine: INTEREQRANK        Evaluates the inter-equation truncation rank
C###  Routine: LCURVE             Evaluates LCURVE function
C###  Routine: LCFUN         (fn) Auxillary routine for LCURVE
C###  Routine: PICARD             Evaluates PICARD criterion
C###  Routine: QUASIOPT           Evaluates Quasi-opt criterion
C###  Routine: QUASIFUN      (fn) Auxillary routeinf ro QUASIOPT
C###  Routine: REOPT              Evaluates relative optimium
C###  Routine: REFUN         (fn) Auxillary routine for REOPT
C###  Routine: REIDFUN       (fn) Evaluates the intrinsic residual
C###  Routine: TGSVD              Evaluates truncated SVD
C###  Routine: TIKHONOV           Evaluates Tikhonov
C###  Routine: ZEROCROSSING       Evalutes zerocrossing function
C###  Routine: ZEROFUN       (fn) Auxillary routine for ZEROCROSSING

C#### Variable: ISIZE_MFI(3,NSSM)
C###  Type: INTEGER
C###  Set_up: EVSOLU
C###  Description:
C###    <HTML>
C###    Stores the sizes of the MFI array.
C###    <UL>
C###      <LI> 1st index stores the (1) number of electrodes,
C###           (2) number of time steps, (3) number of components
C###      <LI> The 2nd index stores the signal set the
C###           1st index refers to.
C###    </UL>
C###    </HTML>
C###  See-Also: MFI

C#### Variable: ISIZE_PHI(2)
C###  Type: INTEGER
C###  Set_up: EVPHI
C###  Description:
C###    <HTML>
C###    <UL>
C###      <LI> ISIZE_PHI(1) is the number of rows of PHI (electrodes)
C###      <LI> ISIZE_PHI(2) is the number of cols of PHI (time)
C###    </UL>
C###    </HTML>
C###  See-Also: PHI, ISIZE_PHIH, PHI_H, PHI_H_EXACT

C#### Variable: ISIZE_PHIH(2)
C###  Type: INTEGER
C###  Set_up: EVPHI
C###  Description:
C###    <HTML>
C###    <UL>
C###      <LI> ISIZE_PHIH(1) is the number of rows of PHIH (electrodes)
C###      <LI> ISIZE_PHIH(2) is the number of cols of PHIH (time)
C###    </UL>
C###    </HTML>
C###  See-Also: PHI, ISIZE_PHI, PHI_H, PHI_H_EXACT

C#### Variable: ISIZE_PHI(2)
C###  Type: INTEGER
C###  Set_up: EVPHI
C###  Description:
C###    <HTML>
C###    The size of the Heart PHI matrix. Not the size of PHIH computed
C###    by APTRSF using the activation transfer matrix.
C###    <UL>
C###      <LI> ISIZE_PHIH(1) is the number of rows of PHIH (electrodes)
C###      <LI> ISIZE_PHIH(2) is the number of cols of PHIH (time)
C###    </UL>
C###    </HTML>
C###  See-Also: PHI, ISIZE_PHIH, PHI_H, PHI_H_EXACT

C#### Variable: ISIZE_TBH(2)
C###  Type: INTEGER
C###  Set_up: EVTRSF
C###  Description:
C###    <HTML>
C###    <UL>
C###    <LI> ISIZE_TBH(1) is the number of rows of T_BH
C###    <LI> ISIZE_TBH(2) is the number of cols of T_BH
C###    </UL>
C###    </HTML>
C###  See-Also: T_BH

C#### Variable: DLL(3,nl)
C###  Type: REAL*8
C###  Set_up: EVAERO
C###  Description:

C#### Variable: FIXP(2,ne)
C###  Type: LOGICAL
C###  Set_up: EVSENS
C###  Description:
C###    FIXP(2,ne) set wiether pressure is incremented when
C###    determining sensitivity of optimisation parameters
C###    wrt geometry pressure bcs

C#### Variable: NEERR(NEM,3)
C###  Type: REAL*8
C###  Set_up: EVERROR
C###  Description:
C###    <HTML>
C###    NEERR(ne,1) contains the relative percentage error and <BR>
C###    NEERR(ne,2) contains the element error norm <BR>
C###    This variable may only be temporary as I experiment with
C###    error estimation techniques CS
C###    <BR> Memory Allocation removed for now - 12-JUL-1999
C###    </HTML>

C#### Variable: NQNP
C###  Type: INTEGER
C###  Set_up: UPNODE
C###  Description:
C###    NQNP(np) gives the number of the global grid point located
C###    at node np.

C#### Variable: PE(2,ne)
C###  Type: REAL*8
C###  Set_up: EVSENS
C###  Description:
C###    PE(1..2,ne) is the pressure increment on element face.

C#### Variable: MFI(nd,nt,3,nss)
C###  Type: REAL*8
C###  Set_up: EVSOLU
C###  Description:
C###    <HTML>
C###    MFI is the matrix of Magnetic Field Intensity (H).
C###    where magnetic flux intensity B=uH.
C###    <BR>
C###    It has dimensions:
C###      MFI(electrodes,timesteps,fieldcomponents,samplesets).
C###    <BR>NOTE: nss=1 is where the known solution is stored while
C###    nss=2 is where the estimated solution is stored.
C###    </HTML>
C###  See-Also: ISIZE_MFI

C#### Variable: PHI(nytr,nts)
C###  Type: REAL*8
C###  Set_up: EVPHI
C###  Description:
C###    PHI(nytr,nts) is the matrix of (measured) signals
C###    (the data matrix).
C###  See-Also: ISIZE_PHI

C#### Variable: PHI_H(nytr,nts)
C###  Type: REAL*8
C###  Set_up: EVINVE
C###  Description:
C###    PHI_H(nytr,nts) is the matrix of inverted heart potentials.
C###    Also used as intermediate storage for torso potentials during
C###    inverse calculation of isochrones.

C#### Variable: PHI_H_EXACT(nytr,nts)
C###  Type: REAL*8
C###  Set_up: EVPHI
C###  Description:
C###    PHI_H_EXACT(nytr,nts) is the matrix of exact (measured) heart
C###    signals.

C#### Variable: REG_PARAMETER(nts)
C###  Tupe: REAL*8
C###  Set_up: EVINVE
C###  Sescription:
C###    REG_PARAMETER(nts) is the vector containing the regularisation
C###    parameters for each time step.

C#### Variable: U_PHI(nytr,nytr)
C###  Type: REAL*8
C###  Set_up: EVPHI
C###  Description:
C###    SVD of PHI is U_PHI*SIGMA_PHI*VT_PHI

C#### Variable: U_T_BH(nytr,nytr)
C###  Type: REAL*8
C###  Set_up: EVZEROXING,EVTRSF
C###  Description:
C###    SVD of T_BH is U_T_BH*SIGMA_T_BH*VT_T_BH

C#### Variable: VT_PHI(nts,nts)
C###  Type: REAL*8
C###  Set_up: EVPHI
C###  Description:
C###    SVD of PHI is U_PHI*SIGMA_PHI*VT_PHI

C#### Variable: VT_T_BH(nytr,nytr)
C###  Type: REAL*8
C###  Set_up: EVZEROXING,EVTRSF
C###  Description:
C###    SVD of T_BH is U_T_BH*SIGMA_T_BH*VT_T_BH

C#### Variable: SIGMA_PHI(nytr)
C###  Type: REAL*8
C###  Set_up: EVPHI
C###  Description:
C###    SVD of PHI is U_PHI*SIGMA_PHI*VT_PHI

C#### Variable: SIGMA_T_BH(nytr)
C###  Type: REAL*8
C###  Set_up: EVZEROXING,EVTRSF
C###  Description:
C###    SVD of T_BH is U_T_BH*SIGMA_T_BH*VT_T_BH

C#### Variable: NTST
C###  Type: INTEGER
C###  Set_up: EVPHI
C###  Description:
C###    Total number of time samples.

C#### Variable: ZCROSSING(nytr,nts)
C###  Type: REAL*8
C###  Set_up: EVZEROXING
C###  Description:
C###    ZCROSSING(nyr,nts) is the zero crossing matrix obtained using
C###    using the algorithm of Greensite and Huiskamp.
C###  See-Also: PHI

C#### Variable: SVD_CUTOFF_RATIO
C###  Type: REAL8*
C###  Set_up: EVPHI
C###  Description:
C###    Ratio of singular values (sigma/sigma_max) to be used in
C###    determining the rank of PHI

C#### Variable: T_BH(nytr,nytr)
C###  Type: REAL*8
C###  Set_up: EVTRSF
C###  Description:
C###    T_BH is the transfer matrix from surface h (heart) to surface
C###    b (body).
C###  See-Also: ISIZE_TBH

C#### Variable: T_BH_INV(nytr,nytr)
C###  Type: REAL*8
C###  Set_up: EVINVE
C###  Description:
C###    T_BH_INV is the inverse of the transfer matrix.  It either
C###    stores the inverse explicity, or the regularised matrix from
C###    which the solution is calculated in APINVE.

C#### Variable: WK1_INV(nytr,nytr)
C###  Type: REAL*8
C###  Set_up: EVINVE,EVTRSF
C###  Description:
C###    A work array used in the construction of the transfer (T_BH) and
C###    inverse transfer (T_BH_INV) arrays

C#### Variable: WK2_INV(nytr,nytr)
C###  Type: REAL*8
C###  Set_up: EVINVE,EVTRSF
C###  Description:
C###    A work array used in the construction of the transfer (T_BH) and
C###    inverse transfer (T_BH_INV) arrays

C#### Variable: WK3_INV(nytr,nytr)
C###  Type: REAL*8
C###  Set_up: EVINVE,EVTRSF
C###  Description:
C###    A work array used in the construction of the transfer (T_BH) and
C###    inverse transfer (T_BH_INV) arrays

C#### Variable: WK4_INV(nytr)
C###  Type: REAL*8
C###  Set_up: EVINVE,EVTRSF
C###  Description:
C###    A work array used in the construction of the transfer (T_BH) and
C###    inverse transfer (T_BH_INV) arrays


C#### Module: FE30
C###  Description:
C###    General pde  routines

C###  Routine: ANS001   calculates 2d diffusion analytical solution
C###  Routine: ANS002   calculates 1d diffusion analytical solution
C###  Routine: BFRONT    to cmiss_gridarchive.f
C###  Routine: BLK30     (block data) set up parameter titles in common block
C###  Routine: BRANCH1   calculates gridpoint values at bifurcations
C###  Routine: BRANCH2  evaluates micro-circulation model parameters
C###  Routine: BRANCH3  calculates transition pts (art-vein) parameters    
C###  Routine: CALCMASS  calculates mass in lung model
C###  Routine: CALC_BBM  calculates regression parameters for BB model
C###  Routine: CALC_DIFFUSION calculate grid diffusion explicitly
C###  Routine: CALC_DTAR handles the addition/removal dynamic tracking
C###  Routine: CALC_FE_GRID_COEF grid-based FE coefficients for nq's
C###  Routine: CALC_FV_GRID_COEF grid-based FV coefficients for nq's
C###  Routine: CALC_FV_GRID_SECFLUX grid-based FV secondary flux for all nq's
C###  Routine: CALC_GRID_BOUND_COEF calculation of flux coefficients
C###  Routine: CALC_GRID_COEF coefficients an positions of nonzero nq's
C###  Routine: CALC_GRID_XI calcs xi position within quadratic elem
C###  Routine: CALC_HEMODYNAMICS calcs hemodynamic properties
C###  Routine: CALC_HEMATOCRIT calcs hematocrit for pulm capillaries
C###  Routine: CALC_LATT_COEF calc lattice matrix coefs for nq's
C###  Routine: CALC_LATT_GDOTN calc the tensor GDOTN 
C###  Routine: CALC_STIMULUS calcs stimulus and pseudo stim. currents
C###  Routine: CALC_STIMULUS2 calcs stimulus and pseudo stim. currents
C###  Routine: CALC_TRANSMEMBRANE_RHS calcs RHS vector for Vm
C###  Routine: CALC_TRANSIT_TIME calcs RBC transit time through network
C###  Routine: CALC_R_SEGMENT calcs resistance in capillary segment
C###  Routine: CALC_WALLTEMP calcs temperature at air-mucus interface
C###  Routine: CELL_ASSIGN_SPATIAL handles spatially varying cell params
C###  Routine: CAP_DIMENSION calcs dimensions in pulm capillaries
C###  Routine: CHECK_CONV check YQ-YP convergence
C###  Routine: CHECK_MESH finds elements involved in lung solution
C###  Routine: GEN_EXT_RHS rhs vector for extracellular implicit grid
C###  Routine: GEN_GRID_BEM_RHS generate RHS vector from BEM YP array
C###  Routine: GEN_GRID_MAP makes grid/gauss and gauss/grid maps
C###  Routine: GEN_INT_RHS to cmiss_gridarchive.f
C###  Routine: GET_EXT_CONTRIB phi-e component of Vm equation
C###  Routine: GET_FD_POINTS generates local quadratic grid schemes
C###  Routine: GETYQYMAP get YQ to Y mapping for grid solutions
C###  Routine: GGRADPHIQDN calulate grid flux
C###  Routine: GRADPHIQRC  grad phi in cartesians for grid points
C###  Routine: IONIC_CURRENT to cmiss_gridarchive.f
C###  Routine: INIT_PULM sets up lung parameters
C###  Routine: INIT_PULM_CIRC Sets up arrays for pulmonary capillary
C###  Routine: IPANA3    input for analytic functions
C###  Routine: IPBAS3    input for basis functions
C###  Routine: IPINI3    input initial conditions and bdry conditions
C###  Routine: IPMAT3    input for material parameters
C###  Routine: IPMAT3_CELL  input for cellular spatial variance
C###  Routine: IPMAT3_CELL1 input for cellular spatial variance
C###  Routine: IPMAT3_COUP input for coupling parameters
C###  Routine: LINE_FROM_2_PTS eqn of a line from 2 points
C###  Routine: MECH_FLOW calculates vessel transmural stress
C###  Routine: MESH_DEFORM
C###  Routine: MESH_FLOW calculates flow in 1D mesh
C###  Routine: NORM30    find normal vector to surface for bdry elems
C###  Routine: NORM31    find normal vector for grid problems
C###  Routine: NORM_LATTICE find normal vector for lattice grid problems
C###  Routine: NQDS      find arc length derivative of a grid point
C###  Routine: OPC30     print solutions and reactions
C###  Routine: OPLUNG    output lung information during solution
C###  Routine: OPMAT3    output material parameters
C###  Routine: OPINI3    output initial and boundary data
C###  Routine: PETROV_PREP multipliers for Petrov-Galerkin FEM
C###  Routine: PLANE_FROM_3_PTS eqn of a plane from 3 points
C###  Routine: SET_LUNG_BC sets boundary conditions for lung model
C###  Routine: SORTAC1   add new pt [SSTA(NSA)] in right place in SSTA
C###  Routine: SORTAC2   sort existing SSTA after one point altered
C###  Routine: TFRONT    to cmiss_gridarchive.f
C###  Routine: LINREGRESS finds linear regression equation
C###  Routine: LPM_EVAL  Method of Characteristics for lung model
C###  Routine: TREE_GEOM geometric info for purkinje trees
C###  Routine: TREE_POINTS random point generation for Purkinje trees
C###  Routine: TREE_START tree start calculation for Purkinje trees
C###  Routine: USER3_CORONARY1 function for coronary balck box eqn
C###  Routine: USER3_CORONARY2 function for coronary balck box eqn
C###  Routine: UPDATE_BBM updates BBA model parameters for lung models
C###  Routine: WBC_BLOCK calcs elems blocked by WBCs
C###  Routine: WBC_TRANSIT calcs WBC transit time in pulm capillaries
C###  Routine: XPEQ30    calc element flux matrices for linear eqs
C###  Routine: XPES30    calc element matrices for linear time-dep eqs
C###  Routine: XPFD30    calc molecule matrices for finite diff eqs
C###  Routine: ZEES30    calc element tangent stiff. mtrx. at curr. sol.
C###  Routine: ZERE30    calc element residuals
C###  Routine: ZFFS30    calc face tangent stiffness matrices
C###  Routine: ZFRF30    calc surface residual components

C**** Data statements for material parameter titles and default values
C**** RDMATE(nm,ityp5,ityp2,ityp3) contains material parameters nm
C**** equiv to DSTATIC_ityp2_ityp3(nm) for ityp5=1 (static analysis)
C****   "      DDYNAM       "            "        2 (time integration)
C****   "      DMODAL       "            "        3 (modal analysis)
C****   "      DQUASIST     "            "        4 (Quasistatic analy)
C****   "      DFRONT       "            "        5 (Front path analy)
C****   "      DBUCKLE      "            "        6 (Buckling analysis)
C**** equiv to CSTATIC_ityp2_ityp3(nm) for ityp5=1 (static analysis)
C****   "      CDYNAM       "            "        2 (time integration)
C****   "      CMODAL       "            "        3 (modal analysis)
C****   "      CQUASIST     "            "        4 (Quasistatic analy)
C****   "      CFRONT       "            "        5 (Front path analy)
C****   "      CBUCKLE      "            "        6 (Buckling analysis)

C#### Variable: AQ(maq,nq)
C###  Type: REAL*8
C###  Set_up: IPINI3
C###  Description:
C###    AQ holds what auxiliary parameters are defined at grid
C###    point nq. The maq values are allocated by MAQ_LOC.
C###    Currently it is used to hold current injection information.

C#### Variable: ILP(il,ie,nr,nx)
C###  Type: INTEGER
C###  Set_up: IPMAT3
C###  Description:
C###    <HTML> <PRE>
C###    ILP(il,ie,nr,nx),il=1,ILTOT(NJT,ie,ITYP2(nr,nx),ITYP3(nr,nx)),
C###    records whether equation nx parameters are:
C###       1 Constant everywhere - value  in CE(il,ne)
C###       2 Piecewise constant - values in CE(il,ne)
C###       3 Piecewise linear - values in CP(il,np)
C###       4 Defined by global points - values in CQ(il,nq)
C###       5 Defined elsewhere (eg in IPMESH)
C###    and if parameter IL is time varying ILP(il,ie,nr,nx) is negative
C###    where:
C###       ILTOT is ILTOT1 for static   analysis equations
C###       ILTOT is ILTOT2 for time  integration equations
C###       ILTOT is ILTOT3 for modal    analysis equations
C###       ILTOT is ILTOT4 for Fourier  analysis equations
C###       ILTOT is ILTOT5 for front path analysis equations
C###       ILTOT is ILTOT6 for buckling analysis equations
C###    NB. ie is always 1 except for linear elasticity problems.
C###    </PRE> </HTML>

C#### Variable: ITHRES(i,ng,ne)
C###  Type: INTEGER
C###  Set_up: TFRONT
C###  Description:
C###    ITHRES(1,ng,ne is 0/1/2 if Gauss point ng of element ne is
C###    inactive/active/surrounded by active points.
C###    ITHRES(2,ng,ne) is 0/1 if Gauss point ng is ordinary
C###    myocardium/Purkinje tissue.

C#### Variable: KTYPMBC
C###  Type: INTEGER
C###  Set_up: IPINI3
C###  Description:
C###    KTYPMBC is 1,2 if the mixed boundary condition flux is an
C###    integrated value/a point value.

C#### Variable: NQT
C###  Type: INTEGER
C###  Set_up: IPGRID
C###  Description:
C###    NQT is the number of collocation points in mesh.

C#### Variable: CQ(nm,nq)
C###  Type: REAL*8
C###  Set_up: IPMAT3
C###  Description:
C###    CQ(nm,nq) is the material parameter at collocation point nq.

C#### Variable: ED(nhs,nhs)
C###  Type: REAL*8
C###  Set_up: XPES30
C###  Description:
C###    ED(nhs,nhs) is the element damping matrix (first order time
C###    derivatives).

C#### Variable: EM(nhs,nhs)
C###  Type: REAL*8
C###  Set_up: XPES30
C###  Description:
C###    EM(nhs,nhs) is the element mass matrix (second order time
C###    derivatives).

C#### Variable: ES(nhs,nhs)
C###  Type: REAL*8
C###  Set_up: XPES30
C###  Description:
C###    ES(nhs,nhs) is the element stiffness matrix (zeroth order
C###    time serivatives).

C#### Variable: NE_OLD(NORM)
C###  Type: INTEGER
C###  Set_up: LPM_EVAL, IPMESH2
C###  Description:
C###    NE_OLD is a temporary array used by IPMESH2 for generating
C###    lung meshes, and by LPM_EVAL for calculating x*

C#### Variable: NE_TEMP(NORM)
C###  Type: INTEGER
C###  Set_up: LPM_EVAL, IPMESH2
C###  Description:
C###    NE_TEMP is a temporary array used by IPMESH2 for generating
C###    lung meshes, and by LPM_EVAL for calculating x*

C#### Variable: NQGP(0:NQGM,nq)
C###  Type: INTEGER
C###  Set_up: XQXE,GET_FD_POINTS
C###  Description:
C###    NQGP holds the grid points surrounding nq in the order
C###    in which they are required to create implicit finite
C###    difference matricies (i.e. ascending).  An entire matrix row is
C###    defined for each nq.
C###  See-Also: NQGP_PIVOT,NQGW

C#### Variable: NQGP_PIVOT(NQGM,nq)
C###  Type: INTEGER
C###  Set_up: GET_FD_POINTS
C###  Description:
C###    NQGP_PIVOT holds the reordering information for NQGW
C###    so the grid points can be extracted in ascending order
C###    into the sparse matrix structures.
C###  See-Also: NQGP,NQGW

C#### Variable: NQGW(NQGM,nq)
C###  Type: REAL*8
C###  Set_up: CALC_GRID_COEF
C###  Description:
C###    NQGW hold the matrix coefficients at the points surrounding
C###    nq as given in NQGW(NQGP_PIVOT(nzero,nq),nq) corresponds to grid
C###    point NQGP(nzero,nq).
C###  See-Also: NQGP,NQGP_PIVOT

C#### Variable: NTMP(0:500,NEM)
C###  Type: INTEGER
C###  Set_up: LPM_EVAL
C###  Description:
C###    NTMP stores the element numbers to begin searching for x*
C###    for the current element

C#### Variable: THRES(i,ng,ne)
C###  Type: REAL*8
C###  Set_up: TFRONT
C###  Description:
C###    THRES(1,ng,ne) is time since Gauss point became active;
C###    THRES(2,ng,ne) is membrane potential at Gauss point;
C###    THRES(3,ng,ne) is recovery variable  at Gauss point.

C#### Variable: PULMAT(nm), when ITYP3(nr,nx).LE.2 - airways
C###  Type: REAL*8
C###  Set_up: IPMAT3_CONST
C###  Description:
C###    PULMAT(1) is Gas diffusion coefficient/water vapour diffusivity
C###    PULMAT(2) is Oxygen loss coefficient / thermal diffusivity
C###    PULMAT(3) is Density of air+water vapour
C###    PULMAT(4) is Specific heat of air+water vapour
C###    PULMAT(5) is Thermal conductivity
C###    PULMAT(6) is Heat of evaporation

C#### Variable: PULMAT(nm), when ITYP3(nr,nx).EQ.3 - capillaries
C###  Type: REAL*8
C###  Set_up: IPMAT3_CONST
C###  Description:
C###    PULMAT(1) is the flux cutoff parameter, r
C###    PULMAT(2) is the preferential flux parameter, b
C###    PULMAT(3) is the capillary wall stiffness, kc (cm H2O)

C#### Variable: XIQ(ni,nq)
C###  Type: REAL*8
C###  Set_up: GEN_GRID_POTE_RHS,CALC_LATTICE_XIQ
C###  Description:
C###    XIQ(ni,nq) is used to store the xi position of grid point
C###    nq. The element number is stored in NENQ.

C#### Variable: YQ(nyq,niq,na,nx)
C###  Type: REAL*8
C###  Set_up: BFRONT,MG_COLLOCATE
C###  Description:
C###    <HTML>
C###    YQ(nyq,niq,na,nx) stores grid point state information.
C###    <PRE>
C###    NOTE: for activation problems nyq=nq
C###    For multigrid applications:
C###      YQ(nyq,1,na,nx) is current solution at grid level na
C###      YQ(nyq,2,na,nx) is residual or restriction
C###      YQ(nyq,3,na,nx) is residual or restriction
C###      YQ(nyq,4,na,nx) is restriction
C###      YQ(nyq,5,na,nx) is source term
C###      YQ(nyq,6,na,nx) is soln at previous time step for dynamic eqn
C###    For myocardial activation process:
C###     Ionic current variable values are
C###     For Panfilov FHN:
C###      YQ(nyq,1,1,nx) is transmembrane potential at grid point nq;
C###      YQ(nyq,2,1,nx) is extracellular potential at grid point nq;
C###      YQ(nyq,3,1,nx) is recovery variable (if needed) at point nq.
C###      YQ(nyq,4,1,nx) is calcium level.
C###      YQ(nyq,5,1,nx) is SAC current.
C###     For VCD with Calif. mods only (for Vdot and Tdot):
C###      YQ(nyq,1,1,nx) is transmembrane potential at grid point nq;
C###      YQ(nyq,2,1,nx) is extracellular potential at grid point nq;
C###      YQ(nyq,3,1,nx) is recovery variable (if needed) at point nq.
C###      YQ(nyq,4,1,nx) is calcium level.
C###      YQ(nyq,5,1,nx) is change in recovery variable.
C###      YQ(nyq,6,1,nx) is change in transmembrane potential;
C###      YQ(nyq,7,1,nx) is SAC current (-ve=inwards)
C###     For Beeler-Reuter:
C###      YQ(nyq,1,1,nx) is transmembrane potential at grid point nq;
C###      YQ(nyq,2,1,nx) is extracellular potential at grid point nq;
C###      YQ(NYQ,3,1,nx) is x1
C###      YQ(nyq,4,1,nx) is Cai
C###      YQ(nyq,5,1,nx) is m
C###      YQ(nyq,6,1,nx) is h
C###      YQ(nyq,7,1,nx) is j
C###      YQ(nyq,8,1,nx) is d
C###      YQ(nyq,9,1,nx) is f
C###     For Luo-Rudy:
C###      YQ(nyq,1,1,nx) is transmembrane potential at grid point nq;
C###      YQ(nyq,2,1,nx) is extracellular potential at grid point nq;
C###      YQ(NYQ,3,1,nx) is x
C###      YQ(nyq,4,1,nx) is Cai
C###      YQ(nyq,5,1,nx) is m
C###      YQ(nyq,6,1,nx) is h
C###      YQ(nyq,7,1,nx) is j
C###      YQ(nyq,8,1,nx) is d
C###      YQ(nyq,9,1,nx) is f
C###     For diFrancesco-Noble:
C###      YQ(nyq,1,1,nx) is transmembrane potential at grid point nq;
C###      YQ(nyq,2,1,nx) is extracellular potential at grid point nq;
C###      YQ(nyq,3,1,nx) is m at grid point nq;
C###      YQ(nyq,4,1,nx) is h at grid point nq.
C###    Miscellaneous other grid point information:
C###      YQ(nyq,1,3,nx) is activation time (in ms) of grid point;
C###      YQ(nyq,2,3,nx) used for temporary storage of phi(m) at nq;
C###      YQ(nyq,3,3,nx) used for temporary storage of phi(e) at nq;
C###      YQ(nyq,4,3,nx) is abs largest value of change in potential.
C###    </PRE>
C###    </HTML>

C#### Variable: FIXQ(nyq,niy_fix,nx)
C###  Type: LOGICAL
C###  Set_up: IPINI3
C###  Description:
C###    <HTML>
C###    FIXQ(nyq,niy_fix,nx) stores what type of boundary condition
C###    (if any) is applied on each nyq.
C###    <PRE>
C###    For activation problems:
C###      niy_fix=1 is a dependent variable boundary condition
C###      niy_fix=2 is a flux boundary condition
C###      niy_fix=3 is an analytic boundary condition
C###        (extracellular only)
C###    </PRE>
C###    </HTML>


C#### Module: FE40
C###  Description:
C###    Linear elasticity routines.

C###  Routine: BLK40     (block data) set up parameter titles in common block
C###  Routine: PHT3      (fn) eval 1D cubic Hermite  basis
C###  Routine: PSI20     (fn) eval tensor product Lagrange & Hermite basis
C###  Routine: ACTIVE    compute active stiffness terms
C###  Routine: AZTG45    calc membrane element stresses
C###  Routine: BEAM      calc element matrix for beam elements
C###  Routine: CGC9      calc coefficents for 3D elasticity
C###  Routine: CGC11     calc coefficents for plane stress,strain
C###  Routine: CGS9      calc stress,strain for 3D elasticity
C###  Routine: CGS11     calc stress,strain for plane stress,strain
C###  Routine: ELAS3D    calc element matrix for 3D elasticity
C###  Routine: GAUS20    eval Gauss point arrays for shell elements
C###  Routine: GGRRM     eval Gauss point arrays for shell elements
C###  Routine: IPANA4    input for analytic functions
C###  Routine: IPBAS4    input for basis functions
C###  Routine: IPINI4    input i.c.s and boundary conditions
C###  Routine: IPMAT4    input material parameters
C###  Routine: MIDMN     compute midsurface shell stress resultants
C###  Routine: NORM40    find normal vector in r.c.coords for shells
C###  Routine: OPC40     prints solution variables and reactions
C###  Routine: OPMAT4    output equation parameters for all elements
C###  Routine: OPMAT41   output equation parameters for one element
C###  Routine: OPST40    output stresses and strains at Gauss points
C###  Routine: OPINI4    output initial and boundary data
C###  Routine: PGGR      finds strain functions for shell elements
C###  Routine: PLANE     calc element matrix for plane stress & strain
C###  Routine: PLATE     calc element matrix for plate elements
C###  Routine: SHELL     calc shell element stiffness matrix
C###  Routine: SSGGX     calc shell quantities
C###  Routine: TRUSS     calc element matrix for truss elements
C###  Routine: X3XG      eval Gauss pt array of 3rd derivs of X wrt Xi
C###  Routine: XPES40    calc elem stiffness matrix for linear elast.
C###  Routine: ZEAZ45    calc metric quantities for membrane elements
C###  Routine: ZERE40    calc element residuals
C###  Routine: ZETG41    calc truss element stresses
C###  Routine: ZETG42    calc batten element stresses
C###  Routine: ZETG43    calc beam element stresses
C###  Routine: ZETG44    calc link element stresses
C###  Routine: ZEZG20    calc Gauss point arrays for shell elements
C###  Routine: ZGCHK     check for underflow in shell elements
C###  Routine: ZPZE20    transfer global to element arrays for shells


C#### Variable: D3PG(i,j,k)
C###  Type: REAL*8
C###  Set_up: GAUS20
C###  Description:
C###    D3PG(i,j,k) is the basis function Gauss point array for
C###    Lagrange/Hermite tensor product type basis function nb.

C#### Variable: GG(2,2,48,*)
C###  Type: REAL*8
C###  Set_up: PGGR
C###  Description:
C###    GG(2,2,48,*) are the strain functions of the basis function
C###    PG for all Gauss points of the current element ne.

C#### Variable: KTYP45
C###  Type: INTEGER
C###  Set_up: UNKNOWN (not known to be initialised 12/02/99)
C###  Description:
C###    KTYP45 is 1..4  for type of beam cross-section
C###    KTYP45 = 1 beam cross-section is rectangular solid
C###      "    = 2 beam cross-section is rectangular tube
C###      "    = 3 beam cross-section is ellipsoidal solid
C###      "    = 4 beam cross-section is ellipsoidal tube

C#### Variable: RE(ns,nh)
C###  Type: REAL*8
C###  Set_up: ZERE40
C###  Description:
C###    RE(ns,nh) is the element residual.

C#### Variable: RR(2,2,48,*)
C###  Type: REAL*8
C###  Set_up: PGGR
C###  Description:
C###    RR(2,2,48,*) is the curvature change function of the basis
C###    functions PG for all Gauss points of the current element ne.

C#### Variable: SMX(i,j)
C###  Type: REAL*8
C###  Set_up: MIDMN
C###  Description:
C###    SMX(i,j) (i,j=1,ni) is the moment resultant.

C#### Variable: SNX(i,j)
C###  Type: REAL*8
C###  Set_up: MIDMN
C###  Description:
C###    SNX(i,j) (i,j=1,ni) is the direct stress resultant.

C#### Variable: X3G(nd,nj,ng)
C###  Type: REAL*8
C###  Set_up: X3XG
C###  Description:
C###    <HTML> <PRE>
C###    X3G(1,nj) is d^3(Xj)/d(Xi1)^3;
C###    X3G(2,nj) is d^3(Xj)/d(Xi1)^2/d(Xi2);
C###    X3G(3,nj) is d^3(Xj)/d(Xi1)/d(Xi2)^2;
C###    X3G(4,nj) is d^3(Xj)/d(Xi2)^3.
C###    </PRE> </HTML>

C#### Variable: ZGG(nu,nh)
C###  Type: REAL*8
C###  Set_up: ZEZG20
C###  Description:
C###    ZGG(nu,nh), nu=1,NUT(nb) is a Gauss point array.

C#### Variable: EKC(ms,ns)
C###  Type: COMPLEX*16
C###  Set_up: ACTIVE
C###  Description:
C###    EKC(ms,ns) are active stiffness components.


C#### Module: FE50
C###  Description:
C###    Finite deformation elasticity routines.

C###  Routine: BLK50     (block data) set up parameters in common block
C###  Routine: AG        (fn) integrand of virtual work equation
C###  Routine: D_AG      (fn) deriv of integrand in virtual work eqn
C###  Routine: D_ENERGY  calc 2nd derivs of strain energy function
C###  Routine: D_PFRE_NE deriv of pres constr resid wrt mat/geom pars
C###  Routine: D_PFRF    deriv of press loading resids
C###  Routine: D_ZERE50  deriv of elem resids wrt mat/geom params
C###  Routine: D_ZGTG53  deriv of stress cpts (3D) wrt mat/geom vars
C###  Routine: D_ZGTG54  deriv of stress cpts (mem) wrt mat/geom vars
C###  Routine: E_PRES    calc equilibrium press for polynom strn energy
C###  Routine: ENERGY    calc derivatives of strain energy function
C###  Routine: INITVC    puts VC array into VC_init array
C###  Routine: IPANA5    input for analytic functions
C###  Routine: IPBAS5    input for basis functions
C###  Routine: IPINI5    input i.c.s and boundary conditions
C###  Routine: IPMAT5    input for material parameters
C###  Routine: IPMAT5_SINGLE  inputs a single material parameter
C###  Routine: OPC50     Gauss pt stress & strain fields in an elem
C###  Routine: OPEG50    o/p strain fields at a Gauss point
C###  Routine: OPEG55    o/p strain cpts wrt fibr/ref coords (string)
C###  Routine: OPINI5    output initial and boundary data
C###  Routine: OPMAT5    output passive muscle model parameters
C###  Routine: OPTG50    output stress fields at a Gauss pt
C###  Routine: PFRE_NE   calc stress constr resid terms for elem vars
C###  Routine: PFRE_NP   calc stress constr resid terms for node vars
C###  Routine: PFRF      calc pressure loading residual terms
C###  Routine: UPPRESSVOROIDEAL update press in voro from ideal gas law
C###  Routine: UPPRESSVOROPVCURVE update press in voro from pv curve
C###  Routine: USER51    user strain energy fn of principal invars
C###  Routine: USER52    user strain energy fn of principal extns
C###  Routine: USER53    user strain energy fn of fibre/trans strains
C###  Routine: USER53A   calc derivatives of strain energy fn
C###  Routine: ZEOP      output diagnostic strain & stress fields
C###  Routine: ZERE50    calc element residual
C###  Routine: CONTACT_RESIDUAL   adds contact residuals to global residual vector
C###  Routine: CONTACT_STIFFNESS  adds contact stiffness to global stiffness matrix
C###  Routine: ZERE55    press-vol eles for ventric cavity loading
C###  Routine: ZETX50    calc stress cpts wrt fibre/ref coords
C###  Routine: ZGTG51    eval stress tensor for pl stress problems
C###  Routine: ZGTG52    eval stress tensor for pl strain problems
C###  Routine: ZGTG53    eval stress tensor for 3D problems
C###  Routine: ZGTG54    eval cpts of stress tensor for membrane
C###  Routine: ZGTG55    eval cpts of stress tensor for string
C###  Routine: ZGTG5A    eval actively developed fibre stress


C#### Variable: IL_fluid_conductivity
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    IL_fluid_conductivity is IL number for through wall fluid
C###    conductivity.

C#### Variable: ILPIN(il,nr,nx)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    <HTML><PRE>
C###    ILPIN(il,nr,nx), il=1,ILTIN(nr,nx) records the whether input
C###    material parameter il is:
C###      1 constant spatially
C###      2 Piecewise constant (varies by elements)
C###      3 Piecewise linear (varies by nodes)
C###    This is then transferred into the ILP array which stores
C###    the corresponding info for the constitutive law parameters.
C###    </PRE></HTML>
C###  See-Also: ILP

C#### Variable: IL_sarcomere
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    IL_sarcomere is IL number for stress-free SL distribution.

C#### Variable: IL_thickness
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    IL_thickness  is IL number for membrane thickness.

C#### Variable: ILTIN(nr,nx)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    ILTIN(nr,nx) stores the total number of input material
C###    parameters. This is then transferred into the ILT array
C###    which stores the total number of constitutive law parameters.
C###  See-Also: ILT

C#### Variable: IL_time_delay
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    IL_time_delay is IL number for time-delay variable.

C#### Variable: KTYP51(nr)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    KTYP51(nr) is plane stress/plane strain/3D/membrane/
C###    string/shell.

C#### Variable: KTYP52(nr)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    KTYP52(nr) is compressible/incomp/incomp+fliud/comp+fluid/incomp+inext
C###    with fluid movement.

C#### Variable: KTYP53(nr)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    KTYP53(nr) is theta coordinates/Nu coordinates/Nu coordinates
C###    + active.

C#### Variable: KTYP54(nr)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    KTYP54(nr) is hyperelastic/viscosity/creep.

C#### Variable: KTYP55(nr)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    KTYP55(nr) is principal strain invariants/extension ratios/
C###    fibre strains.

C#### Variable: KTYP55a(nr)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    KTYP55a(nr) is Exponential law selection  Bloomgarden/Holmes

C#### Variable: KTYP56(nr)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    KTYP56(nr) is constant/poly/specific func/specific func/
C###    specific func/user defined.

C#### Variable: KTYP57(nr)
C###  Type: INTEGER
C###  Set_up: IPINI5
C###  Description:
C###    KTYP57(nr) is no pressure bc/pressure increments entered

C#### Variable: KTYP58(nr)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    KTYP58(nr) is geometric coordinates/displacements for dependent
C###    variable.

C#### Variable: KTYP59(nr)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    KTYP59(nr) is elastance/Hill-type/fading-memory formulation.

C#### Variable: KTYP5A(nr)
C###  Type: INTEGER
C###  Set_up: IPINI5
C###  Description:
C###    KTYP5A(nr) uses Xi3 face constraints which enforce
C###    incompressibility/continuous normal Cauchy stress to
C###    determine free hydrostatic pressure variables.

C#### Variable: KTYP5B
C###  Type: INTEGER
C###  Set_up: IPINI5
C###  Description:
C###    KTYP5B is the material parameter number of time-delay variable.

C#### Variable: KTYP5C
C###  Type: INTEGER
C###  Set_up: IPINI5
C###  Description:
C###    KTYP5C is the material parameter number of through wall fluid
C###    conductivity.

C#### Variable: KTYP5F(nr)
C###  Type: INTEGER
C###  Set_up: IPINI5
C###  Description:
C###    KTYP5F is 1 if the shear parameters in the pole-zero law
C###    are determined from the fibre distribution model, and
C###    2 if they are input.

C#### Variable: NHT50(ktyp52(nr),ktyp51(nr))
C###  Type: INTEGER
C###  Set_up: IPBAS5
C###  Description:
C###    NHT50(ktyp52(nr),ktyp51(nr)) is number of global variables
C###    reqd for problem type specified by ktyp51 and ktyp52 as follows:
C###    Plain stress compress:NHT50(1,1)=NJT;incompress: NHT50(2,1)=NJT
C###    Plain strain      "           2   2       "              2   2
C###    3-dimensional     "           3   3       "              3   4
C###    Membrane          "           4 NJT+1     "              4  NJT
C###    String            "           5 NJT+1     "              5  NJT
C###    Shell                         6   3       "              6   3

C#### Variable: NMBIN(nm,nr,nx)
C###  Type: INTEGER
C###  Set_up: IPMAT5
C###  Description:
C###    NMBIN(nm,nr,nx) stores the basis number for piecewise linear
C###    material parameter variation of entered material parameters.
C###    This is then transferred into the NMB array which stores
C###    the same info for constitutive law parameters.
C###  See-Also: NMB

C#### Variable: AXL(i,j)
C###  Type: REAL*8
C###  Set_up: ZERE50
C###  Description:
C###    AXL(i,j) is a deformed metric tensor component wrt undeformed
C###    Nu-coordinates (KTYP53(nr)>1) or wrt undeformed theta
C###    coordinates (KTYP53(nr)=1).

C#### Variable: AXU(i,j)
C###  Type: REAL*8
C###  Set_up: ZERE50
C###  Description:
C###    AXU(i,j) is a deformed metric tensor component wrt undeformed
C###    Nu-coordinates (KTYP53(nr)>1) or wrt undeformed theta
C###    coordinates (KTYP53(nr)=1).

C#### Variable: AZL(i,j)
C###  Type: REAL*8
C###  Set_up: ZERE50
C###  Description:
C###    AZL(i,j) is a deformed metric tensor component wrt undeformed
C###    Nu-coordinates (KTYP53(nr)>1) or wrt undeformed theta
C###    coordinates (KTYP53(nr)=1).

C#### Variable: AZU(i,j)
C###  Type: REAL*8
C###  Set_up: ZERE50
C###  Description:
C###    AZU(i,j) is a deformed metric tensor component wrt deformed
C###    Nu-coordinates (KTYP53(nr)>1) or wrt deformed theta
C###    coordinates (KTYP53(nr)=1).

C#### Variable: CIN(il,0:ng,nep)
C###  Type: REAL*8
C###  Set_up: IPMAT5
C###  Description:
C###    CIN(il,0:ng,nep) are the input material parameters
C###    (which may be constant or element/nodally varying CIN(il,0,nep)
C###    or gauss point varying CIN(il,ng,nep))
C###    from which constitutive parameters are calculated or copied.

C#### Variable: DXIX(ni,nj)
C###  Type: REAL*8
C###  Set_up: ZERE50
C###  Description:
C###    DXIX(ni,nj) are partial derivatives of Xi wrt Xj coords
C###    (if KTYP53(nr)=1) or wrt Nu coords (if KTYP53(nr)>1).

C#### Variable: D_RE(ns,nh,nop)
C###  Type: REAL*8
C###  Set_up: D_ZERE50
C###  Description:
C###    D_RE(ns,nh,nop) is the derivative of element residuals
C###    wrt material/geometric parameters.

C#### Variable: D_RI3(nhs)
C###  Type: REAL*8
C###  Set_up: D_ZERE50
C###  Description:
C###    D_RI3(nhs) is the derivative of the third strain invariant
C###    wrt material/geometric parameters.

C#### Variable: D_TG(i,j,nhs)
C###  Type: REAL*8
C###  Set_up: D_ZERE50
C###  Description:
C###    D_TG(i,j,nhs) is the derivative of the 2nd Piola-Kirchoff
C###    stress tensor wrt material/geometric parameters.

C#### Variable: D_ZG(nhx,nu,nhs)
C###  Type: REAL*8
C###  Set_up: D_ZERE50
C###  Description:
C###    D_ZG(nhx,nu,nhs) is the derivative of the deformed coordinates
C###    wrt material/geometric parameters.

C#### Variable: GXL(i,j)
C###  Type: REAL*8
C###  Set_up: ZERE50
C###  Description:
C###    GXL(i,j) is an undeformed metric tensor component wrt
C###    Xi-coordinates.

C#### Variable: GXU(i,j)
C###  Type: REAL*8
C###  Set_up: ZERE50
C###  Description:
C###    GXU(i,j) is an undeformed metric tensor component wrt
C###    Xi-coordinates.

C#### Variable: ZG1(nhx,nu)
C###  Type: REAL*8
C###  Set_up: D_ZERE50
C###  Description:
C###    ZG1(nhx,nu) are the Gauss point values (nu=1) and derivatives
C###    (nu=2..) of dependent variables.


C#### Module: FE60
C###  Description:
C###    General fluid dynamics and Delaunay/Voronoi routines

C###  Routine: ANG3P          (fn) Calculates the angle between 3 points
C###  Routine: DEGREE         (fn) Genmesh function
C###  Routine: DETERMINANT    (fn) Determinant of a 3x3 matrix A
C###  Routine: ENCASPULATED   (fn) Determines whether point inside tet
C###  Routine: ENTRIES_MATCH  (fn) Finds whether all of X are in B
C###  Routine: FOUND_LIST     (fn) Constructs a search list for tets
C###  Routine: IND            (fn) Genmesh function.
C###  Routine: INTERNAL       (fn) Tests whether tet is internal
C###  Routine: INTERSECT      (fn) Checks whether 2 faces share nodes
C###  Routine: SEARCH_NEXT    (fn) Finds closest tet face to a node
C###  Routine: SPACE          (fn) FInds whether space available for IN
C###  Routine: ADD_NODE       adds node to seed points, genmesh
C###  Routine: ADD_TETRAHEDRAL adds a tet to triangulation set. Genmesh
C###  Routine: ASSEMBLE_VORO  Genmesh routine.
C###  Routine: ASSERT_VORO    Error handling for a Genmesh routine.
C###  Routine: CALCDIVU       calculates velocity divergence
C###  Routine: CALC_NENP_VORO calculates nenp for voronoi
C###  Routine: CALC_NUMCAP calcs anatomical estimate of # of capillaries
C###  Routine: CALCULATE_FSPACE calculates the 'f-space'. Genmesh.
C###  Routine: CALC_VORO      calculates voronoi diagram
C###  Routine: CALC_VORO_CAPS creates capillary/voronoi mesh
C###  Routine: CALC_VORO_ELEMENTS converts Voronoi to node and elem
C###  Routine: CALC_VORO_FACE2D calculates voronoi face area
C###  Routine: CALC_VORO_FACE3D calculates voronoi face area
C###  Routine: CALC_VORO_LUNG calculates embedded voronoi-lung mesh
C###  Routine: CHECK_INTEGRITY Checks Delaunay mesh integrity
C###  Routine: CHECKMSH       checks delaunay triangulation
C###  Routine: CHECK_TETRAHEDRAL Genmesh routine.
C###  Routine: CHECK_TETRAHEDRAL2 Genmesh routine
C###  Routine: CIRCUM         computes circumcentre
C###  Routine: CIRCUM2     computes circumcentre/or centre of triangle
C###  Routine: CONSTRUCT      constructs intercept of a plane
C###  Routine: CORRVELO       corrects fluid velocities
C###  Routine: CREATE_CAPILLARY creates capillary mesh from Voronoi
C###  Routine: DELAUNAY       generates a 3D Delaunay mesh
C###  Routine: DELETE_TETRAHEDRAL deletes a tet from data set. Genmesh
C###  Routine: EDGEFLIP       flips diagonal of 2 adjoined triangles
C###  Routine: FIND_TETRAHEDRAL_TAIL fins tail of tets. Genmesh
C###  Routine: FIX_TETRAHEDRAL Genmesh routine.
C###  Routine: FLIP2D         flips 2d  delaunay triangulation
C###  Routine: FLIP_EDGE_TO_FACE flips edge to face. Genmesh
C###  Routine: FORMPPEM       forms pressure poisson equation
C###  Routine: FORMPPEM_OUTLET calculates PPE matrix for outlet
C###  Routine: FSOLV3V        solves 3x3 system 'quickly'
C###  Routine: GET_CIRCUMDATA determines tet circumdata. Genmesh
C###  Routine: GET_SPACE      Genmesh routine
C###  Routine: INSERTION      Sorts an array into ascending order.
C###  Routine: INTERNAL_NODES Automatically generates IN nodes. Genmesh
C###  Routine: IPBAS6         inputs basis functions for FE60
C###  Routine: DELAUNAY_NE    writes tet elements to ipelem format
C###  Routine: IPINI6         inputs initial conditions for FE60
C###  Routine: DELAUNAY_NP    writes tet nodes to ipnode format
C###  Routine: MAKE_TETRAHEDRAL creates initial enclsing tet. Genmesh
C###  Routine: MESHFLUX       computes flux from moving mesh
C###  Routine: MVINTNOD       moves internal nodes
C###  Routine: NENXI_VORO     calculates nxi for Voronoi.
C###  Routine: NETTFLUX       computes nett momentum of fluid
C###  Routine: NFNXF          sets up NXF from face information.
C###  Routine: NPNFNP         sets up NFNP and XPNP
C###  Routine: OPINI6         outputs initial/boundary conditions
C###  Routine: OUTLET_MEM     calc memory needed for duct flow obc
C###  Routine: PEEL_AWAY      peels away boundary and initial tets.
C###  Routine: PEEL_TETRAHEDRAL removes peeled tets from list. Genmesh
C###  Routine: RECONNECT      reconnects delaunay mesh by 'flipping'
C###  Routine: READ_GEOMETRY  reads boundary geometry. Genmesh
C###  Routine: REJECT_TRIANGLES retriangulates for boundary.
C###  Routine: RHIECHOW       calculates rhie chow fluxes
C###  Routine: RHIECHOW_IBUFFER buffer to RHIECHOW for parallelisation
C###  Routine: SOLVPPEQ       solves the pressure poisson equation
C###  Routine: SPHERE_CALC_CIRCUM calcs ccentre of 3 points on sphere
C###  Routine: SPHERE_CALC_DIST calcs dist |.| 2 points on sphere
C###  Routine: SPHERE_DELAUNAY calcs Delaunay tri of vertices on sphere
C###  Routine: SPHERE_POINTS  generates points on surface of sphere
C###  Routine: SPHERE_PROJECT projects point onto surface of unit sphere
C###  Routine: SPRSVORO       calculates the Voronoi sparsity structure
C###  Routine: TENTVELO       calculates tentative velocities
C###  Routine: TIMESTEP       calculates max allowable timestep
C###  Routine: TRI_NE_NP      puts Delaunay tri into node and elem
C###  Routine: TWO_THREE      Voronoi routine.
C###  Routine: VORO_XPXE      puts XP info into XE for geometry only
C###  Routine: XINORM         gets point on face nrml to internal point
C###  Routine: XINORM1D       gets XI coords of point on line.
C###  Routine: XNORMXI        gets X coords and nrml of point on face

C#### Variable: nfvl
C###  Type: INTEGER
C###  Set-up: CALC_VORO
C###  Description:
C###    nfvl is the local Voronoi cell face number

C#### Variable: nvc
C###  Type: INTEGER
C###  Set-up: CALC_VORO
C###  Description:
C###    nvc is the global Voronoi cell loop variable

C#### Variable: nfv
C###  Type: INTEGER
C###  Set-up: CALC_VORO
C###  Description:
C###    nfv is the global Voronoi face number

C#### Variable: NENFVC(0:nfvcl,nfl)
C###  Type: INTEGER
C###  Set-up: CALC_VORO_FACE3D
C###  Description:
C###    <HTML>
C###    NENFVC(0,nfl) is the number of vertices in face nfl.
C###    NENFVC(1..NFVCM,nfl) are the elements (ne) that contribute
C###      to the vertices of face nfl.
C###    </HTML>

CALC_VORO_FACE3D

C#### Variable: NEPZ(ne)
C###  Type: INTEGER
C###  Set-up: CALC_VORO_LUNG
C###  Description:
C###    <HTML>
C###    NEPZ is the global host element number for vertex ne of a
C###         Voronoi cell-Lumped parameter acinus.
C###    </HTML>

C#### Variable: NFVC(2,0:nfvl,nvc)
C###  Type: INTEGER
C###  Set-up: CALC_VORO
C###  Description:
C###    <HTML>
C###    NFVC is the connectivity array for the Voronoi cells
C###    <PRE>
C###    NFVC(1,0,nvc)    = Number of nodes connected to Voronoi cell nvc
C###    NFVC(2,0,nvc)    = Currently set to same as above
C###    NFVC(1,nfvl,nvc) = adjacent nonode number of local face nfvl of
C###                       Voronoi cell nvc (nfvl=1..NFVC(1,0,nvc))
C###    NFVC(2,nfvl,nvc) = global face number of local face nfvl of
C###                       Voronoi cell nvc
C###    </PRE>
C###    </HTML>

C#### Variable: NFVCL(0:10,nfv)
C###  Type: INTEGER
C###  Set-up: CALC_VORO_LUNG
C###  Description:
C###    <HTML>
C###    NFVCL maps between Voroni faces and Voroni vertices, for
C###          embedded Voronoi-lung models.
C###    <PRE>
C###    NFVCL(0,nvc)    = Number of Voronoi vertices in face nfv
C###    NFVCL(npl,nvc)  = Voronoi vertex # for face nfv
C###    </PRE>
C###    </HTML>

C#### Variable: NVCBBM(np)
C###  Type: INTEGER
C###  Set-up: CALC_VORO_LUNG
C###  Description:
C###    <HTML>
C###    NVCBBM is the mapping array between Voronoi cells and
C###           Lumped Parameter respiratory airway models
C###    </HTML>

C#### Variable: XNFV(-(NJM+1):NJM,nfv)
C###  Type: REAL*8
C###  Set-up: CALC_VORO_FACE2D,CALC_VORO_FACE3D
C###  Description:
C###    <HTML>
C###    XNFV is the geometric description of global Voronoi face nfv
C###    <PRE>
C###    XNFV(-2..-4,nfv) = face centroid
C###    XNFV(-1,nfv)     = Area of Voronoi face nfv
C###    XNFV(0,nfv)      = Inverted distance between the two nodes that
C###                       span Voronoi face nfv
C###   XNFV(1..NJM,nfv)  = Components of unit normal of face nfv
C###    </PRE>
C###    </HTML>

C#### Variable: VC(0:nvc)
C###  Type: REAL*8
C###  Set-up: CALC_VORO
C###  Description:
C###    <HTML>
C###    VC is the volumetric description of the Voronoi cell
C###    <PRE>
C###    VC(nvc) = Volume of Voronoi cell nvc
C###    VC(0)   = Volume of Entire Voronoi mesh. (Sum of VC(1..NVCT))
C###    </PRE>
C###    </HTML>

C#### Variable: OLDVC(0:nvc)
C###  Type: REAL*8
C###  Set-up: MARCH60
C###  Description:
C###    <HTML>
C###    OLDVC is the volumetric description of the Voronoi cell at the
C###    previous time step (Used for deforming meshes)
C###    <PRE>
C###    OLDVC(nvc) = Volume of Voronoi cell nvc
C###    OLDVC(0)   = Volume of Entire Voronoi mesh.(Sum of VC(1..NVCT))
C###    </PRE>
C###    </HTML>

C#### Variable: NFVT
C###  Type: INTEGER
C###  Set-up: CALC_VORO
C###  Description:
C###    NFVT is the total number of global Voronoi faces

C#### Variable: NVCT
C###  Type: INTEGER
C###  Set-up: CALC_VORO
C###  Description:
C###    NVCT is the total number of Voronoi cells

C#### Variable: NVBT
C###  Type: INTEGER
C###  Set-up: IPMESH11
C###  Description:
C###    NVBT is the number of boundary nodes for Voronoi meshes.

C#### Variable: NVIBT
C###  Type: INTEGER
C###  Set-up: IPMESH11
C###  Description:
C###    NVIBT is the number of Internal boundary Voronoi cells.

C#### Variable: NVIT
C###  Type: INTEGER
C###  Set-up: IPMESH11
C###  Description:
C###    NVIT is the number of Internal Voronoi cells.

C#### Variable: ZA(1,1..4,1,ne)
C###  Type: REAL*8
C###  Set-up: UPDELA
C###  Description:
C###    <HTML>
C###    For Voronoi type problems, ZA is different from the usual
C###    set up. For the triangles/tetrahedra (Dual of the Voronoi mesh),
C###    ZA contains the following:
C###    <PRE>
C###    ZA(1,1..3,1,ne) = x,y,z coordinates of the circumcentre of
C###                      the circle/sphere generated by the
C###                      triangle/tetrahedron. See Subroutine CIRCUM
C###                      as to how this is calculated.
C###    ZA(1,4,1,ne)    = the radius squared of the circle/sphere.
C###
C###    </PRE>
C###    </HTML>

C#### Variable: ZNFV(nfv)
C###  Type: REAL*8
C###  Set-up: MARCH60
C###  Description:
C###    <HTML>
C###    ZNFV is the array that stores Voronoi face based flux values.
C###    <PRE>
C###    ZNFV(nfv) = solenoidal Rhie-Chow face flux
C###    </PRE>
C###    </HTML>

C#### Variable: ZNFVMSH(nfv)
C###  Type: REAL*8
C###  Set-up: MARCH60
C###  Description:
C###    <HTML>
C###    ZNFVMSH is the array that stores Voronoi face based mesh flux
C###    values
C###    <PRE>
C###    ZNFV(nfv) = mesh flux = rate of volume swept by Voronoi face nfv
C###    </PRE>
C###    </HTML>

C#### Variable: NXI(0,0:nei,ne)
C###  Type: INTEGER
C###  Set-up: NENXI_VORO
C###  Description:
C###    <HTML>
C###    NXI_VORO has a slightly different meaning for Voronoi problems
C###    <PRE>
C###    NXI(0,0,ne)  = number of adjoining elements to element ne
C###    NXI(0,nn,ne) = (nn = 1..nei) Elements adjoined to element ne.
C###                   For each local node in a tetrahedra, this
C###                   local node will be opposite a face of the
C###                   tetrahedra. This face will be adjoined to
C###                   another element, or will be adjoined to
C###                   nothing, ie the boundary. Thus the local node nn
C###                   is "opposite" an element.
C###    </PRE>
C###    </HTML>

C#### Variable: NVCB(-1:3,nvcb)
C###  Type: INTEGER
C###  Set-up: IPINI6
C###  Description:
C###    <HTML>
C###    NVCB contains the boundary node and boundary condition types
C###         for fe60 Voronoi fluid dynamics problems
C###    <PRE>
C###    NVCB(-1,nvcb) = Boundary condition type for Voronoi boundary
C###                    node nvcb
C###    NVCB(0,nvcb)  = Number of internal boundary nodes adjoined to
C###                    boundary node nvcb
C###    NVCB(1..3,nvcb)=Internal boundary nodes adjoined to boundary
C###                    node
C###    </PRE>
C###    </HTML>

C#### Variable: NODENVC
C###  Type: INTEGER
C###  Set-up: IPINI6
C###  Description:
C###    <HTML>
C###    NODENVC is used to return the nonode number for Voronoi cell
C###    nvc for the Voronoi region
C###    <PRE>
C###    NODENVC(nvc) = nonode for Voronoi region of Voronoi node nvc
C###    </PRE>
C###    </HTML>

C#### Variable: NODENVCB
C###  Type: INTEGER
C###  Set-up: IPINI6
C###  Description:
C###    <HTML>
C###    NODENVC is used to return the nonode number for Voronoi boundary
C###    node nvcb for the Voronoi region
C###    <PRE>
C###    NODENVCB(nvcb) = nonode for Voronoi region of Voronoi node nvcb
C###    </PRE>
C###    </HTML>

C#### Variable: NVCNODE
C###  Type: INTEGER
C###  Set-up: IPINI6
C###  Description:
C###    <HTML>
C###    Returns either the nvc number (Voronoi cell) or nvcb number
C###    (Voronoi boundary node) given a nonode number in the Voronoi
C###    region.
C###    <PRE>
C###    NVCNODE(1,nonode) = Type of Voronoi node (Boundary = 1)
C###                                             (Internal Boundary = 2)
C###                                             (Internal = 3)
C###                        for a nonode number in the Voronoi region
C###    NVCNODE(2,nonode) = boundary node or Voronoi cell number for
C###                        a nonode number in the Voronoi region
C###    </PRE>
C###    </HTML>


C#### Module: FE70
C###  Description:
C###    Sail design routines.
C###  Routine: D:H3     interpolate deriv. with 1D cubic Hermite basis
C###  Routine: X:H3     interpolates with 1D cubic Hermite basis
C###  Routine: X:LENG   integrates length of cubic in x,y or x,z plane
C###  Routine: IPSAIL   sail parameter input
C###  Routine: OPSAIL   sail parameter output
C###  Routine: SCHORD   shows shape of individual chord
C###  Routine: XHERM    returns cubic Hermite polyline in global coords

C Note: GKS calls here need passing through fe14 routines.

C#### Module: FE90
C###  Description:
C###    Boundary element routines
C###  Routine: ANGLE       calculates the angle between two vectors.
C###  Routine: CALC_CURVCORRECT calculates the curvature correction.
C###  Routine: CALC_INT_SCHEME calculates the integration scheme.
C###  Routine: CHKSOL      compares solution with analytic solutions
C###  Routine: CHEBYSHEV2  *REMOVED* calculates 2nd kind Chebyshev polynomial
C###  Routine: COEFF       calcs diagonal coefficient.
C###  Routine: DIPOLE_EVALUATE evaluates soln for dipole in sphere(s)
C###  Routine: DIPOLE_SOLVE finds coeffs for dipole inside sphere(s)
C###  Routine: DIST        calcs min dist. from node to element
C###  Routine: DIST_LOOSE  calcs min dist. from node to elem, larger tol
C###  Routine: DOMSOL      calculates BE sol at predefined domain points
C###  Routine: DOMSOL2     calculates BE sol at predefined domain points
C###  Routine: DOMSOLOPTI2D 2d Laplace optimised DOMSOL
C###  Routine: DOMSOLOPTI2D 3d Laplace optimised DOMSOL
C###  Routine: ECCENTRIC_DIPOLE calculates pot due to eccentric dipole
C###  Routine: EQTYPE      identifies Green's fn and physical constants
C###  Routine: FINDK       calcs the k values for the Runge-Kutta scheme
C###  Routine: GETDIPOLE   calcs the dipole value at time t
C###  Routine: GRADDOMSOL  calcs grad(phi) at domain points
C###  Routine: GRADDOMSOLOPTI2D calcs grad(phi) at domain points
C###  Routine: GRADDOMSOLOPTI3D calcs grad(phi) at domain points
C###  Routine: IPANA9      input for analytic functions
C###  Routine: IPBAS9      input for basis functions
C###  Routine: IPINI9      inputs initial conds and b.c. for BE probs.
C###  Routine: MAGSOL      calcs magnetic solution from potentials
C###  Routine: MAGSOLOPTI  optimised verson of MAGSOL
C###  Routine: MAGSOLANALYTIC calcs magnetic analytic solution
C###  Routine: NORMAL      finds norm vector to a surface for bdry elems
C###  Routine: OPC90       outputs solution and reactions
C###  Routine: OPINI9      outputs initial and boundary data
C###  Routine: QUADBE      identifies appropriate quadrature scheme
C###  Routine: RUNGE_KUTTA 4th order R-K scheme to track current lines
C###  Routine: TELLES      Calculates Telles adaptive int. params.
C###  Routine: XEGKGQ      calcs bdry element stiffness matrices for ne
C###  Routine: XEGKGQ_3DL  efficient 3d Laplace eqtn XEGKGQ
C###  Routine: XEGKGQ_3DL_BUFFER Buffered efficient 3d Laplace eqtn XEGKGQ
C###  Routine: XEPGKGQ     calcs bdry element integrals
C###  Routine: XEPGKGQDOMSOL calcs bdry elem integrals for domain sols
C###  Routine: XEPGKMAGSOL calcs bdry elem integrls for magnetic H sols
C###  Routine: XEPGKMAGSOL2 calcs bdry elem integrls for magnetic A sols
C###  Routine: XEPGKGQ_3DL efficient 3d Laplace XEPKGKGQ routine
C###  Routine: XEPGKGQHYP  calcs derivative BEM integrals
C###  Routine: XEPGKGQHYP_3DL efficient 3d Laplace eqtn XEPGKGQHYP
C###  Routine: XEPGKGQHYPS calcs singular derivative BEM integrals
C###  Routine: XEPGKGQHYPS_3DL efficient 3d Laplace eqtn XEPGKGQHYPS
C###  Routine: XPGKGQ      calcs bdry element stiffness matrices for np
C###  Routine: XPGKGQ_3DL  efficient 3d Laplace eqtn XPGKGQ
C###  Routine: XPGKGQ_3DL_BUFFER  Buffered efficient 3d Laplace eqtn XPGKGQ
C###  Routine: XPGD        calcs domain integrals and/or source terms
C###  Routine: XPXNO       transfers XP to XNO for derivative BIE
C###  Routine: UPBE        *** ARCHIVED ***

C**** General purpose boundary element routines

C#### Variable: CURVCORRECT(i,j,nn,ne)
C###  Type: REAL*8
C###  Set-up: CALC_CURVCORRECT
C###  Description:
C###    CURVCORRECT(i,j,nn,ne) is the curvature correction applied to
C###    du/ds_i when the normal derivative is differentiated in the
C###    s_j direction at local node nn of element ne.

C#### Variable: DIPOLE_COEFF(2*nspheres,max_mn_index)
C###  Type: REAL*8
C###  Set_up: DIPOLE_SOLVE
C###  Description:
C###    DIPOLE_COEFF(2*nspheres,max_mn_index) contains the coefficients
C###    in the series expansion of the solution for a dipole inside
C###    nspheres spheres. The first index is the (A,B,C,D,...)
C###    coefficients. The second index is to get the coefficients
C###    for the (m,n)th terms in the series expansion and the index
C###    is calculated from n*(max_m+1)+m.

C#### Variable: IGREN(nr)
C###  Type: INTEGER
C###  Set-up: EQTYPE
C###  Description:
C###    <HTML>
C###    IGREN(nr) is the type of Green's function used in the region nr.
C###    <PRE>
C###    IGREN(nr) = 1    Laplace equation, 2D cartesian.
C###              = 2    Laplace equation, 3D cartesian.
C###              = 3    Helmholtz equation, 2D cartesian.
C###              = 4    Helmholtz equation, 3D cartesian.
C###              = 5    Modified Helmholtz (Yukawa) equation, 2D cart.
C###              = 6    Modified Helmholtz (Yukawa) equation, 3D cart.
C###              = 7    Generalised Laplace equation, 2D cartesian.
C###              = 8    Generalised Laplace equation, 3D cartesian.
C###              = 9    Poisson equation, 2D cartesian, special source
C###              = 10   Poisson equation, 3D cartesian, special source
C###              = 11   Helmholtz equation, 3D cyl. sym.
C###                     (about x axis), complex
C###              = 12   Modified Helmholtz(Yukawa) eq., 3D cyl. sym.
C###                     about z axis
C###              = 13   2D linear elasticity - plane strain
C###              = 14   2D linear elasticity - plane stress
C###              = 15   3D linear elasticity
C###    </PRE>
C###    </HTML>


C#### Variable: KTYP93(nc,nr)
C###  Type: INTEGER
C###  Set_up: IPBAS9
C###  Description:
C###    <HTML><P>
C###    KTYP93(nc,nr)=0,1 if cross derivatives/nocross derivatives are
C###    used for the dependent variables in the BE derivative equations.</P>
C###    <P>      
C###    In Jan 2000 a new question was added (for bicubic hermite
C###    BE basis functions) that asked whether cross derivatives
C###    were to be set to zero or not (previously this was asked in
C###    ipequa).  Answering yes to the new question results in
C###    no ny being set up for the cross derivatives and NKH adjusted
C###    accordingly (previously NKH was not adjusted and
C###    KTYP93 was subtracted from NKH in dependent
C###    variable calculations).  This new method was introduced so that
C###    one could apply transfer matrices constructed using bicubic
C###    hermite BE basis functions to fitted fields.</P>
C###    <P>
C###    Thus currently there are two ways to get a bicubic hermite
C###    BE solution (with no cross derivatives) - setting them to zero
C###    in the ipbase file (and KTYP93 will be set to 0 and no further
C###    question should be asked in the ipequa stage) or not setting
C###    them to zero in the ipbas# routine and setting them to zero at
C###    the ipequa stage (ie the old way - this will set KTYP93 to 1).
C###    </P><P>
C###    Results obtained from each should be the same, but the
C###    first method must be used if a transfer matrix is going to be
C###    constructed and applied to fitted bicubic hermite BE fields.
C###    The first method is preferred and KTYP93 is obsolescent.</P>
C###    </HTML>      

C#### Variable: NPB(0:np,5)
C###  Type: INTEGER
C###  Set_up: CALC_INIT_SCHEME
C###  Description:
C###    NPB(0:np,5) stores the list of nodes for each type of
C###    integration scheme. NPB(0,1..5) is the number of nodes for the
C###    split element/adaptive/high/medium/low integration schemes.
C###    NPB(1..NPB(0,i),i) is the list of nodes in the ith integration
C###    scheme.

C#### Variable: DRDN(ng)
C###  Type: REAL*8
C###  Set_up: XEGKGQ
C###  Description:
C###    DRDN(ng) is the derivative of RAD in the normal direction.

C#### Variable: DRDNO(ng,nk)
C###  Type: REAL*8
C###  Set_up: XEPGKGQHYP,XEPGKGQHYP_3DL,XEPGKGQHYPS,XEPGKGQHYPS_3DL
C###  Description:
C###    DRDNO(ng,nk) is the derivative of RAD at the Gauss point ng
C###    in the differentiation direction nk of the boundary integral
C###    equation.

C#### Variable: DIPOLE_F
C###  Type: REAL*8
C###  Set_up: DIPOLE_SOLVE
C###  Description:
C###    DIPOLE_F is the eccentricity of the dipole in the eccentric
C###    spheres, multiple shell analytic solution.

C#### Variable: DIPOLE_STRENGTH(3)
C###  Type: REAL*8
C###  Set_up: DIPOLE_SOLVE
C###  Description:
C###    DIPOLE_STRENGTH contains the components of the dipole, as if
C###    the dipole was located on the z axis.

C#### Variable: DIPOLE_AXIS
C###  Type: INTEGER
C###  Set_up: DIPOLE_SOLVE
C###  Description:
C###    DIPOLE_AXIS is the axis that the dipole is eccentric on.

C#### Variable: RAD(ng)
C###  Type: REAL*8
C###  Set_up: XEGKGQ
C###  Description:
C###    RAD(ng) is the distance from the Gauss point to the current
C###    node point under consideration (BEM).

C#### Variable: RD(ng)
C###  Type: REAL*8
C###  Set_up: XPGKGQ
C###  Description:
C###    RD(ng) is the correction at the Gauss point ng to the element
C###    integrals for symmetric problems.

C#### Variable: XG1(nj,nu,ng)
C###  Type: REAL*8
C###  Set_up:
C###  Description:
C###    XG1(nj,nu,ng) stores XG for several Gauss points.

C#### Variable: XN(nj,ng)
C###  Type: REAL*8
C###  Set_up: XEGKGQ
C###  Description:
C###    XN(nj,ng) is the njth coordinate of the normal vector at the
C###    Gauss point ng.

C#### Variable: XN_GRAD(nj,ng)
C###  Type: REAL*8
C###  Set_up: XPGKGQ
C###  Description:
C###    XN_GRAD is a BEM array and is the dot product of XN with the
C###    gradients of the various arclengths (there are nje-1 arclengths
C###    and the outward normal (stored in the nj=nje position).

C#### Variable: XR(nj,ng)
C###  Type: REAL*8
C###  Set_up: XEGKGQ
C###  Description:
C###    XR(nj,ng) is the njth coordinate of the vector from the Gauss
C###    point ng to the singularity point.

C#### Variable: XR_GRAD(nj,ng)
C###  Type: REAL*8
C###  Set_up: XPGKGQ
C###  Description:
C###    XR_GRAD is a BEM array and is the dot product of XR with the
C###    gradients of the various arclengths (there are nje-1 arclengths
C###    and the outward normal (stored in the nj=nje position).

C#### Variable: YD(nh)
C###  Type: REAL*8
C###  Set_up: DOMSOL
C###  Description:
C###   Contains the domain solution.

C#### Variable: ZF(ns,nh)
C###  Type: REAL*8
C###  Set_up: DOMSOL
C###  Description:
C###    ZF(ns,nh) contains deformed face information.


C#### Module: FECELLML
C###  Description:
C###    Routines which define the interface to the CellML library.

C###  Routine: CELLML_INITIALISE Must be called to initialise CellML processor
C###  Routine: CELLML_TERMINATE Must be called when finished processing CellML
C###  Routine: CELLML_PARSE Parses a given CellML file

C#### Module: FEFUTILS
C###  Description:
C###    Fortran Utility routines.
C###  Routine: FSTRINGLEN    Determines the length of a fortran string
C###  Routine: CREATE_CSTRING allocates a C string from a fortran string
C###  Routine: DESTROY_CSTRING deallocates memory from CREATE_CSTRING
C###  Routine: C2FSTRING     Converts a c string to a fortran string
C###  Routine: F2CSTRING     Converts a fortran string to a c string

C#### Module: fegx.f
C###  Description:
C###    Routines for interfacing with the gx graphics library.
C###  Routine: CHECK_GX_OPEN    opens GX if not already open


C#### Module: FEINTERPRETER
C###  Description:
C###    Routines called by the command interpreter.

C###  Routine: EXECUTE_COMMAND               execute interpreted command
C###  Routine: INTERPRET_COMMAND_LINE        interpret and execute command
C###  Routine: SET_USER_CHARACTER            set a character interp variable
C###  Routine: SET_USER_DOUBLE               set a character interp variable
C###  Routine: SET_USER_INTEGER              set a character interp variable


C#### Module: fesolver.f
C###  Description:
C###    Routines for interfacing with the libsolver linear solver library.
C###  Routine: SOLVE_SYSTEM    Solves a linear system of equations
C###  Routine: IPSOLU_SOLVER   Inputs solution parameters
C###  Routine: OPSOLU_SOLVER   Outputs solution parameters


C#### Module: FEUSER
C###  Description:
C###    Spreadsheet fns and user-defined subroutines

C###  Routine: BR        *REMOVED* called by ODE integrator ADAMS
C###  Routine: BeelerReuter *** ARCHIVED ***
C###  Routine: RBR       *** ARCHIVED ***
C###  Routine: CUBIC     *** ARCHIVED ***
C###  Routine: DUFFING   *** ARCHIVED ***
C###  Routine: HH        *** ARCHIVED ***
C###  Routine: L2SYSTEM  *** ARCHIVED ***
C###  Routine: LINEAR    *** ARCHIVED ***
C###  Routine: LORENZ    *** ARCHIVED ***
C###  Routine: NAGUMO    *** ARCHIVED ***
C###  Routine: POINCARE1 *** ARCHIVED ***
C###  Routine: POINCARE2 *** ARCHIVED ***
C###  Routine: POINCARE3 *** ARCHIVED ***
C###  Routine: ROSSLER   *** ARCHIVED ***
C###  Routine: VANDERPOL *** ARCHIVED ***
C###  Routine: USER_9    *** ARCHIVED ***
C###  Routine: USER_10   Valve operation for fe60 problems
C###  Routine: USER_11   Moves boundaries for fe60 problems

C#### Module: FEUSER_CELL
C###  Description:
C###    User defined cellular models
C###  Routine: USER_CELL1     User defined model 1
C###  Routine: USER_CELL2     User defined model 2
C###  Routine: USER_CELL3     User defined model 3
C###  Routine: USER_CELL4     User defined model 4
C###  Routine: USER_CELL5     User defined model 5

C**** CMISS Module FD16: Dummy Coronary angiography routines

C**** CMISS Module FD18: Dummy Bead analysis routines

C**** GX Graphics Library dummy routines

C**** CMISS Module FZLAPACK.F: Dummy LAPACK routines

C**** MINOS dummy routines

C**** CMISS Module FZSOCKET.F: Dummy SOCKET routines

C**** CMISS Module fzsolver.f: Dummy routines for the solver library interface

C       SUBROUTINE ALLOC_SOLVER()
C       CALL FLAG_ERROR(0,'Link with Solver library')
C       RETURN 1
C       END

C**** CMISS Module fzxblas.f: Dummy routines for the xblas library

