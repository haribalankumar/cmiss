function umfpack_report (Control, Info)
%UMFPACK_REPORT
%
%       umfpack_report (Control, Info) ;
%
%   Prints the current Control settings for umfpack, and the statistical
%   information returned by umfpack in the Info array.  If Control is
%   an empty matrix, then the default control settings are printed.
%
%   Control is 20-by-1, and Info is 90-by-1.  Not all entries are used.
%
%   Alternative usages:
%
%       umfpack_report ([ ], Info) ;    print the default control parameters
%                                       and the Info array.
%       umfpack_report (Control) ;      print the control parameters only.
%       umfpack_report ;                print the default control parameters
%                                       only.
%
%   UMFPACK Version 4.0 (Apr 11, 2002).  Copyright (c) 2002 by Timothy A.
%   Davis.  All Rights Reserved.  Type "help umfpack_details" for License.
%
%   See also umfpack, umfpack_make, umfpack_details,
%   umfpack_demo, and umfpack_simple.

%   The contents of Control and Info are defined in umfpack.h

if (nargin < 1)
    Control = umfpack ;
end
if (isempty (Control))
    Control = umfpack ;
end

%-------------------------------------------------------------------------------
% control settings
%-------------------------------------------------------------------------------

fprintf ('\numfpack version 4.0:  Control settings:\n') ;
fprintf ('    Control (1):  print level: %d\n', Control (1)) ;
drow = Control (2) ;
fprintf ('    Control (2):  dense row parameter:    %g\n', drow) ;
fprintf ('       ("dense" rows have    > max (16, (%g)*16*sqrt(n_col)) entries)\n', drow) ;
dcol = Control (3) ;
fprintf ('    Control (3):  dense column parameter: %g\n', dcol) ;
fprintf ('       ("dense" columns have > max (16, (%g)*16*sqrt(n_row)) entries)\n', dcol) ;
fprintf ('    Control (4):  pivot tolerance: %g\n', Control (4)) ;
fprintf ('    Control (5):  max block size for dense matrix kernels: %d\n', Control (5)) ;
fprintf ('    Control (6):  relaxed amalgamation parameter: %g\n', Control (6)) ;
fprintf ('    Control (7):  initial allocation ratio: %g\n', Control (7)) ;
fprintf ('    Control (8):  max iterative refinement steps: %d\n', Control (8)) ;
fprintf ('    Control (14): relax2 amalgamation parameter: %g\n', Control (14)) ;
fprintf ('    Control (15): relax3 amalgamation parameter: %g\n\n', Control (15)) ;

% compile-time options:

fprintf ('    The following options can only be changed at compile-time:\n') ;

if (Control (9) == 1)
    fprintf ('    Control (9): compiled to use the BLAS\n') ;
else
    fprintf ('    Control (9): compiled without the BLAS\n') ;
    fprintf ('        (you will not get the best possible performance)\n') ;
end

if (Control (10) == 1)
    fprintf ('    Control (10): compiled for MATLAB\n') ;
elseif (Control (10) == 2)
    fprintf ('    Control (10): compiled for MATLAB\n') ;
    fprintf ('        Uses internal utMalloc, utFree, utRealloc, utPrintf\n') ;
    fprintf ('        utDivideComplex, and utFdlibm_hypot routines.\n') ;
else
    fprintf ('    Control (10): not compiled for MATLAB\n') ;
    fprintf ('        Uses ANSI C malloc, free, realloc, and printf\n') ;
    fprintf ('        instead of mxMalloc, mxFree, mxRealloc, and mexPrintf.\n') ;
    fprintf ('        Printing will be in terms of 0-based matrix indexing,\n') ;
    fprintf ('        not 1-based as is expected in MATLAB.  Diary output may\n') ;
    fprintf ('        not be properly recorded.\n') ;
end

if (Control (11) == 1)
    fprintf ('    Control (11): uses getrusage to get CPU time.\n') ;
else
    fprintf ('    Control (11): uses ANSI C clock to get CPU time.\n') ;
    fprintf ('        The CPU time may wrap around, type "help cputime".\n') ;
end

if (Control (12) == 1)
    fprintf ('    Control (12): compiled with debugging enabled\n') ;
    fprintf ('        This will be exceedingly slow!\n') ;
    if (Control (10) == 1)
        fprintf ('        Uses mxAssert.\n') ;
    elseif (Control (10) == 2)
        fprintf ('        Uses utAssert.\n') ;
    else
        fprintf ('        Uses ANSI C assert instead of mxAssert.\n') ;
    end
else
    fprintf ('    Control (12): compiled for normal operation\n') ;
end

%-------------------------------------------------------------------------------
% Info
%-------------------------------------------------------------------------------

if (nargin < 2)
    return
end

if (isempty (Info))
    Info = -ones (90,1) ;
    Info (1) = 0 ;
end

status = Info (1) ;
fprintf ('\numfpack status:  Info (1): %d, ', status) ;

if (status == 0)
    fprintf ('OK\n') ;
elseif (status == -1)
    fprintf ('ERROR    out of memory\n') ;
elseif (status == 1)
    fprintf ('WARNING  matrix is singular\n') ;
elseif (status == -3)
    fprintf ('ERROR    numeric LU factorization is invalid\n') ;
elseif (status == -4)
    fprintf ('ERROR    symbolic LU factorization is invalid\n') ;
elseif (status == -5)
    fprintf ('ERROR    required argument is missing\n') ;
elseif (status == -6)
    fprintf ('ERROR    n <= 0\n') ;
elseif (status <= -7 & status >= -12)
    fprintf ('ERROR    matrix A is corrupted\n') ;
elseif (status == -13)
    fprintf ('ERROR    invalid system\n') ;
elseif (status == -14)
    fprintf ('ERROR    invalid triplet form\n') ;
elseif (status == -15)
    fprintf ('ERROR    invalid permutation\n') ;
elseif (status == -16)
    fprintf ('ERROR    problem too large\n') ;
elseif (status == -911)
    fprintf ('ERROR    internal error!\n') ;
    fprintf ('Please report this error to Tim Davis (davis@cise.ufl.edu)\n') ;
else
    fprintf ('ERROR    unrecognized error.  Info array corrupted\n') ;
end

fprintf ('    Detailed statistics (a -1 means the entry has not been computed):\n') ;
fprintf ('    Info (2):  %d, number of rows of A\n', Info (2)) ;
fprintf ('    Info (17): %d, number of columns of A\n', Info (17)) ;
fprintf ('    Info (3): %d, nonzeros in A\n', Info (3)) ;
fprintf ('    Info (4): %d, Unit size, in bytes, for memory usage reported below\n', Info (4)) ;

fprintf ('\n    Computed in the symbolic analysis:\n') ;
fprintf ('    Info (5): %d, size of int (in bytes)\n', Info (5)) ;
fprintf ('    Info (6): %d, size of long (in bytes)\n', Info (6)) ;
fprintf ('    Info (7): %d, size of pointer (in bytes)\n', Info (7)) ;
fprintf ('    Info (8): %d, size of numerical entry (in bytes)\n', Info (8)) ;
fprintf ('    Info (9): %d, number of "dense" rows\n', Info (9)) ;
fprintf ('    Info (10): %d, number of "empty" rows (entries only in "dense" columns)\n', Info (10)) ;
fprintf ('    Info (11): %d, number of "dense" columns\n', Info (11)) ;
fprintf ('    Info (12): %d, number of "empty" columns (entries only in "dense" rows)\n', Info (12)) ;
fprintf ('    Info (13): %d, defragmentations during symbolic analysis\n', Info (13)) ;
fprintf ('    Info (14): %d, memory used during symbolic analysis (Units)\n', Info (14)) ;
fprintf ('    Info (15): %d, final size of symbolic factors (Units)\n', Info (15)) ;
fprintf ('    Info (16): %d, symbolic analysis time (seconds)\n', Info (16)) ;

fprintf ('    Info (18..20):  unused\n') ;

fprintf ('    Info (21): %d, estimated size of LU factors (Units)\n', Info (21)) ;
fprintf ('    Info (22): %d, estimated total peak memory usage (Units)\n', Info (22)) ;
fprintf ('    Info (23): %d, estimated factorization flop count\n', Info (23)) ;
fprintf ('    Info (24): %d, estimated number of nonzeros in L (incl. diag)\n', Info (24)) ;
fprintf ('    Info (25): %d, estimated number of nonzeros in U (incl. diag)\n', Info (25)) ;
fprintf ('    Info (26): %d, est. initial size of variable-part of LU factors (Units)\n', Info (26)) ;
fprintf ('    Info (27): %d, est. peak size of variable-part of LU factors (Units)\n', Info (27)) ;
fprintf ('    Info (28): %d, est. final size of variable-part of LU factors (Units)\n', Info (28)) ;
fprintf ('    Info (29): %d, est. max frontal matrix size (number of numerical entries)\n', Info (29)) ;

fprintf ('    Info (30..40):  unused\n') ;

fprintf ('\n    Computed in the numeric factorization:\n') ;
fprintf ('    Info (41): %d, size of LU factors (Units)\n', Info (41)) ;
fprintf ('    Info (42): %d, total peak memory usage (Units)\n', Info (42)) ;
fprintf ('    Info (43): %d, factorization flop count\n', Info (43)) ;
fprintf ('    Info (44): %d, number of nonzeros in L (incl. diag)\n', Info (44)) ;
fprintf ('    Info (45): %d, number of nonzeros in U (incl. diag)\n', Info (45)) ;
fprintf ('    Info (46): %d, initial size of variable-part of LU factors (Units)\n', Info (46)) ;
fprintf ('    Info (47): %d, peak size of variable-part of LU factors (Units)\n', Info (47)) ;
fprintf ('    Info (48): %d, final size of variable-part of LU factors (Units)\n', Info (48)) ;
fprintf ('    Info (49): %d, maximum frontal matrix size (number of numerical entries)\n', Info (49)) ;

fprintf ('    Info (50..60):  unused\n') ;

fprintf ('    Info (61): %d, defragmentations during numerical factorization\n', Info (61)) ;
fprintf ('    Info (62): %d, reallocations during numerical factorization\n', Info (62)) ;
fprintf ('    Info (63): %d, costly reallocations during numerical factorization\n', Info (63)) ;
fprintf ('    Info (64): %d, integer indices in compressed pattern of L and U\n', Info (64)) ;
fprintf ('    Info (65): %d, numerical values stored in L and U\n', Info (65)) ;
fprintf ('    Info (66): %d, numeric factorization time (seconds)\n', Info (66)) ;
fprintf ('    Info (67): %d, number of nonzeros on diagonal of U\n', Info (67)) ;
fprintf ('    Info (68): %g, reciprocal condition number estimate\n', Info (68)) ;

fprintf ('    Info (69..80):  unused\n') ;

fprintf ('\n    Computed in the solve step:\n') ;
fprintf ('    Info (81): %d, iterative refinement steps taken\n', Info (81)) ;
fprintf ('    Info (82): %d, iterative refinement steps attempted\n', Info (82)) ;
fprintf ('    Info (83): %g, omega(1), sparse-backward error estimate\n', Info (83)) ;
fprintf ('    Info (84): %g, omega(2), sparse-backward error estimate\n', Info (84)) ;
fprintf ('    Info (85): %d, solve flop count\n', Info (85)) ;
fprintf ('    Info (86): %d, solve time (seconds)\n', Info (86)) ;

fprintf ('    Info (87..90):  unused\n\n') ;

