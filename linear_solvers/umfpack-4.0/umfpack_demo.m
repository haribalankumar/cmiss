%UMFPACK DEMO
%
%   A demo of UMFPACK for MATLAB.
%
%   UMFPACK Version 4.0 (Apr 11, 2002).  Copyright (c) 2002 by Timothy A.
%   Davis.  All Rights Reserved.  Type "help umfpack_details" for License.
%
%   See also umfpack, umfpack_make, umfpack_details, umfpack_report,
%   and umfpack_simple.

%-------------------------------------------------------------------------------
% solve a simple system
%-------------------------------------------------------------------------------

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('Factor and solve a 5-by-5 system, Ax=b, using default parameters\n') ;
fprintf ('(except for verbose printing enabled)\n') ;

fprintf ('A, as a dense matrix:\n') ;
A = [
 2  3  0  0  0
 3  0  4  0  6
 0 -1 -3  2  0
 0  0  1  0  0
 0  4  2  0  1
]

fprintf ('A, as a sparse matrix:\n') ;
A = sparse (A)

b = [8 45 -3 3 19]'

control = umfpack ;
control (1) = 5 ;     % print everything

[xu, info] = umfpack (A, '\', b, control) ;

fprintf ('Solution to Ax=b via UMFPACK:\n') ;
x = xu

fprintf ('Solution to Ax=b via MATLAB:\n') ;
xm = A\b ;
x = xm

fprintf ('Difference between UMFPACK and MATLAB solution:\n') ;
xu - xm

%-------------------------------------------------------------------------------
% spy the results
%-------------------------------------------------------------------------------

figure (1)
clf

subplot (2,3,1)
spy (A)
title ('The matrix A') ;

subplot (2,3,2)
[Ptree, Qtree, Fr, Ch, Info] = umfpack (A, 'symbolic') ;
treeplot (Fr (:,2)') ;
title ('Supernodal column elimination tree') ;

subplot (2,3,3)
spy (A (Ptree, Qtree))
title ('A, with initial row and column order') ;

subplot (2,3,4)
[L, U, P, Q] = umfpack (A) ;
spy (A (P, Q))
title ('A, with final row/column order') ;

subplot (2,3,5)
spy (spones (L) + spones (U))
title ('UMFPACK LU factors') ;

subplot (2,3,6)
fprintf ('If you are using a version of MATLAB prior to V6.0, then the\n') ;
fprintf ('following statement (Q2 = colamd (A)) may fail.  Either download\n');
fprintf ('colamd from http://www.cise.ufl.edu/research/sparse, upgrade to\n') ;
fprintf ('MATLAB V6.0 or later, or replace the statement with\n') ;
fprintf ('Q2 = colmmd (A) ;\n') ;
Q2 = colamd (A) ;
[L2, U2, P2] = lu (A (:,Q2)) ;
spy (spones (L2) + spones (U2))
title ('MATLAB LU factors') ;

%-------------------------------------------------------------------------------
% solve A'x=b
%-------------------------------------------------------------------------------

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('Then solve A''x=b:\n') ;

[xu, info] = umfpack (b', '/', A, control) ;
xu = xu' ;

fprintf ('Solution to A''x=b via UMFPACK:\n') ;
x = xu

fprintf ('Solution to A''x=b via MATLAB:\n') ;
xm = (b'/A)' ;
x = xm

fprintf ('Difference between UMFPACK and MATLAB solution:\n') ;
xu - xm

%-------------------------------------------------------------------------------
% modify A and solve Ax=b
%-------------------------------------------------------------------------------

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('Set A (2,5) to zero and solve Ax=b\n') ;

A (2,5) = 0

[xu, info] = umfpack (A, '\', b, control) ;

fprintf ('Solution to modified Ax=b via UMFPACK:\n') ;
x = xu

fprintf ('Solution to modified Ax=b via MATLAB:\n') ;
xm = A\b ;
x = xm

fprintf ('Difference between UMFPACK and MATLAB solution:\n') ;
xu - xm

%-------------------------------------------------------------------------------
% modify all of A and solve Ax=b
%-------------------------------------------------------------------------------

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('Next, change all of the entries in A, but not the pattern.\n') ;

%     [ 2 13  0  0  0 ]      [  8 ]                  [  8.5012 ]
%     [ 2  0 23  0 39 ]      [ 45 ]                  [ -0.6925 ]
% A = [ 0  7 15 30  0 ], b = [ -3 ]. Solution is x = [  0.1667 ].
%     [ 0  0 18  0  0 ]      [  3 ]                  [ -0.0218 ]
%     [ 0 10 18  0 37 ]      [ 19 ]                  [  0.6196 ]


A = [
 2 13  0  0  0
 2  0 23  0 39
 0  7 15 30  0
 0  0 18  0  0
 0 10 18  0 37
]

A = sparse (A)

fprintf ('Compute the LU factorization of A via UMFPACK:\n') ;

[L, U, P, Q, info] = umfpack (A, control) ;

fprintf ('Here are the LU factors of A:\n') ;
L
U
P
Q

fprintf ('A (P,Q) - L*U should be zero:\n') ;
fprintf ('\nA (P,Q) - L*U =\n') ;
A (P,Q) - L*U

fprintf ('\nsolve Ax=b using the factors of A:\n') ;

xu = U \ (L \ b (P)) ;
xu (Q) = xu ;

fprintf ('    x = U \\ (L \\ b (P)) ;\n') ;
fprintf ('    x (Q) = x ;\n') ;
x = xu

fprintf ('Solve Ax=b via MATLAB:\n') ;
xm = A\b ;
x = xm

fprintf ('Difference between UMFPACK and MATLAB solution:\n') ;
xu - xm

%-------------------------------------------------------------------------------
% factor A' and then solve Ax=b using the factors of A'
%-------------------------------------------------------------------------------

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('Finally, compute B = A'', and compute the LU factorization of B\n') ;
fprintf ('Factorizing A'' can sometimes be better than factorizing A itself\n');
fprintf ('(less work and memory usage).  Solve B''x=b; the solution is the\n') ;
fprintf ('same as the solution to Ax=b for the original A.\n');

B = A'

% factorize B (P,Q) = L*U
[L, U, P, Q, info] = umfpack (B, control) ;

fprintf ('Here are the LU factors of B:\n') ;
L
U
P
Q

fprintf ('B (P,Q) - L*U should be zero:\n') ;
fprintf ('\nB (P,Q) - L*U =\n') ;
B (P,Q) - L*U

fprintf ('Solution to Ax=b via UMFPACK, using the factors of B:\n') ;
fprintf ('    x = L'' \\ (U'' \\ b (Q)) ;\n') ;
fprintf ('    x (P) = x ;\n') ;

x = L' \ (U' \ b (Q)) ;
x (P) = x ;
x

fprintf ('Solution to Ax=b via MATLAB:\n') ;
xm = A\b ;
x = xm

fprintf ('Difference between UMFPACK and MATLAB solution:\n') ;
xu - xm

