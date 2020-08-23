function [L,U,P,Q,Info] = umfpack_factorize (A, Control)
%UMFPACK_FACTORIZE
%
%   [L,U,P,Q] = umfpack_factorize (A) ;
%   [L,U,P,Q,Info] = umfpack_factorize (A, Control) ;
%
%   Computes the LU factorization A(P,Q)=L*U, using umfpack.  Determines
%   the upper bound on the number of floating-point operations required to
%   factorize A and A', and choses the one with the smaller upper bound.  If A
%   is factorized, then L is unit lower-triangular; otherwise U is unit
%   upper-triangular.
%
%   UMFPACK Version 4.0 (Apr 11, 2002).  Copyright (c) 2002 by Timothy A.
%   Davis.  All Rights Reserved.  Type "help umfpack_details" for License,
%   and for details on the Control and Info arguments.
%
%   See also:  umfpack, umfpack_details, umfpack_report, umfpack_demo,
%   umfpack_simple, umfpack_make

if (nargin == 1)
    % get default parameters
    Control = umfpack ;
end

[P1,Q1,Fr1,Ch1,Info1] = umfpack (A,  'symbolic', Control) ;
[P2,Q2,Fr2,Ch2,Info2] = umfpack (A', 'symbolic', Control) ;

if (nargin > 0 & Control (1) > 0)
    fprintf ('Upper bound analysis:\n') ;
    fprintf ('Factorize A:  flop count: %d   memory usage: %d (bytes)\n', ...
        Info1 (23), Info1 (22) * Info1 (4)) ;
    fprintf ('Factorize A'': flop count: %d   memory usage: %d (bytes)\n', ...
        Info2 (23), Info2 (22) * Info2 (4)) ;
end

if (Info1 (23) < Info2 (23))
    [L,U,P,Q,Info] = umfpack (A, Control) ;
else
    [l,u,p,q,Info] = umfpack (A', Control) ;
    L = u' ;
    U = l' ;
    P = q ;
    Q = p ;
end

if (nargin > 0 & Control (1) > 0)
    fprintf ('Actual flop count: %d   actual memory usage: %d (bytes)\n', ...
        Info (43), Info (42) * Info (4)) ;
end
