function [out1, out2, out3, out4, out5] = umfpack (in1, in2, in3, in4, in5)
%UMFPACK
%
%   UMFPACK is a MATLAB mexFunction for solving sparse linear systems for
%   MATLAB V6.0 or later.
%
%   UMFPACK:                            |  MATLAB approximate equivalent:
%   ---------------------------------------------------------------------
%   x = umfpack (A, '\', b) ;           |  x = A \ b
%                                       |
%   x = umfpack (b, '/', A) ;           |  x = b / A
%                                       |
%   [L,U,P,Q] = umfpack (A) ;           |  Q = colamd (A) ;
%                                       |  [L,U,P] = lu (A (:,Q)) ;
%                                       |
%   [P,Q,F,C] = umfpack (A, 'symbolic') |  Q = colamd (A) ;
%                                       |  [count,h,parent,post] = ...
%                                       |  symbfact (A (:,Q), 'col') ;
%
%   A must be sparse.  It can be complex, singular, and/or rectangular.
%   A must be square for '/' or '\'.  b must be a dense real or complex
%   vector.
%
%   UMFPACK Version 4.0 (Apr 11, 2002).  Copyright (c) 2002 by Timothy A.
%   Davis.  All Rights Reserved.  Type "help umfpack_details" for License.
%
%   See also:
%   umfpack_make      to compile umfpack for use in MATLAB
%   umfpack_details   type "help umfpack_details" for more information
%   umfpack_report    prints optional control settings and statistics
%   umfpack_demo      a long demo
%   umfpack_simple    a simple demo

error ('umfpack mexFunction not found!  Use umfpack_make to compile umfpack.') ;

