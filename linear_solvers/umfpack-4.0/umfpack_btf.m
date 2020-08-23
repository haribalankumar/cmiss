function [x, Info] = umfpack_btf (A, b, Control)
%UMFPACK_BTF
%
%  [x, Info] = umfpack_btf (A, b, Control)
%
%  solve Ax=b using umfpack and dmperm
%
%  Info is the same as for umfpack, except:
%
%       20: number of diagonal blocks (before merging 1-by-1's)
%       18: number of diagonal (or upper triangular) blocks (after
%           merging 1-by-1's)
%       19: number of off-diagonal entries
%
%   UMFPACK Version 4.0 (Apr 11, 2002).  Copyright (c) 2002 by Timothy A.
%   Davis.  All Rights Reserved.  Type "help umfpack_details" for License.
%

[m n] = size (A) ;
if (m ~= n | ~issparse (A) | issparse (b))
    error ('umfpack_btf:  A must be square and sparse; b must be full') ;
end

if (nargin < 3)
    Control = umfpack ;
end

%-------------------------------------------------------------------------------
% find the block triangular form
%-------------------------------------------------------------------------------

tdmperm = cputime ;
[p,q,r] = dmperm (A) ;
nblocks = length (r) - 1 ;
tdmperm = cputime - tdmperm ;
Info (20) = nblocks ;

%-------------------------------------------------------------------------------
% solve the system
%-------------------------------------------------------------------------------

tfactor = cputime ;
if (nblocks == 1)

    %---------------------------------------------------------------------------
    % matrix is irreducible
    %---------------------------------------------------------------------------

    [x, Info] = umfpack (A, '\', b, Control) ;
    Info (16) = Info (16) + tdmperm ;
    Info (18) = 1 ;
    Info (19) = 0 ;

else

    %---------------------------------------------------------------------------
    % A (p,q) is in block triangular form
    %---------------------------------------------------------------------------

    % Info
    terms_latest = [1 4 5 6 7 8] ;
    terms_sum    = [9 10:13 15 16 21 23:25 41 43:45 61:65 66 85 86] ;
    terms_max    = [1 26 28 29 46 47 48 49 81:84] ;
    terms_peak   = [14 22 27 42 47] ;
    terms_prev   = [15 21 28 41 48] ;
    Info = zeros (1,90) ;
    Info (2) = n ;
    Info (3) = nnz (A) ;
    Info (16) = tdmperm ;

    c = b (p) ;
    B = A (p,q) ;
    x = zeros (n,1) ;

    %---------------------------------------------------------------------------
    % merge adjacent singletons into a single upper triangular block
    %---------------------------------------------------------------------------

    bsize = r (2:nblocks+1) - r (1:nblocks) ;
    t = [0 (bsize == 1)] ;
    z = (t (1:nblocks) == 0 & t (2:nblocks+1) == 1) | t (2:nblocks+1) == 0 ;
    y = [(find (z)) nblocks+1] ;
    r = r (y) ;
    nblocks = length (y) - 1 ;
    is_triangular = y (2:nblocks+1) - y (1:nblocks) > 1 ;
    % bsize = r (2:nblocks+1) - r (1:nblocks) ;
    clear t z y bsize
    Info (18) = nblocks ;

    %---------------------------------------------------------------------------
    % solve the system: x (q) = B\c
    %---------------------------------------------------------------------------

    for k = nblocks:-1:1

        k1 = r (k) ;
        k2 = r (k+1) - 1 ;
        nk = k2-k1+1 ;
        Bk = B (k1:k2, k1:k2) ;

        if (is_triangular (k))

            % back substitution only
            t = cputime ;
            x (k1:k2) = Bk \ c (k1:k2) ;
            t = cputime - t ;

            % solve time
            nnzB = nnz (Bk) ;
            Info (85) = Info (85) + 2*nnzB ;
            Info (86) = Info (86) + t ;

            % nz in LU:
            Info ([21 41   ]) = Info ([21 41   ]) + nnzB ;
            Info ([24 44   ]) = Info ([24 44   ]) + nk ;
            Info ([25 45 48]) = Info ([25 45 48]) + nnzB ;

        elseif (nk < 4)

            % solve it as a dense linear system
            x (k1:k2) = full (Bk) \ c (k1:k2) ;

            nzl = ((nk^2 - nk) / 2) + nk ;
            Info ([21 41 48]) = Info ([21 41 48]) + nk^2 ;
            Info ([24 44   ]) = Info ([24 44   ]) + nzl ;
            Info ([25 45   ]) = Info ([25 45   ]) + nzl ;

        else

            % solve it as a sparse linear system
            [xk, BInfo] = umfpack (Bk, '\', c (k1:k2), Control) ;
            x (k1:k2) = xk ;
            clear xk

            Info (terms_peak) = max (Info (terms_peak), ...
                 Info (terms_prev) + BInfo (terms_peak)) ;
            Info (terms_latest) = BInfo (terms_latest) ;
            Info (terms_sum) = Info (terms_sum) + BInfo (terms_sum) ;
            Info (terms_max) = max (Info (terms_max), BInfo (terms_max)) ;

        end

        % off-diagonal block back substitution
        Boff = B (1:k1-1, k1:k2) ;
        t = cputime ;
        c (1:k1-1) = c (1:k1-1) - Boff * x (k1:k2) ;
        t = cputime - t ;

        nnzoff = nnz (Boff) ;
        Info (19) = Info (19) + nnzoff ;

        % solve time
        Info (85) = Info (85) + 2*nnzoff ;
        Info (86) = Info (86) + t ;

        clear Bk Boff xk

    end

    x (q) = x ;

end

tfactor = cputime - tfactor ;
tall = tdmperm + tfactor ;

