function umfpack_make
%UMFPACK_MAKE
%
%   Compiles the UMFPACK mexFunction and then runs a simple demo.
%
%   UMFPACK Version 4.0 (Apr 11, 2002).  Copyright (c) 2002 by Timothy A.
%   Davis.  All Rights Reserved.  Type "help umfpack_details" for License.
%
%   See also umfpack, umfpack_details, umfpack_report,
%   umfpack_demo, and umfpack_simple.

fprintf ('\nNow compiling the umfpack mexFunction.\n') ;

obj = 'o' ;
blas_lib = '' ;
if (ispc)
    obj = 'obj' ;
end

% BLAS option
fprintf ('\nUsing the BLAS is faster, but might not compile correctly.\n') ;
fprintf ('\nPlease select one of the following options: \n') ;
fprintf ('   1:  attempt to compile with the BLAS (default)\n') ;
fprintf ('   2:  do not use the BLAS\n') ;
blas = input ('or type control-C if you do not wish to proceed: ') ;
if (isempty (blas))
    blas = 1 ;
end
if (blas == 1)
    % try to link to MATLAB's built-in BLAS
    blas = '' ;
    if (ispc)
        % the default lcc compiler needs this library to access the BLAS
        blas_lib = 'libmwlapack.lib' ;
	msg = [ ...
	'\nTo use the BLAS, you must first copy the libmwlapack.lib file ', ...
	'from the\numfpack\\lcc_lib\\ directory to the ', ...
	'<matlab>\\extern\\lib\\win32\\lcc\\\ndirectory, where <matlab> ', ...
	'is where MATLAB is installed.  Next, type\n\n    mex -setup\n\n', ...
	'at the MATLAB prompt, and ask MATLAB to select the lcc compiler.  ',...
	'You can skip\nall of this if you have already done it, or have ', ...
	'configured mex to use\na different compiler.  If you are using ', ...
	'Norton anti-virus software on Windows\n98SE, then you need to ', ...
	'exit MATLAB, turn off virus checking, and restart MATLAB\n', ...
	'before you can use the mex command or compile UMFPACK.\n', ...
	'\nHit enter to continue, or type control-C if you do not wish to '] ;
	fprintf (msg) ;
	input ('proceed: ') ;
    end
else
    % No BLAS
    blas = '-DNBLAS' ;
end

% mex command
mx = sprintf ('mex -inline -O -DNDEBUG %s ', blas) ;
fprintf ('\nCompile options:\n%s\nNow compiling.  Please wait ...\n\n', mx) ;

umf = { 'analyze', 'apply_order', 'assemble', 'blas3_update', ...
    'build_tuples', 'build_tuples_usage', 'colamd', 'create_element', ...
    'dump', 'extend_front', 'free', 'garbage_collection', 'get_memory', ...
    'init_front', 'is_permutation', 'kernel', 'kernel_init', ...
    'kernel_init_usage', 'kernel_wrapup', 'local_search', 'lsolve', ...
    'ltsolve', 'malloc', 'mem_alloc_element', 'mem_alloc_head_block', ...
    'mem_alloc_tail_block', 'mem_free_tail_block', ...
    'mem_init_memoryspace', 'order_front_tree', 'report_perm', ...
    'realloc', 'report_vector', 'row_search', 'scale_column', ...
    'set_stats', 'solve', 'symbolic_usage', 'transpose', ...
    'tuple_lengths', 'usolve', 'utsolve', 'valid_numeric', ...
    'valid_symbolic' } ;

umfu = { 'col_to_triplet', 'defaults', 'free_numeric', 'free_symbolic', ...
    'get_numeric', 'get_lunz', 'get_symbolic', 'numeric', 'qsymbolic', ...
    'report_control', 'report_info', 'report_matrix', 'report_numeric', ...
    'report_perm', 'report_status', 'report_symbolic', 'report_triplet', ...
    'report_vector', 'solve', 'symbolic', 'transpose', 'triplet_to_col' } ;

M = sprintf ('%s -v -output umfpack umfpackmex.c ', mx) ;

%-------------------------------------------------------------------------------
% compile umf_* routines
%-------------------------------------------------------------------------------

for k = 1:length (umf)

    cmd (sprintf ('%s -DDINT -c umf_%s.c', mx, umf {k})) ;
    src = sprintf ('umf_%s.%s', umf {k}, obj) ;
    dst = sprintf ('umf_md_%s.%s', umf {k}, obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DZINT -c umf_%s.c', mx, umf {k})) ;
    src = sprintf ('umf_%s.%s', umf {k}, obj) ;
    dst = sprintf ('umf_mz_%s.%s', umf {k}, obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

end

%-------------------------------------------------------------------------------
% compile umfpack_* routines
%-------------------------------------------------------------------------------

for k = 1:length (umfu)

    cmd (sprintf ('%s -DDINT -c umfpack_%s.c', mx, umfu {k})) ;
    src = sprintf ('umfpack_%s.%s', umfu {k}, obj) ;
    dst = sprintf ('umfpack_md_%s.%s', umfu {k}, obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DZINT -c umfpack_%s.c', mx, umfu {k})) ;
    src = sprintf ('umfpack_%s.%s', umfu {k}, obj) ;
    dst = sprintf ('umfpack_mz_%s.%s', umfu {k}, obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

end


%-------------------------------------------------------------------------------
% umf_*hsolve:
%-------------------------------------------------------------------------------

    cmd (sprintf ('%s -DDINT -DCONJUGATE_SOLVE -c umf_ltsolve.c', mx)) ;
    src = sprintf ('umf_ltsolve.%s', obj) ;
    dst = sprintf ('umf_md_lhsolve.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DZINT -DCONJUGATE_SOLVE -c umf_ltsolve.c', mx)) ;
    src = sprintf ('umf_ltsolve.%s', obj) ;
    dst = sprintf ('umf_mz_lhsolve.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DDINT -DCONJUGATE_SOLVE -c umf_utsolve.c', mx)) ;
    src = sprintf ('umf_utsolve.%s', obj) ;
    dst = sprintf ('umf_md_uhsolve.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DZINT -DCONJUGATE_SOLVE -c umf_utsolve.c', mx)) ;
    src = sprintf ('umf_utsolve.%s', obj) ;
    dst = sprintf ('umf_mz_uhsolve.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

%-------------------------------------------------------------------------------
% umf_*_triplet_*map_*x, all from the single umf_triplet.c file:
%-------------------------------------------------------------------------------

    cmd (sprintf ('%s -DDINT -DDO_MAP -DDO_VALUES -c umf_triplet.c', mx)) ;
    src = sprintf ('umf_triplet.%s', obj) ;
    dst = sprintf ('umf_md_triplet_map_x.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DZINT -DDO_MAP -DDO_VALUES -c umf_triplet.c', mx)) ;
    src = sprintf ('umf_triplet.%s', obj) ;
    dst = sprintf ('umf_mz_triplet_map_x.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DDINT -DDO_MAP -c umf_triplet.c', mx)) ;
    src = sprintf ('umf_triplet.%s', obj) ;
    dst = sprintf ('umf_md_triplet_map_nox.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DZINT -DDO_MAP -c umf_triplet.c', mx)) ;
    src = sprintf ('umf_triplet.%s', obj) ;
    dst = sprintf ('umf_mz_triplet_map_nox.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DDINT -DDO_VALUES -c umf_triplet.c', mx)) ;
    src = sprintf ('umf_triplet.%s', obj) ;
    dst = sprintf ('umf_md_triplet_nomap_x.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DZINT -DDO_VALUES -c umf_triplet.c', mx)) ;
    src = sprintf ('umf_triplet.%s', obj) ;
    dst = sprintf ('umf_mz_triplet_nomap_x.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DDINT -c umf_triplet.c', mx)) ;
    src = sprintf ('umf_triplet.%s', obj) ;
    dst = sprintf ('umf_md_triplet_nomap_nox.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DZINT -c umf_triplet.c', mx)) ;
    src = sprintf ('umf_triplet.%s', obj) ;
    dst = sprintf ('umf_mz_triplet_nomap_nox.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

%-------------------------------------------------------------------------------
% umfpack_wsolve:
%-------------------------------------------------------------------------------

    cmd (sprintf ('%s -DDINT -DWSOLVE -c umfpack_solve.c', mx)) ;
    src = sprintf ('umfpack_solve.%s', obj) ;
    dst = sprintf ('umfpack_md_wsolve.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

    cmd (sprintf ('%s -DZINT -DWSOLVE -c umfpack_solve.c', mx)) ;
    src = sprintf ('umfpack_solve.%s', obj) ;
    dst = sprintf ('umfpack_mz_wsolve.%s', obj) ;
    mvfile (src, dst) ;
    M = [M, ' ', dst] ;

%-------------------------------------------------------------------------------
% umfpack_timer:
%-------------------------------------------------------------------------------

    cpfile ('umfpack_timer.c', 'umfpack_mtimer.c') ;
    cmd (sprintf ('%s -c umfpack_mtimer.c', mx)) ;
    rmfile ('umfpack_mtimer.c') ;
    dst = sprintf ('umfpack_mtimer.%s', obj) ;
    M = [M, ' ', dst] ;

%-------------------------------------------------------------------------------
% add the BLAS library
%-------------------------------------------------------------------------------

    M = [M, ' ', blas_lib] ;
    
%-------------------------------------------------------------------------------
% Finally, compile the mexFunction:
%-------------------------------------------------------------------------------

    cmd (M) ;

fprintf ('\nCompilation has completed.  Now trying the umfpack_simple demo.\n');
umfpack_simple



%-------------------------------------------------------------------------------
% rmfile:  delete a file, but only if it exists
%-------------------------------------------------------------------------------

function rmfile (file)
if (length (dir (file)) > 0)
    fprintf ('delete %s\n', file) ;
    delete (file) ;
end

%-------------------------------------------------------------------------------
% cpfile:  copy the src file to the filename dst, overwriting dst if it exists
%-------------------------------------------------------------------------------

function cpfile (src, dst)
rmfile (dst)
fprintf ('copy %s to %s\n', src, dst) ;
if (length (dir (src)) == 0)
    error (sprintf ('File does not exist: %s\n', src)) ;
end
copyfile (src, dst) ;

%-------------------------------------------------------------------------------
% mvfile:  move the src file to the filename dst, overwriting dst if it exists
%-------------------------------------------------------------------------------

function mvfile (src, dst)
cpfile (src, dst) ;
rmfile (src) ;

%-------------------------------------------------------------------------------
% cmd:  display and execute a command
%-------------------------------------------------------------------------------

function cmd (s)
fprintf ('%s\n', s) ;
eval (s) ;
