#General CMISS octave functions

#Global variable definitions

global cmiss_diag = 0;

#####

#Global function definitions

#####

function cmiss_flag_warning( warning_message )

  printf(">>WARNING: %s\n",warning_message);

endfunction

#####

function cmiss_flag_error( error_message )

  error(">>ERROR: %s\n",error_message);

endfunction

#####

function cmiss_diag_message( diag_message )

  printf(">>DIAG: %s\n",diag_message);

endfunction

#####

function ret_string = cmiss_chardata2string( chardata, number )

  ret_string = "";
  if( number > 0 )
    for i  = 1:number
      ret_string = sprintf("%s%c",ret_string, chardata(i));
    endfor
  endif

endfunction

#####

function [ fileid, numbertags, filebinvertype, errcode ] = cmiss_open_bin_file( type, version, command, extent, file )

  global cmiss_diag;

  errcode = 0;

  extension = strcat( ".bin", extent );
  filename = strcat( file, extension );

  if( cmiss_diag ) 
    diag_message = sprintf("The command is %s",command);
    cmiss_diag_message( diag_message );
    diag_message = sprintf("The filename is %s",filename);
    cmiss_diag_message( diag_message );
  endif
   
  if( strcmp( command, "read") )

    #Open the binary file for reading

    fileid = fopen( filename, "r+");

    #Read identity head section
     
    [ chardata, count ] = fread( fileid, 2, "char" );
    if( chardata(1) != 0x7 ) 
      cmiss_flag_error("Not a CMISS binary file");
    endif
    filebinvertype = chardata(2);
    if( cmiss_diag )
      cmiss_diag_message("File identity header information:");
      diag_message = sprintf("   Binary file header format version is %d",filebinvertype);
      cmiss_diag_message( diag_message );
    endif
    if( filebinvertype == 0x0 )
      #Identiry format 0 - has 2 extra dummy bytes. Skip them
      cmiss_flag_warning("Old binary file identity format found. Please update file");
      [ chardata, count ] = fread( fileid, 2, "char" );      
    elseif( filebinvertype == 0x1 )
      #Identity format 1 - has machine header section
      cmiss_flag_warning("Old binary file identity format found. Please update file");
      [ chardata, count ] = fread( fileid, 1, "char" );
      nummachineheaderbytes = chardata(1);
      [ chardata, count ] = fread( fileid, nummachineheaderbytes, "char" );
      if( nummachineheaderbytes == 11)
        filemachtype = chardata(1);
        fileostype = chardata(2);
        fileendiantype = chardata(3);
        filespformtype = chardata(4);
        filedpformtype = chardata(5);
        filecharsize = chardata(6);
        fileintsize = chardata(7);
        filesintsize = chardata(8);
        filespsize = chardata(9);
        filedpsize = chardata(10);
        filelogsize = chardata(11);
        if( cmiss_diag ) 
          cmiss_diag_message("File machine header information:");
          diag_message = sprintf("   Machine type: %d",filemachtype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Operating system type: %d",fileostype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Endian type: %d",fileendiantype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Single precision format type: %d",filespformtype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Double precision format type: %d",filedpformtype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Character size: %d",filecharsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Integer size: %d",fileintsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Short integer size: %d",filesintsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Single precision size: %d",filespsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Double precision size: %d",filedpsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Logical size: %d",filelogsize);
          cmiss_diag_message( diag_message );
        endif
      else
        cmiss_flag_error("Invalid number of machine header bytes");
      endif
    elseif( filebinvertype == 0x2 )
      #Identity format 2 - has expanded machine header section and subtags
      [ chardata, count ] = fread( fileid, 1, "char" );
      nummachineheaderbytes = chardata(1);
      [ chardata, count ] = fread( fileid, nummachineheaderbytes, "char" );
      if( nummachineheaderbytes == 16)
        filemachtype = chardata(1);
        fileostype = chardata(2);
        fileendiantype = chardata(3);
        filecharformtype = chardata(4);
        fileintformtype = chardata(5);
        filespformtype = chardata(6);
        filedpformtype = chardata(7);
        filecharsize = chardata(8);
        fileintsize = chardata(9);
        filesintsize = chardata(10);
        filelintsize = chardata(11);
        filespsize = chardata(12);
        filedpsize = chardata(13);
        filelogsize = chardata(14);
        filespcsize = chardata(15);
        filedpcsize = chardata(16);
        if( cmiss_diag ) 
          cmiss_diag_message("File machine header information:");
          diag_message = sprintf("   Machine type: %d",filemachtype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Operating system type: %d",fileostype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Endian type: %d",fileendiantype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Character format type: %d",filecharformtype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Integer format type: %d",fileintformtype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Single precision format type: %d",filespformtype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Double precision format type: %d",filedpformtype);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Character size: %d",filecharsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Integer size: %d",fileintsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Short integer size: %d",filesintsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Long integer size: %d",filelintsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Single precision size: %d",filespsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Double precision size: %d",filedpsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Logical size: %d",filelogsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Single precision complex size: %d",filespcsize);
          cmiss_diag_message( diag_message );
          diag_message = sprintf("   Double precision complex size: %d",filedpcsize);
          cmiss_diag_message( diag_message );
        endif
      else
        cmiss_flag_error("Invalid number of machine header bytes");
      endif
    else 
      cmiss_flag_error("Unknown binary file header identity format");
    endif
 
    #Read file header section
    #Read file type, version and header

    [ intdata, count ] = fread( fileid, 1, "int32" );
    filetype = intdata(1);
    if( filetype != type )
      warning_message = sprintf("File has a different file type\n   File type is %d\n   Expected file type is %d",filetype, type);
      cmiss_flag_warning(warning_message);
    endif

    if( filebinvertype == 0x0 || filebinvertype == 0x1 )
      [ spdata, count ] = fread( fileid, 1, "float" );
      binfileversion(1) = spdata(1);
      if( binfileversion(1) != version(1) ) 
        warning_message = sprintf("File has a different version number\n   File version is %f\n   Expected version is %f",binfileversion(1),version(1));
        cmiss_flag_warning(warning_message);
      endif
    else
      [ intdata, count ] = fread( fileid, 3, "int32" );
      binfileversion(1) = intdata(1);
      binfileversion(2) = intdata(2);
      binfileversion(3) = intdata(3);
      if( binfileversion(1) != version(1) || binfileversion(2) != version(2) || binfileversion(3) || version(3) )
        warning_message = sprintf("File has a different version number\n   File version is %d.%d.%d\n   Expected version is %d.%d.%d",binfileversion(1),binfileversion(2),binfileversion(3),version(1),version(2),version(3));
        cmiss_flag_warning(warning_message);
      endif
    endif

    #Read heading
    [ intdata, count ] = fread( fileid, 1, "int32" );
    headsize = intdata(1);

    [ chardata, count ] = fread( fileid, headsize, "char" );
    if( headsize == 0 )
      fprintf("File heading:\n");
    else
      fprintf("File heading: %s\n",cmiss_chardata2string(chardata,headsize));
    endif

    [ intdata, count ] = fread( fileid, 1, "int32" );
    numbertags = intdata(1);

    if( cmiss_diag )
      cmiss_diag_message("File header information:");
      diag_message = sprintf("  File type: %d",filetype);
      cmiss_diag_message( diag_message );
      diag_message = sprintf("  File version: %d.%d.%d",binfileversion(1),binfileversion(2),binfileversion(3) );
      cmiss_diag_message( diag_message );
      diag_message = sprintf("  Number of header bytes: %d", headsize);
      cmiss_diag_message( diag_message );
      diag_message = sprintf("  Number of tags: %d", numbertags);
      cmiss_diag_message( diag_message );
    endif

  else
    cmiss_flag_error("Opening a cmiss binary file for writing is not implemented");
  endif

endfunction

#####

function [ tagindex, numsubtags, numtagbytes, errcode ] = cmiss_read_bin_tag_header( fileid, filebinvertype )

  global cmiss_diag;

  errcode = 0;

  [ intdata, count ] = fread( fileid, 2, "int32" );
  tagindex = intdata(1);
  numtagheaderbytes = intdata(2);
 
  [ chardata, count ] = fread( fileid, numtagheaderbytes, "char" );

  if( filebinvertype == 0x2 )
    [ intdata, count ] = fread( fileid, 1, "int32" );
    numsubtags = intdata(1);
  else
    numsubtags = 0;
  endif

  if( numsubtags == 0 )
    [ intdata, count ] = fread( fileid, 1, "int32" );
    numtagbytes = intdata(1);
  endif

  if( cmiss_diag )
    cmiss_diag_message( "Tag header information:" );
    diag_message = sprintf("   Tag index = %d", tagindex);
    cmiss_diag_message( diag_message );
    diag_message = sprintf("   Number tag headerbytes = %d", numtagheaderbytes);
    cmiss_diag_message( diag_message );
    diag_message = sprintf("   Tag header = %s", cmiss_chardata2string( chardata, numtagheaderbytes ) );
    cmiss_diag_message( diag_message );
    diag_message = sprintf("   Number subtags = %d", numsubtags );
    cmiss_diag_message( diag_message );
    diag_message = sprintf("   Number tag bytes = %i", numtagbytes );
    cmiss_diag_message( diag_message );
  endif

endfunction

#####

function return_matrix = cmiss_read_bin_matrix( filename )

  #Open binary file

  filetype = 1;
  version = [ 1, 1, 0 ];
  [ fileid, numbertags, filebinvertype, errcode ] = cmiss_open_bin_file( filetype, version, "read", "mat", filename );

  #Read matrices (tags)
  if( numbertags == 1 )
    [ tagindex, numsubtags, numtagbytes, errcode ] = cmiss_read_bin_tag_header( fileid, filebinvertype );
    #Read the number of indicies
    [ intdata, count ] = fread( fileid, 1, "int32" );
    numindicies = intdata(1);
    if( numindicies == 2 )
      [ intdata, count ] = fread( fileid, numindicies, "int32" );
      matrix_m = intdata(1);
      matrix_n = intdata(2);
      return_matrix = zeros( matrix_m, matrix_n );
      #Read in the size of the matrix
      [ intdata, count ] = fread( fileid, 1, "int32" );
      nztot = intdata(1);
      #Read in the number of values stored
      [ intdata, count ] = fread( fileid, 1, "int32" );
      numvalues = intdata(1);
      #Read in the timecode
      [ intdata, count ] = fread( fileid, 1, "int32" );
      timecode = intdata(1);
      if( timecode == 1 ) then
        [ dpdata, count ] = fread( fileid, 1, "double" );
        timecode = dpdata(1);
      endif
      #Read in the sparsity type
      [ intdata, count ] = fread( fileid, 1, "int32" );
      sparsitycode = intdata(1);
      if( sparsitycode > 0 )
        sparsity_type = sparsitycode;
      else
        sparsity_type = 0;
      endif
      if( sparsity_type > 0 )
        if( sparsity_type == 1 ) #Compressed row
          #Read in the row offsets
          [ isr_matrix, count ] = fread( fileid, matrix_m+1, "int32" );
          #Read in the column numbers
          [ isc_matrix, count ] = fread( fileid, isr_matrix(matrix_m+1)-1, "int32" );
          nztot = isr_matrix(matrix_m+1)-1;
        elseif( sparsity_type == 2 ) #Row-column
          #Read the number of non-zeros
          [ intdata, count ] = fread( fileid, 1, "int32" );
          nztot = intdata(1);
          #Read the row numbers
          [ isr_matrix, count ] = fread( fileid, nztot, "int32" );
          #Read the column numbers
          [ isc_matrix, count ] = fread( fileid, nztot, "int32" );
        elseif( sparsity_type == 3 ) #Compressed column
	  #Read the column offsets
          [ isc_matrix, count ] = fread( fileid, matrix_n+1, "int32" );
          #Read in the row numbers
          [ isr_matrix, count ] = fread( fileid, isc_matrix(matrix_n+1)-1, "int32" );
          nztot = isc_matrix(matrix_n+1)-1;          
        elseif( sparsity_type == 5 ) #UMFpack row-column
          #Read in the number of non-zeros
          [ intdata, count ] = fread( fileid, 1, "int32" );
          nztot = intdata(1);
          #Read in the row and column numbers
          [ isc_matrix, count ] = fread( fileid, 2*nztot, "int32" );
        else
          cmiss_flag_error("Unknown sparisty code");
        endif
      endif
      #Read in the values
      if( sparsity_type == -1 ) 
        for column = 1:matrix_n
          for row = 1:matrix_m
     	     [ dpdata, count ] = fread( fileid, 1, "double" );
             return_matrix(row, column) = dpdata(1);
          endfor
        endfor
      elseif( sparsity_type == 0 ) 
        for column = 1:matrix_n
          for row = 1:matrix_m
     	     [ dpdata, count ] = fread( fileid, 1, "double" );
             return_matrix(row, column) = dpdata(1);
          endfor
        endfor
      elseif( sparsity_type == 1 ) 
	cmiss_flag_warning("Not done yet");
      elseif( sparsity_type == 2 ) 
	cmiss_flag_warning("Not done yet");
      elseif( sparsity_type == 3 ) 
	cmiss_flag_warning("Not done yet");
      elseif( sparsity_type == 4 ) 
	cmiss_flag_warning("Not done yet");
      elseif( sparsity_type == 5 ) 
	cmiss_flag_warning("Not done yet");
      else
        cmiss_flag_error("Unknown sparsity type");
      endif
    else
      cmiss_flag_error("Number of indicies !=2 is not implemented");
    endif
  else
    cmiss_flag_error("Matrix files containing more than one matrix is not implemented");
  endif

  fclose( fileid ) ;

endfunction
