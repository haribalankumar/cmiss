#!/bin/tcsh
#ensures that data_page.html exists for all directorys.
if ( ! -d $1/generated ) then
        mkdir $1/generated
	touch $1/generated/data_page.html
endif
