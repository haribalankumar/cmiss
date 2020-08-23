$1 != "name.txt" && $1 != "fort.9" && $1 != "generated" && $1 !~ /cmiss_data/ && $1 !~ /.html/ && $1 !~ /.out/ && $1 !~ /.changes/  && $1 !~ /.dif/ { temp = "dirname";
	gsub("/product/cmiss/data","",temp);
	printf "  <LI> <A HREF=\"%s/%s\">%s</A>\n",temp,$1,$1}
