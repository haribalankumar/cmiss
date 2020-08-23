for file in `find . -name "*.uil"`
do
  i=${file%.*}
  uil -o $i.uid $i.uil
done

for new_file in `find . -name "*.uid"`
do
  j=${new_file%.*}
  perl /hpc/cmiss/cmgui/source/utilities/uid2uidh.pl $j.uid $j.uidh
done

