#!/bin/bash
cncfc_dir=/home/adam/Documents/00.projects/02.python/cncfc

layer_prefix="fus_01_"
counter=0
for FILE in $layer_prefix"_xyuv_"*.knt
do
  name="${FILE%%.knt}"
  prof_numb="${name##$layer_prefix'_xyuv_'}"
  prof_r=$layer_prefix"_b_"$prof_numb.knt
  echo found matching:
    if [ -f "$prof_r" ]; then
      o_file="$layer_prefix""$prof_numb"

      echo "data set:" $counter
      echo found matching rotation data: "$prof_r"
      echo found matching translation data: "$FILE"
      echo output file: "$o_file"

      $cncfc_dir/knots2gcode.py -i $FILE -ir $prof_r -o $o_file -d 418 -cm -sh
      let counter=counter+1
    else
      echo did NOT found matching: "$prof_r" - EXIT
      break
    fi
done

cut_file="cut_""${layer_prefix%%'_'}"".ngc"
touch $cut_file

cat > $cut_file << EOL
(lofted profiles cut test 1)
G21 G90

F200
G0 B0

G0 X60 Y-5 U60 V-5
EOL

for FILE in fus_01*.ngc
do
echo "o<""${FILE%%.ngc}""> call" >> $cut_file
echo "G0 X50 U50" >> $cut_file
done

cat >> $cut_file << EOL
G0 X60 U60
M2
EOL
