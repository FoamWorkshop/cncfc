#!/bin/bash
cncfc_dir="./../../../../cncfc"
dxf_dir="./../../cad"

$cncfc_dir/dxf2knots.py -i $dxf_dir/ex_01.dxf -l 0000d -op 1
$cncfc_dir/dxf2knots.py -i $dxf_dir/ex_01.dxf -l 0100d -op 1
$cncfc_dir/dxf2knots.py -i $dxf_dir/ex_01.dxf -l 0200d -op 1
$cncfc_dir/dxf2knots.py -i $dxf_dir/ex_01.dxf -l 0300d -op 1

$cncfc_dir/knots2gcode.py -i 0000d1.knt -o ex_01_0000 -sh
$cncfc_dir/knots2gcode.py -i 0100d1.knt -o ex_01_0100 -sh
$cncfc_dir/knots2gcode.py -i 0200d1.knt -o ex_01_0200 -sh
$cncfc_dir/knots2gcode.py -i 0300d1.knt -o ex_01_0300 -sh

cut_part='ex_01_cut.ngc'
touch $cut_part
cat > $cut_part << EOL

G21 (Units in millimeters)  G90 (Absolute programming)
F180
o<ex_01_0000> call
o<ex_01_0100> call
o<ex_01_0200> call
o<ex_01_0300> call
M2

EOL
