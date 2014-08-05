#/bin/bash
#
filetype=sym
root_name=Ctrl1
input=${root_name}_uv_plev
infile=../${input}
infilecn=../${root_name}_uv10_pv.gtk
outfile=${root_name}_${filetype}
dhr=1
rm -f ${outfile}.dat
rm -f ${filetype}_tmp1.dat
rm -f ${outfile}.ctl
touch ${outfile}.dat
touch ${outfile}.ctl
ti=3
tf=45
#
for((i=${ti};i<=${tf};i++))
do
t1=$[(${i}-1)*${dhr}+1]
t2=$[${t1}+${dhr}-1]
echo 't1, t2='${t1}' '${t2}
gradsc -bl <<eor
  azi_wh.gs ${infile} ${infilecn} ${t1} ${t2}
eor
cat ${filetype}_tmp1.dat >> ${outfile}.dat
done
#
gradsc -bl <<eor
  write_p_ctl.gs ${input} ${outfile} ${ti} ${tf} ${filetype}
eor
