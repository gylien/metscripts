function write_ctl(args)

************* README (write_ctl.gs) ************ 
*   This file is designed to output ctl-file au-
* tomatically.
*   This is not a flexible as file.
*   For *_sym.dat and *_asym.dat

**************** HOW TO RUN ***************
* This is called by doa*.bash automatically

************************************************

 infile=subwrd(args,1)
 ctlfile=subwrd(args,2)
 ti=subwrd(args,3)
 tf=subwrd(args,4)
 filetype=subwrd(args,5)
 'reinit'
 'open ../'infile'.ctl'
 'set t 'ti''
 t=tf-ti+1
 'q time'
 say result
 time=subwrd(result,3)

var_num=2
'q file'
data=sublin(result,5)
zlev=subwrd(data,9)
******** Output the ctl file ******** 
 outstr='dset ^'ctlfile'.dat'
 rc = write(''ctlfile'.ctl',outstr)
 outstr='byteswapped'
 rc = write(''ctlfile'.ctl',outstr)
 outstr='title MM5 data'
 rc = write(''ctlfile'.ctl',outstr)
 outstr='undef -9999.'
 rc = write(''ctlfile'.ctl',outstr)
 outstr='xdef  41 linear  0. 0.0225'
 rc = write(''ctlfile'.ctl',outstr)
 outstr='ydef  1 linear   0  1.'
 rc = write(''ctlfile'.ctl',outstr)
 outstr='zdef  'zlev' levels'
 rc = write(''ctlfile'.ctl',outstr)
**** vertical level infomation
zz=1
while(zz<=zlev)
 'set z 'zz''
 plev=subwrd(result,4)
 outstr=' 'plev'.00'
 rc = write(''ctlfile'.ctl',outstr)
 zz=zz+1
endwhile
******************************
 outstr='tdef  't' linear 00:'time'   60MN'
 rc = write(''ctlfile'.ctl',outstr)
 outstr='vars 'var_num''
 rc = write(''ctlfile'.ctl',outstr)
 if(filetype='sym')
 outstr='ur   'zlev' 99 the the radial wind speed'
 rc = write(''ctlfile'.ctl',outstr)
 outstr='ut   'zlev' 99 the tangential wind speed'
 rc = write(''ctlfile'.ctl',outstr)
 endif
 if(filetype='asym')
 outstr='uasy   'zlev' 99 the asymmetric part of the horizontal wind'
 rc = write(''ctlfile'.ctl',outstr)
 outstr='vasy   'zlev' 99 the asymmetric part of the moridional wind'
 rc = write(''ctlfile'.ctl',outstr)
 endif
 outstr='endvars'
 rc = write(''ctlfile'.ctl',outstr)
