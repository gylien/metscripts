function write_ctl(args)

************* README (write_ctl.gs) ************
*   This file is designed to output ctl-file au-
* tomatically.
*   This is not a flexible as file.
*   For *_sym.dat and *_asym.dat

**************** HOW TO RUN ***************
* This is called by doa*.bash automatically

******************************(Shiuan,2008,Feb.)
 
 infile=subwrd(args,1)
 ti=subwrd(args,2)
 tf=subwrd(args,3)
 filetype=subwrd(args,4)
 'reinit'
 'open ../../'infile''
 'set t 'ti''
 t=tf-ti+1
 'q time'
 say result
 time=subwrd(result,3)
 ctlfile=''infile'_'filetype''

var_num=2
zlev=31
******** Output the vtl file ******** 
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
 outstr=' 0.99850'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.99600'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.99350'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.99100'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.98750'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.98250'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.97750'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.97250'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.96750'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.96250'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.95500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.94000'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.91000'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.87000'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.82500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.77500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.72500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.67500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.62500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.57500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.52500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.47500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.42500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.37500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.32500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.27500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.22500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.17500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.12500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.07500'
 rc = write(''ctlfile'.ctl',outstr)
 outstr=' 0.02500'
 rc = write(''ctlfile'.ctl',outstr)
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
