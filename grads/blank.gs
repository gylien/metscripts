*** blank.gs *******************************************************************
*
* plot a blank map
*
* Usage: blank lon1 lon2 lat1 lat2
*
********************************************************************************

function blank (args)

msg = 'Usage: blank lon1 lon2 lat1 lat2'

lon1 = subwrd(args, 1)
lon2 = subwrd(args, 2)
lat1 = subwrd(args, 3)
lat2 = subwrd(args, 4)
if (lat2 = '') ; say msg ; return ; endif

'reinit.gs'
'open dummy.ctl'

'set lon 'lon1' 'lon2
'set lat 'lat1' 'lat2

'set clevs 1'
'd dum'

return

********************************************************************************
