*** pmark.gs *******************************************************************
*
* plot a mark with all setting
*
* Usage: pmark x y <marktype <color# <size <thickness>>>>
*
********************************************************************************

function pmark (args)

msg = 'Usage: pmark x y <marktype <color# <size <thickness>>>>'

x     = subwrd(args, 1)
y     = subwrd(args, 2)
type  = subwrd(args, 3)
color = subwrd(args, 4)
size  = subwrd(args, 5)
thick = subwrd(args, 6)
if (y = '') ; say msg ; return ; endif
if (type  = '') ; type  = 1   ; endif
if (color = '') ; color = 1   ; endif
if (size  = '') ; size  = 0.1 ; endif
if (thick = '') ; thick = 3   ; endif

'set line 'color' 1 'thick
'draw mark 'type' 'x' 'y' 'size

return

********************************************************************************
