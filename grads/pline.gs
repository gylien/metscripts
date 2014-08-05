*** pline.gs *******************************************************************
*
* plot a line with all setting
*
* Usage: pline x1 y1 x2 y2 <color# <style <thickness>>>
*
********************************************************************************

function pline (args)

msg = 'Usage: pline x1 y1 x2 y2 <color# <style <thickness>>>'

x1    = subwrd(args, 1)
y1    = subwrd(args, 2)
x2    = subwrd(args, 3)
y2    = subwrd(args, 4)
color = subwrd(args, 5)
style = subwrd(args, 6)
thick = subwrd(args, 7)
if (y2 = '') ; say msg ; return ; endif
if (color = '') ; color = 1 ; endif
if (style = '') ; style = 1 ; endif
if (thick = '') ; thick = 3 ; endif

'set line 'color' 'style' 'thick
'draw line 'x1' 'y1' 'x2' 'y2

return

********************************************************************************
