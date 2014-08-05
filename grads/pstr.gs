*** pstr.gs ********************************************************************
*
* plot a string with all setting
*
* Usage: pstr x y string <color# <size <thick <ref_pt <shift> <vshift>>>>>
*
********************************************************************************

function pstr (args)

msg = 'Usage: pstr x y string <color# <size <thick <ref_pt <shift> <vshift>>>>>'

x      = subwrd(args, 1)
y      = subwrd(args, 2)
str    = subwrd(args, 3)
color  = subwrd(args, 4)
size   = subwrd(args, 5)
thick  = subwrd(args, 6)
ref    = subwrd(args, 7)
arg8   = subwrd(args, 8)
arg9   = subwrd(args, 9)
if (str = '') ; say msg ; return ; endif
if (color = '') ; color = 1   ; endif
if (size  = '') ; size  = 0.2 ; endif
if (thick = '') ; thick = 3   ; endif
if (ref   = '') ; ref   = 'c' ; endif
if (arg8  = '') ; arg8  = 0   ; endif
if (arg9  = '') ; arg9  = 0   ; endif

xs = 0
ys = 0
if (ref = 'l')  ; xs =      arg8 ;                  endif
if (ref = 'r')  ; xs = -1 * arg8 ;                  endif
if (ref = 'bc') ;                  ys =      arg9 ; endif
if (ref = 'tc') ;                  ys = -1 * arg9 ; endif
if (ref = 'bl') ; xs =      arg8 ; ys =      arg9 ; endif
if (ref = 'br') ; xs = -1 * arg8 ; ys =      arg9 ; endif
if (ref = 'tl') ; xs =      arg8 ; ys = -1 * arg9 ; endif
if (ref = 'tr') ; xs = -1 * arg8 ; ys = -1 * arg9 ; endif

x = x + xs
y = y + ys
hsize = size * 0.75

'set string 'color' 'ref' 'thick
'set strsiz 'hsize' 'size
'draw string 'x' 'y' 'str

return

********************************************************************************
