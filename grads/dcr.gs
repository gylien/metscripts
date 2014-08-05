****************************************

function dcr (args)

clev1 = subwrd(args, 1)
clev2 = subwrd(args, 2)
cint  = subwrd(args, 3)
dargs = subwrd(args, 4)

clevsr(clev1, clev2, cint)

'd 'dargs

return

****************************************

function clevsr (clev1, clev2, cint)

levstr=''
cl = clev1
while (cl <= clev2)
  levstr = levstr' 'cl
  cl = cl * cint
endwhile

say 'set clevs'levstr
'set clevs'levstr

return

****************************************
