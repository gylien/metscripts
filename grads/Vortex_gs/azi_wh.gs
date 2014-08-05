function main(arg)
 
  'reinit'
  if(arg=''); say 'infile, inpufile of center, t1, t2'; pull arg; endif
  fname=subwrd(arg,1); fcn=subwrd(arg,2); ti=subwrd(arg,3); tf=subwrd(arg,4)
  fout='sym_tmp1.dat'
* fcn='ctrls.gtk'
  'open 'fname
* 'open 'fname'_u'
* 'open 'fname'_v'
* 'open 'fname'_w'
* 'open 'fname'_t'
* 'open 'fname'_rh'
* 'open 'fname'_th'
* 'open 'fname'_the'
* 'open 'fname'_pv'
* 'open 'fname'_vor'
  
* set parameter **************
  rr   =  0.9
  ndir  = 24
  num  =  41
* ti=1; tf=24; fout=fname%'_sym31.dat'

  'q file'; rec=sublin(result,5); ztop=subwrd(rec,9)
  zbot=1; ztop2=16
  say 'zbot, ztop, ztop2: 'zbot', 'ztop', 'ztop2

* _latc=20.
* lon.0=135.; lon.12=134.95; lon.24=134.9; lon.36=134.8; lon.48=134.6; lon.60=134.55; lon.72=134.5; lon.84=134.5; lon.96=134.5
* start do loop ******************
  'set fwrite -be '%fout

  while(ti<=tf)
    'set t 'ti
    rec=rdcenter(fcn)
*   iii=math_int(ti/12); ti1=iii*12; ti2=ti1+12; _lonc=((ti-ti1)*lon.ti2+(ti2-ti)*lon.ti1)/(ti2-ti1);
    say 't='ti', lonc='_lonc', latc='_latc
    rec=draw_azimu(ndir, rr, num, zbot, ztop, ur, 'u', 'v')
*   rec=draw_azimu(ndir, rr, num, zbot, ztop, ur, 'u.1', 'v.2')
    rec=draw_azimu(ndir, rr, num, zbot, ztop, ut, 'u', 'v')
*   rec=draw_azimu(ndir, rr, num, zbot, ztop, ut, 'u.1', 'v.2')
*   rec=draw_azimu(ndir, rr, num, zbot, ztop, tt, 'w.3')
*   rec=draw_azimu(ndir, rr, num, zbot, ztop, tt, 't.4')
*   rec=draw_azimu(ndir, rr, num, zbot, ztop, tt, 'rh.5')
*   rec=draw_azimu(ndir, rr, num, zbot, ztop, tt, 'th.6')
*   rec=draw_azimu(ndir, rr, num, zbot, ztop, tt, 'the.7')
*   rec=draw_azimu(ndir, rr, num, zbot, ztop, tt, 'pv.8')
*   rec=draw_azimu(ndir, rr, num, zbot, ztop, tt, 'vor.9')
    ti=ti+1
  endwhile

  'disable fwrite'
***** example ***************
* rec=draw_azimu(ndir, rr, num, zbot, ztop, filecenter, pt, th)
* rec=draw_azimu(ndir, rr, num, zbot, ztop, filecenter, tt, rh)
* rec=draw_azimu(ndir, rr, num, zbot, ztop, filecenter, ur, u, v)
* rec=draw_azimu(ndir, rr, num, zbot, ztop, filecenter, ut, u, v)

return

function draw_azimu(ndir, rr, num, zbot, ztop, ind, var1, var2)
*------------------------------------------------------------------
*  GrADS script:  azimu.gs
*
*  Script purpose:  (see below)
*
*  first edited:  17 Apirl 2003  by T.-S. Huang
*
*********************************************************************
* The following lines will display an arbitrary X section
* from one specified point to another.  
*
* input:
*   lonc is the position of Typhoon center at longitude point
*   latc is the position of Typhoon center at latitude point
*   rr   is the radius in degree of latitude
*   ndir is the number of direction to cut 
*           cross-section related to TY center
*        note: ndir must less than 31 for 
*              the limition of collection function
*   num  is the number of point of a cross-section
*   ztop is the number at top level of the cross-section
*   zbot is the number at bottom level of the cross-section
*   filecenter is the file name of the TY-center
*        with the form of fort.66
*   ind  is the index of the variable where:
*        ind='tt' : for any scaler variable var1
*        ind='pt' : for any scaler variable var1 to get a perturbation
*        ind='ur' : for u, v to get ur
*        ind='ut' : for u, v to get ut
* 
*
* the arbitrary cross section.
*
************************************************
******* set parameter **************************

**dr is the interval in degree of latitude
  dr   =  rr/(num-1)
  pi   = 3.1415926
  dth  =  360.*pi/180./ndir

* rec=rdcenter(filecenter)
  latc=_latc
  lonc=_lonc
* say latc' 'lonc

********
  'set z 'zbot' 'ztop
  'set x 1'
  'set y 1'

  if(ind =pt)
    rrr=rr
    'q define'
    if(subwrd(result,1)=varmean); 'undefine 'varmean; endif
    'define 'varmean'=aave('var1',lon='lonc-rrr/cos(latc*pi/180.)',lon='lonc+rrr/cos(latc*pi/180.)',lat='latc-rrr',lat='latc+rrr')'
    say result
  endif

  the  =  0.
  ii   =  1
  while( ii <= ndir)
    dlat = dr*math_cos(the)
*** at f-plane ************************************
    dlon = dr*math_sin(the)/math_cos(latc*pi/180.)
    say 'dr='dr
*   say 'ii='ii', rr='rr', dr='dr', dlat='dlat', dlon='dlon
    if( ind = tt)
      rec=crs(var1, ii, num, latc, lonc, dlat, dlon)
    endif
    if( ind = pt)
      rec= crs_ptr(var1, varmean, ii, num, latc, lonc, dlat, dlon)
    endif
    if( ind = ur)
      rec= crs_ur(var1, var2, ii, num, latc, lonc, dlat, dlon, the)
    endif
    if( ind = ut)
      rec= crs_ut(var1, var2, ii, num, latc, lonc, dlat, dlon, the)
    endif

    ii=ii+1
    the = the+dth
*   pull aa
  endwhile

* 'set xaxis 0 'rr' 'rr/5

  ii=1
  vv='(coll2gr('ii',-u)'
  while(ii<ndir)
    ii=ii+1
    vv=vv'+coll2gr('ii',-u)'
  endwhile
  vv=vv')/'ndir

  rec=draw(ndir, rr, vv, latc, lonc)
* say vv
* 'set zlog on'
* 'set lon 110 120'
* 'set ccolor rainbow'
* 'set cint 1'
* 'set gxout contour'
* 'd 'vv

return

*********************************************************************

function draw(ndir, rr, vv)
* 'set gxout contour'
* say vv
* 'set ccolor rainbow'
  'set ccolor 1'
  'set clopts -1 -1 0.2'
  'set clskip 2'
  'set xaxis 0 'rr' 'rr/5
* 'set zlog on'
  'set lon 110 120'
  'set xlopts 1 -1 0.2'
  'set ylopts 1 -1 0.2'
  'set gxout fwrite'
  'set xlabs 0.| |0.5| |1.| |1.5| |2.| |2.5| |3.| |3.5| |4.'
  'd 'vv
  'set gxout contour'
return


** function crs to get a cross section
function crs(var, numcl, num, latc, lonc, dlat, dlon)

*  say 'ii='numcl', num='num', var='var', dlat='dlat', dlon='dlon', latc='latc', lonc='lonc

********
  lon = lonc
  lat = latc
  'collect 'numcl' free'

  ii=0
  while(ii < num)
*   say 'ii,lat,lon= 'ii', 'lat', 'lon
    'collect 'numcl' gr2stn('var','lon','lat')'
    ii=ii+1
    lat=lat+dlat
    lon=lon+dlon
  endwhile

* 'set lon 110 120'
* 'set gxout shaded'
* 'd coll2gr('numcl')'

return




*********************************************************************
** function crs to get a cross section
function crs_ptr(var, varmean, numcl, num, latc, lonc, dlat, dlon)

  lon = lonc
  lat = latc

  'collect 'numcl' free'

  ii=0
  while(ii < num)
*   say 'ii,lat,lon= 'ii', 'lat', 'lon
    'collect 'numcl' gr2stn('var'-'varmean','lon','lat')'
    ii=ii+1
    lat=lat+dlat
    lon=lon+dlon
  endwhile

* 'set lon 110 120'
* 'set gxout shaded'
* 'd coll2gr('numcl')'
return




*********************************************************************
** function crs to get a cross section
function crs_ur(uu, vv, numcl, num, latc, lonc, dlat, dlon, the)

*  say 'ii='numcl', num='num', var='var', dlat='dlat', dlon='dlon', latc='latc', lonc='lonc

********
  lon = lonc
  lat = latc
  'collect 'numcl' free'

  ii=0
  while(ii < num)
*   say 'ii,lat,lon= 'ii', 'lat', 'lon
    'collect 'numcl' gr2stn('uu','lon','lat')*sin('the')+gr2stn('vv','lon','lat')*cos('the')'
    ii=ii+1
    lat=lat+dlat
    lon=lon+dlon
  endwhile

* 'set lon 110 120'
* 'set gxout shaded'
* 'd coll2gr('numcl')'
return





*********************************************************************
** function crs to get a cross section
function crs_ut(uu, vv, numcl, num, latc, lonc, dlat, dlon, the)

*  say 'ii='numcl', num='num', var='var', dlat='dlat', dlon='dlon', latc='latc', lonc='lonc

********
  lon = lonc
  lat = latc
  'collect 'numcl' free'

  ii=0
  while(ii < num)
*   say 'ii,lat,lon= 'ii', 'lat', 'lon
    'collect 'numcl' gr2stn('vv','lon','lat')*sin('the')-gr2stn('uu','lon','lat')*cos('the')'
    ii=ii+1
    lat=lat+dlat
    lon=lon+dlon
  endwhile

* 'set lon 110 120'
* 'set gxout shaded'
* 'd coll2gr('numcl')'

return

******************************************************************
******************************************************************
*** function to get center from file
*** Remark: there are 2 global varables
***         _latc and _lonc
function rdcenter(fname)

*************
**set data time
  'q dim'
  rec=sublin(result,5)
  rtime=subwrd(rec,6)

************************
* get center from file

  ret=read(fname)
  rc =sublin(ret,1)
  data=sublin(ret,2)
* say rc
* say data

  start=substr(data,7,8)
  end=substr(data,16,8)
  name=substr(data,24,15)
 say name' 'start' 'end
***** Read data *****
  ret=read(fname)
  rc=sublin(ret,1)
  data=sublin(ret,2)
  blank=substr(data,13,1)
  nrec=0
  while (blank =' ')
    nrec=nrec+1
    timenum=substr(data,3,10)
    time.nrec=num2time(timenum)
    latt=substr(data,14,6)
    longg=substr(data,21,7)
    if (substr(latt,1,1)=' ')
      latt=substr(latt,2,5)
    endif
    if (substr(latt,6,1)='N')
        latt=substr(latt,1,5)
    endif
    if (substr(longg,1,1)=' ')
      longg=substr(longg,2,6)
    endif
    if (substr(longg,6,1)='E')
      longg=substr(longg,1,6)
    endif
    lati.nrec=latt
    long.nrec=longg
*    spd=substr(data,36,5)
*    if (substr(spd,1,1)=' ')
*      spd=substr(spd,2,4)
*    endif
*    speed.nrec=spd
*   say 'nrec='nrec', time='time.nrec', lat='lati.nrec', long='long.nrec', speed='speed.nrec
*************************************************
    ret=read(fname)
    rc=sublin(ret,1)
    if(rc=2);break;endif
    data=sublin(ret,2)
    blank=substr(data,13,1)
  endwhile

** check center *****************
  ii=1
  while(ii <= nrec)
*   say ii' 'time.ii' 'rtime 
    if(time.ii = rtime)
      _latc=lati.ii
      _lonc=long.ii
      say 'time='rtime', latc='_latc', lonc='_lonc
      break
    endif
    ii=ii+1
  endwhile
  rec=close(fname)
    
return



*******************************************************************
*******************************************************************
function num2month(num)
  num=num*1
  month.1=JAN; month.2=FEB; month.3=MAR; month.4=APR
  month.5=MAY; month.6=JUE; month.7=JUL; month.8=AUG
  month.9=SEP; month.10=OCT; month.11=NOV; month.12=DEC
* say 'month = 'month.num'('num')'
return(month.num)

*******************************************************************
*******************************************************************
function num2time(timenum)
  yr=substr(timenum,1,2)
  mon=substr(timenum,3,2)
  dy=substr(timenum,5,2)
  hr=substr(timenum,7,2)
  mn=substr(timenum,9,2)

  if(yr<10); yr='0'substr(yr,2,1); endif
  if(yr<=20); yr='20'yr; else; yr='19'yr; endif
  mon=num2month(mon)
* say timenum', yr='yr', mon='mon', day='dy', hr='hr
  if(mn='00' | mn='0')
    time=hr'Z'dy%mon%yr
  else
    time=hr%':'%mn%'Z'dy%mon%yr
  endif
return(time)
