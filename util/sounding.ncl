load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

if (.not.isvar("qcOBS")) then
  qcOBS = "used"
end if
if (.not.isvar("YYYY")) then
  YYYY="*"
end if
if (.not.isvar("MM")) then
  MM="*"
end if
if (.not.isvar("DD")) then
  DD="*"
end if
if (.not.isvar("HH")) then
  HH="*"
end if
if (.not. isvar("outTYPE")) then
  outTYPE = "x11"
end if
DATE = YYYY+"-*"+MM+"-*"+DD+"_*"+HH+"*"

DATADir = "./"
FILES = systemfunc (" ls -1 " + DATADir + "qc_obs_"+qcOBS+"*"+DATE+".nc ")
numFILES = dimsizes(FILES)
print("Plotting Sounding for the following data:")
print("  " + FILES)

do k=0,numFILES-1

  a = addfile(FILES(k),"r")

  datec = stringtochar(FILES(k))
  if(qcOBS.eq."used") then
    date = chartostring(datec(18:30))
  else
    date = chartostring(datec(17:29))
  end if
  print("DATE = " + date)

  wks = gsn_open_wks(outTYPE,"qc_obs_"+qcOBS+".sounding." + date)

  stat = a->stationID
  stats = str_squeeze(chartostring(stat(:,1:6)))
  lat = a->lat
  lon = a->lon

  p = a->pressure/100.
  t = a->temperature-273.15
  sound = a->sounding
  rh = a->rh
  td = a->dew_point-273.15
  sp = a->speed
  dir = a->direction
  u = a->u*1.9
  v = a->v*1.9
  slp = a->slp/100.
  z = a->height

  rep = dimsizes(slp)

  indsound = ind(sound.eq.1 .and. stats.ne."AIREP" .and. stats.ne."SATOB")
  ;print(indsound + " " + stats(indsound))

  if(.not.ismissing(indsound(0))) then
  statn = stats(indsound)
  latn = lat(indsound)
  lonn = lon(indsound)

  nstat = dimsizes(indsound)

  do n=0,nstat-1

    do i=0,rep-1
      if(sound(i).eq.1 .and. stats(i).eq.statn(n)) then
        tn = t(i,:)
        tdn = td(i,:)
        pn = p(i,:)
        zn = z(i,:)
        un = u(i,:)
        vn = v(i,:)
      end if
    end do

    skewtOpts                 = True
    skewtOpts@DrawColAreaFill = True    ; default is False
    skewtOpts@DrawFahrenheit  = False
    skewtOpts@tiMainString    = "Station = " + statn(n) + "  Lat = " + latn(n) + "  Lon = " + lonn(n) + "  Date = " + date
    skewtOpts@tiMainFontHeightF = 0.018

    dataOpts = True
    dataOpts@Parcel = 1
    dataOpts@WspdWdir  = False  ; wind speed and dir [else: u,v]
    dataOpts@HspdHdir  = True   ; wind speed and dir [else: u,v]

    if (.not.all(ismissing(tn)) .and. pn(0).gt.500.) then
      print("STATION = " + statn(n))
      skewt_bkgd = skewT_BackGround (wks,skewtOpts)
      skewt_data = skewT_PlotData(wks,skewt_bkgd,pn,tn,tdn,zn,un,vn,dataOpts)
      draw (skewt_bkgd)
      draw (skewt_data)
      frame(wks)
    end if

  end do

  delete(statn)
  delete(latn)
  delete(lonn)
  delete(nstat)
  delete(tn)
  delete(tdn)
  delete(pn)
  delete(zn)
  delete(un)
  delete(vn)
  delete(skewtOpts)
  delete(dataOpts)

  end if

  delete(a)
  delete(wks)
  delete(datec)
  delete(date)
  delete(stat)
  delete(stats)
  delete(lat)
  delete(lon)
  delete(p)
  delete(t)
  delete(sound)
  delete(rh)
  delete(td)
  delete(sp)
  delete(dir)
  delete(u)
  delete(v)
  delete(slp)
  delete(z)
  delete(rep)
  delete(indsound)

end do

end
