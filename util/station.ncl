load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "./util/gsn_add_text_color.ncl"

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
print("Plotting Stations for the following data:")
print("  " + FILES)

do k=0,numFILES-1

  levels = (/1001,1000,975,950,925,900,875,850,825,800,775,750,725,700,675,650,625,600,\
            575,550,525,500,475,450,425,400,375,350,325,300,275,250,225,200,175,150,125,\
            100,90,80,70,60,50,40,30,20,10,5/)
  nlevels = dimsizes(levels)

  a = addfile(FILES(k),"r")

  datec = stringtochar(FILES(k))
  if(qcOBS.eq."used") then
    date = chartostring(datec(18:30))
  else
    date = chartostring(datec(17:29))
  end if
  print("DATE = " + date)
  stat = a->stationID
  stats = str_squeeze(chartostring(stat(:,1:6)))
  lat = a->lat
  lon = a->lon

  wks = gsn_open_wks(outTYPE,"qc_obs_"+qcOBS+".station." + date)
  gsn_define_colormap(wks,"default")

  p = a->pressure
  t = a->temperature-273.15
  t_qc = a->temperature_QC
  sound = a->sounding
  rh = a->rh
  rh_qc = a->rh_QC
  td = a->dew_point-273.15
  sp = a->speed
  dir = a->direction
  u = a->u*1.9
  u_qc = a->u_QC
  v = a->v*1.9
  v_qc = a->v_QC
  slp = a->slp/100.
  slp_qc = a->slp_QC

  rep = dimsizes(slp)

  res = True
  res@mpFillOn = False
  res@mpDataBaseVersion = "MediumRes"
  ;res@mpOutlineBoundarySets  = "AllBoundaries"
  res@mpOutlineBoundarySets  = "National"
  res@gsnFrame = False
  res@mpGridAndLimbOn = True
;  res@mpGridSpacingF = 45.
  res@pmTickMarkDisplayMode  = "Always"
  res@tmYROn = False
  res@tmXTOn = False

  if(a@MAP_PROJ .eq. 1) then  ;Lambert
    res1 = True
    res1@MAP_PROJ = 1
    res1@TRUELAT1 = a@TRUELAT1
    res1@TRUELAT2 = a@TRUELAT2
    res1@DX = a@DX
    res1@DY = a@DY
    res1@REF_LAT = a@CEN_LAT
    res1@REF_LON = a@CEN_LON
    res1@STAND_LON = a@STAND_LON
    res1@POLE_LAT = 90.0
    res1@POLE_LON = 0.0
    res1@LATINC = 0.0
    res1@LONINC = 0.0
    res1@KNOWNI = (a@i_parent_end)/2
    res1@KNOWNJ = (a@j_parent_end)/2
    loc1 = wrf_ij_to_ll(a@i_parent_start,a@j_parent_start,res1)
    loc2 = wrf_ij_to_ll(a@i_parent_end,a@j_parent_end,res1)

    res@mpLimitMode = "Corners"
    res@mpLeftCornerLatF = loc1(1)
    res@mpLeftCornerLonF = loc1(0)
    res@mpRightCornerLatF = loc2(1)
    res@mpRightCornerLonF = loc2(0)
    res@tfDoNDCOverlay = True
    res@mpProjection = "LambertConformal"
    res@mpLambertParallel1F = a@TRUELAT1
    res@mpLambertParallel2F = a@TRUELAT2
    res@mpLambertMeridianF = a@STAND_LON
  end if

  if(a@MAP_PROJ .eq. 2) then  ;Polar Stereographic
    res1 = True
    res1@MAP_PROJ = 2
    res1@TRUELAT1 = a@TRUELAT1
    res1@TRUELAT2 = a@TRUELAT2
    res1@DX = a@DX
    res1@DY = a@DY
    res1@REF_LAT = a@CEN_LAT
    res1@REF_LON = a@CEN_LON
    res1@STAND_LON = a@STAND_LON
    res1@POLE_LAT = 0.0
    res1@POLE_LON = 0.0
    res1@LATINC = 0.0
    res1@LONINC = 0.0
    res1@KNOWNI = (a@i_parent_end)/2
    res1@KNOWNJ = (a@j_parent_end)/2
    loc1 = wrf_ij_to_ll(a@i_parent_start,a@j_parent_start,res1)
    loc2 = wrf_ij_to_ll(a@i_parent_end,a@j_parent_end,res1)

    res@mpLimitMode = "Corners"
    res@mpLeftCornerLatF = loc1(1)
    res@mpLeftCornerLonF = loc1(0)
    res@mpRightCornerLatF = loc2(1)
    res@mpRightCornerLonF = loc2(0)
    res@tfDoNDCOverlay = True
    res@mpProjection = "Stereographic"
    res@mpCenterLonF = a@STAND_LON
    res@mpCenterLatF = a@CEN_LAT
  end if

  if(a@MAP_PROJ .eq. 3) then  ;Mercator
    res1 = True
    res1@MAP_PROJ = 3
    res1@TRUELAT1 = a@TRUELAT1
    res1@TRUELAT2 = a@TRUELAT2
    res1@DX = a@DX
    res1@DY = a@DY
    res1@REF_LAT = a@CEN_LAT
    res1@REF_LON = a@CEN_LON
    res1@STAND_LON = a@STAND_LON
    res1@POLE_LAT = 90.0
    res1@POLE_LON = 0.0
    res1@LATINC = 0.0
    res1@LONINC = 0.0
    res1@KNOWNI = (a@i_parent_end)/2
    res1@KNOWNJ = (a@j_parent_end)/2
    loc1 = wrf_ij_to_ll(a@i_parent_start,a@j_parent_start,res1)
    loc2 = wrf_ij_to_ll(a@i_parent_end,a@j_parent_end,res1)

    res@mpLimitMode = "Corners"
    res@mpLeftCornerLatF = loc1(1)
    res@mpLeftCornerLonF = loc1(0)
    res@mpRightCornerLatF = loc2(1)
    res@mpRightCornerLonF = loc2(0)
    res@tfDoNDCOverlay = True
    res@mpProjection = "Mercator"
  end if

  txres               = True
  txres@txFontHeightF = 0.01
  txres@txFont        = "helvetica-bold"

  txres_t               = True
  txres_t@txFontHeightF = 0.01
  txres_t@txFont        = "helvetica-bold"

  txres_rh               = True
  txres_rh@txFontHeightF = 0.01
  txres_rh@txFont        = "helvetica-bold"

  txres_slp               = True
  txres_slp@txFontHeightF = 0.01
  txres_slp@txFont        = "helvetica-bold"

  do n=0,nlevels-1

    slpn = new(10000,float)
    tn = new(10000,float)
    rhn = new(10000,float)
    un = new(10000,float)
    vn = new(10000,float)
    latn = new(10000,float)
    lonn = new(10000,float)
    statn = new(10000,string)
    tns = new(10000,string)
    rhns = new(10000,string)
    slpns = new(10000,string)
    colors = new(10000,float)
    colors = 1
    colors@_FillValue = -1
    colors_t = new(10000,float)
    colors_t = -1
    colors_t@_FillValue = -1
    colors_rh = new(10000,float)
    colors_rh = -1
    colors_rh@_FillValue = -1
    colors_slp = new(10000,float)
    colors_slp = -1
    colors_slp@_FillValue = -1

    j = 0

    if(levels(n).eq.1001) then
      do i=0,rep-1
        if(sound(i).eq.0 .and. stats(i).ne."SATOB" .and. stats(i).ne."AIREP" \
           .and. .not.ismissing(p(i,0))) then
          if(.not.ismissing(slp(i)) .or. .not.ismissing(t(i,0)) .or. .not.ismissing(rh(i,0)) \
             .or. .not.ismissing(u(i,0)))
            statn(j) = stats(i)
            latn(j) = lat(i)
            lonn(j) = lon(i)
            slpn(j) = slp(i)
            if(.not.ismissing(slp_qc(i)) .and. slp_qc(i).eq.65536) then
              colors_slp(j) = 2  ;red
            else if(.not.ismissing(slp_qc(i)) .and. slp_qc(i).eq.131072) then
              colors_slp(j) = 3  ;green
            else if(.not.ismissing(slp_qc(i)) .and. slp_qc(i).eq.196608) then
              colors_slp(j) = 4  ;blue
            else
              colors_slp(j) = 1
            end if
            end if
            end if
            tn(j) = t(i,0)
            if(.not.ismissing(t_qc(i,0)) .and. t_qc(i,0).eq.65536) then
              colors_t(j) = 2  ;red
            else if(.not.ismissing(t_qc(i,0)) .and. t_qc(i,0).eq.131072) then
              colors_t(j) = 3  ;green
            else if(.not.ismissing(t_qc(i,0)) .and. t_qc(i,0).eq.196608) then
              colors_t(j) = 4  ;blue
            else
              colors_t(j) = 1
            end if
            end if
            end if
            rhn(j) = rh(i,0)
            if(.not.ismissing(rh_qc(i,0)) .and. rh_qc(i,0).eq.65536) then
              colors_rh(j) = 2  ;red
            else if(.not.ismissing(rh_qc(i,0)) .and. rh_qc(i,0).eq.131072) then
              colors_rh(j) = 3  ;green
            else if(.not.ismissing(rh_qc(i,0)) .and. rh_qc(i,0).eq.196608) then
              colors_rh(j) = 4  ;blue
            else
              colors_rh(j) = 1
            end if
            end if
            end if
            un(j) = u(i,0)
            vn(j) = v(i,0)
            j=j+1
          end if
        end if
      end do

      tns = sprintf("%5.1f",tn)
      rhns = sprintf("%5.1f",rhn)
      slpns = sprintf("%5.0f",slpn)

      indt = ind(.not.ismissing(tn))
      indrh = ind(.not.ismissing(rhn))
      indslp = ind(.not.ismissing(slpn))
      indstat = ind(.not.ismissing(statn))

      dims_slp = dimsizes(indslp)
      colors2_slp = new (dims_slp,float)
      ii = 0
      do jj=0,j
         if ( .not. ismissing(slpn(jj)) ) then
           colors2_slp(ii) = colors_slp(jj)
           ii = ii + 1
         end if
      end do

      dims_t = dimsizes(indt)
      colors2_t = new (dims_t,float)
      ii = 0
      do jj=0,j
         if ( .not. ismissing(tn(jj)) ) then
           colors2_t(ii) = colors_t(jj)
           ii = ii + 1
         end if
      end do

      dims_rh = dimsizes(indrh)
      colors2_rh = new(dims_rh,float)
      ii = 0
      do jj=0,j
        if(.not.ismissing(rhn(jj))) then
          colors2_rh(ii) = colors_rh(jj)
          ii = ii + 1
        end if
      end do

      dims_stat = dimsizes(indstat)
      colors2 = new(dims_stat,float)
      ii = 0
      do jj=0,j
        if(.not.ismissing(statn(jj))) then
          colors2(ii) = colors(jj)
          ii = ii + 1
        end if
      end do

      txres@txFontColors = colors2
      txres@txJust = "TopLeft"
      txres_slp@txFontColors = colors2_slp
      txres_slp@txJust = "BottomLeft"
      txres_t@txFontColors = colors2_t
      txres_t@txJust = "BottomRight"
      txres_rh@txFontColors = colors2_rh
      txres_rh@txJust = "TopRight"

      res@tiMainString = "Level: " + levels(n) + "   Records: " + j + ""
      res@gsnLeftString = "Date: " + date +""

      if (j.ne.0) then
        map = gsn_csm_map(wks,res) ; create map

        print("LEVEL = " + levels(n))

        text1 = gsn_add_text(wks,map," "+"~S10~"+tns(indt) + " ",lonn(indt),latn(indt),txres_t)
        text2 = gsn_add_text(wks,map," "+"~B10~"+rhns(indrh) + " ",lonn(indrh),latn(indrh),txres_rh)
        text3 = gsn_add_text(wks,map," " + "~S10~"+slpns(indslp),lonn(indslp),latn(indslp),txres_slp)
        text4 = gsn_add_text(wks,map," " + "~B10~"+statn(indstat),lonn(indstat),latn(indstat),txres)
      end if

    else

      do i=0,rep-1
        pind = ind(p(i,:).eq.stringtofloat(levels(n)+"00."))
        if(.not.ismissing(pind)) then
          if(.not.ismissing(t(i,pind)) .or. .not.ismissing(rh(i,pind)) \
             .or. .not.ismissing(u(i,pind)))
            statn(j) = stats(i)
            latn(j) = lat(i)
            lonn(j) = lon(i)
            tn(j) = t(i,pind)
            ;print(stats(i) + " " + t_qc(i,pind) + " " + t(i,pind))
            if(.not.ismissing(t_qc(i,pind)) .and. t_qc(i,pind).eq.65536) then
              colors_t(j) = 2  ;red
            else if(.not.ismissing(t_qc(i,pind)) .and. t_qc(i,pind).eq.131072) then
              colors_t(j) = 3  ;green
            else if(.not.ismissing(t_qc(i,pind)) .and. t_qc(i,pind).eq.196608) then
              colors_t(j) = 4  ;blue
            else
              colors_t(j) = 1
            end if
            end if
            end if
            rhn(j) = rh(i,pind)
            if(.not.ismissing(rh_qc(i,pind)) .and. rh_qc(i,pind).eq.65536) then
              colors_rh(j) = 2  ;red
            else if(.not.ismissing(rh_qc(i,pind)) .and. rh_qc(i,pind).eq.131072) then
              colors_rh(j) = 3  ;green
            else if(.not.ismissing(rh_qc(i,pind)) .and. rh_qc(i,pind).eq.196608) then
              colors_rh(j) = 4  ;blue
            else
              colors_rh(j) = 1
            end if
            end if
            end if
            un(j) = u(i,pind)
            vn(j) = v(i,pind)
            j=j+1
          end if
        end if
      end do

      tns = sprintf("%5.1f ",tn)
      rhns = sprintf("%5.1f ",rhn)

      indt = ind(.not.ismissing(tn))
      indrh = ind(.not.ismissing(rhn))
      indstat = ind(.not.ismissing(statn))

      dims_t = dimsizes(indt)
      colors2_t = new (dims_t,float)
      ii = 0
      do jj=0,j
         if ( .not. ismissing(tn(jj)) ) then
           colors2_t(ii) = colors_t(jj)
           ii = ii + 1
         end if 
      end do

      dims_rh = dimsizes(indrh)
      colors2_rh = new(dims_rh,float)
      ii = 0
      do jj=0,j
        if(.not.ismissing(rhn(jj))) then
          colors2_rh(ii) = colors_rh(jj)
          ii = ii + 1
        end if
      end do

      dims_stat = dimsizes(indstat)
      colors2 = new(dims_stat,float)
      ii = 0
      do jj=0,j
        if(.not.ismissing(statn(jj))) then
          colors2(ii) = colors(jj)
          ii = ii + 1
        end if
      end do

      txres@txFontColors = colors2
      txres_t@txFontColors = colors2_t
      txres_rh@txFontColors = colors2_rh

      res@tiMainString = "Level: " + levels(n) + "   Records: " + j + ""

      if (j.ne.0) then
        print("LEVEL = " + levels(n))

        map = gsn_csm_map(wks,res) ; create map

        if (.not.ismissing(indt(0))) then
          text1 = gsn_add_text(wks,map," "+"~S10~"+tns(indt) + " ",lonn(indt),latn(indt),txres_t)
        end if
        if(.not.ismissing(indrh(0))) then
          text2 = gsn_add_text(wks,map," "+"~B10~"+rhns(indrh) + " ",lonn(indrh),latn(indrh),txres_rh)
        end if
        if(.not.ismissing(indstat(0))) then
          text4 = gsn_add_text(wks,map," " + "~B10~"+statn(indstat),lonn(indstat),latn(indstat),txres)
        end if

      end if

    end if

    if(j.ne.0) then

      wmsetp("col",1)
      wmsetp("wbs",0.0315)
      wmbarbmap(wks,latn,lonn,un,vn)

      legend = create "Legend" legendClass wks
       "vpXF" : 0.2
       "vpYF" : 0.165
       "vpWidthF" : 0.4 ; width
       "vpHeightF" : 0.09 ; height
       "lgPerimOn" : False ; perimeter
       "lgItemCount" : 3 ; how many
       "lgLabelStrings" : (/"Failed both errmax and buddy test","Failed buddy test","Failed errmax test"/)
       "lgMonoDashIndex" : True
       "lgLineColors" : (/4,3,2/)
       "lgBoxMinorExtentF" : 0.2
      end create

      draw(map)
      draw(legend)
      frame(wks)

    end if

    delete(statn)
    delete(latn)
    delete(lonn)
    delete(slpn)
    delete(tn)
    delete(rhn)
    delete(un)
    delete(vn)
    delete(slpns)
    delete(tns)
    delete(rhns)
    delete(indt)
    delete(indrh)
    delete(indstat)
    delete(colors2)
    delete(colors2_t)
    delete(colors2_rh)
    delete(txres@txFontColors)
    delete(txres_t@txFontColors)
    delete(txres_rh@txFontColors)
    if(isvar((/"text1"/))) then
      delete(text1)
    end if
    if(isvar((/"text2"/))) then
      delete(text2)
    end if
    if(isvar((/"text4"/))) then
      delete(text4)
    end if

  end do

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
  delete(t_qc)
  delete(rh)
  delete(rh_qc)
  delete(sound)
  delete(td)
  delete(sp)
  delete(dir)
  delete(u)
  delete(u_qc)
  delete(v)
  delete(v_qc)
  delete(slp)
  delete(slp_qc)
  delete(rep)
  delete(indslp)
  delete(txres_slp@txFontColors)
  delete(text3)
  delete(colors2_slp)

end do

end
