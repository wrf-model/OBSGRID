;***********************************************************************;
; Function : gsn_add_text                                               ;
;                   wks: workstation object                             ;
;                plotid: plot object                                    ;
;                  text: array of text strings                          ;
;                     x: n-dimensional array of x data positions        ;
;                     y: n-dimensional array of y data positions        ;
;               resources: optional resources                           ;
;                                                                       ;
; This function adds text strings to the plot "plotid". "x" and "y" are ;
; the x and y locations of each text string, and should be specified in ;
; the same data space as the data space of "plotid". "resources" is an  ;
; optional list of TextItem and AnnoManager resources.                  ;
;                                                                       ;
; This function is different from gsn_text because it actually attaches ;
; the text to the plot. This means that if you resize or move the plot, ;
; the text will stay with the plot.                                     ;
;                                                                       ;
; This special version checks for "txFontColors", which is not a real   ;
; NCL resource. It was added for purposes of allowing multiple colors   ;
; for text strings.                                                     ;
;***********************************************************************;
undef("gsn_add_text")
function gsn_add_text(wks:graphic,plotid:graphic,texto:string, \
                      xo:numeric,yo:numeric, resources:logical )
local txid, txres, amres, just, res2, wksname, am_ids
begin
;
; Any one of xo, yo, and texto can just be one element, but if two or more
; are more than one element, then they must be exactly the same size.
;
  xsizes = dimsizes(xo)
  ysizes = dimsizes(yo)
  tsizes = dimsizes(texto)
  xrank  = dimsizes(xsizes)
  yrank  = dimsizes(ysizes)
  trank  = dimsizes(tsizes)
  if(xrank.gt.1.and.yrank.gt.1.and..not.all(xsizes.eq.ysizes)) then
    print("Error: gsn_add_text: x and y must have the same dimension sizes, or either be a single value.")
    dummy = new(1,graphic)
    return(dummy)
  end if
  if(trank.gt.1.and. \
     (xrank.gt.1.and..not.all(xsizes.eq.tsizes)) .or. \
     (yrank.gt.1.and..not.all(ysizes.eq.tsizes))) then
    print("Error: gsn_add_text: text must be a single string or the same dimension size as x and/or y.")
    dummy = new(1,graphic)
    return(dummy)
  end if
;
; Convert to 1-dimensional arrays of all the same length.
;
  if(xrank.gt.1) then
    x    = ndtooned(new(xsizes, typeof(xo)))
    y    = ndtooned(new(xsizes, typeof(yo)))
    text = ndtooned(new(xsizes, typeof(texto)))
  else
    if(yrank.gt.1) then
      x    = ndtooned(new(ysizes, typeof(xo)))
      y    = ndtooned(new(ysizes, typeof(yo)))
      text = ndtooned(new(ysizes, typeof(texto)))
    else
      x    = new(xsizes > ysizes, typeof(xo))
      y    = new(xsizes > ysizes, typeof(yo))
      text = new(xsizes > ysizes, typeof(texto))
    end if
  end if

  x    = ndtooned(xo)
  y    = ndtooned(yo)
  text = ndtooned(texto)
  len  = dimsizes(x)

  res2 = get_resources(resources)
;
; The "txFuncCode" can't be set during a setvalues call. It must be
; set during the creation of the object.  
;
  wksname = get_res_value_keep(wks,"name","gsnapp")

  txres = get_res_eq(res2,"tx")  ; Get text resources.
  txid  = new(len,graphic)

;---Check for multiple colors
  use_mult_colors = False
  if(isatt(txres,"txFontColors")) then
    colors  = get_res_value(txres,"txFontColors",1)
    ncolors = dimsizes(colors)
    if(ncolors.eq.len) then
      use_mult_colors = True
    else
      color = colors(0)    ; This is in case we only have one color,
                           ; or not enough colors
      print("Warning: gsn_add_text: the wrong number of colors was input via txFontColors.")
      print("         The default foreground color will be used.")
    end if
  end if

  if(res2.and.isatt(res2,"txFuncCode")) then
    do i=0,len-1
      if(use_mult_colors) then
        color = colors(i)
      end if
      txid(i) = create wksname + "_text"+i textItemClass wks
        "txString"    : text(i)
        "txFuncCode"  : res2@txFuncCode
        "txFontColor" : color
      end create
      attsetvalues_check(txid(i),txres)          ; Set text resources.
    end do
  else
    do i=0,len-1
      if(use_mult_colors) then
        color = colors(i)
      end if
      txid(i) = create wksname + "_text"+i textItemClass wks
        "txString"    : text(i)
        "txFontColor" : color
      end create
      attsetvalues_check(txid(i),txres)          ; Set text resources.
    end do
  end if
;
; Get current list of annotations that are already attached to
; the plot.
;
  getvalues plotid
    "pmAnnoViews" : text_ids
  end getvalues
;
; Make sure the next text strings are first in the list.
;
  if(.not.any(ismissing(text_ids)))
    new_text_ids            = new(dimsizes(text_ids)+len,graphic)
    new_text_ids(0:len-1)  = txid
    new_text_ids(len:)     = text_ids
  else
    new_text_ids = txid
  end if
;
; Set the old and new annotations, with the new ones being first.
;
  setvalues plotid
    "pmAnnoViews" : new_text_ids
  end setvalues
;
; Retrieve the id of the AnnoManager object created by the PlotManager and
; then set its location in data coordinate space.
;
  getvalues plotid
    "pmAnnoManagers": am_ids
  end getvalues

  tmp_just  = get_res_value(txres,"txJust","CenterCenter")
  just      = get_res_value(res2,"amJust",tmp_just)

  do i=0,len-1
    setvalues am_ids(i)
      "amDataXF"       : x(i)
      "amDataYF"       : y(i)
      "amResizeNotify" : True
      "amTrackData"    : True
      "amJust"         : just
    end setvalues
  end do

  amres = get_res_eq(res2,"am")           ; Get annomanager resources.
  attsetvalues_check(am_ids(0),amres)     ; Set annomanager resources.

  if(xrank.gt.1) then
    return(onedtond(am_ids(0:len-1),xsizes))
  else
    return(onedtond(am_ids(0:len-1),ysizes))
  end if
end


