function convert_date,idate,base_year=base_year

julian = [[0,31,59,90,120,151,181,212,243,273,304,334], $
          [0,31,60,91,121,152,182,213,244,274,305,335]]
smon=['January','February','March','April','May','June','July', $
      'August','September','October','November','December']
bmon=bytarr(12,9)
for i=0,11 do begin
  for j=0,strlen(smon[i])-1 do begin
    bmon[i,j] = byte(strmid(smon[i],j,1))
  endfor
endfor
date = 0L

if (n_elements(base_year) EQ 0) then base_year=1970
nrec=n_elements(idate)
st=create_struct('Year',              0, $
                 'Month',             0, $
                 'MonthStr', bytarr(10), $
                 'Day',               0, $
                 'JulDay',            0, $
                 'Date',             0L, $
                 'JulDate',          0L, $
                 'ITime',            0L, $
                 'Time',            0.0)
if (nrec GT 1) then begin
  a=replicate(st,nrec)
endif else begin
  a=st
endelse

for n=0,nrec-1 do begin
  if (idate[n] LT 10000000) then begin
    year = idate[n] / 1000
    iyr  = year
    julday = idate[n] MOD 1000
;    if (year GT 50) then begin
;        year = year + 1900
;    endif else begin
;        year = year + 2000
;    endelse
    if (year MOD 4 EQ 0) then ind=1 else ind=0
    imon = 12
    for i=1, 12 do begin
       if (julian[i-1,ind] GT julday-1) then begin
          imon = i-1
          break
       endif
    endfor
    nleap_year = (year - 1) / 4 - (base_year - 1) / 4
    nday = 366 * nleap_year + 365 * (year - base_year - nleap_year) + julday - 1
    itime = 86400L * nday
    iday = julday - julian[imon-1,ind]
    a[n].Year      = year
    a[n].Month     = imon
    a[n].Day       = iday
    a[n].JulDay    = julday
    a[n].MonthStr  = bmon[imon-1,*]
    a[n].Date      = iyr * 10000 + imon * 100 + iday
    a[n].JulDate   = iyr*1000 + julday
    if (ind EQ 0) then begin
      a[n].Time      = float(year) + float(julday) / 365.0
    endif else begin
      a[n].Time      = float(year) + float(julday) / 366.0
    endelse
    a[n].ITime       = itime
  endif else begin
    year = idate[n] / 10000
    iyr  = year
    imon = (idate[n] - year * 10000) / 100
    iday = (idate[n] - year * 10000 - imon * 100)

;    if (year GT 50) then begin
;      year = year + 1900
;    endif else begin
;      year = year + 2000
;    endelse
    if (year MOD 4 EQ 0) then ind=1 else ind=0
    julday = julian[imon-1,ind] + iday
    nleap_year = (year - 1) / 4 - (base_year - 1) / 4
    nday = 366 * nleap_year + 365 * (year - base_year - nleap_year) + julday - 1
    itime = 86400L * nday
    a[n].Year      = year
    a[n].Month     = imon
    a[n].Day       = iday
    a[n].JulDay    = julday
    a[n].MonthStr  = bmon[imon-1,*]
    a[n].Date      = iyr * 10000 + imon * 100 + iday
    a[n].JulDate   = iyr*1000 + julday
    if (ind EQ 0) then begin
      a[n].time      = year + julday / 365.0
    endif else begin
      a[n].time      = year + julday / 366.0
    endelse
    a[n].ITime       = itime
  endelse
endfor
return,a

end
