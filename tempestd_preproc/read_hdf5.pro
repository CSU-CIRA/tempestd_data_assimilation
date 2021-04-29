function read_hdf5,filename,vars=vars,verbose=verbose

c=[' ','-','.','/','\',')','(','+',','] ; Define special characters (replaced with _ in output variable names)

nvars = n_elements(vars)
FileID=H5F_OPEN(filename)
nobjs=H5G_GET_NUM_OBJS(fileid)
if (keyword_set(verbose)) then print,'nobjs = ',nobjs
name=strarr(nobjs)
var1=strarr(1000)
var2=strarr(1000)
var=strarr(1000)
nvar=0
for n=0,nobjs-1 do begin
  result=H5G_GET_OBJ_NAME_BY_IDX(fileID,n)
  name[n] = result
  info=H5G_GET_OBJINFO(fileID,name[n])
  if (keyword_set(verbose)) then help,info
  if (info.type EQ 'GROUP') then begin
    groupID=H5G_OPEN(fileID,name[n])
    nmem=H5G_GET_NMEMBERS(fileID,name[n])
    if (keyword_set(verbose)) then print,'nmem = ',nmem
    for k=0,nmem-1 do begin
      result=H5G_GET_MEMBER_NAME(fileID,name[n],k)
      if (keyword_set(verbose)) then print,k,nmem,name[n],result
      var1[k]='/'+name[n]+'/'+result
      vname = name[n]+'_'+result
      for j=0,n_elements(c)-1 do begin
        pos=strpos(vname,c[j])
        while (pos GE 0) do begin
          if (pos LT strlen(vname)-1) then begin
            strput,vname,'_',pos
            pos=strpos(vname,c[j])
          endif else begin
            strput,vname,' ',pos
            pos=strpos(vname,c[j])
          endelse
        endwhile
        pos=strpos(vname,'__')
        while (pos GE 0) do begin
          strput,vname,'_ ',pos
          pos=strpos(vname,'__')
        endwhile
      endfor
      var[nvar] = strcompress(vname,/remove_all)

      if (keyword_set(verbose)) then print,k,var1[k],var[nvar],format='("Var[",i2.2,"] = ",a,a)'
      info=H5G_GET_OBJINFO(fileID,var1[k])
      if (info.type EQ 'GROUP') then begin
        groupID=H5G_OPEN(fileID,var1[k])
        nmem=H5G_GET_NMEMBERS(fileID,var1[k])
        for m=0,nmem-1 do begin
          result=H5G_GET_MEMBER_NAME(fileID,var1[k],m)
          var2[nvar]=var1[k]+'/'+result
          vname = name[n]+'_'+result
          for j=0,n_elements(c)-1 do begin
            pos=strpos(vname,c[j])
            while (pos GE 0) do begin
              if (pos LT strlen(vname)-1) then begin
                strput,vname,'_',pos
                pos=strpos(vname,c[j])
              endif else begin
                strput,vname,' ',pos
                pos=strpos(vname,c[j])
              endelse
            endwhile
            pos=strpos(vname,'__')
            while (pos GE 0) do begin
              strput,vname,'_ ',pos
              pos=strpos(vname,'__')
            endwhile
          endfor
          var[nvar] = strcompress(vname,/remove_all)

          if (keyword_set(verbose)) then print,nvar,var2[nvar],format='("Var[",i2.2,"] = ",a)'
          info=H5G_GET_OBJINFO(fileID,var2[nvar])
          if (info.type EQ 'DATASET') then begin
            if (nvars GT 0) then ind = where(strupcase(vname) EQ strupcase(vars), cnt) else cnt=0
            if (cnt GT 0 OR nvars EQ 0) then begin
              dataID=h5d_open(fileID,var2[nvar])
              data=H5D_READ(dataID)
              if (keyword_set(verbose)) then print,nvar,var[nvar],format='("Var[",i3.3,"] = ",a)'
              if (nvar EQ 0) then begin
                a=create_struct(var[nvar],data)
              endif else begin
                a=create_struct(a,var[nvar],data)
              endelse
              nvar = nvar + 1
            endif
          endif
        endfor
      endif else begin
        if (info.type EQ 'DATASET') then begin
          var2[nvar] = var1[k]
          if (nvars GT 0) then ind = where(strupcase(vname) EQ strupcase(vars), cnt) else cnt=0
          if (cnt GT 0 OR nvars EQ 0) then begin
            dataID=h5d_open(fileID,var2[nvar])
            data=H5D_READ(dataID)
            ;help,data
            if (keyword_set(verbose)) then print,nvar,var[nvar],format='("Var[",i3.3,"] = ",a)'
            if (nvar EQ 0) then begin
              a=create_struct(var[nvar],data)
            endif else begin
              a=create_struct(a,var[nvar],data)
            endelse
            nvar = nvar + 1
          endif
        endif
      endelse
    endfor
  endif else begin
    if (nvars GT 0) then ind = where(strupcase(name[n]) EQ strupcase(vars), cnt) else cnt=0
    if (cnt GT 0 OR nvars EQ 0) then begin
      dataID=h5d_open(fileID,name[n])
      data=H5D_READ(dataID)
      vname = name[n]
      for j=0,n_elements(c)-1 do begin
        pos=strpos(vname,c[j])
        while (pos GE 0) do begin
          if (pos LT strlen(vname)-1) then begin
            strput,vname,'_',pos
            pos=strpos(vname,c[j])
          endif else begin
            strput,vname,' ',pos
            pos=strpos(vname,c[j])
          endelse
        endwhile
        pos=strpos(vname,'__')
        while (pos GE 0) do begin
          strput,vname,'_ ',pos
          pos=strpos(vname,'__')
        endwhile
      endfor
      name[n] = strcompress(vname,/remove_all)
      if (keyword_set(verbose)) then print,nvar,name[n],format='("Var[",i3.3,"] = ",a)'
      if (nvar EQ 0) then begin
        a=create_struct(name[n],data)
      endif else begin
        a=create_struct(a,name[n],data)
      endelse
      nvar = nvar + 1
    endif
  endelse
endfor

;H5F_CLOSE, FileID
return,a
end

