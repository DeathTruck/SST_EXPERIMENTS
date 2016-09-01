pro sst5
!x.style=1
!y.style =1 
a='sst.nc'
temp = nc_get(a,'var11')
lon = nc_get(a,'lon')
lat = nc_get(a,'lat')
nlons=n_elements(lon)
nlats=n_elements(lat)
print, lon, lat
time = nc_get(a,'time')
contour, temp[*,*,0], lon, lat, levels = 250+findgen(50),/cell_fill
map_world,/over,/continents
help, temp
dum = temp
dumlat = lat
anomp =20

; center of anomaly
;ilon = 190.0
;jlat = 45
;dumlon = abs(lon-ilon)
;dumlat = abs(lat-jlat)
;ix = where(dumlon eq min(dumlon))
;jy = where(dumlat eq min(dumlat))
;anom= reform(temp[*,*,0])* 0.0
;;;;;;; have correlation map
restore, '/home/ccr-01/buenning/climssts/ensemble_cors01.sav'
;;;;;;; place the box in a sea of zeros
nlons2 = n_elements(lon2)
nlats2 = n_elements(lat2)
anoms = fltarr(nlons2,nlats2)
lonbox = where(lons eq lon2)
latbox = where(lats eq lat2)
nlons3 = n_elements(lons)
nlats3 = n_elements(lats)
lonbox3 = intarr(nlons3)
latbox3 = intarr(nlats3)
for i = 0, nlons3 -1 do lonbox3[i] = where(lons[i] eq lon2)
for i = 0, nlats3 -1 do latbox3[i] = where(lats[i] eq lat2)

firstanom = fltarr(nlons2,nlats2)
for i = 0, nlons3-1 do for j= 0, nlats3-1 do begin
   firstanom[lonbox3[i],latbox3[j]] = image2[i,j]
endfor

dumanom = fltarr(nlons2,nlats2-2,2)
dumanom[*,*,0] = firstanom[*,1:nlats2-2]
dumanom[*,*,1] = firstanom[*,1:nlats2-2]
;dumanom[*,*,2] = firstanom[*,1:nlats2-2]
;dumanom[*,*,3] = firstanom[*,1:nlats2-2]
;dumanom[*,*,4] = firstanom[*,1:nlats2-2]
;dumanom[*,*,5] = firstanom[*,1:nlats2-2]
;dumanom[*,*,6] = firstanom[*,1:nlats2-2]
;dumanom[*,*,7] = firstanom[*,1:nlats2-2]
;dumanom[*,*,8] = firstanom[*,1:nlats2-2]
;dumanom[*,*,9] = firstanom[*,1:nlats2-2]
;dumanom[*,*,10] = firstanom[*,1:nlats2-2]
;dumanom[*,*,11] = firstanom[*,1:nlats2-2]


; lat2 = lat2[1:nlats2-2]
; regrid,dumanom,lon2,lat2,nlons,nlats,fnew,lonnew=lonnew,latnew=latnew, wrap=1, gauss = 0

dumanom = congrid(firstanom, nlons, nlats,/interp)

print, max(image2), min(image2)

   




; for i = 0,nlons-1 do for j= 0, nlats-1 do begin
;    if abs(ix - i) lt 6 and abs(jy - j) lt 6 then begin
;      anom[i,j] = anomp; /(abs(ix-i) + abs(jy-j))
;    endif
; endfor

;;;;;;;;;;;;;;;apply anomaly
; for i = 0,nlons-1 do for j= 0, nlats-1 do for k = 0, 11 do begin
;    
;    temp[i,j,k]=temp[i,j,k] - anom[i,j]
; 
; endfor
; dum = temp

;lon = [lon,360+lon]
;print, lon
;newtemp = fltarr(nlons*2,nlats,12)
;newtemp[0:95,*,*] = temp
;newtemp[96:191,*,*]=temp
;temp=newtemp
; temp=shift(temp,48,0,0)
 ; lon=shift(lon,48)
 ; temp2=temp[48:48+95,*,*]
 ; temp = temp[0:95,*,*]
 ; lon2= lon[0:95]
 ; lon= lon[48:48+95]
; lon[48:95]=360+lon[48:95]
print, lon
SST= temp
print, max(SST), min(SST)
;;;;;;;;;;;:::::REGRIDDING
; SST=congrid(temp, 360,180,12,/interp,/center)
;  nx = 360		; number of longitudes
;  ny = 180		; number of latitudes

;  regrid,temp,lon,lat,nx,ny,fnew,lonnew=lonnew,latnew=latnew, wrap=1
;SST=fnew
;  regrid,temp2,lon,lat,nx,ny,fnew,lonnew=lonnew2,latnew=latnew2, wrap=1
;SST[354,*,*]=fnew[174,*,*]
;SST[355,*,*]=fnew[175,*,*]
;SST[356,*,*]=fnew[176,*,*]
;SST[357,*,*]=fnew[177,*,*]
;SST[358,*,*]=fnew[178,*,*]
;SST[359,*,*]=fnew[179,*,*]
;SST[0,*,*]=fnew[180,*,*]
;SST[1,*,*]=fnew[181,*,*]
  ; SST=temp2  
;  lon =0.5+findgen(360)
; lat=-89.5+findgen(180)



window, 1
firstanom = firstanom[*,1:nlats2-2]
lat2 = lat2[1:nlats2-2]
contour,  firstanom, lon2, lat2, levels = -0.8 + findgen(17)/16.0 * 1.6 ,/cell_fill
map_world,/over,/continents
window, 2
 contour,  dumanom, lon, lat, levels = -0.8 + findgen(17)/16.0 * 1.6 ,/cell_fill
map_world,/over,/continents

;;;;;;;;;flip
newanom = dumanom
for i = 0, nlons-1  do for j = 0, nlats-1 do begin
   newanom[i,j] = dumanom[i,nlats-j-1]
endfor 
window, 3
 contour,  newanom, lon, lat, levels = -0.8 + findgen(17)/16.0 * 1.6 ,/cell_fill
map_world,/over,/continents

;;;;;;;;; get rid of low correlations
; newanom = newanom * (abs(newanom) gt 0.1)
print, max(newanom), min(newanom)
newanom = newanom * (-4.0)

SSTold = SST

for k = 0, 11 do begin
   SST[*,*,k] = SST[*,*,k] + newanom
endfor

window, 4
contour, SST[*,*,0], lon, lat, levels = 250+findgen(50),/cell_fill
map_world,/over,/continents

normal
contour, SST[*,*,0]-SSTold, lon, lat, levels = -2.0 + findgen(21)/20.0*4.0,/cell_fill
map_world,/over,/continents
device,/close
stop





filo='sst5.new.nc'
id=ncdf_create(filo,/clobber)
zid=ncdf_dimdef(id,'time',12)
yid=ncdf_dimdef(id,'lat',n_elements(lat))
xid=ncdf_dimdef(id,'lon',n_elements(lon))
;yid=ncdf_dimdef(id,'lat',180)
;xid=ncdf_dimdef(id,'lon',360)
timeid=ncdf_vardef(id,'time',[zid])
latid=ncdf_vardef(id,'lat',[yid])
lonid=ncdf_vardef(id,'lon',[xid])
vid=ncdf_vardef(id,'SST',[xid,yid,zid], /FLOAT)
ncdf_control, id,/endef
ncdf_varput, id, vid, SST
ncdf_varput, id, timeid, time
ncdf_varput, id, lonid, lon
ncdf_varput, id, latid, lat
ncdf_close, id
stop

end
