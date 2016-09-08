pro sst21
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
; map_world,/over,/continents
help, temp
dum = temp
dumlat = lat
anomp =20

;;;;;;; have anomaly  map
b ='/home/ccr-01/buenning/climssts/blob.nc'
lon2 = nc_get(b, 'lon')
lat2 = nc_get(b, 'lat')
anom = nc_get(b, 'skt')
nlons2= n_elements(lon2)
nlats2= n_elements(lat2)


;;;;;;; place the box in a sea of zeros

firstanom = anom
dumanom = congrid(firstanom, nlons, nlats,/interp)

print, lon
SST= temp
print, max(SST), min(SST)


window, 1
firstanom = firstanom[*,1:nlats2-2]
lat2 = lat2[1:nlats2-2]
contour,  firstanom, lon2, lat2, levels = -0.8 + findgen(17)/16.0 * 1.6 ,/cell_fill
; map_world,/over,/continents
window, 2
 contour,  dumanom, lon, lat, levels = -0.8 + findgen(17)/16.0 * 1.6 ,/cell_fill
; map_world,/over,/continents

;;;;;;;;;flip
newanom = dumanom*0.0
minlon = 0
maxlon = 360
minlat = -85
maxlat = 85

latbox = where(lat ge minlat and lat le maxlat)
lonbox = where(lon ge minlon and lon le maxlon)
jy0 =min(latbox)
jy1 =max(latbox)
ix0 =min(lonbox)
ix1 =max(lonbox)
for i = ix0, ix1  do for j = jy0, jy1 do begin
    newanom[i,j] = dumanom[i,j]
endfor 


window, 3
 contour,  newanom, lon, lat, levels = -0.8 + findgen(17)/16.0 * 1.6 ,/cell_fill
; map_world,/over,/continents

;;;;;;;;; get rid of low correlations
; newanom = newanom * (abs(newanom) gt 0.1)
print, max(newanom), min(newanom)
; newanom = newanom * (-4.0)

SSTold = SST

for k = 0, 11 do begin
   SST[*,*,k] = SST[*,*,k] + newanom
endfor

window, 4
contour, SST[*,*,0], lon, lat, levels = 250+findgen(50),/cell_fill
; map_world,/over,/continents

; normal
contour, SST[*,*,0]-SSTold, lon, lat, levels = -2.0 + findgen(21)/20.0*4.0,/cell_fill
; map_world,/over,/continents
device,/close





filo='sst21.new.nc'
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
