import rasterio
import os
import csv
import sys
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from colormaps import*
from osgeo import gdal_array,gdal,gdalconst
from numpy import where,linspace,meshgrid,nan,isnan,ma,log,zeros


def cell_to_coord(row, col, cellx, celly, xmin, ymax):
	lon = cellx*col + xmin
	lat = ymax-celly*row
	return lat,lon


####make raster 0.5 x 0.5 for regions
#from osgeo import gdal_array,gdal,gdalconst,osr,ogr
#from numpy import ceil
#srs = osr.SpatialReference()
#srs.ImportFromEPSG(4326)
#
#vector_fn = "./MEOW/meow_ecos.shp"
### Open the data source and read in the extent
#source_ds = ogr.Open(vector_fn)
#source_layer = source_ds.GetLayer()
#x_min =-180.0
#x_max = 180.0
#y_min = -90.0
#y_max = 90.0
#NoData_value = 0
#pixel_size = 0.5
#x_res = int(ceil((x_max - x_min) / pixel_size))
#y_res = int(ceil((y_max - y_min) / pixel_size))
#raster_fn = './MEOW/meow_05.tif'
#target_ds = gdal.GetDriverByName('GTiff').Create(raster_fn, x_res, y_res, 1, gdal.GDT_Byte, options=["COMPRESS=LZW"])
#target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
#target_ds.SetProjection(srs.ExportToWkt())
#band = target_ds.GetRasterBand(1)
#band.SetNoDataValue(NoData_value)
#target_ds.GetMetadata()
#gdal.RasterizeLayer(target_ds, [1], source_layer, options=['ALL_TOUCHED=TRUE',"ATTRIBUTE=RLM_CODE"])
#target_ds = None
####


meow_fields = csv.reader(open('./MEOW/meow_ecos.csv','r'),delimiter=';')
fields = list(enumerate(next(meow_fields)))
#rlm_code = 4
#rlm = 5
reg_names = dict([])
for i in meow_fields:
	reg_names[int(i[4])] = i[5]



csv.field_size_limit(sys.maxsize)

assis = csv.reader(open('./Macrophytes/database_pruned_assis_et_al.csv','r'),delimiter=';')
head = list(enumerate(next(assis)))

#species = 101
#kingdom = 102
mac_alg = dict([])
for i in assis:
	try:
		if i[102]!='':
			mac_alg[i[101].replace(' ','_')] = i[102]
	except:
		pass




reg_rast = rasterio.open('./MEOW/meow_05.tif').read(1)
regions = [1*(reg_rast==i) for i in range(1,13)]

l_mask = rasterio.open('./light_at_bottom/Present_Benthic_Max_Depth_Light_bottom_Max_05.tif').read(1)
l_mask = 1*(l_mask>50) #check the threshold

area = rasterio.open('land_area_05.tif').read(1)

tot_area_light_ok = (l_mask*area).sum()
###maybe useful, this is the total area where light in principle permits macrophytes' existence

spp = sorted([i[:-4] for i in os.listdir('./macrophyte_ranges_05/')])

out = open('range_changes.csv','w')
out.write('ssp,region,kingdom,species,year,loss,suit,prop\n')

for sce in ['245','370','585']:
	fut_dir = './macrophyte_future_suit_05/ssp'+sce+'/'
	ratio_occ = dict([]) #occupied vs. potential in the present
	pres_range = dict([]) #present range
	for sp in spp:
		ahull = rasterio.open('./macrophyte_ranges_05/'+sp+'.tif').read(1)
		suit = rasterio.open(fut_dir+sp+'_2015.tif').read(1)*l_mask
		pres_range[sp] = ahull*l_mask #EOO
		ratio_occ[sp] = (ahull*suit*area).sum()/(suit*area).sum()
	for year in range(2015,2101):
		for sp in spp:
			suit = rasterio.open(fut_dir+sp+'_'+str(year)+'.tif').read(1)
			range_loss = pres_range[sp]*suit*area
			range_suit = suit*l_mask*area
			range_prop = range_suit*ratio_occ[sp]
			out.write(','.join(map(str,[sce,'global',mac_alg[sp],sp,year,range_loss.sum(),range_suit.sum(),range_prop.sum()]))+'\n')
			for reg_n in range(12):
				out.write(','.join(map(str,[sce,reg_names[reg_n+1],mac_alg[sp],sp,year,
							 (regions[reg_n]*range_loss).sum(),
							 (regions[reg_n]*range_suit).sum(),
							 (regions[reg_n]*range_prop).sum()]))+'\n')
	#	print (sce,year)


out.close()




#####HABITAT
from numpy import arange
out = open('hab_changes.csv','w')
out.write('ssp,region,kingdom,year,hab_ext_full,'+','.join(map(str,arange(0,100,10)))+'\n')
for sce in ['245','370','585']:
	for taxon in ['Chromista','Plantae']:
		fut_dir = './macrophyte_future_suit_05_hab/ssp'+sce+'/'
		for year in range(2015,2101):
			suit = rasterio.open(fut_dir+taxon+'_'+str(year)+'.tif').read(1)
			range_tre = []
			range_suit = suit*l_mask*area
			for tre in arange(0,1,0.1):
				range_tre.append(1*(suit>tre)*l_mask*area)
			out.write(','.join(map(str,[sce,'global',taxon,year,range_suit.sum()]+[i.sum() for i in range_tre]))+'\n')
			for reg_n in range(12):
				out.write(','.join(map(str,[sce,reg_names[reg_n+1],taxon,year,
							 (regions[reg_n]*range_suit).sum()]+[(regions[reg_n]*i).sum() for i in range_tre]))+'\n')
		#	print (sce,year,taxon)


out.close()
########END HABITAT




###Generate rasters for maps
rst_fn = 'ocean_05.tif'
rst = rasterio.open(rst_fn)
meta = rst.meta.copy()
meta['dtype']='float64'


for sce in ['245','370','585']:
	fut_dir = './macrophyte_future_suit_05/ssp'+sce+'/'
	pres_range = dict([]) #present range
	for sp in spp:
		ahull = rasterio.open('./macrophyte_ranges_05/'+sp+'.tif').read(1)
		suit = rasterio.open(fut_dir+sp+'_2015.tif').read(1)*l_mask
		pres_range[sp] = ahull*suit
	for year in range(2015,2101):
		for taxon in ['Chromista', 'Plantae']:
			mat_loss = zeros([360,720])
			mat_suit = zeros([360,720])
			for sp in spp:
				if mac_alg[sp]==taxon:
					suit = rasterio.open(fut_dir+sp+'_'+str(year)+'.tif').read(1)
					pres_range[sp] = 1*(pres_range[sp]>0)*suit
					range_loss = pres_range[sp]
					range_suit = suit*l_mask
					mat_loss+=range_loss
					mat_suit+=range_suit
			out=rasterio.open('./maps/rasters/map_loss'+sce+'_'+str(year)+'_'+taxon+'.tif', 'w', **meta)
			out.write(mat_loss.astype(rasterio.float64),1)
			out.close()
			out=rasterio.open('./maps/rasters/map_suit'+sce+'_'+str(year)+'_'+taxon+'.tif', 'w', **meta)
			out.write(mat_suit.astype(rasterio.float64),1)
			out.close()
	#	print (sce,year)






out = open('lat_trends.csv','w')
out.write('scenario,taxon,year,lat,loss_log,loss\n')
ls = 'loss'
for sce in ['245','370','585']:
	for taxon in ['Chromista', 'Plantae']:
		data_2015 = gdal_array.LoadFile('./maps/rasters/map_'+ls+sce+'_2015_'+taxon+'.tif')
		cur_dist = where(data_2015>0)
		for year in range(2015,2101):#ls in ['loss']:#['loss','suit']:
			data_2100 = gdal_array.LoadFile('./maps/rasters/map_'+ls+sce+'_'+str(year)+'_'+taxon+'.tif')
			data_log = 100*log(data_2100/data_2015)
			data = 100*(data_2100-data_2015)/data_2015
			for i in range(len(cur_dist[0])):
				y,x = cur_dist[0][i],cur_dist[1][i]
				lat,lon = cell_to_coord(y,x,0.5,0.5,-180,90)
				out.write(','.join(map(str,[sce,taxon,year,lat,data_log[y][x],data[y][x]]))+'\n')
	#		print (','.join(map(str,[sce,taxon,year])))



out.close()



###rescale to 2 x 2
referencefile = 'ocean2.tif'#Path to reference file
reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
referenceProj = reference.GetProjection()
referenceTrans = reference.GetGeoTransform()
x = reference.RasterXSize
y = reference.RasterYSize

for sce in ['245','370','585']:
	for year in [2015,2100]:
		for taxon in ['Chromista', 'Plantae']:
			for ff in ['./maps/rasters/map_loss'+sce+'_'+str(year)+'_'+taxon+'.tif','./maps/rasters/map_suit'+sce+'_'+str(year)+'_'+taxon+'.tif']:
				input = gdal.Open(ff, gdalconst.GA_ReadOnly)
				inputProj = input.GetProjection()
				inputTrans = input.GetGeoTransform()
				driver= gdal.GetDriverByName('GTiff')
				bandreference = input.GetRasterBand(1)
				output = driver.Create(ff.replace('.tif','_2x2.tif'),x,y,1,bandreference.DataType)
				output.SetGeoTransform(referenceTrans)
				output.SetProjection(referenceProj)
				gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Max)
				del output




for sce in ['245','370','585']:
	for year in [2015,2100]:
		for taxon in ['Chromista', 'Plantae']:
			for ff in ['map_loss'+sce+'_'+str(year)+'_'+taxon+'_2x2.tif','map_suit'+sce+'_'+str(year)+'_'+taxon+'_2x2.tif']:
				data = gdal_array.LoadFile('./maps/rasters/'+ff)
				data = log(1+data)
				data[where(data==0)] = nan
				plt.rcParams['hatch.linewidth'] = 0.001
				m = Basemap(resolution='l', projection='cyl',llcrnrlat=-90, llcrnrlon=-180, urcrnrlat=90, urcrnrlon= 180) #projection='vandg',lon_0=0)#
				m.drawmapboundary(fill_color='black',zorder=-1)
				m.fillcontinents(color='darkgrey', lake_color='black',zorder=0)  # zorder=0 to paint over continents
				m.drawcoastlines(linewidth=0.1)
				data_m = ma.masked_where(isnan(data),data)
				vmax = data_m.max()
				vmin = data_m.min()
				x = linspace(-180, 180, data.shape[1]+1)
				y = linspace(90, -90, data.shape[0]+1)
				xx, yy = meshgrid(x, y)
				xx,yy = m(xx,yy)
				im = m.pcolormesh(xx, yy, data_m, cmap=magma,vmin=vmin,vmax=vmax)
				plt.colorbar(orientation = 'horizontal',pad=0.01)
				plt.savefig('./maps/pngs/'+ff[:-4]+'.png', dpi=300,bbox_inches='tight')
				clear = [plt.clf() for i in range(10000)]





ls = 'loss'
for sce in ['245','370','585']:
	for taxon in ['Chromista', 'Plantae']:
		data_2015 = gdal_array.LoadFile('./maps/rasters/map_'+ls+sce+'_2015_'+taxon+'_2x2.tif')
		cur_dist = where(data_2015>0)
		data_2100 = gdal_array.LoadFile('./maps/rasters/map_'+ls+sce+'_2100_'+taxon+'_2x2.tif')
		data = 100*log(data_2100/data_2015)
		plt.rcParams['hatch.linewidth'] = 0.001
		m = Basemap(resolution='l', projection='cyl',llcrnrlat=-90, llcrnrlon=-180, urcrnrlat=90, urcrnrlon= 180) #projection='vandg',lon_0=0)#
		m.drawmapboundary(fill_color='black',zorder=-1)
		m.fillcontinents(color='darkgrey', lake_color='black',zorder=0)  # zorder=0 to paint over continents
		m.drawcoastlines(linewidth=0.1)
		data_m = ma.masked_where(isnan(data),data)
		vmax = data_m.max()
		vmin = data_m.min()
		abs_max = max(abs(vmin),abs(vmax))
		if vmin<0:
			vmin = -abs_max
		vmax = abs_max
		cm = plt.get_cmap('seismic_r')
		x = linspace(-180, 180, data.shape[1]+1)
		y = linspace(90, -90, data.shape[0]+1)
		xx, yy = meshgrid(x, y)
		xx,yy = m(xx,yy)
		im = m.pcolormesh(xx, yy, data_m, cmap=cm,vmin=vmin,vmax=vmax)
		plt.colorbar(orientation = 'horizontal',pad=0.01)
		plt.savefig('./maps/pngs/map_'+ls+'_'+sce+'_hotspots_'+taxon+'.png', dpi=300,bbox_inches='tight')
		clear = [plt.clf() for i in range(10000)]


