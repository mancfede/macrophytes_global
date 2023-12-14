#module load geoconda

from numpy import zeros,array,where,float,argsort,arange

from os import listdir
from osgeo import gdal_array,gdal,gdalconst

from rasterio import features
#from rasterio.features import shapes
from scipy.spatial import Delaunay
from shapely.ops import cascaded_union, polygonize
from random import random,sample,shuffle
from sklearn.ensemble import RandomForestClassifier,RandomForestRegressor

import csv
import os
import rasterio
import math
import shapely.geometry as geometry
import collections




####FUNCTIONS
def alpha_shape(points, alpha):
	if len(points) < 4:
		return geometry.MultiPoint(list(points)).convex_hull
	def add_edge(edges, edge_points, coords, i, j):
		if (i, j) in edges or (j, i) in edges:
			return edges.add( (i, j) )
		edge_points.append(coords[ [i, j] ])
	coords =array([point.coords[0] for point in points])
	tri = Delaunay(coords)
	edges = set()
	edge_points = []
	for ia, ib, ic in tri.vertices:
		pa = coords[ia]
		pb = coords[ib]
		pc = coords[ic]
		a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
		b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
		c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
		s = (a + b + c)/2.0
		area = math.sqrt(s*(s-a)*(s-b)*(s-c))
		circum_r = a*b*c/(4.0*area)
		if circum_r < 1.0/alpha:
			add_edge(edges, edge_points, coords, ia, ib)
			add_edge(edges, edge_points, coords, ib, ic)
			add_edge(edges, edge_points, coords, ic, ia)
	m = geometry.MultiLineString(edge_points)
	triangles = list(polygonize(m))
	return cascaded_union(triangles)


def format_coord(x):
	sw = ['\xc2\xb0S','\xc2\xb0W']
	if any(i in x for i in sw):
		return -(float(x.split('\xc2')[0]))
	else:
		return float(x.split('\xc2')[0])



def point_to_cell(point_x, point_y, cellx, celly, xmin, ymax):
    col = int((point_x - xmin) / cellx)
    row = int((point_y - ymax) / -celly)
    return row,col


def cell_to_coord(col, row, cellx, celly, xmin, ymax):
	lon = cellx*col + xmin
	lat = ymax-celly*row
	return lat,lon


def deg_to_km(lat,res,R=6371.0072):
	return (math.sin(math.radians(lat+(res/2.0)))-math.sin(math.radians(lat-(res/2.0))))*math.radians(res)*(R**2)



###generate 05 grid with area of cells
res = 0.5
a_mat = zeros([360,720])
for row in range(360):
	for col in range(720):
		lat,lon = cell_to_coord(col,row,0.5,0.5,-180,90)
		lat-=res
		a_mat[row][col] = deg_to_km(lat,res)




rst_fn = 'ocean_05.tif'
rst = rasterio.open(rst_fn)
meta = rst.meta.copy()
meta['dtype']='float64'

out=rasterio.open('land_area_05.tif', 'w', **meta)
out.write(a_mat.astype(rasterio.float64),1)
out.close()


####harmonize GEBCO bathymetry
##generate 0.05x0.5 maps x hulls
#referencefile = 'ocean_005.tif'#Path to reference file
#reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
#referenceProj = reference.GetProjection()
#referenceTrans = reference.GetGeoTransform()
#x = reference.RasterXSize
#y = reference.RasterYSize
#
#
#inputfile = './GEBCO/gebco.tif'
#input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
#inputProj = input.GetProjection()
#inputTrans = input.GetGeoTransform()
#outputfile = 'bathymetry_005.tif'#Path to output file
#driver= gdal.GetDriverByName('GTiff')
#bandreference = input.GetRasterBand(1)
#output = driver.Create(outputfile,x,y,1,bandreference.DataType)
#output.SetGeoTransform(referenceTrans)
#output.SetProjection(referenceProj)
#gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Max)
#del output
#
#
##generate 1x1 maps
#referencefile = 'ocean.tif'#Path to reference file
#reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
#referenceProj = reference.GetProjection()
#referenceTrans = reference.GetGeoTransform()
#x = reference.RasterXSize
#y = reference.RasterYSize
#
#
#inputfile = './GEBCO/gebco.tif'
#input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
#inputProj = input.GetProjection()
#inputTrans = input.GetGeoTransform()
#outputfile = 'bathymetry.tif'#Path to output file
#driver= gdal.GetDriverByName('GTiff')
#bandreference = input.GetRasterBand(1)
#output = driver.Create(outputfile,x,y,1,bandreference.DataType)
#output.SetGeoTransform(referenceTrans)
#output.SetProjection(referenceProj)
#gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Max)
#del output
#
##generate 0.5x0.5 maps
#referencefile = 'ocean_05.tif'#Path to reference file
#reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
#referenceProj = reference.GetProjection()
#referenceTrans = reference.GetGeoTransform()
#x = reference.RasterXSize
#y = reference.RasterYSize
#
#
#inputfile = './GEBCO/gebco.tif'
#input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
#inputProj = input.GetProjection()
#inputTrans = input.GetGeoTransform()
#outputfile = 'bathymetry_05.tif'#Path to output file
#driver= gdal.GetDriverByName('GTiff')
#bandreference = input.GetRasterBand(1)
#output = driver.Create(outputfile,x,y,1,bandreference.DataType)
#output.SetGeoTransform(referenceTrans)
#output.SetProjection(referenceProj)
#gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Max)
#del output


#####
occs_file = csv.reader(open('./Macrophytes/macrophyte_occs.csv','r'),delimiter=',')
head = list(enumerate(next(occs_file)))

occs = []
errs = 0
for i in occs_file:
	try:
		lat,lon = map(float,i[1:])
		if -90<=lat<=90 and -180<=lon<=180:
			occs.append([i[0],float(i[1]),float(i[2])])
	except:
		errs+=1
		pass




####check number of occurrences per species, to choose a threshold for the min number of occurrences to include a species
from collections import Counter
c = Counter([i[0] for i in occs])
spp = c.keys()
print (len([i for i in spp if c[i]>0]))#	= 380
print(len([i for i in spp if c[i]>4]))#	= 276
print(len([i for i in spp if c[i]>9]))# = 239
print(len([i for i in spp if c[i]>49]))# = 155


##############RANDOM FOREST FOR EACH SPECIES
###This is done at the best resolution for which we have climatic data, i.e. 0.5 degree
#get reference climate
ocean_mat = gdal_array.LoadFile('ocean_05.tif')
ocean_cells = where(ocean_mat>0)
R,C = ocean_mat.shape
cellx=0.5
celly=0.5
xmin =-180.0
ymax = 90.0


rst_fn = 'ocean_05.tif'
rst = rasterio.open(rst_fn)
meta = rst.meta.copy()
meta.update(compress='lzw')
meta['dtype'] = 'float64'
depth05 = gdal_array.LoadFile('bathymetry_05.tif')
solar_05 = [gdal_array.LoadFile('./solar_05/SIS_'+str(i)+'.tif') for i in range(1,13)]

clim_dir = './climate_05_deg_buffered/'
fff = listdir(clim_dir)
scenarios = sorted(list(set([i.split('_')[0] for i in fff])))
clim_vars = sorted(list(set([i.split('_')[1] for i in fff])))

occs_tre = 10
spp_c = collections.Counter([i[0] for i in occs])
spp_k = spp_c.keys()
spp = sorted([sp for sp in spp_k if spp_c[sp]>=occs_tre])
len(spp)
###check variable importance for the random forest
var_imp = open('RF_var_imp.csv','w')
var_imp.write(','.join(['scenario','species','depth']+[clim_var+str(month) for clim_var in ['solar']+clim_vars for month in range(12)])+'\n')
for sce in scenarios:
	ref_clim = [[gdal_array.LoadFile(clim_dir+'_'.join([sce,clim_vars[i],str(month)])+'.tif') for month in range(60)] for i in range(len(clim_vars))]
	monthly_clim = [[] for i in clim_vars]
	for clim_var in range(len(clim_vars)):
		for month in range(12):
			mat = zeros([R,C])
			for year in range(5):
				mat+=ref_clim[clim_var][year*12+month]
			monthly_clim[clim_var].append(mat/5.0)
	occs_and_clim = []
	for sp,lat,lon in occs:
		row,col=point_to_cell(lon, lat, cellx, celly, xmin, ymax)
		if row >= R:
			row = R
		if col >= C:
			col = C
		loc_clim = [monthly_clim[clim_var][month][row][col] for clim_var in range(len(clim_vars)) for month in range(12)]
		loc_depth = depth05[row][col]
		loc_SIS = [solar_05[i][row][col] for i in range(12)]
		if loc_depth>0:
			loc_depth = 0
		occs_and_clim.append([sp,lat,lon]+[loc_depth]+loc_SIS+loc_clim)
	rf_models = []
	for sp in spp:
		occs_n = spp_c[sp]
		rf = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
	               max_depth=100, max_features='sqrt', max_leaf_nodes=None,
	               min_impurity_decrease=0.0,
	               min_samples_leaf=1, min_samples_split=2,
	               min_weight_fraction_leaf=0.0, n_estimators=1000, n_jobs=4,
	               oob_score=True, random_state=0, verbose=0, warm_start=False)
		pres,psab = [],[]
		pres_coords = set([])
		for row in occs_and_clim:
			if row[0]==sp:
				pres.append(row[3:])
				pres_coords.add((round(row[1],0),round(row[2],0)))
			else:
				if (round(row[1],0),round(row[2],0)) not in pres_coords: #check that pseudo absences are defined without including points where the target species actually occurrs
					psab.append(row[3:])
		while len(psab)<len(pres):
			psab+=psab
		psab = sample(psab,len(pres))
		train = [[1]+i for i in pres]+[[0]+i for i in psab]
		shuffle(train)
		train_y = [i[0] for i in train]
		train_x = [i[1:] for i in train]
		rf_sp = rf.fit(train_x,train_y)
		vi = rf_sp.feature_importances_
		var_imp.write(','.join(map(str,[sce,sp]+list(vi)))+'\n')
		print (sce,sp)



var_imp.close()

###now the actual models
depth = gdal_array.LoadFile('./GEBCO/gebco.tif')
R_d,C_d = depth.shape
cellx_d=1/240
celly_d=1/240


if not os.path.exists('macrophyte_suit_05'):
	os.makedirs('macrophyte_suit_05')


if not os.path.exists('macrophyte_future_suit_05'):
	os.makedirs('macrophyte_future_suit_05')

out = open('random_forest_results.csv','w')
out.write('scenarios,species,oob,occs\n')
out.close()
for sce in scenarios:
	ref_clim = [[gdal_array.LoadFile(clim_dir+'_'.join([sce,clim_vars[i],str(month)])+'.tif') for month in range(60)] for i in range(len(clim_vars))]
	monthly_clim = [[] for i in clim_vars]
	for clim_var in range(len(clim_vars)):
		for month in range(12):
			mat = zeros([R,C])
			for year in range(5):
				mat+=ref_clim[clim_var][year*12+month]
			monthly_clim[clim_var].append(mat/5.0)
	occs_and_clim = []
	for sp,lat,lon in occs:
		row,col=point_to_cell(lon, lat, cellx, celly, xmin, ymax)
		if row >= R:
			row = R
		if col >= C:
			col = C
		row_d,col_d=point_to_cell(lon, lat, cellx_d, celly_d, xmin, ymax)
		if row_d >= R_d:
			row_d = R_d
		if col_d >= C_d:
			col_d = C_d
		loc_clim = [monthly_clim[clim_var][month][row][col] for clim_var in range(len(clim_vars)) for month in range(12)]
		loc_depth = depth[row_d][col_d]
		loc_SIS = [solar_05[i][row][col] for i in range(12)]
		if loc_depth>0:
			loc_depth = 0
		occs_and_clim.append([sp,lat,lon]+[loc_depth]+loc_SIS+loc_clim)
	out = open('random_forest_results.csv','a')
	rf_models = []
	for sp in spp:
		occs_n = spp_c[sp]
		pres,psab = [],[]
		pres_coords = set([])
		for row in occs_and_clim:
			if row[0]==sp:
				pres.append(row[3:])
				pres_coords.add((round(row[1],0),round(row[2],0)))
		for row in occs_and_clim:
			if row[0]!=sp:
				if (round(row[1],0),round(row[2],0)) not in pres_coords: #check that pseudo absences are defined without including points where the target species actually occurrs
					psab.append(row[3:])
		while len(psab)<len(pres):
			psab+=psab
		psab = sample(psab,len(pres))
		train = [[1]+i for i in pres]+[[0]+i for i in psab]
		shuffle(train)
		train_y = [i[0] for i in train]
		train_x = [i[1:] for i in train]
		rf_sp = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
	               max_depth=100, max_features='sqrt', max_leaf_nodes=None,
	               min_impurity_decrease=0.0,
	               min_samples_leaf=1, min_samples_split=2,
	               min_weight_fraction_leaf=0.0, n_estimators=1000, n_jobs=4,
	               oob_score=True, random_state=0, verbose=0, warm_start=False).fit(train_x,train_y)
		vi = rf_sp.feature_importances_
		oob = rf_sp.oob_score_
		sel_var = list(range(len(vi)))
		for var_red in arange(60,0,-10):
			sel_var_ = argsort(vi)[::-1][:var_red]
			sel_train_x = [[i[k] for k in sel_var_] for i in train_x]
			sel_rf_sp = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
	               max_depth=100, max_features='sqrt', max_leaf_nodes=None,
	               min_impurity_decrease=0.0,
	               min_samples_leaf=1, min_samples_split=2,
	               min_weight_fraction_leaf=0.0, n_estimators=1000, n_jobs=4,
	               oob_score=True, random_state=0, verbose=0, warm_start=False).fit(sel_train_x,train_y)
			oob_ = sel_rf_sp.oob_score_
			if oob_>=oob:
				oob,sel_var,rf_sp = oob_,sel_var_,sel_rf_sp
			vi = sel_rf_sp.feature_importances_
		rf_models.append([sp,rf_sp,oob,occs_n,sel_var])
		print (sp,oob,occs_n)
		out.write(','.join(map(str,[sce,sp,oob,occs_n]))+'\n')
	out.close()
	###Make maps of suitability at 05x05
	if not os.path.exists('macrophyte_suit_05/'+sce):
		os.makedirs('macrophyte_suit_05/'+sce)
	if not os.path.exists('macrophyte_future_suit_05/'+sce):
		os.makedirs('macrophyte_future_suit_05/'+sce)
	#present suitability
	glob_clim = []
	for i in range(len(ocean_cells[0])):
		row,col = ocean_cells[0][i],ocean_cells[1][i]
		glob_clim.append([depth05[row][col]]+[solar_05[i][row][col] for i in range(12)]+[monthly_clim[clim_var][month][row][col] for clim_var in range(len(clim_vars)) for month in range(12)])
	for mod in rf_models:
		suit = [i[1] for i in mod[1].predict_proba([[i[j] for j in mod[-1]] for i in glob_clim])]
		sp_mat = zeros([R,C])
		sp_mat[ocean_cells] = suit
		out=rasterio.open('./macrophyte_suit_05/'+sce+'/'+mod[0].replace(' ','_')+'.tif', 'w', **meta)
		out.write(sp_mat.astype(rasterio.float64),1)
		out.close()
		print (sce,mod[0],sp_mat.sum())
	#for future
	for year in range(86): #2015 to 2100; 2100 included
		months = [i+year*12 for i in range(12)]
		monthly_clim = [[gdal_array.LoadFile(clim_dir+'_'.join([sce,clim_vars[i],str(month)])+'.tif') for month in months] for i in range(len(clim_vars))]
		glob_clim = []
		for i in range(len(ocean_cells[0])):
			row,col = ocean_cells[0][i],ocean_cells[1][i]
			glob_clim.append([depth05[row][col]]+[solar_05[i][row][col] for i in range(12)]+[monthly_clim[clim_var][month][row][col] for clim_var in range(len(clim_vars)) for month in range(12)])
		for mod in rf_models:
			suit = [i[1] for i in mod[1].predict_proba([[i[j] for j in mod[-1]] for i in glob_clim])]
			sp_mat = zeros([R,C])
			sp_mat[ocean_cells] = suit
			out_rast=rasterio.open('./macrophyte_future_suit_05/'+sce+'/'+mod[0].replace(' ','_')+'_'+str(2015+year)+'.tif', 'w', **meta)
			out_rast.write(sp_mat.astype(rasterio.float64),1)
			out_rast.close()
			print (sce,year,mod[0],sp_mat.sum())


###Alpha hull to obtain geometric ranges
#####RASTERIZE OCCS AT 0.05X0.05 RESOLUTION
cellx=0.05
celly=0.05
xmin =-180.0
ymax = 90.0

ocean_mat = gdal_array.LoadFile('ocean_005.tif')
R,C = ocean_mat.shape

oceans = where(ocean_mat>0)
ocean_list = list(zip(oceans[0],oceans[1]))


mmm = [[collections.Counter() for i in range(C)] for j in range(R)]

sc=0
for sp,lat,lon in occs:
	row,col=point_to_cell(lon, lat, cellx, celly, xmin, ymax)
	if row >= R:
		row = R
	if col >= C:
		col = C
	mmm[row][col][sp]+=1
	sc+=1
	if sc%10000 == 0:
		print (sc)



agg_occs = []
for row,col in ocean_list:
	if len(mmm[row][col])>0:
		lon,lat = (col*cellx)+xmin,row*(-celly)+ymax
		for i in mmm[row][col].items():
			agg_occs.append([lon,lat,i[1],i[0]])
	print (len(agg_occs))



spp = sorted(list(set([i[-1] for i in agg_occs])))
spp_dict = dict([i[::-1] for i in enumerate(spp)])
sp_occ = [[] for i in range(len(spp))]
sc = 0
for i in agg_occs:
	sc+=1
	sp_occ[spp_dict[i[-1]]].append(i[:3])
	print (len(agg_occs)-sc)



###make output folders
if not os.path.exists('macrophyte_ranges'):
	os.makedirs('macrophyte_ranges')


if not os.path.exists('macrophyte_ranges_ocean'):
	os.makedirs('macrophyte_ranges_ocean')


if not os.path.exists('macrophyte_ranges_05'):
	os.makedirs('macrophyte_ranges_05')


###set meta again for presence absence rasters
meta['dtype'] = 'float64'
#COMPUTE ALPHA HULLS AND GRID DATA AT 1X1 RESOLUTION (IN RASTER FILES)
out_rep = open('alpha_hull_report.csv','w')
out_rep.write('final_area,alpha,area_reduction,occurrence_reduction\n')
for sp in range(len(spp)):
	sp_name = spp[sp].replace(' ','_')
	p_tot = float(sum([ooo[2] for ooo in sp_occ[sp]]))
	if len(sp_occ[sp])>=occs_tre: ###threshold for drawing a hull; the higher the more accurate; min >2
		geom_ok = 'no'
		except_count = 0
		while geom_ok == 'no':
			try:
				ppp = [geometry.Point([ooo[0]+random()/1000.0,ooo[1]+random()/1000.0]) for ooo in sp_occ[sp]]
				alpha = 0.001
				ashape = alpha_shape(ppp,alpha)
				while ashape.area == 0:
					ppp = [geometry.Point([ooo[0]+random()/1000.0,ooo[1]+random()/1000.0]) for ooo in sp_occ[sp]]
					alpha/= 1.1
					ashape = alpha_shape(ppp,alpha)
					if ashape.area == 0:
						ashape = ashape.buffer(0.1)
				geom_ok = 'yes'
			except KeyboardInterrupt:
				raise
			except:
				pass
		tradeoff = 10
		aaas = [ashape]
		p_rem = [sum([sp_occ[sp][i][-1] for i in range(len(ppp)) if ppp[i].intersects(aaas[-1])])]
		att = 0
		while tradeoff>=1.5 and att<1000:
			alpha+=0.001
			ashape = alpha_shape(ppp,alpha)
			if ashape.area == 0:
				ashape = ashape.buffer(0.1)
			aaas.append(ashape)
			p_rem.append(sum([sp_occ[sp][i][-1] for i in range(len(ppp)) if ppp[i].intersects(aaas[-1])]))
			pr = (p_rem[-2]-p_rem[-1])/p_tot
			if pr == 0 or aaas[-2].area==0:
				tradeoff = 1.5
			else:
				ar = (aaas[-2].area-aaas[-1].area)/aaas[-2].area #fraction of area reduction
				tradeoff = ar/pr
			att+=1
		aaa = aaas[-2]
		pr = 1-sum([sp_occ[sp][i][-1] for i in range(len(ppp)) if ppp[i].intersects(aaa)])/p_tot
		ar = (aaas[0].area-aaa.area)/aaas[0].area #fraction of area reduction
		out_rep.write(','.join(map(str,[aaa.area,alpha-0.001,ar,pr]))+'\n')
		out_fn='macrophyte_ranges/'+sp_name+'.tif'
		out_fn_ocean='macrophyte_ranges_ocean/'+sp_name+'.tif'
		if aaa.type == 'Polygon':
			geom = [[aaa,1]]
		else:
			geom=[[j,1] for j in aaa.geoms]
	else:
		ppp_multi = [geometry.Point(ooo[:2]) for ooo in sp_occ[sp] if ooo[2]>5]
		geom=[[j.buffer(0.01),1] for j in ppp_multi]
	if geom!=[]:
		out=rasterio.open(out_fn, 'w+', **meta)
		out_arr = out.read(1)
		burned = features.rasterize(shapes=geom, fill=0, out=out_arr, transform=out.transform)
		out.write_band(1, burned)
		out.close()
		mat = gdal_array.LoadFile(out_fn)
		if mat.sum()>0:
			out=rasterio.open(out_fn_ocean, 'w', **meta)
			out.write(mat.astype(rasterio.uint8),1)
			out.close()
		print (sp_name,aaa.area,ar,pr)


out_rep.close()



#generate 0.5x0.5 maps
referencefile = 'ocean_05.tif'#Path to reference file
reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
referenceProj = reference.GetProjection()
referenceTrans = reference.GetGeoTransform()
x = reference.RasterXSize
y = reference.RasterYSize


spp = sorted(list(set([i for i in listdir('macrophyte_ranges_ocean') if i[-4:]=='.tif'])))


for sp in spp:
	inputfile = './macrophyte_ranges_ocean/'+sp
	input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
	inputProj = input.GetProjection()
	inputTrans = input.GetGeoTransform()
	outputfile = './macrophyte_ranges_05/'+sp#Path to output file
	driver= gdal.GetDriverByName('GTiff')
	bandreference = input.GetRasterBand(1)
	output = driver.Create(outputfile,x,y,1,bandreference.DataType)
	output.SetGeoTransform(referenceTrans)
	output.SetProjection(referenceProj)
	gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Average)
	del output
	print (sp)


###habitat
import sys
csv.field_size_limit(sys.maxsize)

ocean_mat = gdal_array.LoadFile('ocean_05.tif')
ocean_cells = where(ocean_mat>0)
R,C = ocean_mat.shape
cellx=0.5
celly=0.5
xmin =-180.0
ymax = 90.0


rst_fn = 'ocean_05.tif'
rst = rasterio.open(rst_fn)
meta = rst.meta.copy()
meta.update(compress='lzw')
meta['dtype'] = 'float64'
depth05 = gdal_array.LoadFile('bathymetry_05.tif')
solar_05 = [gdal_array.LoadFile('./solar_05/SIS_'+str(i)+'.tif') for i in range(1,13)]

clim_dir = './climate_05_deg_buffered/'
fff = listdir(clim_dir)
scenarios = sorted(list(set([i.split('_')[0] for i in fff])))
clim_vars = sorted(list(set([i.split('_')[1] for i in fff])))

occs_tre = 10
spp_c = collections.Counter([i[0] for i in occs])
spp_k = spp_c.keys()
spp = sorted([sp for sp in spp_k if spp_c[sp]>=occs_tre])

assis = csv.reader(open('./Macrophytes/database_pruned_assis_et_al.csv','r'),delimiter=';')
head = list(enumerate(next(assis)))

#species = 101
#kingdom = 102
mac_alg = dict([])
for i in assis:
	try:
		if i[102]!='':
			mac_alg[i[101]] = i[102]
	except:
		pass




spp_used = set([i.replace('.tif','').replace('_',' ') for i in listdir('./macrophyte_ranges_05/')])-set(['Sargassum fluitans','Sargassum natans','Sargassum pusillum'])

mat_chrom = [[set([]) for i in range(int(360/cellx))] for j in range(int(180/celly))]
mat_plant = [[set([]) for i in range(int(360/cellx))] for j in range(int(180/celly))]

for sp,lat,lon in occs:
	if sp in spp_used:
		row,col=point_to_cell(lon, lat, cellx, celly, xmin, ymax)
		if row >= R:
			row = R
		if col >= C:
			col = C
		taxon = mac_alg.get(sp,'none')
		if taxon=='Chromista':
			mat_chrom[row][col].add(sp)
		elif taxon=='Plantae':
			mat_plant[row][col].add(sp)
		else:
			print (sp)


for row in range(int(180/celly)):
	for col in range(int(360/cellx)):
		mat_chrom[row][col] = len(mat_chrom[row][col])
		mat_plant[row][col] = len(mat_plant[row][col])



for row in range(int(180/celly)):
	for col in range(int(360/cellx)):
		if mat_chrom[row][col]>0:
			mat_chrom[row][col] = 1
		if mat_plant[row][col]>0:
			mat_plant[row][col] = 1




l_mask = rasterio.open('./light_at_bottom/Present_Benthic_Max_Depth_Light_bottom_Max_05.tif').read(1)
l_mask = 1*(l_mask>50)

###Make a single prediction for the habitat
out = open('random_forest_results_habitat.csv','w')
out.write('scenarios,taxon,oob\n')
out.close()
for sce in scenarios:
	ref_clim = [[gdal_array.LoadFile(clim_dir+'_'.join([sce,clim_vars[i],str(month)])+'.tif') for month in range(60)] for i in range(len(clim_vars))]
	monthly_clim = [[] for i in clim_vars]
	for clim_var in range(len(clim_vars)):
		for month in range(12):
			mat = zeros([R,C])
			for year in range(5):
				mat+=ref_clim[clim_var][year*12+month]
			monthly_clim[clim_var].append(mat/5.0)
	occs_and_clim = []
	for row in range(int(180/celly)):
		for col in range(int(360/cellx)):
			if l_mask[row][col]==1:
				loc_clim = [monthly_clim[clim_var][month][row][col] for clim_var in range(len(clim_vars)) for month in range(12)]
				if -9999 not in loc_clim:
					loc_depth = depth05[row][col]
					loc_SIS = [solar_05[i][row][col] for i in range(12)]
					if loc_depth>0:
						loc_depth = 0
					occs_and_clim.append(['Chromista',mat_chrom[row][col]]+[loc_depth]+loc_SIS+loc_clim)
					occs_and_clim.append(['Plantae',mat_plant[row][col]]+[loc_depth]+loc_SIS+loc_clim)
	out = open('random_forest_results_habitat.csv','a')
	rf_models = []
	for sp in ['Chromista','Plantae']:
		rf = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
               max_depth=100, max_features='sqrt', max_leaf_nodes=None,
               min_impurity_decrease=0.0, 
               min_samples_leaf=1, min_samples_split=2,
               min_weight_fraction_leaf=0.0, n_estimators=1000, n_jobs=4,
               oob_score=True, random_state=0, verbose=0, warm_start=False)
		pres,psab = [],[]
		for row in occs_and_clim:
			if row[0]==sp:
				if row[1]>0:
					pres.append(row[1:])
				else:
					psab.append(row[1:])
		while len(psab)<len(pres):
			psab+=psab
		psab = sample(psab,len(pres))
		train = [i for i in pres]+[i for i in psab]
		shuffle(train)
		train_y = [i[0] for i in train]
		train_x = [i[1:] for i in train]
		rf_sp = rf.fit(train_x,train_y)
		rf_models.append([sp,rf_sp,rf_sp.oob_score_])
		print (sp,rf_sp.oob_score_)
		out.write(','.join(map(str,[sce,sp,rf_sp.oob_score_]))+'\n')
	out.close()
	###Make maps of suitability at 05x05
	if not os.path.exists('macrophyte_future_suit_05_hab/'+sce):
		os.makedirs('macrophyte_future_suit_05_hab/'+sce)
	for year in range(86): #2015 to 2100; 2100 included
		months = [i+year*12 for i in range(12)]
		monthly_clim = [[gdal_array.LoadFile(clim_dir+'_'.join([sce,clim_vars[i],str(month)])+'.tif') for month in months] for i in range(len(clim_vars))]
		glob_clim = []
		for i in range(len(ocean_cells[0])):
			row,col = ocean_cells[0][i],ocean_cells[1][i]
			glob_clim.append([depth05[row][col]]+[solar_05[i][row][col] for i in range(12)]+[monthly_clim[clim_var][month][row][col] for clim_var in range(len(clim_vars)) for month in range(12)])
		for mod in rf_models:
			suit = [i[1] for i in mod[1].predict_proba(glob_clim)]
			sp_mat = zeros([R,C])
			sp_mat[ocean_cells] = suit
			out_rast=rasterio.open('./macrophyte_future_suit_05_hab/'+sce+'/'+mod[0].replace(' ','_')+'_'+str(2015+year)+'.tif', 'w', **meta)
			out_rast.write(sp_mat.astype(rasterio.float64),1)
			out_rast.close()
			print (sce,year,mod[0],sp_mat.sum())




