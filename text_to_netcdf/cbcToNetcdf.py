import flopy.utils.binaryfile as binaryfile
import netCDF4
import numpy
import osr
import sys

def find_values(cbc, namePosition, layer):

	value = cbc.get_record(namePosition)[layer, :, :]
	return value


def find_average_resolution(fileHandle, numberOfHRUs, numberOfRows, numberOfColumns):

	latitudeValues = []
	longitudeValues = []

	for i in range(numberOfHRUs):
		valuesInLine = fileHandle.next().strip().split()
		longitudeValues.append(float(valuesInLine[1]))
		latitudeValues.append(float(valuesInLine[2]))

	minimumLatitudeValue = min(latitudeValues)
	maximumLatitudeValue = max(latitudeValues)

	minimumLongitudeValue = min(longitudeValues)
	maximumLongitudeValue = max(longitudeValues)

	averageOfLatitudeValues = (maximumLatitudeValue-minimumLatitudeValue)/numberOfRows
	averageOfLongitudeValues = (maximumLongitudeValue-minimumLongitudeValue)/numberOfColumns

	latitudeOfFirstHru = latitudeValues[0]
	longitudeOfFirstHru = longitudeValues[0]

	return averageOfLatitudeValues, averageOfLongitudeValues, latitudeOfFirstHru, longitudeOfFirstHru


def cbc_to_netcdf(cbcFile, disFile, locationFile, fileOutput):

	cbc = binaryfile.CellBudgetFile(cbcFile)
	recordNames = cbc.unique_record_names()
	lengthOfRecordNames = len(recordNames)
	numberOfRecords = cbc.get_nrecords()
	
	fileHandle = open(disFile, 'r')
	for line in fileHandle:
		if "#" not in line.strip().split()[0]:
			values = line.split()
			break

	numberOfLayers = int(values[0])
	numberOfRows = int(values[1])
	numberOfColumns = int(values[2])
	numberOfHRUs = numberOfRows * numberOfColumns
	numberOfStressPeriods = int(values[3])

	timeStepForStressPeriods = []

	for index in range(numberOfStressPeriods, 0, -1):
		fileHandle = open(disFile, 'r')
		lineNumber = index * -1
		values = fileHandle.readlines()[lineNumber].split()
		numberOfTimeSteps = int(values[1])
		timeStepForStressPeriods.append(numberOfTimeSteps)

	fileHandle = open(locationFile, 'r')
	values = find_average_resolution(fileHandle, numberOfHRUs, numberOfRows, numberOfColumns)
	averageOfLatitudeValues = values[0]
	averageOfLongitudeValues = values[1]
	latitudeOfFirstHru = values[2]
	longitudeOfFirstHru = values[3]

	# Initialize new dataset
	ncfile = netCDF4.Dataset(fileOutput, mode='w')

	# Initialize dimensions
	lat_dim = ncfile.createDimension('lat', numberOfRows)
	lon_dim = ncfile.createDimension('lon', numberOfColumns)
	
	latList = []
	latList.append(latitudeOfFirstHru)
	previousValue = latitudeOfFirstHru
	lat = ncfile.createVariable('lat', 'f8', ('lat',))
	lat.long_name = 'latitude'  
	lat.units = 'degrees_north'
	for i in range(numberOfRows - 1):
		newValue = previousValue - averageOfLatitudeValues
		latList.append(newValue)
		previousValue = newValue
	lat[:] = latList

	lonList = []
	lonList.append(longitudeOfFirstHru)
	previousValue = longitudeOfFirstHru
	lon = ncfile.createVariable('lon', 'f8', ('lon',))
	lon.long_name = 'longitude'  
	lon.units = 'degrees_east'
	for i in range(numberOfColumns - 1):
		newValue = previousValue + averageOfLongitudeValues
		lonList.append(newValue)
		previousValue = newValue
	lon[:] = lonList

	sr = osr.SpatialReference()
	sr.ImportFromEPSG(4326)
	crs = ncfile.createVariable('crs', 'S1',)
	crs.spatial_ref = sr.ExportToWkt()

	namePosition = 0

	for i in range(numberOfStressPeriods):
		for j in range(timeStepForStressPeriods[i]):
			for p in recordNames:
				for q in range(numberOfLayers):
					var = ncfile.createVariable(p.strip()+"_"+str(i+1)+"_"+str(j+1)+"_"+str(q+1), 'f8', ('lat', 'lon'))
					var.layer_name = p.strip()+"_"+str(i+1)+"_"+str(j+1)+"_"+str(q+1)
					var.layer_desc = "Record " +p.strip()+" Stress Period "+str(i+1)+ " Time Step " +str(j+1)+ " Layer " +str(q+1)
					#var.layer_units = lengthUnit
					var.grid_mapping = "crs" 
					values = find_values(cbc, namePosition, q)
					var[:] = values
				namePosition += 1
				

if __name__ == "__main__":

	numberOfArgs = len(sys.argv)

	for index in range(numberOfArgs):

		if sys.argv[index] == "-cbc":
			cbcFile = sys.argv[index+1]

		elif sys.argv[index] == "-dis":
			disFile = sys.argv[index+1]

		elif sys.argv[index] == "-loc":
			locationFile = sys.argv[index+1]

	cbc_to_netcdf(cbcFile, disFile, locationFile, 'cbc.nc')