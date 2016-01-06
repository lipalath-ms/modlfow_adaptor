import netCDF4
import osr
import sys

def find_time_unit(timeValue):
	
	if timeValue == 0:
		timeUnit = "undefined"
	elif timeValue == 1:
		timeUnit = "seconds"
	elif timeValue == 2:
		timeUnit = "minutes"
	elif timeValue == 3:
		timeUnit = "hours"
	elif timeValue == 4:
		timeUnit = "days"
	elif timeValue == 5:
		timeUnit = "years"
	return timeUnit			

def find_length_unit(lengthValue):
	
	if lengthValue == 0:
		lengthUnit = "undefined"
	elif lengthValue == 1:
		lengthUnit = "feet"
	elif lengthValue == 2:
		lengthUnit = "meters"
	elif lengthValue == 3:
		lengthUnit = "centimeters"
	return lengthUnit			

def find_values(fileHandle, numberOfHRUs, stressPeriodIndex, layerIndex):

	hruValues = [] 
		
	for line in fileHandle:
		sum = 0
		if "HEAD" in line:
			indexValues = line.strip().split()
			if int(indexValues[1]) == stressPeriodIndex and int(indexValues[7]) == layerIndex:
				while sum != numberOfHRUs:
					values = fileHandle.next().split()
					lengthOfLine = len(values)
					sum += lengthOfLine
					for index in range(lengthOfLine):
						hruValues.append(float(values[index]))

	return hruValues
	

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
	

def fhd_to_netcdf(fhdFile, disFile, locationFile, fileOutput):

	totalNumberOfTimeSteps = 0

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
	timeValue = int(values[4])
	timeUnit = find_time_unit(timeValue)
	lengthValue = int(values[5])
	lengthUnit = find_length_unit(lengthValue)
	
	stressPeriodLengths = []
	timeStepCounts = []
	multiplierCounts = []
	stressPeriodStates = []
	
	for index in range(numberOfStressPeriods, 0, -1):
		fileHandle = open(disFile, 'r')
		lineNumber = index * -1
		values = fileHandle.readlines()[lineNumber].split()
		lengthOfStressPeriod = float(values[0])
		numberOfTimeSteps = int(values[1])
		multiplier = float(values[2])
		stateOfStressPeriod = values[3]

		stressPeriodLengths.append(lengthOfStressPeriod)
		timeStepCounts.append(numberOfTimeSteps)
		multiplierCounts.append(multiplier)
		stressPeriodStates.append(stateOfStressPeriod)

	
	fileHandle = open(locationFile, 'r')
	values = find_average_resolution(fileHandle, numberOfHRUs, numberOfRows, numberOfColumns)
	averageOfLatitudeValues = values[0]
	averageOfLongitudeValues = values[1]
	latitudeOfFirstHru = values[2]
	longitudeOfFirstHru = values[3]

	# Initialize new dataset
	ncfile = netCDF4.Dataset('fhd.nc', mode='w')

	# Initialize dimensions

	for index in range(numberOfStressPeriods):
		timestep = ncfile.createDimension('timestep'+"_"+str(index+1), timeStepCounts[index])

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
	
	
 	for i in range(numberOfStressPeriods):
 		for j in range(numberOfLayers):
 			var = ncfile.createVariable("head_"+str(i+1)+"_"+str(j+1), 'f8', ('timestep'+"_"+str(i+1), 'lat', 'lon'))
			var.description = "Water head for HRUs in Stress Period " +str(i+1)+ " Layer " +str(j+1)
			var.units = lengthUnit
			fileHandle = open(fhdFile, 'r')
			values = find_values(fileHandle, numberOfHRUs, i+1, j+1)
			var[:] = values

	# Global attributes
    ncfile.title = 'Modflow Formatted Head Package'
    ncfile.bands = 1
    ncfile.bands_name = 'nsteps'
    ncfile.bands_desc = 'Variable information for ' + fhdFile
    ncfile.number_of_stress_periods = numberOfStressPeriods
    ncfile.number_of_layers = numberOfLayers
    ncfile.time_unit = timeUnit
    ncfile.length_unit = lengthUnit
    
    # Close the 'ncfile' object
    ncfile.close()

	
if __name__ == "__main__":

	numberOfArgs = len(sys.argv)

	for index in range(numberOfArgs):

		if sys.argv[index] == "-fhd":
			fhdFile = sys.argv[index+1]
		elif sys.argv[index] == "-dis":
			disFile = sys.argv[index+1]
		elif sys.argv[index] == "-loc":
			locationFile = sys.argv[index+1]

	fhd_to_netcdf(fhdFile, disFile, locationFile, 'fhd.nc')
    
    
