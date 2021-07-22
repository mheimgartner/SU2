# Script which computes and outputs area and mass flow averaged variables from Paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`
import numpy as np 
import math

#### Import the simple module from the paraview
from paraview.simple import *
#### Disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'XMLUnstructuredGrid Reader'
flowvtu = XMLUnstructuredGridReader(FileName='/home/hem1dev/Codes/SU2_Github_Mark_Test/SU2/Venturi/flow_verification_davidson_00000749.vtu') 

# Get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1047, 659]

# Show data in view
flowvtuDisplay = Show(flowvtu, renderView1)

# Trace defaults for the display properties.
flowvtuDisplay.Representation = 'Surface'

# Reset view to fit data
renderView1.ResetCamera()

# Update the view to ensure updated data information
renderView1.Update()

# Get layout
layout1 = GetLayout()

# Create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(Input=flowvtu,Source='High Resolution Line Source')

# Init the 'High Resolution Line Source' selected for 'Source'
plotOverLine1.Source.Point1 = [1, 0, 0.0]
plotOverLine1.Source.Point2 = [1, 0.1, 0.0]

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Resolution = 1000

# Update the view to ensure updated data information
renderView1.Update() 

# Get layout
layout1 = GetLayout()

# Retrieve data from Paraview 
plotOverLine1.UpdatePipeline()
data = servermanager.Fetch(plotOverLine1)
pointdata = data.GetPointData()
numCells = pointdata.GetNumberOfTuples()
#print dir(pointdata)

# Allocate arrays
YCH4 = []
MeanMolecularWeight = []
Density = []
LaminarViscosity = []
MassDiffusivity = []
SpecificHeatCp = []
ThermalConductivity = []

#DON'T FORGET TO POSSIBLY SORT THE DATA BECAUSE WHEN YOU LOAD IN THE 30KW_NODE.CASE FILE IN PARAVIEW THE ORDER IS NOT LEXICOGRAPHIC... 
#ALSO CHECK IF YOU CAN DIRECTLY LOAD THE OUTLET DATA INSTEAD OF USING PLOTOVERLINE 

# Loop to retrieve variables from Paraview
for x in range(numCells):
	YCH4.append(pointdata.GetArray("Passive_Scalar").GetValue(x))
	MeanMolecularWeight.append(pointdata.GetArray("Mean_Molecular_Weight").GetValue(x))
	Density.append(pointdata.GetArray("Density").GetValue(x))
	LaminarViscosity.append(pointdata.GetArray("Laminar_Viscosity").GetValue(x))
	MassDiffusivity.append(pointdata.GetArray("Diffusivity").GetValue(x))
	SpecificHeatCp.append(pointdata.GetArray("Specific_Heat_Cp").GetValue(x))
	ThermalConductivity.append(pointdata.GetArray("Conductivity").GetValue(x))
	
# Write out data to .csv file to allow postprocessing in jupyter
data_export = np.column_stack((YCH4,MeanMolecularWeight,Density,LaminarViscosity,MassDiffusivity,SpecificHeatCp,ThermalConductivity))
header = "YCH4,MeanMolecularWeight,Density,LaminarViscosity,MassDiffusivity,SpecificHeatCp,ThermalConductivity"
np.savetxt("smallboxverification_davidson.csv", data_export, header = header, delimiter=",", comments='')
