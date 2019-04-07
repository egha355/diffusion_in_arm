#> \author Elias Ghadam Soltani
#> \brief This is an example program to solve a diffusion equation using OpenCMISS calls.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s):
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. If you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. If you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>





#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,sys,os,time,math,cmath,csv
from shutil import copyfile
from opencmiss.iron import iron
t=time.time()
# Diagnostics
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",[""])
#iron.ErrorHandlingModeSet(iron.ErrorHandlingModes.TRAP_ERROR)
#iron.OutputSetOn("Testing")

# Set the OpenCMISS random seed so that we can test this example by using the
# same parallel decomposition.
numberOfRandomSeeds = iron.RandomSeedsSizeGet()
randomSeeds = [0]*numberOfRandomSeeds
randomSeeds[0] = 100
iron.RandomSeedsSet(randomSeeds)

# Get the computational nodes info
#computationEnvironment = iron.ComputationEnvironment()
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber    = iron.ComputationalNodeNumberGet()

print(' ')

#================================================================================================================================
# Modules
#================================================================================================================================
# Copying current example file into output directory

if computationalNodeNumber==0:
  if os.path.exists("./output"):
    input('\033[1;31m'+'Output directory already exists. RENAME IT first, if required.'+'\033[0m')
  else:
    os.makedirs("./output/Details")
  copyfile(os.path.basename(__file__), './output/Details/file.py')

def module_exists(module_name):
    try:
        __import__(module_name)
    except ImportError:
        return False
    else:
        return True

if computationalNodeNumber==0 and module_exists('DescriptionGenerator'):
  import DescriptionGenerator

#================================================================================================================================
#  Problem Control Panel
#================================================================================================================================
numberOfDimensions     = 1  #(One-dimensional)
#
# derivIdx = 1
#
ProgressDiagnostics = False   # Set to diagnostics
#
# energySolve       = False # Set to solve energy equation
# tissueSolve       = False # Set to solve energy equation
# tissueEnergySolve = True  # Set to solve energy equation
meshOrigin = [0.0,0.0,0.0]

#================================================================================================================================
#  Start Program
#================================================================================================================================

controlLoopNode                        = 0
coordinateSystemUserNumber             = 1
regionUserNumberEnergy                 = 2
regionUserNumberTissue                 = 3
basisUserNumberEnergy                  = 4
BasisUserNumberTissue                  = 5
meshUserNumberEnergy                   = 6
meshUserNumberTissue                   = 7
decompositionUserNumberEnergy          = 8
decompositionUserNumberTissue          = 9
geometricFieldUserNumberEnergy         = 10
geometricFieldUserNumberTissue         = 11
equationsSetUserNumberEnergy           = 12
equationsSetUserNumberTissue           = 13
equationsSetFieldUserNumberEnergy      = 14
equationsSetFieldUserNumberTissue      = 15
dependentFieldUserNumberEnergy         = 16
dependentFieldUserNumberTissue         = 17
materialsFieldUserNumberEnergy         = 18
materialsFieldUserNumberTissue         = 19
sourceFieldUserNumberEnergy            = 20
sourceFieldUserNumberTissue            = 21
independentFieldUserNumber             = 22
problemUserNumber                      = 23

#================================================================================================================================
#  Mesh Reading
#================================================================================================================================

# LENGTH = 0.5 # m
# numberGlobalXElements = 192
numberOfNodesEnergy    = 26
numberOfElementsEnergy = 25

#ELEMENTS
elementNodes= (numberOfElementsEnergy)*[2*[0]]
#Brachial artery
elementNodes[0]=[1,2]
elementNodes[1]=[2,3]
elementNodes[2]=[3,4]
elementNodes[3]=[4,5]
elementNodes[4]=[5,6]

#Radial
elementNodes[5] =[6,7]
elementNodes[6] =[7,8]
elementNodes[7] =[8,9]
elementNodes[8] =[9,10]
elementNodes[9]=[10,11]
elementNodes[10]=[11,12]
elementNodes[11]=[12,13]
elementNodes[12]=[13,14]
elementNodes[13]=[14,15]
elementNodes[14]=[15,16]

#Ulnar
elementNodes[15]=[6,17]
elementNodes[16]=[17,18]
elementNodes[17]=[18,19]
elementNodes[18]=[19,20]
elementNodes[19]=[20,21]
elementNodes[20]=[21,22]
elementNodes[21]=[22,23]
elementNodes[22]=[23,24]
elementNodes[23]=[24,25]
elementNodes[24]=[25,26]

#NODES
xValues = numpy.zeros((numberOfNodesEnergy,1),dtype = numpy.float)
yValues = numpy.zeros((numberOfNodesEnergy,1),dtype = numpy.float)
zValues = numpy.zeros((numberOfNodesEnergy,1),dtype = numpy.float)
radius= (numberOfNodesEnergy)*[1*[0]]
ra=(numberOfElementsEnergy)*[1*[0]]
Tt= (numberOfNodesEnergy)*[1*[0]]
#brachial artery
xValues[0]=152.12 #mm
xValues[1]=164.04
xValues[2]=179
xValues[3]=187.37
xValues[4]=195.89
xValues[5]=215.5

yValues[0]=-80.904
yValues[1]=-76.816
yValues[2]=-76.099
yValues[3]=-82.713
yValues[4]=-86.659
yValues[5]=-85.043

zValues[0]=1243.4
zValues[1]=1214.5
zValues[2]=1171
zValues[3]=1123
zValues[4]=1073.7
zValues[5]=1022.5

radius[0]=3.3062 #mm
radius[1]=3.1344
radius[2]=2.8703
radius[3]=2.5923
radius[4]=2.2908
radius[5]=1.9026

Tt[0]=37.0 #C
#Tt[1]=36.77
Tt[1]=37.0
# Tt[2]=36.67
Tt[2]=37.0
#Tt[3]=36.57
Tt[3]=37.0
# Tt[4]=36.47
Tt[4]=37.0
# Tt[5]=36.37
Tt[5]=37.0

#Radial
# le   = xValues[1][0]
xValues[6] =227.3
xValues[7] =239.54
xValues[8] =250.22
xValues[9] =256.49
xValues[10]=256.77
xValues[11]=258.42
xValues[12]=262.93
xValues[13]=267.01
xValues[14]=270.44
xValues[15]=278.41

yValues[6] =-90.085
yValues[7] =-93.529
yValues[8] =-93.818
yValues[9] =-90.243
yValues[10]=-85.607
yValues[11]=-88.336
yValues[12]=-100.68
yValues[13]=-115.58
yValues[14]=-121.21
yValues[15]=-119.14

zValues[6] =1008.6
zValues[7] =992.83
zValues[8] =975.58
zValues[9] =954.23
zValues[10]=930.66
zValues[11]=905.41
zValues[12]=877.61
zValues[13]=852.92
zValues[14]=823.34
zValues[15]=792.14

radius[6] =1.8158 #mm
radius[7] =1.7652
radius[8] =1.7176
radius[9] =1.6593
radius[10]=1.6072
radius[11]=1.5359
radius[12]=1.4673
radius[13]=1.3942
radius[14]=1.305
radius[15]=0.93851

#Tt[6] =36.26
# Tt[6]=37.0
# Tt[7] =36.2
# Tt[8] =36.14
# Tt[9] =36.08
# Tt[10]=36.02
# Tt[11]=35.96
# Tt[12]=35.90
# Tt[13]=35.84
# Tt[14]=35.76
# Tt[15]=35.70

#Ulnar
xValues[16]=213.85
xValues[17]=212.98
xValues[18]=213.25
xValues[19]=214.23
xValues[20]=215.15
xValues[21]=221.02
xValues[22]=228.77
xValues[23]=235.44
xValues[24]=239.83
xValues[25]=242.49

yValues[16]=-80.244
yValues[17]=-78.828
yValues[18]=-85.493
yValues[19]=-96.217
yValues[20]=-101.56
yValues[21]=-105.22
yValues[22]=-111.98
yValues[23]=-121.22
yValues[24]=-128.52
yValues[25]=-128.92

zValues[16]=1009.5
zValues[17]=999.14
zValues[18]=982.26
zValues[19]=962.45
zValues[20]=928.98
zValues[21]=893.59
zValues[22]=868.46
zValues[23]=844.1
zValues[24]=820.37
zValues[25]=795.86

radius[16]=1.8299
radius[17]=1.797
radius[18]=1.7584
radius[19]=1.6931
radius[20]=1.6111
radius[21]=1.5121
radius[22]=1.4446
radius[23]=1.3709
radius[24]=1.3085
radius[25]=1.1698

# Tt[16]=36.26
# Tt[17]=36.20
# Tt[18]=36.14
# Tt[19]=36.08
# Tt[20]=36.02
# Tt[21]=35.96
# Tt[22]=35.90
# Tt[23]=35.84
# Tt[24]=35.76
# Tt[25]=35.70

for nodeIdx in range(numberOfNodesEnergy):
    Tt[nodeIdx]=37 #C

arteriesElements=range(1,numberOfElementsEnergy+1)

for elemIdx in range(numberOfElementsEnergy):
    ra[elemIdx]=(radius[elemIdx]+radius[elemIdx+1])/2

ra[15]=(radius[5]+radius[16])/2 #ulnar nodes in bifurcation is different

print(ra)
def elementLength(p1,p2):
    length=math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)
    return length

le= (numberOfElementsEnergy)*[1*[0]]
for i in range(numberOfElementsEnergy):
    p1=[xValues[i],yValues[i],zValues[i]]
    p2=[xValues[i+1],yValues[i+1],zValues[i+1]]
    le[i]=elementLength(p1,p2) #mm

brachialElements=[1,2,3,4,5]
radialElements  =[6,7,8,9,10,11,12,13,14,15]
ulnarElements   =[16,17,18,19,20,21,22,23,24,25]


# Tissue mesh
FileName_ele = "MaxVol1000/WholeLeftArm.1.ele"
numberOfNodesTissue    = 40813
numberOfElementsTissue = 237070
numberOfLocalNodes = 4
offset = 1 # offset is 1 if nodes and elements number begin from 0, offset for default = 0

# Printing the head and tail to see if we are reading the information correctly.
printHead = True
printTail = True
head = range(5)
tail = range(numberOfElementsTissue-4,numberOfElementsTissue+1)

# Label for parts
muscleRegionLabel      = 2
leftRadiusRegionLabel  = 3
leftHumerusRegionLabel = 4
leftUlnaRegionLabel    = 5

muscleElements      = []
leftRadiusElements  = []
leftHumerusElements = []
leftUlnaElements    = []

localNodes=[0]*numberOfLocalNodes

elementNumber=0
totalNumberOfElements=1
# elementNodes = False
#with open("upperlimb_refined2.exelem", "r") as f:
#================================================================================================================================
#  Initial Data & Default values
#================================================================================================================================


# Artery =========
k_bl               = 0.5*1e-3      # W/mm.K blood conductivity.
rho_bl             = 1069.0*1e-9   # kg/mm3 blood density
c_bl               = 3650.0        # J/Kg.K blood specific heat
Alpha              = k_bl/(rho_bl*c_bl)       # mm2/s Diffusivity
# r                  = 1.5           # mm, inner radius of the artery
# CArtery            = 4*Alpha/(r*r) # 0.27911 1/s
# Tt                 = 37.0          # C

# Tissue ==========

rho_ms             = 1085.0*1e-9   # kg/mm3   muscle density
c_ms               = 3768.0        # J/Kg.K   muscle specific heat
rho_bn             = 1357.0*1e-9   # kg/mm3    bone density
c_bn               = 1700.0        # J/Kg.K   bone specific heat
rho_sk             = 1085.0*1e-9   # kg/mm3    skin density
c_sk               = 3680.0        # J/Kg.K   skin specific heat

k_ms               = 0.42*1e-3     # W/mm.K muscle conductivity.
k_bn               = 0.75*1e-3     # W/mm.K bone conductivity.
k_sk               = 0.47*1e-3     # W/mm.K skin conductivity.

h_conv             = 2.0*1e-6      # W/mm2.K
#h_conv            = 200.0*1e-6    # W/mm2.K for water
hr_rad             = 5.9*1e-6      # W/mm2.K See example 3.12 Incropera

# R_arm              = 0.03          # m

Tb                 = 37.0          # C blood temeprature
Tair               = 24.0          # C

w                  = 5e-4          # 1/s terminal blood flow per volume of tissue.

cMuscle            = rho_bl*c_bl/(rho_ms*c_ms) *w   # 4.51128e-4 1/s

cBone              = 0.0           #\see equations section

# Set the time parameters
timeIncrement   = 0.5
startTime       = 0.0
stopTime        = 4200

# Set the output parameters
DYNAMIC_SOLVER_DIFFUSION_OUTPUT_FREQUENCY = 100

# Set the solver parameters
#relativeTolerance   = 1.0E-05  # default: 1.0E-05
#absoluteTolerance   = 1.0E-08  # default: 1.0E-10
#DIVERGENCE_TOLERANCE = 1.0e+10  # default: 1.0e+05
MAXIMUM_ITERATIONS   = 1000   # default: 100000
#RESTART_VALUE        = 3000     # default: 30

# Navier-Stokes solver
equationsSetEnergySubtype = iron.EquationsSetSubtypes.ADVECTION_DIFFUSION
equationsSetTissueSubtype = iron.EquationsSetSubtypes.LINEAR_SOURCE_DIFFUSION
ProblemSubtype      = iron.ProblemSubtypes.THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION

#================================================================================================================================
#  Coordinate System
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> COORDINATE SYSTEM << == ")

# Start the creation of RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

print('\033[1;32m'+'Coordinate System COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Region
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> COORDINATE SYSTEM << == ")

# Start the creation of Energy region
regionEnergy = iron.Region()
regionEnergy.CreateStart(regionUserNumberEnergy,iron.WorldRegion)
regionEnergy.LabelSet("Energy")
regionEnergy.CoordinateSystemSet(coordinateSystem)
regionEnergy.CreateFinish()

# Start the creation of Tissue region
regionTissue = iron.Region()
regionTissue.CreateStart(regionUserNumberTissue,iron.WorldRegion)
regionTissue.LabelSet("Tissue")
regionTissue.CoordinateSystemSet(coordinateSystem)
regionTissue.CreateFinish()

print('\033[1;32m'+'Region            COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Bases
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> BASIS << == ")

# Start the creation of bases for blood energy equaiton
basisXiGaussSpace = 3
basisEnergy = iron.Basis()
basisEnergy.CreateStart(basisUserNumberEnergy)
basisEnergy.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basisEnergy.NumberOfXiSet(numberOfDimensions)
basisEnergy.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE])
basisEnergy.QuadratureNumberOfGaussXiSet([basisXiGaussSpace])
basisEnergy.CreateFinish()

# Create a tri-linear Simplex basis
basisTissue = iron.Basis()
basisTissue.CreateStart(BasisUserNumberTissue)
basisTissue.TypeSet(iron.BasisTypes.SIMPLEX)
basisTissue.NumberOfXiSet(3)
basisTissue.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*3)
basisTissue.QuadratureOrderSet(2)
basisTissue.CreateFinish()

print('\033[1;32m'+'Bases             COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Nodes
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> NODES << == ")
#
# # Start the creation of mesh nodes
NodesEnergy = iron.Nodes()
NodesEnergy.CreateStart(regionEnergy,numberOfNodesEnergy)
NodesEnergy.CreateFinish()

# Define nodes for the tissue mesh
nodesTissue = iron.Nodes()
nodesTissue.CreateStart(regionTissue,numberOfNodesTissue)
nodesTissue.CreateFinish()

print('\033[1;32m'+'Nodes             COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Mesh
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> MESH << == ")


# Start the creation of artery mesh for energy =======
# Create a generated mesh
# generatedMesh = iron.GeneratedMesh()
# generatedMesh.CreateStart(meshUserNumberEnergy,regionEnergy)
# generatedMesh.TypeSet(iron.GeneratedMeshTypes.REGULAR)
# generatedMesh.BasisSet([basisEnergy])
# generatedMesh.ExtentSet([0,LENGTH,0])
# generatedMesh.NumberOfElementsSet([numberGlobalXElements])
#
# meshEnergy = iron.Mesh()
# generatedMesh.CreateFinish(meshUserNumberEnergy,meshEnergy)
meshEnergy = iron.Mesh()
meshEnergy.CreateStart(meshUserNumberEnergy,regionEnergy,3)
meshEnergy.NumberOfElementsSet(numberOfElementsEnergy)
meshEnergy.NumberOfComponentsSet(1)

# Start the creation of mesh elements for energy
meshElementsEnergy  = iron.MeshElements()
meshComponentNumber = 1
meshElementsEnergy.CreateStart(meshEnergy,meshComponentNumber,basisEnergy)

for elemIdx in range(0,numberOfElementsEnergy):
    meshElementsEnergy.NodesSet(elemIdx+1,elementNodes[elemIdx])
meshElementsEnergy.CreateFinish()

meshEnergy.CreateFinish()


# Start the creation of tissue mesh
meshTissue = iron.Mesh()
meshTissue.CreateStart(meshUserNumberTissue,regionTissue,3)
meshTissue.origin=meshOrigin
meshTissue.NumberOfComponentsSet(1)
meshTissue.NumberOfElementsSet(numberOfElementsTissue)

# Start the creation of mesh elements for tissue
meshElementsTissue  = iron.MeshElements()
meshComponentNumber = 1
meshElementsTissue.CreateStart(meshTissue, meshComponentNumber, basisTissue)


print( "Elapsed time before reading ele file is: ", time.time()-t)
# reading elements and its local Nodes and setting elements.Nodes
with open(FileName_ele, "r") as f:
    target=f.readlines()
    for lineNum,line in enumerate(target):
        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
    for lineNum,line in enumerate(target):
        if (lineNum !=0 and lineNum<=numberOfElementsTissue):
# reading elements and its localNodes
            elementNumber = int(target[lineNum][0])+offset
            localNodes[0] = int(target[lineNum][1])+offset
            localNodes[1] = int(target[lineNum][2])+offset
            localNodes[2] = int(target[lineNum][3])+offset
            localNodes[3] = int(target[lineNum][4])+offset
            materialType  = int(target[lineNum][5])
# Giving elements material label
            if materialType == muscleRegionLabel:
                muscleElements.append(elementNumber)
            elif materialType == leftRadiusRegionLabel:
                leftRadiusElements.append(elementNumber)
            elif materialType == leftHumerusRegionLabel:
                leftHumerusElements.append(elementNumber)
            elif materialType == leftUlnaRegionLabel:
                leftUlnaElements.append(elementNumber)
# printing the head and the tail of the elements
            if elementNumber in head and printHead:
              print(elementNumber,localNodes,materialType)
            elif elementNumber in tail and printTail:
              print(elementNumber,localNodes,materialType)

            meshElementsTissue.NodesSet(elementNumber,localNodes)

print( "number of muscle Elements = %d\nnumber of left Ulna elements = %d\ntotal number of elements = %d\n"%(len(muscleElements)
,len(leftUlnaElements),len(muscleElements)+len(leftRadiusElements)+len(leftHumerusElements)+len(leftUlnaElements)))

print( "Elapsed time after reading ele file is: ", time.time()-t)
#input("Press Enter to continue...")

meshElementsTissue.CreateFinish()
meshTissue.CreateFinish()

print('\033[1;32m'+'Mesh              COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Decomposition
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> MESH DECOMPOSITION << == ")

# Create a decomposition for artery energy mesh
decompositionEnergy = iron.Decomposition()
decompositionEnergy.CreateStart(decompositionUserNumberEnergy,meshEnergy)
decompositionEnergy.TypeSet(iron.DecompositionTypes.CALCULATED)
decompositionEnergy.NumberOfDomainsSet(numberOfComputationalNodes)
# decompositionEnergy.CalculateFacesSet(True)
decompositionEnergy.CreateFinish()

# Start the creation of Tissue mesh decomposition =======
decompositionTissue = iron.Decomposition()
decompositionTissue.CreateStart(decompositionUserNumberTissue,meshTissue)
decompositionTissue.TypeSet(iron.DecompositionTypes.CALCULATED)
decompositionTissue.NumberOfDomainsSet(numberOfComputationalNodes)
decompositionTissue.CalculateFacesSet(True)
decompositionTissue.CreateFinish()

print('\033[1;32m'+'Decomposition     COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Geometric Field
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> GEOMETRIC FIELD << == ")

# Start the creation of geometric field for artery energy
geometricFieldEnergy = iron.Field()
geometricFieldEnergy.CreateStart(geometricFieldUserNumberEnergy,regionEnergy)
geometricFieldEnergy.LabelSet('Geometric Field')
geometricFieldEnergy.MeshDecompositionSet(decompositionEnergy)
geometricFieldEnergy.NumberOfVariablesSet(1)
geometricFieldEnergy.VariableLabelSet(iron.FieldVariableTypes.U,'Artery Coordinates')
geometricFieldEnergy.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricFieldEnergy.ScalingTypeSet(iron.FieldScalingTypes.NONE)

for componentNumber in range(1,coordinateSystem.dimension+1):
    geometricFieldEnergy.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentNumber,1)
# geometricFieldEnergy.CreateFinish()
# generatedMesh.GeometricParametersCalculate(geometricFieldEnergy)

geometricFieldEnergy.CreateFinish()

# Set the geometric field values for version 1
versionIdx = 1
derivIdx = 1
for nodeIdx in range(0,numberOfNodesEnergy):
    nodeDomain = decompositionEnergy.NodeDomainGet(nodeIdx+1,1)
    if (nodeDomain == computationalNodeNumber):
        geometricFieldEnergy.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
         versionIdx,derivIdx,nodeIdx+1,1,xValues[nodeIdx][0])
        geometricFieldEnergy.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
         versionIdx,derivIdx,nodeIdx+1,2,yValues[nodeIdx][0])
        geometricFieldEnergy.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
         versionIdx,derivIdx,nodeIdx+1,3,zValues[nodeIdx][0])

# Finish the parameter update
geometricFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)


# Start the creation of geometric field for Tissue
geometricFieldTissue = iron.Field()
geometricFieldTissue.CreateStart(geometricFieldUserNumberTissue,regionTissue)
geometricFieldTissue.LabelSet('Geometric Field')
geometricFieldTissue.MeshDecompositionSet(decompositionTissue)
geometricFieldTissue.NumberOfVariablesSet(1)
geometricFieldTissue.VariableLabelSet(iron.FieldVariableTypes.U,'Tissue Coordinates')
geometricFieldTissue.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricFieldTissue.ScalingTypeSet(iron.FieldScalingTypes.NONE)
for componentNumber in range(1,coordinateSystem.dimension+1):
    geometricFieldTissue.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentNumber,1)
geometricFieldTissue.CreateFinish()



# Get nodes
nodes = iron.Nodes()
regionTissue.NodesGet(nodes)
numberOfNodes = nodes.numberOfNodes
print( numberOfNodes)
# Get or calculate geometric Parameters

print( "Elapsed time before reading node file is: ", time.time()-t)

FileName_node = "MaxVol1000/correctedBoundaryMarker.node"

# X input file is mm. if you want to keep it that way you need to multiply k and rho*c by factors of 10^-3 and 10^-9 respectively.
Units = 1e0
kFactor= 1#1e-3
rhoCFactor = 1#1e-9

printHead = True
printTail = True
head = range(5)
tail = range(numberOfNodes-4,numberOfNodes+1)

#nodePositions = False
boundaryMarker=0
nodeNumber = 0

boundaryTissue = []
skinMarker = 2

with open(FileName_node, "r") as f:
    target=f.readlines()
    for lineNum,line in enumerate(target):
        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
    for lineNum,line in enumerate(target):
        if lineNum !=0:
          nodeNumber = int(target[lineNum][0])+offset
          x = float(target[lineNum][1]) * Units
          y = float(target[lineNum][2]) * Units
          z = float(target[lineNum][3]) * Units
          boundaryMarker = int(target[lineNum][4])
          if boundaryMarker == skinMarker:
            boundaryTissue.append(nodeNumber)

          if nodeNumber in head and printHead:
            print( nodeNumber,x,y,z,boundaryMarker)
          elif nodeNumber in tail and printTail:
            print( nodeNumber,x,y,z,boundaryMarker)

          nodeDomain = decompositionTissue.NodeDomainGet(nodeNumber,1)
          if nodeDomain == computationalNodeNumber:
            geometricFieldTissue.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNumber,1,x)
            geometricFieldTissue.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNumber,2,y)
            geometricFieldTissue.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNumber,3,z)

print( "number of boundary nodes = %d"%len(boundaryTissue))

print( "Elapsed time after reading node file is: ", time.time()-t)

# Update the geometric field
geometricFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

print('\033[1;32m'+'Geometric Field   COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Equations Sets
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> EQUATIONS SET << == ")

# Create the equations set for advection diffusion in arteries
# dT/dt+u dT/dx-alpha d2T/dx2-(b-cT)=0
equationsSetEnergy = iron.EquationsSet()
equationsSetFieldEnergy = iron.Field()
# Set the equations set to be a dynamic linear problem
equationsSetEnergySpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
  iron.EquationsSetTypes.ADVECTION_EQUATION,
  equationsSetEnergySubtype]
equationsSetEnergy.CreateStart(equationsSetUserNumberEnergy,regionEnergy,geometricFieldEnergy,
equationsSetEnergySpecification,equationsSetFieldUserNumberEnergy,equationsSetFieldEnergy)
equationsSetEnergy.LabelSet('Advec Diffusion Equation')
equationsSetEnergy.CreateFinish()


# Create standard Diffusion equations set
# dT/dt-div(alpha grad(T))-(b-cT)=0
equationsSetTissue = iron.EquationsSet()
equationsSetFieldTissue = iron.Field()
equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
  iron.EquationsSetTypes.DIFFUSION_EQUATION,
  equationsSetTissueSubtype]
equationsSetTissue.CreateStart(equationsSetUserNumberTissue,regionTissue,geometricFieldTissue,
        equationsSetSpecification,equationsSetFieldUserNumberTissue,equationsSetFieldTissue)
equationsSetTissue.LabelSet('Diffusion Equation')
equationsSetTissue.CreateFinish()

print('\033[1;32m'+'Equations Set     COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Dependent Field
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> DEPENDENT FIELD << == ")

# Create the equations set dependent field variables
dependentFieldEnergy = iron.Field()
equationsSetEnergy.DependentCreateStart(dependentFieldUserNumberEnergy,dependentFieldEnergy)
dependentFieldEnergy.LabelSet('Dependent Field')
dependentFieldEnergy.VariableLabelSet(iron.FieldVariableTypes.U,'Blood Temperature')
dependentFieldEnergy.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'Blood Temperature Gradient')
dependentFieldEnergy.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
dependentFieldEnergy.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
equationsSetEnergy.DependentCreateFinish()

# Initialise dependent field
dependentFieldEnergy.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,37.0)

dependentFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
dependentFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)


# Create dependent field
dependentFieldTissue = iron.Field()
equationsSetTissue.DependentCreateStart(dependentFieldUserNumberTissue,dependentFieldTissue)
dependentFieldTissue.VariableLabelSet(iron.FieldVariableTypes.U,'Tissue Temperature')
dependentFieldTissue.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'Tissue Temperature Gradient')
dependentFieldTissue.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
dependentFieldTissue.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
equationsSetTissue.DependentCreateFinish()

# Initialise dependent field
dependentFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,37.0)

dependentFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
dependentFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

print('\033[1;32m'+'Dependent Field   COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Materials Field
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> MATERIALS FIELD << == ")

# Create the equations set material field variables for artery energy
materialsFieldEnergy = iron.Field()
equationsSetEnergy.MaterialsCreateStart(materialsFieldUserNumberEnergy,materialsFieldEnergy)
materialsFieldEnergy.VariableLabelSet(iron.FieldVariableTypes.U,'Blood Properties')
materialsFieldEnergy.ComponentLabelSet(iron.FieldVariableTypes.U,1,'Blood Diffusivity')
materialsFieldEnergy.ComponentLabelSet(iron.FieldVariableTypes.U,2,'Source Tb coeff.')
equationsSetEnergy.MaterialsCreateFinish()

# Initialise the properties and source values
diffusivity=Alpha #+U*beta*le/2 #U*beta*le/2=0.000416667 almost 3000 times of the real diffusivity Pe=Ule/2a=0.2*0.05/12/2/0.0004=1
materialsFieldEnergy.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,diffusivity)

for elemIdx in arteriesElements:
    elemDomain = decompositionEnergy.ElementDomainGet(elemIdx)
    if elemDomain == computationalNodeNumber:
        # ra=(radius[elemIdx-1]+radius[elemIdx])/2
        cArtery=4*Alpha/ra[elemIdx-1]**2
        materialsFieldEnergy.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elemIdx,2, cArtery)

materialsFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
materialsFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)



# Create the equations set material field variables for tissue
materialsFieldTissue = iron.Field()
equationsSetTissue.MaterialsCreateStart(materialsFieldUserNumberTissue,materialsFieldTissue)
materialsFieldTissue.VariableLabelSet(iron.FieldVariableTypes.U,'Materials')
materialsFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,1,'Diffusivity 1')
materialsFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,2,'Diffusivity 2')
materialsFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,3,'Diffusivity 3')
materialsFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,4,'Source T coeff.')
equationsSetTissue.MaterialsCreateFinish()


for elementNumber in muscleElements:
    elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
    if elementDomain == computationalNodeNumber:
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,1, k_ms*kFactor/(rho_ms*c_ms*rhoCFactor)) #0.42*1e-3
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,2, k_ms*kFactor/(rho_ms*c_ms*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,3, k_ms*kFactor/(rho_ms*c_ms*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,4, cMuscle)  #4.0e6*1e-9 0.42/4.0e6

for elementNumber in leftRadiusElements:
    elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
    if elementDomain == computationalNodeNumber:
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,1, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,2, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,3, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,4, cBone)

for elementNumber in leftHumerusElements:
    elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
    if elementDomain == computationalNodeNumber:
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,1, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,2, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,3, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,4, cBone)

for elementNumber in leftUlnaElements:
    elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
    if elementDomain == computationalNodeNumber:
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,1, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,2, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,3, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
      materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  elementNumber,4, cBone)

#materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,0.0)

materialsFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
materialsFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

print('\033[1;32m'+'Materials Field   COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Source Field
#================================================================================================================================

#------------------------------
#creating ulnar artery
#-----------------------------

#FileName_cell = "MaxVol1000/ulnarArteryFromWholeUpLi"

#ulnarElementList=[]

#with open(FileName_cell, "r") as f:
#    target=f.readlines()
#    for lineNum,line in enumerate(target):
#        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
#    for lineNum,line in enumerate(target):
#        if lineNum !=0:
#          ulnarElementList.append(int(target[lineNum][0])+1) # plus one because we are taking vtk numbers and we need to add 1.

if (ProgressDiagnostics):
    print( " == >> SOURCE FIELD << == ")

# Create source field for artery energy
sourceFieldEnergy = iron.Field()
equationsSetEnergy.SourceCreateStart(sourceFieldUserNumberEnergy,sourceFieldEnergy)
equationsSetEnergy.SourceCreateFinish()

# Source is b-cT; e.g. for my case C(Tt-T), b=CTt, c=C. Because source field is scalar type I cannot define 2 components.
# So it is in 2nd component of the materials field.

#sourceFieldEnergy.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,cArtery*Tt)
for elemIdx in arteriesElements:
    elemDomain = decompositionEnergy.ElementDomainGet(elemIdx)
    if elemDomain == computationalNodeNumber:
        # ra=(radius[elemIdx-1]+radius[elemIdx])/2
        cArtery=4*Alpha/ra[elemIdx-1]**2
        sourceFieldEnergy.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elemIdx,1,cArtery*Tt[elemIdx])

sourceFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
sourceFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Create source field for tissue
sourceFieldTissue = iron.Field()
equationsSetTissue.SourceCreateStart(sourceFieldUserNumberTissue,sourceFieldTissue)
equationsSetTissue.SourceCreateFinish()

for elementNumber in muscleElements:
    elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
    if elementDomain == computationalNodeNumber:
      sourceFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,700.0e-9/(rho_ms*c_ms)+cMuscle*Tb) # source=qm/rhoC+CTb-CT
for elementNumber in leftRadiusElements:
    elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
    if elementDomain == computationalNodeNumber:
      sourceFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,0.0)
for elementNumber in leftHumerusElements:
    elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
    if elementDomain == computationalNodeNumber:
      sourceFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,0.0)
for elementNumber in leftUlnaElements:
    elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
    if elementDomain == computationalNodeNumber:
      sourceFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,0.0)

#sourceField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,700.0/4.0e6)
#for elementNumber in ulnarElementList:
#    elementDomain = decomposition.ElementDomainGet(elementNumber)
#    if elementDomain == computationalNodeNumber:
##      sourceField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,2.5e6/2.0/4.0e6)
#      sourceField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,0.0)
# Modifying source term for the element with the artery passing through.
# for elementNumber in range(1,81):
#   sourceFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,
#     700/(rho_t*c_t)+CWall*Tb*volumeCorrection) # Because element volume here is much smaller than the artery passing it. We could consider more elements as well.

sourceFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
sourceFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

print('\033[1;32m'+'Source Field      COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
# Independent Field
#================================================================================================================================

if (ProgressDiagnostics):
    print (" == >> INDEPENDENT FIELD << == ")

# Create the equations set independent field variables
IndependentFieldEnergy = iron.Field()
#IndependentFieldEnergy.VariableLabelSet(iron.FieldVariableTypes.U,'flow velocity')
# Set the mesh component to be used by the field components.
#IndependentFieldEnergy.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)

# NAVIER-STOKES
equationsSetEnergy.IndependentCreateStart(independentFieldUserNumber,IndependentFieldEnergy)
equationsSetEnergy.IndependentCreateFinish()

# Set the velocity
IndependentFieldEnergy.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  1,1.0)

# Finish the parameter update
IndependentFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
IndependentFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

print('\033[1;32m'+'Independent Field COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Equations
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> EQUATIONS << == ")

# Create equations for artery energy
equationsEnergy = iron.Equations()
equationsSetEnergy.EquationsCreateStart(equationsEnergy)
equationsEnergy.sparsityType = iron.EquationsSparsityTypes.SPARSE
equationsEnergy.outputType = iron.EquationsOutputTypes.NONE
equationsSetEnergy.EquationsCreateFinish()

  # I want to solve this type of equation, dT/dt-div(Sigma grad(T))-(b-cT)=0.
  # Sigma in my case is Sigma=k/rhoC.
  # q=(Sigma grad(T)).n which in my 1D case is q=k/rhoC dT/dx which in the boundary is q=-h/rhoC (T-Tinf)
  # in my case b=qm/rhoC+CTb and c=C

# Create equations for tissue
equationsTissue = iron.Equations()
equationsSetTissue.EquationsCreateStart(equationsTissue)
equationsTissue.sparsityType = iron.EquationsSparsityTypes.SPARSE
equationsTissue.outputType = iron.EquationsOutputTypes.NONE
equationsSetTissue.EquationsCreateFinish()

print('\033[1;32m'+'equations         COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Problems
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> PROBLEM << == ")

# Start the creation of a problem.
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.MULTI_PHYSICS,
                        iron.ProblemTypes.DIFFUSION_ADVECTION_DIFFUSION,
                        ProblemSubtype]
problem.CreateStart(problemUserNumber,problemSpecification)
problem.CreateFinish()

print('\033[1;32m'+'Problems          COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Control Loops
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> PROBLEM CONTROL LOOP << == ")

# Create control loops
problem.ControlLoopCreateStart()
TimeLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],TimeLoop)
TimeLoop.LabelSet('Time Loop')
TimeLoop.TimesSet(startTime,stopTime,timeIncrement)
TimeLoop.TimeOutputSet(DYNAMIC_SOLVER_DIFFUSION_OUTPUT_FREQUENCY)
problem.ControlLoopCreateFinish()

print('\033[1;32m'+'Control Loops     COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Solvers
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> SOLVERS << == ")

# Create problem solver
solverEnergy = iron.Solver()
LinearSolverEnergy = iron.Solver()

solverTissue = iron.Solver()
LinearSolverTissue = iron.Solver()

problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solverEnergy)
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,solverTissue)
solverEnergy.LabelSet('Arterial Energy Solver')
solverTissue.LabelSet('Tissue Energy Solver')
#solver.outputType = iron.SolverOutputTypes.SOLVER
solverEnergy.DynamicLinearSolverGet(LinearSolverEnergy)
solverTissue.DynamicLinearSolverGet(LinearSolverTissue)
#solver.linearType = iron.LinearSolverTypes.ITERATIVE
#solver.linearIterativeAbsoluteTolerance = 1.0E-12
#solver.linearIterativeRelativeTolerance = 1.0E-12
problem.SolversCreateFinish()

print('\033[1;32m'+'Solvers           COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Solver Equations
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> SOLVER EQUATIONS << == ")

# Create solver equations and add equations set to solver equations
solverEnergy = iron.Solver()
solverTissue = iron.Solver()
solverEquationsEnergy = iron.SolverEquations()
solverEquationsTissue = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solverEnergy)
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,solverTissue)
solverEnergy.SolverEquationsGet(solverEquationsEnergy)
solverTissue.SolverEquationsGet(solverEquationsTissue)
solverEquationsEnergy.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex  = solverEquationsEnergy.EquationsSetAdd(equationsSetEnergy)
equationsSetIndex2 = solverEquationsTissue.EquationsSetAdd(equationsSetTissue)
problem.SolverEquationsCreateFinish()

print('\033[1;32m'+'Solver Equations  COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> BOUNDARY CONDITIONS << == ")

boundaryConditionsEnergy = iron.BoundaryConditions()
boundaryConditionsTissue = iron.BoundaryConditions()

solverEquationsEnergy.BoundaryConditionsCreateStart(boundaryConditionsEnergy)

nodeDomain = decompositionEnergy.NodeDomainGet(1,1)
if nodeDomain == computationalNodeNumber:
    boundaryConditionsEnergy.SetNode(dependentFieldEnergy,iron.FieldVariableTypes.U,1,1,1,1,
        iron.BoundaryConditionsTypes.FIXED,[37.0])

dependentFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
dependentFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

solverEquationsEnergy.BoundaryConditionsCreateFinish()

#Tissue boundary conditions =============
solverEquationsTissue.BoundaryConditionsCreateStart(boundaryConditionsTissue)

nodes = iron.Nodes()
regionTissue.NodesGet(nodes)


# Actually in OpenCMISS the equation is divided by rhoC. So q=alpha*gradT.n which alpha=sigma/rhoC.
# So for Robin BCs you need to pass h/rhoC and q_h/rhoC. Also DON'T FORGET ABOUT THE UNITS
q_hUnit = 1#1e-6
hUnit   = 1#1e-6


Rtot=1/((h_conv+hr_rad)*hUnit)+3/k_sk # units terms has mm2 so no 1e6.

for nodeNumber in boundaryTissue:
  nodeDomain = decompositionTissue.NodeDomainGet(nodeNumber,1)
  if nodeDomain == computationalNodeNumber:
    boundaryConditionsTissue.SetNode(dependentFieldTissue,iron.FieldVariableTypes.DELUDELN,1,1,nodeNumber,1,
      iron.BoundaryConditionsTypes.ROBIN,[1.0/(rho_sk*c_sk*rhoCFactor*Rtot),1.0/(rho_sk*c_sk*rhoCFactor*Rtot)* Tair])


dependentFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
dependentFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
#!  !Finish the creation of the equations set boundary conditions
#!  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

solverEquationsTissue.BoundaryConditionsCreateFinish()

print('\033[1;32m'+'Boundary Conditions COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Run Solvers
#================================================================================================================================

#!  !Solve the problem
#!  CALL cmfe_Problem_Solve(Problem,Err)
print( "Solving problem...")
start = time.time()
# Solve the problem
problem.Solve()
end = time.time()
elapsed = end - start
print( "Total Number of Elements = %d " %totalNumberOfElements)
print( "Total Number of nodes = %d " %numberOfNodes)
print( "Calculation Time = %3.4f" %elapsed)
print( "Problem solved!")
print( "#")
print( "number of muscle Elements = %d\nnumber of left radius elements = %d\nHumerus=%d\nUlna=%d\ntotal number of elements = %d\n"%(
len(muscleElements),len(leftRadiusElements),len(leftHumerusElements),len(leftUlnaElements),len(muscleElements)+
len(leftRadiusElements)+len(leftHumerusElements)+len(leftUlnaElements)))
print( "number of boundary nodes = %d"%len(boundary))
print( '\033[1;31m'+"Elapsed time: "+'\033[0m', time.time()-t)
#!# Export results
#!baseName = "laplace"
#!dataFormat = "PLAIN_TEXT"
#!fml = iron.FieldMLIO()
#!fml.OutputCreate(mesh, "", baseName, dataFormat)
#!fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,
#!    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#!fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,
#!    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#!fml.OutputWrite("LaplaceExample.xml")
#!fml.Finalise()


#!  !Output Analytic analysis
#!!  CALL cmfe_AnalyticAnalysis_Output(DependentField,"DiffusionAnalytics_x10_y10_z10_L_T1",Err)
#exportField = True
#if (exportField):
#    fields = iron.Fields()
#    fields.Create(Region)
#    fields.Finalise()
#!  EXPORT_FIELD=.TRUE.
#!  IF(EXPORT_FIELD) THEN
#!    CALL cmfe_Fields_Initialise(Fields,Err)
#!    CALL cmfe_Fields_Create(Region,Fields,Err)
#!!    CALL cmfe_Fields_NodesExport(Fields,"DiffusionConstantSourceAnalytic_x10_y10_z10_L_T1","FORTRAN",Err)
#!!    CALL cmfe_Fields_ElementsExport(Fields,"DiffusionConstantSourceAnalytic_x10_y10_z10_L_T1","FORTRAN",Err)
#!    CALL cmfe_Fields_Finalise(Fields,Err)

#!  ENDIF

#!  !CALL cmfe_Finalise(Err)
#!  WRITE(*,'(A)') "Program successfully completed."

#!  STOP
iron.Finalise()
#!END PROGRAM DIFFUSIONCONSTANTSOURCEEXAMPLE
