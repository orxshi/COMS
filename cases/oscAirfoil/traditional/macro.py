from paraview.simple import *

def createVTK (fileName, mainView):
	vtk = LegacyVTKReader(FileNames=[fileName])
	vtkDis = Show(vtk, mainView)
	Hide(vtk, mainView)
	vtkDis.SetRepresentationType('Surface With Edges')
	return [vtk, vtkDis]

def createThreshold (inpt, var, minVal, maxVal, mainView):
	th = Threshold(Input=inpt)
	th.Scalars = ['CELLS', var]
	th.ThresholdRange = [minVal, maxVal]
	thDis = Show(th, mainView)
	thDis.SetRepresentationType('Surface With Edges')
	return th

paraview.simple._DisableFirstRenderCameraReset()

mainView = GetActiveViewOrCreate('RenderView')
mainView.InteractionMode = '2D'

fileName1 = '/home/orhan/Documents/codeBaku/cases/oscAirfoil/UnifiedGrid/out/10_6_2015_15_41/Grid_0/allVTK_0.vtk' 
fileName2 = '/home/orhan/Documents/codeBaku/cases/oscAirfoil/UnifiedGrid/out/10_6_2015_15_41/Grid_1/allVTK_0.vtk'
fileName3 = '/home/orhan/Documents/codeBaku/cases/oscAirfoil/UnifiedGrid/out/10_6_2015_15_41/tri.vtk'

[vtk1,vtk1Dis] = createVTK(fileName1, mainView)
[vtk2,vtk2Dis] = createVTK(fileName2, mainView)
[vtk3,vtk3Dis] = createVTK(fileName3, mainView)

Show(vtk3, mainView)
vtk3Res = GetDisplayProperties(vtk3, mainView)
ColorBy(vtk3Res, None)
vtk3Res.DiffuseColor = [255, 0, 0]

th1 = createThreshold (vtk1, 'iBlank', 4, 4, mainView)
th2 = createThreshold (th1, 'trim', 0, 0, mainView)
Hide(th1, mainView)

th3 = createThreshold (vtk2, 'iBlank', 4, 4, mainView)
th4 = createThreshold (th3, 'trim', 0, 0, mainView)
#th4 = Threshold (Input = th3, Scalars = 'trim')
#th4.ThresholdRange = [0.0, 0.0]
Hide(th3, mainView)

#th5 = createThreshold (vtk3, 'I', 0, 0, mainView)
#extractEdges1 = ExtractEdges (Input=th5)
#th6 = createThreshold (extractEdges1, 'I', 0, 1, mainView)

#iss = IDSelectionSource (FieldType = 'CELL', InsideOut = False)
#idds = GenerateIds(Input=vtk3, ArrayName = 'pts')
#iss.IDs = idds.pts
#extractSelection1 = ExtractSelection(Input=vtk3, Selection=iss)

mainView.ResetCamera()
