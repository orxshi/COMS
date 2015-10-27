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

fileName1 = '/home/orhan/Documents/codeBaku/cases/oscAirfoil/UnifiedGrid/out/10_3_2015_4_40/Grid_0/allVTK_0.vtk' 
fileName2 = '/home/orhan/Documents/codeBaku/cases/oscAirfoil/UnifiedGrid/out/10_3_2015_4_40/Grid_1/allVTK_0.vtk'

[vtk1,vtk1Dis] = createVTK(fileName1, mainView)
[vtk2,vtk2Dis] = createVTK(fileName2, mainView)

th1 = createThreshold (vtk1, 'iBlank', 4, 4, mainView)
th2 = createThreshold (th1, 'trim', 0, 0, mainView)
Hide(th1, mainView)

th3 = createThreshold (vtk2, 'iBlank', 4, 4, mainView)
th4 = createThreshold (th3, 'trim', 0, 0, mainView)
Hide(th3, mainView)

mainView.ResetCamera()
