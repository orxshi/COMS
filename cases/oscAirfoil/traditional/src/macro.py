>>> from paraview.simple import *
# Create a cone and assign it as the active object
>>> Cone()
<paraview.servermanager.Cone object at 0x2910f0>
# Set a property of the active object
>>> SetProperties(Resolution=32)
# Apply the shrink filter to the active object
# Shrink is now active
>>> Shrink() 
<paraview.servermanager.Shrink object at 0xaf64050>
# Show shrink
>>> Show() 
<paraview.servermanager.UnstructuredGridRepresentation object at 0xaf57f90>
# Render the active view
>>> Render() 
<paraview.servermanager.RenderView object at 0xaf57ff0>
