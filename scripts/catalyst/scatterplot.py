
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.1.2-399-g6324cb6 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.1.2-399-g6324cb6

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1174, 929]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [15.234395240433514, -4.76837158203125e-07, 0.11717867851257324]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [17.27228783074267, 71.63250368971187, -9.285649562306032]
      renderView1.CameraFocalPoint = [15.234395240433512, -4.768371590034197e-07, 0.11717867851257469]
      renderView1.CameraViewUp = [-0.8792420791071403, -0.03732527606443192, -0.474910718023996]
      renderView1.CameraParallelScale = 18.706336660411555
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='electron_scatterplot_%t.png', freq=15, fittoscreen=0, magnification=1, width=1174, height=929, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Partitioned Polydata Reader'
      # create a producer from a simulation input
      electron = coprocessor.CreateProducer(datadescription, 'electron')

      # create a new 'Glyph'
      glyph1 = Glyph(Input=electron,
          GlyphType='Arrow')
      glyph1.Scalars = ['POINTS', 'None']
      glyph1.Vectors = ['POINTS', 'momentum']
      glyph1.ScaleMode = 'scalar'
      glyph1.ScaleFactor = 3.046867900155485
      glyph1.MaximumNumberOfSamplePoints = 30000
      glyph1.Stride = 200
      glyph1.GlyphTransform = 'Transform2'

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'weights'
      weightsLUT = GetColorTransferFunction('weights')
      weightsLUT.RGBPoints = [0.002574920654296875, 0.231373, 0.298039, 0.752941, 0.002582550048828125, 0.865003, 0.865003, 0.865003, 0.002590179443359375, 0.705882, 0.0156863, 0.14902]
      weightsLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'weights'
      weightsPWF = GetOpacityTransferFunction('weights')
      weightsPWF.Points = [0.002574920654296875, 0.0, 0.5, 0.0, 0.002590179443359375, 1.0, 0.5, 0.0]
      weightsPWF.ScalarRangeInitialized = 1

      # get color transfer function/color map for 'GlyphVector'
      glyphVectorLUT = GetColorTransferFunction('GlyphVector')
      glyphVectorLUT.RGBPoints = [0.009482396371128719, 0.231373, 0.298039, 0.752941, 0.5585663030005522, 0.865003, 0.865003, 0.865003, 1.1076502096299758, 0.705882, 0.0156863, 0.14902]
      glyphVectorLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'GlyphVector'
      glyphVectorPWF = GetOpacityTransferFunction('GlyphVector')
      glyphVectorPWF.Points = [0.009482396371128719, 0.0, 0.5, 0.0, 1.1076502096299758, 1.0, 0.5, 0.0]
      glyphVectorPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from electron
      electronDisplay = Show(electron, renderView1)
      # trace defaults for the display properties.
      electronDisplay.ColorArrayName = ['POINTS', 'weights']
      electronDisplay.LookupTable = weightsLUT
      electronDisplay.OSPRayScaleArray = 'weights'
      electronDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      electronDisplay.SelectOrientationVectors = 'momentum'
      electronDisplay.ScaleFactor = 3.046867900155485
      electronDisplay.SelectScaleArray = 'weights'
      electronDisplay.GlyphType = 'Arrow'

      # show color legend
      electronDisplay.SetScalarBarVisibility(renderView1, True)

      # show data from glyph1
      glyph1Display = Show(glyph1, renderView1)
      # trace defaults for the display properties.
      glyph1Display.Representation = 'Points'
      glyph1Display.ColorArrayName = ['POINTS', 'GlyphVector']
      glyph1Display.LookupTable = glyphVectorLUT
      glyph1Display.Opacity = 0.74
      glyph1Display.OSPRayScaleArray = 'weights'
      glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      glyph1Display.SelectOrientationVectors = 'GlyphVector'
      glyph1Display.ScaleFactor = 3.349505639076233
      glyph1Display.SelectScaleArray = 'weights'
      glyph1Display.GlyphType = 'Sphere'

      # show color legend
      glyph1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for weightsLUT in view renderView1
      weightsLUTColorBar = GetScalarBar(weightsLUT, renderView1)
      weightsLUTColorBar.Title = 'weights'
      weightsLUTColorBar.ComponentTitle = ''

      # get color legend/bar for glyphVectorLUT in view renderView1
      glyphVectorLUTColorBar = GetScalarBar(glyphVectorLUT, renderView1)
      glyphVectorLUTColorBar.Position = [0.85, 0.52]
      glyphVectorLUTColorBar.Position2 = [0.12, 0.42999999999999994]
      glyphVectorLUTColorBar.Title = 'GlyphVector'
      glyphVectorLUTColorBar.ComponentTitle = 'Magnitude'

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(glyph1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'electron': [15, 15]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(False, 1)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
