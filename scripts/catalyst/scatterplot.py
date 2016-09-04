
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.1.2 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.1.2

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [791, 512]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [15.234381385147572, -6.9141387939453125e-06, 0.23440146446228027]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [10.644077829285354, 86.09566389824528, -18.05071227374809]
      renderView1.CameraFocalPoint = [15.234381385147582, -6.9141387839099245e-06, 0.23440146446226215]
      renderView1.CameraViewUp = [0.9801362825393347, 0.010192348917000493, -0.198063080036456]
      renderView1.CameraParallelScale = 18.852204262064895
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=791, height=512, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Partitioned Polydata Reader'
      # create a producer from a simulation input
      electron_pvtp = coprocessor.CreateProducer(datadescription, 'electron')

      # create a new 'Glyph'
      glyph1 = Glyph(Input=electron_pvtp,
          GlyphType='Sphere')
      glyph1.Scalars = ['POINTS', 'None']
      glyph1.Vectors = ['POINTS', 'momentum']
      glyph1.ScaleFactor = 3.0468695268034938
      glyph1.MaximumNumberOfSamplePoints = 50000
      glyph1.GlyphTransform = 'Transform2'

      # init the 'Sphere' selected for 'GlyphType'
      glyph1.GlyphType.Radius = 0.05

      # create a new 'Parallel PolyData Writer'
      parallelPolyDataWriter1 = servermanager.writers.XMLPPolyDataWriter(Input=glyph1)

      # register the writer with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the data, etc.
      coprocessor.RegisterWriter(parallelPolyDataWriter1, filename='filename_%t.pvtp', freq=1)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'GlyphVector'
      glyphVectorLUT = GetColorTransferFunction('GlyphVector')
      glyphVectorLUT.RGBPoints = [0.00900907972479836, 0.231373, 0.298039, 0.752941, 0.5441682644121746, 0.865003, 0.865003, 0.865003, 1.0793274490995508, 0.705882, 0.0156863, 0.14902]
      glyphVectorLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'GlyphVector'
      glyphVectorPWF = GetOpacityTransferFunction('GlyphVector')
      glyphVectorPWF.Points = [0.00900907972479836, 0.0, 0.5, 0.0, 1.0793274490995508, 1.0, 0.5, 0.0]
      glyphVectorPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from electron_pvtp
      electron_pvtpDisplay = Show(electron_pvtp, renderView1)
      # trace defaults for the display properties.
      electron_pvtpDisplay.ColorArrayName = [None, '']
      electron_pvtpDisplay.OSPRayScaleArray = 'momentum'
      electron_pvtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      electron_pvtpDisplay.GlyphType = 'Arrow'
      electron_pvtpDisplay.SetScaleArray = [None, '']
      electron_pvtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      electron_pvtpDisplay.OpacityArray = [None, '']
      electron_pvtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'

      # show data from glyph1
      glyph1Display = Show(glyph1, renderView1)
      # trace defaults for the display properties.
      glyph1Display.ColorArrayName = ['POINTS', 'GlyphVector']
      glyph1Display.LookupTable = glyphVectorLUT
      glyph1Display.OSPRayScaleArray = 'GlyphVector'
      glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      glyph1Display.GlyphType = 'Arrow'
      glyph1Display.SetScaleArray = [None, '']
      glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
      glyph1Display.OpacityArray = [None, '']
      glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      glyph1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for glyphVectorLUT in view renderView1
      glyphVectorLUTColorBar = GetScalarBar(glyphVectorLUT, renderView1)
      glyphVectorLUTColorBar.Title = 'GlyphVector'
      glyphVectorLUTColorBar.ComponentTitle = 'Magnitude'

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(parallelPolyDataWriter1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'electron': [5]}
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
