
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

      # Create a new 'Bar Chart View'
      barChartView1 = CreateView('XYBarChartView')
      barChartView1.ViewSize = [598, 942]
      barChartView1.LeftAxisRangeMaximum = 90000.0
      barChartView1.BottomAxisRangeMinimum = -1.0
      barChartView1.BottomAxisRangeMaximum = 1.2000000000000002

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(barChartView1,
          filename='image_1_%t.png', freq=1, fittoscreen=0, magnification=1, width=598, height=942, cinema={})
      barChartView1.ViewTime = datadescription.GetTime()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [599, 942]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [15.234341884031892, 7.867813110351562e-06, 0.23438167572021484]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [15.234341884031892, 7.867813110351562e-06, 73.07345914174965]
      renderView1.CameraFocalPoint = [15.234341884031892, 7.867813110351562e-06, 0.23438167572021484]
      renderView1.CameraParallelScale = 18.852140475906275
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_0_%t.png', freq=1, fittoscreen=0, magnification=1, width=599, height=942, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Partitioned Polydata Reader'
      # create a producer from a simulation input
      ion_10pvtp = coprocessor.CreateProducer(datadescription, 'ion')

      # create a new 'Histogram'
      histogram1 = Histogram(Input=ion_10pvtp)
      histogram1.SelectInputArray = ['POINTS', 'momentum']
      histogram1.BinCount = 30
      histogram1.CustomBinRanges = [0.0102996826171875, 0.0102996826171875]

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'weights'
      weightsLUT = GetColorTransferFunction('weights')
      weightsLUT.RGBPoints = [0.0102996826171875, 0.231373, 0.298039, 0.752941, 0.01029973411611557, 0.865003, 0.865003, 0.865003, 0.01029978561504364, 0.705882, 0.0156863, 0.14902]
      weightsLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'weights'
      weightsPWF = GetOpacityTransferFunction('weights')
      weightsPWF.Points = [0.0102996826171875, 0.0, 0.5, 0.0, 0.01029978561504364, 1.0, 0.5, 0.0]
      weightsPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'barChartView1'
      # ----------------------------------------------------------------

      # show data from histogram1
      histogram1Display = Show(histogram1, barChartView1)
      # trace defaults for the display properties.
      histogram1Display.CompositeDataSetIndex = [0]
      histogram1Display.AttributeType = 'Row Data'
      histogram1Display.UseIndexForXAxis = 0
      histogram1Display.XArrayName = 'bin_extents'
      histogram1Display.SeriesVisibility = ['bin_values']
      histogram1Display.SeriesLabel = ['bin_extents', 'bin_extents', 'bin_values', 'bin_values']
      histogram1Display.SeriesColor = ['bin_extents', '0', '0', '0', 'bin_values', '0.89', '0.1', '0.11']
      histogram1Display.SeriesPlotCorner = ['bin_extents', '0', 'bin_values', '0']

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from ion_10pvtp
      ion_10pvtpDisplay = Show(ion_10pvtp, renderView1)
      # trace defaults for the display properties.
      ion_10pvtpDisplay.ColorArrayName = ['POINTS', 'weights']
      ion_10pvtpDisplay.LookupTable = weightsLUT
      ion_10pvtpDisplay.OSPRayScaleArray = 'weights'
      ion_10pvtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      ion_10pvtpDisplay.GlyphType = 'Arrow'
      ion_10pvtpDisplay.SetScaleArray = ['POINTS', 'weights']
      ion_10pvtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      ion_10pvtpDisplay.OpacityArray = ['POINTS', 'weights']
      ion_10pvtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      ion_10pvtpDisplay.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for weightsLUT in view renderView1
      weightsLUTColorBar = GetScalarBar(weightsLUT, renderView1)
      weightsLUTColorBar.Title = 'weights'
      weightsLUTColorBar.ComponentTitle = ''

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(histogram1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'ion': [1, 1]}
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
