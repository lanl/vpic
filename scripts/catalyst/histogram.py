
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

      # Create a new 'Bar Chart View'
      barChartView1 = CreateView('XYBarChartView')
      barChartView1.ViewSize = [915, 538]
      barChartView1.LeftAxisRangeMaximum = 350000.0
      barChartView1.BottomAxisRangeMinimum = -1.2000000000000002
      barChartView1.BottomAxisRangeMaximum = 1.2000000000000002

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(barChartView1,
          filename='ion_histogram_%t.png', freq=10, fittoscreen=0, magnification=1, width=915, height=538, cinema={})
      barChartView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Partitioned Polydata Reader'
      # create a producer from a simulation input
      ion = coprocessor.CreateProducer(datadescription, 'ion')

      # create a new 'Histogram'
      histogram1 = Histogram(Input=ion)
      histogram1.SelectInputArray = ['POINTS', 'momentum']
      histogram1.BinCount = 30
      histogram1.CustomBinRanges = [0.00257492065429688, 0.00257492065429688]

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
      # finally, restore active source
      SetActiveSource(histogram1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'ion': [10]}
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
