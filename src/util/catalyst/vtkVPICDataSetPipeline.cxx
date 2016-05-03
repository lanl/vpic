#include "vtkVPICDataSetPipeline.h"

#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCameraPass.h"
#include "vtkColorTransferFunction.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkDataSetMapper.h"
#include "vtkExecutive.h"
#include "vtkLightKit.h"
#include "vtkLightsPass.h"
#include "vtkIceTCompositePass.h"
#include "vtkIceTSynchronizedRenderers.h"
#include "vtkInformation.h"
#include "vtkMath.h"
#include "vtkMPIController.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkOpaquePass.h"
#include "vtkOpenGLRenderer.h"
#include "vtkPNGWriter.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderPassCollection.h"
#include "vtkRenderWindow.h"
#include "vtkSequencePass.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkSynchronizedRenderWindows.h"
#include "vtkWindowToImageFilter.h"

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


vtkStandardNewMacro(vtkVPICDataSetPipeline);

vtkVPICDataSetPipeline::vtkVPICDataSetPipeline() 
{
  this->dataSetMapper = vtkDataSetMapper::New();
  this->actor->SetMapper(this->dataSetMapper);
}

vtkVPICDataSetPipeline::~vtkVPICDataSetPipeline()
{
  this->dataSetMapper->Delete();
}

int vtkVPICDataSetPipeline::CoProcess(vtkCPDataDescription* desc)
{
  this->dataSetMapper->SetInputData(vtkDataSet::SafeDownCast(desc->GetInputDescriptionByName("input")->GetGrid()));
  this->dataSetMapper->SetScalarRange(AttributeMinimum, AttributeMaximum);

  vtkVPICPipeline::CoProcess(desc);
  return 1;
}
