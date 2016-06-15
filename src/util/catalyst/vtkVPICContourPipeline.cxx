#include "vtkVPICContourPipeline.h"

#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCameraPass.h"
#include "vtkColorTransferFunction.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataMapper.h"
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


vtkStandardNewMacro(vtkVPICContourPipeline);

vtkVPICContourPipeline::vtkVPICContourPipeline() 
{
  this->contours = vtkContourFilter::New();
  this->polyDataMapper = vtkPolyDataMapper::New();

  this->polyDataMapper->SetInputData(contours->GetOutput());
  this->actor->SetMapper(this->polyDataMapper); 

  this->Isovalue = 0.0;
}

vtkVPICContourPipeline::~vtkVPICContourPipeline()
{
  this->contours->Delete();
  this->polyDataMapper->Delete();
}

int vtkVPICContourPipeline::CoProcess (vtkCPDataDescription* desc)
{
  this->contours->SetInputData(vtkDataSet::SafeDownCast(desc->GetInputDescriptionByName("input")->GetGrid()));
  this->contours->SetValue(0, this->Isovalue);
  this->contours->Update();

  vtkVPICPipeline::CoProcess(desc);
  return 1;
}
