#include "vtkVPICPistonContourPipeline.h"

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
#include "vtkDataSetToPiston.h"
#include "vtkDebugLeaks.h"
#include "vtkPistonContour.h"
#include "vtkPistonMapper.h"
#include "vtkTrivialProducer.h"

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


vtkStandardNewMacro(vtkVPICPistonContourPipeline);

vtkVPICPistonContourPipeline::vtkVPICPistonContourPipeline() 
{
  this->d2p = vtkDataSetToPiston::New();
  this->input = vtkTrivialProducer::New ();
  this->contour = vtkPistonContour::New();
  this->mapper = vtkPistonMapper::New();
  this->Isovalue = 0.0;

  vtkMultiProcessController *ctrl = vtkMultiProcessController::GetGlobalController();
  if (ctrl->GetLocalProcessId() == 0) 
  {
    this->window->Render();
    vtkPistonMapper::InitCudaGL(this->window);
    ctrl->TriggerBreakRMIs();
    ctrl->Barrier();
  } 
  else 
  {
    ctrl->ProcessRMIs();
    ctrl->Barrier();
  }
}

vtkVPICPistonContourPipeline::~vtkVPICPistonContourPipeline()
{
  this->d2p->Delete();
  this->input->Delete();
  this->contour->Delete();
  this->mapper->Delete();
}

int vtkVPICPistonContourPipeline::CoProcess(vtkCPDataDescription* desc)
{
  this->input->SetOutput(desc->GetInputDescriptionByName("input")->GetGrid());
  this->d2p->SetInputConnection(input->GetOutputPort()); 
  this->contour->SetInputConnection(d2p->GetOutputPort());
  this->contour->SetIsoValue(this->Isovalue); 
  this->mapper->SetInputConnection(contour->GetOutputPort());
  this->mapper->Update(); 
  this->actor->SetMapper(this->mapper);

  vtkVPICPipeline::CoProcess(desc);
  return 1;
}
