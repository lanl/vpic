#include "vtkVPICPipeline.h"

#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCameraPass.h"
#include "vtkColorTransferFunction.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
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

vtkStandardNewMacro(vtkVPICPipeline);

vtkVPICPipeline::vtkVPICPipeline() 
{
  this->Filename = 0;

  this->CameraThetaAngle = 0.0;
  this->CameraPhiAngle = 0.0;
  this->CameraDistance = 0.0;

  for (int i=0; i<3; i++)
  {
    this->Bounds[2*i] = -1.0;
    this->Bounds[2*i+1] = 1.0;
  }

  this->AttributeMaximum = 1.0;
  this->AttributeMinimum = 0.0;

  this->actor = vtkActor::New();
  this->lightKit = vtkLightKit::New();
  this->syncRen = vtkIceTSynchronizedRenderers::New();
  this->syncWin = vtkSynchronizedRenderWindows::New();
  this->window = vtkRenderWindow::New();
  this->renderer = vtkRenderer::New();
  this->w2i = vtkWindowToImageFilter::New();
  this->writer = vtkPNGWriter::New();

  vtkMultiProcessController *ctrl = vtkMultiProcessController::GetGlobalController();
 
  this->renderer->SetBackground(0.0, 0.0, 0.0);
  this->renderer->AddActor(actor);

  this->window->AddRenderer(renderer);
  this->window->SetPosition(ctrl->GetLocalProcessId()*(5), 0);
  this->window->SetSize (1024, 1024); 
  this->window->DoubleBufferOn();
  this->window->SwapBuffersOff();

  this->lightKit->AddLightsToRenderer(renderer);

  VTK_CREATE (vtkRenderPassCollection, passes);
  VTK_CREATE (vtkLightsPass, lights);
  passes->AddItem(lights);
  VTK_CREATE (vtkOpaquePass, opaque);
  passes->AddItem(opaque);

  VTK_CREATE (vtkSequencePass, seq);
  seq->SetPasses(passes);  

  VTK_CREATE (vtkIceTCompositePass, iceTPass);
  iceTPass->SetController(ctrl);
  iceTPass->SetRenderPass(seq);
  iceTPass->SetDataReplicatedOnAllProcesses(false);
  iceTPass->SetFixBackground(true);

  VTK_CREATE (vtkCameraPass, cameraP);
  cameraP->SetDelegatePass(iceTPass);
  vtkOpenGLRenderer *glRenderer = vtkOpenGLRenderer::SafeDownCast(this->renderer);
  if( glRenderer != NULL )
  {
    glRenderer->SetPass(cameraP);
  }
  else
  {
    vtkErrorMacro("Cannot cast renderer to vtkOpenGLRenderer!");
    return;
  }

  this->syncWin->SetRenderWindow(window);
  this->syncWin->SetParallelController(ctrl);
  this->syncWin->SetIdentifier(rand()); 

  this->syncRen->SetRenderer(renderer);
  this->syncRen->SetParallelController(ctrl);
  this->syncRen->SetParallelRendering(true);
  this->syncRen->WriteBackImagesOn();
  this->syncRen->SetRootProcessId(0);

  this->w2i->SetInput(window);
  this->w2i->ReadFrontBufferOff();
  this->w2i->ShouldRerenderOff();
  this->w2i->FixBoundaryOn();

  this->writer->SetInputConnection(w2i->GetOutputPort());
}

vtkVPICPipeline::~vtkVPICPipeline()
{
  this->actor->Delete();
  this->lightKit->Delete();
  this->syncRen->Delete();
  this->syncWin->Delete();
  this->window->Delete();
  this->renderer->Delete();
  this->w2i->Delete();
  this->writer->Delete();
}

void vtkVPICPipeline::PrintSelf(ostream& os, vtkIndent indent) 
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Filename: " << this->Filename << endl;
  os << indent << "CameraThetaAngle: " << this->CameraThetaAngle << endl;
  os << indent << "CameraPhiAngle: " << this->CameraPhiAngle << endl;
  os << indent << "CameraDistance: " << this->CameraDistance << endl;
  os << indent << "AttributeMaximum: " << this->AttributeMaximum << endl;
  os << indent << "AttributeMinimum: " << this->AttributeMinimum << endl;
}

int vtkVPICPipeline::RequestDataDescription(vtkCPDataDescription* desc)
{
  if (!desc) 
  {
    vtkWarningMacro("data description is NULL"); 
    return 0;
  }
  return (desc->GetTimeStep() % 10 == 0);
}

int vtkVPICPipeline::CoProcess (vtkCPDataDescription* desc)
{
  int timestep = desc->GetTimeStep();
  vtkMultiProcessController *ctrl = vtkMultiProcessController::GetGlobalController();

  vtkCamera *cam = this->renderer->GetActiveCamera();
  cam->SetFocalPoint(50, 0, 0);
  cam->SetViewUp(0, 1, 0);
  cam->SetPosition(50, 0, 1);
  cam->Azimuth (vtkMath::DegreesFromRadians(this->CameraThetaAngle));
  cam->Elevation (vtkMath::DegreesFromRadians(this->CameraPhiAngle));
  if (this->CameraDistance == 0.0)
  {
    this->renderer->ResetCamera(this->Bounds);
  } 
  else
  {
    cam->Dolly (1.0 / this->CameraDistance);
  }

  if (ctrl->GetLocalProcessId() == 0) 
  {
    this->window->Render();
    ctrl->TriggerBreakRMIs();
    ctrl->Barrier();
  } 
  else 
  {
    ctrl->ProcessRMIs();
    ctrl->Barrier();
  }

  if (ctrl->GetLocalProcessId() == 0) 
  {
    this->w2i->Modified();
    char *outstring = new char[strlen(this->Filename) + 32];
    sprintf (outstring, "%s%d.png", this->Filename, timestep);
    this->writer->SetFileName(outstring);
    this->writer->Write ();
    printf("Wrote image file %s\n", outstring);
    delete [] outstring;
  }

  return 1;
}
