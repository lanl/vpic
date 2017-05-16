#ifndef vtkVPICPipeline_H
#define vtkVPICPipeline_H

#include "vtkCPPipeline.h"

class vtkActor;
class vtkLightKit;
class vtkIceTSynchronizedRenderers;
class vtkSynchronizedRenderWindows;
class vtkRenderWindow;
class vtkRenderer;
class vtkWindowToImageFilter;
class vtkPNGWriter;

class vtkVPICPipeline : public vtkCPPipeline
{
public:
  static vtkVPICPipeline* New();
  vtkTypeMacro (vtkVPICPipeline, vtkCPPipeline);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual int RequestDataDescription(vtkCPDataDescription *desc);  
  virtual int CoProcess(vtkCPDataDescription* desc);

  // Description:
  // name of the image file to output
  vtkSetStringMacro(Filename);
  vtkGetStringMacro(Filename);

  // Description:
  // angle of the camera around the mass of particles in the x-z plane
  vtkSetMacro (CameraThetaAngle, double);
  vtkGetMacro (CameraThetaAngle, double);

  // Description:
  // angle the camera makes with the y axis (relative to particles' center)
  vtkSetMacro (CameraPhiAngle, double);
  vtkGetMacro (CameraPhiAngle, double);

  // Description:
  // distance the camera is from the particles' center
  // zero means to find the optimal viewing distance (default: zero)
  vtkSetMacro (CameraDistance, double);
  vtkGetMacro (CameraDistance, double);

  // Description:
  // bounds of the particle space.  This will be used to set the camera
  // distance (if requested) and the box outline.
  vtkSetVector6Macro (Bounds, double);
  vtkGetVectorMacro (Bounds, double, 6);

  // Description:
  // max and min value for the attributes to use in defining the color lookup
  vtkSetMacro (AttributeMaximum, double);
  vtkGetMacro (AttributeMaximum, double);
  vtkSetMacro (AttributeMinimum, double);
  vtkGetMacro (AttributeMinimum, double);

protected:
  vtkVPICPipeline ();
  virtual ~vtkVPICPipeline ();

  char *Filename;

  double CameraThetaAngle;
  double CameraPhiAngle;
  double CameraDistance;

  double Bounds[6];

  double AttributeMaximum;
  double AttributeMinimum;

  vtkActor* actor;
  vtkLightKit* lightKit;
  vtkIceTSynchronizedRenderers* syncRen;
  vtkSynchronizedRenderWindows* syncWin;
  vtkRenderWindow* window;
  vtkRenderer* renderer;
  vtkWindowToImageFilter* w2i;
  vtkPNGWriter* writer;
};

#endif /* vtkVPICPipeline_H */
