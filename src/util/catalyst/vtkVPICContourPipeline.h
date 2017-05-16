#ifndef vtkVPICContourPipeline_H
#define vtkVPICContourPipeline_H

#include "vtkVPICPipeline.h"

class vtkContourFilter;
class vtkPolyDataMapper;

class vtkVPICContourPipeline : public vtkVPICPipeline
{
public:
  static vtkVPICContourPipeline* New();
  vtkTypeMacro (vtkVPICContourPipeline, vtkVPICPipeline);

  virtual int CoProcess(vtkCPDataDescription* desc);

  // Description:
  // isovalue
  vtkSetMacro (Isovalue, double);
  vtkGetMacro (Isovalue, double);

protected:
  vtkVPICContourPipeline ();
  virtual ~vtkVPICContourPipeline ();

  vtkContourFilter* contours;
  vtkPolyDataMapper* polyDataMapper;

  double Isovalue;
};

#endif /* vtkVPICContourPipeline_H */
