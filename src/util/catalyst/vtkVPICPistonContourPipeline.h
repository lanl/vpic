#ifndef vtkVPICPistonContourPipeline_H
#define vtkVPICPistonContourPipeline_H

#include "vtkVPICPipeline.h"

class vtkDataSetToPiston;
class vtkTrivialProducer;
class vtkPistonContour;
class vtkPistonMapper;

class vtkVPICPistonContourPipeline : public vtkVPICPipeline
{
public:
  static vtkVPICPistonContourPipeline* New();
  vtkTypeMacro (vtkVPICPistonContourPipeline, vtkVPICPipeline);

  virtual int CoProcess(vtkCPDataDescription* desc);

  // Description:
  // isovalue
  vtkSetMacro (Isovalue, double);
  vtkGetMacro (Isovalue, double);

protected:
  vtkVPICPistonContourPipeline ();
  virtual ~vtkVPICPistonContourPipeline ();

  vtkDataSetToPiston *d2p;
  vtkTrivialProducer* input;
  vtkPistonContour *contour;
  vtkPistonMapper *mapper;

  double Isovalue;
};

#endif /* vtkVPICPistonContourPipeline_H */
