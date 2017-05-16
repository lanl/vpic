#ifndef vtkVPICDataSetPipeline_H
#define vtkVPICDataSetPipeline_H

#include "vtkVPICPipeline.h"

class vtkDataSetMapper;
class vtkActor;

class vtkVPICDataSetPipeline : public vtkVPICPipeline
{
public:
  static vtkVPICDataSetPipeline* New();
  vtkTypeMacro (vtkVPICDataSetPipeline, vtkVPICPipeline);

  virtual int CoProcess(vtkCPDataDescription* desc);

protected:
  vtkVPICDataSetPipeline ();
  virtual ~vtkVPICDataSetPipeline ();

  vtkDataSetMapper* dataSetMapper;
};

#endif /* vtkVPICDataSetPipeline_H */
