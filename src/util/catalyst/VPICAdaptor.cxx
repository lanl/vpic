#include "VPICAdaptor.h"

#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkMultiProcessController.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

#include <float.h>

#include "vpic.h"
#include "grid.h"
#include "sf_interface.h"
#include "util_base.h"
#include "dumpnames.h"

namespace
{
  vtkCPProcessor* coProcessor = 0;

  // functions to share VPIC fields and hydros between procs
  enum action_t
  {
    PUT = 0, // put local data into remote array for sending
    GET = 1  // get local data from remote array after receiving
  };

  // ix, iy, iz are this proc's topology index in each direction starting from 0
  void RANK_TO_INDEX(int topology[3], int rank, int& ix, int& iy, int& iz)
  {
    ix  = (rank);  /* ix = ix + gpx*( iy + gpy*iz ) */
    iy  = ix/topology[0]; /* iy = iy + gpy*iz */
    ix -= iy*topology[0]; /* ix = ix */
    iz  = iy/topology[1]; /* iz = iz */
    iy -= iz*topology[1]; /* iy = iy */
  }

  int INDEX_TO_RANK(int topology[3], int ix, int iy, int iz)
  {
    /* Wrap processor index periodically */
    while(ix>=topology[0]) ix-=topology[0];
    while(ix<0) ix+=topology[0];
    while(iy>=topology[1]) iy-=topology[1];
    while(iy<0) iy+=topology[1];
    while(iz>=topology[2]) iz-=topology[2];
    while(iz<0) iz+=topology[2];
    /* Compute the rank */
    return (ix + topology[0]*( iy + topology[1]*iz ));
  }

  void FillData(action_t action, vtkPointData* localArrays, std::vector<float>& remoteData,
                std::vector<vtkIdType>& idmap )
  {
    float vals[9];
    size_t counter = 0;
    for(int i=0;i<localArrays->GetNumberOfArrays();i++)
      {
      vtkFloatArray* localArray = vtkFloatArray::SafeDownCast(localArrays->GetArray(i));
      int numberOfComponents = localArray->GetNumberOfComponents();
      for(size_t j=0;j<idmap.size();j++)
        {
        if(action == PUT)
          {
          localArray->GetTupleValue(idmap[j], vals);
          for(int k=0;k<numberOfComponents;k++)
            {
            remoteData[counter+k] = vals[k];
            }
          counter += numberOfComponents;
          }
        else
          {
          for(int k=0;k<numberOfComponents;k++)
            {
            vals[k] = remoteData[counter+k];
            }
          counter += numberOfComponents;
          localArray->SetTupleValue(idmap[j], vals);
          }
        }
      }
  }

  void GetSendIdMap(int ix, int iy, int iz, int rix, int riy, int riz,
                    vtkImageData* imageData, std::vector<vtkIdType>& idmap)
  {
    // assumes that this process has received all information
    // from it before it sends anything out. otherwise this process
    // may not have the proper info to send
    idmap.clear();
    int dimension[3];
    imageData->GetDimensions(dimension);
    if(ix == rix+1 && iy == riy && iz == riz)
      { // -x direction only
      for(int k=0;k<dimension[2];k++)
        {
        for(int j=0;j<dimension[1];j++)
          {
          idmap.push_back(j*dimension[0]+k*dimension[0]*dimension[1]);
          }
        }
      }
    else if(ix == rix && iy == riy+1 && iz == riz)
      { // -y direction only
      for(int k=0;k<dimension[2];k++)
        {
        for(int i=0;i<dimension[0];i++)
          {
          idmap.push_back(i+k*dimension[0]*dimension[1]);
          }
        }
      }
    else if(ix == rix && iy == riy && iz == riz+1)
      { // -z direction only
      for(int j=0;j<dimension[1];j++)
        {
        for(int i=0;i<dimension[0];i++)
          {
          idmap.push_back(i+j*dimension[0]);
          }
        }
      }
  }

  void GetReceiveIdMap(int ix, int iy, int iz, int rix, int riy, int riz,
                       vtkImageData* imageData, std::vector<vtkIdType>& idmap)
  {
    idmap.clear();
    int dimension[3];
    imageData->GetDimensions(dimension);
    if(ix == rix-1 && iy == riy && iz == riz)
      { // +x direction only
      for(int k=0;k<dimension[2];k++)
        {
        for(int j=0;j<dimension[1];j++)
          {
          idmap.push_back(j*dimension[0]+k*dimension[0]*dimension[1]+dimension[0]-1);
          }
        }
      }
    else if(ix == rix && iy == riy-1 && iz == riz)
      { // +y direction only
      for(int k=0;k<dimension[2];k++)
        {
        for(int i=0;i<dimension[0];i++)
          {
          idmap.push_back(i+k*dimension[0]*dimension[1]+dimension[0]*(dimension[1]-1));
          }
        }
      }
    else if(ix == rix && iy == riy && iz == riz-1)
      { // +z direction only
      for(int j=0;j<dimension[1];j++)
        {
        for(int i=0;i<dimension[0];i++)
          {
          idmap.push_back(i+j*dimension[0]+dimension[0]*dimension[1]*(dimension[2]-1));
          }
        }
      }
  }

  void CommunicateData(int topology[3], vtkImageData* imageData)
  {
    vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
    if(controller->GetNumberOfProcesses() == 1)
      {
      return;
      }
    int ix, iy, iz;
    RANK_TO_INDEX(topology, controller->GetLocalProcessId(), ix, iy, iz);

    int numValuesPerNode = 0;
    for(int i=0;i<imageData->GetPointData()->GetNumberOfArrays();i++)
      {
      numValuesPerNode += imageData->GetPointData()->GetArray(i)->GetNumberOfComponents();
      }

    std::vector<vtkIdType> idmap;
    std::vector<float> data;
    // first do receives to make sure data on this process is correct
    if(ix != topology[0]-1)
      {
      int receiveFromRank = INDEX_TO_RANK(topology, ix+1, iy, iz);
      GetReceiveIdMap(ix, iy, iz, ix+1, iy, iz, imageData, idmap);
      data.resize(idmap.size()*numValuesPerNode);
      controller->Receive(&(data[0]), static_cast<vtkIdType>(data.size()), receiveFromRank, 825);
      FillData(GET, imageData->GetPointData(), data, idmap);
      }
    if(iy != topology[1]-1)
      {
      int receiveFromRank = INDEX_TO_RANK(topology, ix, iy+1, iz);
      GetReceiveIdMap(ix, iy, iz, ix, iy+1, iz, imageData, idmap);
      data.resize(idmap.size()*numValuesPerNode);
      controller->Receive(&(data[0]), static_cast<vtkIdType>(data.size()), receiveFromRank, 825);
      FillData(GET, imageData->GetPointData(), data, idmap);
      }
    if(iz != topology[2]-1)
      {
      int receiveFromRank = INDEX_TO_RANK(topology, ix, iy, iz+1);
      GetReceiveIdMap(ix, iy, iz, ix, iy, iz+1, imageData, idmap);
      data.resize(idmap.size()*numValuesPerNode);
      controller->Receive(&(data[0]), static_cast<vtkIdType>(data.size()), receiveFromRank, 825);
      FillData(GET, imageData->GetPointData(), data, idmap);
      }
    // now do sends with correct data
    if(ix != 0)
      {
      int sendToRank = INDEX_TO_RANK(topology, ix-1, iy, iz);
      GetSendIdMap(ix, iy, iz, ix-1, iy, iz, imageData, idmap);
      data.resize(idmap.size()*numValuesPerNode);
      FillData(PUT, imageData->GetPointData(), data, idmap);
      controller->Send(&(data[0]), static_cast<vtkIdType>(data.size()), sendToRank, 825);
      }
    if(iy != 0)
      {
      int sendToRank = INDEX_TO_RANK(topology, ix, iy-1, iz);
      GetSendIdMap(ix, iy, iz, ix, iy-1, iz, imageData, idmap);
      data.resize(idmap.size()*numValuesPerNode);
      FillData(PUT, imageData->GetPointData(), data, idmap);
      controller->Send(&(data[0]), static_cast<vtkIdType>(data.size()), sendToRank, 825);
      }
    if(iz != 0)
      {
      int sendToRank = INDEX_TO_RANK(topology, ix, iy, iz-1);
      GetSendIdMap(ix, iy, iz, ix, iy, iz-1, imageData, idmap);
      data.resize(idmap.size()*numValuesPerNode);
      FillData(PUT, imageData->GetPointData(), data, idmap);
      controller->Send(&(data[0]), static_cast<vtkIdType>(data.size()), sendToRank, 825);
      }
  }

} // end anonymous namespace

void coprocessorinitialize (std::vector<std::string>& pythonScripts)
{
  if (!coProcessor)
    {
    coProcessor = vtkCPProcessor::New();
    coProcessor->Initialize();
    for(std::vector<std::string>::iterator it=pythonScripts.begin();
        it!=pythonScripts.end();it++)
      {
      cout << "adding in catalyst script " << *it << endl;
      vtkNew<vtkCPPythonScriptPipeline> pipeline;
      pipeline->Initialize(it->c_str());
      coProcessor->AddPipeline(pipeline.GetPointer());
      }
    }
}

// topology is number of blocks in each direction
void coprocessorProcess (long long timestep, double time,
                         vpic_simulation* sim, int topology[3],
                         std::vector<DumpParameters *>& dumpParams)
{
  if (!coProcessor)
    {
    cerr << "CoProcessor has not been properly initialized" << endl;
    return;
    }

  if(vtkMultiProcessController::GetGlobalController()->GetLocalProcessId() == 0)
    {
    cout << "Catalyst working on timestep " << timestep << endl;
    }

  vtkSmartPointer<vtkCPDataDescription> coProcessorData =
    vtkSmartPointer<vtkCPDataDescription>::New();
  coProcessorData->AddInput("input");
  coProcessorData->SetTimeData (time, static_cast<vtkIdType>(timestep));
  if (coProcessor->RequestDataDescription (coProcessorData))
    {
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();

    int ix, iy, iz;
    RANK_TO_INDEX(topology, vtkMultiProcessController::GetGlobalController()->GetLocalProcessId(),
                  ix, iy, iz);
    int extent[6] = {sim->grid->nx*ix, sim->grid->nx*(ix+1),
                     sim->grid->ny*iy, sim->grid->ny*(iy+1),
                     sim->grid->nz*iz, sim->grid->nz*(iz+1)};
    if(ix+1 == sim->px)
      {
      extent[1] -= 1;
      }
    if(iy+1 == sim->py)
      {
      extent[3] -= 1;
      }
    if(iz+1 == sim->pz)
      {
      extent[5] -= 1;
      }
    imageData->SetExtent(extent);
    // x0, y0 and z0 are local quantities
    double origin[3] = {sim->grid->x0, sim->grid->y0, sim->grid->z0};
    vtkMultiProcessController::GetGlobalController()->Broadcast(origin, 3, 0);
    imageData->SetOrigin(origin);
    imageData->SetSpacing(sim->grid->dx, sim->grid->dy, sim->grid->dz);


    // Create a variable list of field variables to output.
    size_t numvars = std::min(dumpParams[0]->output_vars.bitsum(
                                field_indeces, total_field_groups),
                              total_field_groups);
    std::vector<size_t> varlist;

    for(size_t v=0; v<total_field_groups; v++)
      {
      if(dumpParams[0]->output_vars.bitset(field_indeces[v]))
        {
        varlist.push_back(v);
        }
      }

    size_t numvarsscalar = std::min(dumpParams[0]->output_vars.bitsum(),
                                    total_field_variables);
    std::vector<size_t> varlistscalar;

    for(size_t i=0; i<total_field_variables; i++)
      {
      if(dumpParams[0]->output_vars.bitset(i))
        {
        varlistscalar.push_back(i);
        }
      }

    double tuple[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}; // ACB -- not type consistent here
    int varlistscalarcounter = 0; // use to keep track of what arrays I already got out of field_array
    for(size_t v=0;v<numvars;v++)
      {
      vtkSmartPointer<vtkDataArray> field;
      if(strcmp(fieldInfo[varlist[v]].type, "FLOATING_POINT") == 0)
        {
        field = vtkSmartPointer<vtkFloatArray>::New();
        }
      else
        {
        cerr << "VPICAdaptor.cxx: unknown field type " << fieldInfo[varlist[v]].type << endl;
        throw 1;
        }
      field->SetName(fieldInfo[varlist[v]].name);
      int numarrayscalars = 1;
      if(strcmp(fieldInfo[varlist[v]].degree, "VECTOR") == 0)
        {
        field->SetNumberOfComponents(3);
        numarrayscalars = 3;
        }
      else if(strcmp(fieldInfo[varlist[v]].degree, "TENSOR") == 0)
        {
        numarrayscalars = atoi(fieldInfo[varlist[v]].elements);
        field->SetNumberOfComponents(9);
        }
      else if(strcmp(fieldInfo[varlist[v]].degree, "SCALAR") != 0)
        {
        cerr << "VPICAdaptor.cxx: unknown field degree " << fieldInfo[varlist[v]].degree << endl;
        throw 1;
        }
      field->SetNumberOfTuples(imageData->GetNumberOfPoints());
      imageData->GetPointData()->AddArray(field);
      vtkIdType tupleCounter = 0;
      for(int k=1;k<extent[5]-extent[4]+2;k++) // dump.cc has a +2 here
        {
        for(int j=1;j<extent[3]-extent[2]+2;j++) // dump.cc has a +2 here
          {
          for(int i=1;i<extent[1]-extent[0]+2;i++) // dump.cc has a +2 here
            {
            for(int l=0;l<numarrayscalars;l++)
              {
              field_t jj = sim->field_array->f[ VOXEL(i, j, k, sim->grid->nx, sim->grid->ny, sim->grid->nz)];
              float* ff = &(jj.ex);
              tuple[l] = ff[varlistscalar[l+varlistscalarcounter]];
              }
            if(numarrayscalars == 6)
              { // symmetric 2nd order tensor. tuple[0] is unchanged
              tuple[8] = tuple[2];
              tuple[7] = tuple[3];
              tuple[2] = tuple[6] = tuple[4];
              tuple[4] = tuple[1];
              tuple[1] = tuple[3] = tuple[5];
              tuple[5] = tuple[7];
              }
            field->SetTuple(tupleCounter, tuple);
            tupleCounter++;
            }
          }
        }
      varlistscalarcounter += numarrayscalars;
      }


    // Process the species...
    for(size_t dumpParam=1;dumpParam<dumpParams.size();dumpParam++) // starting with dumpParam=1 based on dump.cc
      {
      std::string sp_name;
      if(strcmp(dumpParams[dumpParam]->baseFileName, "ehydro") == 0)
        {
        sp_name = "electron";
        }
      else if(strcmp(dumpParams[dumpParam]->baseFileName, "Hhydro") == 0)
        {
        sp_name = "ion";
        }
      else
        {
        cerr << "VPICAdaptor.cxx: don't know how to process species " <<
          dumpParams[dumpParam]->baseFileName << endl;
        throw 1;
        }
      species_t *sp = find_species_name( sp_name.c_str(), sim->species_list );
      if( !sp )
        {
        ERROR(( "Invalid species \"%s\"", sp_name.c_str() ));
        }

      clear_hydro_array( sim->hydro_array );
      accumulate_hydro_p( sim->hydro_array, sp, sim->interpolator_array );
      synchronize_hydro_array( sim->hydro_array );

      varlistscalarcounter = 0; // use to keep track of what arrays I already got out of this hydro species_array
      numvars = std::min(dumpParams[dumpParam]->output_vars.bitsum(
                           hydro_indeces, total_hydro_groups),
                         total_hydro_groups);
      std::string speciesName = "(";
      speciesName.append(dumpParams[dumpParam]->baseFileName);
      speciesName.append(")");
      varlist.clear();
      for(size_t v=0; v<total_hydro_groups; v++)
        {
        if(dumpParams[dumpParam]->output_vars.bitset(hydro_indeces[v]))
          {
          varlist.push_back(v);
          }
        }

      numvarsscalar = std::min(dumpParams[dumpParam]->output_vars.bitsum(),
                               total_hydro_variables);
      varlistscalar.clear();

      for(size_t i=0; i<total_hydro_variables; i++)
        {
        if(dumpParams[dumpParam]->output_vars.bitset(i))
          {
          varlistscalar.push_back(i);
          }
        }

      for(size_t v=0;v<numvars;v++)
        {
        vtkSmartPointer<vtkDataArray> hydro;
        if(strcmp(hydroInfo[varlist[v]].type, "FLOATING_POINT") == 0)
          {
          hydro = vtkSmartPointer<vtkFloatArray>::New();
          }
        else
          {
          cerr << "VPICAdaptor.cxx: unknown hydro type " << hydroInfo[varlist[v]].type << endl;
          throw 1;
          }
        std::string hydroName = hydroInfo[varlist[v]].name+speciesName;
        hydro->SetName(hydroName.c_str());
        int numarrayscalars = 1;
        if(strcmp(hydroInfo[varlist[v]].degree, "VECTOR") == 0)
          {
          numarrayscalars = 3;
          hydro->SetNumberOfComponents(3);
          }
        else if(strcmp(hydroInfo[varlist[v]].degree, "TENSOR") == 0)
          {
          numarrayscalars = atoi(hydroInfo[varlist[v]].elements);
          hydro->SetNumberOfComponents(9);
          }
        else if(strcmp(hydroInfo[varlist[v]].degree, "SCALAR") != 0)
          {
          cerr << "VPICAdaptor.cxx: unknown hydro degree " << hydroInfo[varlist[v]].degree << endl;
          throw 1;
          }
        hydro->SetNumberOfTuples(imageData->GetNumberOfPoints());
        imageData->GetPointData()->AddArray(hydro);
        vtkIdType tupleCounter = 0;
        for(int k=1;k<extent[5]-extent[4]+2;k++) // dump.cc has a +2 here
          {
          for(int j=1;j<extent[3]-extent[2]+2;j++) // dump.cc has a +2 here
            {
            for(int i=1;i<extent[1]-extent[0]+2;i++) // dump.cc has a +2 here
              {
              for(int l=0;l<numarrayscalars;l++)
                {
                hydro_t jj = sim->hydro_array->h[ VOXEL(i, j, k, sim->grid->nx, sim->grid->ny, sim->grid->nz)];
                float* ff = &(jj.jx);
                tuple[l] = ff[varlistscalar[l+varlistscalarcounter]];
                }
              if(numarrayscalars == 6)
                { // symmetric 2nd order tensor. tuple[0] is unchanged
                tuple[8] = tuple[2];
                tuple[7] = tuple[3];
                tuple[2] = tuple[6] = tuple[4];
                tuple[4] = tuple[1];
                tuple[1] = tuple[3] = tuple[5];
                tuple[5] = tuple[7];
                }
              hydro->SetTuple(tupleCounter, tuple);
              tupleCounter++;
              }
            }
          }
        varlistscalarcounter += numarrayscalars;
        } // iterating over vars
      } // iterating over dumpParams


    // fill in the point data from other processes
    CommunicateData(topology, imageData);

    coProcessorData->GetInputDescriptionByName("input")->SetGrid(imageData);
    int wholeExtent[6] = {0, static_cast<int>(sim->grid->nx*sim->px-1),
                          0, static_cast<int>(sim->grid->ny*sim->py-1),
                          0, static_cast<int>(sim->grid->nz*sim->pz-1)};

    coProcessorData->GetInputDescriptionByName("input")->SetWholeExtent(wholeExtent);
    coProcessor->CoProcess(coProcessorData);
    }
}


void coprocessorfinalize ()
{
  std::cout << "Finalizing in vpicadaptor.cxx\n";
  if (coProcessor)
    {
    coProcessor->Finalize();
    coProcessor->Delete();
    coProcessor = 0;
    }
}
