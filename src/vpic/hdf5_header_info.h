#ifndef VPIC_HDF5_HEAD_INFO
#define VPIC_HDF5_HEAD_INFO

#define FIELD_ARRAY_NAME field_array

// XML header stuff
const char *header = "<?xml version=\"1.0\"?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n\t<Domain>\n";
const char *header_topology = "\t\t<Topology Dimensions=\"%s\" TopologyType=\"3DCoRectMesh\" name=\"topo\"/>\n";
const char *header_geom = "\t\t<Geometry Type=\"ORIGIN_DXDYDZ\" name=\"geo\">\n";
const char *header_origin = "\t\t\t<!-- Origin --> \n\t\t\t<DataItem Dimensions=\"3\" Format=\"XML\">%s</DataItem>\n";
const char *header_dxdydz = "\t\t\t<!-- DxDyDz --> \n\t\t\t<DataItem Dimensions=\"3\" Format=\"XML\">%s</DataItem>\n";
const char *footer_geom = "\t\t</Geometry>\n";
const char *grid_line = "\t\t<Grid CollectionType=\"Temporal\" GridType=\"Collection\" Name=\"TimeSeries\"> \n \
\t\t\t<Time TimeType=\"HyperSlab\"> \n \
\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"XML\" NumberType=\"Float\">";
const char *grid_line_footer = "</DataItem> \n\
\t\t\t</Time>\n";
const char *footer = "\t\t</Grid>\n\t</Domain>\n</Xdmf>\n";

const char *main_body_head = "\t\t\t<Grid GridType=\"Uniform\" Name=\"T%d\"> \n \
\t\t\t\t<Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>   \n \
\t\t\t\t<Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>  \n";
const char *main_body_foot = "\t\t\t</Grid>\n";

const char *main_body_attributeV = "\
        \t\t\t\t <Attribute AttributeType =\"Vector\" Center=\"Node\" Name=\"%s\">  \n \
            \t\t\t\t\t<DataItem Dimensions=\" %s \" Function=\"JOIN($0, $1, $2)\" ItemType=\"Function\">  \n \
                \t\t\t\t\t\t<DataItem ItemType=\"Uniform\" Dimensions=\" %s \" DataType=\"Float\" Precision=\"4\" Format=\"HDF\"> T.%d/%s_%d.h5:/Timestep_%d/%s </DataItem>  \n \
                \t\t\t\t\t\t<DataItem ItemType=\"Uniform\" Dimensions=\" %s \" DataType=\"Float\" Precision=\"4\" Format=\"HDF\"> T.%d/%s_%d.h5:/Timestep_%d/%s </DataItem>  \n \
                \t\t\t\t\t\t<DataItem ItemType=\"Uniform\" Dimensions=\" %s \" DataType=\"Float\" Precision=\"4\" Format=\"HDF\"> T.%d/%s_%d.h5:/Timestep_%d/%s </DataItem>  \n \
            \t\t\t\t\t</DataItem>  \n \
        \t\t\t\t</Attribute>  \n ";

const char *main_body_attributeS = "\
        \t\t\t\t <Attribute AttributeType =\"Scalar\" Center=\"Node\" Name=\"%s\">  \n \
                \t\t\t\t\t\t<DataItem ItemType=\"Uniform\" Dimensions=\" %s \" DataType=\"Float\" Precision=\"4\" Format=\"HDF\"> T.%d/%s_%d.h5:/Timestep_%d/%s </DataItem>  \n \
        \t\t\t\t</Attribute>  \n ";

#define create_file_with_header(xml_file_name, dimensions, orignal, dxdydz, nframes, fields_interval) \
  {                                                                                                   \
    FILE *fp;                                                                                         \
    fp = fopen(xml_file_name, "w");                                                                   \
    fputs(header, fp);                                                                                \
    fprintf(fp, header_topology, dimensions);                                                         \
    fputs(header_geom, fp);                                                                           \
    fprintf(fp, header_origin, orignal);                                                              \
    fprintf(fp, header_dxdydz, dxdydz);                                                               \
    fputs(footer_geom, fp);                                                                           \
    fprintf(fp, grid_line, nframes);                                                                  \
    int i;                                                                                            \
    for (i = 0; i < nframes; i++)                                                                     \
      fprintf(fp, "%d ", i*fields_interval);                                                         \
    fputs(grid_line_footer, fp);                                                                      \
    fclose(fp);                                                                                       \
  }
#define write_main_body_attribute(fpp, main_body_attribute_p, attribute_name, dims_4d_p, dims_3d_p, file_name_pre_p, time_step_p, a1, a2, a3) \
  {                                                                                                                                           \
    fprintf(fpp, main_body_attribute_p, attribute_name, dims_4d_p,                                                                            \
            dims_3d_p, time_step_p, file_name_pre_p, time_step_p, time_step_p, a1,                                                            \
            dims_3d_p, time_step_p, file_name_pre_p, time_step_p, time_step_p, a2,                                                            \
            dims_3d_p, time_step_p, file_name_pre_p, time_step_p, time_step_p, a3);                                                           \
  }

#define invert_field_xml_item(xml_file_name, speciesname_p, time_step, dims_4d, dims_3d, add_footer_flag)                                 \
  {                                                                                                                                       \
    FILE *fp;                                                                                                                             \
    fp = fopen(xml_file_name, "a");                                                                                                       \
    fprintf(fp, main_body_head, time_step);                                                                                               \
    if (field_dump_flag.enabledE())                                                                                                       \
      write_main_body_attribute(fp, main_body_attributeV, "E", dims_4d, dims_3d, speciesname_p, time_step, "ex", "ey", "ez");             \
    if (field_dump_flag.div_e_err)                                                                                                        \
      fprintf(fp, main_body_attributeS, "div_e_err", dims_3d, time_step, speciesname_p, time_step, time_step, "div_e_err");               \
    if (field_dump_flag.enabledCB())                                                                                                      \
      write_main_body_attribute(fp, main_body_attributeV, "B", dims_4d, dims_3d, speciesname_p, time_step, "cbx", "cby", "cbz");          \
    if (field_dump_flag.div_b_err)                                                                                                        \
      fprintf(fp, main_body_attributeS, "div_b_err", dims_3d, time_step, speciesname_p, time_step, time_step, "div_b_err");               \
    if (field_dump_flag.enabledTCA())                                                                                                     \
      write_main_body_attribute(fp, main_body_attributeV, "TCA", dims_4d, dims_3d, speciesname_p, time_step, "tcax", "tcay", "tcaz");     \
    if (field_dump_flag.rhob)                                                                                                             \
      fprintf(fp, main_body_attributeS, "rhob", dims_3d, time_step, speciesname_p, time_step, time_step, "rhob");                         \
    if (field_dump_flag.enabledJF())                                                                                                      \
      write_main_body_attribute(fp, main_body_attributeV, "JF", dims_4d, dims_3d, speciesname_p, time_step, "jfx", "jfy", "jfz");         \
    if (field_dump_flag.rhof)                                                                                                             \
      fprintf(fp, main_body_attributeS, "rhof", dims_3d, time_step, speciesname_p, time_step, time_step, "rhof");                         \
    if (field_dump_flag.enabledEMAT())                                                                                                    \
      write_main_body_attribute(fp, main_body_attributeV, "EMAT", dims_4d, dims_3d, speciesname_p, time_step, "ematx", "ematy", "ematz"); \
    if (field_dump_flag.nmat)                                                                                                             \
      fprintf(fp, main_body_attributeS, "nmat", dims_3d, time_step, speciesname_p, time_step, time_step, "nmat");                         \
    if (field_dump_flag.enabledFMAT())                                                                                                    \
      write_main_body_attribute(fp, main_body_attributeV, "FMAT", dims_4d, dims_3d, speciesname_p, time_step, "fmatx", "fmaty", "fmatz"); \
    if (field_dump_flag.cmat)                                                                                                             \
      fprintf(fp, main_body_attributeS, "cmat", dims_3d, time_step, speciesname_p, time_step, time_step, "cmat");                         \
    fprintf(fp, "%s", main_body_foot);                                                                                                          \
    if (add_footer_flag)                                                                                                                  \
      fputs(footer, fp);                                                                                                                  \
    fclose(fp);                                                                                                                           \
  }
#define invert_hydro_xml_item(xml_file_name, speciesname_p, time_step, dims_4d, dims_3d, add_footer_flag)                          \
  {                                                                                                                                \
    FILE *fp;                                                                                                                      \
    fp = fopen(xml_file_name, "a");                                                                                                \
    fprintf(fp, main_body_head, time_step);                                                                                        \
    if (hydro_dump_flag.enabledJ())                                                                                                \
      write_main_body_attribute(fp, main_body_attributeV, "J", dims_4d, dims_3d, speciesname_p, time_step, "jx", "jy", "jz");      \
    if (hydro_dump_flag.rho)                                                                                                       \
      fprintf(fp, main_body_attributeS, "rho", dims_3d, time_step, speciesname_p, time_step, time_step, "rho");                    \
    if (hydro_dump_flag.enabledP())                                                                                                \
      write_main_body_attribute(fp, main_body_attributeV, "P", dims_4d, dims_3d, speciesname_p, time_step, "px", "py", "pz");      \
    if (hydro_dump_flag.ke)                                                                                                        \
      fprintf(fp, main_body_attributeS, "ke", dims_3d, time_step, speciesname_p, time_step, time_step, "ke");                      \
    if (hydro_dump_flag.enabledTD())                                                                                               \
      write_main_body_attribute(fp, main_body_attributeV, "TD", dims_4d, dims_3d, speciesname_p, time_step, "txx", "tyy", "tzz");  \
    if (hydro_dump_flag.enabledTOD())                                                                                              \
      write_main_body_attribute(fp, main_body_attributeV, "TOD", dims_4d, dims_3d, speciesname_p, time_step, "tyz", "tzx", "txy"); \
    fprintf(fp, "%s", main_body_foot);                                                                                                   \
    if (add_footer_flag)                                                                                                           \
      fputs(footer, fp);                                                                                                           \
    fclose(fp);                                                                                                                    \
  }


#endif // VPIC_HDF5_HEAD_INFO
