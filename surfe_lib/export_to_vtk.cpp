#include <stdio.h>
#include <export_to_vtk.h>
#include <string>
#include <vector>
#include <cstring>

#include <iostream>
using namespace Surfe;
void Surfe::write_to_vtk(Basic_input input, std::string filename) {
    auto points = *input.evaluation_pts;
    FILE* pFile = fopen(filename.c_str(), "w");
    if (pFile!=NULL){
    fprintf(pFile, "%s", "<VTKFile type=\"UnstructuredGrid\"");
    fprintf(pFile, "%s", " version=\"0.1\"");
    fprintf(pFile, "%s\n", " byte_order=\"LittleEndian\">");
    fprintf(pFile, "%s\n", "   <UnstructuredGrid>");
    fprintf(pFile, "%s", "      <Piece  ");
    fprintf(pFile, "%s%i%s", "NumberOfPoints=\"", points.size(), "\"  ");
    fprintf(pFile, "%s%i%s\n", "NumberOfCells=\"", points.size(), "\">");

    // VERTICES
    fprintf(pFile, "%s\n", "         <Points>");
    fprintf(pFile, "%s%s%s\n",
            "            <DataArray type=\"Float32\" NumberOfComponents=\"",
            "3", "\" format=\"ascii\">");
    for (size_t iV = 0; iV < points.size(); iV++) {
        // vertices points (x_{i=0} y_{i=0} z_{i=0} x_{i=1} y_{i=1}
        // z_{i=1}
        //....
        // x_{i=n} y_{i=n} z_{i=n}) i is the node index
        fprintf(pFile, "%f %f %f", points[iV].x(), points[iV].y(),
                points[iV].z());
        fprintf(pFile, "%s\n", "");
    }
    fprintf(pFile, "%s\n", "            </DataArray>");
    fprintf(pFile, "%s\n", "         </Points>");
    //fprintf(pFile, "%s\n", "         <Cells>");
    //fprintf(pFile, "%s\n",
    //        "            <DataArray type=\"Int32\" Name=\"connectivity\" "
    //        "format=\"ascii\">");
    //// Lists all vertices of all nodes
    //for (int iC = 0; iC < points.size(); iC++) {
    //    fprintf(pFile, "%i ", iC);
    //    //fprintf(pFile, "%s\n", "");

    //    //    for (size_t iVl = 0; iVl < c_vertices[iR].size(); iVl++) {
    //    //        fprintf(pFile, "%zu ", c_vertices[iR][iVl]);
    //    //    }
    //}
    //fprintf(pFile, "%s\n", "");
    //fprintf(pFile, "%s\n", "            </DataArray>");

    //fprintf(pFile, "%s\n",
    //        "            <DataArray type=\"Int32\" Name=\"offsets\" "
    //        "format=\"ascii\">");
    //int prev = 0;
    //// probably could be done more efficiently like building fofsets
    //// array when
    //// making cells list but...
    //for (int iC = 0; iC < points.size(); iC++) {
    //    prev += 1;  // vertices.size();  // offsets are the end position of the
    //                // cell in
    //    // the nodes array?

    //    fprintf(pFile, "%i ", prev);
    //}
    //fprintf(pFile, "%s\n", "");

    //fprintf(pFile, "%s\n", "            </DataArray>");
    //fprintf(pFile, "%s\n",
    //        "            <DataArray type=\"Int32\" Name=\"types\" "
    //        "format=\"ascii\">");
    //// probably could be done more efficiently like building fofsets
    //// array when


    //// making cells list but...
    //for (int iC = 0; iC < points.size(); iC++) {
    //    fprintf(pFile, "%i ", 1);
    //}
    //fprintf(pFile, "%s\n", "");
    //fprintf(pFile, "%s\n", "            </DataArray>");
    //fprintf(pFile, "%s\n", "            </Cells>");
    //fprintf(pFile, "%s\n", "         <CellData>");
    //fprintf(pFile, "%s%s%s\n",
    //        "            <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //        "Name=\"",
    //        "prop",
    //        "\" "
    //        "format=\"ascii\">");
    //for (int i = 0; i < points.size(); i++) {
    //    fprintf(pFile, "%f \n",points[i].scalar_field());
    //                  

    //}
    //fprintf(pFile, "%s\n", "");

    //fprintf(pFile, "%s\n", "            </DataArray>");

    //fprintf(pFile, "%s\n", "         </CellData>");
    fprintf(pFile, "%s\n", "      </Piece>  ");
    fprintf(pFile, "%s\n", "   </UnstructuredGrid>");
    fprintf(pFile, "%s\n", "</VTKFile>");
    std::cout<<"Finished"<<std::endl;
}
}
