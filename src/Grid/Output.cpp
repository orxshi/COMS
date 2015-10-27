#include "Grid.h"

void Grid::log()
{
    ofstream out;
    out.open (logDir);
    
    out << "mesh file = " << meshFile << endl;
    out << "number of cells = " << n_in_elm << endl;
    out << "number of bouElms = " << n_bou_elm << endl;
    out << "number of bouElms = " << n_bou_elm << endl;
    out << "total number of elms = " << cell.size() << endl;
    out << "id = " << id << endl;
    out << "nFaces = " << face.size() << endl;
    out << "nPt = " << pt.size() << endl;    
    for (int i=0; i<phys_count; ++i)
    {
        out << "phys_" << phys[i] << " : " << bcVerbose[bc[i]] << endl;
    }
            
    out.close();
}

void Grid::outAllTecplot ()
{
    int len;
    int maxLen = 20;
    ofstream out;
    string dir = outputDir;
    
    string temps = "allTecplot.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);
    
    out.open (dir);

    out << "VARIABLES = ";
    out << "\"X\", "; // 1
    out << "\"Y\", "; // 2
    out << "\"Z\", "; // 3
    out << "\"RHO\", "; // 4
    out << "\"U\", "; // 5
    out << "\"V\", "; // 6
    out << "\"W\", "; // 7
    out << "\"P\", "; // 8
    out << "\"iBlank\", "; // 9
    out << "\"trim\", "; // 10
    out << endl;

    out << "ZONE ";
    out << "ZONETYPE=FEBRICK ";
    out << "N=";
    out << pt.size();
    out << " E=";
    out << n_in_elm;
    out << " DATAPACKING=BLOCK ";
    out << "VARLOCATION=([4,5,6,7,8,9,10]=CELLCENTERED)";
    out << endl;

    // xyz coordinates
    for (int d=0; d<3; ++d)
    {
        len = 0;

        for (const Point& p: pt)
        {
            ++len;
            if (len == maxLen)
            {
                out << endl;
                len = 0;
            }

            out << p.dim[d] << ", ";
        }    

        out << endl;
    }
    
    // prim variables
    for (int n=0; n<N_VAR; ++n)
    {
        len = 0;
    
        for (int e=n_bou_elm; e<cell.size(); ++e)
        {
            ++len;
            if (len == maxLen)
            {
                out << endl;
                len = 0;
            }

            out << cell[e].prim[n] << ", ";
        }
        
        out << endl;
    }
    
    // iBlank    
    len = 0;
    for (int e=n_bou_elm; e<cell.size(); ++e)
    {
        ++len;
        if (len == maxLen)
        {
            out << endl;
            len = 0;
        }

        out << static_cast<int> (cell[e].iBlank) << ", ";
    }

    out << endl;
    
    // trim
    len = 0;
    for (int e=n_bou_elm; e<cell.size(); ++e)
    {
        ++len;
        if (len == maxLen)
        {
            out << endl;
            len = 0;
        }

        out << static_cast<int> (cell[e].trim) << ", ";
    }

    out << endl;

    // connectivity
    for (int ie=n_bou_elm; ie<cell.size(); ++ie)
    {
        const Cell& e = cell[ie];
        
        for (unsigned int v=0; v<e.vtx.size(); ++v)
        {
            if (e.vtx.size()==6 && (v==2 || v==5))
            {
                out << setw(10) << e.vtx[v] + 1;
                //out << setw(10) << e.vtxIndex[v] + 1;
            }
            else if (e.vtx.size()==4)
            {
                if (v==2)
                {
                    out << setw(10) << e.vtx[v] + 1;
                    //out << setw(10) << e.vtxIndex[v] + 1;
                }
                if (v==3)
                {
                    out << setw(10) << e.vtx[v] + 1;
                    out << setw(10) << e.vtx[v] + 1;
                    out << setw(10) << e.vtx[v] + 1;
                    
                    //out << setw(10) << e.vtxIndex[v] + 1;
                    //out << setw(10) << e.vtxIndex[v] + 1;
                    //out << setw(10) << e.vtxIndex[v] + 1;
                }
            }

            out << setw(10) << e.vtx[v] + 1;
            //out << setw(10) << e.vtxIndex[v] + 1;
        }

        out << endl;
    }
    
    out.close();
}

void Grid::outAllUnsteady (string title)
{
    int len;
    int maxLen = 20;
    bool firstTime;    
    
    string dir = outputDir;
    string temps = "unsteady.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);
    
    struct stat buf;
    if (stat(dir.c_str(), &buf) != -1)
    {        
        firstTime = false;
    }
    else
    {        
        firstTime = true;
    }
    
    ofstream out;
    out.open (dir, std::fstream::app);

    out << "VARIABLES = ";
    out << "\"X\", ";
    out << "\"Y\", ";
    out << "\"Z\", ";
    out << "\"RHO\", ";
    out << "\"U\", ";
    out << "\"V\", ";
    out << "\"W\", ";
    out << "\"P\", ";
    out << endl;

    out << "ZONE ";
    out << "T=\"";
    out << title; // alpha
    out << "\" ";
    out << " SOLUTIONTIME=";
    out << time;
    out << " ZONETYPE=FEBRICK ";
    out << "N=";
    out << pt.size();
    out << " E=";
    out << n_in_elm;
    out << " DATAPACKING=BLOCK ";
    out << "VARLOCATION=([4,5,6,7,8]=CELLCENTERED) ";

    if (!firstTime)
    {
        //outUnsteady << "VARSHARELIST=([1,2,3]=1) ";
        out << "CONNECTIVITYSHAREZONE=1";
    }

    out << endl;
    
    // xyz coordinates
    for (int d=0; d<N_DIM; ++d)
    {
        len = 0;

        for (const Point& p: pt)
        {
            ++len;
            if (len == maxLen)
            {
                out << endl;
                len = 0;
            }

            out << p.dim[d] << ", ";
        }    

        out << endl;
    }

    // prim variables
    for (int n=0; n<N_VAR; ++n)
    {
        len = 0;
    
        for (int e=n_bou_elm; e<cell.size(); ++e)
        {
            ++len;
            if (len == maxLen)
            {
                out << endl;
                len = 0;
            }

            out << cell[e].prim[n] << ", ";
        }
        
        out << endl;
    }

    if (firstTime)
    {
        // connectivity
        for (int ie=n_bou_elm; ie<cell.size(); ++ie)
        {
            const Cell& e = cell[ie];

            for (unsigned int v=0; v<e.vtx.size(); ++v)
            {
                if (e.vtx.size()==6 && (v==2 || v==5))
                {
                    out << setw(10) << e.vtx[v] + 1;
                }
                else if (e.vtx.size()==4)
                {
                    if (v==2)
                    {
                        out << setw(10) << e.vtx[v] + 1;
                    }
                    if (v==3)
                    {
                        out << setw(10) << e.vtx[v] + 1;
                        out << setw(10) << e.vtx[v] + 1;
                        out << setw(10) << e.vtx[v] + 1;
                    }
                }

                out << setw(10) << e.vtx[v] + 1;
            }

            out << endl;
        }
    }
    
    out.close();
}

/*void Grid::createOutputDir()
{
    string ids = to_string (id);    
    outputDir = mainDir;
    outputDir.append ("/Grid_");
    outputDir.append (ids);
    outputDir.append ("/");
    
    mkdir (outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}*/


void Grid::outAllVTK (int time)
{
    int len;
    int maxLen = 20;
    int cellListSize = 0;
    ofstream out;
    string dir = outputDir;
    
    string temps = "allVTK_";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);
    
    stringstream ss;
    ss << fixed << time;
    temps = ss.str();
    dir.append (temps);
    temps = ".vtk";
    dir.append (temps);
    
    out.open (dir);

    out << "# vtk DataFile Version 3.0" << endl;
    out << "All in VTK format" << endl;
    out << "ASCII" << endl;
    out << "DATASET UNSTRUCTURED_GRID" << endl;
    out << "POINTS " << pt.size() << " float" << endl;
    
    len = 0;
    
    for (int i=0; i<pt.size(); ++i)
    {
        out << pt[i].dim[0];
        out << " ";
        out << pt[i].dim[1];
        out << " ";
        out << pt[i].dim[2];
        out << endl;
    }
    
    // get cell list size
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        cellListSize += (cell[c].vtx.size() + 1);
    }
    
    out << endl;    
    out << "CELLS " << n_in_elm << " " << cellListSize << endl;
    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        const Cell& cll = cell[c];
        
        out << cll.vtx.size();
        out << " ";
        
        for (int i=0; i<cll.vtx.size(); ++i)
        {
            out << cll.vtx[i];
            out << " ";
        }
        
        out << endl;
    }    
    
    out << "CELL_TYPES " << n_in_elm << endl;
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        const Cell& cll = cell[c];
        
             if (cll.type == elmType_t::TET)
        {
            out << static_cast<int>(vtkCellType_t::TET);
        }
        else if (cll.type == elmType_t::HEX)
        {
            out << static_cast<int>(vtkCellType_t::HEX);
        }
        else if (cll.type == elmType_t::PEN)
        {
            out << static_cast<int>(vtkCellType_t::WEDGE);
        }
        else
        {
            cout << "Error: neither of the specified types in outAllVTK()" << endl;
            cout << "c = " << c << endl;
            exit(-2);
        }
        
        out << endl;
    }
    
    out << "CELL_DATA " << n_in_elm << endl;
    
    out << "SCALARS " << "rho " << "float " << "1" << endl;
    out << "LOOKUP_TABLE default" << endl;    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        const Cell& cll = cell[c];
                     
        out << cll.prim[0] << endl;
    }
    
    out << "SCALARS " << "u " << "float " << "1" << endl;
    out << "LOOKUP_TABLE default" << endl;    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        const Cell& cll = cell[c];
                     
        out << cll.prim[1] << endl;
    }
    
    out << "SCALARS " << "v " << "float " << "1" << endl;
    out << "LOOKUP_TABLE default" << endl;    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        const Cell& cll = cell[c];
                     
        out << cll.prim[2] << endl;
    }
    
    out << "SCALARS " << "w " << "float " << "1" << endl;
    out << "LOOKUP_TABLE default" << endl;    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        const Cell& cll = cell[c];
                     
        out << cll.prim[3] << endl;
    }
    
    out << "SCALARS " << "p " << "float " << "1" << endl;
    out << "LOOKUP_TABLE default" << endl;    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        const Cell& cll = cell[c];
                     
        out << cll.prim[4] << endl;
    }
    
    out << "SCALARS " << "iBlank " << "int " << "1" << endl;
    out << "LOOKUP_TABLE default" << endl;    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        const Cell& cll = cell[c];
                     
        out << static_cast<int>(cll.iBlank) << endl;
    }
    
    out << "SCALARS " << "trim " << "int " << "1" << endl;
    out << "LOOKUP_TABLE default" << endl;    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        const Cell& cll = cell[c];
                     
        out << cll.trim << endl;
    }
    
    out.close();
}