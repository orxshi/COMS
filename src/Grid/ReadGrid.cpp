#include "Grid.h"

void Grid::read_ptSize ()
{
    string temps;
    int tempi;
        
    in.open (meshFile);

    if (in.is_open())
    {
        do
        {
            getline (in, temps);
        }
        while (temps != "$Nodes");

        in >> tempi;
        pt.resize (tempi);
        //ptIndex.resize (tempi);
    }
    else
    {
        cout << "Could not open mesh file" << endl;
    }
}

void Grid::read_pt()
{
    string temps;
    
    for (unsigned int i=0; i<pt.size(); ++i)
    {
        in >> temps;
        in >> pt[i].dim[0];
        in >> pt[i].dim[1];
        in >> pt[i].dim[2];
        
        pt[i].belonging = id;
    }
}

void Grid::read_elmSize ()
{
    std::string temps;
    int tempi;

    do
    {
        getline (in, temps);
    }
    while (temps != "$Elements");

    in >> totalNElms;
}

void Grid::read_elm ()
{
    bool stillBoundary = true;
    int n_tags;
    string temps;
    int tempi;
    int tag_count;
    int n_part_belongs;

    bt.resize (pt.size() + 1);

    for (int e=0; e<totalNElms; ++e)
    {
        in >> temps; // ID
        in >> tempi; // type
        
        Cell tmpElm;
        tmpElm.belonging = id;

        switch (tempi)
        {
            case 2:
                tmpElm.type = elmType_t::TRI;
                break;
            case 3:
                tmpElm.type = elmType_t::QUAD;
                break;
            case 4:
                tmpElm.type = elmType_t::TET;
                break;
            case 5:
                tmpElm.type = elmType_t::HEX;
                break;
            case 6:
                tmpElm.type = elmType_t::PEN;
                break;
            default:
                cout << "unknown boundary element type in read_elm()" << endl;
                exit(-2);
                break;
        }

        if (stillBoundary)
        {
            if (tmpElm.type > elmType_t::QUAD)
            {
                n_bou_elm = e;
                stillBoundary = false;
            }
        }

        in >> n_tags; // num of tags

        tag_count = 0;

        while (n_tags > 0)
        {
            in >> tmpElm.phys;
            ++ tag_count;
            if (tag_count == n_tags) break;
            in >> temps; // geometrical num
            ++ tag_count;
            if (tag_count == n_tags) break;
            in >> n_part_belongs; // num of partitions to which element belongs

            for (int i=0; i<n_part_belongs; ++i)
            {
                in >> temps;
            }

            break;
        }

        switch ( tmpElm.type )
        {
            case elmType_t::TRI:
                //tmpElm.vtx.reserve ( static_cast<int>(nVertices_t::TRI) );
                //tmpElm.face.reserve (static_cast<int>(nFaces_t::TRI));
                tmpElm.nVertices = static_cast<int>(nVertices_t::TRI);
                tmpElm.nFaces = static_cast<int>(nFaces_t::TRI);
                break;

            case elmType_t::HEX:
                //tmpElm.vtx.reserve ( static_cast<int>(nVertices_t::HEX));
                //tmpElm.face.reserve (static_cast<int>(nFaces_t::HEX));
                tmpElm.nVertices = static_cast<int>(nVertices_t::HEX);
                tmpElm.nFaces = static_cast<int>(nFaces_t::HEX);
                break;

            case elmType_t::TET:
                //tmpElm.vtx.reserve ( static_cast<int>(nVertices_t::TET));
                //tmpElm.face.reserve (static_cast<int>(nFaces_t::TET));
                tmpElm.nVertices = static_cast<int>(nVertices_t::TET);
                tmpElm.nFaces = static_cast<int>(nFaces_t::TET);
                break;

            case elmType_t::QUAD:
                //tmpElm.vtx.reserve ( static_cast<int>(nVertices_t::QUAD));
                //tmpElm.face.reserve (static_cast<int>(nFaces_t::QUAD));
                tmpElm.nVertices = static_cast<int>(nVertices_t::QUAD);
                tmpElm.nFaces = static_cast<int>(nFaces_t::QUAD);
                break;

            case elmType_t::PEN:
                //tmpElm.vtx.reserve ( static_cast<int>(nVertices_t::PEN));
                //tmpElm.face.reserve (static_cast<int>(nFaces_t::PEN));
                tmpElm.nVertices = static_cast<int>(nVertices_t::PEN);
                tmpElm.nFaces = static_cast<int>(nFaces_t::PEN);
                break;
            default:
                cout << "elm[" << e << "].type = " << static_cast<int>(tmpElm.type) << endl;
                exit(-2);
                break;
        }

        tmpElm.vtx.reserve (tmpElm.nVertices);
        tmpElm.vtxBelo.reserve (tmpElm.nVertices);
        for (unsigned int i=0; i<tmpElm.nVertices; ++i)
        {
            in >> tempi;
            tmpElm.vtx.push_back(tempi-1);
            tmpElm.vtxBelo.push_back (id);
            //tmpElm.vtx.push_back ( ref(pt[tempi-1]) );
            //tmpElm.vtxIndex.push_back (tempi-1);
            bt[tempi-1].insert(e);
        }
        
        cell.push_back( move(tmpElm) );
        //cell.push_back(tmpElm);
    }
    
    for (int c=0; c<n_bou_elm; ++c)
    {
        cell[c].ghost = true;
    }
    
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        cell[c].ghost = false;
    }

    n_in_elm = totalNElms - n_bou_elm;
}

/*void Grid::read_input()
{
    int tmp;
    string tmps;
    
    ifstream in;
    string inputDir;
    string ids = to_string (id);
    inputDir = "../input_";
    inputDir.append (ids);
    inputDir.append (".dat");
    in.open(inputDir);
    
    in >> tmps; in >> meshFile;
    in >> tmps; in >> steady;
    in >> tmps; in >> cfl;
    in >> tmps; in >> useCFL;
    in >> tmps; in >> dt;
    in >> tmps; in >> tmp;
    
    if (tmp == 1)
    {
        tOrder = tsOrder_t::FIRST;
    }
    else if (tmp == 2)
    {
        tOrder = tsOrder_t::SECOND;
    }
    else
    {
        cout << "temporal order is given as: " << tmp << endl;
        exit(-2);
    }
    
    in >> tmps; in >> tmp;
    
    if (tmp == 1)
    {
        sOrder = tsOrder_t::FIRST;
    }
    else if (tmp == 2)
    {
        sOrder = tsOrder_t::SECOND;
    }
    else
    {
        cout << "spatial order is given as: " << tmp << endl;
        exit(-2);
    }
    
    in >> tmps; in >> implicit;
    in >> tmps; in >> nGaussIter;
    in >> tmps; in >> tol;
    in >> tmps; in >> maxTimeStep;
    in >> tmps; in >> nSolidBoundaries;
    
    for (int i=0; i<nSolidBoundaries; ++i)
    {
        in >> tmpd;
        solidBoundaryPhys.push_back (tmpd);
    }
    
    in >> tmps; in >> solidBoundaryPhys;
    in >> tmps; in >> phys_count;
    
    phys.resize(phys_count);
    bc.resize(phys_count);
    
    for (int i=0; i<phys_count; ++i)
    {
        in >> tmps; in >> phys[i];
    }
    
    for (int i=0; i<phys_count; ++i)
    {
        in >> tmps; in >> bc[i];
    }
    
    in >> tmps; in >> rhoInf;
    in >> tmps; in >> pInf;
    in >> tmps; in >> MachAirfoil;
    in >> tmps; in >> MachAir;
    in >> tmps; in >> aoa;
    in >> tmps; in >> alphaMean;
    in >> tmps; in >> alphaMax;
    in >> tmps; in >> kc;
    
    for (int d=0; d<N_DIM; ++d)
    {
        in >> tmps; in >> centerAirfoil[d];
    }
    
    in.close();
}*/

/*void Grid::printInput()
{
    cout << ">>> Grid: " << id << endl;
    cout << ">>> Mesh file: " << meshFile << endl;
    cout << ">>> CFL: " << cfl << endl;

    if (steady)
    {
        cout << ">>> Steady: true" << endl;
    }
    else
    {
        cout << ">>> Steady: false" << endl;
    }

    if (useCFL)
    {
        cout << ">>> Use CFL: true" << endl;
    }
    else
    {
        cout << ">>> Use CFL: false" << endl;
        cout << ">>> time step: " << dt << endl;
    }

    cout << ">>> Temporal order: " << static_cast<int>(tOrder) << endl;
    cout << ">>> Spatial order: " << static_cast<int>(sOrder) << endl;

    if (implicit)
    {
        cout << ">>> Implicit: true" << endl;
        cout << ">>> nGaussIter: " << nGaussIter << endl;
    }
    else
    {
        cout << ">>> Implicit: false" << endl;
    }

    cout << ">>> Tolerance: " << tol << endl;
    cout << ">>> maxTimeStep: " << maxTimeStep << endl;
    cout << ">>> solidBoundaryPhys: " << solidBoundaryPhys << endl;
    cout << ">>> Number of phys: " << phys_count << endl;
    
    for (int i=0; i<phys_count; ++i)
    {
        cout << ">>> phys_" << phys[i] << " : " << bcVerbose[bc[i]] << endl;
    }
    
    cout << ">>> Ref density: " << rhoInf << endl;
    cout << ">>> Ref pressure: " << pInf << endl;
    cout << ">>> MachAirfoil: " << MachAirfoil << endl;
    cout << ">>> MachAir: " << MachAir << endl;
    cout << ">>> aoa: " << aoa << endl;
    cout << ">>> alphaMean: " << alphaMean << endl;
    cout << ">>> alphaMax: " << alphaMax << endl;
    cout << ">>> kc: " << kc << endl;
    
    for (int d=0; d<N_DIM; ++d)
    {
        cout << ">>> centerAirfoil[" << d << "] : " << centerAirfoil[d] << endl;
    }
}*/

void Grid::printMeshInfo()
{
    cout << ">>> Number of cells: " << n_in_elm << endl;
}

void Grid::readInput()
{
    string tmps;
    
    ifstream in;
    string inputDir;
    string ids = to_string (id);
    inputDir = "./Grid_";
    inputDir.append (ids);
    inputDir.append ("/input.dat");
    in.open(inputDir);
    
    in >> tmps; in >> meshFile;
    in >> tmps; in >> phys_count;
    
    phys.resize(phys_count);
    bc.resize(phys_count);

    for (int i=0; i<phys_count; ++i)
    {
        in >> tmps; in >> phys[i];
    }
    
    for (int i=0; i<phys_count; ++i)
    {
        in >> tmps; in >> bc[i];
    }
    
    in >> tmps; in >> nHoles;
    holes.resize(nHoles);
    
    for (int i=0; i<nHoles; ++i)
    {
        for (int j=0; j<N_DIM; ++j)
        {
            in >> tmps; in >> holes[i].min[j];
            in >> tmps; in >> holes[i].max[j];
        }
    }
    
    in.close();
}

void Grid::printInput()
{    
    cout << meshFile << endl;
}

void Grid::read_grid ()
{
    readInput();
    printInput();
    read_ptSize ();
    read_pt();
    read_elmSize ();
    read_elm ();    
    printMeshInfo();
    in.close();
}



