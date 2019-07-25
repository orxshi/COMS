#include "MovingMesh.h"

MovingGrid::MovingGrid (double dt)
{
    this->dt = dt;
}

int MovingGrid::getEdge (const vector<Point>& pt, const Face& f)
{
    int success = 0;
    vector <int> tempVectices;
    CVector refVec;
    CVector testVec;
    CVector faceNorm;
    
    for (int v: f.vtx)
    {
        if (pt[v].dim[2] == 0.)
        {
            tempVectices.push_back (v);
        }
    }
    
    if (tempVectices.size() == 2)
    {
        //cout << "size should be two in MovingGrid::getEdge(...)" << endl;
        //cout << "size = " << tempVectices.size() << endl;
        //exit(-2);
        
        success = 1;
        
        vertex0 = pt[tempVectices[0]].dim;
        vertex1 = pt[tempVectices[1]].dim;

        faceNorm = norm(f.area);

        refVec[0] = faceNorm[1];
        refVec[1] = -faceNorm[0];
        refVec[2] = 0.;

        testVec = vertex1 - vertex0;

        if (dotP(testVec, refVec) < 0.)
        {
            vertex0 = pt[tempVectices[1]].dim;
            vertex1 = pt[tempVectices[0]].dim;
        }
        else if (dotP(testVec, refVec) == 0.)
        {
            cout << "cannot be zero in getEdge(...)" << endl;
            exit(-2);
        }
    }
            
    
    
    return success;
}

void MovingGrid::getFaceVelocity(Face& f)
{
    double areaSwept = 0.;
    double edgeLength;
    CVector faceNorm;
    CVector vertices[4];
    
    vertices[0] = vertex0;
    vertices[1] = vertex1;
    vertices[2] = vertex2;
    vertices[3] = vertex3;
    
    faceNorm = norm(f.area);
    edgeLength = mag(vertex1 - vertex0);
    
    for (int i=0; i<3; ++i)
    {
        areaSwept += ( vertices[i][0] * vertices[i+1][1] ) - ( vertices[i+1][0] * vertices[i][1] );
    }

    areaSwept += vertices[3][0] * vertices[0][1];
    areaSwept -= vertices[0][0] * vertices[3][1];
    areaSwept *= 0.5;

    f.vb = (areaSwept / dt / edgeLength) * faceNorm;
    /*cout << areaSwept << endl;
    cout << dt << endl;
    cout << edgeLength << endl;
    exit(-2);*/
}

void MovingGrid::checkGCL(const Grid& gr)
{
    double maxErr = BIG_NEG_NUM;
    for (int e=gr.n_bou_elm; e<gr.cell.size(); ++e)
    {
        double sum = 0.;
        
        for (int iF: gr.cell[e].face)
        {
            const Face& f = gr.face[iF];
            
            double dp = dotP( f.vb,f.area );

            if (e == f.nei[0])
            {
                sum += dp;
            }
            else if (e == f.nei[1])
            {
                sum -= dp;
            }
            else
            {
                cout << "neither nei[0] nor nei[1] in getMovingFaceVelocity(...)" << endl;
                exit(-2);
            }
        }

        maxErr = max (fabs(sum),maxErr);

        if (fabs(sum) > 1e-11)
        {
            cout << "GCL not satisfied in getMovingFaceVelocity(...)" << endl;
            cout << fabs(sum) << endl;
            
            for (int iF: gr.cell[e].face)
            {
                const Face& f = gr.face[iF];
                
                cout << "size = " << gr.cell[e].face.size() << endl;
                cout << "e = " << e << endl;
                
                cout << "f.vb[0] = " << f.vb[0] << endl;
                cout << "f.vb[1] = " << f.vb[1] << endl;
                cout << "f.vb[2] = " << f.vb[2] << endl;
                
                cout << "f.area[0] = " << f.area[0] << endl;
                cout << "f.area[1] = " << f.area[1] << endl;
                cout << "f.area[2] = " << f.area[2] << endl;
            }
            
            exit(-2);
        }
    }
}

void MovingGrid::getAllFaceVelocities (Grid& gr)
{
    for (Face& f: gr.face)
    {
        int success = getEdge (gr.pt, f);
        
        if (success == 0)
        {
            f.vb[0] = 0.;
            f.vb[1] = 0.;
            f.vb[2] = 0.;
        }
        else
        {
            moveEdge();
            getFaceVelocity(f);
        }
    }
    
    //checkGCL(gr);
    
    gr.apply_BCs();
}

void MovingGrid::rotateVectorAroundPoint2D (const CVector pivot, const double angleRad, const CVector r, CVector& ans)
{
    ans[0] = pivot[0] + r[0] * cos(angleRad) - r[1]*sin(angleRad);
    ans[1] = pivot[1] + r[0] * sin(angleRad) + r[1]*cos(angleRad);
}

void MovingGrid::displacePoint2D (const CVector orig, const CVector vel, const double dt, CVector& ans)
{
    ans[0] = orig[0] + vel[0] * dt;
    ans[1] = orig[1] + vel[1] * dt;
}



/*void getMovingFaceVelocity(Grid& gr)
{
    CVector r;
    vector <int> tempVectices;
    double areaSwept;
    CVector vertices[4];
    CVector norma;
    CVector refVec;
    CVector testVec;
    CVector faceNorm;

    double delAlphaRad = delAlpha * DEG_TO_RAD;
    double alphaRad = alpha * DEG_TO_RAD;

    for (Face& f: face)
    {
        tempVectices.clear();

        for (int v: f.vtx)
        {
            if (pt[v].dim[2] == 0.)
            {
                tempVectices.push_back (v);
            }
        }

        if (tempVectices.size() == 2)
        {
            vertices[0] = pt[tempVectices[0]].dim;
            vertices[1] = pt[tempVectices[1]].dim;
            
            faceNorm = norm(f.area);
            
            refVec[0] = faceNorm[1];
            refVec[1] = -faceNorm[0];
            refVec[2] = 0.;
            
            testVec = vertices[1].dim - vertices[0].dim;
            
            if (dotP(testVec, refVec) < 0.)
            {
                vertices[0] = pt[tempVectices[1]].dim;
                vertices[1] = pt[tempVectices[0]].dim;
            }
            else if (dotP(testVec, refVec) == 0.)
            {
                cout << "cannot be zero in getMovingFaceVelocity(...)" << endl;
                exit(-2);
            }

            r = vertices[1] - centerAirfoil;
            rotateVectorAroundPoint2D (centerAirfoil, delAlphaRad, r, vertices[2]);

            r = vertices[0] - centerAirfoil;
            rotateVectorAroundPoint2D (centerAirfoil, delAlphaRad, r, vertices[3]);

            CVector vel;
            vel[0] = MachAirfoil * cos (alphaRad);
            vel[1] = MachAirfoil * sin (alphaRad);
            vel[2] = 0.;

            displacePoint2D (vertices[2], vel, dt, vertices[2]);
            displacePoint2D (vertices[3], vel, dt, vertices[3]);

            areaSwept = 0.;
            for (int i=0; i<3; ++i)
            {
                areaSwept += ( vertices[i][0] * vertices[i+1][1] ) - ( vertices[i+1][0] * vertices[i][1] );
            }

            areaSwept += vertices[3][0] * vertices[0][1];
            areaSwept -= vertices[0][0] * vertices[3][1];
            areaSwept *= 0.5;

            double dx = vertices[1][0] - vertices[0][0];
            double dy = vertices[1][1] - vertices[0][1];

            //norma[0] = -dy;
            //norma[1] = dx;
            //norma[2] = 0.;

            f.vb = (areaSwept / dt / mag(norma)) * faceNorm;
            //f.vb = (areaSwept / dt / mag(norma)) * norm(norma);
        }
        else
        {
            f.vb[0] = 0.;
            f.vb[1] = 0.;
            f.vb[2] = 0.;
        }
    }

    double maxErr = BIG_NEG_NUM;
    for (int e=n_bou_elm; e<cell.size(); ++e)
    {
        double sum = 0.;

        for (int iF: cell[e].face)
        {
            Face& f = face[iF];
            
            double dp = dotP( f.vb,f.area );

            if (e == f.nei[0])
            {
                sum += dp;
//                sum += ( mag(grid.face[ID].Af) * grid.face[f].vbn );
            }
            else if (e == f.nei[1])
            {
                sum -= dp;
//                sum -= ( mag(grid.face[ID].Af) * grid.face[f].vbn );
            }
            else
            {
                cout << "neither nei[0] nor nei[1] in getMovingFaceVelocity(...)" << endl;
                exit(-2);
            }
        }

        maxErr = max (fabs(sum),maxErr);

        if (fabs(sum) > 1e-11)
        {
            cout << "GCL not satisfied in getMovingFaceVelocity(...)" << endl;
            cout << fabs(sum) << endl;
            exit(-2);
        }
    }
}*/


