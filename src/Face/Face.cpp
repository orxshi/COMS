#include "Face.h"

Face::Face()
{
    vb[0] = 0.;
    vb[1] = 0.;
    vb[2] = 0.;
}

void Face::set_area (const vector<Point>& pt)
{
    // For area of triangle and quad

    unsigned int tri  = static_cast<int>(nVerticesFace_t::TRI);
    unsigned int quad = static_cast<int>(nVerticesFace_t::QUAD);

    if (vtx.size() == tri || vtx.size() == quad)
    {
        area = crossP (pt[vtx[1]].dim - pt[vtx[0]].dim, pt[vtx.back()].dim - pt[vtx[0]].dim);

        if (vtx.size() == tri)
        {
            area *= 0.5;
        }
    }
}

void Face::set_centroid (const vector<Point>& pt)
{
    for (unsigned int i=0; i<cnt.size(); ++i)
    {
        cnt[i] = 0.;
    }

    for (const int p: vtx)
    {        
        cnt += pt[p].dim;
    }

    cnt /= vtx.size();
}