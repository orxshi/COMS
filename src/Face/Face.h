/* 
 * File:   Face.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 12:30 AM
 */

#ifndef FACE_H
#define	FACE_H

#include <string>
#include <vector>
#include <array>
#include "../Cell/Cell.h"

#define N_NEI_FACE 2

using std::vector;
using std::array;
using std::string;
using std::reference_wrapper;

enum class nVerticesFace_t {UNDEFINED=-1, TRI=3, QUAD=4};
enum class face_t {UNDEFINED=-1, BOUNDARY=0, INTERIOR=1};

struct Face
{
    vector <int> vtx;
    vector <int> nei;
    CVector area;
    CVector cnt;
    CVector vb;
    face_t bouType;
    //array < Matrix5, N_NEI_FACE > M; // left matrix, right matrix
    //array < Matrixd<N_VAR,N_VAR>, N_NEI_FACE > M; // left matrix, right matrix
    
    // Methods
    Face();
    void set_area (const vector<Point>& pt);
    void set_centroid (const vector<Point>& pt);
};

#endif	/* FACE_H */

