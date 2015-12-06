/* 
 * File:   Cell.h
 * Author: Orhan Shibliyev
 *
 * Created on June 13, 2014, 3:40 AM
 */

#ifndef CELL_H
#define	CELL_H

#include <vector>
#include <memory>
#include <array>
#include "../Vector/Vector.h"
#include "../Constants.h"
#include "../Matrix5/Matrix5.h"
#include "../Point/Point.h"

using std::vector;
using std::array;
using std::reference_wrapper;
using std::move;

class Face; // forward declaration

enum class elmType_t {UNDEFINED=-1, TRI=2, QUAD=3, TET=4, HEX=5, PEN=6};
enum class nFaces_t {UNDEFINED=-1, TRI=1, QUAD=1, TET=4, PEN=5, HEX=6};
enum class nVertices_t {UNDEFINED=-1, TRI=3, QUAD=4, TET=4, PEN=6, HEX=8};
enum class iBlank_t {UNDEFINED, NA, HOLE, FRINGE, FIELD};
enum class BC {UNDEFINED=-2, NA=-1, EMPTY=0, SLIP_WALL=1, DIRICHLET=2};
enum class vtkCellType_t {UNDEFINED=-1, TET=10, HEX=12, WEDGE=13, TRI=5};
enum class fringeBou_t {UNDEFINED=-1, NO=0, YES=1};

struct Cell
{
    // Fields
    int phys;
    Cell* donor;
    vector<Cell*> receiver;
    iBlank_t iBlank;
    int belonging;
    int nTrims;
    int nVertices;
    int nFaces;
    vector <int> nei;
    double Mach;
    double sigma;
    double wallDistance;
    double vol;
    double r_11, r_12, r_13, r_22, r_23, r_33;
    Vector<N_VAR> R;
    Vector<N_VAR> dQ, old_dQ;
    Vector<N_VAR> prim, cons, old_cons, oldold_cons;
    Vector<N_VAR> resInner, resOuter;
    Matrixd<N_VAR, N_VAR> D;
    bool trim;
    bool newlyCreated;
    bool ghost;
    Vector<3> cnt;
    Vector2D <N_DIM,N_VAR> grad;
    Vector2D <N_DIM,N_VAR> emin;
    Vector2D <N_DIM,N_VAR> emax;
    vector <int> face;    
    elmType_t type;
    vector <int> vtx;
    vector <int> vtxBelo;
    BC bc;    
    fringeBou_t fringeBou;

    // Constructors
    Cell();
    Cell (const Cell& other);
    Cell& operator= (const Cell& other);
    Cell (Cell&& other);

    // Methods
    void set_centroid (const vector<Point>& pt);
    void prim_to_cons();
    void cons_to_prim();
    void interpolate();
};

#endif	/* ELEMENT_H */

