//=INPUTS=====================================================

H = 30;
D = 1;
mx = 100;  //10 // 100
my = 101;   //6 // 51
mz = 1;

//============================================================

ori = newp; // origin
Point(ori)  = {-H,-H,0};
ul = newp;
Point(ul) = {-H,H,0};
ll = newl;
Line(ll) = {ori,ul};

Transfinite Line{ll} = my;

outA[] = Extrude {2*H,0,0} {Line{ll}; Layers{mx}; Recombine;};

Printf("top line = %g", outA[0]);
Printf("surface = %g", outA[1]);
Printf("side lines = %g and %g", outA[2], outA[3]);

outB[] = Extrude {0,0,D} {Surface{outA[1]}; Layers{mz}; Recombine;};

Printf("top surface = %g", outB[0]);
Printf("volume = %g", outB[1]);
Printf("side surfaces = %g and %g and %g and %g", outB[2], outB[3], outB[4], outB[5]);
Printf("base surfaces = %g", outA[0]);

Physical Surface(1) = {outB[0], outA[1]}; // laterals
Physical Surface(2) = {outB[2], outB[4], outB[3], outB[5]}; // sides
Physical Volume(1) = {outB[1]};


