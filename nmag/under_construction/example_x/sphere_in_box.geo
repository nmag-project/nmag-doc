Include "sphere.geo";
lc2 = 1.;

cs = 5.;  //Cube Size 
c = cs/2.; // distance from centre to corner

Point(100) = {-cs, -cs, -cs, lc2 };
Point(101) = { cs, -cs, -cs, lc2 };
Point(102) = {-cs,  cs, -cs, lc2 };
Point(103) = { cs,  cs, -cs, lc2 };
Point(104) = {-cs, -cs,  cs, lc2 };
Point(105) = { cs, -cs,  cs, lc2 };
Point(106) = {-cs,  cs,  cs, lc2 };
Point(107) = { cs,  cs,  cs, lc2 };

Line(101) = {100,101};
Line(102) = {102,103};
Line(103) = {104,105};
Line(104) = {106,107};



Line(105) = {100,102};
Line(106) = {101,103};
Line(107) = {105,107};
Line(108) = {104,106};

Line Loop(101) = {101,106,-102,-105};
Line Loop(102) = {103,107,-104,-108};

Plane Surface(101) = {101};
Plane Surface(102) = {102};

Line(109) = {103,107};
Line(110) = {102,106};

Line Loop(103) = {-109,-102,110,104};
Plane Surface(103) = {103};

Line(111) = {101,105};

Line Loop(104) = {111,107,-109,-106};
Plane Surface(104) = {104};

Line(112) = {100,104};

Line Loop(105) = {101,111,-103,-112};
Plane Surface(105) = {105};

Line Loop(106) = {112,108,-110,-105};
Plane Surface(106) = {106};

Surface Loop(101) = {101, 102, 103, 104, 105, 106, 7, 8, 9, 10};

Volume(2) = {101};

Physical Volume(2) =  {2};

