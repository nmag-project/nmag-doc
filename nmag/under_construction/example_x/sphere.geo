lc = 0.5;
Point(2) = {2, 0.0, 0, lc };
Point(1) = {0, 0, 0, lc };
Point(3) = {0, 2, 0, lc };
Point(4) = {-2, 0, 0, lc };
Point(5) = {-0, -2, 0, lc };
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};


Rotate {{0, 1, 0}, {0, 0, 0}, Pi/2} {
       Duplicata{ Line{1,2,3,4}; } 
}

Line Loop(11) = {1, 6, 7, 4};
Line Loop(12) = {-3, 7, 6, -2};
Line Loop(13) = {2, 3, 5, 8};
Line Loop(14) = {-4, -1, 5, 8};

Ruled Surface(7) = {11} ;
Ruled Surface(8) = {12} ;
Ruled Surface(9) = {13} ;
Ruled Surface(10) = {14};

Surface Loop(15) = {10, 7, 8, 9};

Volume(1) = {15};

Physical Volume(1) =  {1};
