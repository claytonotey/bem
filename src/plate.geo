
m0 = .4;
m1 = .2;
m2 = .1;

w = 2;
L = 5;
d = .5;
H = .8;


Point(1) = {0, 0, 0, m0};
Point(2) = {w-d, 0, 0, m1};
Point(3) = {w, 0, 0, m2};
Point(4) = {w, L, 0, m2};
Point(5) = {w-d, L, 0, m1};
Point(6) = {0, L, 0, m0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line Loop(1000) = {1, 2, 3, 4, 5, 6};

Plane Surface(1) = { 1000 }; 
surf = Translate {-w, -1, 0} { Surface{1}; };
Physical Volume(1) = Extrude { 0, 0, H } { Surface{surf}; }