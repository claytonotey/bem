
m0 = 5; 
m1 = 5;
m2 = 3;

w = 80;
L = 60;
d = 20;
H = .27;


Point(1) = {0,    0, 0, m0};
Point(2) = {-w+d, 0, 0, m1};
Point(3) = {-w,   0, 0, m2};
Point(4) = {-w,   L, 0, m2};
Point(5) = {-w+d, L, 0, m1};
Point(6) = {0,    L, 0, m0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line Loop(1000) = {1, 2, 3, 4, 5, 6};

Plane Surface(1) = { 1000 }; 
surf = Translate {w, 0, 0} { Surface{1}; };
zdir[] = Extrude { 0, 0, H } { Surface{surf}; Layers{1}; };
