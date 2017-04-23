N = 3;

m0 = .3;
m1 = .15;
m2 = .08;

// radius of endcaps
r = .1;
// spacing between fingers
s = .3;
// width of backing plate
w = .2;
// period
p = (2 * r + s);
// total length
L = N * p;
// length of fingers
d = .4;
// z height
H = .8;

Point(1) = {0, 0, 0, m1};
Point(2) = {w, 0, 0, m0};
Point(3) = {w, L, 0, m0};
Point(4) = {0, L, 0, m1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

i = 0;
Point(5) = {0,L-i*p-s,0,m1};
Point(6) = {-d,L-i*p-s,0,m2};
Point(7) = {-d,L-i*p-s-r,0,m2};
Point(8) = {-d,L-i*p-s-r-r,0,m2};

Line(4) = {4, 5};
Line(5) = {5, 6};
Circle(6) = {6, 7, 8};

i = 1;
Point(9) = {0,L-i*p,0,m1};
Point(10) = {0,L-i*p-s,0,m1};
Point(11) = {-d,L-i*p-s,0,m2};
Point(12) = {-d,L-i*p-s-r,0,m2};
Point(13) = {-d,L-i*p-s-r-r,0,m2};

Line(7) = {8, 9};
Line(8) = {9, 10};
Line(9) = {10, 11};
Circle(10) = {11, 12, 13};

i = 2;
Point(14) = {0,L-i*p,0,m1};
Point(15) = {0,L-i*p-s,0,m1};
Point(16) = {-d,L-i*p-s,0,m2};
Point(17) = {-d,L-i*p-s-r,0,m2};
Point(18) = {-d,L-i*p-s-r-r,0,m2};

Line(11) = {13, 14};
Line(12) = {14, 15};
Line(13) = {15, 16};
Circle(14) = {16, 17, 18};

Line(15) = {18, 1};

Line Loop(1000) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

Plane Surface(1) = { 1000 }; 
surf = Translate {d+r, 0, 0} { Surface{1}; };
Physical Volume(1) = Extrude { 0, 0, H } { Surface{surf}; }