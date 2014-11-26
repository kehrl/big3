lc1 = 0.00294106863992; 
lc2 = 0.00294106863992; 
lc3 = 0.0588213727983; 
Point(1) = {0, 0, 0, lc3}; 
Point(2) = {0.823535881605, 0, 0, lc2}; 
Point(3) = {0.941178627202, 0, 0, lc1}; 
Point(4) = {1, 0, 0, lc1}; 
Line(5) = {1, 2}; 
Line(6) = {2, 3}; 
Line(7) = {3, 4}; 
Extrude {0, 1, 0} {Line{5}; Layers{10}; Recombine;} 
Extrude {0, 1, 0} {Line{6}; Layers{10}; Recombine;} 
Extrude {0, 1, 0} {Line{7}; Layers{10}; Recombine;} 
Physical Surface(1) = {11,15,19}; 
Physical Line(1) = {5,6,7}; 
Physical Line(2) = {18}; 
Physical Line(3) = {16,12,8}; 
Physical Line(4) = {9}; 
