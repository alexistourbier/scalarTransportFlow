// // Point(1) = {0, 0, 0, 500};
// // Point(2) = {0.1, 0, 0, 500};
// // Point(3) = {0.1, 0.1, 0, 500};
// // Point(4) = {0, 0.1, 0, 500};
// // Line(1) = {1, 2};
// // Line(2) = {2, 3};
// // Line(3) = {3, 4};
// // Line(4) = {4, 1};
// // Line Loop(6) = {4, 1, 2, 3};
// // Plane Surface(6) = {6};
// // Physical Volume("internal") = {1};
// // Extrude {0, 0, 0.1} {
// //     Surface{6};
// //     Layers{1};
// //     Recombine;
// // }
// // Physical Surface("front") = {28};
// // Physical Surface("back") = {6};
// // Physical Surface("bottom") = {27};
// // Physical Surface("left") = {15};
// // Physical Surface("top") = {19};
// // Physical Surface("right") = {23};
// // Physical Volume("internal") = {1};


// // Inputs
// squareSide = 0.1; //m
// meshThickness = squareSide / 10; 
// gridsize = squareSide / 500;

// // All numbering counterclockwise from bottom-left corner
// Point(1) = {0, 0, 0, gridsize};
// Point(2) = {squareSide, 0, 0, gridsize};
// Point(3) = {squareSide, squareSide, 0, gridsize};
// Point(4) = {0, squareSide, 0, gridsize};
// Line(1) = {1, 2};				// bottom line
// Line(2) = {2, 3};				// right line
// Line(3) = {3, 4};				// top line
// Line(4) = {4, 1};				// left line
// Line Loop(5) = {1, 2, 3, 4}; 	
// // the order of lines in Line Loop is used again in surfaceVector[]
// Plane Surface(6) = {5};

// surfaceVector[] = Extrude {0, 0, meshThickness} {
//  Surface{6};
//  Layers{1};
//  Recombine;
// };
// /* surfaceVector contains in the following order:
// [0]	- front surface (opposed to source surface)
// [1] - extruded volume
// [2] - bottom surface (belonging to 1st line in "Line Loop (6)")
// [3] - right surface (belonging to 2nd line in "Line Loop (6)")
// [4] - top surface (belonging to 3rd line in "Line Loop (6)")
// [5] - left surface (belonging to 4th line in "Line Loop (6)") */
// Physical Surface("front") = surfaceVector[0];
// Physical Volume("internal") = surfaceVector[1];
// Physical Surface("bottom") = surfaceVector[2];
// Physical Surface("right") = surfaceVector[3];
// Physical Surface("top") = surfaceVector[4];
// Physical Surface("left") = surfaceVector[5];
// Physical Surface("back") = {6}; // from Plane Surface (6) ...
    


Point(1) = {0, 0, 0, 0.01};
Point(2) = {5, 0, 0, 0.01};
Point(3) = {5, 3, 0, 0.01};
Point(4) = {0, 3, 0, 0.01};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Physical Volume("internal") = {1};
Extrude {0, 0, 0.1} {
 Surface{6};
 Layers{1};
 Recombine;
}
Physical Surface("front") = {28};
Physical Surface("back") = {6};
Physical Surface("bottom") = {27};
Physical Surface("left") = {15};
Physical Surface("top") = {19};
Physical Surface("right") = {23};