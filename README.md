# coilpattern
Code for design of coil winding pattern

## Compile
Please install the following libraries:
- Eigen
- boost
- OpenMPI

Run the following command to build the code
`mkdir build`
`cmake <Directory of the code>`
`make`

# Use the code
1. Mesh cutting: geometry must be meshed in grid, otherwise result gives large calculation error so that the smooth winding pattern cannot be obtained. The code read the mesh in `vtk` format. There are examples in the directory of `input/vtk` and the input files to create `vtk` files are saved in `input/gmsh`. The following files are meshed in grid.
- `curvedGeo.vtk`: curved elliptical geometry
- `mri_cylinder.vtk`: cylinder geometry

2. Construction: the code for the calculation of winding pattern is written in the file: `examples/test.cxx`. The following code is to load the vtk file and pass to a mesh handler.
```cpp
XVtkFile* file = new XVtkFile( filename );
XMeshHandle2* handle = new XMeshHandle2;
handle->SetVtkFile( file );
delete file;

// input the mesh handler into matrix
XMatrix3Node* mat = new XMatrix3Node;
mat->SetMeshInfo(handle);
```
Then the target region and the magnetic field that has to be optimized is inputted to class `XMatrix3Node` point by point. The target field should be contained into a `VectorXd` class. The response matrix containing the components of magnetic field is obtained from the class `XMatrix3Node`.
```cpp
VectorXd b = VectorXd :: Ones( numOfPoints );

for (int i=0; i<numOfPoints; i++) {
  mat->SetPointAtTarget(i, x.at(i), y.at(i), z.at(i));
  b(i) = byy.at(i);
}

MatrixXd bx_Mat, by_Mat, bz_Mat;
mat->GetFieldMatrix( bx_Mat, by_Mat, bz_Mat );
```
In final, using the class `XSolveTSVD` to solve the problem as follows.
```cpp
XSolveTSVD* solver = new XSolveTSVD;
solver->SetMeshHandle(handle);
solver->SetTruncatedNumber(120);
solver->SetSaveFile("res_tsvd");
solver->Solve( by_mat, b, x );
```
The truncated number has the maximum number equaling to the total points of target field. To obtain the more precision solution, please increase the value of truncated number.

3. Runing the code: the code can be built by user. There is an example in `examples/test.cxx`. To run the test example, please run the following command.
```
./test <directory>/<vtk file>
```

4. The results are written in vtk file. The result summed upto the truncated number is saved in `res_tsvd_mode_998.vtk` file. Please using Paraview to read the files.


