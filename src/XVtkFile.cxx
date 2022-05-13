#include <iostream>
#include <iomanip>
#include <fstream>
#include "XVtkFile.h"
#include "XLogger.h"

XVtkFile :: XVtkFile(const char* filename)
{
  Info("VTK INPUT FILE: " << filename << "."); 
  Load(filename);
}

void XVtkFile :: Load(const char* filename)
{
  // load vtk file
  vector<vector<string>> dataline;
  read_input_file( filename, dataline );

  Info("CHECK NODE POINTS.");
  fill_points_vector( dataline );

  Info("CHECK ELEMENTS.");
  fill_cells_vector( dataline );

  Info("CHECK TYPE OF ELEMENT.");
  fill_cell_types( dataline );

  Info("CHECK BOUNDARIES AND SURFACES.");
  fill_surface_and_boundary( dataline );

  // check types of cell
  Print();
}

void XVtkFile :: read_input_file(const char* filename, vector<vector<string>>& table)
{
  ifstream inputfile( filename );

  // check file exists or not
  if ( !inputfile.is_open() ) {
    cerr << "ERROR: Could not open the file - " << filename << endl;
    Fatal("CANNOT OPEN THE VTK FILE: " << filename << ".");
    throw invalid_argument("could not open the file.");
  }

  cout << "LOADING VTK MESH FILE - " << filename << endl;
  string eachline;

  while (getline(inputfile, eachline)) {
    // put each line into string stream
    istringstream iss( eachline );
    string word;
    vector<string> elements;

    // read word by word
    while (iss >> word) 
      elements.push_back( word );

    // filled each element array into table
    table.push_back( elements );
  }
}

void XVtkFile :: fill_points_vector(const vector<vector<string>>& dataline)
{
  int start{-99}, npoints{-99};
  double x, y, z;

  // search for the number of points and start point
  for (int i=0; i<dataline.size(); i++) {
    if (dataline.at(i).size()!=0 && dataline.at(i).at(0)=="POINTS") {
      if (start<0) npoints = stoi( dataline.at(i).at(1) );
      start = i + 1;
      break; 
    }
  }

  cout << " - TOTAL NUMBER OF POINTS: " << npoints << endl;
  cout << " - RANGE OF POINTS LINE: (" << start << "," << start+npoints << ")" << endl;
  Info("READING THE POINTS FROM LINE " << start << " TO LINE " << start+npoints << ".");
  Info("TOTAL NUMBER OF POINTS: " << npoints << ".");

  // fill the array of point
  for (int i=start; i<start+npoints; i++) {
    x = stod( dataline.at(i).at(0) );
    y = stod( dataline.at(i).at(1) );
    z = stod( dataline.at(i).at(2) );
    vector<double> pts{ x, y, z }; 

    fPoints.push_back( pts );
  }
}

void XVtkFile :: fill_cells_vector(const vector<vector<string>>& dataline)
{
  int start{-99}, ncells{-99};

  // search for the number of cells
  for (int i=0; i<dataline.size(); i++) {
    if (dataline.at(i).size()!=0 && dataline.at(i).at(0)=="CELLS") {
      if (start<0) ncells = stoi( dataline.at(i).at(1) );
      start = i + 1;
      break; 
    }
  }

  cout << " - TOTAL NUMBER OF CELLS: " << ncells << endl;
  cout << " - RANGE OF CELLS LINE: (" << start << "," << start+ncells << ")" << endl;
  Info("READING THE ELEMENTS FROM LINE " << start << " TO LINE " << start+ncells << ".");
  Info("TOTAL NUMBER OF ELEMENTS: " << ncells << ".");

  // fill the array of cell
  for (int i=start; i<start+ncells; i++) {
    vector<int> ptsId( stoi(dataline.at(i).at(0)) );

/*
    if (i<5725+250) {
      cout << "i = " << i << ",";
      for (int j=0; j<ptsId.size(); j++) cout << " " << dataline.at(i).at(j+1); 
      cout << "\n";
    }
*/

    for (int j=0; j<ptsId.size(); j++)
      ptsId.at(j) = stoi( dataline.at(i).at(j+1) );

    fCells.push_back( ptsId );
  }
}

void XVtkFile :: fill_cell_types(const vector<vector<string>>& dataline)
{
  int start{-99}, ncells{-99};

  // search for the number of cell types
  for (int i=0; i<dataline.size(); i++) {
    if (dataline.at(i).size()!=0 && dataline.at(i).at(0)=="CELL_TYPES") {
      if (start<0) ncells = stoi( dataline.at(i).at(1) );
      start = i + 1;
      break; 
    }
  }

  cout << " - TOTAL NUMBER OF CELL TYPES: " << ncells << endl;
  cout << " - RANGE OF CELL TYPE LINE: (" << start << "," << start+ncells << ")" << endl;
  Info("READING THE TYPES FROM LINE " << start << " TO LINE " << start+ncells << ".");
  Info("TOTAL NUMBER OF TYPES: " << ncells << ".");

  // fill the array of cell type
  for (int i=start; i<start+ncells; i++)
    fTypes.push_back( stoi(dataline.at(i).at(0)) );
}

void XVtkFile :: fill_surface_and_boundary(const vector<vector<string>>& dataline)
{
  int start{-99}, ncells{-99};

  // search for the number of cell types
  for (int i=0; i<dataline.size(); i++) {
    if (dataline.at(i).size()!=0 && dataline.at(i).at(0)=="CELL_DATA") {
      if (start<0) ncells = stoi( dataline.at(i).at(1) );
      start = i + 3;
      break; 
    }
  }

  cout << " - TOTAL NUMBER OF CELL ID: " << ncells << endl;
  cout << " - RANGE OF CELL ID LINE: (" << start << "," << start+ncells << ")" << endl;
  Info("READING THE ELEMENT ID FROM LINE " << start << " TO LINE " << start+ncells << ".");
  Info("TOTAL NUMBER OF ELEMENT ID: " << ncells << ".");

  // fill the array of cell type
  if (start>0) {
    for (int i=start; i<start+ncells; i++)
      fCellId.push_back( stoi(dataline.at(i).at(0)) );
    cout << " - CELL ID IS FILLED." << endl;
  }
  else
    cout << " - CELL IS IS NOT DEFINED." << endl;
}

vector<int> XVtkFile :: ListOfCellType()
{
  vector<int> types{fTypes.at(0)};
  bool is_exist = false;

  for (int i=1; i<fTypes.size(); i++) {
    is_exist = false;

    for (int j=0; j<types.size(); j++) {
      if (fTypes.at(i) == types.at(j)) { is_exist = true; break; }
    }

    if (!is_exist) types.push_back( fTypes.at(i) );
  }

  return types;
}

vector<int> XVtkFile :: ListOfCellId()
{
  vector<int> id;
  bool is_exist = false;

  for (int i=0; i<fCellId.size(); i++) {
    is_exist = false;

    for (int j=0; j<id.size(); j++) {
      if (fCellId.at(i) == id.at(j)) { is_exist = true; break; }
    }

    if (!is_exist) id.push_back( fCellId.at(i) );
  }

  return id;
}

string XVtkFile :: GetNameOfType(const int type) const
{
  string name = "";

  switch(type) {
    case 1 : name = "VTK_VERTEX"        ; break;
    case 2 : name = "VTK_POLY_VERTEX"   ; break;
    case 3 : name = "VTK_LINE"          ; break;
    case 4 : name = "VTK_POLY_LINE"     ; break;
    case 5 : name = "VTK_TRIANGLE"      ; break;
    case 6 : name = "VTK_TRIANGLE_STRIP"; break;
    case 7 : name = "VTK_POLYGON"       ; break;
    case 8 : name = "VTK_PIXEL"         ; break;
    case 9 : name = "VTK_QUAD"          ; break;
    case 10: name = "VTK_TETRA"         ; break;
    case 11: name = "VTK_VOXEL"         ; break;
    case 12: name = "VTK_HEXAHEDRON"    ; break;
    case 13: name = "VTK_WEDGE"         ; break;
    case 14: name = "VTK_PYRAMID"       ; break;
    default: name = "UNKNOWN"           ; break;
  }

  return name;
}

void XVtkFile :: Print()
{
  vector<int> types = ListOfCellType();
  vector<int> id    = ListOfCellId();

  // print out cell types
  cout << "LISTS OF CELL TYPE:" << endl;
  for (int i=0; i<types.size(); i++)
    cout << "TYPE:" << setw(5) << fixed << types.at(i) << setw(16) << fixed << GetNameOfType(types.at(i)) << endl; 

  // print out cell id
  cout << "LISTS OF CELL ID:" << endl;
  for (int i=0; i<id.size(); i++)
    cout << "ID:" << setw(5) << fixed << id.at(i) << endl;

/*
  for (int i=0; i<250; i++) {
    cout << "i = " << i << ", ";
    for (int j=0; j<fCells.at(i).size(); j++) cout << fCells.at(i).at(j) << " "; 
    cout << "\n";
  }
*/
}


void XVtkFile :: GetCellInfo(const int idx, vector<int>& nodes, int& type, int& id)
{
  if (idx<0 && idx>=GetNumOfCells()) {
    cerr << "ERROR: the index of cell is out of range." << endl;
    cerr << "       current index:" << idx << ", number of cells:" << GetNumOfCells() << endl;
    throw out_of_range( "the index of cell is out of range." );
  }

  nodes.clear();

  nodes = fCells .at(idx); 
  type  = fTypes .at(idx);
  id    = fCellId.at(idx);
}

void XVtkFile :: GetPointInfo(const int idx, double& x, double& y, double& z)
{
  if (idx<0 && idx>=GetNumOfPoints()) {
    cerr << "ERROR: the index of point is out of range." << endl;
    cerr << "       current index:" << idx << ", number of points:" << GetNumOfPoints() << endl;
    throw out_of_range( "the index of point is out of range." );
  }

  x = fPoints.at(idx).at(0);
  y = fPoints.at(idx).at(1);
  z = fPoints.at(idx).at(2);
}

void XVtkFile :: PrintNode()
{
  cout << "PRINT OUT THE INFORMATION OF NODES." << endl;
  cout << setw( 6) << fixed << "ID"
       << setw(15) << fixed << "X(m)" 
       << setw(15) << fixed << "Y(m)"
       << setw(15) << fixed << "Z(m)" << endl;

  for (int i=0; i<fPoints.size(); i++)
    cout << setw( 6) << fixed << i
         << setw(15) << setprecision(6) << scientific << fPoints.at(i).at(0)
         << setw(15) << setprecision(6) << scientific << fPoints.at(i).at(1)
         << setw(15) << setprecision(6) << scientific << fPoints.at(i).at(2) << endl;
}

void XVtkFile :: PrintElement()
{
  cout << "PRINT OUT THE INFORMATION OF ELEMENTS." << endl;
  cout << setw( 6) << fixed << "ID"
       << setw(15) << fixed << "TYPE"
       << setw( 6) << fixed << "LABEL"
       << setw( 7) << fixed << "NODE_I" << endl;

  for (int i=0; i<fCells.size(); i++) {
    cout << setw( 6) << fixed << i
         << setw(15) << fixed << GetNameOfType(fTypes.at(i))
         << setw( 6) << fixed << fCellId.at(i);
    for (int j=0; j<fCells.at(i).size(); j++) 
      cout << setw(7) << fixed << fCells.at(i).at(j);
    cout << "\n";
  }
}
