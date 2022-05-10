#ifndef XVTKFILE_HH
#define XVTKFILE_HH

#include <vector>
#include <string>

using namespace std;

// @brief vtk cell types
enum Vtk_Type {
  VTK_VERTEX         = 1,
  VTK_POLY_VERTEX    = 2,
  VTK_LINE           = 3,
  VTK_POLY_LINE      = 4,
  VTK_TRIANGLE       = 5,
  VTK_TRIANGLE_STRIP = 6,
  VTK_POLYGON        = 7,
  VTK_PIXEL          = 8,
  VTK_QUAD           = 9,
  VTK_TETRA          = 10,
  VTK_VOXEL          = 11,
  VTK_HEXAHEDRON     = 12,
  VTK_WEDGE          = 13,
  VTK_PYRAMID        = 14
};

// @file   XVtkFile.h
// @author Y.Yang (QST)
// @date   2022.04.26 (created)

class XVtkFile
{
  public:
    /// @brief constructor
    XVtkFile(const char* filename);

    /// @brief deconstructor
    ~XVtkFile(){}

    /// @brief load input file
    void Load(const char* filename);

    /// @brief return the number of cell type
    vector<int> ListOfCellType();

    /// @brief return the number of cell id
    vector<int> ListOfCellId();

    /// @brief return the name of vtk type
    string GetNameOfType(const int type) const; 

    /// @brief print out information about vtk file
    void Print();

    /// @brief return the node points
    vector<vector<double>> GetPoints();

    /// @brief return the element cells
    vector<vector<int>> GetCells();

    /// @brief return the element types
    vector<int> GetCellTypes();

    /// @brief return the element id
    vector<int> GetCellId();

    /// @brief get the number of cells
    int GetNumOfCells() const { return fCells.size(); }

    /// @brief return types and cells
    void GetCellInfo(const int idx, vector<int>& nodes, int& type, int& id);

    /// @brief get the number of points
    int GetNumOfPoints() const { return fPoints.size(); }

    /// @brief return the point at the given node
    void GetPointInfo(const int idx, double& x, double& y, double& z);

    /// @brief print out information of node
    void PrintNode();

    /// @brief print out information of elements
    void PrintElement();


  protected:
    /// @brief read input file to 2D vector
    void read_input_file(const char* filename, vector<vector<string>>& table);

    /// @brief fill the point into array
    void fill_points_vector(const vector<vector<string>>& dataline);

    /// @brief fill the cell into array
    void fill_cells_vector(const vector<vector<string>>& dataline);

    /// @brief fill the array containing cell types
    void fill_cell_types(const vector<vector<string>>& dataline);

    /// @brief fill the array containing the cell of surface and line of boundaries
    void fill_surface_and_boundary(const vector<vector<string>>& dataline);

  private:
    vector<vector<double>> fPoints;
    vector<vector<int>> fCells;
    vector<int> fTypes;
    vector<int> fCellId;
};

#endif
