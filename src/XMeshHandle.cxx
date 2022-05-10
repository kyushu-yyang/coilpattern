#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "XMeshHandle.h"

XMeshHandle :: XMeshHandle(const char* filename)
{
  Load( filename );
}

XMeshHandle :: ~XMeshHandle()
{
  if (fNodes.size()>0) {
    for (auto eachnode : fNodes) delete eachnode;
    fNodes.clear();
  }

  if (fElements.size()>0) {
    for (auto eachelem : fElements) delete eachelem;
    fElements.clear();
  }
}

void XMeshHandle :: Load(const char* filename)
{
  ifstream inputfile( filename );
  if ( !inputfile.is_open() ) {
    cerr << "Could not open the file - " << filename << endl;
    throw invalid_argument("could not open the file.");
  }

  cout << "loading mesh file - " << filename << endl;

  string eachline;
  int npoints, ncells, ntypes, point0, cell0, type0, cnt{0}, ecnt{0}, ncnt{0}, typecnt{0};
  int bcnt{0};

  while (getline(inputfile, eachline)) {
    // put each line into string stream
    istringstream iss( eachline );
    string word;
    vector<string> items;

    // read word by word
    while (iss >> word) 
      items.push_back( word );

    // check the points of node
    if (items.size()>0 && items.at(0)=="POINTS") {
      npoints = stoi(items.at(1));
      point0  = cnt;
      cout << " - number of points: " << npoints << endl;
      cout << " - range of points : (" << point0+1 << "," << point0+npoints+1 << ")" << endl; 
    }

    // fill the points into array
    if (cnt > point0 && cnt <= point0+npoints ) {
      XNode* eachnode = new XNode;
      eachnode->SetId(ncnt);
      eachnode->SetPoint( stod(items.at(0)), stod(items.at(1)), stod(items.at(2)) );
      fNodes.push_back(eachnode);
      ncnt ++;
    }
    
    // check the cells
    if (items.size()>0 && items.at(0)=="CELLS") {
      ncells = stoi(items.at(1));
      cell0  = cnt;
      cout << " - number of cells : " << ncells << endl;
      cout << " - range of cells  : (" << cell0+1 << "," << cell0+ncells+1 << ")" << endl; 
    }

    // fill the cells into array
    if (cnt > cell0 && cnt <= cell0+ncells) {
      // 3-node element
      if (items.size()-1==3) {
        XElement3Node* each_element = new XElement3Node;
        each_element->SetId(ecnt);

        for (int i=1; i<items.size(); i++ ) 
          each_element->SetNode( i-1, fNodes.at(stoi(items[i])) );

        fElements.push_back( each_element );
        ecnt ++;
      }

      // 2-node line
      /*
      if (items.size()-1==2) {
        XLine* each_bc = new XLine;
        each_bc->SetId(bcnt);
        each_bc->SetNode( 0, fNodes.at( stoi(items[1]) ) );
        each_bc->SetNode( 1, fNodes.at( stoi(items[2]) ) );
        fLines.push_back( each_bc );
        bcnt ++;
      }
      */
    }

    // search for the boundary
    if (items.size()>0 && items.at(0)=="CELL_DATA") {
      ntypes = stoi(items.at(1));
      type0  = cnt+2;
      cout << " - SEARCHING FOR BOUNDARIES" << endl;
    }

    if (cnt > type0 && cnt <= type0+ntypes) {
      // if the type of line (4) is found, then push back
      if (items.at(0)=="4") {
        cout << " - BC ID:" << typecnt << ", TYPE:" << items.at(0) << ", NODE1: " << endl;
      }
      typecnt ++;
    }
      
    cnt ++;
  }
}

void XMeshHandle :: load_vtk_file(const char* filename, vector<vector<string>>& table)
{
  ifstream inputfile( filename );
  if ( !inputfile.is_open() ) {
    cerr << "Could not open the file - " << filename << endl;
    throw invalid_argument("could not open the file.");
  }

  cout << "loading mesh file - " << filename << endl;
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

vector<XElement3Node*> XMeshHandle :: GetSurroundingElements(const int node)
{
  vector<XElement3Node*> element;
  int numOfNodes;

  for (int i=0; i<fElements.size(); i++) {
    numOfNodes = fElements.at(i)->GetNumOfNodes();

    // search the element containing the given node id
    for (int j=0; j<numOfNodes; j++) {
      if (node==fElements.at(i)->GetNode(j)->GetId())
        element.push_back(fElements.at(i));
    }
  }

  /*
  cout << "NODE:" << setw(8) << node << endl;
  cout << " - COORDINATE:" << fNodes.at(node)->GetX() << "," << fNodes.at(node)->GetY() << "," << fNodes.at(node)->GetZ() << endl;
  cout << " - SURROUNDING ELEMENTS:" << setw(4) << fixed << element.size() << endl;

  for (int i=0; i<element.size(); i++)
    cout << " - ELEMENT_" << i << ":" << setw(4) << fixed << element.at(i)->GetId() << endl;
  */

  return element;
}

vector<XNode*> XMeshHandle :: GetSurroundingNodes(const int element)
{
  vector<XNode*> node;
  const int numOfNodes = fElements.at(element)->GetNumOfNodes();
  for (int i=0; i<numOfNodes; i++)
    node.push_back( fElements.at(element)->GetNode(i) );

  /*
  cout << "ELEMENT:" << setw(8) << element << endl;
  cout << " - NUMBER OF NODES:" << setw(8) << numOfNodes << endl;
  for (int i=0; i<numOfNodes; i++)
    cout << " - NODE_" << i << ": (" << node.at(i)->GetX() << "," << node.at(i)->GetY() << "," << node.at(i)->GetZ() << ")" <<endl;
  */
  
  return node;
}
