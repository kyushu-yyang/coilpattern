#include <iostream>
#include <iomanip>
#include "XMeshHandle2.h"
#include "XLogger.h"

XMeshHandle2 :: XMeshHandle2()
{}

XMeshHandle2 :: ~XMeshHandle2()
{
  if (fNode.size()>0) {
    for (auto eachnode : fNode) delete eachnode;
    fNode.clear();
  }

  if (fElement.size()>0) {
    for (auto eachelem : fElement) delete eachelem;
    fElement.clear();
  }
}

void XMeshHandle2 :: SetVtkFile(XVtkFile* file)
{
  setup_nodes(file);
  setup_boundaries(file);
  setup_elements(file);
}

vector<XElement3Node*> XMeshHandle2 :: GetSurroundingElements(const int node)
{
  vector<XElement3Node*> element;
  int numOfNodes;

  for (int i=0; i<fElement.size(); i++) {
    numOfNodes = fElement.at(i)->GetNumOfNodes();

    // search the element containing the given node id
    for (int j=0; j<numOfNodes; j++) {
      if (node==fElement.at(i)->GetNode(j)->GetId()) 
      { element.push_back(fElement.at(i)); break; }
    }
  }

  //sum_current_vector(node, element);

  /*
  cout << "NODE:" << setw(8) << node << endl;
  cout << " - COORDINATE:" << fNodes.at(node)->GetX() << "," << fNodes.at(node)->GetY() << "," << fNodes.at(node)->GetZ() << endl;
  cout << " - SURROUNDING ELEMENTS:" << setw(4) << fixed << element.size() << endl;

  for (int i=0; i<element.size(); i++)
    cout << " - ELEMENT_" << i << ":" << setw(4) << fixed << element.at(i)->GetId() << endl;
  */

  return element;
}

vector<XNode*> XMeshHandle2 :: GetSurroundingNodes(const int element)
{
  vector<XNode*> node;
  const int numOfNodes = fElement.at(element)->GetNumOfNodes();
  for (int i=0; i<numOfNodes; i++)
    node.push_back( fElement.at(element)->GetNode(i) );

  /*
  cout << "ELEMENT:" << setw(8) << element << endl;
  cout << " - NUMBER OF NODES:" << setw(8) << numOfNodes << endl;
  for (int i=0; i<numOfNodes; i++)
    cout << " - NODE_" << i << ": (" << node.at(i)->GetX() << "," << node.at(i)->GetY() << "," << node.at(i)->GetZ() << ")" <<endl;
  */
  
  return node;
}

void XMeshHandle2 :: setup_nodes(XVtkFile* file)
{
  const int numOfNodes = file->GetNumOfPoints();
  double x, y, z;

  for (int i=0; i<numOfNodes; i++) {
    file->GetPointInfo(i, x, y, z);

    // push back the setup of node 
    XNode* each_node = new XNode;
    each_node->SetId(i);
    each_node->SetPoint(x,y,z);
    fNode.push_back(each_node);
  }
}

void XMeshHandle2 :: setup_boundaries(XVtkFile* file)
{
  const int numOfCells = file->GetNumOfCells();
  vector<int> node;
  int type, label;
  int cnt = 0;

  for (int i=0; i<numOfCells; i++) {
    file->GetCellInfo( i, node, type, label );

    // search for the boundary
    if (label>0 && node.size()==2) { 
      fNode.at(node.at(0))->AtBoundary(label); 
      fNode.at(node.at(1))->AtBoundary(label);

      cout << "FOUND BOUNDARY AT CELL" << setw(5) << fixed << cnt << ":"
           << setw(6) << fixed << node.at(0)
           << setw(6) << fixed << node.at(1) << endl; 
      Debug( "BOUNDARY CELL" << setw(5) << fixed << cnt << ": ("
                             << setw(5) << fixed << node.at(0) << ") --- ("
                             << setw(5) << fixed << node.at(1) << "), BC LABEL:"
                             << setw(4) << fixed << label );
      cnt ++;
    }
  }
}

void XMeshHandle2 :: setup_elements(XVtkFile* file)
{
  const int numOfCells = file->GetNumOfCells();
  vector<int> node;
  int type, id;
  int cnt = 0;

  for (int i=0; i<numOfCells; i++) {
    file->GetCellInfo( i, node, type, id );
    
    // fill the element vectors
    if (id>0 && node.size()==3) {
      XElement3Node* each_element = new XElement3Node;
      each_element->SetId(i);

      // setup node pointer
      for (int j=0; j<node.size(); j++)
        each_element->SetNode( j, fNode.at(node.at(j)) );

      fElement.push_back( each_element );
      cnt ++;
    }
  }
}

bool XMeshHandle2 :: sum_current_vector(const int nid, const vector<XElement3Node*>& elements)
{
  Vector3d jsum(0.,0.,0.);
  bool is_loop = true;
  int  n_idx = -999;

  for (int i=0; i<elements.size(); i++) {
    jsum += elements.at(i)->GetLoopBasisVector( nid ); 
    /*
    n_idx = elements.at(i)->GetBoundaryNodeIndex();
    if (n_idx>0)
      jsum += elements.at(i)->GetLoopBasisVector( elements.at(i)->GetNode(n_idx)->GetId() );
    */
  }

  for (int i=0; i<3; i++) {
    if (fabs(jsum(i))>1e-9) {is_loop = false; break;} 
  }

  for (int i=0; i<elements.size(); i++) {
    for (int j=0; j<3; j++) {
      if (elements.at(i)->GetNode(j)->Is_AtBoundary()) { is_loop = true; break; }
    }
  }

  if (!is_loop) {
    cout << "ERROR: THE OPENED LOOP IS DETECTED." << endl;
    cout << " - NODE ID:" << setw(5) << fixed << nid << endl;
    cout << " - ELEMENT ID SURROUNDING THE NODE:";
    for (int i=0; i<elements.size(); i++)
      cout << setw(6) << fixed << elements.at(i)->GetId() << ",";
    cout << "\n";
    cout << "JSUM: " << setw(15) << setprecision(6) << scientific << jsum(0);
    cout << setw(15) << setprecision(6) << scientific << jsum(1);
    cout << setw(15) << setprecision(6) << scientific << jsum(2) << endl;
  }

  return is_loop;
}

void XMeshHandle2 :: Print()
{
  print_elements_sur();
}

void XMeshHandle2 :: print_elements_sur()
{
  vector<XElement3Node*> elem;
  Vector3d jn;
  Vector3d ln;

  cout << "PRINT OUT ELEMENTS SURROUNDING THE GIVEN NODE." << endl;
  cout << setw( 5) << fixed << "ID"
       << setw(14) << fixed << "JX"
       << setw(14) << fixed << "JY"
       << setw(14) << fixed << "JZ"
       << setw(14) << fixed << "VX"
       << setw(14) << fixed << "VY"
       << setw(14) << fixed << "VZ"
       << setw(10) << fixed << "ELEMENT_I" << endl;

  for (int i=0; i<fNode.size(); i++) {
    elem = GetSurroundingElements(i);
    cout << setw(5) << fixed << i;

    jn(0) = 0.; jn(1) = 0.; jn(2) = 0.;
    ln(0) = 0.; ln(1) = 0.; ln(2) = 0.;

    for (int j=0; j<elem.size(); j++) {
      jn += elem.at(j)->GetCurrentBasisVector(i);
      ln += elem.at(j)->GetLoopBasisVector(i);
    }

    cout << setw(14) << setprecision(6) << scientific << jn(0)
         << setw(14) << setprecision(6) << scientific << jn(1)
         << setw(14) << setprecision(6) << scientific << jn(2);
    cout << setw(14) << setprecision(6) << scientific << ln(0)
         << setw(14) << setprecision(6) << scientific << ln(1)
         << setw(14) << setprecision(6) << scientific << ln(2);

    for (int j=0; j<elem.size(); j++)
      cout << setw(10) << fixed << elem.at(j)->GetId();

    cout << "\n";
  }
}

vector<int> XMeshHandle2 :: search_node_at_boundary()
{
  vector<int> node_id;
  
  for (int i=0; i<fNode.size(); i++) {
    if (fNode.at(i)->Is_AtBoundary())
      node_id.push_back( fNode.at(i)->GetId() );
  }

  return node_id;
}

void XMeshHandle2 :: GetListOfBoundary(VectorXi& nodeId, VectorXi& bcLabel)
{
  // search for the node at boundary
  vector<int> nodes = search_node_at_boundary();

  // initialize the vectors
  nodeId = VectorXi :: Zero( nodes.size() );
  bcLabel= VectorXi :: Zero( nodes.size() );

  // fill the vectors
  for (int i=0; i<nodes.size(); i++) {
    nodeId (i) = nodes.at(i);
    bcLabel(i) = fNode.at(nodes.at(i))->GetBoundaryId();
  }
}

vector<int> XMeshHandle2 :: GetBoundaryLabels()
{
  bool is_exist;
  VectorXi nodeId, bcLabel;
  vector<int> nodes;

  // search for the node id and boundary label at bc
  GetListOfBoundary( nodeId, bcLabel );

  // push back the labels
  for (int i=0; i<bcLabel.size(); i++) {
    is_exist = false;

    for (int j=0; j<nodes.size(); j++) {
      if (nodes.at(j)==bcLabel(i)) {
        is_exist = true;
        break;
      }
    }

    if (!is_exist)
      nodes.push_back( bcLabel(i) );
  }

  return nodes;
}

vector<int> XMeshHandle2 :: GetNodeIdAtBoundary(const int label)
{
  vector<int> nodes;
  VectorXi nodeId, bcLabel;

  // search for the node id and boundary label at bc
  GetListOfBoundary( nodeId, bcLabel );

  // push back the nodes
  for (int i=0; i<nodeId.size(); i++) {
    if (label==bcLabel(i))
      nodes.push_back( nodeId(i) );
  }

  return nodes;
}
