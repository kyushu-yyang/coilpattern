#include <iostream>
#include <iomanip>
#include "XVtkFile.h"
#include "XMeshHandle2.h"

using namespace std;

void test(const char* filename)
{
  XVtkFile* file = new XVtkFile( filename );
  XMeshHandle2* msh = new XMeshHandle2;
  msh->SetVtkFile( file );
  delete file;

  vector<XNode*> node = msh->GetSurroundingNodes(2);
  for (int i=0; i<node.size(); i++) cout << node.at(i)->GetId() << endl;
}

int main(int argc, char** argv)
{
  try {
    test(argv[1]);
  }
  catch (invalid_argument& except) {
    cerr << except.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS; 
}
