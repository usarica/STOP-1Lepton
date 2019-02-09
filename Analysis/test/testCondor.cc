#include <iostream>
#include "MELAStreamHelpers.hh"
#include "MELAStreamHelpers.hh"
#include "FileTransferHelpers.h"


using namespace std;
using namespace MELAStreamHelpers;
using namespace FileTransferHelpers;


void testCondor(TString strarg){
  MELAout.open("test_output.txt");
  MELAout << strarg << endl;
  MELAout.close();

  // Transfer the file
  InitiateCondorFileTransfer("./", "test_output.txt", "t2.ucsd.edu", "/hadoop/cms/store/user/usarica/STOP_1L/output/testCondor");
}
