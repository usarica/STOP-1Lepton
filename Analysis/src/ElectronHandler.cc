#include "ElectronHandler.h"


using namespace std;


bool ElectronHandler::constructElectrons(){
  if (!currentTree) return false;

  vector<int> const* charge;
  vector<CMSLorentzVector> const* momentum;

  vector<bool> const* conv_vtx_flag;

  vector<int> const* expectedMissingInnerHits;

  vector<float> const* energySC;
  vector<float> const* etaSC;
  vector<float> const* etaSeedSC;
  vector<float> const* rho;
  vector<float> const* sigmaIEtaIEta_full5x5;
  vector<float> const* dEtaIn;
  vector<float> const* dPhiIn;
  vector<float> const* hOverE;
  vector<float> const* ecalEnergy;
  vector<float> const* eOverPIn;
  vector<float> const* dxyPV;
  vector<float> const* dzPV;
  vector<float> const* miniIso_ch;
  vector<float> const* miniIso_nh;
  vector<float> const* miniIso_em;

  if (charge->empty()) return false;

  unsigned int nProducts = charge->size();
  productList.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
    productList.push_back(new ElectronObject(-11*(charge->at(ip)>0 ? 1 : -1), momentum->at(ip)));
    ElectronObject*& obj = productList.back();

    obj->extras.conv_vtx_flag = conv_vtx_flag->at(ip);

    obj->extras.expectedMissingInnerHits = expectedMissingInnerHits->at(ip);

    obj->extras.energySC = energySC->at(ip);
    obj->extras.etaSC = etaSC->at(ip);
    obj->extras.etaSeedSC = etaSeedSC->at(ip);
    obj->extras.rho = rho->at(ip);
    obj->extras.sigmaIEtaIEta_full5x5 = sigmaIEtaIEta_full5x5->at(ip);
    obj->extras.dEtaIn = dEtaIn->at(ip);
    obj->extras.dPhiIn = dPhiIn->at(ip);
    obj->extras.hOverE = hOverE->at(ip);
    obj->extras.ecalEnergy = ecalEnergy->at(ip);
    obj->extras.eOverPIn = eOverPIn->at(ip);
    obj->extras.dxyPV = dxyPV->at(ip);
    obj->extras.dzPV = dzPV->at(ip);
    obj->extras.miniIso_ch = miniIso_ch->at(ip);
    obj->extras.miniIso_nh = miniIso_nh->at(ip);
    obj->extras.miniIso_em = miniIso_em->at(ip);
  }

  return true;
}

bool ElectronHandler::wrapTree(BaseTree* tree){
  if (!IvyBase::wrapTree(tree)) return false;

  return true;
}

void ElectronHandler::bookBranches(BaseTree* tree){

}

