#include "common_includes.h"


class EventAnalyzer : public FrameworkTreeLooperBase{
protected:

  bool runEvent(FrameworkTree* tree, float const& externalWgt, SimpleEntry& product);

public:
  EventAnalyzer();
  EventAnalyzer(FrameworkTree* inTree);
  EventAnalyzer(std::vector<FrameworkTree*> const& inTreeList);
  EventAnalyzer(FrameworkSet const* inTreeSet);

};

EventAnalyzer::EventAnalyzer() : FrameworkTreeLooperBase() {}
EventAnalyzer::EventAnalyzer(FrameworkTree* inTree) : FrameworkTreeLooperBase(inTree) {}
EventAnalyzer::EventAnalyzer(std::vector<FrameworkTree*> const& inTreeList) : FrameworkTreeLooperBase(inTreeList) {}
EventAnalyzer::EventAnalyzer(FrameworkSet const* inTreeSet) : FrameworkTreeLooperBase(inTreeSet) {}

bool EventAnalyzer::runEvent(FrameworkTree* tree, float const& externalWgt, SimpleEntry& product){
  bool validProducts = (tree!=nullptr);
  if (!validProducts) return validProducts;

  float wgt = externalWgt;

  WeightsHandler* wgtHandler=nullptr;
  for (auto it_ivy=this->externalIvyObjects.begin(); it_ivy!=this->externalIvyObjects.end(); it_ivy++){
    WeightsHandler* tmp_ivy = dynamic_cast<WeightsHandler*>(it_ivy->second);
    if (tmp_ivy){ wgtHandler = tmp_ivy; break; }
  }
  validProducts &= (wgtHandler!=nullptr);
  if (!validProducts){
    MELAerr << "EventAnalyzer::runEvent: Weight handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }

  if (wgtHandler){
    validProducts &= wgtHandler->constructWeights();
    if (!validProducts){
      MELAerr << "EventAnalyzer::runEvent: Weight product could not be constructed (Tree: " << tree->sampleIdentifier << ")." << endl;
      return validProducts;
    }

    auto const* wgtProduct = wgtHandler->getProduct();
    validProducts &= (wgtProduct!=nullptr);
    if (!validProducts){
      MELAerr << "EventAnalyzer::runEvent: Weight product is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
      return validProducts;
    }
    wgt *= wgtProduct->extras.wgt_central;
  }

  product.setNamedVal("weight", wgt);
  if (std::isnan(wgt) || std::isinf(wgt) || wgt==0.f){
    if (wgt!=0.f){
      MELAerr << "EventAnalyzer::runEvent: Invalid weight " << wgt << " is being discarded (Tree " << tree->sampleIdentifier << ")." << endl;
      exit(1);
    }
    validProducts=false;
  }
  if (!validProducts) return validProducts;

  return validProducts;
}


void testWeights(){
  std::string stropts = "indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./ outfile=WZZ_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=100 ismc=true";
  FrameworkOptionParser opts(stropts);
  FrameworkSet theSet(opts, CMS4_EVENTS_TREE_NAME);

  //EventAnalyzer analyzer(theSet);
}