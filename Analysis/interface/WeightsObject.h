#ifndef WEIGHTSOBJECT_H
#define WEIGHTSOBJECT_H


class WeightVariables{
public:
  float wgt_central;
  float wgt_muR2;
  float wgt_muR0p5;
  float wgt_muF2;
  float wgt_muF0p5;
  float wgt_PDFVariationUp;
  float wgt_PDFVariationDn;
  float wgt_AsMZUp;
  float wgt_AsMZDn;
  float wgt_PSUp;
  float wgt_PSDn;

  WeightVariables();
  WeightVariables(WeightVariables const& other);
  WeightVariables& operator=(const WeightVariables& other);

  void swap(WeightVariables& other);

};

class WeightsObject{
public:
  WeightVariables extras;

  WeightsObject();
  WeightsObject(const WeightsObject& other);
  WeightsObject& operator=(const WeightsObject& other);
  ~WeightsObject();

  void swap(WeightsObject& other);

};

#endif
