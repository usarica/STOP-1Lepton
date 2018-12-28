// Special PDF cases

struct ExceptionalCases{
  bool specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1;
  bool specialPDF_NNPDF31_NNLO_as_0118_nf_4;
  bool specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1;

  ExceptionalCases() :
    specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1(false),
    specialPDF_NNPDF31_NNLO_as_0118_nf_4(false),
    specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1(false)
  {}
};
ExceptionalCases exceptionalCases;

// madgraph_1000offset special case with POWHEG-like numbering and, unfortunately, NNPDF30_nlo_nf_4_pdfas as the PDF choice
void set_specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1(bool flag){ exceptionalCases.specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1=flag; }
// madgraph_1000offset special case with some numbering differences and MC-type PDF variations
void set_specialPDF_NNPDF31_NNLO_as_0118_nf_4(bool flag){ exceptionalCases.specialPDF_NNPDF31_NNLO_as_0118_nf_4=flag; }
// madgraph_1000offset special case with 45 muR, muF variations with 5 completely different central value choices
void set_specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1(bool flag){ exceptionalCases.specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1=flag; }
