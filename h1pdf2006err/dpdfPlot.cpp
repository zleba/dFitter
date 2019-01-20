extern "C" {
      void qcd_2006_(double *Z, double Q2, IFIT,XPQ,F2,FL,C2,CL)
extern "C" {
   //   void diffpdf_(double* xpom, double*  zpom, double*  Q2, double *pdfs);
   void diffpdferr_(double* xpom, double*  zpom, double*  Q2, int* ipdf, int* ierr, double *pdfs, double *tmax);
}
