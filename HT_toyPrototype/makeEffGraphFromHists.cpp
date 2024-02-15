//_________________________________________________
//
// function for binomial error
//

void BinomialError(Double_t n,Double_t N,Double_t* p,Double_t* err_low,Double_t* err_high){

  *p = n/N;
  *err_high = (sqrt(n + 0.25 - n*n/N) + 0.5 - n/N)/(N + 1.);
  *err_low  = (sqrt(n + 0.25 - n*n/N) - 0.5 + n/N)/(N + 1.);
  return;
}


/******************************************
Divides two histograms and makes graph with asymmetric
binomial Punzi errors.
Histograms are taken from root files
******************************************/


TGraph* makeEffGraphFromHists(TH1D * num_hist, TH1D * den_hist) {


  // kludge: allocate mem for 1000 bins maximum
  // --- correct size should be dynamically allocated ---

  Double_t eff[1000];
  Double_t erry_high[1000];
  Double_t erry_low[1000];
  Double_t errx_high[1000];
  Double_t errx_low[1000];
  Double_t x_value[1000];

  
  // find number of bins to plot

  Int_t nbins = num_hist->GetNbinsX();
  if(den_hist->GetNbinsX() < nbins) nbins = den_hist->GetNbinsX();// minimum nbins
  if(nbins > 1000) nbins=1000;// max 1000 bins 

  // fill arrays with eff and errors

   for(Int_t ibin=0;ibin<nbins;++ibin){
     
     Double_t x1 = num_hist->GetBinContent(ibin+1);// skip undf & ovf
     Double_t x2 = den_hist->GetBinContent(ibin+1);

     x_value[ibin]= num_hist->GetBinCenter(ibin+1);
     errx_high[ibin]=num_hist->GetBinWidth(ibin+1)/2.;
     errx_low[ibin]=num_hist->GetBinWidth(ibin+1)/2.;

     if(x2 != 0.){
       // estimate efficiency and asymmetric binomial errors
       BinomialError(x1,x2,&eff[ibin],&erry_low[ibin],&erry_high[ibin]);
     }
     else{
       // set eff to 1/2 +/- 1/2
      eff[ibin]=0.5;
      erry_high[ibin]=0.5;
      erry_low[ibin]=0.5;
     }
   }

 // create graph of efficiency

   TGraph* gr  = new TGraphAsymmErrors(nbins,x_value,eff,errx_low,errx_high,erry_low,erry_high);

   return gr;
}