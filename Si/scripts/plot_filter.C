void plot_filter(){
TF1 *filter_wiener = new TF1("wiener","exp(-[0]*pow(x,[1]))",0,5);
  filter_wiener->SetParameter(0,1588.);
  filter_wiener->SetParameter(1,15.80);

  filter_wiener->Draw();

}
