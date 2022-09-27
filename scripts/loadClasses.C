{
    TString include = ".include ";
    TString load = ".L ";
    TString pwd = gSystem->pwd();
    pwd = pwd + '/';

    TString prefix = "../Si";
    gROOT->ProcessLine( include + pwd + prefix );
    gROOT->ProcessLine( load + prefix + "/SiProperties.cc+" );

}