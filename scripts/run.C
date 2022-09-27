{
    // gROOT->ProcessLine(".x loadClasses.C" ); -> not working in root 6
    // root -l loadClasses.C run.C

    x = SiProperties();
    x.PrintInfo();
}