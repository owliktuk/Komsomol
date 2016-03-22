void Spectrometer::drawStructure(Medium &Osrodek, QTextBrowser *Browser) {

    double a=Osrodek.itsBasis.getLatticeConstant();

    int RecVec=itsRecVectors;
    double krok=a/50;

    double Density;
    double gx;
    double gy;
    double gx_prim;
    double gy_prim;

    gsl_complex rho;

    //plik z danymi
    std::ofstream plik;
    plik.open("structure.dat");

    Browser->insertPlainText("Counting...");
    QApplication::processEvents();

    for (double x=0; x <= a; x=x+krok) {
        for(double y=0; y <= a; y=y+krok) {

            rho=gsl_complex_rect(0, 0);

            for(int Nx=-RecVec, S=0; Nx<=RecVec; Nx++) {
                for(int Ny=-RecVec; Ny<=RecVec; Ny++) {
                    for(int Nx_prim=-RecVec; Nx_prim<=RecVec; Nx_prim++) {
                        for(int Ny_prim=-RecVec; Ny_prim<=RecVec; Ny_prim++, S++) {


                            double Elasticity[6][6];
                            itsRecBasisSubstance[S].getElasticity(Elasticity);
                            //Density=itsRecBasisSubstance[S].getDensity();
                            gx=2*M_PI*Nx/a;
                            gy=2*M_PI*Ny/a;
                            gx_prim=2*M_PI*Nx_prim/a;
                            gy_prim=2*M_PI*Ny_prim/a;

                            rho = gsl_complex_add(rho, gsl_complex_polar(Elasticity[0][0], (gx-gx_prim)*x+(gy-gy_prim)*y));

                        }
                    }
                }
            }

            //zapisanie wartości do pliku
            plik << x << "\t" << y << "\t" << GSL_REAL(rho)/pow(2*RecVec+1, 2) << "\n";

        }

        plik << "\n";
    }

    plik.close();

}
