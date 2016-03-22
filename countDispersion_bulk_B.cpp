void Spectrometer::countDispersion_bulk_B(Medium &Osrodek, QString DataName, QProgressBar *Progress) {
    //lsfgerhla

    double a=Osrodek.itsBasis.getLatticeConstant();

    int RecVec=itsRecVectors;
    int k_prec=itsK_Precision;
    double wi=itsFrequencies[0];
    double w_prec=itsFrequencies[1];
    double wf=itsFrequencies[2];

    int dimension=3*(2*RecVec+1)*(2*RecVec+1);

    std::ofstream plik;
    DataName="results/"+ DataName + ".dat";
    QByteArray   bytes  = DataName.toAscii();
    const char * CDataName = bytes.data();
    plik.open(CDataName);

    //inicjalizacje wektorów i macierzy

    gsl_matrix *gamma=gsl_matrix_calloc(dimension, dimension);

    gsl_permutation *permutation= gsl_permutation_alloc(dimension);
    int signum;
    gsl_eigen_nonsymmv_workspace *wspce=gsl_eigen_nonsymmv_alloc(dimension);

    double kx, ky, krokx, kroky, boundary_x, boundary_y;
    double k_zred, k_zred0;

    double krok = M_PI/(k_prec*a);

    for (int droga=0; droga<3; droga++) {
        if (droga==0)           //droga M->Gamma
        {
            kx=-M_PI/a;
            krokx=krok;
            boundary_x=0;
            ky=-M_PI/a;
            kroky=krok;
            boundary_y=0;

            k_zred0=-1*sqrt(pow(kx/(2*M_PI/a), 2)+pow(ky/(2*M_PI/a), 2));
        }
        else if (droga==1)
        {
            kx=0;               //droga Gamma->X
            krokx=krok;
            boundary_x=M_PI/a;
            ky=0;
            kroky=0;
            boundary_y=0;

            k_zred0=sqrt(2)/2;
        }
        else if (droga==2)
        {
            kx=M_PI/a;          //Droga X->M
            krokx=0;
            boundary_x=M_PI/a;
            ky=0;
            kroky=krok;
            boundary_y=M_PI/a;

            k_zred0=sqrt(2)/2;
        }

       //petla dla wektorów falowych
       for (; kx <= boundary_x && ky <= boundary_y; kx=kx+krokx, ky=ky+kroky)
       {
           if (droga==0) {
               k_zred = abs(k_zred0 + sqrt( pow(kx/(2*M_PI/a), 2)+pow(ky/(2*M_PI/a), 2)));
           } else {
               k_zred = k_zred0 + kx/(2*M_PI/a) + ky/(2*M_PI/a);
           }


           int postep=int(100*k_zred/1.7);
           Progress->setValue(postep);
           Progress->update();
           QApplication::processEvents();

            //pętla dla częstości w
           for (double w=wi; w<wf; w=w+w_prec)
            {

                //gamma dla Podłoża
                /*
                S - numeruje transformaty tensora sprężystoci i gestosci
                i - numeruje wiersze macierzy
                j - numeruje kolumny macierzy
                */
                for(int Nx=-RecVec, i=0, S=0; Nx<=RecVec; Nx++) {
                    for(int Ny=-RecVec; Ny<=RecVec; Ny++, i=i+3) {
                        for(int Nx_prim=-RecVec, j=0; Nx_prim<=RecVec; Nx_prim++) {
                            for(int Ny_prim=-RecVec; Ny_prim<=RecVec; Ny_prim++, j=j+3, S++) {

                                double Elasticity[6][6];
                                itsRecBasisSubstance[S].getElasticity(Elasticity);
                                double Density=itsRecBasisSubstance[S].getDensity();
                                double gx=2*M_PI*Nx/a;
                                double gy=2*M_PI*Ny/a;
                                double gx_prim=2*M_PI*Nx_prim/a;
                                double gy_prim=2*M_PI*Ny_prim/a;

                                gsl_matrix_set(gamma, i, j, Elasticity[0][0]*(kx+gx)*(kx+gx_prim)+Elasticity[3][3]*(ky+gy)*(ky+gy_prim)-Density*w*w);
                                gsl_matrix_set(gamma, i+1, j, Elasticity[0][1]*(kx+gx_prim)*(ky+gy)+Elasticity[3][3]*(ky+gy_prim)*(kx+gx));
                                gsl_matrix_set(gamma, i, j+1, Elasticity[0][1]*(ky+gy_prim)*(kx+gx)+Elasticity[3][3]*(kx+gx_prim)*(ky+gy));
                                gsl_matrix_set(gamma, i+1, j+1, Elasticity[0][0]*(ky+gy_prim)*(ky+gy)+Elasticity[3][3]*(kx+gx_prim)*(kx+gx)-Density*w*w);
                                gsl_matrix_set(gamma, i+2, j+2, Elasticity[3][3]*(ky+gy_prim)*(ky+gy)+Elasticity[3][3]*(kx+gx_prim)*(kx+gx)-Density*w*w);

                                //if(Browser->insertPlainText("Counting...");
                            }
                        }
                    }
                }



    //            for (int c=0; c < 54; c++) {
    //                //for (int d=0; d < EigValNumber; d++) {
    //                    Check[c]=GSL_REAL(gsl_vector_complex_get(StrEigenVal, c));
    //               // }
    //            }


                //obliczenie wyznacznika z gamma
                gsl_linalg_LU_decomp(gamma, permutation, &signum);
                double DetVal = gsl_linalg_LU_lndet(gamma);

                //zapisanie wartości do pliku
                plik << k_zred << "\t" << w << "\t" << DetVal << "\n";

            }
            plik << "\n";
        }

    }

    plik.close();
    gsl_matrix_free(gamma);
    gsl_eigen_nonsymmv_free(wspce);

}


