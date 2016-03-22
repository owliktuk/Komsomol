void Spectrometer::countBCDet_B(Medium &Osrodek, QString DataName, QProgressBar *Progress, int factor, double kx, double ky) {
    //lsfgerhla

    double a=Osrodek.itsBasis.getLatticeConstant();

    int RecVec=itsRecVectors;
    double wi=itsFrequencies[0];
    double w_prec=itsFrequencies[1];
    double wf=itsFrequencies[2];

    int const half=3*(2*RecVec+1)*(2*RecVec+1); //połowa wymiaru macierzy
    int const dimension=6*(2*RecVec+1)*(2*RecVec+1); //wymiar macierzy do zagadnienia własnego
    int const EigValNumber=3*(2*RecVec+1)*(2*RecVec+1); //ilość wybieranych wartości własnych

    //plik z danymi
    std::ofstream plik;
    DataName="results/"+ DataName + ".dat";
    QByteArray   bytes  = DataName.toAscii();
    const char * CDataName = bytes.data();
    plik.open(CDataName);

    //inicjalizacje wektorów i macierzy
    gsl_matrix *gammaC=gsl_matrix_calloc(dimension, dimension);
    gsl_matrix *gammaD=gsl_matrix_calloc(dimension, dimension);

    gsl_eigen_genv_workspace * wspce=gsl_eigen_genv_alloc(dimension);

    gsl_vector_complex *alpha = gsl_vector_complex_calloc(dimension);
    gsl_vector *beta = gsl_vector_calloc(dimension);

    gsl_matrix_complex *BasisEigenVec=gsl_matrix_complex_calloc(dimension, dimension);

    gsl_matrix_complex *ChosenVectors = gsl_matrix_complex_calloc(half, EigValNumber);
    gsl_vector_complex *ChosenValues = gsl_vector_complex_calloc(EigValNumber);
    gsl_matrix_complex *Boundaries=gsl_matrix_complex_calloc(EigValNumber, EigValNumber);


            //pętla dla częstości w
           for (double w=wi; w<wf; w=w+w_prec)
            {
               //pasek postepu
               int postep=int(100*(w-wi)/(wf-wi));
               Progress->setValue(postep);
               Progress->update();
               QApplication::processEvents();

                gsl_matrix_complex_set_all(Boundaries, gsl_complex_rect (0,0)); //ustawienie wartosci wyznacznika na 0
                gsl_matrix_set_all(gammaC, 0);
                gsl_matrix_set_all(gammaD, 0);
                gsl_vector_complex_set_all(alpha, gsl_complex_rect (0,0));
                gsl_vector_set_all(beta, 0);
                gsl_matrix_complex_set_all(BasisEigenVec, gsl_complex_rect (0,0));
                gsl_matrix_complex_set_all(ChosenVectors, gsl_complex_rect (0,0));
                gsl_vector_complex_set_all(ChosenValues, gsl_complex_rect (0,0));

                //gammaC,D dla Podłoża
                /*
                S - numeruje transformaty tensora sprężystoci i gestosci
                i - numeruje wiersze macierzy
                j - numeruje kolumny macierzy
                half - druga polowa macierzy
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

                                gsl_matrix_set(gammaC, i, j, Elasticity[3][3]);
                                gsl_matrix_set(gammaC, i+1, j+1, Elasticity[3][3]);
                                gsl_matrix_set(gammaC, i+2, j+2, Elasticity[0][0]);

                                gsl_matrix_set(gammaD, i+2, j, -Elasticity[0][1]*(kx+gx_prim)-Elasticity[3][3]*(kx+gx));
                                gsl_matrix_set(gammaD, i+2, j+1, -Elasticity[0][1]*(ky+gy_prim)-Elasticity[3][3]*(ky+gy));
                                gsl_matrix_set(gammaD, i, j+2, -Elasticity[0][1]*(kx+gx)-Elasticity[3][3]*(kx+gx_prim));
                                gsl_matrix_set(gammaD, i+1, j+2, -Elasticity[0][1]*(ky+gy)-Elasticity[3][3]*(ky+gy_prim));

                                gsl_matrix_set(gammaD, i, j+half, -Elasticity[0][0]*(kx+gx)*(kx+gx_prim)-Elasticity[3][3]*(ky+gy)*(ky+gy_prim)+Density*w*w);
                                gsl_matrix_set(gammaD, i+1, j+half, -Elasticity[0][1]*(kx+gx_prim)*(ky+gy)-Elasticity[3][3]*(ky+gy_prim)*(kx+gx));
                                gsl_matrix_set(gammaD, i, j+half+1, -Elasticity[0][1]*(ky+gy_prim)*(kx+gx)-Elasticity[3][3]*(kx+gx_prim)*(ky+gy));
                                gsl_matrix_set(gammaD, i+1, j+half+1, -Elasticity[0][0]*(ky+gy_prim)*(ky+gy)-Elasticity[3][3]*(kx+gx_prim)*(kx+gx)+Density*w*w);
                                gsl_matrix_set(gammaD, i+2, j+half+2, -Elasticity[3][3]*(ky+gy_prim)*(ky+gy)-Elasticity[3][3]*(kx+gx_prim)*(kx+gx)+Density*w*w);

                                if (i==j) {
                                    gsl_matrix_set(gammaC, i+half, j+half, 1);
                                    gsl_matrix_set(gammaC, i+half+1, j+half+1, 1);
                                    gsl_matrix_set(gammaC, i+half+2, j+half+2, 1);

                                    gsl_matrix_set(gammaD, i+half, j, 1);
                                    gsl_matrix_set(gammaD, i+half+1, j+1, 1);
                                    gsl_matrix_set(gammaD, i+half+2, j+2, 1);
                                }


                            }
                        }
                    }
                }

                //rozwiazanie zagadnienienia własnego
                gsl_eigen_genv(gammaD, gammaC, alpha, beta, BasisEigenVec, wspce);

                double imagL, realL; //części Re i Im wartości własnych

                int n=0;
                for (int i = 0; i<dimension; i++)
                {
                    //wybieranie odpowiednich wartości i wektorów własnych dla podłoża i przepisanie do ChosenValues, ChosenVectors
                    gsl_complex BValue;
                    //BValue= gsl_vector_complex_get(BasisEigenVal, i);
                    BValue= gsl_complex_div_real(gsl_vector_complex_get(alpha, i), gsl_vector_get(beta,i));
                    imagL=GSL_IMAG(BValue);
                    realL=GSL_REAL(BValue);

                    //Browser->insertPlainText("\n" + QString::number(imagL));


                    if (imagL > 0.00001) //warunek na wartości własne && żeby nie było ich więcej niż połowa
                    {
                        if (n<EigValNumber)
                        {
                            gsl_vector_complex_set(ChosenValues, n, BValue); //wybranie wartości własnej

                            for (int j = half, m=0; j < dimension; j++, m++)
                            {
                               gsl_matrix_complex_set(ChosenVectors, m, n, gsl_matrix_complex_get(BasisEigenVec, j, i));  //wybranie drugiej połowy wektora własnego
                            }
                        }
                        n++;
                    }

                }

                //wybranie rozwiązań imagL==0
                for (int i = 0; i < dimension && n < EigValNumber; i++)
                {
                    gsl_complex BValue;
                    BValue= gsl_complex_div_real(gsl_vector_complex_get(alpha, i), gsl_vector_get(beta,i));
                    imagL=GSL_IMAG(BValue);
                    realL=GSL_REAL(BValue);

                    //Browser->insertPlainText("\n" + QString::number(realL));


                    if (imagL > -0.00001 && imagL < 0.00001 && realL < -0.00001) //warunek na wartości własne && żeby nie było ich więcej niż połowa
                    {
                        gsl_vector_complex_set(ChosenValues, n, BValue); //wybranie wartości własnej

                        //Browser->insertPlainText("s");

                        for (int j = half, m=0; j < dimension; j++, m++)
                        {
                           gsl_matrix_complex_set(ChosenVectors, m, n, gsl_matrix_complex_get(BasisEigenVec, j, i));  //wybranie drugiej połowy wektora własnego
                        }

                        n++;
                    }
                }

//                if(n==EigValNumber)
//                {
//                    Browser->insertPlainText("o");
//                } else {
//                    Browser->insertPlainText("W");
//                }

                //wyznacznik warunków brzegowych - konstrukcja
                /*
                S, S' - numerujš transformaty tensora sprężystoci
                j - numeruje warunki brzegowe dla kolejnych wektorow odwrotnych G
                k - numeruje wektory własne A w pętli dla G'
                L - numeruje wartoci własne
                */

                for (int Nx=-RecVec, S=0, S_prim=0, j=0; Nx <= RecVec; Nx++) {
                    for (int Ny=-RecVec; Ny <= RecVec; Ny++, j=j+3) {

                        for (int L=0; L < EigValNumber; L++) {

                            S_prim = S;
                            for (int Nx_prim=-RecVec, k=0; Nx_prim <= RecVec; Nx_prim++) {
                                for (int Ny_prim=-RecVec; Ny_prim <= RecVec; Ny_prim++, S_prim++, k=k+3) {

                                    double BasisElasticity[6][6];
                                    itsRecBasisSubstance[S_prim].getElasticity(BasisElasticity);


                                    double gx_prim=2*M_PI*Nx_prim/a;
                                    double gy_prim=2*M_PI*Ny_prim/a;


                                    //warunki zerowania naprężenia na powierzchni
                                    gsl_complex w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+2, L), (kx+gx_prim)); //EigenVec(z,L)*(kx+gx_prim)
                                    gsl_complex w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k, L));  // EigenValue(L)*EigenVec(x,L)
                                    gsl_complex BCjL = gsl_complex_add(gsl_complex_mul_real(gsl_complex_add(w2, w3), BasisElasticity[3][3]), gsl_matrix_complex_get(Boundaries, j, L)); // (w2+w3)*C44 + Boundaries(j,L)
                                    gsl_matrix_complex_set(Boundaries, j, L, BCjL);

                                    w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+2, L), (ky+gy_prim)); ////EigenVec(z,L)*(ky+gy_prim)
                                    w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k+1, L)); // EigenValue(L)*EigenVec(y,L)
                                    BCjL = gsl_complex_add(gsl_complex_mul_real(gsl_complex_add(w2, w3), BasisElasticity[3][3]), gsl_matrix_complex_get(Boundaries, j+1, L)); // (w2+w3)*C44 + Boundaries(j+1,L)
                                    gsl_matrix_complex_set(Boundaries, j+1, L, BCjL);

                                    w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k, L), (kx+gx_prim)); //EigenVec(x,L)*(kx+gx_prim)
                                    gsl_complex w22 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+1, L), (ky+gy_prim)); //EigenVec(y,L)*(ky+gy_prim)
                                    w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k+2, L)); // EigenValue(L)*EigenVec(z,L)
                                    gsl_complex w4 = gsl_complex_add(gsl_complex_mul_real(gsl_complex_add(w2, w22), BasisElasticity[0][1]), gsl_complex_mul_real(w3, BasisElasticity[0][0])); // (w2+w22)*C12 + w3*C11
                                    BCjL = gsl_complex_add(w4, gsl_matrix_complex_get(Boundaries, j+2, L)); // w4 + Boundaries(j+2,L)
                                    gsl_matrix_complex_set(Boundaries, j+2, L, BCjL);


                                }
                            }

                        }
                        S=S_prim;
                    }
                }

                //skalowanie macierzy Boundaries
                gsl_complex scale=gsl_complex_rect(pow(10, factor), 0);
                gsl_matrix_complex_scale(Boundaries, scale);

                //obliczenie wyznacznika z Boundaries
                gsl_permutation *Bpermutation = gsl_permutation_alloc(EigValNumber);
                int Bsignum;
                gsl_linalg_complex_LU_decomp(Boundaries, Bpermutation, &Bsignum);
                double DetVal = gsl_linalg_complex_LU_lndet(Boundaries);

                //zapisanie wartości do pliku
                plik <<  w << "\t" << DetVal << "\n";

            }

    plik.close();
    gsl_matrix_free(gammaC);
    gsl_matrix_free(gammaD);
    gsl_eigen_genv_free(wspce);
    gsl_matrix_complex_free(Boundaries);
    gsl_matrix_complex_free(ChosenVectors);
    gsl_vector_complex_free(ChosenValues);

}


