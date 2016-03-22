void Spectrometer::countBCDet_BS(Medium &Osrodek, QString DataName, QProgressBar *Progress, int factor, double kx, double ky) {
    //lsfgerhla

    double a=Osrodek.itsBasis.getLatticeConstant();
    double thickness=Osrodek.itsStructure.getThickness();

    int RecVec=itsRecVectors;
    double wi=itsFrequencies[0];
    double w_prec=itsFrequencies[1];
    double wf=itsFrequencies[2];

    int const half=3*(2*RecVec+1)*(2*RecVec+1); //połowa wymiaru macierzy
    int const dimension=6*(2*RecVec+1)*(2*RecVec+1); //wymiar macierzy do zagadnienia własnego
    int const EigValNumber=9*(2*RecVec+1)*(2*RecVec+1); //ilość wybieranych wartości własnych
    int EigValStrNumber=6*(2*RecVec+1)*(2*RecVec+1);

    //plik z danymi
    std::ofstream plik;
    DataName="results/"+ DataName + ".dat";
    QByteArray   bytes  = DataName.toAscii();
    const char * CDataName = bytes.data();
    plik.open(CDataName);

    //inicjalizacje wektorów i macierzy   

    gsl_matrix *gammaA=gsl_matrix_calloc(dimension, dimension);
    gsl_matrix *gammaB=gsl_matrix_calloc(dimension, dimension);

    gsl_matrix *gammaC=gsl_matrix_calloc(dimension, dimension);
    gsl_matrix *gammaD=gsl_matrix_calloc(dimension, dimension);

    gsl_eigen_genv_workspace *wspce=gsl_eigen_genv_alloc(dimension);
    gsl_eigen_genv_workspace *wspce2=gsl_eigen_genv_alloc(dimension);

    gsl_vector_complex *StrAlpha =gsl_vector_complex_alloc(dimension);
    gsl_vector *StrBeta = gsl_vector_alloc(dimension);
    gsl_matrix_complex *StrEigenVec=gsl_matrix_complex_calloc(dimension, dimension);

    gsl_vector_complex *BAlpha =gsl_vector_complex_alloc(dimension);
    gsl_vector *BBeta = gsl_vector_alloc(dimension);
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
               gsl_matrix_set_all(gammaA, 0);
               gsl_matrix_set_all(gammaB, 0);
               gsl_matrix_set_all(gammaC, 0);
               gsl_matrix_set_all(gammaD, 0);
               gsl_vector_complex_set_all(StrAlpha, gsl_complex_rect (0,0));
               gsl_vector_set_all(BBeta, 0);
               gsl_vector_complex_set_all(BAlpha, gsl_complex_rect (0,0));
               gsl_vector_set_all(StrBeta, 0);
               gsl_matrix_complex_set_all(BasisEigenVec, gsl_complex_rect (0,0));
               gsl_matrix_complex_set_all(StrEigenVec, gsl_complex_rect (0,0));
               gsl_matrix_complex_set_all(ChosenVectors, gsl_complex_rect (0,0));
               gsl_vector_complex_set_all(ChosenValues, gsl_complex_rect (0,0));

               //gammaA,B dla struktury
                           //gammaA,B dla struktury
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
                               itsRecStructureSubstance[S].getElasticity(Elasticity);
                               double Density=itsRecStructureSubstance[S].getDensity();
                               double gx=2*M_PI*Nx/a;
                               double gy=2*M_PI*Ny/a;
                               double gx_prim=2*M_PI*Nx_prim/a;
                               double gy_prim=2*M_PI*Ny_prim/a;

                               gsl_matrix_set(gammaA, i, j, Elasticity[3][3]);
                               gsl_matrix_set(gammaA, i+1, j+1, Elasticity[3][3]);
                               gsl_matrix_set(gammaA, i+2, j+2, Elasticity[0][0]);

                               gsl_matrix_set(gammaB, i+2, j, -Elasticity[0][1]*(kx+gx_prim)-Elasticity[3][3]*(kx+gx));
                               gsl_matrix_set(gammaB, i+2, j+1, -Elasticity[0][1]*(ky+gy_prim)-Elasticity[3][3]*(ky+gy));
                               gsl_matrix_set(gammaB, i, j+2, -Elasticity[0][1]*(kx+gx)-Elasticity[3][3]*(kx+gx_prim));
                               gsl_matrix_set(gammaB, i+1, j+2, -Elasticity[0][1]*(ky+gy)-Elasticity[3][3]*(ky+gy_prim));

                               gsl_matrix_set(gammaB, i, j+half, -Elasticity[0][0]*(kx+gx)*(kx+gx_prim)-Elasticity[3][3]*(ky+gy)*(ky+gy_prim)+Density*w*w);
                               gsl_matrix_set(gammaB, i+1, j+half, -Elasticity[0][1]*(kx+gx_prim)*(ky+gy)-Elasticity[3][3]*(ky+gy_prim)*(kx+gx));
                               gsl_matrix_set(gammaB, i, j+half+1, -Elasticity[0][1]*(ky+gy_prim)*(kx+gx)-Elasticity[3][3]*(kx+gx_prim)*(ky+gy));
                               gsl_matrix_set(gammaB, i+1, j+half+1, -Elasticity[0][0]*(ky+gy_prim)*(ky+gy)-Elasticity[3][3]*(kx+gx_prim)*(kx+gx)+Density*w*w);
                               gsl_matrix_set(gammaB, i+2, j+half+2, -Elasticity[3][3]*(ky+gy_prim)*(ky+gy)-Elasticity[3][3]*(kx+gx_prim)*(kx+gx)+Density*w*w);

                               if (i==j) {
                                   gsl_matrix_set(gammaA, i+half, j+half, 1);
                                   gsl_matrix_set(gammaA, i+half+1, j+half+1, 1);
                                   gsl_matrix_set(gammaA, i+half+2, j+half+2, 1);

                                   gsl_matrix_set(gammaB, i+half, j, 1);
                                   gsl_matrix_set(gammaB, i+half+1, j+1, 1);
                                   gsl_matrix_set(gammaB, i+half+2, j+2, 1);
                               }


                           }
                       }
                   }
               }

               //rozwiazanie zagadnienienia własnego
               gsl_eigen_genv(gammaB, gammaA, StrAlpha, StrBeta, StrEigenVec, wspce);


               //gammaC,D dla Podłoża
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
               gsl_eigen_genv(gammaD, gammaC, BAlpha, BBeta, BasisEigenVec, wspce2);

               double imagL, realL; //części Re i Im wartości własnych

               int n=0;
               for (int i = 0; i<dimension; i++)
               {
                   //przepisanie wartości i wektorów własnych struktury do macierzy Chosen*
                   gsl_complex StrValue;
                   StrValue= gsl_complex_div_real(gsl_vector_complex_get(StrAlpha, i), gsl_vector_get(StrBeta,i));
                   gsl_vector_complex_set(ChosenValues, i, StrValue);
                   for (int j = 0, m=0; j < half; j++, m++)
                   {
                       gsl_matrix_complex_set(ChosenVectors, m, i, gsl_matrix_complex_get(StrEigenVec, j, i));
                   }

                   //wybieranie odpowiednich wartości i wektorów własnych dla podłoża i przepisanie do macierzy Chosen*
                   gsl_complex BValue;
                   BValue= gsl_complex_div_real(gsl_vector_complex_get(BAlpha, i), gsl_vector_get(BBeta,i));
                   imagL=GSL_IMAG(BValue);
                   realL=GSL_REAL(BValue);

                   if (imagL > 0.00001 && n+EigValStrNumber<EigValNumber) //warunek na wartości własne && żeby nie było ich więcej niż połowa
                   {
                       gsl_vector_complex_set(ChosenValues, n+EigValStrNumber, BValue); //wybranie wartości własnej

                       for (int j = 0, m=0; j < half; j++, m++)
                       {
                          gsl_matrix_complex_set(ChosenVectors, m, n+EigValStrNumber, gsl_complex_mul_real(gsl_matrix_complex_get(BasisEigenVec, j, i), -1));  //wybranie drugiej połowy wektora własnego
                       }

                       n++;
                   }

               }

               if (n+EigValStrNumber<EigValNumber)
               {
                   for (int i = 0; i<dimension; i++)
                   {
                       gsl_complex BValue;
                       BValue= gsl_complex_div_real(gsl_vector_complex_get(BAlpha, i), gsl_vector_get(BBeta,i));
                       imagL=GSL_IMAG(BValue);
                       realL=GSL_REAL(BValue);

                       if (imagL < 0.00001 && imagL > -0.00001 && realL < -0.00001 && n+EigValStrNumber<EigValNumber) //warunek na wartości własne && żeby nie było ich więcej niż połowa
                       {
                           gsl_vector_complex_set(ChosenValues, n+EigValStrNumber, BValue); //wybranie wartości własnej

                           for (int j = 0, m=0; j < half; j++, m++)
                           {
                              gsl_matrix_complex_set(ChosenVectors, m, n+EigValStrNumber, gsl_complex_mul_real(gsl_matrix_complex_get(BasisEigenVec, j, i), -1));  //wybranie drugiej połowy wektora własnego
                           }

                           n++;
                       }
                   }
               }

//               gsl_vector_complex *normal = gsl_vector_complex_calloc(3);
//               //unormowanie wektorów ChosenVectors
//               for (int i = 0; i < half; i=i+3) {
//                   for(int j = 0; j < EigValNumber; j++) {

//                       gsl_vector_complex_set(normal, 0, gsl_matrix_complex_get(ChosenVectors, i, j));
//                       gsl_vector_complex_set(normal, 1, gsl_matrix_complex_get(ChosenVectors, i+1, j));
//                       gsl_vector_complex_set(normal, 2, gsl_matrix_complex_get(ChosenVectors, i+2, j));

//                       double norm = 1/gsl_blas_dznrm2(normal);
//                       gsl_vector_complex_scale(normal, gsl_complex_rect(norm, 0));

//                       gsl_matrix_complex_set(ChosenVectors, i, j, gsl_vector_complex_get(normal, 0));
//                       gsl_matrix_complex_set(ChosenVectors, i+1, j, gsl_vector_complex_get(normal, 1));
//                       gsl_matrix_complex_set(ChosenVectors, i+2, j, gsl_vector_complex_get(normal, 2));

//                   }
//               }

               //wyznacznik warunków brzegowych - konstrukcja
               /*
               S, S' - numerujš transformaty tensora sprężystoci
               i - numeruje wektory własne A w pętli dla G
               j - numeruje warunki brzegowe dla kolejnych wektorow odwrotnych G
               k - numeruje wektory własne A w pętli dla G'
               L - numeruje wartoci własne
               */
               for (int Nx=-RecVec, S=0, S_prim=0, i=0, j=0; Nx <= RecVec; Nx++) {
                   for (int Ny=-RecVec; Ny <= RecVec; Ny++, j=j+9, i=i+3) {

                       for (int L=0; L < EigValNumber; L++) {

                           S_prim = S;
                           for (int Nx_prim=-RecVec, k=0; Nx_prim <= RecVec; Nx_prim++) {
                               for (int Ny_prim=-RecVec; Ny_prim <= RecVec; Ny_prim++, S_prim++, k=k+3) {

                                   double StrElasticity[6][6];
                                   itsRecStructureSubstance[S_prim].getElasticity(StrElasticity);
                                   double BasisElasticity[6][6];
                                   itsRecBasisSubstance[S_prim].getElasticity(BasisElasticity);

                                   double gx_prim=2*M_PI*Nx_prim/a;
                                   double gy_prim=2*M_PI*Ny_prim/a;

                                   if (L < EigValStrNumber)
                                   {

                                       //eksponens
                                       gsl_complex exponent = gsl_complex_polar(exp(GSL_IMAG(gsl_vector_complex_get(ChosenValues, L))*thickness), -1*GSL_REAL(gsl_vector_complex_get(ChosenValues, L))*thickness);

                                       //warunki zerowania się naprężenia na powierzchni
                                       gsl_complex w1 = gsl_complex_mul_real(exponent, StrElasticity[3][3]);
                                       gsl_complex w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+2, L), (kx+gx_prim));
                                       gsl_complex w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k, L));
                                       gsl_complex BCjL = gsl_complex_add(gsl_complex_mul(gsl_complex_add(w2, w3), w1), gsl_matrix_complex_get(Boundaries, j, L));
                                       gsl_matrix_complex_set(Boundaries, j, L, BCjL);

                                       w1 = gsl_complex_mul_real(exponent, StrElasticity[3][3]);
                                       w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+2, L), (ky+gy_prim));
                                       w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k+1, L));
                                       BCjL = gsl_complex_add(gsl_complex_mul(gsl_complex_add(w2, w3), w1), gsl_matrix_complex_get(Boundaries, j+1, L));
                                       gsl_matrix_complex_set(Boundaries, j+1, L, BCjL);

                                       w1 = gsl_complex_mul_real(exponent, StrElasticity[0][0]);
                                       gsl_complex w11 = gsl_complex_mul_real(exponent, StrElasticity[0][1]);
                                       w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k, L), (kx+gx_prim));
                                       gsl_complex w22 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+1, L), (ky+gy_prim));
                                       w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k+2, L));
                                       gsl_complex w4 = gsl_complex_add(gsl_complex_mul(gsl_complex_add(w2, w22), w11), gsl_complex_mul(w3, w1));
                                       BCjL = gsl_complex_add(w4, gsl_matrix_complex_get(Boundaries, j+2, L));
                                       gsl_matrix_complex_set(Boundaries, j+2, L, BCjL);

                                       //warunki równości naprężeń na granicy ośrodków - część dla struktury
                                       w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+2, L), (kx+gx_prim));
                                       w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k, L));
                                       BCjL = gsl_complex_add(gsl_complex_mul_real(gsl_complex_add(w2, w3), StrElasticity[3][3]), gsl_matrix_complex_get(Boundaries, j+3, L));
                                       gsl_matrix_complex_set(Boundaries, j+3, L, BCjL);

                                       w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+2, L), (ky+gy_prim));
                                       w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k+1, L));
                                       BCjL = gsl_complex_add(gsl_complex_mul_real(gsl_complex_add(w2, w3), StrElasticity[3][3]), gsl_matrix_complex_get(Boundaries, j+4, L));
                                       gsl_matrix_complex_set(Boundaries, j+4, L, BCjL);

                                       w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k, L), (kx+gx_prim));
                                       w22 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+1, L), (ky+gy_prim));
                                       w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k+2, L));
                                       w4 = gsl_complex_add(gsl_complex_mul_real(gsl_complex_add(w2, w22), StrElasticity[0][1]), gsl_complex_mul_real(w3, StrElasticity[0][0]));
                                       BCjL = gsl_complex_add(w4, gsl_matrix_complex_get(Boundaries, j+5, L));
                                       gsl_matrix_complex_set(Boundaries, j+5, L, BCjL);

                                   } else {

                                       //warunki równości naprężeń na granicy ośrodków - część dla podłoża
                                       gsl_complex w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+2, L), (kx+gx_prim));
                                       gsl_complex w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k, L));
                                       gsl_complex BCjL = gsl_complex_add(gsl_complex_mul_real(gsl_complex_add(w2, w3), BasisElasticity[3][3]), gsl_matrix_complex_get(Boundaries, j+3, L));
                                       gsl_matrix_complex_set(Boundaries, j+3, L, BCjL);

                                       w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+2, L), (ky+gy_prim));
                                       w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k+1, L));
                                       BCjL = gsl_complex_add(gsl_complex_mul_real(gsl_complex_add(w2, w3), BasisElasticity[3][3]), gsl_matrix_complex_get(Boundaries, j+4, L));
                                       gsl_matrix_complex_set(Boundaries, j+4, L, BCjL);

                                       w2 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k, L), (kx+gx_prim));
                                       gsl_complex w22 = gsl_complex_mul_real(gsl_matrix_complex_get(ChosenVectors, k+1, L), (ky+gy_prim));
                                       w3 = gsl_complex_mul(gsl_vector_complex_get(ChosenValues, L), gsl_matrix_complex_get(ChosenVectors, k+2, L));
                                       gsl_complex w4 = gsl_complex_add(gsl_complex_mul_real(gsl_complex_add(w2, w22), BasisElasticity[0][1]), gsl_complex_mul_real(w3, BasisElasticity[0][0]));
                                       BCjL = gsl_complex_add(w4, gsl_matrix_complex_get(Boundaries, j+5, L));
                                       gsl_matrix_complex_set(Boundaries, j+5, L, BCjL);

                                   }
                               }
                           }

                           // warunek równości wychyleń na granicy ośrodków
                           gsl_matrix_complex_set(Boundaries, j+6, L, gsl_matrix_complex_get(ChosenVectors, i, L));
                           gsl_matrix_complex_set(Boundaries, j+7, L, gsl_matrix_complex_get(ChosenVectors, i+1, L));
                           gsl_matrix_complex_set(Boundaries, j+8, L, gsl_matrix_complex_get(ChosenVectors, i+2, L));

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

               //usuwanie NaN
               if(DetVal != DetVal) DetVal = 0;

               //usuwanie inf
               if(isanyinf(DetVal)) DetVal = 0;

               //zapisanie wartości do pliku
               plik << w << "\t" << DetVal << "\n";

            }

    plik.close();
    gsl_matrix_free(gammaA);
    gsl_matrix_free(gammaB);
    gsl_vector_free(BBeta);
    gsl_vector_free(StrBeta);
    gsl_matrix_free(gammaC);
    gsl_matrix_free(gammaD);
    gsl_vector_complex_free(StrAlpha);
    gsl_vector_complex_free(BAlpha);
    gsl_eigen_genv_free(wspce);
    gsl_eigen_genv_free(wspce2);
    gsl_matrix_complex_free(Boundaries);
    gsl_matrix_complex_free(ChosenVectors);
    gsl_vector_complex_free(ChosenValues);

}


