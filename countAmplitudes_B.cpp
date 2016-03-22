void Spectrometer::countAmplitudes_B(Medium &Osrodek, QString DataName, QProgressBar *Progress, double kx, double ky, double w, int polarisation, double x_lenght, double z_lenght, double precision) {


    double a=Osrodek.itsBasis.getLatticeConstant();

    int RecVec=itsRecVectors;

    int half=3*(2*RecVec+1)*(2*RecVec+1); //połowa wymiaru macierzy
    int dimension=6*(2*RecVec+1)*(2*RecVec+1); //wymiar macierzy do zagadnienia własnego
    int EigValNumber=3*(2*RecVec+1)*(2*RecVec+1); //ilość wybieranych wartości własnych

    //plik z danymi
    std::ofstream plik;
    DataName="results/amplitudes/"+ DataName + ".dat";
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
        BValue= gsl_complex_div_real(gsl_vector_complex_get(alpha, i), gsl_vector_get(beta,i));
        imagL=GSL_IMAG(BValue);
        realL=GSL_REAL(BValue);



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



        if (imagL > -0.00001 && imagL < 0.00001 && realL < -0.00001) //warunek na wartości własne && żeby nie było ich więcej niż połowa
        {
            gsl_vector_complex_set(ChosenValues, n, BValue); //wybranie wartości własnej

            for (int j = half, m=0; j < dimension; j++, m++)
            {
               gsl_matrix_complex_set(ChosenVectors, m, n, gsl_matrix_complex_get(BasisEigenVec, j, i));  //wybranie drugiej połowy wektora własnego
            }

            n++;
        }
    }

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

    //rozwiązanie układu równań
    complex <double> Boundaries_svd[EigValNumber*EigValNumber];
    complex <double> e[EigValNumber+EigValNumber];

    int lda=EigValNumber;
    int ldu=lda;
    int ldv=lda;
    int job=01;

    complex <double> s[EigValNumber+EigValNumber];
    complex <double> u_svd[EigValNumber*EigValNumber];
    complex <double> v[EigValNumber*EigValNumber];

    //przepisanie macierzy Boundaries
    for(int i=0; i < EigValNumber; i++) {
        for(int j=0; j < EigValNumber; j++) {
            Boundaries_svd[i+j*lda] = complex <double> ( GSL_REAL(gsl_matrix_complex_get(Boundaries, i, j)),  GSL_IMAG(gsl_matrix_complex_get(Boundaries, i, j)));
        }
    }

    zsvdc( Boundaries_svd, lda, EigValNumber, EigValNumber, s, e, u_svd, ldu, v, ldv, job );

    gsl_complex u, ux, uy, uz;

    if(polarisation==3) {

        for(double x=0; x < x_lenght; x=x+precision) {

            double y=x*ky/kx;

            for(double z=0; z < z_lenght; z=z+precision) {
                ux = gsl_complex_rect(0, 0);
                uy = gsl_complex_rect(0, 0);
                uz = gsl_complex_rect(0, 0);

                //pasek postepu
                int postep=int(100*x/x_lenght);
                Progress->setValue(postep);
                Progress->update();
                QApplication::processEvents();

                for (int Nx=-RecVec, i=0; Nx <= RecVec; Nx++) {
                    for (int Ny=-RecVec; Ny <= RecVec; Ny++, i=i+3) {

                        double Gx=2*M_PI*Nx/a;
                        double Gy=2*M_PI*Ny/a;

                        for (int L=0; L < EigValNumber; L++) {

                            gsl_complex EigenValue = gsl_vector_complex_get(ChosenValues, L);
                            gsl_complex Ax = gsl_matrix_complex_get(ChosenVectors, i, L);
                            gsl_complex Ay = gsl_matrix_complex_get(ChosenVectors, i+1, L);
                            gsl_complex Az = gsl_matrix_complex_get(ChosenVectors, i+2, L);
                            gsl_complex Cl = gsl_complex_rect(real(v[L+(EigValNumber-1)*ldv]), imag(v[L+(EigValNumber-1)*ldv]));

                            gsl_complex expin1 = gsl_complex_rect(0, (kx+Gx)*x+(ky+Gy)*y);
                            gsl_complex exp1 = gsl_complex_exp(expin1);

                            gsl_complex exp2 = gsl_complex_exp(gsl_complex_mul(gsl_complex_rect(0,1), gsl_complex_mul_real(EigenValue, z)));

                            gsl_complex multiply = gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(Cl, Ax), exp2), exp1);
                            ux = gsl_complex_add(ux, multiply);

                            multiply = gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(Cl, Ay), exp2), exp1);
                            uy = gsl_complex_add(uy, multiply);

                            multiply = gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(Cl, Az), exp2), exp1);
                            uz = gsl_complex_add(uz, multiply);
                        }
                    }

                }

                double U = sqrt(gsl_complex_abs(ux)*gsl_complex_abs(ux)+gsl_complex_abs(uy)*gsl_complex_abs(uy)+gsl_complex_abs(uz)*gsl_complex_abs(uz));



                //zapisanie wartości do pliku
                QString nazwa = Osrodek.itsBasis.getFillingSubstance().getName();
                if(nazwa=="Air")
                {
                    double rad=Osrodek.itsBasis.getRadius();
                    if(sqrt(x*x+y*y) > rad && sqrt((a-x)*(a-x)) > rad)
                    {
                        plik << sqrt(x*x+y*y) << "\t" << -z << "\t" << U << "\n";
                    } else {
                        plik << sqrt(x*x+y*y) << "\t" << -z << "\t" << 0 << "\n";
                    }
                } else {
                    plik << sqrt(x*x+y*y) << "\t" << -z << "\t" << U << "\n";
                }

            }

            plik << "\n";
        }

    } else if(polarisation==4) {

        gsl_complex uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz;
        double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;

            for(double x=0; x < x_lenght; x=x+precision) {

                double y=x*ky/kx;

                for(double z=0; z < z_lenght; z=z+precision) {
                    uxx = gsl_complex_rect(0, 0);
                    uxy = gsl_complex_rect(0, 0);
                    uxz = gsl_complex_rect(0, 0);
                    uyx = gsl_complex_rect(0, 0);
                    uyy = gsl_complex_rect(0, 0);
                    uyz = gsl_complex_rect(0, 0);
                    uzx = gsl_complex_rect(0, 0);
                    uzy = gsl_complex_rect(0, 0);
                    uzz = gsl_complex_rect(0, 0);

                    //pasek postepu
                    int postep=int(100*x/x_lenght);
                    Progress->setValue(postep);
                    Progress->update();
                    QApplication::processEvents();

                    double rad=Osrodek.itsBasis.getRadius();
                    double Ela[6][6];
                    if(sqrt(x*x+y*y) > rad && sqrt((a-x)*(a-x)) > rad)
                    {
                        Osrodek.itsBasis.getSubstance().getElasticity(Ela);
                    } else {
                        Osrodek.itsBasis.getFillingSubstance().getElasticity(Ela);
                    }

                    for (int Nx=-RecVec, i=0; Nx <= RecVec; Nx++) {
                        for (int Ny=-RecVec; Ny <= RecVec; Ny++, i=i+3) {

                            double Gx=2*M_PI*Nx/a;
                            double Gy=2*M_PI*Ny/a;

                            for (int L=0; L < EigValNumber; L++) {

                                gsl_complex EigenValue = gsl_vector_complex_get(ChosenValues, L);
                                gsl_complex Ax = gsl_matrix_complex_get(ChosenVectors, i, L);
                                gsl_complex Ay = gsl_matrix_complex_get(ChosenVectors, i+1, L);
                                gsl_complex Az = gsl_matrix_complex_get(ChosenVectors, i+2, L);
                                gsl_complex Cl = gsl_complex_rect(real(v[L+(EigValNumber-1)*ldv]), imag(v[L+(EigValNumber-1)*ldv]));

                                gsl_complex expin1 = gsl_complex_rect(0, (kx+Gx)*x+(ky+Gy)*y);
                                gsl_complex exp1 = gsl_complex_exp(expin1);

                                gsl_complex exp2 = gsl_complex_exp(gsl_complex_mul(gsl_complex_rect(0,1), gsl_complex_mul_real(EigenValue, z)));

                                gsl_complex multiply = gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(Cl, Ax), exp2), exp1);

                                gsl_complex multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), kx+Gx));
                                uxx = gsl_complex_add(uxx, multidiff);

                                multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), ky+Gy));
                                uxy = gsl_complex_add(uxy, multidiff);

                                multidiff = gsl_complex_mul(multiply, gsl_complex_mul(gsl_complex_rect(0,1), EigenValue));
                                uxz = gsl_complex_add(uxz, multidiff);

                                //uy
                                multiply = gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(Cl, Ay), exp2), exp1);

                                multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), kx+Gx));
                                uyx = gsl_complex_add(uyx, multidiff);

                                multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), ky+Gy));
                                uyy = gsl_complex_add(uyy, multidiff);

                                multidiff = gsl_complex_mul(multiply, gsl_complex_mul(gsl_complex_rect(0,1), EigenValue));
                                uyz = gsl_complex_add(uyz, multidiff);

                                //uz
                                multiply = gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(Cl, Az), exp2), exp1);

                                multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), kx+Gx));
                                uzx = gsl_complex_add(uzx, multidiff);

                                multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), ky+Gy));
                                uzy = gsl_complex_add(uzy, multidiff);

                                multidiff = gsl_complex_mul(multiply, gsl_complex_mul(gsl_complex_rect(0,1), EigenValue));
                                uzz = gsl_complex_add(uzz, multidiff);
                            }
                        }

                    }

                    Uxx=gsl_complex_abs(uxx);
                    Uxy=gsl_complex_abs(uxy);
                    Uxz=gsl_complex_abs(uxz);
                    Uyx=gsl_complex_abs(uyx);
                    Uyy=gsl_complex_abs(uyy);
                    Uyz=gsl_complex_abs(uyz);
                    Uzx=gsl_complex_abs(uzx);
                    Uzy=gsl_complex_abs(uzy);
                    Uzz=gsl_complex_abs(uzz);

                    double U = Ela[0][0]*(Uxx*Uxx+Uyy*Uyy+Uzz*Uzz)+2*Ela[0][1]*(Uxx*Uyy+Uxx*Uzz+Uyy*Uzz)+0.25*Ela[3][3]*((Uyz+Uzy)*(Uyz+Uzy)+(Uxz+Uzx)*(Uxz+Uzx)+(Uxy+Uyx)*(Uxy+Uyx));


                    //zapisanie wartości do pliku
                    QString nazwa = Osrodek.itsBasis.getFillingSubstance().getName();
                    if(nazwa=="Air")
                    {
                        double rad=Osrodek.itsBasis.getRadius();
                        if(sqrt(x*x+y*y) > rad && sqrt((a-x)*(a-x)) > rad)
                        {
                            plik << sqrt(x*x+y*y) << "\t" << -z << "\t" << U << "\n";
                        } else {
                            plik << sqrt(x*x+y*y) << "\t" << -z << "\t" << 0 << "\n";
                        }
                    } else {
                        plik << sqrt(x*x+y*y) << "\t" << -z << "\t" << U << "\n";
                    }

                }

                plik << "\n";
            }

    } else {


        for(double x=0; x < x_lenght; x=x+precision) {

            double y=x*ky/kx;

            for(double z=0; z < z_lenght; z=z+precision) {
                u = gsl_complex_rect(0, 0);


                //pasek postepu
                int postep=int(100*x/x_lenght);
                Progress->setValue(postep);
                Progress->update();
                QApplication::processEvents();

                for (int Nx=-RecVec, i=0; Nx <= RecVec; Nx++) {
                    for (int Ny=-RecVec; Ny <= RecVec; Ny++, i=i+3) {

                        double Gx=2*M_PI*Nx/a;
                        double Gy=2*M_PI*Ny/a;

                        for (int L=0; L < EigValNumber; L++) {

                            gsl_complex EigenValue = gsl_vector_complex_get(ChosenValues, L);
                            gsl_complex EigenVector = gsl_matrix_complex_get(ChosenVectors, i+polarisation, L);
                            gsl_complex Cl = gsl_complex_rect(real(v[L+(EigValNumber-1)*ldv]), imag(v[L+(EigValNumber-1)*ldv]));

                            gsl_complex expin1 = gsl_complex_rect(0, (kx+Gx)*x+(ky+Gy)*y);
                            gsl_complex exp1 = gsl_complex_exp(expin1);

                            gsl_complex exp2 = gsl_complex_exp(gsl_complex_mul(gsl_complex_rect(0,1), gsl_complex_mul_real(EigenValue, z)));

                            gsl_complex multiply = gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(Cl, EigenVector), exp2), exp1);
                            u = gsl_complex_add(u, multiply);
                        }
                    }

                }

                double U = gsl_complex_abs(u);

                //zapisanie wartości do pliku
                plik << sqrt(x*x+y*y) << "\t" << -z << "\t" << U << "\n";

            }

            plik << "\n";
        }

    }


    plik.close();
    gsl_matrix_free(gammaC);
    gsl_matrix_free(gammaD);
    gsl_eigen_genv_free(wspce);
    gsl_matrix_complex_free(Boundaries);
    gsl_matrix_complex_free(ChosenVectors);
    gsl_vector_complex_free(ChosenValues);

}
