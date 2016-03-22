void Spectrometer::countAmplitudes_bulk_B(Medium &Osrodek, QString DataName, QProgressBar *Progress, QTextBrowser *Browser, double kx, double ky, double w, int polarisation, double x_lenght, double y_lenght, double precision) {

    double a=Osrodek.itsBasis.getLatticeConstant();

    int RecVec=itsRecVectors;

    int dimension=3*(2*RecVec+1)*(2*RecVec+1);

    std::ofstream plik;
    DataName="results/amplitudes/"+ DataName + ".dat";
    QByteArray   bytes  = DataName.toAscii();
    const char * CDataName = bytes.data();
    plik.open(CDataName);

    //inicjalizacje wektorów i macierzy
    gsl_matrix *gamma=gsl_matrix_calloc(dimension, dimension);
    gsl_matrix *V=gsl_matrix_calloc(dimension, dimension);
    gsl_vector *S=gsl_vector_calloc(dimension);
    gsl_vector *work=gsl_vector_calloc(dimension);


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

                }
            }
        }
    }

    //rozwiązanie układu równań
    gsl_linalg_SV_decomp(gamma, V, S, work);

    gsl_complex u;
    double Gx, Gy;

    if (polarisation==3) {

        gsl_complex ux, uy, uz;

        for(double x=0; x < x_lenght; x=x+precision) {
            for(double y=0; y < y_lenght; y=y+precision) {
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

                        Gx=2*M_PI*Nx/a;
                        Gy=2*M_PI*Ny/a;

                        double Ax = gsl_matrix_get(V, i, dimension-1);
                        double Ay = gsl_matrix_get(V, i+1, dimension-1);
                        double Az = gsl_matrix_get(V, i+2, dimension-1);

                        gsl_complex expin1 = gsl_complex_rect(0, (kx+Gx)*x+(ky+Gy)*y);
                        gsl_complex exp1 = gsl_complex_exp(expin1);

                        gsl_complex multiply = gsl_complex_mul_real(exp1, Ax);
                        ux = gsl_complex_add(ux, multiply);

                        expin1 = gsl_complex_rect(0, (kx+Gx)*x+(ky+Gy)*y);
                        exp1 = gsl_complex_exp(expin1);

                        multiply = gsl_complex_mul_real(exp1, Ay);
                        uy = gsl_complex_add(uy, multiply);

                        expin1 = gsl_complex_rect(0, (kx+Gx)*x+(ky+Gy)*y);
                        exp1 = gsl_complex_exp(expin1);

                        multiply = gsl_complex_mul_real(exp1, Az);
                        uz = gsl_complex_add(uz, multiply);
                    }
                }

                double U = sqrt(gsl_complex_abs(ux)*gsl_complex_abs(ux)+gsl_complex_abs(uy)*gsl_complex_abs(uy)+gsl_complex_abs(uz)*gsl_complex_abs(uz));

                //zapisanie wartości do pliku
                QString nazwa = Osrodek.itsBasis.getFillingSubstance().getName();
                if(nazwa=="Air")
                {
                    double rad=Osrodek.itsBasis.getRadius();
                    if(sqrt(x*x+y*y) > rad && sqrt((a-x)*(a-x)+y*y) > rad && sqrt((a-x)*(a-x)+(a-y)*(a-y)) > rad && sqrt((a-y)*(a-y)+x*x) > rad)
                    {
                        plik << x << "\t" << y << "\t" << U << "\n";
                    }
                } else {
                    plik << x << "\t" << y << "\t" << U << "\n";
                }

            }



            plik << "\n";
        }

    } else if (polarisation==4) {

            gsl_complex uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz;
            double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;

            for(double x=0; x < x_lenght; x=x+precision) {
                for(double y=0; y < y_lenght; y=y+precision) {

                    uxx = gsl_complex_rect(0, 0);
                    uxy = gsl_complex_rect(0, 0);
                    uxz = gsl_complex_rect(0, 0);
                    uyx = gsl_complex_rect(0, 0);
                    uyy = gsl_complex_rect(0, 0);
                    uyz = gsl_complex_rect(0, 0);
                    uzx = gsl_complex_rect(0, 0);
                    uzy = gsl_complex_rect(0, 0);
                    uzz = gsl_complex_rect(0, 0);

                    double rad=Osrodek.itsBasis.getRadius();
                    double Ela[6][6];
                    if(sqrt(x*x+y*y) > rad && sqrt((a-x)*(a-x)+y*y) > rad && sqrt((a-x)*(a-x)+(a-y)*(a-y)) > rad && sqrt((a-y)*(a-y)+x*x) > rad)
                    {
                        Osrodek.itsBasis.getSubstance().getElasticity(Ela);
                    } else {
                        Osrodek.itsBasis.getFillingSubstance().getElasticity(Ela);
                    }

                    //pasek postepu
                    int postep=int(100*x/x_lenght);
                    Progress->setValue(postep);
                    Progress->update();
                    QApplication::processEvents();

                    for (int Nx=-RecVec, i=0; Nx <= RecVec; Nx++) {
                        for (int Ny=-RecVec; Ny <= RecVec; Ny++, i=i+3) {

                            Gx=2*M_PI*Nx/a;
                            Gy=2*M_PI*Ny/a;

                            double Ax = gsl_matrix_get(V, i, dimension-1);
                            double Ay = gsl_matrix_get(V, i+1, dimension-1);
                            double Az = gsl_matrix_get(V, i+2, dimension-1);

                            gsl_complex expin1 = gsl_complex_rect(0, (kx+Gx)*x+(ky+Gy)*y);
                            gsl_complex exp1 = gsl_complex_exp(expin1);

                            gsl_complex multiply = gsl_complex_mul_real(exp1, Ax);

                            //uxx
                            gsl_complex multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), kx+Gx));
                            uxx = gsl_complex_add(uxx, multidiff);

                            //uxy
                            multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), ky+Gy));
                            uxy = gsl_complex_add(uxy, multidiff);


                            expin1 = gsl_complex_rect(0, (kx+Gx)*x+(ky+Gy)*y);
                            exp1 = gsl_complex_exp(expin1);

                            //uy
                            multiply = gsl_complex_mul_real(exp1, Ay);

                            //uyx
                            multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), kx+Gx));
                            uyx = gsl_complex_add(uyx, multidiff);

                            //uyy
                            multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), ky+Gy));
                            uyy = gsl_complex_add(uyy, multidiff);

                            expin1 = gsl_complex_rect(0, (kx+Gx)*x+(ky+Gy)*y);
                            exp1 = gsl_complex_exp(expin1);

                            multiply = gsl_complex_mul_real(exp1, Az);

                            //uzx
                            multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), kx+Gx));
                            uzx = gsl_complex_add(uzx, multidiff);

                            //uzy
                            multidiff = gsl_complex_mul(multiply, gsl_complex_mul_real(gsl_complex_rect(0,1), ky+Gy));
                            uzy = gsl_complex_add(uzy, multidiff);
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
                        if(sqrt(x*x+y*y) > rad && sqrt((a-x)*(a-x)+y*y) > rad && sqrt((a-x)*(a-x)+(a-y)*(a-y)) > rad && sqrt((a-y)*(a-y)+x*x) > rad)
                        {
                            plik << x << "\t" << y << "\t" << U << "\n";
                        }
                    } else {
                        plik << x << "\t" << y << "\t" << U << "\n";
                    }

                }



                plik << "\n";
            }


    } else {

        for(double x=0; x < x_lenght; x=x+precision) {
            for(double y=0; y < y_lenght; y=y+precision) {
                u = gsl_complex_rect(0, 0);

                //pasek postepu
                int postep=int(100*x/x_lenght);
                Progress->setValue(postep);
                Progress->update();
                QApplication::processEvents();

                for (int Nx=-RecVec, i=0; Nx <= RecVec; Nx++) {
                    for (int Ny=-RecVec; Ny <= RecVec; Ny++, i=i+3) {

                        Gx=2*M_PI*Nx/a;
                        Gy=2*M_PI*Ny/a;

                        double A = gsl_matrix_get(V, i+polarisation, dimension-1);

                        gsl_complex expin1 = gsl_complex_rect(0, (kx+Gx)*x+(ky+Gy)*y);
                        gsl_complex exp1 = gsl_complex_exp(expin1);

                        gsl_complex multiply = gsl_complex_mul_real(exp1, A);
                        u = gsl_complex_add(u, multiply);
                    }
                }

                double U = gsl_complex_abs(u);

                //zapisanie wartości do pliku
                plik << x << "\t" << y << "\t" << U << "\n";

            }



            plik << "\n";
        }

    }


    plik.close();
    gsl_matrix_free(gamma);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
}



