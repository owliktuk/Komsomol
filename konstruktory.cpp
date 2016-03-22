//***********************************************************************
//konstruktory klas

Medium::Medium() {
}

Basis::Basis() {
}

Layer::Layer() {
}

Substance::Substance() {
    itsName="Nickel";
    int i, j;
    for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
            itsElasticity[i][j]=0;
        }
    }
    itsElasticity[0][0]=itsElasticity[1][1]=itsElasticity[2][2]=20000000000;
    itsElasticity[0][1]=itsElasticity[1][0]=itsElasticity[2][1]=itsElasticity[1][2]=itsElasticity[0][2]=itsElasticity[2][0]=5300000000;
    itsElasticity[3][3]=itsElasticity[4][4]=itsElasticity[5][5]=50000000;
    itsDensity=1000;
}

Substance::Substance(const Substance & rhs) {
    itsName = rhs.getName();
    itsDensity = rhs.getDensity();
    double Ela[6][6];
    rhs.getElasticity(Ela);
    for (int i=0; i<6; i++) {
        for (int j=0; j<6; j++) {
            itsElasticity[i][j]=Ela[i][j];
        }
    }
}

//konstruktor klasy Spectrometer
Spectrometer::Spectrometer(Medium &Osrodek, double (&rVel)[3], int k_prec, int RecVec, bool surface) {

    itsFrequencies[0]=rVel[0];
    itsFrequencies[1]=rVel[1];
    itsFrequencies[2]=rVel[2];
    itsK_Precision=k_prec;
    itsRecVectors=RecVec;
    isSurface=surface;

    //stałe struktury
    double R=Osrodek.itsBasis.getRadius();
    double a=Osrodek.itsBasis.getLatticeConstant();
    double ff=M_PI*R*R/(a*a);

    //pobranie substancji
    Substance BasisMatrixSubstance=Osrodek.itsBasis.getSubstance();
    Substance BasisFillingSubstance=Osrodek.itsBasis.getFillingSubstance();
    Substance StructureMatrixSubstance=Osrodek.itsStructure.getSubstance();
    Substance StructureFillingSubstance=Osrodek.itsStructure.getFillingSubstance();

    //parametry substancji w strukturze
    double StructureMatrixElasticity[6][6];
    StructureMatrixSubstance.getElasticity(StructureMatrixElasticity);
    double StructureMatrixDensity=StructureMatrixSubstance.getDensity();

    double StructureFillingElasticity[6][6];
    StructureFillingSubstance.getElasticity(StructureFillingElasticity);
    double StructureFillingDensity=StructureFillingSubstance.getDensity();


    //parametry substancji w podłożu
    double BasisMatrixElasticity[6][6];
    BasisMatrixSubstance.getElasticity(BasisMatrixElasticity);
    double BasisMatrixDensity=BasisMatrixSubstance.getDensity();

    double BasisFillingElasticity[6][6];
    BasisFillingSubstance.getElasticity(BasisFillingElasticity);
    double BasisFillingDensity=BasisFillingSubstance.getDensity();

    //sprężystości i gęstości w przestrzeni odwrotnej
    double BasisRecElasticity[6][6];
    double BasisRecDensity;
    Substance RecBasisSubstance;

    double StructureRecElasticity[6][6];
    double StructureRecDensity;
    Substance RecStructureSubstance;

    double G_lenght, g1x, g1y, g2x, g2y;

    itsRecBasisSubstance.resize(pow(2*RecVec+1, 4));
    itsRecStructureSubstance.resize(pow(2*RecVec+1, 4));

    //obliczenie transformaty Fouriera Cij oraz gęstości
    for(int N1x=-RecVec, c=0; N1x<=RecVec; N1x++) {
        for(int N1y=-RecVec; N1y<=RecVec; N1y++) {
            for(int N2x=-RecVec; N2x<=RecVec; N2x++) {
                for(int N2y=-RecVec; N2y<=RecVec; N2y++, c++) {
                    g1x=2*M_PI*N1x/a;
                    g1y=2*M_PI*N1y/a;
                    g2x=2*M_PI*N2x/a;
                    g2y=2*M_PI*N2y/a;
                    G_lenght=sqrt(pow(g1x-g2x,2)+pow(g1y-g2y,2));

                    for(int i=0; i<6; i++) {
                        for(int j=0; j<6; j++) {
                            if (N1x-N2x==0 && N1y-N2y==0) {
                                if (isSurface) StructureRecElasticity[i][j]=ff*StructureFillingElasticity[i][j]+(1-ff)*StructureMatrixElasticity[i][j];
                                BasisRecElasticity[i][j]=ff*BasisFillingElasticity[i][j]+(1-ff)*BasisMatrixElasticity[i][j];
                            } else {
                                if (isSurface) StructureRecElasticity[i][j]=2*ff*gsl_sf_bessel_J1(G_lenght*R)*(StructureFillingElasticity[i][j]-StructureMatrixElasticity[i][j])/(G_lenght*R);
                                BasisRecElasticity[i][j]=2*ff*gsl_sf_bessel_J1(G_lenght*R)*(BasisFillingElasticity[i][j]-BasisMatrixElasticity[i][j])/(G_lenght*R);
                            }
                        }
                    }
                    if (N1x-N2x==0 && N1y-N2y==0) {
                        if (isSurface) StructureRecDensity=ff*StructureFillingDensity+(1-ff)*StructureMatrixDensity;
                        BasisRecDensity=ff*BasisFillingDensity+(1-ff)*BasisMatrixDensity;
                    } else {
                        if (isSurface) StructureRecDensity=2*ff*gsl_sf_bessel_J1(G_lenght*R)*(StructureFillingDensity-StructureMatrixDensity)/(G_lenght*R);
                        BasisRecDensity=2*ff*gsl_sf_bessel_J1(G_lenght*R)*(BasisFillingDensity-BasisMatrixDensity)/(G_lenght*R);
                    }
                    if (isSurface) RecStructureSubstance.setElasticity(StructureRecElasticity);
                    if (isSurface) RecStructureSubstance.setDensity(StructureRecDensity);
                    if (isSurface) itsRecStructureSubstance[c]=RecStructureSubstance;

                    RecBasisSubstance.setElasticity(BasisRecElasticity);
                    RecBasisSubstance.setDensity(BasisRecDensity);
                    itsRecBasisSubstance[c]=RecBasisSubstance;
                }
            }
        }
    }


};


//***********************************************************************
//akcesory

// *****************************Medium*******************************

double Medium::getThickness(QString element) const {
    if (element=="Layer") return itsLayer.getThickness(); 
    if (element=="Structure") return itsStructure.getThickness();
    }

void Medium::setThickness(QString element, int value) {
    if (element=="Layer") itsLayer.setThickness(value);
    if (element=="Structure") itsStructure.setThickness(value);
    }
    
// *****************************Spectrometer*******************************

void Spectrometer::getFrequencies(double (&rVelocity)[3]) const {
     rVelocity[0]=itsFrequencies[0];
     rVelocity[1]=itsFrequencies[1];
     rVelocity[2]=itsFrequencies[2];
     }

void Spectrometer::setFrequencies(double (&rVelocity)[3]) {
     itsFrequencies[0]=rVelocity[0];
     itsFrequencies[1]=rVelocity[1];
     itsFrequencies[2]=rVelocity[2];
     }
     

// *****************************Substance*******************************

void Substance::getElasticity(double (&Elasticity)[6][6]) const {
    int i, j;
    for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
            Elasticity[i][j]=itsElasticity[i][j];
        }
    }
}

void Substance::setElasticity(double (&Elasticity)[6][6]) {
    int i, j;
    for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
            itsElasticity[i][j]=Elasticity[i][j];
        }
    }
}

ostream& Substance::Save(ostream &os)
        {
                 char Cname[20];
                 strcpy(Cname, itsName.toAscii().constData());
                 os.write((char*) &Cname, sizeof(char[20]));
                 os.write((char*) &itsDensity, sizeof(itsDensity));
                 int i, j;
                                  for (i=0; i<6; i++) {
                                      for (j=0; j<6; j++) {
                                 os.write((char*) &itsElasticity[i][j], sizeof(double));
                                 }
                                 }

                return os;
        }

istream& Substance::Load(istream &is)
        {
                char Cname[20];
                is.read((char*) &Cname, sizeof(char[20]));
                QString Snazwa(Cname);
                itsName=Snazwa;
                is.read((char*) &itsDensity, sizeof(double));
                int i, j;
                                 for (i=0; i<6; i++) {
                                     for (j=0; j<6; j++) {
                                is.read((char*) &itsElasticity[i][j], sizeof(double));
                                }
                                }

                return is;
        }


//***********************************************************************
//inne funkcje

//zamiany jednostek
void ChangeToGPa(double (&Ela)[6][6]) {
    for (int i=0; i<6; i++) {
        for (int j=0; j<6; j++) {
   Ela[i][j]=Ela[i][j]*0.000000001;
   }
   }
}

void ChangeToPa(double (&Ela)[6][6]) {
    for (int i=0; i<6; i++) {
        for (int j=0; j<6; j++) {
   Ela[i][j]=Ela[i][j]*1000000000;
   }
   }
}

//wywołanie gnójplota
void gnuplot(const char *gnucommand)
{
  char syscommand[1024];
  sprintf(syscommand, "echo \"%s\" | gnuplot -persist", gnucommand);
  system(syscommand);
}
