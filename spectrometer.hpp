class Spectrometer
{
        public:
                Spectrometer(Medium &Osrodek, double (&Vel)[3], int k_prec, int RecVec, bool surface);
                ~Spectrometer() {}
                int getRecVectors() { return itsRecVectors; }
                void setRecVectors(int R) { itsRecVectors=R; }
                void getFrequencies(double (&Velocity)[3]) const;
                void setFrequencies(double (&Velocity)[3]);
                bool checkSurface() { return isSurface; }
                void enableSurface(bool sur) { isSurface=sur; }
                void setRecBasisSubstance(Substance S) { itsRecBasisSubstance.push_back(S); }
                Substance getRecBasisSubstance(int i) { return itsRecBasisSubstance[i]; }
                void setRecLayerSubstance(Substance S) { itsRecLayerSubstance.push_back(S); }
                Substance getRecLayerSubstance(int i) { return itsRecLayerSubstance[i]; }
                void setRecStructureSubstance(Substance S) { itsRecStructureSubstance.push_back(S); }
                Substance getRecStructureSubstance(int i) { return itsRecStructureSubstance[i]; }
                void countBCDet();
                void countDispersion_BS(Medium &Osrodek, QString DataName, QProgressBar *Progress, int factor);
                void countDispersion_B(Medium &Osrodek, QString DataName, QProgressBar *Progress, int factor, QTextBrowser *Browser);
                void countDispersion_bulk_B(Medium &Osrodek, QString DataName, QProgressBar *Progress);
                void drawStructure(Medium &Osrodek, QTextBrowser *Browser);
                void countAmplitudes_B(Medium &Osrodek, QString DataName, QProgressBar *Progress, double kx, double ky, double w, int polarisation, double x_lenght, double z_lenght, double precision);
                void countAmplitudes_bulk_B(Medium &Osrodek, QString DataName, QProgressBar *Progress, QTextBrowser *Browser, double kx, double ky, double w, int polarisation, double x_lenght, double y_lenght, double precision);
                void countAmplitudes_BS(Medium &Osrodek, QString DataName, QProgressBar *Progress, int factor, double kx, double ky, double w, int polarisation, double x_lenght, double z_lenght, double precision);
                void countBCDet_BS(Medium &Osrodek, QString DataName, QProgressBar *Progress, int factor, double kx, double ky);
                void countBCDet_B(Medium &Osrodek, QString DataName, QProgressBar *Progress, int factor, double kx, double ky);
                void countBCDet_bulk_B(Medium &Osrodek, QString DataName, QProgressBar *Progress, int factor, double kx, double ky);

        private:
                int itsRecVectors;
                double itsFrequencies[3];
                int itsK_Precision;
                bool isSurface;
                vector<Substance> itsRecBasisSubstance;
                vector<Substance> itsRecLayerSubstance;
                vector<Substance> itsRecStructureSubstance;

};
