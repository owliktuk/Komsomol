//deklaracje klas
using namespace std;

class Substance
{
      public:
            Substance();
            Substance (const Substance &);
            ~Substance() {}
             QString getName() const {return itsName; }
             void getElasticity(double (&Elasticity)[6][6]) const;
             double getDensity() const { return itsDensity; }
             void setName(QString N) { itsName = N; }
             void setElasticity(double (&Elasticity)[6][6]);
             void setDensity(double D) { itsDensity = D; }
             ostream& Save(ostream &os);
             istream& Load(istream &is);
      public:
              QString itsName;
              double itsElasticity[6][6];
              double itsDensity;
};

class Basis
{
      public:
             Basis();
             ~Basis() {}
             void setSubstance(Substance S) { itsSubstance=S; }
             Substance getSubstance() { return itsSubstance; }
             void setFillingSubstance(Substance S) { itsFillingSubstance=S; }
             void setLatticeConstant(double L) { itsLatticeConstant=L; }
             void setRadius(double R) { itsRadius=R; }
             Substance getFillingSubstance() { return itsFillingSubstance; }
             double getLatticeConstant() { return itsLatticeConstant; }
             double getRadius() { return itsRadius; }
      protected:
              Substance itsSubstance;
              Substance itsFillingSubstance;
              double itsLatticeConstant;
              double itsRadius;
};

class Layer : public Basis
{
      public:
             Layer();
             ~Layer() {}
             double getThickness() const { return itsThickness; }
             void setThickness(double Thickness) { itsThickness=Thickness; }
      protected:
              double itsThickness;
};

class Medium
{
    friend class Spectrometer;
      public:
             Medium();
             ~Medium() {}
             Basis getBasis() const { return itsBasis; }
             Layer getLayer() const { return itsLayer; }
             Layer getStructure() const { return itsStructure; }
             void setBasis(Basis B) { itsBasis=B; }
             void setLayer(Layer L) { itsLayer=L; }
             void setStructure(Layer S) { itsStructure=S; }
             double getThickness(QString element) const;
             void setThickness(QString element, int value);
             void getBCDet() const;
             void getDispersion() const;
             void propagate();
             bool checkLayer() const { return isLayer; }
             bool checkStructure() const { return isStructure; }
             void enableLayer(bool value) { isLayer=value; }
             void enableStructure(bool value) { isStructure=value; }
    public:
              Basis itsBasis;
              Layer itsLayer;
              Layer itsStructure;
              bool isLayer;
              bool isStructure;
              double BCDet[1];
              double Dispersion[1];
              
};

//class Spectrometer
//{
//	public:
//                Spectrometer(Medium &Osrodek, double (&Vel)[3], int RecVec);
//		~Spectrometer() {}
//                int getRecVectors() { return itsRecVectors; }
//                void setRecVectors(int R) { itsRecVectors=R; }
//                void getVelocities(double (&Velocity)[3]) const;
//                void setVelocities(double (&Velocity)[3]);
//                void setRecBasisSubstance(Substance S) { itsRecBasisSubstance.push_back(S); }
//                Substance getRecBasisSubstance(int i) { return itsRecBasisSubstance[i]; }
//                void setRecLayerSubstance(Substance S) { itsRecLayerSubstance.push_back(S); }
//                Substance getRecLayerSubstance(int i) { return itsRecLayerSubstance[i]; }
//                void setRecStructureSubstance(Substance S) { itsRecStructureSubstance.push_back(S); }
//                Substance getRecStructureSubstance(int i) { return itsRecStructureSubstance[i]; }
//                void countBCDet();
//                void countDispersion(Medium &Osrodek, double (&Eigenv)[54][2]);
//        private:
//                int itsRecVectors;
//                double itsVelocities[3];
//                vector<Substance> itsRecBasisSubstance;
//                vector<Substance> itsRecLayerSubstance;
//                vector<Substance> itsRecStructureSubstance;
//};
