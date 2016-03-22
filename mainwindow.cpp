#include "mainwindow.h"
#include "cmath"
#include "fstream"
#include "iostream"
#include "newitem.hpp"
#include "agree.hpp"
#include "ui_mainwindow.h"
#include <QtCore/QTextStream>
#include "newitem.cpp"
#include "agree.cpp"
#include "spectrometer.hpp"
#include "konstruktory.cpp"
#include "countDispersion_BS.cpp"
#include "countDispersion_B.cpp"
#include "countDispersion_bulk_B.cpp"
#include "drawStructure.cpp"
#include "countAmplitudes_B.cpp"
#include "countAmplitudes_bulk_B.cpp"
#include "countAmplitudes_BS.cpp"
#include "countBCDet_B.cpp"
#include "countBCDet_BS.cpp"
#include "countBCDet_bulk_B.cpp"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    //współczynnik wypełnienia
    countFillingFraction();

    //wczytanie listy substancji
        ifstream plik("materials.dat",ios::binary);
        if (!plik) {
            Substance Subs;

            SubstanceList.push_back(Subs);

            ui->comboBox_5->addItem(SubstanceList[0].getName());
            ui->comboBox_2->addItem(SubstanceList[0].getName());
            ui->comboBox_3->addItem(SubstanceList[0].getName());
            ui->comboBox_6->addItem(SubstanceList[0].getName());
            ui->comboBox_7->addItem(SubstanceList[0].getName());


        } else {

          ifstream is("materials.dat", ios::binary);
            int m;
          if(is.read((char*) &m, sizeof(int)))
            {
                   SubstanceList.resize(m);

                   for(int i = 0; i < m && SubstanceList[i].Load(is); i++);
            }
      is.close();

    ui->comboBox_5->clear();
    ui->comboBox_2->clear();
    ui->comboBox_3->clear();
    ui->comboBox_6->clear();
    ui->comboBox_7->clear();
    for (int i=0; i<m; i++) {
        ui->comboBox_2->addItem(SubstanceList[i].getName());
         ui->comboBox_5->addItem(SubstanceList[i].getName());
         ui->comboBox_3->addItem(SubstanceList[i].getName());
         ui->comboBox_6->addItem(SubstanceList[i].getName());
         ui->comboBox_7->addItem(SubstanceList[i].getName());
     }
        }
}

MainWindow::~MainWindow()
{
    //zapisanie listy substancji
    ofstream os ("materials.dat", ios::binary);
    int n = SubstanceList.size();
    os.write((char*) &n, sizeof(int));
    for(int i = 0; i < n; i++) SubstanceList[i].Save(os);
    os.close();

    delete ui;
}


void MainWindow::on_actionProduct_Help_activated()
{
   ui->textBrowser->setPlainText("Product Help");
}

void MainWindow::on_actionAbout_P_Surve_activated()
{
   ui->textBrowser->setPlainText("About...");
}

void MainWindow::updateData()
{
//    if (ui->checkBox->checkState()==2) {
//        Osrodek.enableLayer(true);
//        Osrodek.setThickness("Layer", ui->doubleSpinBox->value());
//        Layer Warstwa;
 //       Warstwa.setSubstance(SubstanceList[]);
//    } else {
//        Osrodek.enableLayer(false);
//    }
    if (ui->checkBox_2->isChecked()) {

        Osrodek.enableStructure(true);

        Osrodek.itsStructure.setThickness(ui->doubleSpinBox_2->value());//*1E-6);
        Osrodek.itsStructure.setSubstance(SubstanceList[ui->comboBox_2->currentIndex()]);
        Osrodek.itsStructure.setFillingSubstance(SubstanceList[ui->comboBox_3->currentIndex()]);

        Osrodek.itsBasis.setSubstance(SubstanceList[ui->comboBox_6->currentIndex()]);
        Osrodek.itsBasis.setFillingSubstance(SubstanceList[ui->comboBox_7->currentIndex()]);
        Osrodek.itsBasis.setLatticeConstant(ui->doubleSpinBox_3->value());//*1E-6);
        Osrodek.itsBasis.setRadius(ui->doubleSpinBox_4->value());//*1E-6);

    } else {
        Osrodek.enableStructure(false);

        Osrodek.itsBasis.setSubstance(SubstanceList[ui->comboBox_6->currentIndex()]);
        Osrodek.itsBasis.setFillingSubstance(SubstanceList[ui->comboBox_7->currentIndex()]);
        Osrodek.itsBasis.setLatticeConstant(ui->doubleSpinBox_3->value());//*1E-6);
        Osrodek.itsBasis.setRadius(ui->doubleSpinBox_4->value());//*1E-6);
    }


//kontrola
//double getWel[3];
//WPack.getVelocities(getWel);
//QString qStr = QString::number(getWel[1]);
//ui->textBrowser->setPlainText(qStr);

}

void MainWindow::countFillingFraction()
{
    double a = ui->doubleSpinBox_3->value();
    double r = ui->doubleSpinBox_4->value();
    double ff=M_PI*r*r/(a*a);
    ui->doubleSpinBox_8->setValue(ff);
}

void MainWindow::StructureActive()
{
    if(ui->radioButton->isChecked()) {
        ui->checkBox_2->setChecked(false);
    }
    if (ui->checkBox_2->isChecked() && !ui->radioButton->isChecked()) {
        ui->groupBox_5->setEnabled(true);
    } else {
        ui->checkBox_2->setChecked(false);
        ui->groupBox_5->setEnabled(false);
}

}


//************************Dodawanie substancji****************************
void MainWindow::NewItem()
{
//wywołanie okna
newitem *dial =new newitem;
if (dial->exec()==QDialog::Accepted) {
    dial->setAttribute(Qt::WA_DeleteOnClose);
    ui->textBrowser->setText("O kurwa");

    //pobranie danych z okna
    double Ela[6][6];
    Ela[0][0]=dial->ui->doubleSpinBox->value();
    Ela[0][1]=Ela[1][0]=dial->ui->doubleSpinBox_2->value();
    Ela[0][2]=Ela[2][0]=dial->ui->doubleSpinBox_3->value();
    Ela[0][3]=Ela[3][0]=dial->ui->doubleSpinBox_4->value();
    Ela[0][4]=Ela[4][0]=dial->ui->doubleSpinBox_5->value();
    Ela[0][5]=Ela[5][0]=dial->ui->doubleSpinBox_6->value();
    Ela[1][1]=dial->ui->doubleSpinBox_7->value();
    Ela[2][1]=Ela[1][2]=dial->ui->doubleSpinBox_8->value();
    Ela[1][3]=Ela[3][1]=dial->ui->doubleSpinBox_9->value();
    Ela[1][4]=Ela[4][1]=dial->ui->doubleSpinBox_10->value();
    Ela[1][5]=Ela[5][1]=dial->ui->doubleSpinBox_11->value();
    Ela[2][2]=dial->ui->doubleSpinBox_12->value();
    Ela[2][3]=Ela[3][2]=dial->ui->doubleSpinBox_13->value();
    Ela[2][4]=Ela[4][2]=dial->ui->doubleSpinBox_14->value();
    Ela[2][5]=Ela[5][2]=dial->ui->doubleSpinBox_15->value();
    Ela[3][3]=dial->ui->doubleSpinBox_16->value();
    Ela[3][4]=Ela[4][3]=dial->ui->doubleSpinBox_17->value();
    Ela[3][5]=Ela[5][3]=dial->ui->doubleSpinBox_18->value();
    Ela[4][4]=dial->ui->doubleSpinBox_19->value();
    Ela[4][5]=Ela[5][4]=dial->ui->doubleSpinBox_20->value();
    Ela[5][5]=dial->ui->doubleSpinBox_21->value();

    ChangeToPa(Ela);

    //zapisanie materiału
    Substance Subs;
    Subs.setDensity(dial->ui->doubleSpinBox_37->value());
    Subs.setElasticity(Ela);
    Subs.setName(dial->ui->lineEdit->text());

    SubstanceList.push_back(Subs);

    //odświeżenie listy
    int n = SubstanceList.size();

    ui->comboBox_5->addItem(SubstanceList[n-1].getName());
    ui->comboBox_2->addItem(SubstanceList[n-1].getName());
    ui->comboBox_3->addItem(SubstanceList[n-1].getName());
    ui->comboBox_6->addItem(SubstanceList[n-1].getName());
    ui->comboBox_7->addItem(SubstanceList[n-1].getName());
}
}


//************************Edycja substancji****************************

void MainWindow::EditItem()
{
//inicjalizacja okna
newitem *dial =new newitem;
dial->ui->label_4->setText("Edit Substance");

//pobranie danych do okna
double Ela[6][6];
Substance &CurrentSubs=SubstanceList[ui->comboBox_5->currentIndex()];
CurrentSubs.getElasticity(Ela);
QString name=CurrentSubs.getName();
double density=CurrentSubs.getDensity();

ChangeToGPa(Ela);
ui->textBrowser->setText("Tekx");

dial->ui->lineEdit->setText(name);
dial->ui->doubleSpinBox_37->setValue(density);
dial->ui->doubleSpinBox->setValue(Ela[0][0]);
dial->ui->doubleSpinBox_2->setValue(Ela[0][1]);
dial->ui->doubleSpinBox_3->setValue(Ela[0][2]);
dial->ui->doubleSpinBox_4->setValue(Ela[0][3]);
dial->ui->doubleSpinBox_5->setValue(Ela[0][4]);
dial->ui->doubleSpinBox_6->setValue(Ela[0][5]);
dial->ui->doubleSpinBox_7->setValue(Ela[1][1]);
dial->ui->doubleSpinBox_8->setValue(Ela[2][1]);
dial->ui->doubleSpinBox_9->setValue(Ela[3][1]);
dial->ui->doubleSpinBox_10->setValue(Ela[4][1]);
dial->ui->doubleSpinBox_11->setValue(Ela[5][1]);
dial->ui->doubleSpinBox_12->setValue(Ela[2][2]);
dial->ui->doubleSpinBox_13->setValue(Ela[3][2]);
dial->ui->doubleSpinBox_14->setValue(Ela[4][2]);
dial->ui->doubleSpinBox_15->setValue(Ela[5][2]);
dial->ui->doubleSpinBox_16->setValue(Ela[3][3]);
dial->ui->doubleSpinBox_17->setValue(Ela[3][4]);
dial->ui->doubleSpinBox_18->setValue(Ela[3][5]);
dial->ui->doubleSpinBox_19->setValue(Ela[4][4]);
dial->ui->doubleSpinBox_20->setValue(Ela[4][5]);
dial->ui->doubleSpinBox_21->setValue(Ela[5][5]);

//wywołanie okna
if (dial->exec()==QDialog::Accepted) {

    //pobranie danych z okna
    Ela[0][0]=dial->ui->doubleSpinBox->value();
    Ela[0][1]=Ela[1][0]=dial->ui->doubleSpinBox_2->value();
    Ela[0][2]=Ela[2][0]=dial->ui->doubleSpinBox_3->value();
    Ela[0][3]=Ela[3][0]=dial->ui->doubleSpinBox_4->value();
    Ela[0][4]=Ela[4][0]=dial->ui->doubleSpinBox_5->value();
    Ela[0][5]=Ela[5][0]=dial->ui->doubleSpinBox_6->value();
    Ela[1][1]=dial->ui->doubleSpinBox_7->value();
    Ela[2][1]=Ela[1][2]=dial->ui->doubleSpinBox_8->value();
    Ela[1][3]=Ela[3][1]=dial->ui->doubleSpinBox_9->value();
    Ela[1][4]=Ela[4][1]=dial->ui->doubleSpinBox_10->value();
    Ela[1][5]=Ela[5][1]=dial->ui->doubleSpinBox_11->value();
    Ela[2][2]=dial->ui->doubleSpinBox_12->value();
    Ela[2][3]=Ela[3][2]=dial->ui->doubleSpinBox_13->value();
    Ela[2][4]=Ela[4][2]=dial->ui->doubleSpinBox_14->value();
    Ela[2][5]=Ela[5][2]=dial->ui->doubleSpinBox_15->value();
    Ela[3][3]=dial->ui->doubleSpinBox_16->value();
    Ela[3][4]=Ela[4][3]=dial->ui->doubleSpinBox_17->value();
    Ela[3][5]=Ela[5][3]=dial->ui->doubleSpinBox_18->value();
    Ela[4][4]=dial->ui->doubleSpinBox_19->value();
    Ela[4][5]=Ela[5][4]=dial->ui->doubleSpinBox_20->value();
    Ela[5][5]=dial->ui->doubleSpinBox_21->value();

    ChangeToPa(Ela);

    //zapisanie materiału
    CurrentSubs.setDensity(dial->ui->doubleSpinBox_37->value());
    CurrentSubs.setElasticity(Ela);
    CurrentSubs.setName(dial->ui->lineEdit->text());

    //odświeżenie listy
    ui->comboBox_5->setItemText(ui->comboBox_5->currentIndex(), dial->ui->lineEdit->text());
    ui->comboBox_2->setItemText(ui->comboBox_5->currentIndex(), dial->ui->lineEdit->text());
    ui->comboBox_3->setItemText(ui->comboBox_5->currentIndex(), dial->ui->lineEdit->text());
    ui->comboBox_6->setItemText(ui->comboBox_5->currentIndex(), dial->ui->lineEdit->text());
    ui->comboBox_7->setItemText(ui->comboBox_5->currentIndex(), dial->ui->lineEdit->text());
}
}

//************************Usuwanie substancji****************************
void MainWindow::DeleteItem() {
    agree *dial = new agree;
    int sindex=ui->comboBox_5->currentIndex();
    QString sname=SubstanceList[sindex].getName();
    dial->ui->label->setText("Are you sure to remove substance " + sname + "\nfrom database?");


    if (dial->exec()==QDialog::Accepted) {
        ui->comboBox_5->clear();
        ui->comboBox_2->clear();
        ui->comboBox_3->clear();
        ui->comboBox_6->clear();
        ui->comboBox_7->clear();

        SubstanceList.erase(SubstanceList.begin()+sindex);
        int n = SubstanceList.size();

        for(int i = 0; i < n; i++) {
            ui->comboBox_5->addItem(SubstanceList[i].getName());
            ui->comboBox_2->addItem(SubstanceList[i].getName());
            ui->comboBox_3->addItem(SubstanceList[i].getName());
            ui->comboBox_6->addItem(SubstanceList[i].getName());
            ui->comboBox_7->addItem(SubstanceList[i].getName());
        }
}
}

void MainWindow::on_actionCount_BC_Determinant_triggered()
{
    ui->tabWidget->setCurrentIndex(0);
    updateData();

    ui->progressBar->setEnabled(true);
    ui->tab_3->setEnabled(false);
    QApplication::processEvents();

    double V[3];
    V[0]=ui->doubleSpinBox_5->value();//*1E9;
    V[1]=ui->doubleSpinBox_6->value();//*1E9;
    V[2]=ui->doubleSpinBox_7->value();//*1E9;
    int RecVec=ui->spinBox_5->value();
    int k_prec=ui->spinBox_4->value();
    int factor=ui->spinBox->value();
    bool Surface=ui->radioButton_2->isChecked();

    QString MatrixName=Osrodek.itsStructure.getSubstance().getName();
    QString FillingName=Osrodek.itsStructure.getFillingSubstance().getName();
    QString BasisMatrixName=Osrodek.itsBasis.getSubstance().getName();
    QString BasisFillingName=Osrodek.itsBasis.getFillingSubstance().getName();
    QString LatticeConst=QString::number(Osrodek.itsBasis.getLatticeConstant());
    QString radius=QString::number(Osrodek.itsBasis.getRadius());
    QString thick=QString::number(Osrodek.itsStructure.getThickness());
    QString rec=QString::number(RecVec);

Spectrometer *Detector=new Spectrometer(Osrodek, V, k_prec, RecVec, Surface);

QProgressBar *Progress = ui->progressBar;
QTextBrowser *Browser = ui->textBrowser;
QString DataName;

if (Detector->checkSurface())
{

    if (Osrodek.checkStructure())
    {
        DataName="Su_BS_" + FillingName + MatrixName + "_on_" + BasisMatrixName + BasisFillingName + "_a" + LatticeConst + "_r" + radius + "_t" + thick + "_REC" + rec;
        ui->textBrowser->insertPlainText("Waves: surface; \n Basis: yes; \n Filling:" + BasisFillingName + "; Matrix: " + BasisMatrixName + "\n Lattice Constant: " + LatticeConst + "; Radius of cylinders: " + radius + "\n Structure: yes \n " + "Structure thickness: " + thick + "Filling: " + FillingName + "; Matrix: " + MatrixName);
        ui->textBrowser->insertPlainText("\n Counting Dispersion relation...");

        Detector->countDispersion_BS(Osrodek, DataName, Progress, factor);
    } else {
        DataName="Su_B_" + BasisFillingName + "_in_" + BasisMatrixName + "_a" + LatticeConst + "_r" + radius + "_REC" + rec;
        ui->textBrowser->insertPlainText("Waves: surface; \n Basis: yes \n Filling:" + BasisFillingName + "; Matrix: " + BasisMatrixName + "\n Lattice Constant: " + LatticeConst + "; Radius of cylinders: " + radius + "\n Structure: no \n");
        ui->textBrowser->insertPlainText("Counting Dispersion relation...");
        Detector->countDispersion_B(Osrodek, DataName, Progress, factor, Browser);
    }

} else {

    DataName="bulk_B_" + BasisFillingName + "_in_" + BasisMatrixName + "_a" + LatticeConst + "_r" + radius + "_REC" + rec;
    ui->textBrowser->insertPlainText("Waves: bulk; \n Filling:" + BasisFillingName + "; Matrix: " + BasisMatrixName + "\n Lattice Constant: " + LatticeConst + "; Radius of cylinders: " + radius + "\n");
    ui->textBrowser->insertPlainText("Counting Dispersion relation...");
    Detector->countDispersion_bulk_B(Osrodek, DataName, Progress);
}

QString command= "set terminal png size 730,550 \n set tmargin -200 \n set bmargin 0 \n set title 'Dispersion relation' \n set xlabel 'Wave vector' \n set ylabel 'Frequency' \n set output 'results/" + DataName + ".png' \n set pm3d map \n set palette rgb 33,13,10 \n splot 'results/" + DataName + ".dat'";
const char * Ccommand = command.toAscii().data();
gnuplot(Ccommand);

QGraphicsScene *scene = new QGraphicsScene;
QPixmap lPixmap;
lPixmap.load("results/" + DataName + ".png");
scene->addPixmap(lPixmap);
ui->graphicsView->setScene(scene);

//for (int i=0; i <9; i++) {
//    //for (int j=0; j <81; j++) {
//        qStr = QString::number(Eigenv[i]);
//        ui->textBrowser->insertPlainText(qStr);
//        //qStr = QString::number(Eigenv[i][1]);
//        //ui->textBrowser->insertPlainText(" + i " + qStr);
//        //ui->textBrowser->insertPlainText(",\t");
//   // }
//    ui->textBrowser->insertPlainText("; \n");
//}

delete Detector;

ui->textBrowser->insertPlainText("Done.");
ui->tab_3->setEnabled(true);
ui->progressBar->setEnabled(false);
ui->tabWidget->setCurrentIndex(2);

}

void MainWindow::on_actionFinish_triggered()
{
    QApplication::quit();
}

void MainWindow::on_actionDraw_Structure_triggered()
{
    ui->tabWidget->setCurrentIndex(0);
    updateData();

    double V[3];
    V[0]=ui->doubleSpinBox_5->value();
    V[1]=ui->doubleSpinBox_6->value();
    V[2]=ui->doubleSpinBox_7->value();
    int RecVec=ui->spinBox_5->value();
    int k_prec=ui->spinBox_4->value();
    bool Surface=ui->radioButton_2->isChecked();

Spectrometer *Detector=new Spectrometer(Osrodek, V, k_prec, RecVec, Surface);
QTextBrowser *Browser = ui->textBrowser;
Detector->drawStructure(Osrodek, Browser);

QString command= "set terminal png size 730,550 \n set tmargin -200 \n set bmargin 0 \n set title 'Structure' \n set output 'results/structure.png' \n set pm3d map \n set palette rgb 33,13,10 \n splot 'structure.dat'";
QByteArray   bytes  = command.toAscii();
const char * Ccommand = bytes.data();

gnuplot(Ccommand);

QGraphicsScene *scene = new QGraphicsScene;
QPixmap lPixmap;
lPixmap.load("results/structure.png");
scene->addPixmap(lPixmap);
ui->graphicsView->setScene(scene);

delete Detector;

ui->textBrowser->insertPlainText("Done.");
ui->tab_3->setEnabled(true);
ui->progressBar->setEnabled(false);
ui->tabWidget->setCurrentIndex(2);
}

void MainWindow::on_action2D_triggered()
{

    double zbocze = ui->doubleSpinBox->value();

    QDir directory;
    QString path = QFileDialog::getOpenFileName(this,tr("Open File"),
                                                directory.path(),
                                                tr("Dat files (*.dat)"));

    directory.setPath(path);
    QString DataName = directory.dirName();

    QByteArray   bytes  = path.toAscii();
    const char * CDataName = bytes.data();
    ifstream plik(CDataName);

    std::ofstream wykres;

    DataName.chop(4);

    QString SData = "results/2D_" + DataName + ".dat";
    bytes  = SData.toAscii();
    CDataName = bytes.data();
    wykres.open(CDataName);

    char ch;
    int a = 0;
    QString Qk, Qw, Qdet, Qk_prev, Qw_prev, Qdet_prev, Qdet_prev2 = "";
    double k_prev, w, w_prev, det, det_prev, det_prev2;
    double det_max = 0;
    bool enter = true;

    while(plik.get(ch))
    {
        if(ch=='\t')
        {
            a++;

        } else if(ch=='\n' && enter)
        {

            QString inf = "-inf";
            if(Qdet==inf)   Qdet='0';

            if(!(Qdet_prev.isEmpty()) && !(Qdet_prev2.isEmpty()))
            {

                det = Qdet.toDouble();
                det_prev = Qdet_prev.toDouble();
                det_prev2 = Qdet_prev2.toDouble();

                k_prev = Qk_prev.toDouble();
                w_prev = Qw_prev.toDouble();
                w = Qw.toDouble();

                if(det_prev >= det_prev2 && det_prev >= det && det_prev > det_max)
                {
                    det_max = det_prev;
                }


                if(det_prev < det && det_prev < det_prev2 &&  det_max-det_prev > zbocze && k_prev != 0.707107 && k_prev != 1.70711)
                {
                    wykres << k_prev << "\t" << w_prev << "\n";
                    det_max=-2000;
                }
            }

            a=0;
            Qdet_prev2 = Qdet_prev;
            Qdet_prev = Qdet;
            Qdet.clear();

            Qk_prev = Qk;
            Qk.clear();

            Qw_prev = Qw;
            Qw.clear();
            enter=false;

            if(w < w_prev) det_max = det;

        } else {
            enter=true;
            if(a==0)
            {
                QString Tk(ch);
                Qk += Tk;
            } else if(a==1)
            {
                QString Tw(ch);
                Qw += Tw;
            } else if(a==2)
            {
                QString Td(ch);
                Qdet += Td;
            }
        }

    }
    wykres.close();
    plik.close();


    QString command= "set terminal png size 730,450 \n set title 'Dispersion relation' \n set xlabel 'Wave vector' \n set ylabel 'Frequency' \n set output 'results/2D_" + DataName + ".png' \n plot 'results/2D_" + DataName + ".dat'";
    bytes  = command.toAscii();
    const char * Ccommand = bytes.data();
    gnuplot(Ccommand);

    QGraphicsScene *scene = new QGraphicsScene;
    QPixmap lPixmap;
    lPixmap.load("results/2D_" + DataName + ".png");
    scene->addPixmap(lPixmap);
    ui->graphicsView->setScene(scene);

    ui->tabWidget->setCurrentIndex(2);

}

void MainWindow::on_actionCount_Amplitudes_triggered()
{
    double kx = ui->doubleSpinBox_9->value();
    double ky = ui->doubleSpinBox_11->value();
    double w = ui->doubleSpinBox_10->value();

    updateData();

    ui->progressBar->setEnabled(true);
    ui->tab_3->setEnabled(false);
    QApplication::processEvents();

    double V[3];
    V[0]=ui->doubleSpinBox_5->value();
    V[1]=ui->doubleSpinBox_6->value();
    V[2]=ui->doubleSpinBox_7->value();
    int RecVec=ui->spinBox_5->value();
    int k_prec=ui->spinBox_4->value();
    int polarisation=ui->comboBox->currentIndex();
    double x_lenght = ui->doubleSpinBox_12->value();
    double z_lenght = ui->doubleSpinBox_13->value();
    double precision = ui->doubleSpinBox_14->value();
    bool Surface=ui->radioButton_2->isChecked();
    int factor=ui->spinBox->value();

    QString MatrixName=Osrodek.itsStructure.getSubstance().getName();
    QString FillingName=Osrodek.itsStructure.getFillingSubstance().getName();
    QString BasisMatrixName=Osrodek.itsBasis.getSubstance().getName();
    QString BasisFillingName=Osrodek.itsBasis.getFillingSubstance().getName();
    QString LatticeConst=QString::number(Osrodek.itsBasis.getLatticeConstant());
    QString radius=QString::number(Osrodek.itsBasis.getRadius());
    QString thick=QString::number(Osrodek.itsStructure.getThickness());
    QString rec=QString::number(RecVec);
    QString polar=ui->comboBox->currentText();
    QString DataName;

    QProgressBar *Progress = ui->progressBar;
    QTextBrowser *Browser = ui->textBrowser;

    Spectrometer *Detector=new Spectrometer(Osrodek, V, k_prec, RecVec, Surface);

    if (Detector->checkSurface())
    {

        if (Osrodek.checkStructure())
        {
            DataName="Amp_" + polar +"_kx_y" + QString::number(kx) + "_" + QString::number(ky) + "_w" + QString::number(w) + "_Su_BS_" + FillingName + MatrixName + "_on_" + BasisMatrixName + BasisFillingName + "_a" + LatticeConst + "_r" + radius + "_t" + thick + "_REC" + rec;
            Detector->countAmplitudes_BS(Osrodek, DataName, Progress, factor, kx, ky, w, polarisation, x_lenght, z_lenght, precision);
        } else {
            DataName="Amp_" + polar + "_kx_y" + QString::number(kx) + "_" + QString::number(ky) + "_w" + QString::number(w) + "_Su_B_" + BasisFillingName + "_in_" + BasisMatrixName + "_a" + LatticeConst + "_r" + radius + "_REC" + rec;
            Detector->countAmplitudes_B(Osrodek, DataName, Progress, kx, ky, w, polarisation, x_lenght, z_lenght, precision);
        }

    } else {
        DataName="Amp_" + polar + "_kx_y" + QString::number(kx) + "_" + QString::number(ky) + "_w" + QString::number(w) + "_bulk_B_" + BasisFillingName + "_in_" + BasisMatrixName + "_a" + LatticeConst + "_r" + radius + "_REC" + rec;
        Detector->countAmplitudes_bulk_B(Osrodek, DataName, Progress, Browser, kx, ky, w, polarisation, x_lenght, z_lenght, precision);
    }




    QString command= "set terminal png size 730,550 \n set tmargin -200 \n set bmargin 0 \n set title 'Amplitudes' \n set xlabel 'x' \n set ylabel '-z or y' \n set output 'results/amplitudes/" + DataName + ".png' \n set pm3d map \n set palette rgb 33,13,10 \n splot 'results/amplitudes/" + DataName + ".dat'";
    QByteArray   bytes  = command.toAscii();
    const char * Ccommand = bytes.data();

    gnuplot(Ccommand);

    QGraphicsScene *scene = new QGraphicsScene;
    QPixmap lPixmap;
    lPixmap.load("results/amplitudes/" + DataName + ".png");
    scene->addPixmap(lPixmap);
    ui->graphicsView->setScene(scene);

    delete Detector;

    ui->tab_3->setEnabled(true);
    ui->progressBar->setEnabled(false);
    ui->tabWidget->setCurrentIndex(2);
}

void MainWindow::on_actionCount_BCDet_triggered()
{
    ui->tabWidget->setCurrentIndex(0);
    updateData();

    ui->progressBar->setEnabled(true);
    ui->tab_3->setEnabled(false);
    QApplication::processEvents();

    double V[3];
    V[0]=ui->doubleSpinBox_5->value();//*1E9;
    V[1]=ui->doubleSpinBox_6->value();//*1E9;
    V[2]=ui->doubleSpinBox_7->value();//*1E9;
    int k_prec=ui->spinBox_4->value();
    int RecVec=ui->spinBox_5->value();
    int factor=ui->spinBox->value();
    bool Surface=ui->radioButton_2->isChecked();
    double kx=ui->doubleSpinBox_15->value();
    double ky=ui->doubleSpinBox_16->value();

    QString MatrixName=Osrodek.itsStructure.getSubstance().getName();
    QString FillingName=Osrodek.itsStructure.getFillingSubstance().getName();
    QString BasisMatrixName=Osrodek.itsBasis.getSubstance().getName();
    QString BasisFillingName=Osrodek.itsBasis.getFillingSubstance().getName();
    QString LatticeConst=QString::number(Osrodek.itsBasis.getLatticeConstant());
    QString radius=QString::number(Osrodek.itsBasis.getRadius());
    QString thick=QString::number(Osrodek.itsStructure.getThickness());
    QString rec=QString::number(RecVec);
    QString Kx=QString::number(kx);
    QString Ky=QString::number(ky);

Spectrometer *Detector=new Spectrometer(Osrodek, V, k_prec, RecVec, Surface);

QProgressBar *Progress = ui->progressBar;
QString DataName;

if (Detector->checkSurface())
{

    if (Osrodek.checkStructure())
    {
        DataName="Su_BS_BCDet_kx" + Kx + "_ky" + Ky + "_" + FillingName + MatrixName + "_on_" + BasisMatrixName + BasisFillingName + "_a" + LatticeConst + "_r" + radius + "_t" + thick + "_REC" + rec;
        ui->textBrowser->insertPlainText("Waves: surface; \n Basis: yes; \n Filling:" + BasisFillingName + "; Matrix: " + BasisMatrixName + "\n Lattice Constant: " + LatticeConst + "; Radius of cylinders: " + radius + "\n Structure: yes \n " + "Structure thickness: " + thick + "Filling: " + FillingName + "; Matrix: " + MatrixName);
        ui->textBrowser->insertPlainText("Counting Boundary condition determinant...");
        Detector->countBCDet_BS(Osrodek, DataName, Progress, factor, kx, ky);

    } else {
        DataName="Su_BCDet_kx" + Kx + "_ky" +Ky + "_" + BasisFillingName + "_in_" + BasisMatrixName + "_a" + LatticeConst + "_r" + radius + "_REC" + rec;
        ui->textBrowser->insertPlainText("Waves: surface; \n Basis: yes \n Filling:" + BasisFillingName + "; Matrix: " + BasisMatrixName + "\n Lattice Constant: " + LatticeConst + "; Radius of cylinders: " + radius + "\n Structure: no \n");
        ui->textBrowser->insertPlainText("Counting Boundary condition determinant...");
        Detector->countBCDet_B(Osrodek, DataName, Progress, factor, kx, ky);
    }

} else {

    DataName="bulk_BCDet_kx" + Kx + "_ky" +Ky + "_" + BasisFillingName + "_in_" + BasisMatrixName + "_a" + LatticeConst + "_r" + radius + "_REC" + rec;
    ui->textBrowser->insertPlainText("Waves: bulk; \n Filling:" + BasisFillingName + "; Matrix: " + BasisMatrixName + "\n Lattice Constant: " + LatticeConst + "; Radius of cylinders: " + radius + "\n");
    ui->textBrowser->insertPlainText("Counting Boundary condition determinant...");
    Detector->countBCDet_bulk_B(Osrodek, DataName, Progress, factor, kx, ky);
}

QString command= "set terminal png size 730,450 \n set title 'Dispersion relation' \n set xlabel 'Frequency' \n set ylabel 'BCDet' \n set output 'results/" + DataName + ".png' \n plot 'results/" + DataName + ".dat'";
const char * Ccommand = command.toAscii().data();
gnuplot(Ccommand);

QGraphicsScene *scene = new QGraphicsScene;
QPixmap lPixmap;
lPixmap.load("results/" + DataName + ".png");
scene->addPixmap(lPixmap);
ui->graphicsView->setScene(scene);

delete Detector;

ui->textBrowser->insertPlainText("Done.");
ui->tab_3->setEnabled(true);
ui->progressBar->setEnabled(false);
ui->tabWidget->setCurrentIndex(2);
}
