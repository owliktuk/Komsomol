#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include <QDir>
#include <QtConcurrentRun>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <cstdlib>
#include <complex>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "medium.hpp"
#include  "linpack_z.hpp"
#include  "blas1_z.hpp"


template<typename T>
inline bool isanyinf(T value)
{
return value >= std::numeric_limits<T>::min() && value <=
std::numeric_limits<T>::max();
}

namespace Ui
{
    class MainWindow;
}

class newitem;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void NewItem();
    void EditItem();
    void DeleteItem();
    void countFillingFraction();

    void on_actionProduct_Help_activated();

    void on_actionAbout_P_Surve_activated();

    void updateData();
    void StructureActive();

    void on_actionCount_BC_Determinant_triggered();

    void on_actionFinish_triggered();

    void on_actionDraw_Structure_triggered();

    void on_action2D_triggered();

    void on_actionCount_Amplitudes_triggered();

    void on_actionCount_BCDet_triggered();

private:
    Ui::MainWindow *ui;
    Medium Osrodek;
    vector<Substance> SubstanceList;
};

#endif // MAINWINDOW_H
