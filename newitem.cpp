#include "newitem.hpp"
#include "mainwindow.h"
#include "ui_newitem.h"
using namespace std;

newitem::newitem(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::newitem)
{
    ui->setupUi(this);
}

newitem::~newitem()
{
    delete ui;
}

