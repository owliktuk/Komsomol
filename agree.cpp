#include "agree.hpp"
#include "ui_agree.h"

agree::agree(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::agree)
{
    ui->setupUi(this);
}

agree::~agree()
{
    delete ui;
}
