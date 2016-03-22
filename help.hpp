#ifndef HELP_HPP
#define HELP_HPP

#include <QDialog>

namespace Ui {
    class Help;
}

class Help : public QDialog
{
    Q_OBJECT

public:
    explicit Help(QWidget *parent = 0);
    ~Help();

private:
    Ui::Help *ui;
};

#endif // HELP_HPP
