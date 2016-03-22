#ifndef AGREE_HPP
#define AGREE_HPP

#include <QDialog>

namespace Ui {
    class agree;
}

class agree : public QDialog
{
    Q_OBJECT
    friend class MainWindow;

public:
    explicit agree(QWidget *parent = 0);
    ~agree();

private:
    Ui::agree *ui;
};

#endif // AGREE_HPP
