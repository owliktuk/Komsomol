#ifndef NEWITEM_HPP
#define NEWITEM_HPP

#include <QDialog>

namespace Ui {
    class newitem;
}

class newitem : public QDialog
{
    Q_OBJECT

public:
    explicit newitem(QWidget *parent = 0);
    ~newitem();
     Ui::newitem *ui;
};

#endif // NEWITEM_HPP
