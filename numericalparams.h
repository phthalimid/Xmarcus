#ifndef NUMERICALPARAMS_H
#define NUMERICALPARAMS_H

#include <QDialog>
#include <QDoubleValidator>

#include <iostream>

namespace Ui {
class NumericalParams;
}

class NumericalParams : public QDialog
{
    Q_OBJECT

public:
    explicit NumericalParams(QWidget *parent = nullptr);
    ~NumericalParams();

    void setBeta(double);
    double getBeta(void);
    void setDm(double); double getDm();
    void setNLerror(double);
    double getNLerror(void);
    void setMaxiter(int);
    int getMaxiter(void);
    void setNLbox(bool);
    bool getNLbox(void);

public slots:

private:
    Ui::NumericalParams *ui;
};

#endif // NUMERICALPARAMS_H
