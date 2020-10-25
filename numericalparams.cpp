#include "numericalparams.h"
#include "ui_numericalparams.h"

NumericalParams::NumericalParams(QWidget *parent) : QDialog(parent), ui(new Ui::NumericalParams) {
    ui->setupUi(this);

    ui->lineEdit_beta->setValidator(new QDoubleValidator(0.0, 1.0, 2));
    ui->lineEdit_Dm->setValidator(new QDoubleValidator(0.0, 100.0, 2));
    ui->lineEdit_nlerror->setValidator(new QDoubleValidator(1.0e-10, 1.0, 2));
    ui->lineEdit_maxiter->setValidator(new QDoubleValidator(1, 1.0e8, 2));
}

NumericalParams::~NumericalParams() {
    delete ui;
}

void NumericalParams::setBeta(double beta) {ui->lineEdit_beta->setText(QString::number(beta));}
double NumericalParams::getBeta(void) {return ui->lineEdit_beta->text().toDouble();}
void NumericalParams::setDm(double Dm) {ui->lineEdit_Dm->setText(QString::number(Dm));}
double NumericalParams::getDm(void) {return ui->lineEdit_Dm->text().toDouble();}
void NumericalParams::setNLerror(double nlerror) {ui->lineEdit_nlerror->setText(QString::number(nlerror));}
double NumericalParams::getNLerror(void) {return ui->lineEdit_nlerror->text().toDouble();}
void NumericalParams::setMaxiter(int maxiter) {ui->lineEdit_maxiter->setText(QString::number(maxiter));}
int NumericalParams::getMaxiter(void) {return ui->lineEdit_maxiter->text().toInt();}
void NumericalParams::setNLbox(bool check) {ui->checkBox_NL->setChecked(check);}
bool NumericalParams::getNLbox(void) {return ui->checkBox_NL->isChecked();}
