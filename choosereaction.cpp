#include "choosereaction.h"
#include "ui_choosereaction.h"

ChooseReaction::ChooseReaction(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ChooseReaction) {

    ui->setupUi(this);

    ui->label->setStyleSheet("QLabel{color : red;}");

    QObject::connect(ui->comboBox_type, SIGNAL(currentIndexChanged(int)), this, SLOT(chooseType(int)));

    model = new QStringListModel();

    chooseType(ui->comboBox_type->currentIndex());
}

ChooseReaction::~ChooseReaction()
{
    delete ui;
}

int ChooseReaction::getType(void) {

    switch(ui->comboBox_type->currentIndex()) {
    case 0:
        return RCTFLAG_E;
    case 1:
        return RCTFLAG_C;
    case 2:
        return RCTFLAG_C2;
    case 3:
        return RCTFLAG_CD;
    case 4:
        return RCTFLAG_CC;
    default:
        return -1;
    }
}

QString ChooseReaction::getStrReaction(void) {

    QString strReaction;

    switch(ui->comboBox_type->currentIndex()) {
    case 0:
        strReaction = ui->comboBox_A->currentText() + " + e ⇄ " + ui->comboBox_C->currentText();
        break;
    case 1:
        strReaction = ui->comboBox_B->currentText() + " ⇄ " + ui->comboBox_C->currentText();
        break;
    case 2:
        strReaction = ui->comboBox_A->currentText() + " + " + ui->comboBox_B->currentText() + " ⇄ " + ui->comboBox_C->currentText() + " + " + ui->comboBox_D->currentText();
        break;
    case 3:
        strReaction = ui->comboBox_B->currentText() + " ⇄ " + ui->comboBox_C->currentText() + " + " + ui->comboBox_D->currentText();
        break;
    case 4:
        strReaction = ui->comboBox_A->currentText() + " + " + ui->comboBox_B->currentText() + " ⇄ " + ui->comboBox_C->currentText();
        break;
    default:
        strReaction = "Error!";
    }

    return strReaction;
}

/*
 * Shows and hides the choices for different reaction types.
 */
void ChooseReaction::chooseType(int idx) {

    ui->comboBox_A->setModel(model);
    ui->comboBox_B->setModel(model);
    ui->comboBox_C->setModel(model);
    ui->comboBox_D->setModel(model);

    switch(idx) {
    case 0:
        ui->comboBox_A->show(); ui->label_p1->show(); ui->label_electron->show(); ui->comboBox_C->show();
        ui->comboBox_B->hide(); ui->label_p2->hide(); ui->comboBox_D->hide();
        break;
    case 1:
        ui->comboBox_B->show(); ui->comboBox_C->show();
        ui->comboBox_A->hide(); ui->label_p1->hide(); ui->label_electron->hide(); ui->label_p2->hide(); ui->comboBox_D->hide();
        break;
    case 2:
        ui->comboBox_A->show(); ui->label_p1->show(); ui->comboBox_B->show(); ui->comboBox_C->show(); ui->label_p2->show(); ui->comboBox_D->show();
        ui->label_electron->hide();
        break;
    case 3:
        ui->comboBox_B->show(); ui->comboBox_C->show(); ui->label_p2->show(); ui->comboBox_D->show();
        ui->comboBox_A->hide(); ui->label_p1->hide(); ui->label_electron->hide();
        break;
    case 4:
        ui->comboBox_A->show(); ui->label_p1->show(); ui->comboBox_B->show(); ui->comboBox_C->show();
        ui->label_electron->hide(); ui->label_p2->hide(); ui->comboBox_D->hide();
        break;
    default:
        std::cout << "This should never happen!\n";
    }
    // qDebug() << "index = " << idx;
}
