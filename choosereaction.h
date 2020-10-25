#ifndef CHOOSEREACTION_H
#define CHOOSEREACTION_H

#include <QDialog>
#include <QStringListModel>
#include <iostream>

#include <QDebug>

#include "reaction.h"

namespace Ui {
class ChooseReaction;
}

class ChooseReaction : public QDialog
{
    Q_OBJECT

public:
    explicit ChooseReaction(QWidget *parent = nullptr);
    ~ChooseReaction();

    void setSpeciesList(QStringList _list) {model->setStringList(_list);}
    int getType(void);
    QString getStrReaction(void);

public slots:
    void chooseType(int);

private:
    Ui::ChooseReaction *ui;

    QStringListModel *model;

};

#endif // CHOOSEREACTION_H
