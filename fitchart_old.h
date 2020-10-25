#ifndef FITCHART_H
#define FITCHART_H

#include "xmarcus.h"
#include "ecdata.h"

#include <QWidget>
#include <QMainWindow>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>

class Xmarcus;

namespace Ui {
class fitChart;
}

class fitChart : public QMainWindow {
    Q_OBJECT

public:
    explicit fitChart(QWidget *parent = nullptr);

    Ui::fitChart *ui;

    Xmarcus *mainWindow;

    void plotGraph(QVector<double>, QVector<double>, box);

signals:
    void autoFillClicked(void);

private slots:
    void emitAutoFillSignal(void);

private:


};

#endif // FITCHART_H
