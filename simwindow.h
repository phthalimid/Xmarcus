#ifndef SIMWINDOW_H
#define SIMWINDOW_H

#include <QWidget>

#include "xmarcus.h"

class ECData;

namespace Ui {
class SimWindow;
}

class SimWindow : public QWidget
{
    Q_OBJECT

public:
    explicit SimWindow(QWidget *parent = nullptr);
    ~SimWindow();

    void setSimData(ECData *_simData) {simData = _simData;}
    void setFitData(ECData *_fitData) {fitData = _fitData;}

    void plotGraph(QVector<double>, QVector<double>, box);

public slots:
    void saveSim();
    void updateSimWindow();

signals:
    void autoFillClicked(void);

private slots:
    void emitAutoFillSignal(void);

private:
    Ui::SimWindow *ui;

    ECData *simData = nullptr;
    ECData *fitData = nullptr;

};

#endif // SIMWINDOW_H
