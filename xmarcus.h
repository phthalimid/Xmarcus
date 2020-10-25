#ifndef XMARCUS_H
#define XMARCUS_H

#define VERSION 0.7

#include <QMainWindow>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QStandardItemModel>
#include <QFlags>
//#include <QFuture>
#include <QtConcurrent/QtConcurrent>

#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/*
 * Little helpers.
 */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct box {
    double xmin, xmax, ymin, ymax;
};

#include "ecdata.h"
#include "species.h"
#include "reaction.h"
//#include "fitchart.h"
#include "numericalparams.h"
#include "choosereaction.h"
#include "simwindow.h"
#include "workersimcv.h"

#include "doubledelegate.h"
#include "speciestabledelegate.h"
#include "reactionstabledelegate.h"

//#define XMFLAG_NOSIMUPD 0b0000000001
//#define XMFLAG_NOCLRFIT 0b0000000010 // Do not erase fitData when loading fitfile.
#define XMFLAG_FIT      0b0000000100
#define XMFLAG_IRC      0b0000001000
#define XMFLAG_AIRCUPD  0b0000010000
#define XMFLAG_BV       0b0000100000
#define XMFLAG_MH       0b0001000000
#define XMFLAG_ADS      0b0010000000 // At lease one substance adsorbing...
#define XMFLAG_NL       0b0100000000 // Non-linear problem
#define XMFLAG_FRCNL    0b1000000000 // Force non-linear solver

class ECData;
class SimWindow;
class WorkerSimCV;

namespace Ui {
class Xmarcus;
}

class Xmarcus : public QMainWindow {
    Q_OBJECT

public:
    explicit Xmarcus(QWidget *parent = nullptr);
    ~Xmarcus();

    void setFlags(int _flags) {flags = _flags;} int getFlags() {return flags;}
    QVector<Species*> getAllSpecies(void) {return allSpecies;}
    QVector<Reaction*> getAllReactions(void) {return allReactions;}
    ECData *getSimData() {return simData;}
    ECData *getFitData() {return fitData;}
    gsl_matrix_int *getM() {return M;}
    WorkerSimCV *getWorker() {return worker;}
    double getResistance(int);
    double getDLCapacitance(int);
    SimWindow *getSimWindow() {return simWindow;}

    NumericalParams *getNumParams() {return &numParams;}
    void setQmax(double _qmax) {qmax = _qmax;} double getQmax(void) {return qmax;}
    void setQmax(void);

    bool isAreaFitChecked();
    bool isResFitChecked();
    bool isDLCapFitChecked();
    bool isReactFitChecked(int);
    void setUiArea(double);
    void setUiResistance(double);
    void setUiDLCapacitance(double);
    void setUiReactParam(int, int, double);

public slots:
    void openFile();
    void saveFile();
    void saveAsFile();
    void startSimulation();
    void modelChanged(int);
    void loadFitFile();
    void autoFillFitData();
    void fittingChanged(bool);
    void sweepsChanged(int);
    void scanRateChanged(double);
    void tablePotChanged(void);
    void iRCChanged(bool);
    void checkNumParams(void);
    void tableSpeciesChanged(void);
    void speciesChanged(int);
    void tableReactionsClicked(const QModelIndex &);
    void tableReactionsChanged(void);
    void reactionsChanged(int);
    void checkReactions(void);
    void errorString(QString);
    void receiveProgressSim(int);
    void receiveFinishedCV(ECData *);

protected:

private:
    int flags = 0;

    Ui::Xmarcus *ui;
    QThread *thread;    // Thread to do the simulation.

    void setupConnections();

    NumericalParams numParams;
    ChooseReaction chooseReaction;

    SimWindow *simWindow = nullptr;

    QString confFile;
    QString workingDir;

    ECData *fitData = nullptr;
    ECData *simData = nullptr;

    void updatePotPlot(double t = -1.0, double E = 0.0);

    void cpRPots2Table();
    void updateSimData(void);
    void checkNonLinear(void);

    QVector<Species*> allSpecies;
    QVector<Reaction*> allReactions;

    // Reactions matrix
    gsl_matrix_int *M = nullptr;

    WorkerSimCV *worker= nullptr;

    double qmax;

signals:
    void startSimThread();     // Activate simulation

};

#endif // XMARCUS_H
