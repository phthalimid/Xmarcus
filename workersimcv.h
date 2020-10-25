#ifndef WORKERSIMCV_H
#define WORKERSIMCV_H

#include "xmarcus.h"
#include "ecdata.h"
#include "species.h"
#include "reaction.h"

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>

#include <QObject>

#define LDEGREE 7	// Degree of the Lagrange polynomial for current calculation.
#define INITSS 10   // Fraction of the total range for fitting: initial step size.
#define IRCLIMIT 3  // Multiples of stepsize to calculate lower and higher limit in case of iRC.

#define FITFLAG_AREA  0b00000000001  // Fit area
#define FITFLAG_RES   0b00000000010  // Fit resistance
#define FITFLAG_DLC   0b00000000100  // Fit double-layer capacitance
#define FITFLAG_REACT 0b00000001000  // Fit reaction
#define FITFLAG_QMAX  0b00000010000  // Fit q_max
#define FITFLAG_ADS   0b00000100000  // Fit adsorption

struct mh_params { // Parameters needed for the Marcus-Hush calculation (Feldberg method).
    double	T;
    double	Estar;
    double	lambdastar;
};

class ECData;
class Xmarcus;

/*
 * Followed these instructions: https://mayaposch.wordpress.com/2011/11/01/how-to-really-truly-use-qthreads-the-full-explanation/
 */

class WorkerSimCV : public QObject {
    Q_OBJECT

public:
    explicit WorkerSimCV(QObject *parent = nullptr);
    ~WorkerSimCV();

    void setActive(bool _active) {active = _active;}

    void setXmParent(Xmarcus *_xmParent) {xmParent = _xmParent;}
    void setKillThread(bool _killThread) {killThread = _killThread;}
    gsl_matrix *getLinm() {return linm;}
    gsl_vector *getC() {return c;}
    int getNop() {return nop;}
    double getDt() {return dt;}
    gsl_matrix_int *getM() {return M;}
    void setWorkAmount(int _value) {workAmount = _value;}

signals:
    void finished();
    void error(QString);
    void sendThreadStatus(QString);
    void sendProgressSim(int);      // Emits the progress of a single CV simulation (datapoint).
    void sendFinishedCV(ECData *);

public slots:
    void process();
    void startSim();    // Activates simulation thread.

private:
    bool active;
    int workAmount;
    int workDone;
    QMutex lock;

    void setThreadStatus(QString);

    WorkerSimCV *whoami;
    bool killThread = false;

    void initialiseAll();
    void updateNumParams();
    void resetConcVector();

    // Calculates the x position according to: A. W. Bott, Current Separations 2000, 19, 45–48.
    double xi(int i) { return dx*(exp(beta*(i+0.5))-1.0)/(exp(beta)-1.0); } // i starts at 1 in the paper.

    static double intgl1 (double, void *);
    static double intgl2 (double, void *);
    void calcHeteroKinetics(double);
    double calcFCurrent();
    double calcCCurrent();
    double calcAdsCurrent();
    void solveSSE(double);
    static int CV_f(const gsl_vector *, void *, gsl_vector *);
    static double iRC_f(double, void *);

    void procCV();      // Simulate CV
    void procFitCV();   // Fit CV

    Xmarcus *xmParent = nullptr;
    ECData *simData;
    ECData *fitData;
    QVector<Species *> allSpecies;
    QVector<Reaction *> allReactions;

    gsl_matrix_int *M = nullptr;

    // Stuff related to the grid.
    double Dmax, xmax, dx, dt;
    int nop;                        // Number of grid points.

    // Stuff we need to solve the matrix.
    gsl_vector *c = nullptr;        // Concentrations vector.
    //gsl_vector *chelp = nullptr;    // Vector to rescue concentrations in case of iRC-drop.
    gsl_vector *alphai = nullptr;
    gsl_vector *betai = nullptr;
    gsl_matrix *m = nullptr;
    gsl_matrix *linm = nullptr;

    gsl_vector *oldThetai = nullptr;    // Stores "old" theta_i for current calculation.

    gsl_vector *zetaj = nullptr;    // Lagrange factors.

    // Stuff we need for fitting.
    int fitFlags = 0;
    int fitNoVariables = 0;

    static double FitCV_f(const gsl_vector *, void *);
    gsl_vector *calcFitVariables();
    void updateFitVariables(const gsl_vector *);
    gsl_vector *setInitStepSize();
    void resultFitVariables(gsl_vector *);

    // Numerical variables.
    double beta = 0.2;          //	0.2		# exponential grid parameter [usually 0.2]
    double Dm = 10.0;           //	10		# funny DigiSim parameter to calculate delta x
    double nlerror = 1.0e-8;    //	1.0e-8	# non-linear: absolute error bound
    int maxiter = 50000;        //	50000	# non-linear: maximum number of iterations

};

#endif // WORKERSIMCV_H
