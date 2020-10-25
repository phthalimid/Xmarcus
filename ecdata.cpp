#include "ecdata.h"
//#include "fitchart.h"

#include <QFile>
#include <QTextStream>
#include <QMessageBox>

#include <iostream>

// Helper to round variables.
/*
inline double round(double val) {
    if(val < 0) return ceil(val - 0.5);
    return floor(val + 0.5);
}
*/

ECData::ECData(QWidget *parent) : QWidget(parent) {
    noSweeps = 0;
    scanRate = 0.0;
    dataFile.clear();

}

/*
 * Reads in data file.
 */
void ECData::readFile(QString fileName, int dataSource) {

    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::information(this, tr("Unable to open file"), file.errorString());
        return;
    }

    // Read data...
    int nols = 0;
    QTextStream in(&file);

    switch (dataSource) {
    case 0: // CSV
        while (!in.atEnd()) {
            QString line = in.readLine();
            QStringList fields = line.split(',');
            time.append(fields.at(0).toDouble());
            potential.append(round(fields.at(1).toDouble()*10000.0)/10000.0);
            current.append(fields.at(2).toDouble());
            nols++;
        }
        break;
    case 1: // Zahner

        qDebug() << "Zahner not implemented yet!";

        break;
    case 2: // Autolab
        in.readLine();             // Dump first line.
        while (!in.atEnd()) {
            QString line = in.readLine();
            QStringList fields = line.split(',');
            potential.append(fields[0].toDouble());
            time.append(fields[1].toDouble());
            current.append(fields[2].toDouble());
            nols++;
        }
        nols--; // minus first line.
        break;
    default:
        break;
    }
    file.close();

    // need to clean the data here!
    cleanData();

    // Calculate potential programme parameters.
    int ppsCounter = 2;
    rPotentials.append(potential[0]);
    rTimes.append(time[0]);
    int slope = sgn(potential[1] - potential[0]);
    noSweeps = 1;
    for(int i=2; i<potential.size(); i++) {
        if(sgn(potential[i]-potential[i-1]) != slope) {
            slope *= -1;
            rPotentials.append(potential[i-1]);
            rTimes.append(time[i-1]);
            ppSweep.append(ppsCounter);
            noSweeps++;
            ppsCounter = 0;
        }
        ppsCounter++;
    }
    rPotentials.append(potential.last());
    rTimes.append(time.last());
    ppSweep.append(ppsCounter);
    scanRate = fabs(rPotentials[1] - rPotentials[0]);

    for(int i=2; i<rPotentials.size(); i++)
        scanRate += fabs(rPotentials[i] - rPotentials[i-1]);
    potentialSteps = scanRate/nols;
    scanRate /= (time.last()-time[0]);

    // Round the return potentials.
    for(int i=0; i<rPotentials.size(); i++)
        rPotentials[i] = rPotentials[i] < 0 ? ceil(rPotentials[i]*RDIGS + 0.5)/RDIGS : floor(rPotentials[i]*RDIGS - 0.5)/RDIGS;

}

/*
 * We need to remove all the rubbish from the data.
 * We always take the "last" point in case of multiples.
 * Set first point t = 0.
 */
void ECData::cleanData() {
    for(QVector<double>::iterator it=(potential.begin()+1); it!=potential.end();) {
        if(qFuzzyCompare(*(it-1), *it)) {
            time.remove(int((it-1)-potential.begin()));
            current.remove(int((it-1)-potential.begin()));
            it = potential.erase(it);
        } else
            it++;
    }

    for(int i=(time.size()-1); i>=0; i--)
        time[i] -= time[0];
}

/*
 * Initially fills the 'time' and 'potential' vectors.
 * Was 'CVSweeps()' in old version.
 */
void ECData::calcCVTimePot(void) {

    time.clear();
    potential.clear();
    effPotential.clear();
    double dt = potentialSteps/scanRate;

    qDebug() << "dE = " << potentialSteps << "\t dt = " << dt;

    // Starting time and potential
    time.append(0.0);
    potential.append(rPotentials[0]);

    double slope = sgn(rPotentials[1] - rPotentials[0]);
    for(int i=1; i<=noSweeps; i++) {
        if(slope > 0.0) {
            while((potential.last()+potentialSteps) <= rPotentials[i]) {
                time.append(time.last()+dt);
                potential.append(potential.last()+potentialSteps);
            }
        } else {
            while((potential.last()-potentialSteps) >= rPotentials[i]) {
                time.append(time.last()+dt);
                potential.append(potential.last()-potentialSteps);
            }
        }
        rPotentials[i] = potential.last();
        rTimes[i] = time.last();
        slope *= -1.0;
    }

    effPotential = potential;

    // Read iRC parameters.
    if(xmParent->getFlags() & XMFLAG_IRC) {
        resistance.clear();
        dlcapacitance.clear();
        resistance.append(xmParent->getResistance(0));
        dlcapacitance.append(xmParent->getDLCapacitance(0));
        if(xmParent->getFlags() & XMFLAG_FIT) {
            resistance.append(xmParent->getResistance(1));
            dlcapacitance.append(xmParent->getDLCapacitance(1));
            resistance.append(xmParent->getResistance(2));
            dlcapacitance.append(xmParent->getDLCapacitance(2));
        }
    }

}

/*
 * Interpolate data to match simulation.
 * IMPORTANT: the interpolated data must lie between the endpoints of the experimental data!
 */
void ECData::splineData(QVector<double> simt, QVector<double> simE) {

    splineCurrent.clear();

    double *E = nullptr;
    double *curr = nullptr;

    QVectorIterator<double> itPot(potential);
    QVectorIterator<double> itCurr(current);
    QVectorIterator<double> itSimE(simE);
    QVectorIterator<double> itSimt(simt);

    gsl_spline *spline = nullptr;
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

//    std::ofstream splFile;
//    splFile.open("/Users/nannth/Desktop/splinefile.csv");

    for(int i=1; i<(noSweeps+1); i++) {

//        qDebug() << "ppsweep = " << ppSweep[i-1];

        E = new double[static_cast<unsigned long>(ppSweep[i-1])];
        curr = new double[static_cast<unsigned long>(ppSweep[i-1])];
        for(int j=0; j<ppSweep[i-1]; j++) {
            E[j] = itPot.next();
            curr[j] = itCurr.next();
        }

        // Save sweep direction in sign and reverse if necessary for splineing.
        int sign = sgn(E[1] - E[0]);
        if(sign == -1) {
            std::reverse(E, &E[ppSweep[i-1]]);
            std::reverse(curr, &curr[ppSweep[i-1]]);
        }

        spline = gsl_spline_alloc(gsl_interp_cspline, static_cast<size_t>(ppSweep[i-1]));
        gsl_spline_init(spline, E, curr, static_cast<size_t>(ppSweep[i-1]));

//        qDebug() << "rtimes = " << rTimes[i];

        // first one by hand.
        splineCurrent.append(gsl_spline_eval(spline, itSimE.next(), acc));
//        splFile << itSimE.peekPrevious() << ",\t" << splineCurrent.last() << "\n";

        while((sgn(itSimE.peekNext()-itSimE.peekPrevious()) == sign) && itSimE.hasNext()) {
            splineCurrent.append(gsl_spline_eval(spline, itSimE.next(), acc));
//            qDebug() << "E = " << itSimE.peekPrevious() << "\t i = " << splineCurrent.last();
//            splFile << itSimE.peekPrevious() << ",\t" << splineCurrent.last() << "\n";
        }

        free(E);
        free(curr);
        gsl_spline_free(spline);
    }
    gsl_interp_accel_free(acc);

//    splFile.close();
}

void ECData::clearData(void) {
    time.clear();
    potential.clear();
    effPotential.clear();
    current.clear();
    splineCurrent.clear();
    noSweeps = 1;
    scanRate = potentialSteps = 0.0;
    rPotentials.clear();
    rTimes.clear();
    ppSweep.clear();
    temperature = 298.15;
    area.clear();
    resistance.clear();
    dlcapacitance.clear();
    dataFile.clear();
}
