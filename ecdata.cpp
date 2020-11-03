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
    case 3: // EClab, BioLogic
        in.readLine();             // Dump first line.
        while (!in.atEnd()) {
            QString line = in.readLine();
            QStringList fields = line.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
            potential.append(fields[0].toDouble());
            current.append(fields[1].toDouble()/1000.0);
            time.append(fields[2].toDouble());
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
/*
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
*/

}

/*
 * We need to remove all the rubbish from the data.
 *
 * We always take the "last" point in case of multiples.
 * Set first point t = 0.
 *
 * 1) Remove unwanted observations
 * 2) Fix structural errors
 * 3) Manage outliers
 * 4) Handle missing data
 *
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

    // 1) remove unwanted observations (data with no slope at beginning and end)

    _findSegments(50);

}


/*
 * Find segments in the provided potential data.
 * Use this idea: https://www.codeproject.com/Articles/5282014/Segmented-Linear-Regression
 *
 */
void ECData::_findSegments(int window) {

    // Make sure window is odd number (symetric).
    if (window % 2 == 0)
        window += 1;

    gsl_vector *gslPotential = gsl_vector_alloc (potential.size());
    gsl_vector *gslSmaDiffPot = gsl_vector_alloc (potential.size());
    for (int i = 0; i < potential.size(); ++i)
        gsl_vector_set (gslPotential, i, potential.at(i));

    // Calculate (symmetric) moving average from potentials vector
    double sma = 0.0;
    int halfwindow = floor(window/2);
    for (int i = 0; i < halfwindow; i++)
        sma += potential.at(i);
    for (int i = 0; i < potential.size(); ++i) {
        if (i <= halfwindow) {   // Head of data.
            sma += potential.at(i+halfwindow);
            gsl_vector_set (gslSmaDiffPot, i, sma/(i+1+halfwindow));
        } else if ((i > halfwindow) && (i < potential.size()-halfwindow)) {
            sma -= potential.at(i-halfwindow-1);
            sma += potential.at(i+halfwindow);
            gsl_vector_set (gslSmaDiffPot, i, sma/window);
        } else if (i >= potential.size()-halfwindow) {
            sma -= potential.at(i-halfwindow-1);
            gsl_vector_set (gslSmaDiffPot, i, sma/(potential.size()-i+halfwindow));
        }
    }

    // Calculate difference between SMA and data vectors
    gsl_vector_sub(gslSmaDiffPot, gslPotential);

    FILE *f = fopen ("/Users/tn438/Desktop/test.txt", "wb");
    gsl_vector_fprintf (f, gslSmaDiffPot, "%g");
    fclose (f);

    // Find turning points.
    int extreme = 0;
    gsl_vector_view subv;
    QVector<int> idxTurns;
    idxTurns.append(0); // add frist point
    noSweeps = 0;
    for (int i = 0; i < potential.size(); ++i) {
        if (abs(gsl_vector_get (gslSmaDiffPot, i)) >= 0.003) {      // turning point found.
            subv = gsl_vector_subvector(gslSmaDiffPot, i, window);
            if (gsl_vector_get (gslSmaDiffPot, i) > 0.0) {          // minimum found.
                extreme = gsl_vector_max_index(&subv.vector) + i;
            } else if (gsl_vector_get (gslSmaDiffPot, i) < 0.0) {   // maximum found.
                extreme = gsl_vector_min_index(&subv.vector) + i;
            }
            idxTurns.append(extreme);
            noSweeps++;
            i += window;
        }
    }
    idxTurns.append(potential.size()-1);    // add last data point
    noSweeps++;

    // Calculate segments
    vector<double> vc0;
    vector<double> vc1;

    for (int i=0; i<noSweeps; i++) {
        vector<double> x;

        int idxcut = floor((idxTurns[i+1] - idxTurns[i])/20.0);  // use 90% of segment for linear regression.
        int idxstart = idxTurns[i] + idxcut;
        int idxstop = idxTurns[i+1] - idxcut;
        int seglength = idxstop - idxstart;

        x.resize(seglength);
        iota(x.begin(), x.end(), idxstart); // range()
        double *y = gsl_vector_ptr(gslPotential, idxstart);

        double c0, c1, cov00, cov01, cov11, sumsq;
        gsl_fit_linear(x.data(), 1, y, 1, seglength, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
        vc0.insert(vc0.end(), c0);
        vc1.insert(vc1.end(), c1);

        // Correct idx turns
        if (i>0) {
            if (isnan(vc0[i-1]))  // For some odd reason the first gsl_fit_linear is sometimes returning NaN
                idxTurns[i] = floor((vc0[i] - potential.at(i-1))/(-vc1[i]) - 1);
            else
                idxTurns[i] = floor((vc0[i] - vc0[i-1])/(vc1[i-1] - vc1[i]));
        }
        x.clear();
    }

    // Smooth out potentials using the linear regression
    int segment = 0;
    for (int i=0; i<potential.size(); i++) {
        if (i>idxTurns[segment+1])
            segment++;
        potential[i] = vc0[segment] + vc1[segment] * i;
    }

    // Cut head off if necessary
    if (idxTurns[1] <= window) {
        potential.remove(potential.front(), idxTurns[1]);
        time.remove(time.front(), idxTurns[1]);
        double time_zero = time[0];
        for(int i=0; i<time.size(); i++)
            time[i] -= time_zero;
        current.remove(current.front(), idxTurns[1]);
        for(int i=2; i<idxTurns.size(); i++)
            idxTurns[i] -= idxTurns[1];
        idxTurns.removeFirst();
        idxTurns[0] = 0;
        vc0.erase(vc0.begin());
        vc1.erase(vc1.begin());
        noSweeps--;
    }

    // Cut tail off if necessary
    if (vc1[vc1.size()-1] < 1.0e-5) {
        potential.remove(idxTurns.end()[-2]+1, idxTurns.end()[-1]-idxTurns.end()[-2]);
        time.remove(idxTurns.end()[-2]+1, idxTurns.end()[-1]-idxTurns.end()[-2]);
        current.remove(idxTurns.end()[-2]+1, idxTurns.end()[-1]-idxTurns.end()[-2]);
        idxTurns.removeLast();
        vc0.pop_back();
        vc1.pop_back();
        noSweeps--;
    }

    // Correct turning points
    for (int i=0; i<noSweeps-1; i++) {
        int count = idxTurns[i+1]-halfwindow;
        int max_pot = count;
        if (potential[count+1] - potential[count] > 0.0) {  // maximum
            for (int j=1; j<window; j++) {
                if (potential[count+j]>potential[count+j-1])
                    max_pot = count+j;
            }
        } else if (potential[count+1] - potential[count] < 0.0) {  // minimum
            for (int j=1; j<window; j++) {
                if (potential[count+j]<potential[count+j-1])
                    max_pot = count+j;
            }
        }
        idxTurns[i+1] = max_pot;
    }



/*
    for(int i=0; i < vc0.size(); i++)
       cout << vc0.at(i) << ' ';
    cout << endl;

    for(int i=0; i < vc1.size(); i++)
       cout << vc1.at(i) << ' ';
    cout << endl;
*/

    // Fill in the return potentials and points per sweep
    rPotentials.clear();
    rTimes.clear();
    ppSweep.clear();

    rPotentials.append(potential[0]);
    rTimes.append(time[0]);
    for (int i=1; i<idxTurns.size(); i++) {
        rPotentials.append(potential[idxTurns[i]]);
        rTimes.append(time[idxTurns[i]]);
        ppSweep.append(idxTurns[i]-idxTurns[i-1]);
    }

    // Round the return potentials.
    for(int i=0; i<rPotentials.size(); i++)
        rPotentials[i] = rPotentials[i] < 0 ? ceil(rPotentials[i]*RDIGS + 0.5)/RDIGS : floor(rPotentials[i]*RDIGS - 0.5)/RDIGS;

    // Calculate average scanrate
    scanRate = 0.0;
    for (int i=0; i<noSweeps; i++)
        scanRate += abs(vc1[i]);
    scanRate /= noSweeps;   // this is actually dE, because our timestep was 1

    potentialSteps = round(scanRate*RDIGS)/RDIGS;
    scanRate = round(potentialSteps/(time.end()[-1]/time.size())*RDIGS)/RDIGS;

    int count = 0;
    for (int i=0; i<noSweeps; i++) {
        int direction = rPotentials[i+1] > rPotentials[i] ? 1 : -1;
        for (int j=0; j<ppSweep[i]-1; j++) {
            if (direction*(potential[count+1]-potential[count]) <= 0.0)
                cout << "problem!!!" << endl;
            count++;
        }
        count++;
    }

    gsl_vector_free (gslPotential);
    gsl_vector_free (gslSmaDiffPot);

}


/*
 * Initially fills the 'time' and 'potential' vectors.
 * Was 'CVSweeps()' in old version.
 */
void ECData::calcCVTimePot(void) {

    time.clear();

//    potential not same number of points than below

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
void ECData::splineData(QVector<double> simE) {
//void ECData::splineData(QVector<double> simt, QVector<double> simE) {

    splineCurrent.clear();

    cout << simE.size() << endl << flush;

/*
    double *E = nullptr;
    double *curr = nullptr;

    QVectorIterator<double> itPot(potential);
    QVectorIterator<double> itCurr(current);
    QVectorIterator<double> itSimE(simE);
    QVectorIterator<double> itSimt(simt);
*/
    gsl_spline *spline = nullptr;
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

//    std::ofstream splFile;
//    splFile.open("/Users/nannth/Desktop/splinefile.csv");

    int data_pos = 0;
    for (int i=0; i<noSweeps; i++) {

        vector<double> x(&potential.data()[data_pos], &potential.data()[data_pos]+ppSweep[i]+1);
        vector<double> y(&current.data()[data_pos], &current.data()[data_pos]+ppSweep[i]+1);

//        cout << "size x = " << x.size() << "\t ppSweep = " << ppSweep[i] << endl << flush;

        // Reverse if necessary
        if (x[1] - x[0] < 0.0) {
            reverse(x.begin(), x.end());
            reverse(y.begin(), y.end());
        }

        spline = gsl_spline_alloc(gsl_interp_cspline, x.size());
        gsl_spline_init(spline, x.data(), y.data(), x.size());

        // Add to spline current vector
        if (i==0)   // First datapoint
            splineCurrent.append(gsl_spline_eval(spline, simE[0], acc));
        for (int j=1; j<int(x.size()); j++)
            splineCurrent.append(gsl_spline_eval(spline, simE[data_pos+j], acc));

        data_pos += ppSweep[i];

        gsl_spline_free(spline);
    }
    gsl_interp_accel_free(acc);


/*
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
*/

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
