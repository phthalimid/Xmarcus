#ifndef ECDATA_H
#define ECDATA_H

#include "xmarcus.h"

#include <QWidget>
#include <QString>
#include <QVector>

#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <gsl/gsl_spline.h>

//#define ECDFLAG_AREA 0b00001    // Fit area
//#define ECDFLAG_RES  0b00010    // Fit resistance
//#define ECDFLAG_DLC  0b00100    // Fit capacitance

#define RDIGS 10000.0 // Determines digits for rounding the potential

class Xmarcus;

class ECData: public QWidget {
    Q_OBJECT

public:
    ECData(QWidget *parent = nullptr);

    void setTemperature(double _T) {temperature = _T;} double getTemperature(void) {return temperature;}

    void setXMParent(Xmarcus *_xmParent) {xmParent = _xmParent;}
    QString getDataFile() {return dataFile;} void setDataFile(QString _file) {dataFile = _file;}
    box getDataRange(void) {return dataRange;}
    int getIData(void) {return iData;} void setIData(int _iData) {iData = _iData;}
    double getPotentialSteps() {return potentialSteps;} void setPotentialSteps(double _value) {potentialSteps = _value;}
    double getScanRate() {return scanRate;} void setScanRate(double _scanRate) {scanRate = _scanRate;}
    double getTime(int i) {return time.at(i);} QVector<double> getTime() {return time;}
    double getPotential(int i) {return potential.at(i);} QVector<double> getPotential() {return potential;}
    double getCurrent(int i) {return current.at(i);} QVector<double> getCurrent() {return current;} void appendCurrent(double _value) {current.append(_value);} void clearCurrent() {current.clear();}
    double getSplineCurrent(int i) {return splineCurrent[i];}
    double getEffPotential(int i) {return effPotential.at(i);} void setEffPotential(int _i, double _Eff) {effPotential.replace(_i, _Eff);}
    double getArea(int i) {return area.at(i);} QVector<double> getArea() {return area;} void appendArea(double _value) {area.append(_value);} void setArea(int _i, double _value) {area.replace(_i, _value);}
    double getResistance(int i) {return resistance.at(i);} QVector<double> getResistance() {return resistance;} void appendResistance(double _value) {resistance.append(_value);} void setResistance(int _i, double _value) {resistance.replace(_i, _value);}
    double getDLCapacitance(int i) {return dlcapacitance.at(i);} QVector<double> getDLCapacitance() {return dlcapacitance;} void appendDLCapacitance(double _value) {dlcapacitance.append(_value);} void setDLCapacitance(int _i, double _value) {dlcapacitance.replace(_i, _value);}
    int getNoSweeps() {return noSweeps;} void setNoSweeps(int _n) {noSweeps = _n;}
    double getRPotentials(int i) {return rPotentials[i];} QVector<double> getRPotentials() {return rPotentials;} void appendRPotentials(double _value) {rPotentials.append(_value);} void setRPotentials(QVector<double> _vector) {rPotentials = _vector;}
    double getRTimes(int i) {return rTimes[i];} QVector<double> getRTimes() {return rTimes;} void appendRTimes(double _value) {rTimes.append(_value);} void setRTimes(QVector<double> _vector) {rTimes = _vector;}

    void readFile(QString, int);
    void clearData(void);
    void calcCVTimePot(void);
    void splineData(QVector<double>, QVector<double>);

private:
    void cleanData();

    Xmarcus *xmParent = nullptr;

    int iData = 0;       // current datapoint.
    QString dataFile;

    double potentialSteps;
    double scanRate;
    QVector<double> time;
    QVector<double> potential;
    QVector<double> effPotential;
    QVector<double> current;
    QVector<double> splineCurrent;

    int noSweeps;
    QVector<int> ppSweep;
    QVector<double> rTimes;
    QVector<double> rPotentials;

    box dataRange;

    double temperature;

    QVector<double> area;
    QVector<double> resistance;
    QVector<double> dlcapacitance;

};

#endif // ECDATA_H
