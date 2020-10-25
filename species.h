#ifndef SPECIES_H
#define SPECIES_H

#include <QString>

class Species
{
public:
    Species();
    Species(QString);
    Species(double, double, int, double, double);

    void setAdsorb(bool _adsorb) {adsorb = _adsorb;}
    bool getAdsorb(void) {return adsorb;}
    void setSymbol(QString _symbol) {symbol = _symbol;}
    QString getSymbol(void) {return symbol;}
    void setCini(double _cini) {cini = _cini;}
    double getCini(void) {return cini;}
    void setD(double _D) {D = _D;}
    double getD(void) {return D;}
    void setCharge(int _charge) {charge = _charge;}
    int getCharge(void) {return charge;}
    void setKa(double _ka) {ka = _ka;}
    double getKa(void) {return ka;}
    void setKd(double _kd) {kd = _kd;}
    double getKd(void) {return kd;}

private:
    bool adsorb;
    QString symbol;
    double cini;
    double D;
    int charge;
    double ka, kd;

};

#endif // SPECIES_H
