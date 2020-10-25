#include "species.h"

Species::Species() {
    adsorb = false;
    symbol = "xyz";
    cini = 0.0;
    D = 1.0e-9;
    charge = 1;
    ka = 0.0;
    kd = 0.0;
}

Species::Species(QString _symbol) {
    adsorb = false;
    symbol = _symbol;
    cini = 0.0;
    D = 1.0e-9;
    charge = 1;
    ka = 0.0;
    kd = 0.0;
}

Species::Species(double c, double diff, int chrg, double kads, double kdes) {
    cini = c;
    D = diff;
    charge = chrg;
    ka = kads;
    kd = kdes;
}
