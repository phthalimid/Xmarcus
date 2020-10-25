#ifndef REACTION_H
#define REACTION_H

#include <QString>
#include <QVector>

#define RCTFLAG_E  0b0000000001   // Electrochem. reaction
#define RCTFLAG_C  0b0000000010   // First order homogeneous
#define RCTFLAG_C2 0b0000000100   // Second order homogeneous
#define RCTFLAG_CD 0b0000001000   // Homogeneous disproportionation
#define RCTFLAG_CC 0b0000010000   // Homogeneous comproportionation

class Reaction {

public:
    Reaction();

    void setFit(bool _fit) {fit = _fit;} bool getFit(void) {return fit;}
    void setType(int _type) {type = _type;} int getType(void) {return type;}
    void setStrReaction(QString _strReaction) {strReaction = _strReaction;} QString getStrReaction(void) {return strReaction;}

    void setRctEPara(int _i, double _value) {rctEPara[_i] = _value;} double getRctEPara(int _i) {return rctEPara[_i];}
    void setRctCPara(int _i, double _value) {rctCPara[_i] = _value;} double getRctCPara(int _i) {return rctCPara[_i];}

private:
    int type = 0;
    QString strReaction;

    bool fit = false;         // fit parameters

    QVector<double> rctEPara = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // NEW reaction parameters.
    QVector<double> rctCPara = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

};

#endif // REACTION_H
