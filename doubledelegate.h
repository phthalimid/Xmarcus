#ifndef DOUBLEDELEGATE_H
#define DOUBLEDELEGATE_H

#include <QItemDelegate>
#include <QLineEdit>
#include <QDoubleValidator>

#include "xmarcus.h"

class DoubleDelegate : public QItemDelegate {

    Q_OBJECT

public:
    DoubleDelegate();

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &,
                             const QModelIndex &index) const override;

};

#endif // DOUBLEDELEGATE_H
