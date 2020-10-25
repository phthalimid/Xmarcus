#include "speciestabledelegate.h"

SpeciesTableDelegate::SpeciesTableDelegate(QObject *parent) : QStyledItemDelegate(parent) {

}

QWidget *SpeciesTableDelegate::createEditor(QWidget *parent,
    const QStyleOptionViewItem &/* option */,
    const QModelIndex &index) const {

    QLineEdit *lineEdit = new QLineEdit(parent);

    if(index.column() > 1) { // only beyond second column
        // Set validator
        QDoubleValidator *validator = new QDoubleValidator(lineEdit);
        lineEdit->setValidator(validator);
    }

    return lineEdit;
}
