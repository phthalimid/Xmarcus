#include "reactionstabledelegate.h"

ReactionsTableDelegate::ReactionsTableDelegate(QObject *parent) : QStyledItemDelegate(parent) {

}

QWidget *ReactionsTableDelegate::createEditor(QWidget *parent,
    const QStyleOptionViewItem &/* option */,
    const QModelIndex &index) const {

    QLineEdit *lineEdit = new QLineEdit(parent);

    if(index.column() > 2) { // only beyond second column
        // Set validator
        QDoubleValidator *validator = new QDoubleValidator(lineEdit);
        lineEdit->setValidator(validator);
    }

    return lineEdit;
}
