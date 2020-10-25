#ifndef SPECIESTABLEDELEGATE_H
#define SPECIESTABLEDELEGATE_H

#include <QStyledItemDelegate>
#include <QLineEdit>

class SpeciesTableDelegate : public QStyledItemDelegate {

    Q_OBJECT

public:
    SpeciesTableDelegate(QObject *parent = nullptr);

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const override;

};

#endif // SPECIESTABLEDELEGATE_H
