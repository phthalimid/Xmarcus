#ifndef REACTIONSTABLEDELEGATE_H
#define REACTIONSTABLEDELEGATE_H

#include <QStyledItemDelegate>
#include <QLineEdit>

class ReactionsTableDelegate : public QStyledItemDelegate {

    Q_OBJECT

public:
    ReactionsTableDelegate(QObject *parent = nullptr);

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const override;

};

#endif // REACTIONSTABLEDELEGATE_H
