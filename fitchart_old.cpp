#include "fitchart.h"
#include "ui_fitchart.h"

fitChart::fitChart(QWidget *parent) : QMainWindow(parent), ui(new Ui::fitChart) {

    ui->setupUi(this);

    QObject::connect(ui->pushButton, SIGNAL(clicked(bool)), this, SLOT(emitAutoFillSignal()));

}

void fitChart::emitAutoFillSignal() {
    emit autoFillClicked();
}

void fitChart::plotGraph(QVector<double> x, QVector<double> y, box wbox) {

    ui->customPlot->addGraph();
    ui->customPlot->graph(0)->setData(x, y);
    ui->customPlot->graph(0)->setPen(QColor(102, 153, 51, 255));
    ui->customPlot->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->customPlot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));

    ui->customPlot->xAxis->setRange(wbox.xmin, wbox.xmax);
    ui->customPlot->xAxis->setLabel("E / [V]");
    ui->customPlot->yAxis->setRange(wbox.ymin, wbox.ymax);
    ui->customPlot->yAxis->setLabel("i / [A]");

    ui->customPlot->replot();

}
