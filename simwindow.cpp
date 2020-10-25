#include "simwindow.h"
#include "ui_simwindow.h"

SimWindow::SimWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SimWindow)
{
    ui->setupUi(this);

    connect(ui->pushButton_SaveSim, SIGNAL(clicked()), this, SLOT(saveSim()));
    connect(ui->pushButton_CancelSim, SIGNAL(clicked()), this, SLOT(close()));
    connect(ui->pushButton_autoFill, SIGNAL(clicked(bool)), this, SLOT(emitAutoFillSignal()));

}

SimWindow::~SimWindow() {
    delete ui;
}

void SimWindow::emitAutoFillSignal() {
    emit autoFillClicked();
}

/*
 * Updates the simulation window
 */
void SimWindow::updateSimWindow() {

    ui->customPlot_CVs->addGraph();
    if(fitData) {
        ui->customPlot_CVs->graph(0)->setPen(QPen(Qt::blue));
        ui->customPlot_CVs->graph(0)->setLineStyle(QCPGraph::lsNone);
        ui->customPlot_CVs->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
        ui->customPlot_CVs->graph(0)->setData(fitData->getPotential(), fitData->getCurrent());
   }

    ui->customPlot_CVs->addGraph();
    if(simData) {
        ui->customPlot_CVs->graph(1)->setPen(QPen(Qt::red));
        ui->customPlot_CVs->graph(1)->setLineStyle(QCPGraph::lsNone);
        ui->customPlot_CVs->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
        ui->customPlot_CVs->graph(1)->setData(simData->getPotential(), simData->getCurrent());
    }

    ui->customPlot_CVs->xAxis->setLabel("E / [V]");
    ui->customPlot_CVs->yAxis->setLabel("i / [A]");
    ui->customPlot_CVs->rescaleAxes();
    ui->customPlot_CVs->replot();
}


void SimWindow::plotGraph(QVector<double> x, QVector<double> y, box wbox) {

    ui->customPlot_CVs->addGraph();
    ui->customPlot_CVs->graph(0)->setData(x, y);
    ui->customPlot_CVs->graph(0)->setPen(QColor(204, 0, 0, 255));
    ui->customPlot_CVs->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->customPlot_CVs->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));

    ui->customPlot_CVs->xAxis->setRange(wbox.xmin, wbox.xmax);
    ui->customPlot_CVs->xAxis->setLabel("E / [V]");
    ui->customPlot_CVs->yAxis->setRange(wbox.ymin, wbox.ymax);
    ui->customPlot_CVs->yAxis->setLabel("i / [A]");

    ui->customPlot_CVs->replot();
}

void SimWindow::saveSim() {

    QString simFile = QFileDialog::getSaveFileName(this, tr("Save Simulation"), "", tr("comma separated values (*.csv);;All Files (*)"));

    if (simFile.isEmpty())
        return;
    else {
        QFile file(simFile);
        if (!file.open(QIODevice::WriteOnly)) {
            QMessageBox::information(this, tr("Unable to open file"), file.errorString());
            return;
        }

        QTextStream out(&file);

        for(int i=0; i<simData->getTime().size(); i++)
            out << simData->getTime(i) << ", " << simData->getPotential(i) << ", " << simData->getCurrent(i) << Qt::endl;

        file.close();
    }

    qDebug() << "saveSim(): ... was here ...";

}
