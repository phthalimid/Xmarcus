class ECData;

#include "xmarcus.h"
#include "ui_xmarcus.h"

//fitChart *fitWindow;

Xmarcus::Xmarcus(QWidget *parent) : QMainWindow(parent), ui(new Ui::Xmarcus) {
    ui->setupUi(this);

    workingDir.clear();
    fitData = nullptr;
    simData = new ECData(); // this is where the simulation goes.
    simData->setXMParent(this);

    simWindow = new SimWindow();

    // Disable fitting (initially)
    ui->checkBox->setChecked(false);

    // Customise potentials table
    ui->tableWidget->setColumnWidth(0, 80);
    ui->tableWidget->setItemDelegate(new DoubleDelegate());

    // Initial settings for AiRC table.
    ui->tableWidget_AiRC->setColumnWidth(0, 25);
    QWidget *pWidget = new QWidget();
    QCheckBox *pCheckBox = new QCheckBox();
    QHBoxLayout *pLayout = new QHBoxLayout(pWidget);
    pLayout->addWidget(pCheckBox); pLayout->setAlignment(Qt::AlignCenter); pLayout->setContentsMargins(0,0,0,0); pWidget->setLayout(pLayout);
    ui->tableWidget_AiRC->setCellWidget(0, 0, pWidget);
    pWidget = new QWidget(); pCheckBox = new QCheckBox(); pLayout = new QHBoxLayout(pWidget);
    pLayout->addWidget(pCheckBox); pLayout->setAlignment(Qt::AlignCenter); pLayout->setContentsMargins(0,0,0,0); pWidget->setLayout(pLayout);
    ui->tableWidget_AiRC->setCellWidget(1, 0, pWidget);
    pWidget = new QWidget(); pCheckBox = new QCheckBox(); pLayout = new QHBoxLayout(pWidget);
    pLayout->addWidget(pCheckBox); pLayout->setAlignment(Qt::AlignCenter); pLayout->setContentsMargins(0,0,0,0); pWidget->setLayout(pLayout);
    ui->tableWidget_AiRC->setCellWidget(2, 0, pWidget);
    ui->tableWidget_AiRC->setItemDelegate(new DoubleDelegate());

    QDoubleValidator *dv = new QDoubleValidator(0.0, 5000.0, 3);
    ui->lineEdit_temp->setValidator(dv);

    // Species sub-window.
    ui->tableWidget_Species->setColumnWidth(0, 25);
    ui->tableWidget_Species->setColumnWidth(1, 50);
    pWidget = new QWidget(); pCheckBox = new QCheckBox(); pLayout = new QHBoxLayout(pWidget);
    pLayout->addWidget(pCheckBox); pLayout->setAlignment(Qt::AlignCenter); pLayout->setContentsMargins(0,0,0,0); pWidget->setLayout(pLayout);
    ui->tableWidget_Species->setCellWidget(0, 0, pWidget);
    connect(pCheckBox, SIGNAL(clicked(bool)), this, SLOT(tableSpeciesChanged()));
    pWidget = new QWidget(); pCheckBox = new QCheckBox(); pLayout = new QHBoxLayout(pWidget);
    pLayout->addWidget(pCheckBox); pLayout->setAlignment(Qt::AlignCenter); pLayout->setContentsMargins(0,0,0,0); pWidget->setLayout(pLayout);
    ui->tableWidget_Species->setCellWidget(1, 0, pWidget);
    connect(pCheckBox, SIGNAL(clicked(bool)), this, SLOT(tableSpeciesChanged()));
    ui->tableWidget_Species->setItemDelegate(new SpeciesTableDelegate());

    // q_max box
    dv = new QDoubleValidator(0.0, 1.0, 3);
    ui->lineEdit_qmax->setValidator(dv);
    ui->lineEdit_qmax->setEnabled(false);

    // Reactions panel.
    ui->tableWidget_Reactions->setColumnWidth(0, 25);
    ui->tableWidget_Reactions->setColumnWidth(1, 30);
    ui->tableWidget_Reactions->item(0, 1)->setFlags(ui->tableWidget_Reactions->item(0, 1)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
    pWidget = new QWidget(); pCheckBox = new QCheckBox(); pLayout = new QHBoxLayout(pWidget);
    pLayout->addWidget(pCheckBox); pLayout->setAlignment(Qt::AlignCenter); pLayout->setContentsMargins(0,0,0,0); pWidget->setLayout(pLayout);
    ui->tableWidget_Reactions->setCellWidget(0, 0, pWidget);
    connect(pCheckBox, SIGNAL(clicked(bool)), this, SLOT(tableReactionsChanged()));
    ui->tableWidget_Reactions->setItemDelegate(new ReactionsTableDelegate());

    // Butler-Volmer to start with...
    flags |= XMFLAG_BV;

    setupConnections();

    // Initialise stuff...
    ui->progressBar->setRange(0, 100);
    ui->progressBar->setValue(0);

    updateSimData();
    updatePotPlot();
    fittingChanged(ui->checkBox->isChecked());
    tableSpeciesChanged();
    tableReactionsChanged();
}

Xmarcus::~Xmarcus() {

    // Quit thread
    thread->quit();
    // Wait for it to be closed properly
    while(!thread->isFinished());
    // Delete thread and UI
    delete thread;

    delete ui;
}

void Xmarcus::errorString(QString error) {QMessageBox::information(this, tr("Error"), error);}
double Xmarcus::getResistance(int i) {return ui->tableWidget_AiRC->item(1, i+1)->text().toDouble();}
double Xmarcus::getDLCapacitance(int i) {return ui->tableWidget_AiRC->item(2, i+1)->text().toDouble();}
void Xmarcus::setQmax(void) {qmax = ui->lineEdit_qmax->text().toDouble()*10000.0;}
bool Xmarcus::isAreaFitChecked() {return ui->tableWidget_AiRC->cellWidget(0, 0)->findChild<QCheckBox *>()->isChecked();}
void Xmarcus::setUiArea(double _value) {ui->tableWidget_AiRC->item(0, 1)->setText(QString::number(_value*10000.0));}
bool Xmarcus::isResFitChecked() {return ui->tableWidget_AiRC->cellWidget(1, 0)->findChild<QCheckBox *>()->isChecked();}
void Xmarcus::setUiResistance(double _value) {ui->tableWidget_AiRC->item(1, 1)->setText(QString::number(_value));}
bool Xmarcus::isDLCapFitChecked() {return ui->tableWidget_AiRC->cellWidget(2, 0)->findChild<QCheckBox *>()->isChecked();}
void Xmarcus::setUiDLCapacitance(double _value) {ui->tableWidget_AiRC->item(2, 1)->setText(QString::number(_value));}
bool Xmarcus::isReactFitChecked(int i) {return ui->tableWidget_Reactions->cellWidget(i, 0)->findChild<QCheckBox *>()->isChecked();}
void Xmarcus::setUiReactParam(int i, int j, double _value) {ui->tableWidget_Reactions->item(i, j)->setText(QString::number(_value));}

/*
 * Set up connections for main window.
 */
void Xmarcus::setupConnections(void) {

    connect(simWindow, SIGNAL(autoFillClicked()), this, SLOT(autoFillFitData()));

    // Connect widgets on the main window.
    connect(ui->actionOpen_File, SIGNAL(triggered(bool)), this, SLOT(openFile()));
    connect(ui->actionSave, SIGNAL(triggered(bool)), this, SLOT(saveFile()));
    connect(ui->actionSave_As, SIGNAL(triggered(bool)), this, SLOT(saveAsFile()));
    connect(ui->actionStart_Simulation, SIGNAL(triggered(bool)), this, SLOT(startSimulation()));
    connect(ui->comboBox_model, SIGNAL(currentIndexChanged(int)), this, SLOT(modelChanged(int)));
    connect(ui->checkBox, SIGNAL(clicked(bool)), this, SLOT(fittingChanged(bool)));
    connect(ui->pushButton_Load, SIGNAL(clicked()), this, SLOT(loadFitFile()));
    connect(ui->spinBox, SIGNAL(valueChanged(int)), this, SLOT(sweepsChanged(int)));
    connect(ui->doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(scanRateChanged(double)));
    connect(ui->tableWidget, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tablePotChanged()));
    connect(ui->checkBox_iRC, SIGNAL(clicked(bool)), this, SLOT(iRCChanged(bool)));
    connect(ui->pushButton_NumParams, SIGNAL(clicked()), this, SLOT(checkNumParams()));
    connect(ui->tableWidget_Species, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableSpeciesChanged()));
    connect(ui->spinBox_species, SIGNAL(valueChanged(int)), this, SLOT(speciesChanged(int)));
    connect(ui->tableWidget_Reactions, SIGNAL(clicked(const QModelIndex &)), this, SLOT(tableReactionsClicked(const QModelIndex &)));
    connect(ui->tableWidget_Reactions, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableReactionsChanged()));
    connect(ui->spinBox_reactions, SIGNAL(valueChanged(int)), this, SLOT(reactionsChanged(int)));
    connect(ui->pushButton_checkReactions, SIGNAL(clicked()), this, SLOT(checkReactions()));

    // creating a new Worker instance and putting it on a QThread instance
    thread = new QThread();
    worker = new WorkerSimCV();
    QTimer *timer = new QTimer();
    timer->setInterval(0);  // Timer's inteveral set to 0 means that timer will trigger an event as soon as there are no other events to be processed

    connect(worker, SIGNAL(error(QString)), this, SLOT(errorString(QString)));
    connect(worker, SIGNAL(finished()), thread, SLOT(quit()));
    connect(worker, SIGNAL(sendProgressSim(int)), this, SLOT(receiveProgressSim(int)));
    connect(worker, SIGNAL(sendFinishedCV(ECData *)), this, SLOT(receiveFinishedCV(ECData *)), Qt::QueuedConnection);

    connect(this, SIGNAL(startSimThread()), worker, SLOT(startSim()));

    connect(timer, SIGNAL(timeout()), worker, SLOT(process()));
    connect(thread, SIGNAL(started()), timer, SLOT(start()));

    connect(thread, SIGNAL(finished()), worker, SLOT(deleteLater()));
    connect(thread, SIGNAL(finished()), timer, SLOT(deleteLater()));
    connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));

    // Start timer and move to thread
    timer->moveToThread(thread);

    // Move worker to thread
    worker->moveToThread(thread);

    worker->setXmParent(this);

    thread->start();

}

/*
 * Updates progress bar
 */
void Xmarcus::receiveProgressSim(int i) {
    ui->progressBar->setValue(i);

//    updatePotPlot(simData->getTime(i), simData->getPotential(i));

}

/*
 * Starts the simulation.
 */
void Xmarcus::startSimulation() {

    updateSimData();

    simData->calcCVTimePot();
    worker->setWorkAmount(simData->getPotential().size());

    checkReactions();
    checkNonLinear();

    emit startSimThread();
}

/*
 * Plot result when finished simulation.
 */
void Xmarcus::receiveFinishedCV(ECData *plotData) {

    ui->progressBar->setValue(100);

    simWindow->setSimData(plotData);
    simWindow->updateSimWindow();
    simWindow->show();

}

/*
 * Reads configuration file.
 */
void Xmarcus::openFile() {

    disconnect(ui->doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(scanRateChanged(double)));
    disconnect(ui->tableWidget, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tablePotChanged()));
    disconnect(ui->spinBox, SIGNAL(valueChanged(int)), this, SLOT(sweepsChanged(int)));

    int nos = 0; // number of species.
    int nor = 0; // number of reactions.
    QCheckBox *checkbox;

    QString inputFile = QFileDialog::getOpenFileName(this,
            tr("Open Configuration File"), "",
            tr("marcus (*.mrc);;marcus - old style (*.conf);;All Files (*)"));

    if (inputFile.isEmpty())
        return;
    else {
        QFile file(inputFile);
        if (!file.open(QIODevice::ReadOnly)) {
            QMessageBox::information(this, tr("Unable to open file"), file.errorString());
            return;
        }

        // Clear potentially full data structure.
        simData->clearData();

        // File info ...
        QFileInfo fi(inputFile);
        workingDir = fi.path();
        confFile = fi.fileName();

        if(!QDir::setCurrent(workingDir))
            qDebug() << "Could not change the current working directory";
        qDebug() << QDir::currentPath();

        // First go ....
        QTextStream in(&file);
        while (!in.atEnd()) {
            // Read line and check if empty or comment.
            QString line = in.readLine();
            QStringList fields = line.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
            if(fields.isEmpty()) continue;
            else if((fields[0][0] == '#') || fields[0].isEmpty()) continue;

            // Read-in pretty much copied from 'marcus'
            if(fields[0] == "type") {
                if(fields[1] == "CV") ui->comboBox->setCurrentIndex(ui->comboBox->findText("Cyclic Voltammetry"));
                else if(fields[1] == "TF") ui->comboBox->setCurrentIndex(ui->comboBox->findText("Tafel Plot"));
                continue;
            }

            if(fields[0] == "fit") {
                if(fields[1] == "ON") ui->checkBox->setChecked(true); //ui->pushButton->setEnabled(true); ui->comboBox_2->setEnabled(true); }
                else ui->checkBox->setChecked(false); // ui->pushButton->setEnabled(false); ui->comboBox_2->setEnabled(false); }
                continue;
            }
            if(fields[0] == "fitdata") {
                if(fields[1] == "ZAHNER") ui->comboBox_DataTypes->setCurrentIndex(ui->comboBox_DataTypes->findText("Zahner"));
                else if(fields[1] == "CSV") ui->comboBox_DataTypes->setCurrentIndex(ui->comboBox_DataTypes->findText("Comma Separated Values (CSV)"));
                else if(fields[1] == "ALB") ui->comboBox_DataTypes->setCurrentIndex(ui->comboBox_DataTypes->findText("Metrohm Autolab"));
                else if(fields[1] == "ECL") ui->comboBox_DataTypes->setCurrentIndex(ui->comboBox_DataTypes->findText("EClab, BioLogic"));
                if(fields.size() <= 2)
                    ui->label_FitFile->setText("no_file");
                else
                    ui->label_FitFile->setText(fields[2]);
                continue;
            }

            // Stuff for the Potential sub-frame...
            if(fields[0] == "sweeps") {simData->setNoSweeps(static_cast<int>(fields[1].toDouble()+0.5)); ui->spinBox->setValue(simData->getNoSweeps()); continue; }
            if(fields[0] == "rate") {simData->setScanRate(fields[1].toDouble()); ui->doubleSpinBox->setValue(simData->getScanRate()); continue; }
            if(fields[0] == "steps") {simData->setPotentialSteps(fields[1].toDouble()); ui->doubleSpinBox_2->setValue(simData->getPotentialSteps()); continue; }

            if(fields.at(0) == "temp") {
                simData->setTemperature(fields.at(1).toDouble());
                ui->lineEdit_temp->setText(fields.at(1));
                continue;
            }

            if(fields.at(0) == "area") {
                simData->appendArea(fields.at(1).toDouble()/10000.0); // convert into m^2
                ui->tableWidget_AiRC->item(0, 1)->setText(fields.at(1));
                continue;
            }
            if(fields.at(0) == "*area") {
                //simData->fitFlags |= ECDFLAG_AREA;
                checkbox = ui->tableWidget_AiRC->cellWidget(0, 0)->findChild<QCheckBox *>();
                checkbox->setEnabled(true); checkbox->setChecked(true);
                for(int i=0; i<3; i++) {
                    simData->appendArea(fields.at(i+1).toDouble()/10000.0); // convert into m^2
                    ui->tableWidget_AiRC->item(0, i+1)->setText(fields.at(i+1));
                }
                continue;
            }

            if(fields.at(0) == "ircdrop") {
                if(fields.at(1) == "ON") { ui->checkBox_iRC->setChecked(true); flags |= XMFLAG_IRC; }
                else { ui->checkBox_iRC->setChecked(false); flags &= ~XMFLAG_IRC; }
                continue;
            }

            if(fields.at(0) == "resist") {
                simData->appendResistance(fields.at(1).toDouble());
                ui->tableWidget_AiRC->item(1, 1)->setText(fields.at(1));
                continue;
            }
            if(fields.at(0) == "*resist") {
                //simData->fitFlags |= ECDFLAG_RES;
                checkbox = ui->tableWidget_AiRC->cellWidget(1, 0)->findChild<QCheckBox *>();
                checkbox->setEnabled(true); checkbox->setChecked(true);
                for(int i=0; i<3; i++) {
                    simData->appendResistance(fields.at(i+1).toDouble());
                    ui->tableWidget_AiRC->item(1, i+1)->setText(fields.at(i+1));
                }
                continue;
            }

            if(fields.at(0) == "dlcap") {
                simData->appendDLCapacitance(fields.at(1).toDouble());
                ui->tableWidget_AiRC->item(2, 1)->setText(fields.at(1));
                continue;
            }
            if(fields.at(0) == "*dlcap") {
                //simData->fitFlags |= ECDFLAG_DLC;
                checkbox = ui->tableWidget_AiRC->cellWidget(2, 0)->findChild<QCheckBox *>();
                checkbox->setEnabled(true); checkbox->setChecked(true);
                for(int i=0; i<3; i++) {
                    simData->appendDLCapacitance(fields.at(i+1).toDouble());
                    ui->tableWidget_AiRC->item(2, i+1)->setText(fields.at(i+1));
                }
                continue;
            }

            // Species parameters
            if(fields.at(0) == "species") {
                nos = fields.at(1).toInt();
                ui->spinBox_species->setValue(nos);
                continue;
            }

            if(fields.at(0) == "qmax") {
                qmax = fields.at(1).toDouble() * 10000.0; // Convert into m^2
                ui->lineEdit_qmax->setText(fields.at(1));
                continue;
            }

            // Electrochemical model
            if(fields.at(0) == "model") {
                flags &= ~XMFLAG_BV; flags &= ~XMFLAG_MH;
                if(fields.at(1) == "BV") {
                    flags |= XMFLAG_BV;
                    ui->comboBox_model->setCurrentIndex(0);
                } else if(fields.at(1) == "MH") {
                    flags |= XMFLAG_MH;
                    ui->comboBox_model->setCurrentIndex(1);
                } else
                    std::cout << "openFile(): unknown model.\n";
                continue;
            }

            // We might not need this! Just count the number on reactions in the first round...
/*            if(fields.at(0) == "react") {
                nor = fields.at(1).toInt();
                ui->spinBox_reactions->setValue(nor);
                continue;
            }
*/
            // Count number of reactions.
            if((fields.at(0) == "r") || (fields.at(0) == "*r"))
                nor++;

            // Numerical parameters
            if(fields.at(0) == "beta") {numParams.setBeta(fields.at(1).toDouble()); continue;}
            if(fields.at(0) == "Dm") {numParams.setDm(fields.at(1).toDouble()); continue;}
            if(fields.at(0) == "nlerror") {numParams.setNLerror(fields.at(1).toDouble()); continue;}
            if(fields.at(0) == "maxiter") {numParams.setMaxiter(fields.at(1).toInt()); continue;}
            if(fields.at(0) == "nlflag") {
                if(fields.at(1).toInt() == 0) flags &= ~XMFLAG_FRCNL;
                else flags |= XMFLAG_FRCNL;
                continue;
            }

        }

        // Second go ....
//        simData->getRPotentials().clear(); simData->getRTimes().clear();
        reactionsChanged(nor);
        int reactionCounter = 0;

//        qDebug() << "openFile(): nor = " << nor;

        in.seek(0);
        while (!in.atEnd()) {
            // Read line and check if empty or comment.
            QString line = in.readLine();
            QStringList fields = line.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
            if(fields.isEmpty()) continue;
            else if((fields.at(0).at(0) == '#') || fields.at(0).isEmpty()) continue;

            if(fields[0] == "pots") {
                double lastPot = fields[1].toDouble();
                simData->appendRPotentials(lastPot);
                simData->appendRTimes(0.0);
                double absPot = 0.0;

                for(int i=1; i<(simData->getNoSweeps()+1); i++) {
                    simData->appendRPotentials(fields[i+1].toDouble()); // first one is "pots"
                    absPot += fabs(fields[i+1].toDouble() - lastPot);
                    lastPot = fields[i+1].toDouble();
                    simData->appendRTimes(absPot/simData->getScanRate());
                }

                cpRPots2Table();

                continue;
            }

            // Read in species parameters
            QObject::disconnect(ui->tableWidget_Species, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableSpeciesChanged()));

            if(fields.at(0) == "symbols") {
                for(int i=0; i<nos; i++) {
                    allSpecies.at(i)->setSymbol(fields.at(i+1));
                    ui->tableWidget_Species->item(i, 1)->setText(fields.at(i+1));
                }
                continue;
            }

            if(fields.at(0) == "conc") {
                for(int i=0; i<nos; i++) {
                    allSpecies.at(i)->setCini(fields.at(i+1).toDouble() * 1000.0); // Convert concentrations into mol/m^3
                    ui->tableWidget_Species->item(i, 2)->setText(fields.at(i+1));
                }
                continue;
            }

            if(fields.at(0) == "diff") {
                for(int i=0; i<nos; i++) {
                    allSpecies.at(i)->setD(fields.at(i+1).toDouble() / 10000.0); // Convert into m^2/s
                    ui->tableWidget_Species->item(i, 3)->setText(fields.at(i+1));
                }
                continue;
            }

            if(fields.at(0) == "charge") {
                for(int i=0; i<nos; i++) {
                    allSpecies.at(i)->setCharge(fields.at(i+1).toInt());
                    ui->tableWidget_Species->item(i, 4)->setText(fields.at(i+1));
                }
                continue;
            }

            // Units!!!
            if(fields.at(0) == "adsco") {
                for(int i=0; i<nos; i++) {
                    allSpecies.at(i)->setKa(fields.at(2*i+1).toDouble());
                    allSpecies.at(i)->setKd(fields.at(2*i+2).toDouble());
                    ui->tableWidget_Species->item(i, 5)->setText(fields.at(2*i+1));
                    ui->tableWidget_Species->item(i, 6)->setText(fields.at(2*i+2));
                    if((allSpecies.at(i)->getKa() != 0.0) || (allSpecies.at(i)->getKd() != 0.0)) {
                        allSpecies.at(i)->setAdsorb(true);
                        checkbox = ui->tableWidget_Species->cellWidget(i, 0)->findChild<QCheckBox *>();
                        checkbox->setChecked(true);
                    } else
                        allSpecies.at(i)->setAdsorb(false);
                }
                continue;
            }

            QObject::connect(ui->tableWidget_Species, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableSpeciesChanged()));

            // Read in reaction parameters. MUST COME AFTER THE SPECIES STUFF!
            QObject::disconnect(ui->tableWidget_Reactions, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableReactionsChanged()));
            if((fields.at(0) == "r") || (fields.at(0) == "*r")) {
                QString strReaction;
                int counter = 0;
                bool disp = false;

                if(fields.at(0) == "*r") {
                    checkbox = ui->tableWidget_Reactions->cellWidget(reactionCounter, 0)->findChild<QCheckBox *>();
                    checkbox->setEnabled(true); checkbox->setChecked(true);
                }

                switch(fields.at(1).toInt()) {
                case 1: // E
                    ui->tableWidget_Reactions->item(reactionCounter, 1)->setText("E");
                    for(int i=0; i<nos; i++) if(fields.at(i+2) == "-1")
                        strReaction = ui->tableWidget_Species->item(i, 1)->text();
                    strReaction += " + e ⇄ ";
                    for(int i=0; i<nos; i++) if(fields.at(i+2) == "1")
                        strReaction += ui->tableWidget_Species->item(i, 1)->text();
                    for(int i=0; i<9; i++)
                        ui->tableWidget_Reactions->item(reactionCounter, i+3)->setText(fields.at(i+nos+2));
                    break;
                case 2: // C
                    ui->tableWidget_Reactions->item(reactionCounter, 1)->setText("C");
                    for(int i=0; i<nos; i++) if(fields.at(i+2) == "-1")
                        strReaction = ui->tableWidget_Species->item(i, 1)->text();
                    strReaction += " ⇄ ";
                    for(int i=0; i<nos; i++) if(fields.at(i+2) == "1")
                        strReaction += ui->tableWidget_Species->item(i, 1)->text();
                    ui->tableWidget_Reactions->item(reactionCounter, 3)->setText(fields.at(nos+2));
                    ui->tableWidget_Reactions->item(reactionCounter, 4)->setText(fields.at(nos+3));
                    ui->tableWidget_Reactions->item(reactionCounter, 5)->setText("-");
                    ui->tableWidget_Reactions->item(reactionCounter, 6)->setText(fields.at(nos+4));
                    ui->tableWidget_Reactions->item(reactionCounter, 7)->setText(fields.at(nos+5));
                    ui->tableWidget_Reactions->item(reactionCounter, 8)->setText("-");
                    ui->tableWidget_Reactions->item(reactionCounter, 9)->setText(fields.at(nos+6));
                    ui->tableWidget_Reactions->item(reactionCounter, 10)->setText(fields.at(nos+7));
                    ui->tableWidget_Reactions->item(reactionCounter, 11)->setText("-");
                    break;
                case 3: // higer-order C
                    counter = 0;
                    disp = false;
                    for(int i=0; i<nos; i++) {
                        if(fields.at(i+2) == "-1")
                            counter++;
                        else if(fields.at(i+2) == "-2")
                            counter += 2;
                    }
                    if(counter == 1) {
                        ui->tableWidget_Reactions->item(reactionCounter, 1)->setText("Cd");
                        for(int i=0; i<nos; i++) if(fields.at(i+2) == "-1")
                            strReaction = ui->tableWidget_Species->item(i, 1)->text() + " ⇄ ";
                        disp = true;
                    } else {
                        int intSpecies[2], j = 0; // placeholder for two species.
                        for(int i=0; i<nos; i++) {
                            if(fields.at(i+2) == "-1") {
                                intSpecies[j] = i;
                                j++;
                            } else if(fields.at(i+2) == "-2") {
                                intSpecies[j] = intSpecies[j+1] = i;
                                continue;
                            }
                        }
                        strReaction = ui->tableWidget_Species->item(intSpecies[0], 1)->text() + " + " +
                                ui->tableWidget_Species->item(intSpecies[1], 1)->text() + " ⇄ ";
                    }
                    counter = 0;
                    for(int i=0; i<nos; i++) {
                        if(fields.at(i+2) == "1")
                            counter++;
                        else if(fields.at(i+2) == "2")
                            counter += 2;
                    }
                    if(counter == 1) {
                        ui->tableWidget_Reactions->item(reactionCounter, 1)->setText("Cc");
                        for(int i=0; i<nos; i++) if(fields.at(i+2) == "1")
                            strReaction += ui->tableWidget_Species->item(i, 1)->text();
                    } else {
                        if(!disp)
                            ui->tableWidget_Reactions->item(reactionCounter, 1)->setText("C2");
                        int intSpecies[2], j = 0; // placeholder for two species.
                        for(int i=0; i<nos; i++) {
                            if(fields.at(i+2) == "1") {
                                intSpecies[j] = i;
                                j++;
                            } else if(fields.at(i+2) == "2") {
                                intSpecies[j] = intSpecies[j+1] = i;
                                continue;
                            }
                        }
                        strReaction += ui->tableWidget_Species->item(intSpecies[0], 1)->text() + " + " +
                                ui->tableWidget_Species->item(intSpecies[1], 1)->text();
                    }
                    ui->tableWidget_Reactions->item(reactionCounter, 3)->setText(fields.at(nos+2));
                    ui->tableWidget_Reactions->item(reactionCounter, 4)->setText(fields.at(nos+3));
                    ui->tableWidget_Reactions->item(reactionCounter, 5)->setText("-");
                    ui->tableWidget_Reactions->item(reactionCounter, 6)->setText(fields.at(nos+4));
                    ui->tableWidget_Reactions->item(reactionCounter, 7)->setText(fields.at(nos+5));
                    ui->tableWidget_Reactions->item(reactionCounter, 8)->setText("-");
                    ui->tableWidget_Reactions->item(reactionCounter, 9)->setText(fields.at(nos+6));
                    ui->tableWidget_Reactions->item(reactionCounter, 10)->setText(fields.at(nos+7));
                    ui->tableWidget_Reactions->item(reactionCounter, 11)->setText("-");
                    break;
                default:
                    std::cout << "openFile(): error in reading reactions.";
                }
                ui->tableWidget_Reactions->item(reactionCounter, 2)->setText(strReaction);
                reactionCounter++;
            }
            QObject::connect(ui->tableWidget_Reactions, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableReactionsChanged()));

            // Done, reaction parameters.

        }
        file.close();
    }

    // Load fitting data if possible.
    if(ui->label_FitFile->text().compare("no_file")) {
        if(!fitData) fitData = new ECData(this);
        else fitData->clearData();
        fitData->setDataFile(ui->label_FitFile->text());
        loadFitFile();
    }

    fittingChanged(ui->checkBox->isChecked());
    updatePotPlot();
    tableSpeciesChanged();
    tableReactionsChanged();

    connect(ui->doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(scanRateChanged(double)));
    connect(ui->tableWidget, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tablePotChanged()));
    connect(ui->spinBox, SIGNAL(valueChanged(int)), this, SLOT(sweepsChanged(int)));

}

/*
 * Opens 'Save as ...' window.
 */
void Xmarcus::saveAsFile() {
    confFile = QFileDialog::getSaveFileName(this, tr("Save Config File"), "", tr("marcus (*.mrc);;All Files (*)"));
    saveFile();
}

/*
 * Saves configuration file.
 */
void Xmarcus::saveFile() {
    if(confFile.isEmpty())
        confFile = QFileDialog::getSaveFileName(this, tr("Save Config File"), "", tr("marcus (*.mrc);;All Files (*)"));

    if (confFile.isEmpty())
        return;
    else {
        QFile file(confFile);
        if (!file.open(QIODevice::WriteOnly)) {
            QMessageBox::information(this, tr("Unable to open file"), file.errorString());
            return;
        }

        updateSimData();

        // Read everything from simData instead of widgets

        QTextStream out(&file);

        switch(ui->comboBox->currentIndex()) {
        case 0: out << "type\tCV" << Qt::endl; break;
        case 1: out << "type\tTF" << Qt::endl; break;
        }

        if(ui->checkBox->isChecked()) out << "fit\tON" << Qt::endl;
        else out << "fit\tOFF" << Qt::endl;

        switch(ui->comboBox_DataTypes->currentIndex()) { // Check that file exists!!!
        case 0: out << "fitdata\t" << "CSV\t" << ui->label_FitFile->text() << Qt::endl; break;
        case 1: out << "fitdata\t" << "ZAHNER\t" << ui->label_FitFile->text() << Qt::endl; break;
        case 2: out << "fitdata\t" << "ALB\t" << ui->label_FitFile->text() << Qt::endl; break;
        case 3: out << "fitdata\t" << "ECL\t" << ui->label_FitFile->text() << Qt::endl; break;
        }

        out << "sweeps\t" << ui->spinBox->value() << Qt::endl;
        out << "rate\t" << ui->doubleSpinBox->value() << Qt::endl;
        out << "steps\t" << ui->doubleSpinBox_2->value() << Qt::endl;
        out << "pots";
        int row = 0, col = 0;
        for(int i=0; i<=ui->spinBox->value(); i++) {
            out << '\t' << ui->tableWidget->item(row, col)->text();
            row = row ? 0 : 1; col = row ? col : col + 1;
        }
        out << Qt::endl;

        out << "temp\t" << simData->getTemperature() << Qt::endl;

        if((flags & XMFLAG_FIT) && isAreaFitChecked()) {
                //(simData->fitFlags & ECDFLAG_AREA)) {
            out << "*area";
            for(int i=0; i<simData->getArea().size(); i++)
                out << '\t' << simData->getArea(i) * 10000.0;
        } else
            out << "area\t" << simData->getArea(0) * 10000.0;
        out << Qt::endl;

        if(flags & XMFLAG_IRC) {
            out << "ircdrop\tON" << Qt::endl;
            if((flags & XMFLAG_FIT) && isResFitChecked()) {
                    //(simData->fitFlags & ECDFLAG_RES)) {
                out << "*resist";
                for(int i=0; i<simData->getResistance().size(); i++)
                    out << '\t' << simData->getResistance(i);
            } else
                out << "resist\t" << simData->getResistance(0);
            out << Qt::endl;
            if((flags & XMFLAG_FIT) && isDLCapFitChecked()) {
                    //(simData->fitFlags & ECDFLAG_DLC)) {
                out << "*dlcap";
                for(int i=0; i<simData->getDLCapacitance().size(); i++)
                    out << '\t' << simData->getDLCapacitance(i);
            } else
                out << "resist\t" << simData->getDLCapacitance(0);
            out << Qt::endl;
        } else
            out << "ircdrop\tOFF" << Qt::endl;

        // Species parameters.
        out << "species\t" << allSpecies.size() << Qt::endl;
        out << "symbols"; for(int i=0; i<allSpecies.size(); i++) out << "\t" << allSpecies.at(i)->getSymbol(); out << Qt::endl;
        out << "conc"; for(int i=0; i<allSpecies.size(); i++) out << "\t" << allSpecies.at(i)->getCini()/1000.0; out << Qt::endl;
        out << "diff"; for(int i=0; i<allSpecies.size(); i++) out << "\t" << allSpecies.at(i)->getD()*10000.0; out << Qt::endl;
        out << "charge"; for(int i=0; i<allSpecies.size(); i++) out << "\t" << allSpecies.at(i)->getCharge(); out << Qt::endl;
        out << "adsco";
        for(int i=0; i<allSpecies.size(); i++)
            out << "\t" << allSpecies.at(i)->getKa() << "\t" << allSpecies.at(i)->getKd();
        out << Qt::endl;
        out << "qmax\t" << ui->lineEdit_qmax->text() << Qt::endl;

        // Reaction parameter (and electrochemical model).
        if(ui->comboBox_model->currentIndex() == 0)
            out << "model\tBV" << Qt::endl;
        else if(ui->comboBox_model->currentIndex() == 1)
            out << "model\tMH" << Qt::endl;
        else
            out << "model\tunknown model" << Qt::endl;

        checkReactions();
        for(int i=0; i<allReactions.size(); i++) {
            if(allReactions.at(i)->getFit()) out << "*r";
            else out << "r";
            switch(allReactions.at(i)->getType()) {
            case RCTFLAG_E:
                out << "\t1";
                break;
            case RCTFLAG_C:
                out << "\t2";
                break;
            default:
                out << "\t3";
                break;
            }
            for(int j=0; j<allSpecies.size(); j++)
                out << "\t" << gsl_matrix_int_get(M, static_cast<size_t>(i), static_cast<size_t>(j));
            switch (allReactions.at(i)->getType()) {
            case RCTFLAG_E:
                out << "\t" << allReactions[i]->getRctEPara(0) << "\t" << allReactions[i]->getRctEPara(1) << "\t" << allReactions[i]->getRctEPara(2) <<
                       "\t" << allReactions[i]->getRctEPara(3) << "\t" << allReactions[i]->getRctEPara(4) << "\t" << allReactions[i]->getRctEPara(5) <<
                       "\t" << allReactions[i]->getRctEPara(6) << "\t" << allReactions[i]->getRctEPara(7) << "\t" << allReactions[i]->getRctEPara(8);
                break;
            default:
                out << "\t" << allReactions[i]->getRctCPara(0) << "\t" << allReactions[i]->getRctCPara(1) << // kf, kb
                       "\t" << allReactions[i]->getRctCPara(2) << "\t" << allReactions[i]->getRctCPara(3) << // Lkf, Lkb
                       "\t" << allReactions[i]->getRctCPara(4) << "\t" << allReactions[i]->getRctCPara(5); // Ukf, Ukb
                break;
            }
            out << Qt::endl;
        }

        // Numerical parameters.
        out << "beta\t" << numParams.getBeta() << Qt::endl;
        out << "Dm\t" << numParams.getDm() << Qt::endl;
        out << "nlerror\t" << numParams.getNLerror() << Qt::endl;
        out << "maxiter\t" << numParams.getMaxiter() << Qt::endl;
        if(flags & XMFLAG_FRCNL) out << "nlflag\t1" << Qt::endl;
        else out << "nlflag\t0" << Qt::endl;

        file.close();
    }

}

/*
 * Rebuilds the potentials table using rPotentials from simData.
 */
void Xmarcus::cpRPots2Table() {

    // Add additional columns.
    int n = static_cast<int>(ceil((simData->getNoSweeps()+1.0)/2.0));
    ui->tableWidget->setColumnCount(n);
    for(int i=0; i<n; i++) ui->tableWidget->setColumnWidth(i, 80);

    // Fill in the return potentials
    int row = 0;
    int col = 0;
    Qt::ItemFlags flags;
    for(int i=0; i<simData->getRPotentials().size(); i++) {
        if(!ui->tableWidget->item(row, col)) {
            ui->tableWidget->setItem(row, col, new QTableWidgetItem(QString::number(simData->getRPotentials(i))));
        } else if(ui->tableWidget->item(row, col)->text() == "-") {
            flags = ui->tableWidget->item(row, col)->flags();
            flags |= Qt::ItemIsSelectable | Qt::ItemIsEditable;
            ui->tableWidget->item(row, col)->setFlags(flags);
        } else
            ui->tableWidget->item(row, col)->setText(QString::number(simData->getRPotentials(i)));
        row = row ? 0 : 1; col = row ? col : col + 1;
    }

    // Disable unused items
    if(simData->getNoSweeps() % 2 == 0) { // even
        if(!ui->tableWidget->item(1, n-1))
            ui->tableWidget->setItem(1, n-1, new QTableWidgetItem("-"));
        else
            ui->tableWidget->item(1, n-1)->setText("-");
        ui->tableWidget->item(1, n-1)->setFlags(ui->tableWidget->item(1, n-1)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
        ui->tableWidget->item(1, n-1)->setBackground(QColor(200, 200, 200, 255));
        ui->tableWidget->item(1, n-1)->setForeground(QColor(100, 100, 100, 255));
    }
}

/*
 * BV or MH have changed.
 */
void Xmarcus::modelChanged(int model) {

    flags &= ~XMFLAG_BV; flags &= ~XMFLAG_MH;

    if(model == 0)
        flags |= XMFLAG_BV;
    else
        flags |= XMFLAG_MH;
}

/*
 * Load data file for fitting.
 */
void Xmarcus::loadFitFile() {
    QString fileName;

    if(!fitData) fitData = new ECData(this);
    else {
        fileName = fitData->getDataFile();
        fitData->clearData();
        fitData->setDataFile(fileName);
    }

    if(fitData->getDataFile().isEmpty())
        fileName = QFileDialog::getOpenFileName(this, tr("Open CV File"), "", tr("Comma Separated Values (*.csv);;Text File (*.txt);;All Files (*)"));
//    else
//        fileName = fitData->getDataFile();

    QFileInfo fileInfo(fileName);
    ui->label_FitFile->setText(fileInfo.fileName());

    if (fileName.isEmpty())
        return;
    else {

        fitData->readFile(fileName, ui->comboBox_DataTypes->currentIndex());

        simWindow->setFitData(fitData);
        simWindow->updateSimWindow();
        simWindow->show();
    }

//    updateSimData();
    updatePotPlot();
}

void Xmarcus::autoFillFitData() {

    disconnect(ui->doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(scanRateChanged(double)));
    disconnect(ui->tableWidget, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tablePotChanged()));

    ui->spinBox->setValue(fitData->getNoSweeps());
    ui->doubleSpinBox->setValue(fitData->getScanRate());
    ui->doubleSpinBox_2->setValue(fitData->getPotentialSteps());

    simData->setRPotentials(fitData->getRPotentials());
    simData->setRTimes(fitData->getRTimes());

    cpRPots2Table();

    connect(ui->doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(scanRateChanged(double)));
    connect(ui->tableWidget, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tablePotChanged()));

    updateSimData();
    updatePotPlot();
}

void Xmarcus::sweepsChanged(int sweeps) {
    int n = static_cast<int>(ceil((sweeps+1.0)/2.0));
    ui->tableWidget->setColumnCount(n);
    for(int i=0; i<n; i++) ui->tableWidget->setColumnWidth(i, 80);

    // Fill new cells with numbers?
    for(int i=0; i<ui->tableWidget->columnCount(); i++) {
        for(int j=0; j<ui->tableWidget->rowCount(); j++) {
            if(!ui->tableWidget->item(j, i)) {
                ui->tableWidget->setItem(j, i, new QTableWidgetItem(ui->tableWidget->item(j, i-1)->text()));
            } else if(ui->tableWidget->item(j, i)->text() == "-") {
                ui->tableWidget->item(j, i)->setFlags(ui->tableWidget->item(j, i)->flags() | Qt::ItemIsSelectable | Qt::ItemIsEditable);
                ui->tableWidget->item(j, i)->setText(ui->tableWidget->item(j, i-1)->text());
                ui->tableWidget->item(j, i)->setBackground(QColor(255, 255, 255, 255));
                ui->tableWidget->item(j, i)->setForeground(QColor(0, 0, 0, 255));
            }
        }
    }

    // Disable unused items
    if(sweeps % 2 == 0) { // even
        ui->tableWidget->item(1, n-1)->setText("-");
        ui->tableWidget->item(1, n-1)->setFlags(ui->tableWidget->item(1, n-1)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
        ui->tableWidget->item(1, n-1)->setBackground(QColor(200, 200, 200, 255));
        ui->tableWidget->item(1, n-1)->setForeground(QColor(100, 100, 100, 255));
    }

    updateSimData();
    updatePotPlot();
}

void Xmarcus::scanRateChanged(double) {
    updateSimData();
    updatePotPlot();
}

void Xmarcus::tablePotChanged(void) {
    updateSimData();
    updatePotPlot();
}

/*
 * Updates the simData object with the data from the UI.
 */
void Xmarcus::updateSimData(void) {

    simData->clearData();

    simData->setNoSweeps(ui->spinBox->value());
    simData->setScanRate(ui->doubleSpinBox->value());
    simData->setPotentialSteps(ui->doubleSpinBox_2->value());

    // Fill first one "by hand" ...
    QTableWidgetItem *item = ui->tableWidget->item(0, 0);
    double lastPot = item->text().toDouble();
    simData->appendRPotentials(lastPot);
    simData->appendRTimes(0.0);
    // ... and now the rest.
    int row = 1; int col = 0;
    double absPot = 0.0;
    for(int i = 1; i<=simData->getNoSweeps(); i++) {
        item = ui->tableWidget->item(row, col);
        simData->appendRPotentials(item->text().toDouble());
        absPot += fabs(item->text().toDouble() - lastPot);
        lastPot = item->text().toDouble();
        simData->appendRTimes(absPot/simData->getScanRate());
        row = row ? 0 : 1;
        col = row ? col : col + 1;
    }

    simData->setTemperature(ui->lineEdit_temp->text().toDouble());

    // Electrode area
    simData->appendArea(ui->tableWidget_AiRC->item(0, 1)->text().toDouble()/10000.0);
    if((flags & XMFLAG_FIT) && isAreaFitChecked()) {
        simData->appendArea(ui->tableWidget_AiRC->item(0, 2)->text().toDouble()/10000.0);
        simData->appendArea(ui->tableWidget_AiRC->item(0, 3)->text().toDouble()/10000.0);
    }
    if(flags & XMFLAG_IRC) {
        simData->appendResistance(ui->tableWidget_AiRC->item(1, 1)->text().toDouble());
        simData->appendDLCapacitance(ui->tableWidget_AiRC->item(2, 1)->text().toDouble());
        if(flags & XMFLAG_FIT) {
            if(isResFitChecked()) {
                simData->appendResistance(ui->tableWidget_AiRC->item(1, 2)->text().toDouble());
                simData->appendResistance(ui->tableWidget_AiRC->item(1, 3)->text().toDouble());
            }
            if(isDLCapFitChecked()) {
                simData->appendDLCapacitance(ui->tableWidget_AiRC->item(2, 2)->text().toDouble());
                simData->appendDLCapacitance(ui->tableWidget_AiRC->item(2, 3)->text().toDouble());

            }
        }
    }
}

/*
 * Update the potential vs time plot.
 */
void Xmarcus::updatePotPlot(double t, double E) {
    box dataRange;

    dataRange.xmin = dataRange.xmax = simData->getRTimes(0);
    dataRange.ymin = dataRange.ymax = simData->getRPotentials(0);
    for(int i=1; i<simData->getRPotentials().size(); i++) {
        dataRange.xmin = dataRange.xmin < simData->getRTimes(i) ? dataRange.xmin : simData->getRTimes(i);
        dataRange.xmax = dataRange.xmax > simData->getRTimes(i) ? dataRange.xmax : simData->getRTimes(i);
        dataRange.ymin = dataRange.ymin < simData->getRPotentials(i) ? dataRange.ymin : simData->getRPotentials(i);
        dataRange.ymax = dataRange.ymax > simData->getRPotentials(i) ? dataRange.ymax : simData->getRPotentials(i);
    }
    double offset = fabs(dataRange.xmin * 0.05) + fabs(dataRange.xmax * 0.05);
    dataRange.xmin -= offset;
    dataRange.xmax += offset;
    offset = fabs(dataRange.ymin * 0.05) + fabs(dataRange.ymax * 0.05);
    dataRange.ymin -= offset;
    dataRange.ymax += offset;

    ui->customPlot->addGraph();
    ui->customPlot->graph(0)->setLineStyle(QCPGraph::lsLine);
    ui->customPlot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->customPlot->graph(0)->setData(simData->getRTimes(), simData->getRPotentials());

    if(t >= 0.0) {
        QVector<double> thelp; thelp.append(t);
        QVector<double> Ehelp; Ehelp.append(E);
        ui->customPlot->addGraph();
        ui->customPlot->graph(1)->setPen(QPen(Qt::red));
        ui->customPlot->graph(1)->setLineStyle(QCPGraph::lsNone);
        ui->customPlot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 8));
        ui->customPlot->graph(1)->setData(thelp, Ehelp);
    }

    ui->customPlot->xAxis->setRange(dataRange.xmin, dataRange.xmax);
    ui->customPlot->xAxis->setLabel("time / [s]");
    ui->customPlot->yAxis->setRange(dataRange.ymin, dataRange.ymax);
    ui->customPlot->yAxis->setLabel("E / [V]");

    ui->customPlot->replot();
}

/*
 * Fitting checkbox has been clicked.
 */
void Xmarcus::fittingChanged(bool checked) {
    QCheckBox *checkbox;

    if(checked) {
        flags |= XMFLAG_FIT;
//        ui->pushButton->setEnabled(true);
//        ui->comboBox_2->setEnabled(true);
        checkbox = ui->tableWidget_AiRC->cellWidget(0, 0)->findChild<QCheckBox *>();
        checkbox->setEnabled(true); checkbox->setChecked(true);
        ui->tableWidget_AiRC->item(0, 0)->setBackground(QColor(255, 255, 255, 255));
        ui->tableWidget_AiRC->item(0, 2)->setFlags(ui->tableWidget_AiRC->item(0, 2)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
        ui->tableWidget_AiRC->item(0, 2)->setBackground(QColor(255, 255, 255, 255));
        ui->tableWidget_AiRC->item(0, 2)->setForeground(QColor(0, 0, 0, 255));
        ui->tableWidget_AiRC->item(0, 3)->setFlags(ui->tableWidget_AiRC->item(0, 3)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
        ui->tableWidget_AiRC->item(0, 3)->setBackground(QColor(255, 255, 255, 255));
        ui->tableWidget_AiRC->item(0, 3)->setForeground(QColor(0, 0, 0, 255));
        if(flags & XMFLAG_IRC) {
            for(int i=1; i<3; i++) { // Enable fitting in AiRC table
                checkbox = ui->tableWidget_AiRC->cellWidget(i, 0)->findChild<QCheckBox *>();
                checkbox->setEnabled(true); checkbox->setChecked(true);
                ui->tableWidget_AiRC->item(i, 0)->setBackground(QColor(255, 255, 255, 255));
                ui->tableWidget_AiRC->item(i, 2)->setFlags(ui->tableWidget_AiRC->item(i, 2)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
                ui->tableWidget_AiRC->item(i, 2)->setBackground(QColor(255, 255, 255, 255));
                ui->tableWidget_AiRC->item(i, 2)->setForeground(QColor(0, 0, 0, 255));
                ui->tableWidget_AiRC->item(i, 3)->setFlags(ui->tableWidget_AiRC->item(i, 3)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
                ui->tableWidget_AiRC->item(i, 3)->setBackground(QColor(255, 255, 255, 255));
                ui->tableWidget_AiRC->item(i, 3)->setForeground(QColor(0, 0, 0, 255));
            }
        }
    } else {
        flags &= ~XMFLAG_FIT;
//        ui->pushButton->setEnabled(false);
//        ui->comboBox_2->setEnabled(false);
        for(int i=0; i<3; i++) { // Disable fitting in AiRC table.
            checkbox = ui->tableWidget_AiRC->cellWidget(i, 0)->findChild<QCheckBox *>();
            checkbox->setEnabled(false); checkbox->setChecked(false);
            ui->tableWidget_AiRC->item(i, 0)->setBackground(QColor(200, 200, 200, 255));
//            ui->tableWidget_AiRC->item(i, 0)->setForeground(QColor(0, 0, 0, 255));
            ui->tableWidget_AiRC->item(i, 2)->setFlags(ui->tableWidget_AiRC->item(i, 2)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_AiRC->item(i, 2)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_AiRC->item(i, 2)->setForeground(QColor(100, 100, 100, 255));
            ui->tableWidget_AiRC->item(i, 3)->setFlags(ui->tableWidget_AiRC->item(i, 3)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_AiRC->item(i, 3)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_AiRC->item(i, 3)->setForeground(QColor(100, 100, 100, 255));
        }
    }

    if(flags & XMFLAG_AIRCUPD) {
        flags &= ~XMFLAG_AIRCUPD;
        ui->tableWidget_AiRC->viewport()->update();
    } else {
        flags |= XMFLAG_AIRCUPD;
        iRCChanged(ui->checkBox_iRC->isChecked());
    }

    tableReactionsChanged();

}

/*
 * iRC checkbox has been changed.
 */
void Xmarcus::iRCChanged(bool checked) {
    QCheckBox *checkbox;

    if(checked) {
        flags |= XMFLAG_IRC;
         if(flags & XMFLAG_FIT) {
             checkbox = ui->tableWidget_AiRC->cellWidget(1, 0)->findChild<QCheckBox *>();
             checkbox->setEnabled(true); checkbox->setChecked(true);
             checkbox = ui->tableWidget_AiRC->cellWidget(2, 0)->findChild<QCheckBox *>();
             checkbox->setEnabled(true); checkbox->setChecked(true);
             ui->tableWidget_AiRC->item(1, 0)->setBackground(QColor(255, 255, 255, 255));
             ui->tableWidget_AiRC->item(1, 0)->setForeground(QColor(0, 0, 0, 255));
             ui->tableWidget_AiRC->item(2, 0)->setBackground(QColor(255, 255, 255, 255));
             ui->tableWidget_AiRC->item(2, 0)->setForeground(QColor(0, 0, 0, 255));
             for(int i=1; i<4; i++) {
                 ui->tableWidget_AiRC->item(1, i)->setFlags(ui->tableWidget_AiRC->item(1, i)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
                 ui->tableWidget_AiRC->item(1, i)->setBackground(QColor(255, 255, 255, 255));
                 ui->tableWidget_AiRC->item(1, i)->setForeground(QColor(0, 0, 0, 255));
                 ui->tableWidget_AiRC->item(2, i)->setFlags(ui->tableWidget_AiRC->item(2, i)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
                 ui->tableWidget_AiRC->item(2, i)->setBackground(QColor(255, 255, 255, 255));
                 ui->tableWidget_AiRC->item(2, i)->setForeground(QColor(0, 0, 0, 255));
             }
         } else {
             ui->tableWidget_AiRC->item(1, 1)->setFlags(ui->tableWidget_AiRC->item(1, 1)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
             ui->tableWidget_AiRC->item(1, 1)->setBackground(QColor(255, 255, 255, 255));
             ui->tableWidget_AiRC->item(1, 1)->setForeground(QColor(0, 0, 0, 255));
             ui->tableWidget_AiRC->item(2, 1)->setFlags(ui->tableWidget_AiRC->item(2, 1)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
             ui->tableWidget_AiRC->item(2, 1)->setBackground(QColor(255, 255, 255, 255));
             ui->tableWidget_AiRC->item(2, 1)->setForeground(QColor(0, 0, 0, 255));
         }
    } else {
        flags &= ~XMFLAG_IRC;
        checkbox = ui->tableWidget_AiRC->cellWidget(1, 0)->findChild<QCheckBox *>();
        checkbox->setEnabled(false); checkbox->setChecked(false);
        checkbox = ui->tableWidget_AiRC->cellWidget(2, 0)->findChild<QCheckBox *>();
        checkbox->setEnabled(false); checkbox->setChecked(false);
        ui->tableWidget_AiRC->item(1, 0)->setBackground(QColor(200, 200, 200, 255));
        ui->tableWidget_AiRC->item(1, 0)->setForeground(QColor(100, 100, 100, 255));
        ui->tableWidget_AiRC->item(2, 0)->setBackground(QColor(200, 200, 200, 255));
        ui->tableWidget_AiRC->item(2, 0)->setForeground(QColor(100, 100, 100, 255));
        for(int i=1; i<4; i++) {
            ui->tableWidget_AiRC->item(1, i)->setFlags(ui->tableWidget_AiRC->item(1, i)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_AiRC->item(1, i)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_AiRC->item(1, i)->setForeground(QColor(100, 100, 100, 255));
            ui->tableWidget_AiRC->item(2, i)->setFlags(ui->tableWidget_AiRC->item(2, i)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_AiRC->item(2, i)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_AiRC->item(2, i)->setForeground(QColor(100, 100, 100, 255));
        }
    }

    if(flags & XMFLAG_AIRCUPD) {
        flags &= ~XMFLAG_AIRCUPD;
        ui->tableWidget_AiRC->viewport()->update();
    } else {
        flags |= XMFLAG_AIRCUPD;
        fittingChanged(ui->checkBox->isChecked());
    }

}

/*
 * Opens the numerical parameters window.
 */
void Xmarcus::checkNumParams(void) {
/*
    numParams.setBeta(thePlaza->getBeta());
    numParams.setDm(thePlaza->getDm());
    numParams.setNlerror(thePlaza->getNlerror());
    numParams.setMaxiter(thePlaza->getMaxiter());
    if(flags & XMFLAG_FRCNL)
        numParams.setNLbox(true);
    else
        numParams.setNLbox(false);
*/
    if(numParams.exec() == QDialog::Accepted) {

        // Need to update the variables of NumParams.

//        thePlaza->setBeta(numParams.getBeta());
//        thePlaza->setDm(numParams.getDM());
//        thePlaza->setNLerror(numParams.getNLerror());
//        thePlaza->setMaxiter(numParams.getMaxiter());
        if(numParams.getNLbox())
            flags |= XMFLAG_FRCNL;
        else
            flags &= ~XMFLAG_FRCNL;
    }

}

/*
 * Runs if there is a change in the Species table. Update Species.
 */
void Xmarcus::tableSpeciesChanged(void) {
    QCheckBox *checkbox;

    disconnect(ui->tableWidget_Species, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableSpeciesChanged()));

    allSpecies.clear();
    flags &= ~XMFLAG_ADS;

    for(int i=0; i<ui->tableWidget_Species->rowCount(); i++) {
        allSpecies.append(new Species());
        checkbox = ui->tableWidget_Species->cellWidget(i, 0)->findChild<QCheckBox *>();
        allSpecies.last()->setAdsorb(checkbox->isChecked());
        allSpecies.last()->setSymbol(ui->tableWidget_Species->item(i, 1)->text());
        allSpecies.last()->setCini(ui->tableWidget_Species->item(i, 2)->text().toDouble()*1000.0); //Convert concentrations into mol/m^3
        allSpecies.last()->setD(ui->tableWidget_Species->item(i, 3)->text().toDouble()/10000.0); // Convert into m^2/s
        allSpecies.last()->setCharge(ui->tableWidget_Species->item(i, 4)->text().toInt());
        if(checkbox->isChecked()) { // ToDo: make sure ka and kd are not zero!
            flags |= XMFLAG_ADS;
            allSpecies.last()->setKa(ui->tableWidget_Species->item(i, 5)->text().toDouble()); // Unit!!!!!!
            allSpecies.last()->setKd(ui->tableWidget_Species->item(i, 6)->text().toDouble());

            ui->tableWidget_Species->item(i, 5)->setFlags(ui->tableWidget_Species->item(i, 5)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Species->item(i, 5)->setBackground(QColor(255, 255, 255, 255));
            ui->tableWidget_Species->item(i, 5)->setForeground(QColor(0, 0, 0, 255));
            ui->tableWidget_Species->item(i, 6)->setFlags(ui->tableWidget_Species->item(i, 6)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Species->item(i, 6)->setBackground(QColor(255, 255, 255, 255));
            ui->tableWidget_Species->item(i, 6)->setForeground(QColor(0, 0, 0, 255));
        } else {
            allSpecies.last()->setKa(0.0); ui->tableWidget_Species->item(i, 5)->setText("1.0");
            allSpecies.last()->setKd(0.0); ui->tableWidget_Species->item(i, 6)->setText("1.0");

            ui->tableWidget_Species->item(i, 5)->setFlags(ui->tableWidget_Species->item(i, 5)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Species->item(i, 5)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_Species->item(i, 5)->setForeground(QColor(100, 100, 100, 255));
            ui->tableWidget_Species->item(i, 6)->setFlags(ui->tableWidget_Species->item(i, 6)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Species->item(i, 6)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_Species->item(i, 6)->setForeground(QColor(100, 100, 100, 255));
        }
    }

    // Enable/disable q_max input.
    ui->lineEdit_qmax->setEnabled(false);
    for(int i=0; i<allSpecies.size(); i++)
        if(allSpecies.at(i)->getAdsorb())
            ui->lineEdit_qmax->setEnabled(true);

    connect(ui->tableWidget_Species, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableSpeciesChanged()));

//    qDebug() << "was here..." << allSpecies.last()->getSymbol() << "..." << allSpecies.last()->getKd();
}

/*
 * Number of species has changed.
 */
void Xmarcus::speciesChanged(int species) {

    if(species == allSpecies.size())
        return;
    else if(species > allSpecies.size()) {
        QObject::disconnect(ui->tableWidget_Species, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableSpeciesChanged()));
        for(int i=ui->tableWidget_Species->rowCount(); i<species; i++) {
            allSpecies.append(new Species(QString(QChar(static_cast<char>(65+i)))));

            ui->tableWidget_Species->setRowCount(ui->tableWidget_Species->rowCount()+1);

            QWidget *pWidget = new QWidget();
            QCheckBox *pCheckBox = new QCheckBox();
            QHBoxLayout *pLayout = new QHBoxLayout(pWidget);
            pLayout->addWidget(pCheckBox); pLayout->setAlignment(Qt::AlignCenter); pLayout->setContentsMargins(0,0,0,0); pWidget->setLayout(pLayout);
            ui->tableWidget_Species->setCellWidget(i, 0, pWidget);
            QObject::connect(pCheckBox, SIGNAL(clicked(bool)), this, SLOT(tableSpeciesChanged()));
            ui->tableWidget_Species->setItem(i, 1, new QTableWidgetItem(allSpecies.last()->getSymbol()));
            ui->tableWidget_Species->setItem(i, 2, new QTableWidgetItem(QString::number(allSpecies.last()->getCini()/1000.0)));
            ui->tableWidget_Species->setItem(i, 3, new QTableWidgetItem(QString::number(allSpecies.last()->getD()*10000.0)));
            ui->tableWidget_Species->setItem(i, 4, new QTableWidgetItem(QString::number(allSpecies.last()->getCharge())));
            ui->tableWidget_Species->setItem(i, 5, new QTableWidgetItem(QString::number(allSpecies.last()->getKa())));
            ui->tableWidget_Species->setItem(i, 6, new QTableWidgetItem(QString::number(allSpecies.last()->getKd())));
        }
        QObject::connect(ui->tableWidget_Species, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableSpeciesChanged()));
    } else if(species < allSpecies.size()) {
        QObject::disconnect(ui->tableWidget_Species, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableSpeciesChanged()));
        for(int i=ui->tableWidget_Species->rowCount()-1; i>species-1; i--) {
            allSpecies.remove(i);
            ui->tableWidget_Species->removeRow(i);
        }
        QObject::connect(ui->tableWidget_Species, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableSpeciesChanged()));
    }

    tableSpeciesChanged();
}

/*
 * Reaction table clicked. Open choose reaction window.
 */
void Xmarcus::tableReactionsClicked(const QModelIndex &index) {

    int row = index.row();
    int column = index.column();

    if(column == 2) {
        QStringList speciesList;
        for(int i=0; i<allSpecies.size(); i++)
            speciesList << allSpecies.at(i)->getSymbol();
        chooseReaction.setSpeciesList(speciesList);

        if(chooseReaction.exec() == QDialog::Accepted) {
            allReactions.at(row)->setType(chooseReaction.getType());
            switch(chooseReaction.getType()) {
            case RCTFLAG_E:
                ui->tableWidget_Reactions->item(row, 1)->setText("E");
                break;
            case RCTFLAG_C:
                ui->tableWidget_Reactions->item(row, 1)->setText("C");
                break;
            case RCTFLAG_C2:
                ui->tableWidget_Reactions->item(row, 1)->setText("C2");
                break;
            case RCTFLAG_CD:
                ui->tableWidget_Reactions->item(row, 1)->setText("Cd");
                break;
            case RCTFLAG_CC:
                ui->tableWidget_Reactions->item(row, 1)->setText("Cc");
                break;
            default:
                qDebug() << "\'tableReactionsClicked()\' error here!";
            }
            if(chooseReaction.getType() != RCTFLAG_E) {
                ui->tableWidget_Reactions->item(row, 5)->setText("-");
                ui->tableWidget_Reactions->item(row, 8)->setText("-");
                ui->tableWidget_Reactions->item(row, 11)->setText("-");
                ui->tableWidget_Reactions->item(row, 5)->setFlags(ui->tableWidget_Reactions->item(row, 5)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
                ui->tableWidget_Reactions->item(row, 5)->setBackground(QColor(200, 200, 200, 255));
                ui->tableWidget_Reactions->item(row, 5)->setForeground(QColor(100, 100, 100, 255));
                ui->tableWidget_Reactions->item(row, 8)->setFlags(ui->tableWidget_Reactions->item(row, 8)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
                ui->tableWidget_Reactions->item(row, 8)->setBackground(QColor(200, 200, 200, 255));
                ui->tableWidget_Reactions->item(row, 8)->setForeground(QColor(100, 100, 100, 255));
                ui->tableWidget_Reactions->item(row, 11)->setFlags(ui->tableWidget_Reactions->item(row, 11)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
                ui->tableWidget_Reactions->item(row, 11)->setBackground(QColor(200, 200, 200, 255));
                ui->tableWidget_Reactions->item(row, 11)->setForeground(QColor(100, 100, 100, 255));
            } else {
                ui->tableWidget_Reactions->item(row, 5)->setText("1.0e3");
                ui->tableWidget_Reactions->item(row, 8)->setText("1.0e3");
                ui->tableWidget_Reactions->item(row, 11)->setText("1.0e3");
            }
            allReactions.at(row)->setStrReaction(chooseReaction.getStrReaction());
            ui->tableWidget_Reactions->item(row, column)->setText(chooseReaction.getStrReaction());
        }
    }
//    qDebug() << "row = " << row << "\tcol = " << column;
}

/*
 * Reaction table changed. Rebuild reactions.
 */
void Xmarcus::tableReactionsChanged(void) {
    QCheckBox *checkbox;

    disconnect(ui->tableWidget_Reactions, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableReactionsChanged()));

    allReactions.clear();

    for(int i=0; i<ui->tableWidget_Reactions->rowCount(); i++) {
        allReactions.append(new Reaction());
        checkbox = ui->tableWidget_Reactions->cellWidget(i, 0)->findChild<QCheckBox *>();

        disconnect(checkbox, SIGNAL(clicked(bool)), this, SLOT(tableReactionsChanged()));

        if(flags & XMFLAG_FIT) { // Checkbox enable/disable
            checkbox->setEnabled(true); //checkbox->setChecked(true);
            ui->tableWidget_Reactions->item(i, 0)->setBackground(QColor(255, 255, 255, 255));
        } else {
            checkbox->setEnabled(false); checkbox->setChecked(false);
            ui->tableWidget_Reactions->item(i, 0)->setBackground(QColor(200, 200, 200, 255));
        }

        allReactions.last()->setFit(checkbox->isChecked());
        allReactions.last()->setStrReaction(ui->tableWidget_Reactions->item(i, 2)->text());
        if(ui->tableWidget_Reactions->item(i, 1)->text() == "E") {
            allReactions.last()->setType(RCTFLAG_E);
            allReactions.last()->setRctEPara(0, ui->tableWidget_Reactions->item(i, 3)->text().toDouble()); // E0
            allReactions.last()->setRctEPara(1, ui->tableWidget_Reactions->item(i, 4)->text().toDouble()); // alpha, lambda
            allReactions.last()->setRctEPara(2, ui->tableWidget_Reactions->item(i, 5)->text().toDouble()); // k0
            allReactions.last()->setRctEPara(3, ui->tableWidget_Reactions->item(i, 6)->text().toDouble()); // lE0
            allReactions.last()->setRctEPara(4, ui->tableWidget_Reactions->item(i, 7)->text().toDouble()); // lalpha, llambda
            allReactions.last()->setRctEPara(5, ui->tableWidget_Reactions->item(i, 8)->text().toDouble()); // lk0
            allReactions.last()->setRctEPara(6, ui->tableWidget_Reactions->item(i, 9)->text().toDouble()); // uE0
            allReactions.last()->setRctEPara(7, ui->tableWidget_Reactions->item(i, 10)->text().toDouble()); // ualpha, ulambda
            allReactions.last()->setRctEPara(8, ui->tableWidget_Reactions->item(i, 11)->text().toDouble()); // uk0
        } else { // ... any other reaction.
            if(ui->tableWidget_Reactions->item(i, 1)->text() == "C")
                allReactions.last()->setType(RCTFLAG_C);
            else if(ui->tableWidget_Reactions->item(i, 1)->text() == "C2")
                allReactions.last()->setType(RCTFLAG_C2);
            else if(ui->tableWidget_Reactions->item(i, 1)->text() == "Cd")
                allReactions.last()->setType(RCTFLAG_CD);
            else if(ui->tableWidget_Reactions->item(i, 1)->text() == "Cc")
                allReactions.last()->setType(RCTFLAG_CC);

            allReactions.last()->setRctCPara(0, ui->tableWidget_Reactions->item(i, 3)->text().toDouble()); // kf
            allReactions.last()->setRctCPara(1, ui->tableWidget_Reactions->item(i, 4)->text().toDouble()); // kb
            allReactions.last()->setRctCPara(2, ui->tableWidget_Reactions->item(i, 6)->text().toDouble()); // Lkf
            allReactions.last()->setRctCPara(3, ui->tableWidget_Reactions->item(i, 7)->text().toDouble()); // Lkb
            allReactions.last()->setRctCPara(4, ui->tableWidget_Reactions->item(i, 9)->text().toDouble()); // Ukf
            allReactions.last()->setRctCPara(5, ui->tableWidget_Reactions->item(i, 10)->text().toDouble()); // Ukb
        }

        // Enables/disables items for fitting.
        if(checkbox->isChecked() && (flags & XMFLAG_FIT)) {
            ui->tableWidget_Reactions->item(i, 6)->setFlags(ui->tableWidget_Reactions->item(i, 6)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->item(i, 6)->setBackground(QColor(255, 255, 255, 255));
            ui->tableWidget_Reactions->item(i, 6)->setForeground(QColor(0, 0, 0, 255));
            ui->tableWidget_Reactions->item(i, 7)->setFlags(ui->tableWidget_Reactions->item(i, 7)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->item(i, 7)->setBackground(QColor(255, 255, 255, 255));
            ui->tableWidget_Reactions->item(i, 7)->setForeground(QColor(0, 0, 0, 255));
            ui->tableWidget_Reactions->item(i, 9)->setFlags(ui->tableWidget_Reactions->item(i, 9)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->item(i, 9)->setBackground(QColor(255, 255, 255, 255));
            ui->tableWidget_Reactions->item(i, 9)->setForeground(QColor(0, 0, 0, 255));
            ui->tableWidget_Reactions->item(i, 10)->setFlags(ui->tableWidget_Reactions->item(i, 10)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->item(i, 10)->setBackground(QColor(255, 255, 255, 255));
            ui->tableWidget_Reactions->item(i, 10)->setForeground(QColor(0, 0, 0, 255));
            if(allReactions.last()->getType() == RCTFLAG_E) {
                ui->tableWidget_Reactions->item(i, 5)->setFlags(ui->tableWidget_Reactions->item(i, 5)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
                ui->tableWidget_Reactions->item(i, 5)->setBackground(QColor(255, 255, 255, 255));
                ui->tableWidget_Reactions->item(i, 5)->setForeground(QColor(0, 0, 0, 255));
                ui->tableWidget_Reactions->item(i, 8)->setFlags(ui->tableWidget_Reactions->item(i, 8)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
                ui->tableWidget_Reactions->item(i, 8)->setBackground(QColor(255, 255, 255, 255));
                ui->tableWidget_Reactions->item(i, 8)->setForeground(QColor(0, 0, 0, 255));
                ui->tableWidget_Reactions->item(i, 11)->setFlags(ui->tableWidget_Reactions->item(i, 11)->flags() | (Qt::ItemIsSelectable | Qt::ItemIsEditable));
                ui->tableWidget_Reactions->item(i, 11)->setBackground(QColor(255, 255, 255, 255));
                ui->tableWidget_Reactions->item(i, 11)->setForeground(QColor(0, 0, 0, 255));
            }
        } else { // grey out fitting stuff.
            ui->tableWidget_Reactions->item(i, 6)->setFlags(ui->tableWidget_Reactions->item(i, 6)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->item(i, 6)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_Reactions->item(i, 6)->setForeground(QColor(100, 100, 100, 255));
            ui->tableWidget_Reactions->item(i, 7)->setFlags(ui->tableWidget_Reactions->item(i, 7)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->item(i, 7)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_Reactions->item(i, 7)->setForeground(QColor(100, 100, 100, 255));
            ui->tableWidget_Reactions->item(i, 8)->setFlags(ui->tableWidget_Reactions->item(i, 8)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->item(i, 8)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_Reactions->item(i, 8)->setForeground(QColor(100, 100, 100, 255));
            ui->tableWidget_Reactions->item(i, 9)->setFlags(ui->tableWidget_Reactions->item(i, 9)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->item(i, 9)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_Reactions->item(i, 9)->setForeground(QColor(100, 100, 100, 255));
            ui->tableWidget_Reactions->item(i, 10)->setFlags(ui->tableWidget_Reactions->item(i, 10)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->item(i, 10)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_Reactions->item(i, 10)->setForeground(QColor(100, 100, 100, 255));
            ui->tableWidget_Reactions->item(i, 11)->setFlags(ui->tableWidget_Reactions->item(i, 11)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->item(i, 11)->setBackground(QColor(200, 200, 200, 255));
            ui->tableWidget_Reactions->item(i, 11)->setForeground(QColor(100, 100, 100, 255));
        }

        connect(checkbox, SIGNAL(clicked(bool)), this, SLOT(tableReactionsChanged()));

    }

//    qDebug() << "tableReactionsChanged(): was here! <- " << allReactions.last()->getk0();

    connect(ui->tableWidget_Reactions, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableReactionsChanged()));
}

/*
 * Number of reactions has changed.
 */
void Xmarcus::reactionsChanged(int reactions) {

    if(reactions == allReactions.size())
        return;
    else if(reactions > allReactions.size()) {
        disconnect(ui->tableWidget_Reactions, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableReactionsChanged()));

        for(int i=ui->tableWidget_Reactions->rowCount(); i<reactions; i++) {
            ui->tableWidget_Reactions->setRowCount(ui->tableWidget_Reactions->rowCount()+1);

            QWidget *pWidget = new QWidget(); QCheckBox *pCheckBox = new QCheckBox(); QHBoxLayout *pLayout = new QHBoxLayout(pWidget);
            pLayout->addWidget(pCheckBox); pLayout->setAlignment(Qt::AlignCenter); pLayout->setContentsMargins(0,0,0,0); pWidget->setLayout(pLayout);
            ui->tableWidget_Reactions->setCellWidget(i, 0, pWidget);
            QObject::connect(pCheckBox, SIGNAL(clicked(bool)), this, SLOT(tableReactionsChanged()));
            ui->tableWidget_Reactions->setItem(i, 0, new QTableWidgetItem(" "));
            ui->tableWidget_Reactions->setItem(i, 1, new QTableWidgetItem("E"));
            ui->tableWidget_Reactions->item(i, 1)->setFlags(ui->tableWidget_Reactions->item(i, 1)->flags() & ~(Qt::ItemIsSelectable | Qt::ItemIsEditable));
            ui->tableWidget_Reactions->setItem(i, 2, new QTableWidgetItem("A + e ⇄ B"));
            ui->tableWidget_Reactions->setItem(i, 3, new QTableWidgetItem("0.2"));
            ui->tableWidget_Reactions->setItem(i, 4, new QTableWidgetItem("0.5"));
            ui->tableWidget_Reactions->setItem(i, 5, new QTableWidgetItem("1.0e3"));
            ui->tableWidget_Reactions->setItem(i, 6, new QTableWidgetItem("0.1"));
            ui->tableWidget_Reactions->setItem(i, 7, new QTableWidgetItem("0.3"));
            ui->tableWidget_Reactions->setItem(i, 8, new QTableWidgetItem("1.0e-3"));
            ui->tableWidget_Reactions->setItem(i, 9, new QTableWidgetItem("0.3"));
            ui->tableWidget_Reactions->setItem(i, 10, new QTableWidgetItem("0.7"));
            ui->tableWidget_Reactions->setItem(i, 11, new QTableWidgetItem("1.0e6"));
            ui->tableWidget_Reactions->setItemDelegate(new ReactionsTableDelegate());
        }

        connect(ui->tableWidget_Reactions, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableReactionsChanged()));
    } else if(reactions < allReactions.size()) {
        disconnect(ui->tableWidget_Reactions, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableReactionsChanged()));

        for(int i=ui->tableWidget_Reactions->rowCount()-1; i>reactions-1; i--)
            ui->tableWidget_Reactions->removeRow(i);

        connect(ui->tableWidget_Reactions, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(tableReactionsChanged()));
    }

//    qDebug() << "reactionsChanged(" << reactions << ") ...";

    tableReactionsChanged();
}

/*
 * Build reactions matrix and check for TSR
 *
 * Source: W. Luo, S. W. Feldberg, M. Rudolph, Journal of Electroanalytical Chemistry 1994, 368, 109–113.
 */
void Xmarcus::checkReactions(void) {

    gsl_matrix_int_free(M); M = gsl_matrix_int_calloc(size_t(allReactions.size()), size_t(allSpecies.size()));

    for(int i=0; i<allReactions.size(); i++) {
        QStringList fields = allReactions.at(i)->getStrReaction().split(QRegExp("\\s+"), Qt::SkipEmptyParts);
        for(int j=0; j<allSpecies.size(); j++) {
            switch(allReactions.at(i)->getType()) {
            case RCTFLAG_E:
                if(fields.at(0) == allSpecies.at(j)->getSymbol())   // Reagent: Ox
                    gsl_matrix_int_set(M, size_t(i), size_t(j), -1);
                if(fields.at(4) == allSpecies.at(j)->getSymbol())   // Product: Red
                    gsl_matrix_int_set(M, size_t(i), size_t(j), 1);
                break;
            case RCTFLAG_C:
                if(fields.at(0) == allSpecies.at(j)->getSymbol())       // Reagent
                    gsl_matrix_int_set(M, size_t(i), size_t(j), -1);
                if(fields.at(2) == allSpecies.at(j)->getSymbol())  // Product
                    gsl_matrix_int_set(M, size_t(i), size_t(j), 1);
                break;
            case RCTFLAG_C2:
                if(fields.at(0) == allSpecies.at(j)->getSymbol())       // Reagent
                    gsl_matrix_int_set(M, size_t(i), size_t(j), gsl_matrix_int_get(M, size_t(i), size_t(j))-1);
                if(fields.at(2) == allSpecies.at(j)->getSymbol())  // Reagent 2
                    gsl_matrix_int_set(M, size_t(i), size_t(j), gsl_matrix_int_get(M, size_t(i), size_t(j))-1);
                if(fields.at(4) == allSpecies.at(j)->getSymbol())  // Product
                    gsl_matrix_int_set(M, size_t(i), size_t(j), gsl_matrix_int_get(M, size_t(i), size_t(j))+1);
                if(fields.at(6) == allSpecies.at(j)->getSymbol())  // Product 2
                    gsl_matrix_int_set(M, size_t(i), size_t(j), gsl_matrix_int_get(M, size_t(i), size_t(j))+1);
                break;
            case RCTFLAG_CC:
                if(fields.at(0) == allSpecies.at(j)->getSymbol())       // Reagent
                    gsl_matrix_int_set(M, size_t(i), size_t(j), gsl_matrix_int_get(M, size_t(i), size_t(j))-1);
                if(fields.at(2) == allSpecies.at(j)->getSymbol())  // Reagent 2
                    gsl_matrix_int_set(M, size_t(i), size_t(j), gsl_matrix_int_get(M, size_t(i), size_t(j))-1);
                if(fields.at(4) == allSpecies.at(j)->getSymbol())  // Product
                    gsl_matrix_int_set(M, size_t(i), size_t(j), 1);
                break;
            case RCTFLAG_CD:
                if(fields.at(0) == allSpecies.at(j)->getSymbol())       // Reagent
                    gsl_matrix_int_set(M, size_t(i), size_t(j), -1);
                if(fields.at(2) == allSpecies.at(j)->getSymbol())  // Product
                    gsl_matrix_int_set(M, size_t(i), size_t(j), gsl_matrix_int_get(M, size_t(i), size_t(j))+1);
                if(fields.at(4) == allSpecies.at(j)->getSymbol())  // Product 2
                    gsl_matrix_int_set(M, size_t(i), size_t(j), gsl_matrix_int_get(M, size_t(i), size_t(j))+1);
                break;
            default:
                std::cout << "checkReactions(): we should never get here!\n";
            }
        }
    }

//    qDebug() << "checkReactions() ...";

}

/*
 * Check if we need to use non-linear solver.
 */
void Xmarcus::checkNonLinear(void) {

    flags &= ~XMFLAG_NL;

    if((flags & XMFLAG_FRCNL) || (flags & XMFLAG_ADS))
        flags |= XMFLAG_NL;

    for(int i=0; i<allReactions.size(); i++) {
        if((allReactions.at(i)->getType() == RCTFLAG_C2) || (allReactions.at(i)->getType() == RCTFLAG_CC) || (allReactions.at(i)->getType() == RCTFLAG_CD))
            flags |= XMFLAG_NL;
    }

//    qDebug() << "checkNonLinear(): was here...";

}
