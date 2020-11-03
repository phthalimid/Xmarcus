#include "workersimcv.h"

WorkerSimCV::WorkerSimCV(QObject *parent) : QObject(parent), active(false) {

    whoami = this;

    // you could copy data from constructor arguments to internal variables here.

}

// --- DECONSTRUCTOR ---
WorkerSimCV::~WorkerSimCV() {
    // free resources

    gsl_vector_free(zetaj);
    gsl_vector_free(alphai);
    gsl_vector_free(betai);
    gsl_vector_free(c);
    gsl_matrix_free(m);
    gsl_matrix_free(linm);
    gsl_vector_free(oldThetai);

}

void WorkerSimCV::setThreadStatus(QString status) {
    emit sendThreadStatus(status);
}

/*
 * Activates the simulation.
 */
void WorkerSimCV::startSim() {

    QMutexLocker locker(&lock);
    active = !active;

    if(workDone == workAmount) {
        workDone = 0;
        emit sendProgressSim(workDone);
    }
}

// --- PROCESS ---
// Start processing data.
void WorkerSimCV::process() {

    if(active) {
        setThreadStatus("Working");

        simData = xmParent->getSimData();
        fitData = xmParent->getFitData();
        allSpecies = xmParent->getAllSpecies();
        allReactions = xmParent->getAllReactions();
        M = xmParent->getM();

        initialiseAll();

        if(xmParent->getFlags() & XMFLAG_FIT) {

            fitData->splineData(simData->getPotential());
//           fitData->splineData(simData->getTime(), simData->getPotential());
            procFitCV();

        } else {

            procCV();
        }

//        emit finished();
        active = false;

    } else setThreadStatus("    Idle    ");

}

/*
 * f(x) = f(c) = M * x - c = 0
 */
int WorkerSimCV::CV_f (const gsl_vector *x, void *params, gsl_vector *f) {

    Xmarcus *parent = static_cast<Xmarcus*>(params);
    gsl_matrix *linm = parent->getWorker()->getLinm();
    gsl_vector *c = parent->getWorker()->getC();
    QVector<Reaction *> allReactions = parent->getAllReactions();
    QVector<Species *> allSpecies = parent->getAllSpecies();

    int nos = allSpecies.size();
    int nop = parent->getWorker()->getNop();
    double dt = parent->getWorker()->getDt();
    gsl_matrix_int *M = parent->getWorker()->getM();
    double dHelp;

    gsl_vector_set_zero (f);

/*
    for(int i=0; i<nop*nos/5; i++)
        qDebug() << "x[" << i << "] = " << gsl_vector_get(x, i);
    for(int i=(nos*nop+nos)-nop*nos/5; i<(nos*nop+nos); i++)
        qDebug() << "x[" << i << "] = " << gsl_vector_get(x, i);
*/

    /*
     * Heterogeneous boundary
     */
    int react = -1; int prod = -1;
    for (int i=0; i<nos; i++) {
        dHelp = gsl_matrix_get(linm, size_t(i), size_t(i)) * gsl_vector_get(x, size_t(i)) +
                gsl_matrix_get(linm, size_t(i), size_t(nos+i)) * gsl_vector_get(x, size_t(nos+i));
        gsl_vector_set (f, size_t(i), dHelp);
    }
    for (int i=0; i<allReactions.size(); i++) if(allReactions.at(i)->getType() == RCTFLAG_E) {
        for (int j=0; j<nos; j++) {
            if(gsl_matrix_int_get(M, size_t(i), size_t(j)) == -1) { // reactant found
                react = j;
                continue;
            } else if(gsl_matrix_int_get(M, size_t(i), size_t(j)) == 1) { // product found
                prod = j;
                continue;
            }
        }        
        gsl_vector_set(f, size_t(react), gsl_vector_get(f, size_t(react)) +
                       (allReactions.at(i)->getRctCPara(0) * gsl_vector_get(x, size_t(react)) -
                       allReactions.at(i)->getRctCPara(1) * gsl_vector_get(x, size_t(prod)))*dt);
        gsl_vector_set(f, size_t(prod), gsl_vector_get(f, size_t(prod)) +
                       (allReactions.at(i)->getRctCPara(1) * gsl_vector_get(x, size_t(prod)) -
                       allReactions.at(i)->getRctCPara(0) * gsl_vector_get(x, size_t(react)))*dt);
    }

    /*
     * Adsorption: surface concentrations are at the end of the concentrations vector
     */
    if(parent->getFlags() & XMFLAG_ADS) {

        // First calculate theta.
        double theta = 0.0;
        for(int i=0; i<nos; i++) if(allSpecies.at(i)->getAdsorb()) { // Calculate theta-i and theta
            theta += gsl_vector_get(x, size_t(nop*nos+i)) / parent->getQmax();
        }
        theta = theta > 1.0 ? 1.0 : theta;
        theta = theta < 0.0 ? 0.0 : theta;

        /*
         * theta_i-s are at the bottom right of the matrix -> end of the vector.
         *
         * d q_i/dt = k_ai c_i (1-theta) - k_di theta_i
         *
         * dt??????
         *
         * ALL with q_i now...
         */
        for(int i=0; i<nos; i++) if(allSpecies.at(i)->getAdsorb()) {
            gsl_vector_set(f, size_t(nop*nos+i), gsl_vector_get(f, size_t(nop*nos+i)) +
                           (allSpecies.at(i)->getKa() * gsl_vector_get(x, size_t(i)) * (1.0-theta) -
                           allSpecies.at(i)->getKd() * gsl_vector_get(x, size_t(nop*nos+i))/parent->getQmax())*dt);
            gsl_vector_set(f, size_t(i), gsl_vector_get(f, size_t(i)) +
                           (allSpecies.at(i)->getKd() * gsl_vector_get(x, size_t(nop*nos+i))/parent->getQmax() -
                            allSpecies.at(i)->getKa() * gsl_vector_get(x, size_t(i)) * (1.0-theta))*dt);
        } else
            gsl_vector_set(f, size_t(nop*nos+i), gsl_vector_get(x, size_t(nop*nos+i)));

        // Now we are looking for redox-active adsorbed species.
        for (int i=0; i<allReactions.size(); i++) if(allReactions.at(i)->getType() == RCTFLAG_E) {
            for (int j=0; j<nos; j++) {
                if(gsl_matrix_int_get(M, size_t(i), size_t(j)) == -1) { // reactant found
                    react = j;
                    continue;
                } else if(gsl_matrix_int_get(M, size_t(i), size_t(j)) == 1) { // product found
                    prod = j;
                    continue;
                }
            }
            if(allSpecies.at(react)->getAdsorb() && allSpecies.at(prod)->getAdsorb()) {
                gsl_vector_set(f, size_t(nop*nos+react), gsl_vector_get(f, size_t(nop*nos+react)) +
                               (-dt*allReactions.at(i)->getRctCPara(0) * gsl_vector_get(x, size_t(nop*nos+react)) +
                               dt*allReactions.at(i)->getRctCPara(1) * gsl_vector_get(x, size_t(nop*nos+prod))));
                gsl_vector_set(f, size_t(nop*nos+prod), gsl_vector_get(f, size_t(nop*nos+prod)) +
                               (-dt*allReactions.at(i)->getRctCPara(1) * gsl_vector_get(x, size_t(nop*nos+prod)) +
                               dt*allReactions.at(i)->getRctCPara(0) * gsl_vector_get(x, size_t(nop*nos+react))));
            } else if(allSpecies.at(react)->getAdsorb() || allSpecies.at(prod)->getAdsorb()) {
                emit parent->getWorker()->error("If an adsorbed species is electro-active, the reagent/product of the redox reaction must also be adsorbed!");
                parent->getWorker()->setKillThread(true);
                return GSL_FAILURE;
                // Kill thread here.
            }
        }

//        qDebug() << "theta = " << theta << "\tt1 = " << gsl_vector_get(x, 0) << "\tt2 = "
  //               << gsl_vector_get(x, 1) << "\tt3 = " << gsl_vector_get(x, 2);


    }



    /*
     * Diffusion and 1st order homogeneous reactions.
     */
    for(int i=1; i<(nop-1); i++) for(int j=0; j<nos; j++) {
        dHelp = gsl_matrix_get(linm, size_t(i*nos+j), size_t((i-1)*nos+j)) * gsl_vector_get(x, size_t((i-1)*nos+j));
        for(int k=0; k<nos; k++)
            dHelp += gsl_matrix_get(linm, size_t(i*nos+j), size_t(i*nos+k)) * gsl_vector_get(x, size_t(i*nos+k));
        dHelp += gsl_matrix_get(linm, size_t(i*nos+j), size_t((i+1)*nos+j)) * gsl_vector_get(x, size_t((i+1)*nos+j));
        dHelp -= gsl_vector_get(c, size_t(i*nos+j));
        gsl_vector_set(f, size_t(i*nos+j), gsl_vector_get(f, size_t(i*nos+j)) + dHelp);
    }

    /*
     * 2nd order homogeneous reactions:
     * A + B = C + D
     * dc_A/dt = kf c_A c_B - kb c_C c_D
     * A = 2B
     * dc_A/dt = kf c_A - 2 kb c_B c_B
     */
    double kfhelp, kbhelp;

    for(int j=0; j<allReactions.size(); j++) {
        if(allReactions.at(j)->getType() == RCTFLAG_C2 || allReactions.at(j)->getType() == RCTFLAG_CC || allReactions.at(j)->getType() == RCTFLAG_CD) {
            for(int i=1; i<(nop-1); i++) {
                kfhelp = allReactions.at(j)->getRctCPara(0) * dt; // kf
                kbhelp = allReactions.at(j)->getRctCPara(1) * dt; // kb
                for(int k=0; k<nos; k++) { // multiply rate constants with concentrations...
                    if(gsl_matrix_int_get(M, size_t(j), size_t(k)) < 0)
                        kfhelp *= pow(gsl_vector_get(x, size_t(i*nos+k)), -gsl_matrix_int_get(M, size_t(j), size_t(k)));
                    else if(gsl_matrix_int_get(M, size_t(j), size_t(k)) > 0)
                        kbhelp *= pow(gsl_vector_get(x, size_t(i*nos+k)), gsl_matrix_int_get(M, size_t(j), size_t(k)));
                }
                for (int k=0; k<nos; k++) { // add concentration changes...
                    if(gsl_matrix_int_get(M, size_t(j), size_t(k)) < 0)
                        gsl_vector_set(f, size_t(i*nos+k), gsl_vector_get(f, size_t(i*nos+k)) - gsl_matrix_int_get(M, size_t(j), size_t(k))*(kfhelp-kbhelp));
                    else if(gsl_matrix_int_get(M, size_t(j), size_t(k)) > 0)
                        gsl_vector_set(f, size_t(i*nos+k), gsl_vector_get(f, size_t(i*nos+k)) + gsl_matrix_int_get(M, size_t(j), size_t(k))*(-kfhelp+kbhelp));
                }
            }
        }
    }

    /*
     * Semiinfinite boundary
     */
    for(int j=0; j<nos; j++)
        gsl_vector_set(f, size_t((nop-1)*nos+j), gsl_vector_get(x, size_t((nop-2)*nos+j)) - gsl_vector_get(x, size_t((nop-1)*nos+j)));

    return GSL_SUCCESS;
}

/*
 * Function that returns 0 if effective potential is found (for root finding).
 */
double WorkerSimCV::iRC_f(double E, void *params) {

    Xmarcus *parent = static_cast<Xmarcus*>(params);
    gsl_vector *c = parent->getWorker()->getC();
    ECData *simData = parent->getSimData();
    int iData = simData->getIData();

    //struct iRC_params *p = (struct iRC_params *) params;

    double result = 0.0;

    // rescue the concentration vectors.
    gsl_vector *chelp = gsl_vector_calloc(c->size);
    gsl_vector_memcpy(chelp, c);

    parent->getWorker()->solveSSE(E);

    // There is an error in my PhD thesis. Stoerzbach's publication is correct.
    if (iData != 0)
        result = E - (simData->getPotential(iData) + simData->getResistance(0) * (-parent->getWorker()->calcFCurrent()
                 + simData->getDLCapacitance(0) * simData->getEffPotential(iData-1) / (parent->getWorker()->getDt()
                 * simData->getArea(0) * GSL_CONST_MKSA_FARADAY))) / (1.0 + simData->getResistance(0) * simData->getDLCapacitance(0)
                 / (parent->getWorker()->getDt() * simData->getArea(0) * GSL_CONST_MKSA_FARADAY));
    else	// In case i=0, we have to set something else for Eeff[i-1]...
        result = E - (simData->getPotential(iData) + simData->getResistance(0) * (-parent->getWorker()->calcFCurrent()
                 + simData->getDLCapacitance(0) * simData->getPotential(iData) / (parent->getWorker()->getDt()
                 * simData->getArea(0) * GSL_CONST_MKSA_FARADAY))) / (1.0 + simData->getResistance(0) * simData->getDLCapacitance(0)
                 / (parent->getWorker()->getDt() * simData->getArea(0) * GSL_CONST_MKSA_FARADAY));

    gsl_vector_memcpy(c, chelp);
    free (chelp);

    return result;
}

/*
 * Solve the simultaneous system of equations (SSE) for one potential.
 */
void WorkerSimCV::solveSSE(double E) {

    int n = int(c->size);
    int nos = xmParent->getAllSpecies().size();

    // For linear & multiroot solvers.
    gsl_permutation *p = gsl_permutation_alloc (size_t(n));
    gsl_vector *x = gsl_vector_alloc (size_t(n));
    int status;

    calcHeteroKinetics(E);

    if(xmParent->getFlags() & XMFLAG_NL) { // non-linear solver.

        const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
        gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (T, size_t(n));

        size_t iter = 0;

        gsl_multiroot_function f = {&CV_f, size_t(n), static_cast<void *>(xmParent)};

        gsl_multiroot_fsolver_set (s, &f, c);

        if(killThread) {
            emit finished();
            return;
        }

        do {
            iter++;

            if ((status = gsl_multiroot_fsolver_iterate(s)))
                break;

            status = gsl_multiroot_test_residual(s->f, nlerror);
        } while (status == GSL_CONTINUE && iter < size_t(maxiter));

        gsl_vector_memcpy (c, s->x);

//            qDebug() << "status =" << gsl_strerror(status);

        gsl_multiroot_fsolver_free (s);

    } else { // linear problem.
        gsl_matrix_memcpy (m, linm); // do we really need m?

        // Not sure why we need to do this...
        for (int j=0; j<nos; j++) {
            gsl_vector_set (c, size_t(j), 0.0);
            gsl_vector_set (c, size_t((nop-1)*nos+j), 0.0); // Semiinfinite boundary condition
        }

        // Add the kf and kbs
        int react = -1; int prod = -1;
        for(int i=0; i<allReactions.size(); i++) if(allReactions.at(i)->getType() == RCTFLAG_E) {
            for (int j=0; j<nos; j++) {
                if(gsl_matrix_int_get(M, size_t(i), size_t(j)) == -1) { /* reactant found */
                    react = j;
                    continue;
                } else if(gsl_matrix_int_get(M, size_t(i), size_t(j)) == 1) { /* product found */
                    prod = j;
                    continue;
                }
            }
            gsl_matrix_set(m, size_t(react), size_t(react), gsl_matrix_get(m, size_t(react), size_t(react)) + allReactions.at(i)->getRctCPara(0)*dt);
            gsl_matrix_set(m, size_t(react), size_t(prod), gsl_matrix_get(m, size_t(react), size_t(prod)) - allReactions.at(i)->getRctCPara(1)*dt);
            gsl_matrix_set(m, size_t(prod), size_t(react), gsl_matrix_get(m, size_t(prod), size_t(react)) - allReactions.at(i)->getRctCPara(0)*dt);
            gsl_matrix_set(m, size_t(prod), size_t(prod), gsl_matrix_get(m, size_t(prod), size_t(prod)) + allReactions.at(i)->getRctCPara(1)*dt);
        }

        // Solve simultaneous equations...
        gsl_linalg_LU_decomp (m, p, &status);
        gsl_linalg_LU_solve (m, p, c, x);
        gsl_vector_memcpy (c, x);
    }

    gsl_vector_free (x);
    gsl_permutation_free (p);

}

/*
 * CV least squares: we return the sum of the squared residuals S = sum (y_i - f_i(params))^2
 * In our "function", we want to minimise the difference between experiement and simulation!
 */
double WorkerSimCV::FitCV_f(const gsl_vector *x, void *params) {

    double S = 0.0;

    Xmarcus *parent = static_cast<Xmarcus*>(params);

    parent->getWorker()->updateFitVariables(x);
    parent->getWorker()->procCV();
/*
    std::ofstream tstFile;
    tstFile.open("/Users/nannth/Desktop/testfile.csv");

    for(int i=0; i<parent->getSimData()->getCurrent().size(); i++)
        tstFile << parent->getSimData()->getPotential(i) << ",\t"
                << parent->getFitData()->getCurrent(i) << ",\t"
                << parent->getFitData()->getSplineCurrent(i) << ",\t"
                << parent->getSimData()->getCurrent(i) << "\n";

    tstFile.close();
*/
    for(int i=0; i<parent->getSimData()->getCurrent().size(); i++)
        S += (parent->getFitData()->getSplineCurrent(i) - parent->getSimData()->getCurrent(i))*(parent->getFitData()->getSplineCurrent(i) - parent->getSimData()->getCurrent(i));

    qDebug() << "S = " << S;

    return S;
}

/*
 * Main function for CV simulation.
 */
void WorkerSimCV::procCV() {

    double Eeff, current;
//    int counter = 0;

    // Reset stuff in case of repeated simulations.
    simData->clearCurrent();
    resetConcVector();

    auto start_time = std::chrono::steady_clock::now(); // Start time.
    workAmount = simData->getPotential().size();
    workDone = 0;

    for(int i=0; i<simData->getPotential().size(); i++) {

        if(xmParent->getFlags() & XMFLAG_IRC) {
            int iter = 0, max_iter = 200, status;
            double x_hi, x_lo;

            // Pass on current datapoint to class.
            simData->setIData(i);

            gsl_function F;
            F.function = &iRC_f;
            F.params = xmParent;

            const gsl_root_fsolver_type *T;
            gsl_root_fsolver *s;
            T = gsl_root_fsolver_brent;
            s = gsl_root_fsolver_alloc(T);

            // Set the boundaries for the root search
            if (i != 0) {
                x_lo = simData->getEffPotential(i-1);
                x_hi = x_lo + IRCLIMIT * simData->getPotentialSteps();

                while(iRC_f(x_lo, xmParent) >= 0.0)
                    x_lo -= simData->getPotentialSteps();
                while(iRC_f(x_hi, xmParent) <= 0.0)
                    x_hi += simData->getPotentialSteps();
            } else {	// i = 0
                double dSign = sgn(simData->getPotential(1) - simData->getPotential(0)); // >= 0.0 ? 1.0 : -1.0;

                if(dSign > 0.0) { // pos slope
                    x_lo = simData->getPotential(0) - IRCLIMIT * simData->getPotentialSteps();
                    x_hi = simData->getPotential(1) + IRCLIMIT * simData->getPotentialSteps();
                } else { // neg slope
                    x_lo = simData->getPotential(1) - IRCLIMIT * simData->getPotentialSteps();
                    x_hi = simData->getPotential(0) + IRCLIMIT * simData->getPotentialSteps();
                }
            }

            gsl_root_fsolver_set (s, &F, x_lo, x_hi);

            do {
                iter++;

                status = gsl_root_fsolver_iterate(s);
                Eeff = gsl_root_fsolver_root(s);
                x_lo = gsl_root_fsolver_x_lower(s);
                x_hi = gsl_root_fsolver_x_upper(s);
                status = gsl_root_test_interval(x_lo, x_hi, 0, 1.0e-4);
/*
                if (status == GSL_SUCCESS)
                    qDebug() << "Converged:\n";
                qDebug() << iter << x_lo << x_hi << simData->getPotential(i) << "--->" << Eeff << x_hi - x_lo;
*/
            } while (status == GSL_CONTINUE && iter < max_iter);

            simData->setEffPotential(i, Eeff);

            gsl_root_fsolver_free(s);
        } else
            Eeff = simData->getPotential(i);

        solveSSE(Eeff);

        // Calculate and add up currents
        current = calcFCurrent();
        if(xmParent->getFlags() & XMFLAG_ADS)
            current += calcAdsCurrent();
        if((xmParent->getFlags() & XMFLAG_IRC) && (simData->getDLCapacitance(0) != 0.0)) {
            if(i==0)
                current += sgn(simData->getPotential(1) - simData->getPotential(0)) * calcCCurrent();
            else
                current += sgn(simData->getPotential(i) - simData->getPotential(i-1)) * calcCCurrent();
        }
        simData->appendCurrent(current);

        // Emit progress signal every second.
        if(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_time).count() > 0.2) {
            emit sendProgressSim(i*100/simData->getPotential().size());
            start_time = std::chrono::steady_clock::now();
        }

        workDone++;
    }

    emit sendFinishedCV(simData);
}

/*
 * Main function for CV fitting.
 */
void WorkerSimCV::procFitCV() {

    // Fitting constants
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = nullptr;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    double size;

    // !!!!!!! iData somewhere  !!!!!!!1


    QVector<double> qvSize;
    QVector<double> qvFunc;

    int status;
    int iter = 0;

    // Allocate and fill initial variable vector.
    x = calcFitVariables();

    // Allocate and set initial step size.
    ss = setInitStepSize();

    // Initialize method and iterate
    minex_func.n = static_cast<size_t>(fitNoVariables);
    minex_func.f = &FitCV_f;
    minex_func.params = static_cast<void *>(xmParent);

    s = gsl_multimin_fminimizer_alloc (T, static_cast<size_t>(fitNoVariables));
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do {
        iter++;

        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
          break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-9);

        if (status == GSL_SUCCESS)
            qDebug() << "converged to minimum at";
        qDebug() << "iteration: " << iter << ", x[0] = " << gsl_vector_get(s->x, 0) << ", f() = " << s->fval << ", size = " << size;

        qvSize.append(size);
        qvFunc.append(s->fval);

        // Plot size and f() here.


    } while (status == GSL_CONTINUE && iter < 1000);

//    qDebug() << "procFitCV(): .... was here ...";

    resultFitVariables(s->x);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
}

/*
 * Creates a GSL vector for the fitting parameters and fills it.
 */
gsl_vector *WorkerSimCV::calcFitVariables() {

    fitNoVariables = 0;
    fitFlags = 0;

    if(xmParent->isAreaFitChecked()) {fitFlags |= FITFLAG_AREA; fitNoVariables++;}
    if(xmParent->isResFitChecked()) {fitFlags |= FITFLAG_RES; fitNoVariables++;}
    if(xmParent->isDLCapFitChecked()) {fitFlags |= FITFLAG_DLC; fitNoVariables++;}
    for(int i=0; i<allReactions.size(); i++) if(xmParent->isReactFitChecked(i)) {
        fitFlags |= FITFLAG_REACT;
        if(allReactions[i]->getType() == RCTFLAG_E)
            fitNoVariables += 3;
        else
            fitNoVariables += 2;
    }

    // Now we are going to allocate the memory and fill in initial values.
    gsl_vector *x = gsl_vector_alloc(static_cast<size_t>(fitNoVariables));

    size_t it = 0;
    if(fitFlags & FITFLAG_AREA) {gsl_vector_set(x, it, simData->getArea(0)); it++;}
    if(fitFlags & FITFLAG_RES) {gsl_vector_set(x, it, simData->getResistance(0)); it++;}
    if(fitFlags & FITFLAG_DLC) {gsl_vector_set(x, it, simData->getDLCapacitance(0)); it++;}
    if(fitFlags & FITFLAG_REACT) {
        for(int i=0; i<allReactions.size(); i++) if(xmParent->isReactFitChecked(i)) {
            if(allReactions[i]->getType() == RCTFLAG_E) {
                gsl_vector_set(x, it, allReactions[i]->getRctEPara(0)); it++; // E0
                gsl_vector_set(x, it, allReactions[i]->getRctEPara(1)); it++; // alpha
                gsl_vector_set(x, it, allReactions[i]->getRctEPara(2)); it++; // k0
            } else {
                gsl_vector_set(x, it, allReactions[i]->getRctCPara(0)); it++; // kf
                gsl_vector_set(x, it, allReactions[i]->getRctCPara(1)); it++; // kb
            }
        }
    }

    return x;
}

/*
 * Update simulation parameters with fit variables "x".
 */
void WorkerSimCV::updateFitVariables(const gsl_vector *x) {

    size_t it = 0;
    if(fitFlags & FITFLAG_AREA) {
        if(gsl_vector_get(x, it) < simData->getArea(1)) gsl_vector_set(const_cast<gsl_vector *>(x), it, simData->getArea(1));
        if(gsl_vector_get(x, it) > simData->getArea(2)) gsl_vector_set(const_cast<gsl_vector *>(x), it, simData->getArea(2));
        simData->setArea(0, gsl_vector_get(x, it));
        qDebug() << "A(rea) = " << simData->getArea(0) * 10000.0;
        it++;
    }
    if(fitFlags & FITFLAG_RES) {
        if(gsl_vector_get(x, it) < simData->getResistance(1)) gsl_vector_set(const_cast<gsl_vector *>(x), it, simData->getResistance(1));
        if(gsl_vector_get(x, it) > simData->getResistance(2)) gsl_vector_set(const_cast<gsl_vector *>(x), it, simData->getResistance(2));
        simData->setResistance(0, gsl_vector_get(x, it));
        qDebug() << "R = " << simData->getResistance(0);
        it++;
    }
    if(fitFlags & FITFLAG_DLC) {
        if(gsl_vector_get(x, it) < simData->getDLCapacitance(1)) gsl_vector_set(const_cast<gsl_vector *>(x), it, simData->getDLCapacitance(1));
        if(gsl_vector_get(x, it) > simData->getDLCapacitance(2)) gsl_vector_set(const_cast<gsl_vector *>(x), it, simData->getDLCapacitance(2));
        simData->setDLCapacitance(0, gsl_vector_get(x, it));
        qDebug() << "C = " << simData->getDLCapacitance(0);
        it++;
    }
    if(fitFlags & FITFLAG_REACT) {
        for(int i=0; i<allReactions.size(); i++) if(xmParent->isReactFitChecked(i)) {
            if(allReactions[i]->getType() == RCTFLAG_E) {
                if(gsl_vector_get(x, it) < allReactions[i]->getRctEPara(3)) gsl_vector_set(const_cast<gsl_vector *>(x), it, allReactions[i]->getRctEPara(3));
                if(gsl_vector_get(x, it) > allReactions[i]->getRctEPara(6)) gsl_vector_set(const_cast<gsl_vector *>(x), it, allReactions[i]->getRctEPara(6));
                allReactions[i]->setRctEPara(0, gsl_vector_get(x, it));
                qDebug() << "E0 = " << allReactions[i]->getRctEPara(0);
                it++;
                if(gsl_vector_get(x, it) < allReactions[i]->getRctEPara(4)) gsl_vector_set(const_cast<gsl_vector *>(x), it, allReactions[i]->getRctEPara(4));
                if(gsl_vector_get(x, it) > allReactions[i]->getRctEPara(7)) gsl_vector_set(const_cast<gsl_vector *>(x), it, allReactions[i]->getRctEPara(7));
                allReactions[i]->setRctEPara(1, gsl_vector_get(x, it));
                if(xmParent->getFlags() & XMFLAG_BV) qDebug() << "α = " << allReactions[i]->getRctEPara(1);
                else qDebug() << "λ = " << allReactions[i]->getRctEPara(1);
                it++;
                if(gsl_vector_get(x, it) < allReactions[i]->getRctEPara(5)) gsl_vector_set(const_cast<gsl_vector *>(x), it, allReactions[i]->getRctEPara(5));
                if(gsl_vector_get(x, it) > allReactions[i]->getRctEPara(8)) gsl_vector_set(const_cast<gsl_vector *>(x), it, allReactions[i]->getRctEPara(8));
                allReactions[i]->setRctEPara(2, gsl_vector_get(x, it));
                qDebug() << "k0 = " << allReactions[i]->getRctEPara(2);
                it++;
            } else {
                if(gsl_vector_get(x, it) < allReactions[i]->getRctCPara(2)) gsl_vector_set(const_cast<gsl_vector *>(x), it, allReactions[i]->getRctCPara(2)); // Lkf
                if(gsl_vector_get(x, it) > allReactions[i]->getRctCPara(4)) gsl_vector_set(const_cast<gsl_vector *>(x), it, allReactions[i]->getRctCPara(4)); // Ukf
                allReactions[i]->setRctCPara(0, gsl_vector_get(x, it));
                qDebug() << "kf = " << allReactions[i]->getRctCPara(0);
                it++;
                if(gsl_vector_get(x, it) < allReactions[i]->getRctCPara(3)) gsl_vector_set(const_cast<gsl_vector *>(x), it, allReactions[i]->getRctCPara(3)); // Lkb
                if(gsl_vector_get(x, it) > allReactions[i]->getRctCPara(5)) gsl_vector_set(const_cast<gsl_vector *>(x), it, allReactions[i]->getRctCPara(5)); // Ukb
                allReactions[i]->setRctCPara(1, gsl_vector_get(x, it));
                qDebug() << "kb = " << allReactions[i]->getRctCPara(1);
                it++;
            }
        }
    }
}

/*
 * Allocates memory and sets initial step-size (INITSS-th fraction of the total range).
 */
gsl_vector *WorkerSimCV::setInitStepSize() {

    gsl_vector *ss = gsl_vector_alloc(static_cast<size_t>(fitNoVariables));

    size_t it = 0;
    if(fitFlags & FITFLAG_AREA) {gsl_vector_set(ss, it, (simData->getArea(2)-simData->getArea(1))/INITSS); it++;}
    if(fitFlags & FITFLAG_RES) {gsl_vector_set(ss, it, (simData->getResistance(2)-simData->getResistance(1))/INITSS); it++;}
    if(fitFlags & FITFLAG_DLC) {gsl_vector_set(ss, it, (simData->getDLCapacitance(2)-simData->getDLCapacitance(1))/INITSS); it++;}
    if(fitFlags & FITFLAG_REACT) {
        for(int i=0; i<allReactions.size(); i++) if(xmParent->isReactFitChecked(i)) {
            if(allReactions[i]->getType() == RCTFLAG_E) {
                gsl_vector_set(ss, it, (allReactions[i]->getRctEPara(6)-allReactions[i]->getRctEPara(3))/INITSS); it++; // E0
                gsl_vector_set(ss, it, (allReactions[i]->getRctEPara(7)-allReactions[i]->getRctEPara(4))/INITSS); it++; // alpha, lambda
                gsl_vector_set(ss, it, pow(10, log10(allReactions[i]->getRctEPara(2)) +  // k0
                               (log10(allReactions[i]->getRctEPara(8))-log10(allReactions[i]->getRctEPara(5)))/INITSS) - allReactions[i]->getRctEPara(2)); it++;
            } else {
                gsl_vector_set(ss, it, pow(10, log10(allReactions[i]->getRctCPara(0)) +  // kf
                                           (log10(allReactions[i]->getRctCPara(4))-log10(allReactions[i]->getRctCPara(2)))/INITSS) - allReactions[i]->getRctCPara(0)); it++;
                gsl_vector_set(ss, it, pow(10, log10(allReactions[i]->getRctCPara(1)) +  // kb
                                           (log10(allReactions[i]->getRctCPara(5))-log10(allReactions[i]->getRctCPara(3)))/INITSS) - allReactions[i]->getRctCPara(1)); it++;
            }
        }
    }

    return ss;
}

/*
 * updates the parameters with the fitting results.
 */
void WorkerSimCV::resultFitVariables(gsl_vector *x) {

    size_t it = 0;
    if(fitFlags & FITFLAG_AREA) {xmParent->setUiArea(gsl_vector_get(x, it)); it++;}
    if(fitFlags & FITFLAG_RES) {xmParent->setUiResistance(gsl_vector_get(x, it)); it++;}
    if(fitFlags & FITFLAG_DLC) {xmParent->setUiDLCapacitance(gsl_vector_get(x, it)); it++;}
    if(fitFlags & FITFLAG_REACT) {
        for(int i=0; i<allReactions.size(); i++) if(xmParent->isReactFitChecked(i)) {
            xmParent->setUiReactParam(i, 3, gsl_vector_get(x, it)); it++;
            xmParent->setUiReactParam(i, 4, gsl_vector_get(x, it)); it++;
            if(allReactions[i]->getType() == RCTFLAG_E) {xmParent->setUiReactParam(i, 5, gsl_vector_get(x, it)); it++;}
        }
    }
}

/*
 * Update numerical parameters.
 */
void WorkerSimCV::updateNumParams() {

    beta = xmParent->getNumParams()->getBeta();
    Dm = xmParent->getNumParams()->getDm();
    nlerror = xmParent->getNumParams()->getNLerror();
    maxiter = xmParent->getNumParams()->getMaxiter();
    if(xmParent->getNumParams()->getNLbox())
        xmParent->setFlags(xmParent->getFlags() | XMFLAG_FRCNL);
    else
        xmParent->setFlags(xmParent->getFlags() & ~XMFLAG_FRCNL);

}

/*
 * Reset the concentration vector.
 */
void WorkerSimCV::resetConcVector() {

    int nos = allSpecies.size();

    for (int j=0; j<nos; j++) {
        for (int i=0; i<nop; i++)
            gsl_vector_set(c, static_cast<size_t>(i*nos+j), allSpecies.at(j)->getCini());
    }
}

/*
 * Get ready for simulation: initialise everything.
 */
void WorkerSimCV::initialiseAll() {

    // Update the numerical parameters.
    updateNumParams();

    int nos = allSpecies.size(); // just a little helper.

    /*
     * Find Dmax, xmax and dx according to:
     *
     * [1] A. W. Bott, Current Separations 2000, 19, 45-48.
     */
    Dmax = allSpecies[0]->getD();
    for(int i=1; i<nos;i++)
        Dmax = allSpecies[i]->getD() > Dmax ? allSpecies[i]->getD() : Dmax;

    xmax = 6.0 * sqrt (Dmax * simData->getRTimes().last());

    // Now, dx:
    double mu = 0.0;
    for (int i=0; i<allReactions.size(); i++)
        if(allReactions.at(i)->getType() != RCTFLAG_E)
            mu = mu > sqrt(Dmax/(allReactions.at(i)->getRctCPara(0) + allReactions.at(i)->getRctCPara(1))) ?
                        mu : sqrt(Dmax/(allReactions.at(i)->getRctCPara(0) + allReactions.at(i)->getRctCPara(1)));

    dt = simData->getPotentialSteps()/simData->getScanRate();
    if (mu != 0.0)
        dx = sqrt(mu * mu * 0.02);
    else
        dx = sqrt(Dmax/Dm * dt);

    /*
     * Calculate number of grid points according to:
     *
     * [1] S. W. Feldberg, Journal of Electroanalytical Chemistry 1981, 127, 1-10.
     */
    if(beta == 0.0)
        nop = int(xmax/dx) + 2;
    else
        nop = int(log(xmax/dx * (exp(beta) - 1.0) + 1.0)/beta + 2.0);

    // Find out if there is anything adsorbing.
    int n = nos * nop;  // number of matrix elements.
    for(int i=0; i<nos; i++) if(allSpecies.at(i)->getAdsorb()) {
        xmParent->setFlags(xmParent->getFlags() | XMFLAG_ADS);
        n += nos;
        break;
    }

    /*
     * Calculate Lagrange factors for current calculation:
     *
     * [1] M. Störzbach, J. Heinze, Journal of Electroanalytical Chemistry 1993, 346, 1–27.
     *
     * Used the better method by Neil Dodgson. Attention: derivative at x = 0 NOT x[0]!
     */
    gsl_vector_free(zetaj);
    zetaj = gsl_vector_calloc (LDEGREE);

    double dpij, pij;
    for (int j=0; j<LDEGREE; j++) {
        dpij = 0.0;
        for (int k=0; k<LDEGREE; k++) {
            if (k!=j) {
                pij = 1.0;
                for (int i=0; i<LDEGREE; i++)
                    if ((i!=k) && (i!=j))
                        pij *= -xi(i);
                dpij += pij;
            }
        }

        pij = 1.0;
        for (int i=0; i<LDEGREE; i++)
            if (i!=j)
                pij *= (xi(j) - xi(i));
        gsl_vector_set(zetaj, size_t(j), dpij/pij);
    }

    /*
     * Initialise diffusion coefficient vectors according to:
     *
     * [1] M. Rudolph, J. Electroanal. Chem. 1991, 314, 13-22.
     */
    gsl_vector_free(alphai); gsl_vector_free(betai);
    alphai = gsl_vector_calloc(size_t(n)); // bit too long, because it contains the ads species...
    betai = gsl_vector_calloc(size_t(n));

    for (int i=1; i<nop; i++) { // why start at 1?? Electrode?
        for (int j=0; j<nos; j++) {
                if (i==1)
                    gsl_vector_set(betai, size_t(nos+j), allSpecies.at(j)->getD()*dt/(dx*dx) * (exp(beta)-1.0)/(exp(beta*0.5)-1.0));
                else
                    gsl_vector_set(betai, size_t(i*nos+j), allSpecies.at(j)->getD()*dt/(dx*dx) * exp(2.0*beta*(1.25-i)));
                gsl_vector_set(alphai, size_t(i*nos+j), allSpecies.at(j)->getD()*dt/(dx*dx) * exp(2.0*beta*(0.75-i)));
        }
    }

    /*
     * Set up concentration vector with initial concentrations.
     */
    gsl_vector_free(c);
    c = gsl_vector_alloc(size_t(n)); // again... has the q_is for adsorption at the end.

    for (int j=0; j<nos; j++) {
        for (int i=0; i<nop; i++)
            gsl_vector_set(c, size_t(i*nos+j), allSpecies.at(j)->getCini());
    }

    /*
     * Allocate and initialise linear elements of matrices for both, linear and non-linear.
     */
    m = gsl_matrix_calloc (size_t(n), size_t(n));
    linm = gsl_matrix_calloc (size_t(n), size_t(n));

    /*
     * In case of adsorption: initialise q_i. All species are added to vector/matrix.
     *
     * q_i = theta_i * qmax
     */
    // CAREFUL!!! kd must not be zero!

 //   qDebug() << "qmax = " << xmParent->getQmax();

    double dHelp = 0.0;
    if(xmParent->getFlags() & XMFLAG_ADS) {
        gsl_vector_free(oldThetai); oldThetai = gsl_vector_calloc(size_t(nos));
        for(int i=0; i<nos; i++) if(allSpecies.at(i)->getAdsorb())
            dHelp += allSpecies.at(i)->getKa()/allSpecies.at(i)->getKd() * allSpecies.at(i)->getCini();
        for(int i=0; i<nos; i++) {
            if(allSpecies.at(i)->getAdsorb()) {
                gsl_vector_set(c, size_t(nop*nos+i), xmParent->getQmax() * (allSpecies.at(i)->getKa()/allSpecies.at(i)->getKd() * allSpecies.at(i)->getCini()/(1.0+dHelp)));
            } else {
                gsl_vector_set(c, size_t(nop*nos+i), 0.0);
                gsl_matrix_set(linm, size_t(nop*nos+i), size_t(nop*nos+i), 1.0); // Initialise thetas with 1.0.
            }
            gsl_vector_set(oldThetai, size_t(i), gsl_vector_get(c, size_t(nop*nos+i)));
        }
    }

    xmParent->setQmax();

    /*
     * Set up heterogeneous boundary conditions (just D/dx part)
     */
    double dx1 = dx * (exp(beta*0.5) - 1.0)/(exp(beta) - 1.0);

    for (int j=0; j<nos; j++) {
        gsl_matrix_set(linm, size_t(j), size_t(j), allSpecies.at(j)->getD()/dx1);
        gsl_matrix_set(linm, size_t(j), size_t(j+nos), -allSpecies.at(j)->getD()/dx1);
    }

    /*
     * Set up diffusion matrix according to:
     *
     * [1] M. Rudolph, Journal of Electroanalytical Chemistry 1991, 314, 13-22.
     */
    for(int i=1; i<(nop-1); i++) {
        for(int j=0; j<nos; j++) {
            gsl_matrix_set(linm, size_t(i*nos+j), size_t((i-1)*nos+j), -gsl_vector_get(betai, size_t(i*nos+j)));
            gsl_matrix_set(linm, size_t(i*nos+j), size_t(i*nos+j), 1.0 + gsl_vector_get(alphai, size_t(i*nos+j)) +
                           gsl_vector_get(betai, size_t(i*nos+j)));
            gsl_matrix_set(linm, size_t(i*nos+j), size_t((i+1)*nos+j), -gsl_vector_get(alphai, size_t(i*nos+j)));
        }
    }

    /*
     * Set up first order homogeneous reactions.
     */
    int react = 0; int prod = 0;
    for (int i=0; i<allReactions.size(); i++) if(allReactions.at(i)->getType() == RCTFLAG_C) { /* 1st order homogeneous reaction found */
        for (int j=0; j<nos; j++) {
            if(gsl_matrix_int_get(M, size_t(i), size_t(j)) == -1) { /* reactant found */
                react = j;
                continue;
            } else if(gsl_matrix_int_get(M, size_t(i), size_t(j)) == 1) { /* product found */
                prod = j;
                continue;
            }
        }
        for (int k=1; k<(nop-1); k++) {
            gsl_matrix_set(linm, size_t(k*nos+react), size_t(k*nos+react),
                           gsl_matrix_get(linm, size_t(k*nos+react), size_t(k*nos+react)) + allReactions.at(i)->getRctCPara(0)*dt);
            gsl_matrix_set(linm, size_t(k*nos+react), size_t(k*nos+prod),
                           gsl_matrix_get(linm, size_t(k*nos+react), size_t(k*nos+prod)) - allReactions.at(i)->getRctCPara(1)*dt);
            gsl_matrix_set(linm, size_t(k*nos+prod), size_t(k*nos+prod),
                           gsl_matrix_get(linm, size_t(k*nos+prod), size_t(k*nos+prod)) + allReactions.at(i)->getRctCPara(1)*dt);
            gsl_matrix_set(linm, size_t(k*nos+prod), size_t(k*nos+react),
                           gsl_matrix_get(linm, size_t(k*nos+prod), size_t(k*nos+react)) - allReactions.at(i)->getRctCPara(0)*dt);
        }
    }

    /*
     * Set semiinfinite boundary conditions.
     */
    for (int j=0; j<nos; j++) {
        gsl_matrix_set(linm, size_t((nop-1)*nos+j), size_t((nop-1)*nos+j), -1.0);
        gsl_matrix_set(linm, size_t((nop-1)*nos+j), size_t((nop-2)*nos+j), 1.0);
    }

}

/*
 * Integral function 1 for Marcus-Hush.
 */
double WorkerSimCV::intgl1 (double x, void *p) {

    struct mh_params *params = static_cast<struct mh_params*>(p);
    double help;

    help = exp(-(x-params->Estar)*(x-params->Estar)/(4.0*params->lambdastar));

    return help/(2.0*cosh(x/2.0));
}

/*
 * Integral function 2 for Marcus-Hush.
 */
double WorkerSimCV::intgl2 (double x, void *p) {

    struct mh_params *params = static_cast<struct mh_params*>(p);
    double help;

    help = exp(-x*x/(4.0*params->lambdastar));

    return help/(2.0*cosh(x/2.0));
}

/*
 * Calculates the heterogeneous rate constants at a given potential.
 */
void WorkerSimCV::calcHeteroKinetics(double E) {

    QVector<Reaction *> allReactions = xmParent->getAllReactions();
    double T = xmParent->getSimData()->getTemperature();

    struct mh_params params;
    gsl_function func;
    func.params = &params;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
    double	result, result2, error;

    for(int i=0; i<allReactions.size(); i++) if(allReactions.at(i)->getType() == RCTFLAG_E) {
        if(xmParent->getFlags() & XMFLAG_BV) {
            allReactions[i]->setRctCPara(0, allReactions[i]->getRctEPara(2)/100.0 * exp(-allReactions[i]->getRctEPara(1)*GSL_CONST_MKSA_FARADAY/(GSL_CONST_MKSA_MOLAR_GAS*T) * (E-allReactions[i]->getRctEPara(0))));
            allReactions[i]->setRctCPara(1, allReactions[i]->getRctEPara(2)/100.0 * exp((1.0-allReactions[i]->getRctEPara(1))*GSL_CONST_MKSA_FARADAY/(GSL_CONST_MKSA_MOLAR_GAS*T) * (E-allReactions[i]->getRctEPara(0))));
        } else if(xmParent->getFlags() & XMFLAG_MH) {
            /*
             * We are using the Feldberg method for calculating the kinetics:
             * [1] S. W. Feldberg, Analytical Chemistry 2010, 82, 5176–5183.
             */
            params.Estar = (E - allReactions[i]->getRctEPara(0)) * (GSL_CONST_MKSA_ELECTRON_VOLT/(GSL_CONST_MKSA_BOLTZMANN * T));
            params.lambdastar = allReactions[i]->getRctEPara(1) / (GSL_CONST_MKSA_BOLTZMANN * T);

            func.function = &intgl1;
            gsl_integration_qagi (&func, 0, 1e-7, 1000, w, &result, &error);
            func.function = &intgl2;
            gsl_integration_qagi (&func, 0, 1e-7, 1000, w, &result2, &error);

            allReactions[i]->setRctCPara(0, allReactions[i]->getRctEPara(2) * exp(-params.Estar/2.0) * result/result2);
            allReactions[i]->setRctCPara(1, allReactions[i]->getRctEPara(2) * exp(params.Estar/2.0) * result/result2);
        } else {
            qDebug() << "calcHeteroKinetics(): we should never get here!";
        }
    }
}

/*
 * Calculate Faradayic current.
 */
double WorkerSimCV::calcFCurrent(void) {

    QVector<Species *> allSpecies = xmParent->getAllSpecies();
    double 	current = 0.0;

    for (int i=0; i<allSpecies.size(); i++)
        for (int j=0; j<LDEGREE; j++)
            current -= allSpecies.at(i)->getCharge() * allSpecies.at(i)->getD() *
                    gsl_vector_get(zetaj, size_t(j)) * gsl_vector_get(c, size_t(j*allSpecies.size()+i));

    return -xmParent->getSimData()->getArea(0) * current * GSL_CONST_MKSA_FARADAY;
}

/*
 * Calculate capacitive current.
 */
double WorkerSimCV::calcCCurrent() {
    double current = simData->getScanRate() * simData->getDLCapacitance(0);
    double tau = simData->getResistance(0) * simData->getDLCapacitance(0);
    //double thelp; // = simData->getTime(simData->getIData());

    int j = simData->getRTimes().size() - 1;
    for(; simData->getRTimes(j)>simData->getTime(simData->getIData()); j--) {}
    double thelp = simData->getRTimes(j);

//    qDebug() << "thelp = " << thelp << "\t t - thelp = " << simData->getTime(simData->getIData())-thelp;

    if (simData->getTime(simData->getIData()) > simData->getRTimes(1))
        current *= 1.0 - 2.0 * exp(-(simData->getTime(simData->getIData())-thelp)/tau);
    else
        current *= 1.0 - exp(-simData->getTime(simData->getIData())/tau);

    return current; // * GSL_CONST_MKSA_FARADAY * ecdata->A;
}

/*
 * Calculate adsorption current.
 */
double WorkerSimCV::calcAdsCurrent(void) {

    QVector<Species *> allSpecies = xmParent->getAllSpecies();
    int nos = allSpecies.size();
    double 	current = 0.0;

    for(int i=0; i<nos; i++) if(allSpecies.at(i)->getAdsorb()) {
        current += allSpecies.at(i)->getCharge() * (gsl_vector_get(c, size_t(nop*nos+i)) - gsl_vector_get(oldThetai, size_t(i)))/dt;
        gsl_vector_set(oldThetai, size_t(i), gsl_vector_get(c, size_t(nop*nos+i)));
    }

    return xmParent->getSimData()->getArea(0) * GSL_CONST_MKSA_FARADAY * current;// * xmParent->getQmax();
}
