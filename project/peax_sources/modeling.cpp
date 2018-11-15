#include "modeling.h"

bool intersection(int* A, int* B){
    // A's bottom below B's
    if (A[0] > B[2]) return false;
    // A's top above B's bottom
    if (A[2] < B[0]) return false;
    // A's Left Edge to left of B's right edge
    if (A[1] > B[3]) return false;
    // A's right edge to right of B's left edge
    if (A[3] < B[1]) return false;

    return true;
}




void em_estimate_2D(IMSMatrix* ims_matrix, uint C, double** datenArr, int* box, double min_peak_signal){
    uint r, t, rj, tj, j, l, tr, tt;

    uint R = box[2] - box[0] + 1;
    uint T = box[3] - box[1] + 1;
    
    matrix_t* matrix = ims_matrix->data;
    

    // -1 wegen HGM
    double HGM = R * T;
    int TC = T * C;
    int RC = R * C;
    double tmpWeight[1];
    double*** mbs = 0;


    // lineares Array für Punkte
    double sumVals = 0;
    double *matrixArr = new double[(int)HGM];
    int hgm = 0;
    
    bool stop_modeling = true;
    
    for_loop (r, R){
        for_loop (t, T){
            sumVals += matrixArr[hgm++] = matrix->at(box[0] + r)->at(box[1] + t);
            if (matrix->at(box[0] + r)->at(box[1] + t) > min_peak_signal) stop_modeling = false;
        }
    }
    
    
    double *sumRet = new double[TC];
    double *sumDrift = new double[RC];
    
    
    // lineares Array für Parameter
    for (j = 0; j < C; j++){
        datenArr[j][2] -= box[0];
        datenArr[j][5] -= box[1];
    }
            
    // alte Werte zur Konvergenzüberprüfung
    double **alt = new double*[C];
    double e, s, m;
    
    for (j = 0; j < C; j++){
            alt[j] = new double[6];
    }
    
    bool fertig = false;    
    double tmpnenner = 0, tmp = 0, tmpzaehler = 0, point;
    int reply = 0;
    
    double newWeight[C + 1];
    
    
    double *thetaR = new double[RC];
    double *thetaT = new double[TC];
    double normalize[C];
    double t1, t2, t3, altw, rep;
    
    while (!fertig && reply++ < REPLY){
        // E-Schritt

        rj = 0;
        tj = 0;
        for (j = 0; j < C; j++){
            newWeight[j] = 0;
            for (r = 0; r < R; r++, rj++){
                thetaR[rj] = datenArr[j][6] ? datenArr[j][6] * ig(r, datenArr[j][0], datenArr[j][1], datenArr[j][2]) : 0;
                sumDrift[rj] = 0;
            }
            for (t = 0; t < T; ++t, tj++){
                thetaT[tj] = datenArr[j][6] ? ig(t, datenArr[j][3], datenArr[j][4], datenArr[j][5]) : 0;
                sumRet[tj] = 0;
            }  
        }
        newWeight[C] = 0;


        tmpWeight[0] = datenArr[C][6] / HGM;

        hgm = 0;
        for (r = 0; r < R; r++){
            for (t = 0; t < T; t++){
                
                for (j = 0, rj = r, tj = t, tmpnenner = tmpWeight[0]; j < C; j++, rj += R, tj += T){
                    tmpnenner += normalize[j] = thetaR[rj] * thetaT[tj];
                        
                }
                point = matrixArr[hgm++] / tmpnenner;
                
                for (j = 0, rj = r, tj = t; j < C; j++, rj += R, tj += T){
                    tmp = point * normalize[j];
                    sumDrift[rj] += tmp;
                    sumRet[tj] += tmp;
                }
            }
        }
        
        

        
        // Membership für HGM berechnen
        newWeight[C] = sumVals;
        int end = 0;
        for (j = 0; j < C; j++){
            tmp = 0;
            end += R;
            for (r = end - R; r < end; r++){
                tmp += sumDrift[r];
            }
            newWeight[C] -= newWeight[j] = tmp;
        }

            
        // Alte Werte abspeichern für Überprüfung
        // von Konvergenzverhalten
        for (j = 0; j < C; j++){
            getDescriptors(datenArr[j][0], datenArr[j][1], datenArr[j][2], &e, &s, &m);
            alt[j][0] = e;
            alt[j][1] = s;
            alt[j][2] = m;
            getDescriptors(datenArr[j][3], datenArr[j][4], datenArr[j][5], &e, &s, &m);
            alt[j][3] = e;
            alt[j][4] = s;
            alt[j][5] = m;
        }
            
            
        // M-Schritt
        // z berechnen
        for (j = 0; j < C + 1; j++){
            tmp = newWeight[j] / sumVals;
            if (isnan(tmp) || isinf(tmp) || tmp < 0 || tmp > 1 ){
                datenArr[j][6] = 0;
            }
            else
                datenArr[j][6] = tmp;
        }
            
        // Mue_r berechnen
        rj = 0; tj = 0;
        for (j = 0; j < C; j++){
            tmpzaehler = 0;
            for (tr = 0; tr < R; tr++){
                tmpzaehler += sumDrift[rj++] * tr;
            }
            tmp = (tmpzaehler / (datenArr[j][6] * sumVals)) - datenArr[j][2];            
            if (isnan(tmp) || isinf(tmp) || tmp < 0 ){
                datenArr[j][6] = 0;
            }
            else
                datenArr[j][0] = tmp;
        }
            
        // Mue_d berechnen
        for (j = 0; j < C; j++){
            tmpzaehler = 0;
            for(tt = 0; tt < T; tt++){
                tmpzaehler += sumRet[tj++] * tt;
            }
            tmp = (tmpzaehler / (datenArr[j][6] * sumVals)) - datenArr[j][5];
            if (isnan(tmp) || isinf(tmp) || tmp <= 0 ){
                datenArr[j][6] = 0;
            }
            else
                datenArr[j][3] = tmp;
        }
            
        // Lambda_r berechnen
        rj = 0; tj = 0;
        for (j = 0; j < C; j++){
            tmp = 1.0 / datenArr[j][0], tmpzaehler = 0;
            for (r = 0; r < R; r++){
                if ((r - datenArr[j][2]) != 0){
                    tmpzaehler += sumDrift[rj++] * (1.0 / (r - datenArr[j][2]) - tmp);
                }
            }
            tmp = (datenArr[j][6] * sumVals) / tmpzaehler;
            if (isnan(tmp) || isinf(tmp) || tmp <= 0 ){
                datenArr[j][6] = 0;
            }
            else
                datenArr[j][1] = tmp;
        }
            
        // Lambda_d berechnen
        for (j = 0; j < C; j++){
            tmp = (1.0 / datenArr[j][3]), tmpzaehler = 0;
            for(t = 0; t < T; t++){
                if ((t - datenArr[j][5]) != 0){
                    tmpzaehler += sumRet[tj++] * (1.0 / (t - datenArr[j][5]) - tmp);
                }
            }

            tmp = (datenArr[j][6] * sumVals) / tmpzaehler;
            if (isnan(tmp) || isinf(tmp) || tmp <= 0 ){
                datenArr[j][6] = 0;
            }
            else
                datenArr[j][4] = tmp;
        }
            
            
        // Offset_r berechnen
        rj = 0; tj = 0;
        for (j = 0; j < C; j++){
            tmp = datenArr[j][2], tmpzaehler = 0, tmpnenner = 0;
            altw = tmp + 1000;
            rep = 0;
            t1 = 3.0 / datenArr[j][1], t2 = 1.0 / square(datenArr[j][0]);
            while (!convergence(tmp, altw, CONVERGENCE_THRESHOLD) && rep++ < 15){
                rj = j * R;
                altw = tmp;
                tmpzaehler = 0, tmpnenner = 0;
                for (r = 0; r < R; r++){
                    t3 = r - tmp;
                    if (t3 != 0){
                        tmpzaehler += sumDrift[rj] * ((t1 / t3) + t2 - (1. / (t3 * t3)));
                        tmpnenner += sumDrift[rj++] * ((t1 / (t3 * t3)) - (1. / (t3 * t3 * t3)));
                    }
                }

                tmp = tmp - (tmpzaehler / tmpnenner);
                if (isnan(tmp) || isinf(tmp)){
                    datenArr[j][6] = 0;
                    break;
                }
            }
            datenArr[j][2] = tmp;
            
        }
            
            
        // Offset_d berechnen
        for (j = 0; j < C; j++){
            tmp = datenArr[j][5], tmpzaehler = 0, tmpnenner = 0;
            altw = tmp + 1000;
            rep = 0;
            t1 = 3.0 / datenArr[j][4], t2 = 1.0 / square(datenArr[j][3]);
            while (!convergence(tmp, altw, CONVERGENCE_THRESHOLD) && rep++ < 15){
                altw = tmp;
                tj = j * T;
                tmpzaehler = 0, tmpnenner = 0;
                for(t = 0; t < T; t++){
                    t3 = t - tmp;
                    if (t3 != 0){
                        tmpzaehler += sumRet[tj] * ((t1 * t3) + t2 - (1. / (t3 * t3)));
                        tmpnenner += sumRet[tj++] * ((t1 * t3 * t3) - (1. / (t3 * t3 * t3)));
                    }
                }
                
                tmp = tmp - (tmpzaehler / tmpnenner);
                if (isnan(tmp) || isinf(tmp)){
                    datenArr[j][6] = 0;
                    break;
                }
            }
            datenArr[j][5] = tmp;
        }

        fertig = true;
        // Berechnung, ob Funktion konvergiert
        for (j = 0; j < C; j++){
            getDescriptors(datenArr[j][0], datenArr[j][1], datenArr[j][2], &e, &s, &m);
            if (!convergence(e, alt[j][0], CONVERGENCE_THRESHOLD)){
                fertig = false;
                break;
            }
            if (!convergence(s, alt[j][1], CONVERGENCE_THRESHOLD)){
                fertig = false;
                break;
            }
            if (!convergence(m, alt[j][2], CONVERGENCE_THRESHOLD)){
                fertig = false;
                break;
            }                       
            
            getDescriptors(datenArr[j][3], datenArr[j][4], datenArr[j][5], &e, &s, &m);
            if (!convergence(e, alt[j][3], CONVERGENCE_THRESHOLD)){
                fertig = false;
                break;
            }
            if (!convergence(s, alt[j][4], CONVERGENCE_THRESHOLD)){
                fertig = false;
                break;
            }
            if (!convergence(m, alt[j][5], CONVERGENCE_THRESHOLD)){
                fertig = false;
                break;
            }

        }
    }
    
    for (j = 0; j < C; j++){
        datenArr[j][2] += box[0];
        datenArr[j][5] += box[1];
    }
    datenArr[C][0] = sumVals;
    
    for (j = 0; j < C; j++){
        delete []alt[j];
    }
    delete []alt;

    delete []thetaT;
    delete []thetaR;

    delete []sumRet;
    delete []sumDrift;
    delete []matrixArr;

}


/*!
 * \brief Modeling peaks with statistical distributions using EM algorithm
 * @param measurement Instance of an IMSMeasurement
 * @param peak_list Peaklist (with final peaks)
 * @return Peaklist (with modeled peaks and additional parameters)
 */
modeling(pme_modeling){
    // First phase for accelereration is a segmentation of the measurement in indipendent
    // peak areas. Boxes are build, enclosing every single peak. When two boxes opverlap
    // they'll be joined
    
    
    assert(pipeline_parameters->find("expansion_size") != pipeline_parameters->end());    
    double expansion_size = atof(pipeline_parameters->at("expansion_size").c_str());
    
    assert(pipeline_parameters->find("intensity_threshold") != pipeline_parameters->end());    
    double min_peak_signal = atof(pipeline_parameters->at("intensity_threshold").c_str());
    
    assert(pipeline_parameters->find("tol_rt") != pipeline_parameters->end());    
    double tol_rt = atof(pipeline_parameters->at("tol_rt").c_str());
    
    assert(pipeline_parameters->find("tol_rt_percent") != pipeline_parameters->end());    
    double tol_rt_procent = atof(pipeline_parameters->at("tol_rt_percent").c_str());
    

    // double[0] retention-value of lower left rect point (index not time)
    // double[1] rim-value of lower left rect point (index not time)
    // double[2] retention-value of upper right rect point (index not time)
    // double[3] rim-value of upper right rect point (index not time)
    vector< int* >* rects = new vector< int* >;
    int len_list = peak_list->ims_peak_list->size();
    
    if (!len_list) return new IMSPeakList(peak_list);
    
    int R = measurement->data->n_rows;
    int T = measurement->data->n_cols;
    
    assert(R >= expansion_size);
    assert(T >= expansion_size);

    int** mx = new int*[R];
    for_loop (r, R){
        mx[r] = new int[T];
        for_loop (t, T){
            mx[r][t] = -1;
        }
    }
    
    deque< int* >* queue = new deque< int* >;
    vector< vector< int >* >* indizes = new vector< vector< int >* >;
    
    int r, t, index = -1, tm, tp, rm, rp;
    for_loop (d, len_list){
        
        if (mx[peak_list->ims_peak_list->at(d)->index_r][peak_list->ims_peak_list->at(d)->index_t] != -1){
            indizes->at(mx[peak_list->ims_peak_list->at(d)->index_r][peak_list->ims_peak_list->at(d)->index_t])->push_back(d);
            continue;
        }
        else {
            index++;
            indizes->push_back(new vector<int>);
            indizes->back()->push_back(d);
        }
        
        
        queue->push_back(new int[2]);
        queue->back()[0] = peak_list->ims_peak_list->at(d)->index_r;
        queue->back()[1] = peak_list->ims_peak_list->at(d)->index_t;
        mx[queue->back()[0]][queue->back()[1]] = index;
        rm = queue->back()[0];
        rp = queue->back()[0];
        tm = queue->back()[1];
        tp = queue->back()[1];
        
        while (!queue->empty()){
            
            r = queue->front()[0];
            t = queue->front()[1];
            delete [](queue->front());
            queue->pop_front();
            
            for (int dr = r - expansion_size; dr <= r + expansion_size; dr++){
                for (int dt = t - expansion_size; dt <= t + expansion_size; dt++){
                        
                    if (dr >= 0 && dt >= 0 && dr < R && dt < T){
                        if (mx[dr][dt] == -1){
                            if (measurement->data->data->at(dr)->at(dt) >= MINHEIGHT_SEGMENTATION){
                                queue->push_back(new int[2]);
                                queue->back()[0] = dr;
                                queue->back()[1] = dt;
                            }
                            mx[dr][dt] = index;
                            if (rm > dr) rm = dr;
                            if (rp < dr) rp = dr;
                            if (tm > dt) tm = dt;
                            if (tp < dt) tp = dt;
                        }
                    }
                }
            }
        }
        int* tmp_rect = new int[4];
        tmp_rect[0] = rm;
        tmp_rect[1] = tm;
        tmp_rect[2] = rp;
        tmp_rect[3] = tp;
        rects->push_back(tmp_rect);
    }
    
    for_loop (r, R){
        delete []mx[r];
    }
    delete []mx;
    delete queue;
   
    
    for (int i = 0; i < rects->size() - 1; i++){
        for (int j = i + 1; j < rects->size(); j++){
            if (intersection(rects->at(i), rects->at(j))){
                
                rects->at(i)[0] = min(rects->at(i)[0], rects->at(j)[0]);
                rects->at(i)[1] = min(rects->at(i)[1], rects->at(j)[1]);
                rects->at(i)[2] = max(rects->at(i)[2], rects->at(j)[2]);
                rects->at(i)[3] = max(rects->at(i)[3], rects->at(j)[3]);
                
                delete []rects->at(j);
                rects->erase(rects->begin() + j);
                
                for (int k = 0; k < indizes->at(j)->size(); k++){
                    indizes->at(i)->push_back(indizes->at(j)->at(k));
                }
                delete indizes->at(j);
                indizes->erase(indizes->begin() + j);
                
                
                assert(rects->at(i)[2] - rects->at(i)[0] + 1 > expansion_size);
                assert(rects->at(i)[3] - rects->at(i)[1] + 1 > expansion_size);
                
                i = -1;
                break;
            }
        }
    }
    
        
        
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Rects and indexes for segmented EM algorithm computed
    
    // Copy existing 
    IMSPeakList* new_peak_list = new IMSPeakList(peak_list);
    new_peak_list->parameter_names->push_back("mue_t");
    new_peak_list->parameter_names->push_back("sigma_t");
    new_peak_list->parameter_names->push_back("mue_r");
    new_peak_list->parameter_names->push_back("sigma_r");
    
    
    if (rects->size()){
        //#pragma omp parallel for
        for_loop (i, rects->size()){
            // Creating parameter table (+1 for background model)
            double** parameters = new double*[indizes->at(i)->size() + 1];
            double distrib = 1. / (indizes->at(i)->size() + 1.);
            for_loop (j, indizes->at(i)->size()){
                uint peak_index = indizes->at(i)->at(j);
                parameters[j] = new double[7];
                
                double r, t, mue, lambda, offset;
                r = peak_list->ims_peak_list->at(peak_index)->index_r;
                t = peak_list->ims_peak_list->at(peak_index)->index_t;
                //double sig = min(2., (r * tol_rt_procent + tol_rt) / 3.);
                getModelParams(r + 0.001, 5, r, &mue, &lambda, &offset);
                parameters[j][0] = mue;
                parameters[j][1] = lambda;
                parameters[j][2] = offset;
                
                getModelParams(t + 0.001, 5, t, &mue, &lambda, &offset);
                parameters[j][3] = mue;
                parameters[j][4] = lambda;
                parameters[j][5] = offset;
                parameters[j][6] = distrib;
                //cout << peak_list->ims_peak_list->at(peak_index)->r << " " << peak_list->ims_peak_list->at(peak_index)->t << endl;
                //cout << r << " " << t << endl;
            }
            // background model
            parameters[indizes->at(i)->size()] = new double[7];
            parameters[indizes->at(i)->size()][0] = 1;
            parameters[indizes->at(i)->size()][1] = 1;
            parameters[indizes->at(i)->size()][2] = 1;
            parameters[indizes->at(i)->size()][3] = 1;
            parameters[indizes->at(i)->size()][4] = 1;
            parameters[indizes->at(i)->size()][5] = 1;
            parameters[indizes->at(i)->size()][6] = distrib;
            
            em_estimate_2D(measurement->data, indizes->at(i)->size(), parameters, rects->at(i), min_peak_signal);
            
            #pragma omp critical(dataupdate)
            for_loop (j, indizes->at(i)->size()){
                uint peak_index = indizes->at(i)->at(j);
                if (parameters[j][6]){
                    double mue_t, sigma_t, mode_t, mue_r, sigma_r, mode_r;
                    getDescriptors(parameters[j][3], parameters[j][4], parameters[j][5], &mue_t, &sigma_t, &mode_t);
                    getDescriptors(parameters[j][0], parameters[j][1], parameters[j][2], &mue_r, &sigma_r, &mode_r);

                
                    new_peak_list->ims_peak_list->at(peak_index)->r = measurement->computeRetention(mode_r);
                    new_peak_list->ims_peak_list->at(peak_index)->t = measurement->computeInverseMobility(mode_t);
                    
                    new_peak_list->ims_peak_list->at(peak_index)->signal = parameters[indizes->at(i)->size()][0] * ig(mode_t, parameters[j][3], parameters[j][4], parameters[j][5]) * ig(mode_r, parameters[j][0], parameters[j][1], parameters[j][2]);
                    new_peak_list->ims_peak_list->at(peak_index)->volume = parameters[j][6] * parameters[indizes->at(i)->size()][0];

                    new_peak_list->ims_peak_list->at(peak_index)->input_parameter("mue_t", measurement->computeInverseMobility(mue_t, true));
                    new_peak_list->ims_peak_list->at(peak_index)->input_parameter("sigma_t", measurement->computeInverseMobility(sigma_t, true));
                    new_peak_list->ims_peak_list->at(peak_index)->input_parameter("mue_r", measurement->computeRetention(mue_r));
                    new_peak_list->ims_peak_list->at(peak_index)->input_parameter("sigma_r", measurement->computeRetention(sigma_r));
                }
                else {
                    new_peak_list->ims_peak_list->at(peak_index)->input_parameter("mue_t", -1);
                    new_peak_list->ims_peak_list->at(peak_index)->input_parameter("sigma_t", -1);
                    new_peak_list->ims_peak_list->at(peak_index)->input_parameter("mue_r", -1);
                    new_peak_list->ims_peak_list->at(peak_index)->input_parameter("sigma_r", -1);
                }
                delete []parameters[j];
            }
            delete []parameters[indizes->at(i)->size()];
            delete []parameters;
        }

    }

    
    for_loop (i, rects->size()){
        delete []rects->at(i);
    }
    delete rects;

    for_loop (i, indizes->size()){
        delete indizes->at(i);
    }
    delete indizes;
    
    return new_peak_list;
}

/*!
 * \brief Dummy empty function that only produces a hard copy of the input peaklist
 * @param measurement Instance of an IMSMeasurement
 * @param ims_peak_list Peaklist (with final peaks)
 * @return Peaklist (with modeled peaks and no additional parameters)
 */
modeling(no_modeling){
    return new IMSPeakList(peak_list);
}
