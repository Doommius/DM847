#include "picking.h"

bool compare_signals_IMSPeaks (IMSPeak* peak1, IMSPeak* peak2){
	double signal1 = peak1->signal; 
	double signal2 = peak2->signal;
	return (signal1 >= signal2);
}

/*!
 * \brief Checks whether two peaks are too close 
 * @param current_peak peak that is checked
 * @param maxp considered peak
 * @param tRT distance for the retention time
 * @param tolDT ditance for the inverse reduced mobility
 * @return true if current_peak is too close to maxp
 */
bool too_near(IMSPeak* current_peak, IMSPeak* maxp, double tRT, double tolDT) {
	double diff_DT = fabs(maxp->t - current_peak->t);
	double diff_RT = fabs(maxp->r - current_peak->r);
	return ( (diff_DT < tolDT) && ( diff_RT < tRT) ); 
}


/*!
 * \brief Calculates the score between two peaks. Note that a low distance between the peaks yields a high similarity.
 * @param pointOne first peak
 * @param pointTwo second peak
 * @param pipeline_parameters map of parameters (map<Key = string, Value = string>)
 * @return similarity score between the two peaks
 */
double sim_score(IMSPeak* pointOne, IMSPeak* pointTwo, pmap* pipeline_parameters){

    assert(pipeline_parameters->find("tol_rt_percent") != pipeline_parameters->end());
    assert(pipeline_parameters->find("tol_rt") != pipeline_parameters->end());
    assert(pipeline_parameters->find("tol_rim") != pipeline_parameters->end());
	assert(pipeline_parameters->find("ce_weight_exponent") != pipeline_parameters->end());

	int ce_weight_exponent = atoi( pipeline_parameters->at("ce_weight_exponent").c_str() );
    double tol_rt_percent = atof(pipeline_parameters->at("tol_rt_percent").c_str());
    double tol_rt = atof(pipeline_parameters->at("tol_rt").c_str());

    double tolRT = 0.0;

    if(pointOne->r >= pointTwo->r){
        tolRT = (tol_rt_percent * pointOne->r  + tol_rt) * 1.0;
    }else{
        tolRT = (tol_rt_percent * pointTwo->r  + tol_rt) * 1.0;
    }

    double tolT = atof(pipeline_parameters->at("tol_rim").c_str());
    double deltaX = pointOne->t - pointTwo->t;
    double deltaY = pointOne->r - pointTwo->r;

    double distance = ( pow(deltaX/tolT,2) + pow(deltaY/tolRT, 2) ) / 2.0;

	if (distance <= 1.0){
		return (double)( pow(2,ce_weight_exponent - ce_weight_exponent * distance ) - 1.0);
	} else{
		return (double) ( - (distance - 1.0) );
    }
}

/*!
 * \brief Merges candidate peaks into a final peak list, also called merging by intensity
 * @param peak_list ims_peak_list (with candidate peaks)
 * @param pipeline_parameters map of parameters (map<Key = string, Value = string>)
 * @return Peaklist (with final peaks)
 */
picking(merging_by_signal){
	
    assert(pipeline_parameters->find("tol_rt_percent") != pipeline_parameters->end());
    assert(pipeline_parameters->find("tol_rt") != pipeline_parameters->end());
    assert(pipeline_parameters->find("tol_rim") != pipeline_parameters->end());
    
    double tolDT = atof(pipeline_parameters->at("tol_rim").c_str());
    double tolRT = atof(pipeline_parameters->at("tol_rt").c_str());
    double tolRTproz = atof(pipeline_parameters->at("tol_rt_percent").c_str());
    
    //sort the incoming peaks by intensity/signal descending
    IMSPeakList* tmp_list = new IMSPeakList(peak_list);
    vector< IMSPeak* >* sorted_vector = tmp_list->ims_peak_list;
    sort(sorted_vector->begin(), sorted_vector->end(), compare_signals_IMSPeaks);
	
    //////////////Merging process////////////////////////////////////////////////////////////////////////////// 
    for ( size_t i = 0; i < sorted_vector->size(); i++ ){
        IMSPeak* maxp = sorted_vector->at(i); // peak with the highest intensity/signal
        double tRT = tolRT + maxp->r * tolRTproz;
		for ( size_t j = i + 1; j < sorted_vector->size(); j++ ) {            
            if ( too_near( sorted_vector->at(j), maxp, tRT, tolDT) ) {
                //delete the smaller one, because of the descending order that's always j
                sorted_vector->erase(sorted_vector->begin() + j);
                j--;
            }
        }
    }
    
    ////////////////// Constructing the IMSPeaklist/////////////////////
    IMSPeakList* picked_peaklist = new IMSPeakList;
    picked_peaklist->measurement_source = peak_list->measurement_source;
    
	for(size_t i = 0; i < sorted_vector->size(); i++){
        IMSPeak* maxp = sorted_vector->at(i);
        picked_peaklist->ims_peak_list->push_back(new IMSPeak(maxp));        
    }
    
    delete tmp_list;
    return picked_peaklist;
}


/*!
 * \brief Merges candidate peaks into a final peak list by using cluster editing
 * @param peak_list ims_peak_list (with candidate peaks)
 * @param pipeline_parameters map of parameters (map<Key = string, Value = string>)
 * @return Peaklist (with final peaks)
 */
picking(cluster_editing){
	int incoming_peaks =  peak_list->ims_peak_list->size();
	assert(incoming_peaks != 0);
	
	map< string, IMSPeak* >* name_peak_map = new map< string, IMSPeak* >;
    for ( int i = 0; i < peak_list->ims_peak_list->size(); i++ ){
        name_peak_map->insert(pair<string, IMSPeak*>(peak_list->ims_peak_list->at(i)->peak_name.c_str(), peak_list->ims_peak_list->at(i)));            
    }
    
    /* Calculate the score_matrix for cluster editing
    Matrix for peaks p1, p2, p3, p4 with score of px,py as s(x,y)
    s(1,2) s(1,3) s(1,4)
    s(2,3) s(2,4)
    s(3,4)    
    */
    matrix_t* score_matrix = new matrix_t;
    score_matrix->resize(peak_list->ims_peak_list->size() - 1);
    #pragma omp parallel for
    for ( int i = 0; i < peak_list->ims_peak_list->size() - 1; i++ ){
        vector< double >* row = new vector< double >;
        for ( int j = i + 1; j < peak_list->ims_peak_list->size(); j++ ){
            double score = sim_score(peak_list->ims_peak_list->at(i),peak_list->ims_peak_list->at(j), pipeline_parameters);
            row->push_back(score);
        }
        score_matrix->at(i) = row;
    }
    
    /////Clustering method used here is yoshiko/////////////////////
    int z = system("rm -r cluster");
	int y = system("mkdir cluster");

    ofstream file;
	file.open("./cluster/score_matrix.cm");
	file << peak_list->ims_peak_list->size() << endl;
	for ( int i = 0; i < peak_list->ims_peak_list->size(); i++ ){
	    file << peak_list->ims_peak_list->at(i)->peak_name << endl;
	}

	for ( int i = 0; i < score_matrix->size(); i++ ){
	    for ( int j = 0; j < score_matrix->at(i)->size(); j++ ){
	        file << score_matrix->at(i)->at(j) << "\t";
	    }
	    file << endl;
	}
    
	vector< vector< IMSPeak* >* >* clusters = new vector< vector< IMSPeak* >* >;
    
	//Reading the Clusters from files generated by yoshiko.
	int c = system("./yoshiko -F 0 -O 0 -r 111111 -f cluster/score_matrix.cm -o cluster/cluster_results");
	if( c != 0){
		cout << "Please provide the software yoshiko. The executable of yoshiko should be in the same folder as the executable of peax. Visit http://www.cwi.nl/research/planet-lisa to download the yoshiko software." << endl;
		exit(0);
	}
	fstream cluster_file("./cluster/cluster_results_0.csv", ios::in);

	string line;
	while ( !cluster_file.eof() ){
		vector< IMSPeak* >* a_cluster = new vector< IMSPeak* >;
		getline (cluster_file,line);
		stringlist peaks = split(line, '\t');
		for (int c = 0; c < peaks.size(); ++c){
	    	a_cluster->push_back(name_peak_map->at(peaks.at(c)));
	    }
		if(a_cluster->size() > 0 ){
	        clusters->push_back(a_cluster);
	    }
	}
	cluster_file.close();
    

    ////Choose an representant: Peak in the cluster with highest intensity////
    IMSPeakList* clustered_list = new IMSPeakList;
    clustered_list->measurement_source = peak_list->measurement_source;
    for (int i = 0; i < clusters->size(); i++){
        double intensity = clusters->at(i)->at(0)->signal;
        IMSPeak* peak = clusters->at(i)->at(0);
        for (int j = 1; j < clusters->at(i)->size(); j++){
            if (intensity < clusters->at(i)->at(j)->signal){
                intensity = clusters->at(i)->at(j)->signal;
                peak = clusters->at(i)->at(j);
            }
        }
        
		clustered_list->ims_peak_list->push_back(new IMSPeak(peak)); 
    }

    //////////////////Deleting pointer//////////////////////////////////////////////// 
    for_loop(i,score_matrix->size()){
        delete score_matrix->at(i);
    }
    delete score_matrix;
    
    delete name_peak_map;
	
    for_loop(i, clusters->size()){
        delete clusters->at(i);
    }
    delete clusters;
    //////////////////////////////////////////////////////////////////////////////////
	
	int outgoing_peaks = clustered_list->ims_peak_list->size();
	assert(incoming_peaks >= outgoing_peaks);
	
	
    return clustered_list;
}

/*!
 * \brief Merges candidate peaks into a final peak list by using a modified EM algorithm
 * @param peak_list ims_peak_list (with candidate peaks)
 * @param pipeline_parameters map of parameters (map<Key = string, Value = string>)
 * @return Peaklist (with final peaks)
 */
picking(em_clustering){
	int incoming_peaks =  peak_list->ims_peak_list->size();
	assert(incoming_peaks != 0);
	
    assert(pipeline_parameters->find("tol_rt_percent") != pipeline_parameters->end());
    assert(pipeline_parameters->find("tol_rt") != pipeline_parameters->end());
    assert(pipeline_parameters->find("tol_rim") != pipeline_parameters->end());
    
    double tol_rim = atof(pipeline_parameters->at("tol_rim").c_str());
    double tol_rt_percent = atof(pipeline_parameters->at("tol_rt_percent").c_str());
    double tol_rt = atof(pipeline_parameters->at("tol_rt").c_str());
    
    
    
    bool multiple_merging = false;
    if (pipeline_parameters->find("multiple_merging") != pipeline_parameters->end()){
        multiple_merging = atoi(pipeline_parameters->at("multiple_merging").c_str());
    }
    uint len_field = peak_list->ims_peak_list->size(); 
    
    vector<double*>* data = new vector<double*>;
    vector< vector < int >* >* cluster = 0;
    
    if (multiple_merging){
        cluster = new vector< vector < int >* >;
    }
    
    
    double equal_dist = 1. / double(len_field);
    for_loop (i, len_field){
        double* tmp = new double[6];
        tmp[0] = peak_list->ims_peak_list->at(i)->t;
        tmp[1] = peak_list->ims_peak_list->at(i)->r;
        tmp[2] = tol_rim / 2.;
        tmp[3] = (tol_rt_percent * tmp[1] + tol_rt) / 2.;
        tmp[4] = equal_dist;
        tmp[5] = i;
        data->push_back(tmp);
        
        if (multiple_merging){
            vector <int>* tmp_cluster = new vector <int>;
            tmp_cluster->push_back(i);
            cluster->push_back(tmp_cluster);
        }
    }
    
    // Hidden variables describing normalized membership of data point to model
    // Don't need to be initialized reasonable, because in first expectation step
    // immediately computed
    vector<double*>* member = new vector<double*>;
    for_loop (i, len_field){
        double* tmp = new double[len_field];
        member->push_back(tmp);
    }
    
    
    for_loop (a, 20){
        //E Step
        #pragma omp parallel for
        for_loop (i, len_field){
            double* p = new double[data->size()];
            double n = 1e-20;
            
            
            for_loop (j, data->size()){
                double rr = g(peak_list->ims_peak_list->at(i)->t, data->at(j)[0], data->at(j)[2]);
                n += p[j] = data->at(j)[4] * rr * g(peak_list->ims_peak_list->at(i)->r, data->at(j)[1], data->at(j)[3]);
            }
                
            for_loop (j, data->size()){
                #pragma omp critical(dataupdate)
                member->at(j)[i] = p[j] / n;
            }
            delete []p;
        } 
        //M Step
        //estimate weight
        #pragma omp parallel for
        for_loop (i, member->size()){
            double summe = 0;
            for_loop (j, len_field){
                summe += member->at(i)[j];
            }
            #pragma omp critical(dataupdate)
            data->at(i)[4] = summe / double(len_field);
        }
        
        //estimate mu_x, mu_y
        #pragma omp parallel for
        for_loop (j, data->size()){
            double mean_t = 0, mean_r = 0;
            double n = 0;
            for_loop (i, len_field){
                mean_t += member->at(j)[i] * peak_list->ims_peak_list->at(i)->t;
                mean_r += member->at(j)[i] * peak_list->ims_peak_list->at(i)->r;
                n += member->at(j)[i];
            }
            #pragma omp critical(dataupdate)
            {
                data->at(j)[0] = mean_t / n;
                data->at(j)[1] = mean_r / n;
            }
        }
            
        //estimate sigma_x, sigma_y
        #pragma omp parallel for
        for_loop (j, data->size()){
            double sta_dev_t = 0, sta_dev_r = 0;
            double n = 0;
            for_loop (i, len_field){
                sta_dev_t += member->at(j)[i] * square(peak_list->ims_peak_list->at(i)->t - data->at(j)[0]);
                sta_dev_r += member->at(j)[i] * square(peak_list->ims_peak_list->at(i)->r - data->at(j)[1]);
                n += member->at(j)[i];
            }
            sta_dev_t /= n;
            sta_dev_r /= n;
            
            
            
            #pragma omp critical(dataupdate)
            {
                data->at(j)[2] = (sta_dev_t > 1e-10) ? sqrt(sta_dev_t) : 1e-5;
                data->at(j)[3] = (sta_dev_r > 1e-8) ? sqrt(sta_dev_r) : 1e-4;
            }
        }
        
        // check all peaks against all peaks, if two peaks are too close to each other, the smaller one
        // will be deleted and weight of smaller (deleted) peak is added to remaining (merged) peak
        for (int i = 0; i < data->size() - 1; i++){
            for (int j = data->size() - 1; j > i; j--){
                if (fabs(data->at(i)[0] - data->at(j)[0]) < atof(pipeline_parameters->at("tol_rim").c_str()) && fabs(data->at(i)[1] - data->at(j)[1]) < atof(pipeline_parameters->at("tol_rt_percent").c_str()) * data->at(i)[1] + atof(pipeline_parameters->at("tol_rt").c_str())){
                    
                    if (peak_list->ims_peak_list->at(data->at(i)[5])->signal < peak_list->ims_peak_list->at(data->at(j)[5])->signal){
                        swap(data->at(i), data->at(j));
                        swap(member->at(i), member->at(j));
                        if (multiple_merging){
                            swap(cluster->at(i), cluster->at(j));
                            }
                    }
                    if (multiple_merging){
                        for_loop(m, cluster->at(j)->size()){
                            cluster->at(i)->push_back(cluster->at(j)->at(m));
                        }
                        delete cluster->at(j);
                        cluster->erase(cluster->begin() + j);
                    }
                    
                    data->at(i)[4] += data->at(j)[4];
                    delete []data->at(j);
                    data->erase(data->begin() + j);
                    for_loop (r, len_field){
                        member->at(i)[r] += member->at(j)[r];
                    }
                    delete []member->at(j);
                    member->erase(member->begin() + j);
                }
            }
        }
    }
    
    for_loop (i, member->size()){
        delete []member->at(i);
    }
    delete member;
    
    
    // Copy original peak list to modify entries
    IMSPeakList* ipl = 0;
    
    if (!multiple_merging){
        ipl = new IMSPeakList(peak_list);
        
        // set is used to check whick peak will be discarded
        set<int> invalid_peaks;
        for_loop (i, ipl->ims_peak_list->size()){
            invalid_peaks.insert(i);
        }
        
        
        
        // modifying peaks in new list and deleting indexes from set that determine that
        // modified peak shall stay in new list
        for_loop (i, data->size()){
            invalid_peaks.erase(data->at(i)[5]);
            // taking mean of models as position of peaks
            /*ipl->ims_peak_list->at(data->at(i)[5])->r = data->at(i)[1];
            ipl->ims_peak_list->at(data->at(i)[5])->t = data->at(i)[0];
            if (ipl->measurement_source){
                ipl->ims_peak_list->at(data->at(i)[5])->index_r = ipl->measurement_source->getRetentionIndex(data->at(i)[1]);
                ipl->ims_peak_list->at(data->at(i)[5])->index_t = ipl->measurement_source->getRIMIndex(data->at(i)[0]);
            }*/
            delete []data->at(i);
        }
        delete data;
        
        // invalid peaks will be deleted in new list by considering all remaining indexes in set
        for (set<int>::iterator it = invalid_peaks.begin(); it != invalid_peaks.end(); it++){
            delete ipl->ims_peak_list->at(*it);
            ipl->ims_peak_list->at(*it) = 0;
        }
        // all entries in new list where enrty is NULL are deleted
        for (int i = ipl->ims_peak_list->size() - 1; i >= 0; i--){
            if (!ipl->ims_peak_list->at(i)){
                ipl->ims_peak_list->erase(ipl->ims_peak_list->begin() + i);
            }
        }
    }
    else {
        assert(pipeline_parameters->find("peak_per_cluster") != pipeline_parameters->end());
        int min_peak_pro_cluster = atoi(pipeline_parameters->at("peak_per_cluster").c_str());
        ipl = new IMSPeakList();
        
        ipl->parameter_names->push_back("cluster");
        ipl->parameter_names->push_back("t_tol");
        ipl->parameter_names->push_back("r_tol");
        for_loop(i, cluster->size()){
            
                
            for_loop(j, cluster->at(i)->size()){
                peak_list->ims_peak_list->at(cluster->at(i)->at(j))->input_parameter("cluster", i);
                peak_list->ims_peak_list->at(cluster->at(i)->at(j))->input_parameter("r_tol", -1);
                peak_list->ims_peak_list->at(cluster->at(i)->at(j))->input_parameter("t_tol", -1);
            }
            
            IMSPeak* peak = new IMSPeak;
            peak->measurement_name = "consensus_layer";
            peak->peak_name = "C" + i2s(i + 1);
            peak->r = data->at(i)[1];
            peak->t = data->at(i)[0];
            peak->signal = 0;
            peak->volume = 0;
            if (cluster->at(i)->size() >= min_peak_pro_cluster){
                peak->input_parameter("cluster", -1);
            }
            else {
                peak->input_parameter("cluster", -2);
            }
            
            double rt = tol_rt_percent * peak->r + tol_rt;
            peak->input_parameter("r_tol", (data->at(i)[3] < rt ? rt : data->at(i)[3]));
            peak->input_parameter("t_tol", (data->at(i)[2] <  tol_rim ?  tol_rim : data->at(i)[2]));
            
            ipl->ims_peak_list->push_back(peak);
            delete cluster->at(i);
        }
        delete cluster;
        
        //delete null peaks
        for (int i = ipl->ims_peak_list->size() - 1; i >= 0; i--){
            if (!ipl->ims_peak_list->at(i)){
                ipl->ims_peak_list->erase(ipl->ims_peak_list->begin() + i);
            }
        }
    }
    int outgoing_peaks = ipl->ims_peak_list->size();
    if (!multiple_merging){
        assert(incoming_peaks >= outgoing_peaks);
    }
    return ipl;
}

/*!
 * \brief Dummy empty function that only produces a hard copy of the input peaklist
 * @param ims_peak_list Peaklist (with candidate peaks)
 * @return Peaklist (with final peaks)
 */
picking(no_picking){
    return new IMSPeakList(peak_list);
}
