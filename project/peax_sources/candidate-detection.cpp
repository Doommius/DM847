#include "candidate-detection.h"

int area_size;
int th;

/*!
 * \brief Finds peak areas.
 *  A peak area is a set of points (t,r) that are contigiuous in the given matrix.
 *  A peak area contains at least area_size points. 
 * @param possiblePoints all points with signal intensity > threshold
 * @param matrix entire matrix
 * @return A vector of peak areas as vectors of pairs
 */
vector< vector< pair<int,int>* >* >* find_peak_areas(vector< pair<int,int>* >* possiblePoints, matrix_t* matrix){
   
    vector< vector< pair<int,int>* >* >* peak_areas = new vector< vector< pair<int,int>* >* >;
    while (!possiblePoints->empty()){
        pair<int,int>* start = possiblePoints->back();
		possiblePoints->pop_back();
        vector< pair<int,int>* >* a_peak_area = new vector< pair<int,int>* >;
        a_peak_area->push_back(start);         
        
        int a = a_peak_area->size();    
        for(int i = 0; i < a; i++){
            int x = a_peak_area->at(i)->first;
            int y = a_peak_area->at(i)->second;
            
            int p = possiblePoints->size();
            for(int j= 0; j < p; j++){
                //cout << "Area size " << a_peak_area->size() << endl;
                //cout << "Point " << x << " " << y << " has neighbor: " << endl;
                //cout << "Pointer at pos " << j<< endl;
                int x_new = possiblePoints->at(j)->first;
                int y_new = possiblePoints->at(j)->second;
                //cout << "Point1 " << x << " " << y << endl;
                //cout << "Point2 " << x_new << " " << y_new << endl;
                //bool x_ab = abs(x_new - x) <= 1;
                //bool y_ab = abs(y_new - y) <= 1;
                //cout << "x abstand " << x_ab << "y abstand " << y_ab << endl;
                if(abs(x_new - x) <= 1 && abs(y_new - y) <= 1){
                    //cout << "Point1 " << x << " " << y << endl;
                    //cout << "Point " << x_new << " " << y_new << endl;
                    a_peak_area->push_back(possiblePoints->at(j));
                    //a_peak_area->push_back(new pair<int,int>(x_new, y_new));
                    //delete possiblePoints->at(j);
                    possiblePoints->erase(possiblePoints->begin() + j);
                    --j;
                    ++a;
                    --p;
                }
                
            }
        }
        //only areas containing min. area_size points are considered peak candidates
        if((signed int)a_peak_area->size() >= area_size){
        
            peak_areas->push_back(a_peak_area);
        
        }else{
            ///////Deleting pointers of unused area///////////
            for_loop(i, a_peak_area->size()){
                delete a_peak_area->at(i);
            }
            delete a_peak_area;
            //////////////////////////////////////////////////
        }  
    }
    
    return peak_areas;
}

/*!
* \brief Checks if all 8 neighbor are in the area.
* @param label identifes the considered area
* @param area considered area
* @param x coordinate of considered points
* @param y coordinate of considered points
* @param area_matrix labeled matrix
* return True if all 8 neighbor have intensity > threshold
*/
bool points_in_area(int label, vector< pair<int,int>* >* area, int x, int y, vector < vector < short >* >* area_matrix){
    pair<int,int>* one = new pair<int,int>(x-1, y);
    pair<int,int>* two = new pair<int,int>(x+1, y);
    pair<int,int>* three = new pair<int,int>(x, y-1); 
    pair<int,int>* four = new pair<int,int>(x, y+1); 
    pair<int,int>* five = new pair<int,int>(x+1, y+1);
    pair<int,int>* six = new pair<int,int>(x-1, y-1);
    pair<int,int>* seven = new pair<int,int>(x+1, y-1);
    pair<int,int>* eight = new pair<int,int>(x-1, y+1);
    

    vector< pair<int,int>* >* neighbors = new vector< pair<int,int>* >;
    neighbors->push_back(one);
    neighbors->push_back(two);
    neighbors->push_back(three);
    neighbors->push_back(four);
    neighbors->push_back(five);
    neighbors->push_back(six);
    neighbors->push_back(seven);
    neighbors->push_back(eight);
    
    
    for_loop(i, neighbors->size()){
        if(area_matrix->at(neighbors->at(i)->second)->at(neighbors->at(i)->first) != label){
			for_loop(n, neighbors->size()){
				delete neighbors->at(n); 
			}
			delete neighbors;
            return false;
        }
    }
    
    for_loop(n, neighbors->size()){
	    delete neighbors->at(n); 
	}
    delete neighbors;
    return true;    
}

/*!
 * \brief Finds peak candidates by considering neighborhood connectivity
 * @param measurement Instance of an IMSMeasurement
 * @param pipeline_parameters map of parameters (map<Key = string, Value = string>) 
 * @return IMSPeakList of IMSPeaks
 */
candidate_detection(local_maxima){
    assert(pipeline_parameters->find("area_size") != pipeline_parameters->end());
    assert(pipeline_parameters->find("intensity_threshold") != pipeline_parameters->end());
    
    area_size = atoi(pipeline_parameters->at("area_size").c_str());
    th = atoi(pipeline_parameters->at("intensity_threshold").c_str());
    
    size_t sizeT =  measurement->data->n_cols;
    size_t sizeR =  measurement->data->n_rows;  

    //////Generating the possible points, all points have intensity > threshold///////
    vector< pair<int,int>* >* possiblePoints= new vector< pair<int,int>* >;
	double value = 0.0;
	for(size_t i = 1; i < sizeT-1; i++){ //excluding the margin
		for(size_t j = 1; j < sizeR-1; j++){ //excluding the margin
			value = measurement->data->data->at(j)->at(i);
			if (value > th) {
				possiblePoints->push_back(new pair<int,int>(i,j));
			}
		}
	}
	//cout << "# of possible Peak: " << possiblePoints->size() << endl;
	/////////////////////////////////////////////////////////////////////////////////	
	
	////////Find local maxima containing more/equal points than area_size////////////////
    vector< vector< pair<int,int>* >* >* peak_areas = find_peak_areas(possiblePoints, measurement->data->data);
    //cout << "# of Peak Areas: " << peak_areas->size() << endl;
    
    vector < vector < short >* >* area_matrix = new vector < vector < short >* >;
    for_loop(a, sizeR){
        vector < short >* row_T = new vector < short >; 
        row_T->resize(sizeT, -1);
        area_matrix->push_back(row_T);        
    }
    
    for_loop(i, peak_areas->size()){
        for_loop(j, peak_areas->at(i)->size()){
            area_matrix->at(peak_areas->at(i)->at(j)->second)->at(peak_areas->at(i)->at(j)->first) = i;
        }
    }
    
	/* 
    ofstream am("area-matrix.dat");
    for_loop(i, area_matrix->size()){
        for_loop(j, area_matrix->at(i)->size()){
            am << area_matrix->at(i)->at(j) << " " ;
        }
        am << endl;
    }
    cout << "Matrix erstellt!" << endl;
    */


    vector< pair<int,int>* >* local_maxima = new vector< pair<int,int>* >;
    #pragma omp parallel for
    for_loop(i, peak_areas->size()){
		#pragma omp parallel for
        for_loop(j, peak_areas->at(i)->size()){
            int x = peak_areas->at(i)->at(j)->first;
            int y = peak_areas->at(i)->at(j)->second;
            double value = measurement->data->data->at(y)->at(x);     
            if( peak_areas->at(i)->size() >= area_size && points_in_area(i, peak_areas->at(i), x, y, area_matrix)){
				//look at values with position (x-1, y), (x+1, y), (x, y-1), (x, y+1), (x+1, y+1), (x-1, y-1), (x+1, y-1), (x-1, y+1)
				double value1 = measurement->data->data->at(y)->at(x-1);
				double value2 = measurement->data->data->at(y)->at(x+1);
				double value3 = measurement->data->data->at(y-1)->at(x);
				double value4 = measurement->data->data->at(y+1)->at(x);
				double value5 = measurement->data->data->at(y+1)->at(x+1);
				double value6 = measurement->data->data->at(y-1)->at(x-1);
				double value7 = measurement->data->data->at(y+1)->at(x-1);            
				double value8 = measurement->data->data->at(y-1)->at(x+1);
				if(value >= value1 && value >= value2 && value >= value3 && value >= value4 && 
					value >= value5 && value >= value6 && value >= value7 && value >= value8){
					#pragma omp critical(dataupdate)
					local_maxima->push_back(peak_areas->at(i)->at(j)); 
				}
            }
        }  
    }
    //cout << "# of Peak: " << local_maxima->size() << endl;
    ///////////////////////////////////////////////////////////////////////////////////
    
    ////Constructing new IMSPeaks and an IMSPeaklist////////////////////////////////////
    IMSPeakList* candidate_peaklist = new IMSPeakList;
    candidate_peaklist->measurement_source = measurement;
	for_loop(i,local_maxima->size()){
		IMSPeak* candidate_peak = new IMSPeak;
		candidate_peak->measurement_name = measurement->filename;
		candidate_peak->peak_name = "p" + i2s(i);
		candidate_peak->index_r = local_maxima->at(i)->second;
		candidate_peak->index_t = local_maxima->at(i)->first;
		candidate_peak->r = measurement->retention_times->at(local_maxima->at(i)->second);
		candidate_peak->t = measurement->reduced_inversed_mobilities->at(local_maxima->at(i)->first);
		candidate_peak->signal = measurement->data->data->at(local_maxima->at(i)->second)->at(local_maxima->at(i)->first);
		
		candidate_peaklist->ims_peak_list->push_back(candidate_peak);
	}
	//////////////////////////////////////////////////////////////////////////////////
	
	//////////////////Deleting pointer////////////////////////////////////////////////
	for_loop(i, possiblePoints->size()){
	    delete possiblePoints->at(i); 
	}
	delete possiblePoints;
	for_loop(i, peak_areas->size()){
        for_loop(j,peak_areas->at(i)->size()){
            delete peak_areas->at(i)->at(j);
        }
        delete peak_areas->at(i);
    }
    for_loop(i, area_matrix->size()){
        delete area_matrix->at(i);
    }
    delete area_matrix;
	delete peak_areas;
	delete local_maxima;
	//////////////////////////////////////////////////////////////////////////////////
	
	return candidate_peaklist;
}


candidate_detection(bit_parallel){
    int R = measurement->data->n_rows;
    int T = measurement->data->n_cols;
    
    uint* data = new uint[(T >> 5) + 1];
    uint pos = -1;
    int r = 189;
    
    for_loop (t, T - 1){
        if (!(t & 31)){
            data[++pos] = 0;
        }
        
        data[pos] |= (measurement->data->data->at(r)->at(t) < measurement->data->data->at(r)->at(t + 1)) << (t & 31);
    }
    
    const uint pattern = 31;
    uint text = 0;
    uint match;
    uint t = 0;
    while (t < T - 1){
        text = (data[t >> 5] >> (t & 31)) | (data[(t >> 5) + 1] << (31 - (t & 31)));
        
        int hits = __builtin_popcount((~(text ^ pattern)) & 511);
        pos = t + 3;
        if (hits >= 8 && measurement->data->data->at(r)->at(pos) >= 5){
            for (int i = t + 3; i < t + 8; i++){
                if (measurement->data->data->at(r)->at(pos) < measurement->data->data->at(r)->at(i)) pos = i;
            }
            cout << "hit: " << pos << endl;
            t += 8;
        }
        else {
            int w = 1;
            for (; w <= 29; w++){
                if (((text >> w) & 7) == 7) break;
            }
            t += w;
        }
    }
    
    
    delete []data;
    
    return 0;
}



void surroundByZero(IMSMatrix* matrix, int iterations) {

    // Insert 0s in all exisitng row at front and end
    for_loop (i, matrix->n_rows){
        for_loop (count, iterations) {
            matrix->data->at(i)->insert(matrix->data->at(i)->begin(), 0.0);
            matrix->data->at(i)->push_back(0.0);
        }
    }

    // Zeroline at beginning
    for_loop (count, iterations) {
        ims_spectrum* zeroList = new ims_spectrum;
        for_loop (i, matrix->n_cols + 2 * iterations)
            zeroList->push_back(0.0);
        matrix->data->insert(matrix->data->begin(), zeroList);
    }

    // Zeroline at end
    for_loop (count, iterations) {
        ims_spectrum* zeroList = new ims_spectrum;
        for_loop (i, matrix->n_cols + 2 * iterations)
            zeroList->push_back(0.0);
        matrix->data->push_back(zeroList);
    }
    
    matrix->n_cols += 2 * iterations;
    matrix->n_rows += 2 * iterations;
}


void shrinkByZero(IMSMatrix* matrix, int iterations) {
    matrix->n_cols -= 2 * iterations;
    matrix->n_rows -= 2 * iterations;
    for_loop (count, iterations) {
        delete matrix->data->front();
        matrix->data->erase(matrix->data->begin());
        
        delete matrix->data->back();
        matrix->data->pop_back();
    }
    
    for_loop (i, matrix->n_rows){
        for_loop (count, iterations) {
            matrix->data->at(i)->erase(matrix->data->at(i)->begin());
            matrix->data->at(i)->pop_back();
        }
    }
}


inline double get_score(int current, int new_p){
    return (fabs(current - new_p) > 5) ? 0 : 1. / (1. + fabs(current - new_p));
}


vector< align* >* global_alignment(vector< line* >* current, vector<int>* new_points, int dimension){
    vector< align* >* alignment = new vector< align* >;
    int n_current = current->size();
    int n_new = new_points->size();
    
    int c = 0, d = 0;
    
    // erstes Spektrum oder keine aktiven Peaks
    if (!current || !current->size()){
       
        for_loop (j, n_new){
            alignment->push_back(new align);
            alignment->back()->old_index = -1;
            alignment->back()->new_index = j;
            c++;
        }
    }
    else {
        // Liste "Daten" aufsteigend nach Modalwert sortieren
        for (int i = 0; i < current->size(); i++){
            for (int j = 1; j < current->size() - i; j++){
                if (dimension == X_ALIGN){
                    if (current->at(j - 1)->points->back()->x > current->at(j)->points->back()->x){
                        swap(current->at(j - 1), current->at(j));
                    }
                }
                if (dimension == Y_ALIGN){
                    if (current->at(j - 1)->points->back()->y > current->at(j)->points->back()->y){
                        swap(current->at(j - 1), current->at(j));
                    }
                }
            }
        }
        double **matrix = new double*[current->size() + 1];
        for (int i = 0; i <= n_current; i++){
                matrix[i] = new double[new_points->size() + 1];
                matrix[i][0] = 0;
        }
        for (int i = 0; i <= n_new; i++){
                matrix[0][i] = 0;
        }
                            

        align*** trace = new align**[n_current + 1];
        for (int i = 0; i <= n_current; i++){
                trace[i] = new align*[n_new + 1];
                trace[i][0] = new align;
                trace[i][0]->old_index = i - 1;
                trace[i][0]->new_index = 0;
        }               
        
        for (int i = 1; i <= n_new; i++){
                trace[0][i] = new align;
                trace[0][i]->old_index = 0;
                trace[0][i]->new_index = i - 1;
        }
        trace[0][0]->old_index = 0;
        trace[0][0]->new_index = 0;



        double score = 0;
        int type = 0;
        double best = 0;
        for (int i = 1; i <= n_current; i++){
            for (int j = 1; j <= n_new; j++){
                if (dimension == X_ALIGN){
                    score = get_score(current->at(i - 1)->points->back()->x, new_points->at(j - 1));
                }
                else if (dimension == Y_ALIGN){
                    score = get_score(current->at(i - 1)->points->back()->y, new_points->at(j - 1));
                }
                best = matrix[i - 1][j - 1] + score;
                type = 0;
                
                if (best < matrix[i][j - 1] + GAP){
                        best = matrix[i][j - 1] + GAP;
                        type = 1;
                }
                
                if (best < matrix[i - 1][j] + GAP){
                        best = matrix[i - 1][j] + GAP;
                        type = 2;
                }
                
                matrix[i][j] = best;
                trace[i][j] = new align;
                if (type == 0){
                        trace[i][j]->old_index = i - 1;
                        trace[i][j]->new_index = j - 1;
                }
                else if (type == 1){
                        trace[i][j]->old_index = i;
                        trace[i][j]->new_index = j - 1;
                }
                else if (type == 2){
                        trace[i][j]->old_index = i - 1;
                        trace[i][j]->new_index = j;
                }
            }

        }
            
        // Trace
        int ti = n_current;
        int tj = n_new;

        while (!(ti == 0 && tj == 0)){
            if (trace[ti][tj]->old_index != ti && trace[ti][tj]->new_index != tj){
                    alignment->insert(alignment->begin(), new align);
                    alignment->at(0)->old_index = --ti;
                    alignment->at(0)->new_index = --tj;
                    c++;
            }
            else if (trace[ti][tj]->old_index == ti){
                    alignment->insert(alignment->begin(), new align);
                    alignment->at(0)->old_index = -1;
                    alignment->at(0)->new_index = --tj;
                    c++;
            }
            else if (trace[ti][tj]->new_index == tj){
                    alignment->insert(alignment->begin(), new align);
                    alignment->at(0)->old_index = --ti;
                    alignment->at(0)->new_index = -1;
                    c++;
            }
        }
        

        for_loop (i, n_current + 1){
            for_loop (j, n_new + 1){
                delete trace[i][j];
            }
            delete []trace[i];
            d++;
        }
        delete []trace;
                        
        for_loop (i, n_current + 1){
            delete []matrix[i];
        }
        delete []matrix;
    }
    
    return alignment;    
}



/*!
 * \brief Finds peak candidates by searching the rootes of both derivated spectra and contour lines
 * @param measurement Instance of an IMSMeasurement
 * @param pipeline_parameters map of parameters (map<Key = string, Value = string>)
 * @return A list of peak candidates
 */
candidate_detection(cross_finding){
    // take parameters from parameter dictionary, if containing
    assert(pipeline_parameters->find("intensity_threshold") != pipeline_parameters->end());
    double min_peak_signal = atoi(pipeline_parameters->at("intensity_threshold").c_str());
    
    // above 2 no significant changes happen
    int margin_size = 2;
    
    IMSMatrix* matrix = measurement->data;
    surroundByZero(matrix, margin_size);
    

    
    uint R = matrix->n_rows;
    uint T = matrix->n_cols;

    // Create matrices D and E for derivative
    IMSMatrix* D = new IMSMatrix(matrix);
    IMSMatrix* E = new IMSMatrix(matrix);

    
    // Derive
    #pragma omp parallel for
    for_loop (r, R) {
        for (uint t = 1; t < T; t++){
            D->data->at(r)->at(t - 1) = D->data->at(r)->at(t) - D->data->at(r)->at(t - 1);
        }
        D->data->at(r)->back() = 0;
    }

    // Derive
    #pragma omp parallel for
    for_loop (t, T) {
        for (uint r = 1; r < R; r++){
            E->data->at(r - 1)->at(t) = E->data->at(r)->at(t) - E->data->at(r - 1)->at(t);
        }
        E->data->back()->at(t) = 0;
    }
 
    
    vector< line* >* valid_d = new vector< line* >;
    vector< line* >* valid_e = new vector< line* >;
    
    #pragma omp parallel sections
    {{
    vector< line* >* vectors_d = new vector< line* >;
    vector< int >* new_points_d = new vector< int >;
    vector< int > dels_d;
    for_loop (r, R){
        new_points_d->clear();
        dels_d.clear();
        for (uint t = 1; t < T - 1; t++){
            if (D->data->at(r)->at(t - 1) > D->data->at(r)->at(t) &&
                D->data->at(r)->at(t - 1) >= 0 &&
                D->data->at(r)->at(t) < 0 &&
                matrix->data->at(r)->at(t) > 1 &&
                (matrix->data->at(r)->at(t - 1) > 1 || matrix->data->at(r)->at(t + 1) > 1)){
                new_points_d->push_back(t);
            }
        }
        vector< align* >* alignment = global_alignment(vectors_d, new_points_d, X_ALIGN);
        for_loop(i, alignment->size()){
            if (alignment->at(i)->old_index != -1 && alignment->at(i)->new_index != -1){
                int le = min(vectors_d->at(alignment->at(i)->old_index)->points->back()->x, new_points_d->at(alignment->at(i)->new_index));
                int re = max(vectors_d->at(alignment->at(i)->old_index)->points->back()->x, new_points_d->at(alignment->at(i)->new_index));
                
                for (int j = le; j <= re; j++){
                    vectors_d->at(alignment->at(i)->old_index)->points->push_back(new point);
                    vectors_d->at(alignment->at(i)->old_index)->points->back()->x = j;
                    vectors_d->at(alignment->at(i)->old_index)->points->back()->y = r - 1;
                    
                    
                    vectors_d->at(alignment->at(i)->old_index)->ll.x = min((int)vectors_d->at(alignment->at(i)->old_index)->ll.x, j);
                    vectors_d->at(alignment->at(i)->old_index)->ll.y = min((int)vectors_d->at(alignment->at(i)->old_index)->ll.y, (int)(r - 1));
                    vectors_d->at(alignment->at(i)->old_index)->ur.x = max((int)vectors_d->at(alignment->at(i)->old_index)->ur.x, j);
                    vectors_d->at(alignment->at(i)->old_index)->ur.y = max((int)vectors_d->at(alignment->at(i)->old_index)->ur.y, (int)(r));
                    
                    vectors_d->at(alignment->at(i)->old_index)->points->push_back(new point);
                    vectors_d->at(alignment->at(i)->old_index)->points->back()->x = j;
                    vectors_d->at(alignment->at(i)->old_index)->points->back()->y = r;
                    
                }
                
                vectors_d->at(alignment->at(i)->old_index)->points->push_back(new point);
                vectors_d->at(alignment->at(i)->old_index)->points->back()->x = new_points_d->at(alignment->at(i)->new_index);
                vectors_d->at(alignment->at(i)->old_index)->points->back()->y = r;
                
                
            }
            else if (alignment->at(i)->old_index == -1){
                vectors_d->push_back(new line);
                vectors_d->back()->points = new vector< point* >;
                vectors_d->back()->points->push_back(new point);
                vectors_d->back()->points->back()->x = new_points_d->at(alignment->at(i)->new_index);
                vectors_d->back()->points->back()->y = r;
                vectors_d->back()->ll.x = new_points_d->at(alignment->at(i)->new_index);
                vectors_d->back()->ll.y = r;
                vectors_d->back()->ur.x = new_points_d->at(alignment->at(i)->new_index);
                vectors_d->back()->ur.y = r;
                
            }
            else if (alignment->at(i)->new_index == -1){
                int index = alignment->at(i)->old_index;
                if (vectors_d->at(index)->points->back()->y - vectors_d->at(index)->points->front()->y >= 4){
                    valid_d->push_back(vectors_d->at(index));
                }
                else {
                    for_loop(j, vectors_d->at(index)->points->size()){
                        delete vectors_d->at(index)->points->at(j);
                    }
                    delete vectors_d->at(index)->points;
                    delete vectors_d->at(index);
                }
                dels_d.push_back(index);
            }
            delete alignment->at(i);
        }
        delete alignment;
        for (int i = dels_d.size() - 1; i >= 0; i--){
            vectors_d->erase(vectors_d->begin() + dels_d.at(i));
            
        }
    }
    for (int i = vectors_d->size() - 1; i >= 0; i--){
        if (vectors_d->at(i)->points->back()->y - vectors_d->at(i)->points->front()->y >= 4){
            valid_d->push_back(vectors_d->at(i));
        }
        else {
            for_loop(j, vectors_d->at(i)->points->size()){
                delete vectors_d->at(i)->points->at(j);
            }
            delete vectors_d->at(i)->points;
            delete vectors_d->at(i);
        }
        vectors_d->erase(vectors_d->begin() + i);
    }
    delete vectors_d;
    delete new_points_d;
    } // omp
    
    #pragma omp section
    { // omp
    vector< line* >* vectors_e = new vector< line* >;
    vector< int >* new_points_e = new vector< int >;
    vector< int > dels_e;
    for_loop (t, T){
        new_points_e->clear();
        dels_e.clear();
        for (uint r = 1; r < R - 1; r++){
            if (E->data->at(r - 1)->at(t) > E->data->at(r)->at(t) &&
                E->data->at(r - 1)->at(t) >= 0 &&
                E->data->at(r)->at(t) < 0 &&
                matrix->data->at(r)->at(t) > 1 &&
                (matrix->data->at(r - 1)->at(t) > 1 || matrix->data->at(r + 1)->at(t) > 1)){
                new_points_e->push_back(r);
            }
        }
        vector< align* >* alignment = global_alignment(vectors_e, new_points_e, Y_ALIGN);
        
        for_loop(i, alignment->size()){
            if (alignment->at(i)->old_index != -1 && alignment->at(i)->new_index != -1){
                int le = min(vectors_e->at(alignment->at(i)->old_index)->points->back()->y, new_points_e->at(alignment->at(i)->new_index));
                int re = max(vectors_e->at(alignment->at(i)->old_index)->points->back()->y, new_points_e->at(alignment->at(i)->new_index));
                
                for (int j = le; j <= re; j++){
                    vectors_e->at(alignment->at(i)->old_index)->points->push_back(new point);
                    vectors_e->at(alignment->at(i)->old_index)->points->back()->x = t - 1;
                    vectors_e->at(alignment->at(i)->old_index)->points->back()->y = j;
                    vectors_e->at(alignment->at(i)->old_index)->points->push_back(new point);
                    vectors_e->at(alignment->at(i)->old_index)->points->back()->x = t;
                    vectors_e->at(alignment->at(i)->old_index)->points->back()->y = j;
                    
                    vectors_e->at(alignment->at(i)->old_index)->ll.x = min((int)vectors_e->at(alignment->at(i)->old_index)->ll.x, (int)(t - 1));
                    vectors_e->at(alignment->at(i)->old_index)->ll.y = min((int)vectors_e->at(alignment->at(i)->old_index)->ll.y, (int)(j));
                    vectors_e->at(alignment->at(i)->old_index)->ur.x = max((int)vectors_e->at(alignment->at(i)->old_index)->ur.x, (int)(t));
                    vectors_e->at(alignment->at(i)->old_index)->ur.y = max((int)vectors_e->at(alignment->at(i)->old_index)->ur.y, (int)(j));
                    
                }
                
                vectors_e->at(alignment->at(i)->old_index)->points->push_back(new point);
                vectors_e->at(alignment->at(i)->old_index)->points->back()->x = t;
                vectors_e->at(alignment->at(i)->old_index)->points->back()->y = new_points_e->at(alignment->at(i)->new_index);
                
            }
            else if (alignment->at(i)->old_index == -1){
                vectors_e->push_back(new line);
                vectors_e->back()->points = new vector< point* >;
                vectors_e->back()->points->push_back(new point);
                vectors_e->back()->points->back()->x = t;
                vectors_e->back()->points->back()->y = new_points_e->at(alignment->at(i)->new_index);
                vectors_e->back()->ll.x = t;
                vectors_e->back()->ll.y = new_points_e->at(alignment->at(i)->new_index);
                vectors_e->back()->ur.x = t;
                vectors_e->back()->ur.y = new_points_e->at(alignment->at(i)->new_index);
            }
            else if (alignment->at(i)->new_index == -1){
                int index = alignment->at(i)->old_index;
                if (vectors_e->at(index)->points->back()->x - vectors_e->at(index)->points->front()->x >= 4){
                    valid_e->push_back(vectors_e->at(index));
                }
                else {
                    for_loop(j, vectors_e->at(index)->points->size()){
                        delete vectors_e->at(index)->points->at(j);
                    }
                    delete vectors_e->at(index)->points;
                    delete vectors_e->at(index);
                }
                dels_e.push_back(index);
            }
            delete alignment->at(i);
        }
        delete alignment;
        for (int i = dels_e.size() - 1; i >= 0; i--){
            vectors_e->erase(vectors_e->begin() + dels_e.at(i));
        }
    }
    for (int i = vectors_e->size() - 1; i >= 0; i--){
        if (vectors_e->at(i)->points->back()->x - vectors_e->at(i)->points->front()->x >= 4){
            valid_d->push_back(vectors_e->at(i));
        }
        else {
            for_loop(j, vectors_e->at(i)->points->size()){
                delete vectors_e->at(i)->points->at(j);
            }
            delete vectors_e->at(i)->points;
            delete vectors_e->at(i);
        }
        vectors_e->erase(vectors_e->begin() + i);
    }
    delete vectors_e;
    delete new_points_e;
    }} // omp
    
    delete D;
    delete E;
    
    vector < point* >* a_points = new vector< point* >;
    
    #pragma omp parallel for collapse(2)
    for_loop(i, valid_d->size()){
        for_loop(j, valid_e->size()){
            line* l1 = valid_d->at(i);
            if (l1->ur.y - l1->ll.y <= 3) continue;
            line* l2 = valid_e->at(j);
            if (l2->ur.x - l2->ll.x <= 3) continue;
            if (!(l1->ll.x > l2->ur.x || l1->ur.x < l2->ll.x || l1->ll.y > l2->ur.y || l1->ur.y < l2->ll.y)){
                int x1 = -1, x2 = -1, y1 = -1, y2 = -1;
                int hx = -1, hy = -1;
                double hs = -1;
                for(uint m = 0; m < l1->points->size(); m++){
                    x1 = l1->points->at(m)->x;
                    y1 = l1->points->at(m)->y;
                    for_loop(l, l2->points->size()){
                        x2 = l2->points->at(l)->x;
                        y2 = l2->points->at(l)->y;
                        if (x1 == x2 && y1 == y2 && measurement->data->data->at(y1)->at(x1) >= min_peak_signal && measurement->data->data->at(y1)->at(x1) > hs)
                        {
                            hs = measurement->data->data->at(y1)->at(x1);
                            hx = x1;
                            hy = y1;
                        }
                    }
                }
                if (hs > 0 && hx >= 0 && hy >= 0){
                    #pragma omp critical(dataupdate)
                    {
                        a_points->push_back(new point);
                        a_points->back()->x = hx;
                        a_points->back()->y = hy;
                    }
                }
            }
        }
    }
    for_loop (i, valid_d->size()){
        for_loop(j, valid_d->at(i)->points->size()){
            delete valid_d->at(i)->points->at(j);
        }
        delete valid_d->at(i)->points;
        delete valid_d->at(i);
    }
    delete valid_d;
    for_loop (i, valid_e->size()){
        for_loop(j, valid_e->at(i)->points->size()){
            delete valid_e->at(i)->points->at(j);
        }
        delete valid_e->at(i)->points;
        delete valid_e->at(i);
    }
    delete valid_e;
        
    shrinkByZero(matrix, margin_size);
    
    IMSPeakList* peaklist = new IMSPeakList;
    peaklist->measurement_source = measurement;
    for_loop (i, a_points->size()){
        if (a_points->at(i)->x - margin_size >= 0 && a_points->at(i)->y - margin_size >= 0){
            IMSPeak* peak = new IMSPeak;
            peak->measurement_name = measurement->filename;
            stringstream ss; 
            ss << (i + 1);
            peak->peak_name = "p" + ss.str();
            peak->r = a_points->at(i)->y - margin_size;
            peak->t = a_points->at(i)->x - margin_size;
            peak->signal = matrix->data->at(peak->r)->at(peak->t);
            peak->index_r = peak->r;
            peak->index_t = peak->t;
            peak->r = measurement->retention_times->at(peak->r);
            peak->t = measurement->reduced_inversed_mobilities->at(peak->t);
            peak->volume = peak->signal;
            peaklist->ims_peak_list->push_back(peak);
        }
        delete a_points->at(i);
    }
    delete a_points;
    
    return peaklist;
}
