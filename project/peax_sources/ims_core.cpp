#include "ims_core.h"


/*!
 * \brief Own implementation of atof which is faster than atof, needed for loading files
 * @param string s
 * @return float
 */
float s2f(string s) {
    const char* p = s.c_str();
    int len = s.length();
    int r = 0.0;
    bool neg = false;
    float n = 1;
    int i = 0;
    // check for leading sign
    if (*p == '-') {
        neg = true;
        ++p;
    }
    // process digits
    while (i < len && *p >= '0' && *p <= '9') {
        r = (r * 10) + *p - 48;
        p++;
        i++;
    }
    // if next letter is a decimal point, go ahead
    if (i < len && *p != '.') {
        return neg ? -r : r;
    }
    p++;
    i++;
    // process digits after decimal point
    while (i < len && *p >= '0' && *p <= '9') { 
        r = (r * 10) + *p - 48;
        p++;
        n *= 10;
        i++;
    }
    // assemble number
    return neg ? -float(r) / float(n) : float(r) / float(n);
}

/*!
 * \brief Own implementation of integer -> string which is faster than stringstream, needed for storing files
 * @param integer s
 * @return string
 */
string i2s(int val){
        if (val == 0) return "0";
        string str = "";
        for(int i = floor(log10(val) + 1); i > 0; i--) {
                str += (char)((int)((val % (int)pow(10, i)) / (double)pow(10, i - 1)) + 48);
        }
        return str;  
}

/*!
 * \brief Constructor for class IMSPeak
 */
IMSPeak::IMSPeak(){
    measurement_name = "";
    peak_name = "";
    r = 0;
    t = 0;
    signal = 0;
    volume = 0;
    index_r = 0;
    index_t = 0;
    peak_parameters = new map<string, double>;
}

/*!
 * \brief Copy constructor for class IMSPeak
 */
IMSPeak::IMSPeak(const IMSPeak* source){
    measurement_name = source->measurement_name;
    peak_name = source->peak_name;
    r = source->r;
    t = source->t;
    signal = source->signal;
    volume = source->volume;
    index_r = source->index_r;
    index_t = source->index_t;
    peak_parameters = new map<string, double>;

    for(map<string,double>::iterator it = source->peak_parameters->begin(); it != source->peak_parameters->end(); ++it){
        peak_parameters->insert(pair<string,double>(it->first, it->second));
    }
}

/*!
 * \brief Destructor for class IMSPeak
 */
IMSPeak::~IMSPeak(){
    delete peak_parameters;
}


void IMSPeak::input_parameter(string s, double d){
    peak_parameters->insert(pair<string,double>(s, d));
}


/////////////////////////////////////////////////////////////////////////////////


/*!
 * \brief Constructor for class IMSPeakList
 */
IMSPeakList::IMSPeakList(){
    ims_peak_list = new vector<IMSPeak*>;
    parameter_names = new stringlist;
    measurement_source = 0;
}

/*!
 * \brief Copy constructor for class IMSPeakList
 */
IMSPeakList::IMSPeakList(const IMSPeakList* source){
    ims_peak_list = new vector<IMSPeak*>;
    parameter_names = new stringlist;
    for_loop(i, source->ims_peak_list->size()){
        ims_peak_list->push_back(new IMSPeak(source->ims_peak_list->at(i)));
    }
    for_loop (s, source->parameter_names->size()){
        parameter_names->push_back(source->parameter_names->at(s));
    }
    
    measurement_source = source->measurement_source;
}

/*!
 * \brief Destructor for class IMSPeakList
 */
IMSPeakList::~IMSPeakList(){
    for_loop(i, ims_peak_list->size()){
        delete ims_peak_list->at(i);
    }
    delete ims_peak_list;
    delete parameter_names;
}

/*!
 * \brief Parse peak list file
 * @param string input_filename
 */
void IMSPeakList::load(string input_filename){
    // Opening file
    ifstream peak_list_file(input_filename.c_str());
    
    // check wheather file exists
    if (peak_list_file.is_open()){
        
        string line = "";
        
        // check for additional parameters and store
        getline(peak_list_file, line);
        stringlist column_names = split(line, '\t');
        if (column_names.size() > 7){
            for (uint i = 7; i < column_names.size(); i++){
                parameter_names->push_back(column_names.at(i));
            }
        }
        
        // reading all lines
        uint n_lines = 1;
        if (column_names.size() >= 7){
            while (!peak_list_file.eof()){
                n_lines++;
                getline(peak_list_file, line);
                string test_line = strip(line);
                // empty line check
                if (test_line.length() == 0) continue;
                // only line will be read when it contains at least 5 entries
                stringlist values = split(line, '\t');
                if (values.size() >= 7){
                    IMSPeak* peak = new IMSPeak;
                    peak->measurement_name = values.at(0);
                    peak->peak_name = values.at(1);
                    peak->t = atof(values.at(2).c_str());
                    peak->r = atof(values.at(3).c_str());
                    peak->signal = atof(values.at(4).c_str());
                    peak->volume = atof(values.at(4).c_str());
                    peak->index_t = atoi(values.at(5).c_str());
                    peak->index_r = atoi(values.at(6).c_str());
            
                    // additional entries will be stored as peak parameters
                    if (column_names.size() > 7){
                        for (uint i = 7; i < column_names.size(); i++){
                            peak->input_parameter(column_names.at(i), atof(values.at(i).c_str()));
                        }
                    }
                    ims_peak_list->push_back(peak);
                }
                else {
                    cout << "Error: line " << n_lines << " in file " << input_filename << " has incompatible format" << endl;
                }
            }
        }
        else {
            cout << "Error: " << input_filename << " has incompatible format" << endl;
        }
        peak_list_file.close();
    }
    else {
        cout << "Error: " << input_filename << " could not be opened" << endl;
    }
}

/*!
 * \brief Store peak list in file
 * @param string output_filename
 */
void IMSPeakList::save(string output_filename){
    // opening output stream
    ofstream output_list(output_filename.c_str());
    // creating first row with column names
    output_list << "measurement_name\tpeak_name\tt\tr\tsignal\tindex_t\tindex_r";
    if (parameter_names->size() > 0){
        for_loop(i, parameter_names->size()){
            output_list << "\t" << parameter_names->at(i);
        }
    }
    output_list << endl;
    
    // writing all peak values out
    for_loop(i, ims_peak_list->size()){
        output_list << ims_peak_list->at(i)->measurement_name << "\t";
        output_list << ims_peak_list->at(i)->peak_name << "\t";
        output_list << ims_peak_list->at(i)->t << "\t";
        output_list << ims_peak_list->at(i)->r << "\t";
        output_list << ims_peak_list->at(i)->signal << "\t";
        output_list << ims_peak_list->at(i)->index_t << "\t";
        output_list << ims_peak_list->at(i)->index_r;
        
        // writing if existing all additional peak parameters out
        if (parameter_names->size() > 0){
            for_loop(j, parameter_names->size()){
                output_list << "\t" << ims_peak_list->at(i)->peak_parameters->at(parameter_names->at(j));
            }
        }
        output_list << endl;
    }
    output_list.close();
}



/////////////////////////////////////////////////////////////////////////////////


/*!
 * \brief Constructor for class IMSMatrix
 */
IMSMatrix::IMSMatrix(){
    n_rows = 0;
    n_cols = 0;
    data = new matrix_t;
}

/*!
 * \brief Copy constructor for class IMSMatrix
 */
IMSMatrix::IMSMatrix(const IMSMatrix* source){
    n_rows = source->n_rows;
    n_cols = source->n_cols;
    data = new matrix_t;
    
    for_loop(i, n_rows){
        data->push_back(new ims_spectrum);
        for_loop(j, n_cols){
            data->at(i)->push_back(source->data->at(i)->at(j));
        }
    }
}

/*!
 * \brief Destructor for class IMSMatrix
 */
IMSMatrix::~IMSMatrix(){
    for_loop(i, n_rows){
        delete data->at(i);
    }
    delete data;
}


/////////////////////////////////////////////////////////////////////////////////


/*!
 * \brief Constructor for class IMSMeasurement
 */
IMSMeasurement::IMSMeasurement(string _filename){
    stringlist dirs = split(_filename, '/');
    filename = dirs.back();
    retention_times = new times;
    drift_times = new times;
    reduced_inversed_mobilities = new times;
    data = new IMSMatrix;
    measurement_parameters = new map<string, string>;
}

/*!
 * \brief Copy constructor for class IMSMeasurement
 */
IMSMeasurement::IMSMeasurement(const IMSMeasurement* source) : filename(source->filename){
    retention_times = new times;
    drift_times = new times;
    reduced_inversed_mobilities = new times;
    measurement_parameters = new map<string, string>;
    
    for_loop(r, source->retention_times->size()){
        retention_times->push_back(source->retention_times->at(r));
    }
    for_loop(d, source->drift_times->size()){
        drift_times->push_back(source->drift_times->at(d));
    }
    for_loop(t, source->reduced_inversed_mobilities->size()){
        reduced_inversed_mobilities->push_back(source->reduced_inversed_mobilities->at(t));
    }
    
    data = new IMSMatrix(source->data);
    
    for(map<string,string>::iterator it = source->measurement_parameters->begin(); it != source->measurement_parameters->end(); ++it){
        measurement_parameters->insert(pair<string,string>(it->first, it->second));
    }
}

/*!
 * \brief Destructor for class IMSMeasurement
 */
IMSMeasurement::~IMSMeasurement(){
    delete retention_times;
    delete drift_times;
    delete reduced_inversed_mobilities;
    delete data;
    delete measurement_parameters;
}

void IMSMeasurement::input_parameter(string s, string p){
    measurement_parameters->insert(pair<string,string>(s, p));
}

/*!
 * \brief Computes the matrix retention index in retention axes having a particular measurement
 * @param double r
 * @return int
 */
int IMSMeasurement::getRetentionIndex(double r){
    assert (!isnan(r));
    assert (!isinf(r));
    assert (r >= retention_times->front());
    assert (retention_times->back() >= r);
    int left = 0;
    int right = data->n_rows;
    int curr;
    
    while (true){
        curr = (left + right) >> 1;
        if (retention_times->at(curr) <= r){
            left = curr;
        }
        else {
            right = curr;
        }
        if (right - left <= 1){
            break;
        }
    }
    return curr;
}

/*!
 * \brief Computes the matrix drift index in drift axes having a particular measurement
 * @param double d
 * @return int
 */
int IMSMeasurement::getDriftIndex(double d){
    assert (!isnan(d));
    assert (!isinf(d));
    assert (d >= drift_times->front());
    assert (drift_times->back() >= d);
    int left = 0;
    int right = data->n_cols;
    int curr;
    
    while (true){
        curr = (left + right) >> 1;
        if (drift_times->at(curr) <= d){
            left = curr;
        }
        else {
            right = curr;
        }
        if (right - left <= 1){
            break;
        }
    }
    return curr;
}

/*!
 * \brief Computes the matrix RIM index in drift axes having a particular measurement
 * @param double t
 * @return int
 */
int IMSMeasurement::getRIMIndex(double t){
    assert (!isnan(t));
    assert (!isinf(t));
    assert (t >= reduced_inversed_mobilities->front());
    assert (reduced_inversed_mobilities->back() >= t);
    
    int left = 0;
    int right = data->n_cols;
    int curr;
    
    while (true){
        curr = (left + right) >> 1;
        if (reduced_inversed_mobilities->at(curr) <= t){
            left = curr;
        }
        else {
            right = curr;
        }
        if (right - left <= 1){
            break;
        }
    }
    return curr;
}

/*!
 * \brief Computes retention time having a particular measurement and matrix retention index
 * @param double index_r
 * @return double
 */
double IMSMeasurement::computeRetention(double index_r){
    
    int i_index_r = index_r;
    assert (!isnan(index_r));
    assert (!isinf(index_r));
    assert (i_index_r >= 0);
    assert (i_index_r < retention_times->size());
    
    double r = retention_times->at(i_index_r);
    if (i_index_r + 1 < data->n_rows){
        r += (retention_times->at(i_index_r + 1) - retention_times->at(i_index_r)) * (index_r - i_index_r);
    }
    return r;
}

/*!
 * \brief Computes RIM having a particular measurement and matrix RIM index
 * @param double index_rim
 * @param bool use_start_offset
 * @return double
 */
double IMSMeasurement::computeInverseMobility(double index_t, bool use_start_offset){
    int i_index_t = index_t;
    assert (!isnan(index_t));
    assert (!isinf(index_t));
    assert (i_index_t >= 0);
    assert (i_index_t < reduced_inversed_mobilities->size());
    
    double t = reduced_inversed_mobilities->at(i_index_t);
    if (i_index_t + 1 < (int)data->n_cols){
        t += (reduced_inversed_mobilities->at(i_index_t + 1) - reduced_inversed_mobilities->at(i_index_t)) * (index_t - i_index_t);
    }
    if (use_start_offset){
        t -= reduced_inversed_mobilities->front();
    }
    return t;
}

/////////////////////////////////////////////////////////////////////////////////


/*!
 * \brief Parse MCC/IMS measurement in csv format from file
 * @param string filename
 * @return IMSMeasurement*
 */
IMSMeasurement* IMSFileCSV::load(string filename){
    ifstream ims_file(filename.c_str());
    if (!ims_file.is_open()){
        cout << "Error: " << filename << " could not be opened" << endl;
        return 0;
    }
    
    int param_lines[] = {0, 1, 2, 3, 4, 6, 7, 8, 13, 14, 15, 16, 17, 24, 25, 26, 30, 49, 51, 53, 55, 56, 60, 67, 75, 80, 88, 89, 99, 18, 19, 28, 29, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 47, 50, 52, 54, 57, 58, 61, 62, 63, 64, 65, 66, 68, 69, 70, 71, 72, 73, 81, 82, 83, 84, 85, 86, 87, 90, 91, 92, 93, 94, 100, 101, 102, 103, 104, 105, 107, 108, 109, 110, 111, 112, 114, 115, 116, 119};
    set<int> set_param_lines(param_lines, param_lines + 93);
   
    string line;
    stringlist* lineList = new stringlist;
    uint n_lines = 0;
    IMSMeasurement* measurement = new IMSMeasurement(filename);
    
    uint T = 0;
    while(!ims_file.eof()){
        getline(ims_file, line);
        lineList->push_back(line);
        if (n_lines >= 132 && line.length() && line.at(0) != ' '){
            T++;
        }
        n_lines++;
    }
    ims_file.close();
    
    if (lineList->at(0)[0] == 'i' && lineList->at(0)[1] == 'm' && lineList->at(0)[2] == 's'){
        delete lineList;
        return 0;
    }
    
    // Meta data
    for_loop(i, 131){  
        line = lineList->at(i);
        stringlist qsp = split(line, ',');
        while(qsp.size() && qsp.back() == "") qsp.pop_back();
        if (line.length() < 2 || qsp.size() < 1) continue;
        
        if (i == 130){
            if (qsp.size() < 2){
                cout << "Error: line " << i << " incomplete"  << endl;
                return 0;
            }
            for (uint ii = 2; ii < qsp.size(); ii++){
                measurement->retention_times->push_back(atof(qsp.at(ii).c_str()));
                measurement->data->data->push_back(new ims_spectrum(T, 0));
            }
        }
        else {
            if (qsp.size() > 2){
                if (set_param_lines.find(i) != set_param_lines.end()){
                    measurement->input_parameter(qsp.at(1), qsp.at(2));
                }
            }
        }
    }
    
    bool error = false;
    measurement->reduced_inversed_mobilities->resize(T, 0);
    measurement->drift_times->resize(T, 0);
    
    #pragma omp parallel for
    for (uint i = 132; i < n_lines; i++){
        
        string line = lineList->at(i);
        stringlist qsp = split(line, ',');
        while(qsp.size() && qsp.back() == "") qsp.pop_back();
        if (line.length() < 2 || qsp.size() < 1) continue;
        
        if (qsp.size() < 2){
            cout << "Error: line " << i << " incomplete"  << endl;
            error = true;
        }
        
        measurement->reduced_inversed_mobilities->at(i - 132) = s2f(qsp.at(0).c_str());
        measurement->drift_times->at(i - 132) = s2f(qsp.at(1).c_str());
        for (uint ii = 2; ii < qsp.size(); ii++){
            measurement->data->data->at(ii - 2)->at(i - 132) = -s2f(qsp.at(ii).c_str());
        }
    }
    
    if(error) return 0;
    
    delete lineList;
    
    measurement->data->n_rows = measurement->data->data->size();
    measurement->data->n_cols = measurement->data->data->at(0)->size();
    return measurement;
}

/*!
 * \brief Store MCC/IMS measurement in csv format in file
 * @param IMSMeasurement* measurement
 * @param string filename
 */
void IMSFileCSV::save(IMSMeasurement* measurement, string filename){
    ofstream file(filename.c_str());
    file << "#,data type," << measurement->measurement_parameters->at("data type");
    for_loop (i, measurement->data->n_rows) file << ",";
    file << endl;
    
    file << "#,version," << measurement->measurement_parameters->at("version");
    for_loop (i, measurement->data->n_rows) file << ",";
    file << endl;
    
    file << "#,template version," << measurement->measurement_parameters->at("template version");
    for_loop (i, measurement->data->n_rows) file << ",";
    file << endl;
    
    file << "#,AD-board type," << measurement->measurement_parameters->at("AD-board type");
    for_loop (i, measurement->data->n_rows) file << ",";
    file << endl;
    
    file << "#,ser.-no.," << measurement->measurement_parameters->at("ser.-no.");
    for_loop (i, measurement->data->n_rows) file << ",";
    file << endl;
    
    
    file << "#" << endl;
    file << "#,date," << measurement->measurement_parameters->at("date") << endl;
    file << "#,time," << measurement->measurement_parameters->at("time") << endl;
    file << "#,file," << measurement->measurement_parameters->at("file") << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#,SAMPLE INFORMATION," << endl;
    file << "#" << endl;
    file << "#,sample type," << measurement->measurement_parameters->at("sample type") << endl;
    file << "#,sample ID," << measurement->measurement_parameters->at("sample ID") << endl;
    file << "#,comment," << measurement->measurement_parameters->at("comment") << endl;
    file << "#,location," << measurement->measurement_parameters->at("location") << endl;
    file << "#,location name," << measurement->measurement_parameters->at("location name") << endl;
    file << "#,height ASL / m," << measurement->measurement_parameters->at("height ASL / m") << endl;
    file << "#,total data acquisition time / s," << measurement->measurement_parameters->at("total data acquisition time / s") << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#,IMS - INFORMATION," << endl;
    file << "#" << endl;
    file << "#,operator," << measurement->measurement_parameters->at("operator") << endl;
    file << "#,operator name," << measurement->measurement_parameters->at("operator name") << endl;
    file << "#,IMS," << measurement->measurement_parameters->at("IMS") << endl;
    file << "#" << endl;
    file << "#,K0 RIP positive / cm^2/Vs," << measurement->measurement_parameters->at("K0 RIP positive / cm^2/Vs") << endl;
    file << "#,K0 RIP negative / cm^2/Vs," << measurement->measurement_parameters->at("K0 RIP negative / cm^2/Vs") << endl;
    file << "#,polarity," << measurement->measurement_parameters->at("polarity") << endl;
    file << "#,grid opening time / us," << measurement->measurement_parameters->at("grid opening time / us") << endl;
    file << "#" << endl;
    file << "#,pause / s," << measurement->measurement_parameters->at("pause / s") << endl;
    file << "#,tD interval (corr.) / ms from," << measurement->measurement_parameters->at("tD interval (corr.) / ms from") << endl;
    file << "#,tD interval (corr.) / ms to," << measurement->measurement_parameters->at("tD interval (corr.) / ms to") << endl;
    file << "#,1/K0 interval / Vs/cm^2 from," << measurement->measurement_parameters->at("1/K0 interval / Vs/cm^2 from") << endl;
    file << "#,1/K0 interval / Vs/cm^2 to," << measurement->measurement_parameters->at("1/K0 interval / Vs/cm^2 to") << endl;
    file << "#,no. of data points per spectra," << measurement->measurement_parameters->at("no. of data points per spectra") << endl;
    file << "#,no. of spectra," << measurement->measurement_parameters->at("no. of spectra") << endl;
    file << "#,no. averaged spectra," << measurement->measurement_parameters->at("no. averaged spectra") << endl;
    file << "#,baseline / signal units," << measurement->measurement_parameters->at("baseline / signal units") << endl;
    file << "#,baseline / V," << measurement->measurement_parameters->at("baseline / V") << endl;
    file << "#,V / signal unit," << measurement->measurement_parameters->at("V / signal unit") << endl;
    file << "#" << endl;
    file << "#,drift length / mm," << measurement->measurement_parameters->at("drift length / mm") << endl;
    file << "#,HV / kV," << measurement->measurement_parameters->at("HV / kV") << endl;
    file << "#,amplification / V/nA," << measurement->measurement_parameters->at("amplification / V/nA") << endl;
    file << "#" << endl;
    
    
    file << "#,drift gas," << measurement->measurement_parameters->at("drift gas") << endl;
    file << "#,drift gas flow / mL/min," << measurement->measurement_parameters->at("drift gas flow / mL/min") << endl;
    file << "#,sample gas," << measurement->measurement_parameters->at("sample gas") << endl;
    file << "#,sample gas flow / mL/min," << measurement->measurement_parameters->at("sample gas flow / mL/min") << endl;
    file << "#,carrier gas," << measurement->measurement_parameters->at("carrier gas") << endl;
    file << "#,carrier gas flow / mL/min," << measurement->measurement_parameters->at("carrier gas flow / mL/min") << endl;
    file << "#,pre-separation type," << measurement->measurement_parameters->at("pre-separation type") << endl;
    file << "#,pre-separation T / deg C," << measurement->measurement_parameters->at("pre-separation T / deg C") << endl;
    file << "#,sample loop T / deg C," << measurement->measurement_parameters->at("sample loop T / deg C") << endl;
    file << "#,sample loop volume / mL," << measurement->measurement_parameters->at("sample loop volume / mL") << endl;
    file << "#" << endl;
    
    file << "#,ambient T source," << measurement->measurement_parameters->at("ambient T source") << endl;
    file << "#,ambient T / deg C," << measurement->measurement_parameters->at("ambient T / deg C") << endl;
    file << "#,ambient T x^2," << measurement->measurement_parameters->at("ambient T x^2") << endl;
    file << "#,ambient T x^1," << measurement->measurement_parameters->at("ambient T x^1") << endl;
    file << "#,ambient T x^0," << measurement->measurement_parameters->at("ambient T x^0") << endl;
    file << "#,ambient T x^-1," << measurement->measurement_parameters->at("ambient T x^-1") << endl;
    file << "#,ambient T x^-2," << measurement->measurement_parameters->at("ambient T x^-2") << endl;
    file << "#,ambient p source," << measurement->measurement_parameters->at("ambient p source") << endl;
    file << "#,ambient p / hPa," << measurement->measurement_parameters->at("ambient p / hPa") << endl;
    file << "#,ambient p x^2," << measurement->measurement_parameters->at("ambient p x^2") << endl;
    file << "#,ambient p x^1," << measurement->measurement_parameters->at("ambient p x^1") << endl;
    file << "#,ambient p x^0," << measurement->measurement_parameters->at("ambient p x^0") << endl;
    file << "#,ambient p x^-1," << measurement->measurement_parameters->at("ambient p x^-1") << endl;
    file << "#,ambient p x^-2," << measurement->measurement_parameters->at("ambient p x^-2") << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#,6-way valve," << measurement->measurement_parameters->at("6-way valve") << endl;
    file << "#" << endl;
    file << "#,EXTERNAL SAMPLING CONTROL," << endl;
    file << "#" << endl;
    file << "#,control status," << measurement->measurement_parameters->at("control status") << endl;
    file << "#,control zero / signal units," << measurement->measurement_parameters->at("control zero / signal units") << endl;
    file << "#,control zero / V," << measurement->measurement_parameters->at("control zero / V") << endl;
    file << "#,control threshold / signal units," << measurement->measurement_parameters->at("control threshold / signal units") << endl;
    file << "#,control threshold / V," << measurement->measurement_parameters->at("control threshold / V") << endl;
    file << "#,control threshold2 / signal units," << measurement->measurement_parameters->at("control threshold2 / signal units") << endl;
    file << "#,control threshold2 / V," << measurement->measurement_parameters->at("control threshold2 / V") << endl;
    file << "#,control sampling time / s," << measurement->measurement_parameters->at("control sampling time / s") << endl;
    file << "#,control variable," << measurement->measurement_parameters->at("control variable") << endl;
    file << "#,control dimension," << measurement->measurement_parameters->at("control dimension") << endl;
    file << "#,control x^2," << measurement->measurement_parameters->at("control x^2") << endl;
    file << "#,control x^1," << measurement->measurement_parameters->at("control x^1") << endl;
    file << "#,control x^0," << measurement->measurement_parameters->at("control x^0") << endl;
    file << "#,control x^-1," << measurement->measurement_parameters->at("control x^-1") << endl;
    file << "#,control x^-2," << measurement->measurement_parameters->at("control x^-2") << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#,STATISTICS," << endl;
    file << "#" << endl;
    
    file << "#,RIP detection," << measurement->measurement_parameters->at("RIP detection") << endl;
    file << "#,tD  (RIP corr.) / ms," << measurement->measurement_parameters->at("tD  (RIP corr.) / ms") << endl;
    file << "#,1/K0 (RIP) / Vs/cm^2," << measurement->measurement_parameters->at("1/K0 (RIP) / Vs/cm^2") << endl;
    file << "#,K0 (RIP) / cm^2/Vs," << measurement->measurement_parameters->at("K0 (RIP) / cm^2/Vs") << endl;
    file << "#,SNR (RIP)," << measurement->measurement_parameters->at("SNR (RIP)") << endl;
    file << "#,WHM (RIP) / Vs/cm^2," << measurement->measurement_parameters->at("WHM (RIP) / Vs/cm^2") << endl;
    file << "#,res. power (RIP)," << measurement->measurement_parameters->at("res. power (RIP)") << endl;
    file << "#" << endl;
    file << "#,tD  (preRIP corr.) / ms," << measurement->measurement_parameters->at("tD  (preRIP corr.) / ms") << endl;
    file << "#,1/K0 (preRIP) / Vs/cm^2," << measurement->measurement_parameters->at("1/K0 (preRIP) / Vs/cm^2") << endl;
    file << "#,K0 (preRIP) / cm^2/Vs," << measurement->measurement_parameters->at("K0 (preRIP) / cm^2/Vs") << endl;
    file << "#,SNR (preRIP)," << measurement->measurement_parameters->at("SNR (preRIP)") << endl;
    file << "#,WHM (preRIP) / Vs/cm^2," << measurement->measurement_parameters->at("WHM (preRIP) / Vs/cm^2") << endl;
    file << "#,res. power (preRIP)," << measurement->measurement_parameters->at("res. power (preRIP)") << endl;
    file << "#" << endl;
    file << "#,signal RIP / V," << measurement->measurement_parameters->at("signal RIP / V") << endl;
    file << "#,signal preRIP / V," << measurement->measurement_parameters->at("signal preRIP / V") << endl;
    file << "#,RIP / preRIP," << measurement->measurement_parameters->at("RIP / preRIP") << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#,Fims / cm^2/kV," << measurement->measurement_parameters->at("Fims / cm^2/kV") << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "#" << endl;
    file << "\\   , tR";
    
    
    for_loop(i, measurement->data->n_rows) file << ", " << measurement->retention_times->at(i);
    file << endl;
    file << "1/K0, tDcorr.\\SNr";
    for_loop(i, measurement->data->n_rows) file << ", " << i;
    file << endl;
    
    
    
    for_loop(i, measurement->data->n_cols){
            file.precision(6);
            if (measurement->reduced_inversed_mobilities->at(i) > 1e-4){
                    file << measurement->reduced_inversed_mobilities->at(i);
            }
            else {
                    //file << fixed << measurement->reduced_inversed_mobilities->at(i);
                    file << measurement->reduced_inversed_mobilities->at(i);
            }
            file << ", " << measurement->drift_times->at(i);
            file.precision(4);
            
            for_loop(j, measurement->data->n_rows){
                    //file << fixed << ", " << -measurement->data->data->at(j)->at(i);
                    file << fixed << ", " << -measurement->data->data->at(j)->at(i);
            }
            file << endl;
    }
    file << endl;
    file.close();
    
}
  


/////////////////////////////////////////////////////////////////////////////////

/*!
 * \brief Parse MCC/IMS measurement in ims format from file
 * @param string filename
 * @return IMSMeasurement*
 */  
IMSMeasurement* IMSFileIMS::load(string filename){
    ifstream is;
    is.open (filename.c_str(), ios::binary);
    if (!is.good()){
        cout << "Error: " << filename << " could not be opened" << endl;
        return 0;
    }
    
    IMSMeasurement* measurement = new IMSMeasurement(filename);
    is.seekg (0, ios::end);
    long length = is.tellg();
    is.seekg (0, ios::beg);
    char *t = new char[length];
    is.read (t, length);    
    is.close();
    
    short ret, drift;
    uint R, T;
    unsigned long offset = 0;
    
    if (!(t[0] == 'i' && t[1] == 'm' && t[2] == 's')){
        delete []t;
        delete measurement;
        return 0;
    }
    offset += 3;
    
    // Alle String Metainformationen auslesen
    string str_names[] = {"data type", "version", "template version", "AD-board type", "ser.-no.", "date", "time", "file", "sample type", "sample ID", "comment", "location", "location name", "operator", "operator name", "IMS", "polarity", "drift gas", "sample gas", "carrier gas", "pre-separation type", "pre-separation T / deg C", "ambient T source", "ambient p source", "6-way valve", "control status", "control variable", "control dimension", "RIP detection"};
    for_loop (i, 29){
            short zahl = *( short *) (t + offset);
            offset += 2;
            string t_str = "";
            
            for_loop (j, (uint)zahl){
                    t_str += t[offset++];
            }
            measurement->input_parameter(str_names[i], t_str);
    }
    
    // Alle Float Metainformationen auslesen
    string flt_names[] = {"height ASL / m", "total data acquisition time / s",
    "K0 RIP positive / cm^2/Vs", "K0 RIP negative / cm^2/Vs", "grid opening time / us", "pause / s",
    "tD interval (corr.) / ms from", "tD interval (corr.) / ms to", "1/K0 interval / Vs/cm^2 from",
    "1/K0 interval / Vs/cm^2 to", "no. of data points per spectra", "no. of spectra", "no. averaged spectra",
    "baseline / signal units", "baseline / V", "V / signal unit", "drift length / mm", "HV / kV",
    "amplification / V/nA", "drift gas flow / mL/min", "sample gas flow / mL/min", "carrier gas flow / mL/min",
    "sample loop T / deg C", "sample loop volume / mL", "ambient T / deg C", "ambient T x^2", "ambient T x^1",
    "ambient T x^0", "ambient T x^-1", "ambient T x^-2", "ambient p / hPa", "ambient p x^2", "ambient p x^1",
    "ambient p x^0", "ambient p x^-1", "ambient p x^-2", "control zero / signal units", "control zero / V", 
    "control threshold / signal units", "control threshold / V", "control threshold2 / signal units", 
    "control threshold2 / V", "control sampling time / s", "control x^2", "control x^1", "control x^0", 
    "control x^-1", "control x^-2", "tD  (RIP corr.) / ms", "1/K0 (RIP) / Vs/cm^2", "K0 (RIP) / cm^2/Vs", 
    "SNR (RIP)", "WHM (RIP) / Vs/cm^2", "res. power (RIP)", "tD  (preRIP corr.) / ms", "1/K0 (preRIP) / Vs/cm^2", 
    "K0 (preRIP) / cm^2/Vs", "SNR (preRIP)", "WHM (preRIP) / Vs/cm^2", "res. power (preRIP)", "signal RIP / V", 
    "signal preRIP / V", "RIP / preRIP", "Fims / cm^2/kV"};
    for_loop (i, 64){
            
            stringstream number_string;
            number_string << *( float *) (t + offset);   
            offset += 4;
            string t_str = number_string.str();
            measurement->input_parameter(flt_names[i], t_str);
    }

    // Retentionszeit und Driftzeit auslesen
    ret = *( short *) (t + offset);
    measurement->data->n_rows = ret;
    offset += 2;
    drift = *( short *) (t + offset);
    measurement->data->n_cols = drift;
    offset += 2;    
    R = ret;
    T = drift;
    
    
    // write whether intensities are stored as 2 or 4 byte floats
    int bytesize = *( char *) (t + offset);
    if (!(bytesize == HALF || bytesize == FLOAT) ){
        delete []t;
        delete measurement;
        return 0;
    }
    offset++; 
    
                  
    // Read reduced inverse mobilities
    for_loop(i, T){
        measurement->reduced_inversed_mobilities->push_back( *( float *) (t + offset) );
        offset += 4;
    }
          
    // Read drift times
    for_loop (i, T){
        measurement->drift_times->push_back( *( float *) (t + offset) );
        offset += 4;
    }
    
    measurement->retention_times->resize(ret, 0.);
    measurement->data->data->resize(ret, 0);
    

    
    // Alle Spektren auslesen
    // erster Wert = Retentionszeit
    // alle weiteren Werte = Intensitäten
    #pragma omp parallel for
    for_loop (i, R){
        unsigned long offset2 = offset + i * (drift * bytesize + 4);
        measurement->retention_times->at(i) = *( float *) (t + offset2);
        measurement->data->data->at(i) = new ims_spectrum(drift, 0.);
        
        if (bytesize == HALF){
            for_loop (j, T){
                measurement->data->data->at(i)->at(j) = -(float)( *( half *) (t + (offset2 + 4 + j * bytesize) ) );
            } 
        }
        else {
            for_loop (j, T){
                measurement->data->data->at(i)->at(j) = -( *( float *) (t + (offset2 + 4 + j * bytesize) ) );
            } 
        }
        
         
    }
    delete []t;
    
    return measurement;
}

/*!
 * \brief Store MCC/IMS measurement in csv format in file ommiting bytesize of floats
 * @param IMSMeasurement* measurement
 * @param string filename
 */
void IMSFileIMS::save(IMSMeasurement* measurement, string filename){
    save(measurement, filename, FLOAT);
}

/*!
 * \brief Store MCC/IMS measurement in csv format in file
 * @param IMSMeasurement* measurement
 * @param string filename
 * @param int bytesize (for storing floats, 2 or 4 bytes can be chosen)
 */
void IMSFileIMS::save(IMSMeasurement* measurement, string filename, int bytesize){
    int num_chars = 0;
    string str_names[] = {"data type", "version", "template version", "AD-board type", "ser.-no.", "date", "time",
    "file", "sample type", "sample ID", "comment", "location", "location name", "operator", "operator name", "IMS",
    "polarity", "drift gas", "sample gas", "carrier gas", "pre-separation type", "pre-separation T / deg C",
    "ambient T source", "ambient p source", "6-way valve", "control status", "control variable", "control dimension",
    "RIP detection"};
    for_loop (i, 29){
        num_chars += measurement->measurement_parameters->at(str_names[i]).length();
    }
    
    unsigned long offset = 0;
    short ret = measurement->data->n_rows;
    short drift = measurement->data->n_cols;
    uint R = ret, T = drift;
    
    /* 318 = 64 * 4 + 29 * 2 + 3 + 1 */
    uint checksum = num_chars + 318 + ((measurement->retention_times->size() + measurement->reduced_inversed_mobilities->size() + measurement->drift_times->size()) * 4) + (ret * drift) * bytesize + 2 + 2;
    char *content = new char[checksum];
            
    // first three magic letters are "ims"
    content[0] = 'i';
    content[1] = 'm';
    content[2] = 's';
    
    offset += 3;
    
    // Alle String Values schreiben
    for_loop (i, 29){
        string t_str = measurement->measurement_parameters->at(str_names[i]);
        short length = t_str.length();
        memcpy(content + offset, &length, 2);
        offset += 2;
        memcpy(content + offset, t_str.c_str(), length);
        offset += length;
    }
    
    
    // Alle Float Values schreiben
    float *fp = (float*)(content + offset);
    string flt_names[] = {"height ASL / m", "total data acquisition time / s", "K0 RIP positive / cm^2/Vs",
    "K0 RIP negative / cm^2/Vs", "grid opening time / us", "pause / s", "tD interval (corr.) / ms from",
    "tD interval (corr.) / ms to", "1/K0 interval / Vs/cm^2 from", "1/K0 interval / Vs/cm^2 to",
    "no. of data points per spectra", "no. of spectra", "no. averaged spectra", "baseline / signal units",
    "baseline / V", "V / signal unit", "drift length / mm", "HV / kV", "amplification / V/nA",
    "drift gas flow / mL/min", "sample gas flow / mL/min", "carrier gas flow / mL/min", "sample loop T / deg C",
    "sample loop volume / mL", "ambient T / deg C", "ambient T x^2", "ambient T x^1", "ambient T x^0",
    "ambient T x^-1", "ambient T x^-2", "ambient p / hPa", "ambient p x^2", "ambient p x^1", "ambient p x^0",
    "ambient p x^-1", "ambient p x^-2", "control zero / signal units", "control zero / V",
    "control threshold / signal units", "control threshold / V", "control threshold2 / signal units",
    "control threshold2 / V", "control sampling time / s", "control x^2", "control x^1", "control x^0",
    "control x^-1", "control x^-2", "tD  (RIP corr.) / ms", "1/K0 (RIP) / Vs/cm^2", "K0 (RIP) / cm^2/Vs",
    "SNR (RIP)", "WHM (RIP) / Vs/cm^2", "res. power (RIP)", "tD  (preRIP corr.) / ms", "1/K0 (preRIP) / Vs/cm^2",
    "K0 (preRIP) / cm^2/Vs", "SNR (preRIP)", "WHM (preRIP) / Vs/cm^2", "res. power (preRIP)", "signal RIP / V",
    "signal preRIP / V", "RIP / preRIP", "Fims / cm^2/kV"};
    for_loop (i, 64){
        fp[i] = atof(measurement->measurement_parameters->at(flt_names[i]).c_str());
    }
    offset += 256;
    
    *((short*)(content + offset)) = ret;
    offset += 2;
    *((short*)(content + offset)) = drift;
    offset += 2;
    
    
    
    // write whether intensities are stored as 2 or 4 byte floats
    *((char*)(content + offset)) = char(bytesize);
    offset++;
    
    
    // write RIM scale
    fp = (float*)(content + offset);
    for_loop (i, T){
        fp[i] = measurement->reduced_inversed_mobilities->at(i); 
    }
    offset += T << 2;
    
    // write drift times scale
    fp = (float*)(content + offset);
    for_loop (i, T){
        fp[i] = measurement->drift_times->at(i);
    }
    offset += T << 2;
    ulong o3 = offset;
    
    
    #pragma omp parallel for
    for_loop (i, R){
        // Scala Correction schreiben
        long offset2 = o3 + (4 + drift * bytesize) * i;
        *((float*)(content + offset2)) = measurement->retention_times->at(i);
        offset2 += 4;
        
        // Signale in Array schreien
        if (bytesize == HALF){
            half *h = (half*)(content + offset2);
            for_loop (j, T){
                h[j] = (half)-measurement->data->data->at(i)->at(j);
            }
        }
        else {
            float *f = (float*)(content + offset2);
            for_loop (j, T){
                f[j] = -measurement->data->data->at(i)->at(j);
            }
        }
    }
    offset += (4 + drift * bytesize) * R;

    assert(checksum == offset);
    
    // write out content
    ofstream wrt(filename.c_str(), ios::binary);
    wrt.write(content, offset);
    wrt.close();
    
    delete []content;
}
/////////////////////////////////////////////////////////////////////////////////

#ifdef usexls
/*!
 * \brief Parse peak list in xls format from file
 * @param IMSMeasurement* measurement
 * @param string filename
 * @param string layer_name
 */
IMSPeakList* IMSLayerXLS::load(string filename){
    IMSPeakList* peaklist = new IMSPeakList;
    peaklist->file_name = filename;
    xlsWorkBook* workbook = xls_open((char*) filename.c_str() , "UTF-8"); 
    xlsWorkSheet* worksheet = xls_getWorkSheet(workbook, 0); // first sheet;
    xls_parseWorkSheet(worksheet);
    
    WORD cellRow;
    vector<IMSPeak*>* loaded_peaklist = new vector<IMSPeak*>;
    
    // Rows 0-4: Meta-information
    // Row 4: first peak
    for (cellRow = 4; cellRow <= worksheet->rows.lastrow; cellRow++) {
        IMSPeak* current_peak = new IMSPeak;
        map<string, double>* current_peakparameters = new map<string, double>;
        xlsCell *cell_name = xls_cell(worksheet, cellRow, 0);
        xlsCell *cell_t = xls_cell(worksheet, cellRow, 2);
        xlsCell *cell_r = xls_cell(worksheet, cellRow, 3);
        xlsCell *cell_tol_t = xls_cell(worksheet, cellRow, 4);
        xlsCell *cell_tol_r = xls_cell(worksheet, cellRow, 5);
        
        current_peak->peak_name = ((char *)cell_name->str);
        current_peak->t = cell_t->d;
        current_peak->r = cell_r->d;
        current_peakparameters->insert(pair<string, double>("r_tol", cell_tol_r->d));
        current_peakparameters->insert(pair<string, double>("t_tol", cell_tol_t->d));
        current_peak->peak_parameters = current_peakparameters; 
        
        loaded_peaklist->push_back(current_peak);
    }
    xls_close(workbook); // close the sheet	
    peaklist->ims_peak_list = loaded_peaklist;
    return peaklist;
}



/*!
 * \brief Store peak list in xls format in file
 * @param IMSMeasurement* measurement
 * @param string filename
 * @param string layer_name
 */
void IMSLayerXLS::save(IMSPeakList* peaklist, string filename, string layer_name){
    workbook* wb = new workbook();     
    worksheet* sheet = wb->sheet(layer_name.c_str());
    
    sheet->label(0, 0, "ID");
    sheet->label(0, 1, "Type");
    sheet->label(0, 2, "default Color");
    sheet->label(0, 3, "Properties");
        
    sheet->label(1, 0, "IMSLayer.xls");
    sheet->label(1, 1, "ovaldot");
    sheet->number(1, 2, -1.7e7);
    sheet->label(1, 3, "aced000570");
        
    sheet->label(3, 0, "Name");
    sheet->label(3, 1, "Comment");
    sheet->label(3, 2, "1/K0");
    sheet->label(3, 3, "RT");
    sheet->label(3, 4, "1/K0 radius");
    sheet->label(3, 5, "RT radius");
    sheet->label(3, 6, "Color");
    int g = 0;
    for (size_t i = 0; i < peaklist->ims_peak_list->size(); i++){            
        sheet->label(4 + g, 0, peaklist->ims_peak_list->at(i)->peak_name);
        sheet->label(4 + g, 1, peaklist->ims_peak_list->at(i)->peak_name); 
        sheet->number(4 + g, 2, peaklist->ims_peak_list->at(i)->t);
        sheet->number(4 + g, 3, peaklist->ims_peak_list->at(i)->r);
        sheet->number(4 + g, 4, peaklist->ims_peak_list->at(i)->peak_parameters->at("t_tol"));
        sheet->number(4 + g, 5, peaklist->ims_peak_list->at(i)->peak_parameters->at("r_tol"));
        sheet->number(4 + g, 6, -16777216);
        g++;
    }
    wb->Dump(filename.c_str());
    delete wb;
}
#endif

/*!
 * \brief 
 * @param 
 * @return 
 */
IMSPeakList* IMSLayerCSV::load(string filename){
    IMSPeakList* peaklist = new IMSPeakList;
    cout << "Not implemented yet!" << endl;
    return peaklist;
}

/*!
 * \brief 
 * @param 
 * @return 
 */
void IMSLayerCSV::save(IMSPeakList* peaklist, string filename, string layer_name){
    cout << "Not implemented yet!" << endl;
}

/*!
 * \brief 
 * @param 
 * @return 
 */
pmap* read_in_config(string path){
    pmap* pipeline_parameters = new pmap;
    
    fstream conf_file( path.c_str(), ios::in );
    if (!conf_file.is_open()){
        cout << "Error: " << path << " could not be opened" << endl;
        return 0;
    }
    string line;
    while ( !conf_file.eof() ){
        getline (conf_file,line);
        if (line[0] == '#') continue;
        if (line.find('\t') == line.npos) continue;
        stringlist parameters = split(line, '\t');
        pipeline_parameters->parameter_insert(parameters[0], parameters[1]);
    }
    return pipeline_parameters;
}

/*!
 * \brief 
 * @param 
 * @return 
 */
string strip(string s, char strip){
    while (s.length() > 0 && s.at(0) == strip) s = s.substr(1, s.length() - 1);
    while (s.length() > 0 && s.at(s.length() - 1) == strip) s = s.substr(0, s.length() - 1);
    return s;
}

stringlist &split(const string &s, char delim, stringlist &elems) {
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(strip(item));
    }
    return elems;
}

/*!
 * \brief String tokenizer
 * @param const string &s
 * @param char delim
 * @return stringlist
 */
stringlist split(const string &s, char delim) {
    stringlist elems;
    return split(s, delim, elems);
}

/*!
 * \brief Computes mode position for inverse Gaussian distribution
 * @param double mue
 * @param double lambda
 * @param double offset
 * @return 
 */
inline double get_mode(double mue, double lambda, double offset){
    return mue * (sqrt(1.0 + square(3.0 * mue) / square(2. * lambda)) - (3.0 * mue) / (2. * lambda)) + offset;
}


/*!
 * \brief Computes inverse Gaussian descriptors from inverse Gaussian parameters
 * @param double mue
 * @param double lambda
 * @param double offset
 * @param double* e
 * @param double* sigma
 * @param double* mode
 */
void getDescriptors(double mue, double lambda, double offset, double* e, double* sigma, double* mode){
    *(e) = mue + offset;
    *(sigma) = sqrt(mue * mue * mue / lambda);
    *(mode) = get_mode(mue, lambda, offset);
}

/*!
 * \brief Computes inverse Gaussian parameters from inverse Gaussian descriptors
 * @param double e
 * @param double sigma
 * @param double mode
 * @param double* mue
 * @param double* lambda
 * @param double* offset
 */
void getModelParams(double e, double sigma, double mode, double* mue, double* lambda, double* offset){
    double p = 0.5 * ((-mode * (2. * e + mode) + 3.0 * (square(e) - square(sigma))) / (2. * (mode - e)));
    double q = (mode * (3.0 * square(sigma) + e * mode) - (e * e * e)) / (2. * (mode - e));

    *(offset) = - p - sqrt(square(p) - q);
    *(mue) = e - *(offset);
    *(lambda) = square(*(mue)) * (*(mue)) / square(sigma);
}

/*!
 * \brief Checks if string fullString has particular ending
 * @param string const &fullString
 * @param string const &ending
 * @return bool
 */
bool ends_with (string const &fullString, string const &ending){
    if (fullString.length() >= ending.length()) {
        return (!fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    }
    else {
        return false;
    }
}
