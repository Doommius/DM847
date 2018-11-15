/*!
 * @author Marianna D'Addario, Dominik Kopczynski
 * @e-mail {Marianna.Daddario, Dominik.Kopczynski}@tu-dortmund.de
 * Copyright (c) 2013 Marianna D'Addario, Dominik Kopczynski
 *
 * LICENSE:
 * PEAX is a free software for academic use only.
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "ims_core.h"
#include "preprocessing.h"
#include "candidate-detection.h"
#include "picking.h"
#include "modeling.h"


class Configuration {
 public:
   string parameters_file;
   string intensity_threshold;

   Configuration() :
     parameters_file(""), intensity_threshold( "" )
     {}
 };


int main(int argc, char* argv[]){
    //Parse argv
    string usage = "Usage: source_measurement target_peak_list -i<int: value for intensity threshold> -p<string: name of your parameters.cfg file> \n Note the parameter i always overrules parameter p.";
    if( argc < 3 ){
        cerr << usage << endl;
        exit(1);
    }

    /* Initialize global configuarations */
    Configuration config;
    
    string source_dir = argv[1];
    string target_dir = argv[2];
    
    int opt = 0;
    while( (opt = getopt (argc, argv, "p:i:")) != -1 ) {
        switch( opt ){
            case 'i':
                config.intensity_threshold = optarg;
                break;

            case 'p':
                config.parameters_file = optarg;
                break;

            case '?':
                cout << usage << endl;
                break;

            default:
                cout << usage << endl;
                break;
        }
		
    }

    pmap* pipeline_parameters;
    if( config.parameters_file != ""){
        pipeline_parameters = read_in_config(config.parameters_file);
    }else{
        pipeline_parameters = read_in_config("parameters.cfg");
    }
    if (!pipeline_parameters) exit(0);

    if( config.intensity_threshold != "" ){
        pipeline_parameters->erase("intensity_threshold");
        pipeline_parameters->insert(pair<string, string>("intensity_threshold", config.intensity_threshold));
    }  
    
    
    // check whether all necessary parameters are set
    if (pipeline_parameters->find("tol_rt_percent") == pipeline_parameters->end()){
        cout << "Error, 'tol_rt_percent' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("preprocessing_first") == pipeline_parameters->end()){
        cout << "Error, 'preprocessing_first' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("preprocessing_second") == pipeline_parameters->end()){
        cout << "Error, 'preprocessing_second' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("preprocessing_third") == pipeline_parameters->end()){
        cout << "Error, 'preprocessing_third' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("candidate_detection") == pipeline_parameters->end()){
        cout << "Error, 'candidate_detection' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("picking") == pipeline_parameters->end()){
        cout << "Error, 'picking' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("modeling") == pipeline_parameters->end()){
        cout << "Error, 'modeling' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("tol_rt") == pipeline_parameters->end()){
        cout << "Error, 'tol_rt' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("tol_rim") == pipeline_parameters->end()){
        cout << "Error, 'tol_rim' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("intensity_threshold") == pipeline_parameters->end()){
        cout << "Error, 'intensity_threshold' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("area_size") == pipeline_parameters->end()){
        cout << "Error, 'area_size' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("fftcutoff") == pipeline_parameters->end()){
        cout << "Error, 'fftcutoff' parameter is not set. exit" << endl;
        exit(0);
    }
    if (pipeline_parameters->find("expansion_size") == pipeline_parameters->end()){
        cout << "Error, 'expansion_size' parameter is not set. exit" << endl;
        exit(0);
    }
    
    

	
    // Initializing all function pointers    
    // Preprocessing: Baseline Correction, mixed_smoothing and de_nosing for permutation
    preprocessing_p pre[3];
    if (pipeline_parameters->at("preprocessing_first") == "dn"){
        pre[0] = de_noising;
    }
    else if (pipeline_parameters->at("preprocessing_first") == "bc"){
        pre[0] = baseline_correction;
    }
    else if (pipeline_parameters->at("preprocessing_first") == "s"){
        pre[0] = mixed_smoothing;
    }
    else {
        cout << "Error, unknown preprocessing modul: '" << pipeline_parameters->at("preprocessing_first") << "'. exit" << endl;
        exit(0);
    }
    
    if (pipeline_parameters->at("preprocessing_second") == "dn"){
        pre[1] = de_noising;
    }
    else if (pipeline_parameters->at("preprocessing_second") == "bc"){
        pre[1] = baseline_correction;
    }
    else if (pipeline_parameters->at("preprocessing_second") == "s"){
        pre[1] = mixed_smoothing;
    }
    else {
        cout << "Error, unknown preprocessing modul: '" << pipeline_parameters->at("preprocessing_second") << "'. exit" << endl;
        exit(0);
    }
    
    if (pipeline_parameters->at("preprocessing_third") == "dn"){
        pre[2] = de_noising;
    }
    else if (pipeline_parameters->at("preprocessing_third") == "bc"){
        pre[2] = baseline_correction;
    }
    else if (pipeline_parameters->at("preprocessing_third") == "s"){
        pre[2] = mixed_smoothing;
    }
    else {
        cout << "Error, unknown preprocessing modul: '" << pipeline_parameters->at("preprocessing_third") << "'. exit" << endl;
        exit(0);
    }
    
    if ((pre[0] == pre[1]) || (pre[0] == pre[2]) || (pre[1] == pre[2])){
        cout << "Error, preprocessing moduls multiple defined. exit" << endl;
        exit(0);
    }

    // Peak candidate detection - local maxima, crossfinding
    candidate_detection_p candidate;
    if (pipeline_parameters->at("candidate_detection") == "lm"){
        candidate = local_maxima;
    }
    else if (pipeline_parameters->at("candidate_detection") == "cf"){
        candidate = cross_finding;
    }
    else {
        cout << "Error, unknown candidate detection modul: '" << pipeline_parameters->at("candidate_detection") << "'. exit" << endl;
        exit(0);
    }
    
    // Peak picking - merging by signal, peace (cluster editing), em-clustering 
    picking_p pick;
    if (pipeline_parameters->at("picking") == "ms"){
        pick = merging_by_signal;
    }
    else if (pipeline_parameters->at("picking") == "ce"){
        pick = cluster_editing;
    }
    else if (pipeline_parameters->at("picking") == "emc"){
        pick = em_clustering;
    }
    else {
        cout << "Error, unknown picking modul: '" << pipeline_parameters->at("picking") << "'. exit" << endl;
        exit(0);
    }

    
    // Peak model estimation
    modeling_p model;
    if (pipeline_parameters->at("modeling") == "pme"){
        model = pme_modeling;
    }
    else if (pipeline_parameters->at("modeling") == "e"){
        model = no_modeling;
    }
    else {
        cout << "Error, unknown modeling modul: '" << pipeline_parameters->at("modeling") << "'. exit" << endl;
        exit(0);
    }
    
    
    // loading the measurement
    IMSMeasurement* measurement = IMSFileCSV::load(argv[1]);
    if (!measurement){
        measurement = IMSFileIMS::load(argv[1]);
        if (!measurement){
            cout << "Error, source measurement could not be opened. exit" << endl;
            exit(0);
        }
    }
    cout << "Measurement file opened: " << argv[1] << endl << endl;
    
    
    // doing the pipeline
    pre[0](measurement, pipeline_parameters);
    cout << "First preprocessing step: "<< pipeline_parameters->at("preprocessing_first") << endl<< endl;
    pre[1](measurement, pipeline_parameters);
    cout << "Second preprocessing step: "<< pipeline_parameters->at("preprocessing_second") << endl<< endl;
    pre[2](measurement, pipeline_parameters);
    cout << "Third preprocessing step: "<< pipeline_parameters->at("preprocessing_third") <<endl<< endl;
    just_positives(measurement, pipeline_parameters);
    IMSPeakList* list1 = candidate(measurement, pipeline_parameters);
    cout << "Candidate detection module " << pipeline_parameters->at("candidate_detection") << " found " << list1->ims_peak_list->size() << " peaks." << endl<< endl;
    IMSPeakList* list2 = pick(list1, pipeline_parameters);
    cout << "Picking module " << pipeline_parameters->at("picking") << " found " << list2->ims_peak_list->size() << " peaks remained." << endl<< endl;
    IMSPeakList* list3 = model(measurement, list2, pipeline_parameters);
    cout << "Modeling module " << pipeline_parameters->at("modeling") << " chosen."  << endl<< endl;
    list3->save(argv[2]);
    cout << "Peak list saved: "<< argv[2] << endl << endl;
    
    
    // cleaning the main memory
    delete list3;
    delete list2;
    delete list1;
    delete measurement;
    delete pipeline_parameters;
    
}
