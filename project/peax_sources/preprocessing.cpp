#include "preprocessing.h"

/*!
 * \brief Corrects the baseline of MCC/IMS measurement using median correction of chromatograms
 * @param matrix instance of IMSMeasurement class
 */


preprocessing(baseline_correction){
    IMSMatrix* matrix = measurement->data;
    uint R = matrix->n_rows;
    uint T = matrix->n_cols;
    
    
    
    for_loop (t, T){
        // correction over all contour lines
        int maxheight = 0;
        int minheight = 0;
        for_loop(r, R){
            maxheight = max(double(maxheight), matrix->data->at(r)->at(t));
            minheight = min(double(minheight), matrix->data->at(r)->at(t));
        }
        
        
        if (maxheight++){
            uint *hist = new uint[maxheight - minheight];
            double* member = new double[maxheight - minheight];
            
            for_loop (i, maxheight - minheight){
                hist[i] = 0;
            }
            for_loop (r, R){
                hist[int(matrix->data->at(r)->at(t)) - minheight]++;
            }
            int index = 0;
            for_loop (i, maxheight - minheight){
                if (hist[index] < hist[i]){
                    index = i;
                }
            }
            
  
            
            
            double mue = index, sigma = 2, tmp = 0;
            double old_mue, old_sigma, old_weight;
            bool run = true;
            int max_em_loops = 15, loops = 0;
            double weight = 0;
            for_loop(i, 5){
                if (0 <= index + i - 2 && index + i - 2 < maxheight - minheight){ 
                    weight += hist[index + i - 2];
                }
            }
            
            
            
            weight = min(weight / double(R), 0.999);
            
            double new_mue;
            while(run && loops++ < max_em_loops){
                old_mue = mue;
                old_sigma = sigma;
                old_weight = weight;
                new_mue = 0;
        
                // Expectation
                tmp = 0;
                #pragma omp parallel for reduction(+:tmp) reduction(+:new_mue)
                for_loop(r, maxheight - minheight){
                    double pdf = weight * g(r, mue, sigma);
                    double mod_pdf = (1. - weight) / double(maxheight - minheight);
                    member[r] = pdf / (pdf + mod_pdf);
                    assert(-0. <= member[r]);
                    assert(member[r] <= 1.);
                    tmp += member[r] * hist[r];
                    new_mue += member[r] * r * hist[r];
                }
                
        
                // Maximization
                mue = new_mue / tmp;
                weight = tmp / double(R);
        
                sigma = 0;
                #pragma omp parallel for reduction(+:sigma)
                for_loop(r, maxheight - minheight){
                    sigma += member[r] * hist[r] * square(mue - r);
                }
                sigma = max(sqrt(sigma / tmp), 1e-6);
        
                assert(!isnan(sigma));
                assert(!isinf(sigma)); 
                assert(sigma > 0.);
                assert(!isnan(mue));
                assert(!isinf(mue));
                assert(!isnan(weight));
                assert(!isinf(weight));
                assert(weight >= 0.);
                assert(weight <= 1.);
                
               
                run = false;
                if (!convergence(mue, old_mue, 1e-3)){
                    run = true;
                }
                if (!convergence(sigma, old_sigma, 1e-3)){
                    run = true;
                }
                if (!convergence(weight, old_weight, 1e-3)){
                    run = true;
                }
            }
            
            double lower = mue + minheight + 2. * sigma;
            
            
            
            //#pragma omp parallel for
            for_loop(r, R){
                double old = matrix->data->at(r)->at(t);
                if (old > 0){
                    matrix->data->at(r)->at(t) = (old - lower >= 0) ? old - lower : 0;
                }
            }
    
            delete []member;
            delete []hist;
        }
    }
}






/*
preprocessing(baseline_correction){
    IMSMatrix* matrix = measurement->data;
    uint R = matrix->n_rows;
    uint T = matrix->n_cols;
    
    double** A = new double*[R];
    for_loop(r, R){
        A[r] = new double[T];
    }
    double* member = new double[R];
    
    #pragma omp parallel for collapse(2)
    for_loop(r, R){
        for_loop(t, T){
            double A_z = 0, A_n = 0;
            for (int i = -4; i <= 4; i++){
                for (int j = -4; j <= 4; j++){
                    if (0 <= r + i && r + i < R && 0 <= t + j && t + j < T){
                        A_z += matrix->data->at(r + i)->at(t + j);
                        A_n++;
                    }
                }
            }
            A[r][t] = A_z / A_n;
        }
    }
    for_loop (t, T){
        // correction over all contour lines
        int maxheight = 0;
        int minheight = 0;
        for_loop(r, R){
            maxheight = max(double(maxheight), matrix->data->at(r)->at(t));
            minheight = min(double(minheight), matrix->data->at(r)->at(t));
        }
        
        
        if (maxheight++){
            uint *hist = new uint[maxheight - minheight];
            //double* member = new double[maxheight - minheight];
            
            for_loop (i, maxheight - minheight){
                hist[i] = 0;
            }
            for_loop (r, R){
                hist[int(matrix->data->at(r)->at(t)) - minheight]++;
            }
            int index = 0;
            for_loop (i, maxheight - minheight){
                if (hist[index] < hist[i]){
                    index = i;
                }
            }
            
  
            
            
            double mue = index + minheight, sigma = 2, tmp = 0;
            double old_mue, old_sigma, old_weight;
            bool run = true;
            int max_em_loops = 15, loops = 0;
            double weight = 0;
            for_loop(i, 5){
                if (0 <= index + i - 2 && index + i - 2 < maxheight - minheight){ 
                    weight += hist[index + i - 2];
                }
            }
            
            
            
            weight = min(weight / double(R), 0.999);
            
            double new_mue;
            while(run && loops++ < max_em_loops){
                old_mue = mue;
                old_sigma = sigma;
                old_weight = weight;
                new_mue = 0;
        
                // Expectation
                tmp = 0;
                #pragma omp parallel for reduction(+:tmp) reduction(+:new_mue)
                for_loop(r, R){
                    double pdf = weight * g(A[r][t], mue, sigma);
                    double mod_pdf = (1. - weight) / double(maxheight - minheight);
                    member[r] = pdf / (pdf + mod_pdf);
                    assert(-0. <= member[r]);
                    assert(member[r] <= 1.);
                    tmp += member[r];
                    new_mue += member[r] * A[r][t];
                }
                
        
                // Maximization
                mue = new_mue / tmp;
                weight = tmp / double(R);
        
                sigma = 0;
                #pragma omp parallel for reduction(+:sigma)
                for_loop(r, R){
                    sigma += member[r] * square(mue - A[r][t]);
                }
                sigma = max(sqrt(sigma / tmp), 1e-6);
        
                assert(!isnan(sigma));
                assert(!isinf(sigma)); 
                assert(sigma > 0.);
                assert(!isnan(mue));
                assert(!isinf(mue));
                assert(!isnan(weight));
                assert(!isinf(weight));
                assert(weight >= 0.);
                assert(weight <= 1.);
                
               
                run = false;
                if (!convergence(mue, old_mue, 1e-3)){
                    run = true;
                }
                if (!convergence(sigma, old_sigma, 1e-3)){
                    run = true;
                }
                if (!convergence(weight, old_weight, 1e-3)){
                    run = true;
                }
            }
            
            #pragma omp parallel for
            for_loop(r, R){
                matrix->data->at(r)->at(t) = (matrix->data->at(r)->at(t) - mue) * (1. - member[r]) + mue;
            }
            
            #pragma omp parallel for
            for_loop(r, R){
                double old = matrix->data->at(r)->at(t);
                if (old > 0){
                    matrix->data->at(r)->at(t) = (old - mue >= 0) ? old - mue : 0;
                }
            }
    
            delete []hist;
        }
    }
    
    for_loop(r, R){
        delete []A[r];
    }
    delete []A;
    delete []member;
}
*/

/*!
 * \brief Detailing MCC/IMS measurement using device function
 * @param matrix instance of IMSMeasurement class
 */
preprocessing(de_tailing){
    // not yet implemented
}

/*!
 * \brief Crops matrix by given cropping point
 * @param matrix instance of IMSMeasurement class
 */
preprocessing(crop){
    IMSMatrix* matrix = measurement->data;
    float R = measurement->retention_times->back();
    float T = measurement->reduced_inversed_mobilities->back();
    
    float T_up = T, R_up = R;
    float T_low = 0, R_low = 0;
    
    if (pipeline_parameters->find("crop_T_high") != pipeline_parameters->end()){
        T_up = atof(pipeline_parameters->at("crop_T_high").c_str());
    }
    if (pipeline_parameters->find("crop_T_low") != pipeline_parameters->end()){
        T_low = atof(pipeline_parameters->at("crop_T_low").c_str());
    }
    if (pipeline_parameters->find("crop_R_high") != pipeline_parameters->end()){
        R_up = atof(pipeline_parameters->at("crop_R_high").c_str());
    }
    if (pipeline_parameters->find("crop_R_low") != pipeline_parameters->end()){
        R_low = atof(pipeline_parameters->at("crop_R_low").c_str());
    }
    
    
    while (measurement->retention_times->back() > R_up && matrix->data->size()) {
        delete matrix->data->back();
        matrix->data->pop_back();
        measurement->retention_times->pop_back();
    }
    
    while (measurement->retention_times->front() < R_low && matrix->data->size()) {
        delete matrix->data->front();
        matrix->data->erase(matrix->data->begin());
        measurement->retention_times->erase(measurement->retention_times->begin());
    }
    
    
    while (measurement->reduced_inversed_mobilities->back() > T_up && matrix->data->front()->size()) {
        for_loop(r, matrix->data->size()){
            matrix->data->at(r)->pop_back();
        }
        measurement->reduced_inversed_mobilities->pop_back();
        measurement->drift_times->pop_back();
    }
        
    uint t_index = 0;
    while (measurement->reduced_inversed_mobilities->size() && measurement->reduced_inversed_mobilities->front() < T_low){
        t_index++;
        measurement->reduced_inversed_mobilities->erase(measurement->reduced_inversed_mobilities->begin());
        measurement->drift_times->erase(measurement->drift_times->begin());
    }
    
    for_loop(r, matrix->data->size()){
        for (uint t = t_index; t < matrix->data->at(r)->size(); t++){
            swap(matrix->data->at(r)->at(t), matrix->data->at(r)->at(t - t_index));
        }
        for_loop (t, t_index){
            matrix->data->at(r)->pop_back();
        }
    }
    matrix->n_rows = matrix->data->size();
    matrix->n_cols = matrix->data->at(0)->size();
}


/*!
 * \brief Baseline correction MCC/IMS measurement by cutting all values below 0 to 0
 * @param matrix instance of IMSMeasurement class
 */
preprocessing(just_positives){
    IMSMatrix* matrix = measurement->data;
    uint R = matrix->n_rows;
    uint T = matrix->n_cols;
//     #pragma omp parallel for collapse(2)
    for_loop(r, R){
        for_loop(t, T){
            if (matrix->data->at(r)->at(t) < 0){
                matrix->data->at(r)->at(t) = 0;
            }
        }
    }
}

/*!
 * \brief Dummy empty function for no preprocessing
 * @param matrix instance of IMSMeasurement class
 */
preprocessing(no_preprocessing){

}


/*!
 * \brief Helper function for mixed_smoothing
 * @param c_points Matrix of complex points of original matrix
 * @param dir Direction of Fourier Transform
 */
void FFT2D(complex** c_points, uint R, uint T, int dir){
    CFFT fft;
    
    #pragma omp parallel for
    for_loop (r, R) {
        dir ? fft.Forward(c_points[r], T) : fft.Inverse(c_points[r], T);
    }
    
    #pragma omp parallel for
    for_loop (t, T) {
        complex* c2 = new complex[R];
        for_loop (r, R) {
            c2[r] = c_points[r][t];
        }
        dir ? fft.Forward(c2, R) : fft.Inverse(c2, R);
        for_loop (r, R) {
            c_points[r][t] = c2[r];
        }
        delete []c2;
    }
}



/*!
 * \brief Several methods to smooth the MCC/IMS measurement and writes a preprocessed MCC/IMS measurement
 * @param matrix instance of IMSMeasurement class
 */
preprocessing(mixed_smoothing){
    assert(pipeline_parameters->find("fftcutoff") != pipeline_parameters->end());
    
    IMSMatrix* matrix = measurement->data;
    uint fftcutoff = atoi(pipeline_parameters->at("fftcutoff").c_str());
    
    uint R = matrix->n_rows;
    uint T = matrix->n_cols;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // low-pass filter
    
    
    
    uint T_up = 1 << (int)ceil(log2(T));
    uint R_up = 1 << (int)ceil(log2(R));
    uint T_up_half = T_up >> 1;
    uint R_up_half = R_up >> 1;

    complex** c_points = new complex*[R_up];
    
    // setting up complex matrix for FFT
    #pragma omp parallel for
    for_loop (r, R_up){
        complex* tmp = new complex[T_up];
        for_loop (t, T_up){
            if (r < R && t < T){
                tmp[t] = matrix->data->at(r)->at(t);
            }
        }
        c_points[r] = tmp;
    }
    
    
    // FFT in forward direction
    FFT2D(c_points, R_up, T_up, FFTFORWARD);
    
    
    // FFT Shift
    #pragma omp parallel for
    for_loop (r, R_up_half){
        swap(c_points[r], c_points[r + R_up_half]);
    }
    #pragma omp parallel for collapse(2)
    for_loop (r, R_up){
        for_loop (t, T_up_half){
            swap(c_points[r][t], c_points[r][t + T_up_half]);
        }
    }


    // cutting all values exceeding threshold
    complex n(0, 0);
    #pragma omp parallel for collapse(2)
    for_loop (r, R){
        for_loop (t, T){
            double val = sqrt((R_up_half - r) * (R_up_half - r) + (T_up_half - t) * (T_up_half - t));
            
            if (val > fftcutoff) {
                c_points[r][t] *= n;
            }
        }
    }


    // FFT Shift
    #pragma omp parallel for
    for_loop (r, R_up_half){
        swap(c_points[r], c_points[r + R_up_half]);
    }
    #pragma omp parallel for collapse(2)
    for_loop (r, R_up){
        for_loop (t, T_up_half){
            swap(c_points[r][t], c_points[r][t + T_up_half]);
        }
    }
    
    
    // FFT in backward direction
    FFT2D(c_points, R_up, T_up, FFTBACKWARD);
    
    
    // writing back smoothed points in matrix
    for_loop (r, R){
        for_loop (t, T){
            matrix->data->at(r)->at(t) = c_points[r][t].re();
        }
    }

    
    // deleting auxiliary data structures
    for_loop (r, R_up){
        delete []c_points[r];
    }
    delete []c_points;
    

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Gaussian blur
    /*
    
    double** gaussbell = new double*[5];
    for_loop(i, 5) gaussbell[i] = new double[5];
    gaussbell[0][0] = gaussbell[0][4] = gaussbell[4][0] = gaussbell[4][4] = 1.;
    gaussbell[0][1] = gaussbell[0][3] = gaussbell[1][0] = gaussbell[1][4] = 4.;
    gaussbell[3][0] = gaussbell[3][4] = gaussbell[4][1] = gaussbell[4][3] = 4.;
    gaussbell[0][2] = gaussbell[2][0] = gaussbell[2][4] = gaussbell[4][2] = 7.;
    gaussbell[1][1] = gaussbell[1][3] = gaussbell[3][1] = gaussbell[3][3] = 20.;
    gaussbell[1][2] = gaussbell[2][1] = gaussbell[2][3] = gaussbell[3][2] = 33.;
    gaussbell[2][2] = 55.;
    
    
    double signal_sum = 0, blur_sum = 0;
    for_loop (r, R){
        for_loop (t, T){
            signal_sum += matrix->data->at(r)->at(t);
        }
    }
    
    // setting up matrix for blured values
    matrix_t *blured = new matrix_t;
    for_loop (r, R){
        blured->push_back(new ims_spectrum);
        for_loop (t, T){
            blured->at(r)->push_back(0);
        }
    }


    
    // computing all blured values by summing all weighted values within an area
    #pragma omp parallel for collapse(2) reduction(+:blur_sum)
    for_loop (r, R){
        for_loop (t, T){
            double gauss_val = 0, weight = 0;
            for (int m = -2; m <= 2; m++){
                for (int n = -2; n <= 2; n++){
                    if (!(((r + m) < 0) || ((t + n) < 0) || ((r + m) >= R) || ((t + n) >= T))){
                        gauss_val += gaussbell[2 + m][2 + n] * matrix->data->at(r + m)->at(t + n);
                        weight += gaussbell[2 + m][2 + n];
                    }
                }
            }
            gauss_val /= weight;
            blured->at(r)->at(t) = gauss_val;
            blur_sum += blured->at(r)->at(t);
        }
    }

    // normalize
    #pragma omp parallel for collapse(2)
    for_loop (r, R){
        for_loop (t, T){
            matrix->data->at(r)->at(t) = (blured->at(r)->at(t) * signal_sum / blur_sum);
        }
    }
    
    for_loop (r, R) delete blured->at(r);
    delete blured;
    
    for_loop (i, 5) delete []gaussbell[i];
    delete []gaussbell;
    
    */
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Savitzky-Golay filter with 9x9 smoothing window
    
    
    
    double coefficients[] = {-21.0/231.0, 14.0/231.0, 39.0/231.0, 54.0/231.0, 59.0/231.0, 54.0/231.0, 39.0/231.0, 14.0/231.0, -21.0/231.0};
    double signal_sum = 0, blur_sum = 0;
    for_loop (r, R){
        for_loop (t, T){
            signal_sum += matrix->data->at(r)->at(t);
        }
    }
    
    // setting up matrix for blured values
    matrix_t *blured = new matrix_t;
    for_loop (r, R){
        blured->push_back(new ims_spectrum);
        for_loop (t, T){
            blured->at(r)->push_back(0);
        }
    }


    
    // computing all blured values by summing all weighted values within an area
    #pragma omp parallel for collapse(2) reduction(+:blur_sum)
    for_loop (r, R){
        for_loop (t, T){
            double blurred_val = 0, weight = 0;
            for (int m = -4; m <= 4; m++){
                for (int n = -4; n <= 4; n++){
                    if (!(((r + m) < 0) || ((t + n) < 0) || ((r + m) >= R) || ((t + n) >= T))){
                        blurred_val += coefficients[4 + m] * coefficients[4 + n] * matrix->data->at(r + m)->at(t + n);
                        weight += coefficients[4 + m] * coefficients[4 + n];
                    }
                }
            }
            blur_sum += blured->at(r)->at(t) = blurred_val/ weight;
        }
    }
    
    blur_sum = 1. / blur_sum;

    // normalize
    #pragma omp parallel for collapse(2)
    for_loop (r, R){
        for_loop (t, T){
            matrix->data->at(r)->at(t) = blured->at(r)->at(t) * signal_sum * blur_sum;
        }
    }
    
    for_loop (r, R) delete blured->at(r);
    delete blured;
}



/*!
 * \brief Eliminates background noise in an MCC/IMS measurement by considering the neighbors of a point wheather they behave normal random distributed in signal dimension
 * @param matrix instance of IMSMeasurement class
 */
preprocessing(de_noising){
    IMSMatrix* matrix = measurement->data;
    int max_em_loops = DEF_MAX_EM_LOOPS;
    if (pipeline_parameters->find("max_em_loops") != pipeline_parameters->end()){
        max_em_loops = atoi(pipeline_parameters->at("max_em_loops").c_str());
    }

    
    uint R = matrix->n_rows;
    uint T = matrix->n_cols;  
    double** member = new double*[R];
    double** member_bg = new double*[R];
    double** A = new double*[R];
    for_loop(r, R){
        member[r] = new double[T];
        member_bg[r] = new double[T];
        A[r] = new double[T];
    }
    
    
    #pragma omp parallel for collapse(2)
    for_loop(r, R){
        for_loop(t, T){
            double A_z = 0, A_n = 0;
            for (int i = -4; i <= 4; i++){
                for (int j = -4; j <= 4; j++){
                    if (0 <= r + i && r + i < R && 0 <= t + j && t + j < T){
                        A_z += matrix->data->at(r + i)->at(t + j);
                        A_n++;
                    }
                }
            }
            A[r][t] = A_z / A_n;
        }
    }
    
    double mue = 0, sigma = 0, tmp = 0, tmp_bg = 0;
    double old_mue, old_sigma, old_weight_noise;


    bool run = true;
    
    double min_val = 0;
    double max_val = 0;

    
    double weight_noise = 0, weight_bg = 0;
    int number_noise = 0, number_data = 0;
    double inverse_mue = 0, lambda = 0;
    
    // searching for mue, sigma for normal distributed noise model
    // and inverse_mue, lambda for inverse normal distributed data model
    for_loop (r, R){
        for_loop (t, T){
            min_val = min(min_val, A[r][t]);
            max_val = max(max_val, A[r][t]);
            if (A[r][t] < 10){
                number_noise++;
            }
            if (A[r][t] >= 1){
                inverse_mue += A[r][t];
                number_data++;
            }
        }
    }
    inverse_mue /= double(number_data);
    weight_noise = double(number_noise) / double(R * T);
    weight_bg = (1. - weight_noise) * 0.01;
    
    for (int r = int(0.1 * R); r < R; r++){
        for (int t = int(0.1 * T); t < int(0.2 * T); t++){
            mue += A[r][t];
        }
        for (int t = int(0.9 * T); t < T; t++){
            mue += A[r][t];
        }
    }
    uint number_mue = (int(R) - int(0.1 * R)) * ((int(T) - int(0.9 * T)) + (int(0.2 * T) - int(0.1 * T)));
    mue /= double(number_mue); 
    
    

    for_loop (r, R){
        for_loop (t, T){
            if (A[r][t] >= 1)
                lambda += 1. / A[r][t] - 1. / inverse_mue;
        }
    }
    lambda = double(number_data) / lambda;
    

    
    for (int r = int(0.1 * R); r < R; r++){
        for (int t = int(0.1 * T); t < int(0.2 * T); t++){
            sigma += square(mue - A[r][t]);
        }
        for (int t = int(0.9 * T); t < T; t++){
            sigma += square(mue - A[r][t]);
        }
    }
    sigma = sqrt(sigma / double(number_mue));
    

    int loops = 0;
    double new_mue, new_inverse_mue;
    double RT = R * T;
    
    while(run && loops++ < max_em_loops){
        old_mue = mue;
        old_sigma = sigma;
        old_weight_noise = weight_noise;
        new_mue = 0, new_inverse_mue = 0, tmp = 0, tmp_bg = 0;
        
        // Expectation
        #pragma omp parallel for collapse(2) reduction(+:tmp) reduction(+:tmp_bg) reduction(+:new_mue) reduction(+:new_inverse_mue)
        for_loop(r, R){
            for_loop(t, T){
                double pdf = weight_noise * g(A[r][t], mue, sigma);
                double pdf_bg = weight_bg / RT;
                double pdf_data = (1. - weight_noise - weight_bg) * ig(A[r][t], inverse_mue, lambda, 0);
                tmp += member[r][t] = pdf / (pdf + pdf_bg + pdf_data);
                assert(0. <= member[r][t]);
                assert(member[r][t] <= 1.);
                
                tmp_bg += member_bg[r][t] = pdf_bg / (pdf + pdf_bg + pdf_data);
                assert(0. <= member_bg[r][t]);
                assert(member_bg[r][t] <= 1.);
                
                new_mue += member[r][t] * A[r][t];
                new_inverse_mue += (1. - member[r][t] - member_bg[r][t]) * A[r][t];
            }
        }
        
        
        // Maximization
        weight_noise = tmp / RT;
        weight_bg = tmp_bg / RT;
        
        mue = new_mue / tmp;
        inverse_mue = new_inverse_mue / (RT - tmp - tmp_bg);
        
        
        sigma = 0;
        lambda = 0;
        double i_inv_mue = 1. / inverse_mue;
        #pragma omp parallel for collapse(2) reduction(+:sigma) reduction(+:lambda)
        for_loop(r, R){
            for_loop(t, T){
                sigma += member[r][t] * square(mue - A[r][t]);
                if (A[r][t] != 0){
                    lambda += (1. - member[r][t] - member_bg[r][t]) * (1. / A[r][t] - i_inv_mue);
                }
            }
        }
        sigma = sqrt(sigma / tmp);
        lambda = (RT - tmp - tmp_bg) / lambda;
        
        assert(!isnan(weight_noise));       
        assert(!isinf(weight_noise));
        assert(!isnan(sigma));
        assert(!isinf(sigma));
        assert(sigma > 0);
        assert(!isnan(mue));
        assert(!isinf(mue));
        
        
        assert(!isnan(inverse_mue));
        assert(!isinf(inverse_mue));
        assert(!isnan(lambda));
        assert(!isinf(lambda));
        assert(lambda > 0);
        assert(inverse_mue > 0);
        
        
        
        run = false;
        if (!convergence(mue, old_mue, 1e-3)){
            run = true;
        }
        if (!convergence(sigma, old_sigma, 1e-3)){
            run = true;
        }
        if (!convergence(weight_noise, old_weight_noise, 1e-3)){
            run = true;
        }
    }
    
    
    #pragma omp parallel for collapse(2)
    for_loop(r, R){
        for_loop(t, T){
            matrix->data->at(r)->at(t) *= (1. - member[r][t]);
        }
    }
    
    for_loop(r, R){
        delete []member[r];
        delete []member_bg[r];
        delete []A[r];
    }
    delete []member;
    delete []member_bg;
    delete []A;
    
}
