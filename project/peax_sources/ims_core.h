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

#ifndef IMS_CORE_H
#define IMS_CORE_H

#include <string>
#include <string.h>
#include <vector>
#include <list>
#include <utility>
#include <queue>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include "sys/time.h"
#include <assert.h>
#include <algorithm>
#ifdef usexls
    #include "xlslib.h"
    #include "libxls/xls.h"
#endif
#include <iomanip>
#include <getopt.h>

#define COMPUTE_INDEX false
#define COMPUTE_TIME true

#define HALF 2
#define FLOAT 4

#define uchar unsigned char
#define ushort unsigned short
#define uint unsigned int
#define ulong unsigned long
#define for_loop(var, max) for(uint var = 0; var < max; var++)
#define parameter_insert(a, b) insert(pair<string,string>(a,b))

#define mmcchar(a) *((char*)((void*)&a))
#define mmcuchar(a) *((uchar*)((void*)&a))
#define mmcshort(a) *((int*)(void*)&a)
#define mmcushort(a) *((int*)(void*)&a)
#define mmcint(a) *((int*)((void*)&a))
#define mmcuint(a) *((uint*)((void*)&a))
#define mmclong(a) *((long*)((void*)&a))
#define mmculong(a) *((ulong*)((void*)&a))
#define mmcfloat(a) *((float*)((void*)&a))
#define mmcdouble(a) *((double*)((void*)&a))

using namespace std;
#ifdef usexls
using namespace xlslib_core;
using namespace xls;
#endif


typedef vector< vector< double >* > matrix_t;
typedef vector< double > ims_spectrum;
typedef vector< double > times;
typedef vector< string > stringlist;
typedef map<string, string> pmap;

struct detection_parameters {
    double tol_rt;
    double tol_rt_procent;
    double tol_rim;
    double min_peak_intensity;
    int area_size;
    int margin_size;
};

class IMSMatrix {
    public:
        IMSMatrix();
        IMSMatrix(const IMSMatrix* source);
        ~IMSMatrix();
        
        matrix_t* data;
        uint n_rows;
        uint n_cols;
};

class IMSMeasurement {
    public:
        IMSMeasurement(string filename);
        IMSMeasurement(const IMSMeasurement* source);
        ~IMSMeasurement();
        void input_parameter(string s, string p);
        int getRetentionIndex(double r);
        int getDriftIndex(double d);
        int getRIMIndex(double t);
        double computeRetention(double x);
//         double computeInverseMobility(double x, bool indexTime, bool dt = false);
        double computeInverseMobility(double x, bool dt = false);
        
        string filename;
        times* retention_times;
        times* drift_times;
        times* reduced_inversed_mobilities;
        IMSMatrix* data;
        map<string, string>* measurement_parameters;
        
};

// In peak_parameters we store the tolerances, as pairs <string, double>, named "t_tol = double" and "r_tol = double". 
class IMSPeak {
    public:
        IMSPeak();
        IMSPeak(const IMSPeak* peak);
        ~IMSPeak();
        void input_parameter(string, double);
        
        // Reason that name of measurement is stored in peaks is
        // because a merged layer can consist of several measurements
        string measurement_name;
        string peak_name;
        double r;
        double t;
        uint index_r;
        uint index_t;
        double signal;
        double volume;
        map<string, double>* peak_parameters;
};

class IMSPeakList {
    public:
        IMSPeakList();
        IMSPeakList(const IMSPeakList* source);
        ~IMSPeakList();
        void load(string input_filename);
        void save(string output_filename);
        
		string file_name;
        vector<IMSPeak*>* ims_peak_list;
        stringlist* parameter_names;
        IMSMeasurement* measurement_source;
};

class IMSFile {
    public:
        static IMSMeasurement* load(string filename);
        static void save(IMSMeasurement* m);
};

class IMSFileCSV : public IMSFile {
    public:
        static IMSMeasurement* load(string filename);
        static void save(IMSMeasurement* m, string filename);
};

class IMSFileIMS : public IMSFile {
    public:
        static IMSMeasurement* load(string filename);
        static void save(IMSMeasurement* m, string filename);
        static void save(IMSMeasurement* m, string filename, int bytesize);
};

class IMSLayer {
    public:
        static IMSPeakList* load(string filename);
        static void save(IMSPeakList* peaklist, string filename);
};

class IMSLayerCSV : public IMSLayer {
    public:
        static IMSPeakList* load(string filename);
        static void save(IMSPeakList* peaklist, string filename, string layer_name = "layer");
};

#ifdef usexls
class IMSLayerXLS : public IMSLayer {
    public:
        static IMSPeakList* load(string filename);
        static void save(IMSPeakList* peaklist, string filename, string layer_name = "layer");
};
#endif

// Additional functions
pmap* read_in_config(string path);
stringlist &split(const string &s, char delim, stringlist &elems);
stringlist split(const string &s, char delim);
string strip(string s, char strip = ' ');
double square(double x);
double get_mode(double mue, double lambda, double o);
void getDescriptors(double mue, double lambda, double offset, double* e, double* sigma, double* mode);
void getModelParams(double e, double sigma, double mode, double* mue, double* lambda, double* offset);
string i2s(int);
float s2f(string);
bool ends_with (string const &fullString, string const &ending);

inline double square(double x){
    return x * x;
}

/*
inline float div2(float x){
    if(!x) return 0.;
    int tmp = ((*((int*)(void*)&x)) & 2155872255) | (((((*((int*)(void*)&x)) & 2139095040) >> 23) - 1) << 23);
    return *((float*)(void*)&tmp);
}

inline float mdiv2(float x){
    if(!x) return 0.;
    int tmp = ((*((int*)(void*)&x)) & 2155872255) | (((((*((int*)(void*)&x)) & 2139095040) >> 23) - 1) << 23);
    tmp = ((~tmp) & 2147483648) | (tmp & 2147483647);
    return *((float*)(void*)&tmp);
}


inline float mul2(float x){
    if(!x) return 0.;
    int tmp = ((*((int*)(void*)&x)) & 2155872255) | (((((*((int*)(void*)&x)) & 2139095040) >> 23) + 1) << 23);
    return *((float*)(void*)&tmp);
}
*/



// function for 3x faster exp approximation
// standard error rate less than 1.2e-05
const int cishift = 52;
const int ciadd = 1022;
const int cdth = -745.133;
const double cda = 4.330968;
const double cde = 28.061659048536896 / (2. * cda * cda);
const double cdg = sqrt(2. * cde);
const double cdb = 0.4712336270551024;
const double cdd = (1. - 0.04303566602796716) * 2.;
const double cdf = cda - cdb;
const double cdilogt = 1. / log(2.);
double inline exp_less_0(double x){
    if (x < cdth)
       return 0;
    x *= cdilogt;
    long y = x;
    
    // approximation
    x = x - y;
    x = (x + cdd - square(cdg * (x + cdb)) / (x - cdf));
    y += ciadd; 
    if (y < 0){
        x /= double(1ull << -(y - 1));
        y = 1;
    }
    y <<= cishift;
    return x * (*(double*)&y);
}



// function to compute normal distributed pdf
const double nrm = sqrt(2. * M_PI);
const double hf = -0.5;
inline double g(double x, double m, double s){
    return exp(hf * square((x - m) / s)) / (s * nrm);
}



// function to compute inverse normal distributed pdf
const double tau = 2. * M_PI;
inline double ig(int x, double mue, double lambda, double off){
    return ((x > off) && (mue > 0) && (lambda > 0)) ? (exp( 0.5 * (log(lambda / (tau * ((double(x) - off) * (double(x) - off) * (double(x) - off)))) - (lambda * square(double(x) - off - mue)) / (square(mue) * (double(x) - off))))) : 0;
}



inline bool convergence(double newVal, double oldVal, double threshold){
        return fabs(newVal - oldVal) / (max(fabs(newVal), fabs(oldVal)) + 1.0) < threshold;
}


// 16 bits floating point, big endian, ieee 754-2008 implementation
class half {
    
    
    public:
        void convert_from(float val){
            value = 0;
            if (fabs(val) >= 0.00390625){
                if (val >= 2048){
                    val = 2047;
                }
                else if (val <= -2048){
                    val = -2047;
                }
                unsigned int i = mmcuint(val);
                value = (i & 2147483648) >> 16;
                value |= (i & 8388607) >> 13;
                unsigned short e = (((i & 2139095040) >> 23) - 112) & 31;
                value |= e << 10;
            }
        }
        
        float convert_to() const{
            if (!value) return 0;
            unsigned int val = ((unsigned int)(value & 32768)) << 16;
            val |= ((((unsigned int)(value) >> 10) & 31) + 112) << 23;
            val |= ((unsigned int)(value) & 1023) << 13;
            return mmcfloat(val);
        }
        
        
        // Read castings for int, long, float, double
        explicit half(int val){
            convert_from(val);
        }
        
        explicit half(long val){
            convert_from(val);
        }        
        
        explicit half(float val){
            convert_from(val);
        }
            
        explicit half(double val){
            convert_from(val);
        }
        
        explicit half(){
            value = 0;
        }
        
        
        // Write castings for int, long, float, double
        operator int() {
            return convert_to();
        }
        
        operator long() {
            return (long)convert_to();
        }
        
        operator float() {
            return convert_to();
        }
        
        operator double() {
            return (double)convert_to();
        }
        
        
        // ofstream operator
        friend ostream& operator<<(ostream &os, const half &p){
            os << p.convert_to();
            return os;
        }
        
    private:
        unsigned short value;
    
};


  


#endif /* IMS_CORE_H */
