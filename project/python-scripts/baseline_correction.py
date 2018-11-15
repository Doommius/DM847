from pymixem import *
from sys import argv
from ims_core import *
        
if __name__ == "__main__":
    
    if len(argv) < 3:
        print("usage:", argv[0], " input_measurement_file output_measurement_file")
        exit()
        
        
    # define all probability density function
    def pdf_g(parameters, data):
        return exp(-0.5 * ((data - parameters[1]) / parameters[2])**2) / (sqrt(2 * pi) * parameters[2])
        
    def pdf_bg(parameters, data):
        return parameters[1]
        
    
    # define all maximum likelihood estimators
    def mle_weight(memberships, parameters, data_set):
        return sum(memberships) / len(data_set)
        
    def mle_g_mue(memberships, parameters, data_set):
        return sum(memberships[n] * data for n, data in enumerate(data_set)) / sum(memberships)
        
    def mle_g_sigma(memberships, parameters, data_set):
        return sqrt(sum(memberships[n] * (data - parameters[1])**2 for n, data in enumerate(data_set)) / sum(memberships))
        
    pdf_set = [pdf_g, pdf_bg]    
    mle_set = [[mle_weight, mle_g_mue, mle_g_sigma], [mle_weight]]
    
    # load measurement
    measurement = ims(argv[1])
    print("file loaded")
    
    R = len(measurement.points)
    T = len(measurement.points[0])
    print(R, T)
    
    for t in range(T):
        min_val, max_val = 0, 0
        data_set = []
        for r in range(R):
            data_set.append(measurement.points[r][t])
            min_val = min(min_val, measurement.points[r][t])
            max_val = max(max_val, measurement.points[r][t])
        histogram = [0] * int(max_val - min_val + 1)
        
        for data in data_set:
            histogram[int(data - min_val)] += 1
            
        mue = 0
        for i, data in enumerate(histogram):
            if histogram[mue] < data:
                mue = i
                
        sigma = 0
        for i, data in enumerate(data_set):
            sigma += (mue - data)**2
        sigma = sqrt(sigma / (R))
        bg_prob = 1. / (max_val - min_val)
        
        parameter_set = [[0.9, mue, sigma], [0.1, bg_prob]]
        
        # start EM baseline correction
        pymixem(parameter_set, data_set, pdf_set, mle_set)
        
        # print some output
        print(t, "of", T, "chromatograms processed")
        
        # correction of chromatogram
        lower = mue + 2. * sigma
        for r in range(R):
            old = measurement.points[r][t]
            if old > 0:
                measurement.points[r][t] = old - lower if old - lower >= 0 else 0
        
    measurement.save(argv[2])