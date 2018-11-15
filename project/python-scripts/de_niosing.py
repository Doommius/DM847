from pymixem import *
from ims_core import ims
from sys import argv

if __name__ == "__main__":
    
    if len(argv) < 3:
        print("usage:", argv[0], " input_measurement_file output_measurement_file")
        exit()
    
    
    # define all probability density function
    def pdf_g(parameters, data):
        return exp(-0.5 * ((data[2] - parameters[1]) / parameters[2])**2) / (sqrt(2 * pi) * parameters[2])
        
    def pdf_ig(parameters, data):
        if data[2] <= 0: return 0
        return exp(-0.5 * parameters[2] * (data[2] - parameters[1])**2 / (parameters[1]**2 * data[2])) * sqrt(parameters[2] / (2 * pi * data[2]**3))

    def pdf_bg(parameters, data):
        return parameters[1]
        
    
    # define all maximum likelihood estimators
    def mle_weight(memberships, parameters, data_set):
        return sum(memberships) / len(data_set)
        
    def mle_g_mue(memberships, parameters, data_set):
        return sum(memberships[n] * data[2] for n, data in enumerate(data_set)) / sum(memberships)
        
    def mle_g_sigma(memberships, parameters, data_set):
        return sqrt(sum(memberships[n] * (data[2] - parameters[1])**2 for n, data in enumerate(data_set)) / sum(memberships))

    def mle_ig_lambda(memberships, parameters, data_set):
        return sum(memberships) / sum(memberships[n] * (1. / data[2] - 1. / parameters[1]) for n, data in enumerate(data_set) if data[2] != 0)

    def mle_bg(memberships, parameters, data_set):
        return parameters[1]
        
    # setup pdf and mle set
    pdf_set = [pdf_g, pdf_ig, pdf_bg]
    mle_set = [[mle_weight, mle_g_mue, mle_g_sigma], [mle_weight, mle_g_mue, mle_ig_lambda], [mle_weight]]
    
    # load measurement
    measurement = ims(argv[1])
    print("file loaded")
    
    # setup data set
    data_set = []
    R = len(measurement.points)
    T = len(measurement.points[0])
    min_val, max_val = 0, 0
    for r in range(R):
        for t in range(T):
            val, num = 0, 0
            for i in range(-1, 2):
                for j in range(-1, 2):
                    if 0 <= r + i < R and 0 <= t + j < T:
                        val += measurement.points[r + i][t + j]
                        num += 1
            val = float(val) / float(num)
            data_set.append([r, t, val])
            min_val = min(min_val, val)
            max_val = max(max_val, val)
    
    #initial parameters
    gauss_weight, gauss_mu, gauss_sigma, inverse_weight, inverse_mu, inverse_lambda, bg_weight, bg_prob = 0., 0., 0., 0., 0., 0., 0., 1. / float(max_val - min_val)
    
    for r in range(int(0.1 * R), R):
        for t in range(int(0.1 * T), int(0.2 * T)):
            gauss_mu += data_set[r * T + t][2]
        for t in range(int(0.9 * T), T):
            gauss_mu += data_set[r * T + t][2]
    num_gauss = ((int(R) - int(0.1 * R)) * ((int(T) - int(0.9 * T)) + (int(0.2 * T) - int(0.1 * T))))
    gauss_mu /= num_gauss
    
    
    
    
    for r in range(int(0.1 * R), R):
        for t in range(int(0.1 * T), int(0.2 * T)):
            gauss_sigma += (data_set[r * T + t][2] - gauss_mu)**2
        for t in range(int(0.9 * T), T):
            gauss_sigma += (data_set[r * T + t][2] - gauss_mu)**2
    gauss_sigma = sqrt(gauss_sigma / num_gauss)
    
    
    num_mu = 0
    for r in range(R):
        for t in range(T):
            if data_set[r * T + t][2] >= gauss_mu + 2. * gauss_sigma:
                inverse_mu += data_set[r * T + t][2] 
                num_mu += 1
    inverse_mu /= num_mu
    
    for r in range(R):
        for t in range(T):
            if data_set[r * T + t][2] >= gauss_mu + 2. * gauss_sigma:
                inverse_lambda += 1. / data_set[r * T + t][2] - 1. / inverse_mu
    inverse_lambda = num_mu / inverse_lambda
    
    gauss_weight = ((R * T) - num_mu) / (R * T)
    inverse_weight = (1. - gauss_weight) * 0.999
    bg_weight = 1. - gauss_weight - inverse_weight
    
    parameter_set = [[gauss_weight, gauss_mu, gauss_sigma], [inverse_weight, inverse_mu, inverse_lambda], [bg_weight, bg_prob]]
    
    
    # setup membership matrix needed for correction
    memberships = [[0 for n in data_set] for c in parameter_set]
    
    # starting EM denoising
    pymixem(parameter_set, data_set, pdf_set, mle_set, memberships = memberships, convergence_threshold = 5e-2)
    
    
    for r in range(R):
        for t in range(T):
            measurement.points[r][t] *= (1. - memberships[0][r * T + t])
    
    measurement.save(argv[2])