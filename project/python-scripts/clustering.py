from pymixem import *
from ims_core import *
from sys import argv

if __name__ == "__main__":
    
    if len(argv) < 4:
        print("usage:", argv[0], "input_peak_list output_peak_list output_cluster_list")
        exit()
    
    # define all probability density function
    def pdf_g(parameters, data):
        return exp(-0.5 * (((data.t - parameters[1]) / parameters[3])**2 + ((data.r - parameters[2]) / parameters[4])**2)) / (2. * pi * parameters[3] * parameters[4])
    
    
    # define all maximum likelihood estimators
    def mle_weight(memberships, parameters, data_set):
        return sum(memberships) / len(data_set)
        
    def mle_g_mue_t(memberships, parameters, data_set):
        return sum(memberships[n] * data.t for n, data in enumerate(data_set)) / sum(memberships)
        
    def mle_g_mue_r(memberships, parameters, data_set):
        return sum(memberships[n] * data.r for n, data in enumerate(data_set)) / sum(memberships)
    
    # to ensure that singleton clusters do not cause a devision by zero,
    # sigma_t and sigma_r are checked not to underrun their minimal values
    def mle_g_sigma_t(memberships, parameters, data_set):
        return max(0.003, sqrt(sum(memberships[n] * (data.t - parameters[1])**2 for n, data in enumerate(data_set)) / sum(memberships)))
    def mle_g_sigma_r(memberships, parameters, data_set):
        return max((0.1 * parameters[2] + 3.0) / 3., sqrt(sum(memberships[n] * (data.r - parameters[2])**2 for n, data in enumerate(data_set)) / sum(memberships)))
    
    # intermediate step between expectation and maximization step to merge
    # all peaks that are to close to each other
    def merging(all_data):
        if all_data[5] > 0:
            parameter_set = all_data[0]
            data_set = all_data[1]
            pdf_set = all_data[2]
            mle_set = all_data[3]
            memberships = all_data[4]
            for i in range(len(parameter_set) - 1):
                for j in range(len(parameter_set) - 1, i, -1):
                    if fabs(parameter_set[i][1] - parameter_set[j][1]) < 0.003 and fabs(parameter_set[i][2] - parameter_set[j][2]) < 0.1 * parameter_set[i][2] + 3.:
                        if data_set[parameter_set[i][5][0]].signal < data_set[parameter_set[j][5][0]].signal:
                            parameter_set[i], parameter_set[j] = parameter_set[j], parameter_set[i]
                            memberships[i], memberships[j] = memberships[j], memberships[i]
                        
                        parameter_set[i][0] += parameter_set[j][0]
                        parameter_set[i][5] += parameter_set[j][5]
                        parameter_set.pop(j)
                        for d in range(len(data_set)):
                            memberships[i][d] += memberships[j][d]
                        memberships.pop(j)
                        pdf_set.pop(j)
                        mle_set.pop(j)
    
    
    # loading peak list containing all peaks from different measurements
    pl = peak_list(argv[1])
    
    # setup pdf and mle set
    pdf_set = [pdf_g] * len(pl.ims_peak_list)
    mle_set = [[mle_weight, mle_g_mue_t, mle_g_mue_r, mle_g_sigma_t, mle_g_sigma_r]] * len(pl.ims_peak_list)
    
    
    # setup parameter set
    parameter_set = []
    for i, p in enumerate(pl.ims_peak_list):
        parameter_set.append([1. / len(pl.ims_peak_list), p.t, p.r, 0.003, (0.1 * p.r + 3.) / 3, [i]])
    
    # start EM clustering
    pymixem(parameter_set, pl.ims_peak_list, pdf_set, mle_set, convergence_threshold = 1e-3, post_expectation = merging)
    
    # expand the original peak list for storing connectivity to a cluster
    pl.parameter_names.append("cluster")
    for p in pl.ims_peak_list:
        p.peak_parameters["cluster"] = -1
        
        
    # create the cluster peak list
    cluster_list = peak_list()
    cluster_list.parameter_names.append("t_tol")
    cluster_list.parameter_names.append("r_tol")
    for i, cluster in enumerate(parameter_set):
        for j in cluster[5]:
            pl.ims_peak_list[j].peak_parameters["cluster"] = i
        p = peak()
        p.measurement_name = "consensus_layer"
        p.peak_name = "C" + str(i + 1)
        p.t = cluster[1]
        p.r = cluster[2]
        p.signal = 0
        p.volume = 0
        p.index_t = 0
        p.index_r = 0
        p.peak_parameters["t_tol"] = max(0.003, cluster[3])
        p.peak_parameters["r_tol"] = max(0.1 * cluster[2] + 3.0, cluster[4])
        cluster_list.ims_peak_list.append(p)
        
    # save the modified peak list and new created cluster list
    pl.save(argv[2])
    cluster_list.save(argv[3])
    