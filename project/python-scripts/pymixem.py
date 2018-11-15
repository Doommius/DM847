from math import *
from copy import deepcopy

def pymixem(parameter_set, data_set, pdf_set, mle_set, memberships = 0, max_iterations = 100, convergence_threshold = 1e-3, pre_expectation = 0, post_expectation = 0, post_maximization = 0):

    if len(parameter_set) < 1:
        print("error: parameter_set is empty")
        return
        
    
    if len(parameter_set) != len(pdf_set):
        print("error: len of parameter_set unequals len of pdf_set")
        return
        
    if len(parameter_set) != len(mle_set):
        print("error: len of parameter_set unequals len of mle_set")
        return
        
    
    for j, parameters in enumerate(parameter_set):
        if len(parameters) < len(mle_set[j]):
            print("error: len of ",j, "th list in parameter_set lower than len of ", j, "th list in mle_set", sep = "")
            return
            
        
    if memberships == 0: memberships = [[0 for n in data_set] for c in parameter_set]
    
    stop, iterations = False, 0
    
    while not stop and iterations < max_iterations:
        iterations += 1
        old_parameter_set = deepcopy(parameter_set)
        
        if pre_expectation is not 0: pre_expectation([parameter_set, data_set, pdf_set, mle_set, memberships, iterations])
        
        # expectation step
        for n, data in enumerate(data_set):
            probabilities = [parameters[0] * pdf_set[c](parameters, data) for c, parameters in enumerate(parameter_set)]
            
            for c, parameters in enumerate(parameter_set):
                memberships[c][n] = probabilities[c] / sum(probabilities)
        
        if post_expectation is not 0: post_expectation([parameter_set, data_set, pdf_set, mle_set, memberships, iterations])
        
        #maximization step computeing parameters
        for c, parameters in enumerate(parameter_set):
            for l, estimator in enumerate(mle_set[c]):
                parameters[l] = estimator(memberships[c], parameters, data_set)
            
        if post_maximization is not 0: post_maximization([parameter_set, data_set, pdf_set, mle_set, memberships, iterations])
        
        #check convergence
        stop = True
        for c in range(len(mle_set)):
            parameters = parameter_set[c]
            for l in range(len(mle_set[c])):
                parameter = parameters[l]
                if fabs(parameter_set[c][l] - old_parameter_set[c][l]) / max(fabs(parameter_set[c][l]), fabs(old_parameter_set[c][l])) >= convergence_threshold:
                    stop = False
                    break
                
            if not stop: break