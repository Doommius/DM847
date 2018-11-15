How to use the python scripts:

Denoising:
- input: a noisy measurement file (consider data archive)
- output: a denoised measurement file (will be created)
- command: python3 de_niosing.py input-measurement.csv output-measurement.csv

Baseline correction:
- input: a measurement file containing a baseline (consider data archive)
- output: a baseline corrected measurement file (will be created)
- command: python3 baseline_correction.py input-measurement.csv output-measurement.csv

Clustering:
- input: a peak list containing consecutively all peak lists extracted from files  (consider peaklist-raw.csv data archive)
- ouput: the same peak list with an additional column containing information connection for every peak to one particular cluster (will be created);
         a cluster list containing information about all clusters (will be created)
- command: python3 clustering.py input-peak-list.csv output-peak-list.csv output-cluster-list.csv

(Visualization):
- input: a measurement file (consider data archive);
         [additional a peak list]
- command: python3 py_visualize.py input-measurement.csv [input-peak-list.csv]