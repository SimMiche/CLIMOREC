CLIMOREC is composed of three main files:
  - climorec.r
  - params.txt
  - run.sh
  
1. climorec.r
  This file contains the main code of the device CLIMOREC and is commented such that someone with a few R knowledge can modify it as he wants.
  
2. params.txt
  This file has to be filled by the user of CLIMOREC for the parameters needed to compute climate index reconstruction, please let "=" as separator between the name of the parameter and the user input.
    - workdir: Enter the directory where CLIMOREC files are and the outputs will be placed, please provide an absolute path
    - name_reconstruction: Name of your experiment, a sub-directory with this name will be created in the workdir for outputs
    - path_database: Path to the proxy records database, first column be the time steps and other columns the proxy records. Only txt and csv files are allowded. Please provide an absolute path
    - path_mode: Path to the file containing the climate mode, the first column must contain the time steps while the second must be the observations of the climate mode. Only txt and csv files are allowded. Please provide an absolute path
    - y_start: Year where the reconstruction starts
    - y_stop: Year where the reconstruction stops
    - R: Number of training/testing splits
    - method: Regression method to be used, must be "rf" for Random Forest, "enet" for Elastic Net, "pls" for Partial Least Squares or "pcr" for Principal Component Regression 
    - freq_train: Relative size of the training samples given in frequency of the length of the total learning period. Must be between 0 and 1, both excluded
    - tests: Does the method apply a correlation tests on the training samples to perform a proxy records selection, set T for True and F for False
    - conf: Confidence level of the correlation tests. Must be between 0 and 1, 1 exluded. Ignored if tests is set to F
    - seed: Set the random seed (facultative). Default is 3
    - trace: Does the user wants do display the execution progress in % of time (facultative). Set T for True and F for False
   
3. run.sh
  When the params.txt file has been filled by the users, he now needs to run this file with the following command: 
  ./run.sh
  
Output
  When run.sh is finished to run, a directory named "name_reconstruction" will be created in "workdir" following the inputs in params.txt
  In this directory, the user can find the following output files:
    - <name_reconstruction>_NSCE_scores.csv: A Rx1 array that contains Nash-Sutcliffe validation scores for the individual reconstruction
    - <name_reconstruction>_RMSE_scores.csv: A Rx1 array that contains Root Mean Squared Error validation scores for the individual reconstruction
    - <name_reconstruction>_correlation_scores.csv: A Rx1 array that contains correlation validation scores for the individual reconstruction
    - <name_reconstruction>_ShapiroWilk_pvalues_residuals.csv: A Rx1 array that contains the p-values of a Shapiro-Wilk normality test on the residuals of each individual reconstruction
    - <name_reconstruction>_individual_reconstructions.csv: A Rx(y2-y1) array that contains the individual reconstructions
    - <name_reconstruction>_individual_se.csv: A Rx1 array that contains the regression standard-errors for the individual reconstruction
    - <name_reconstruction>_name_proxies.csv: Names of the proxy records used for the reconstruction
    - <name_reconstruction>_nb_records.csv: A Rx1 array that contains the number of proxy records used for each individual reconstruction
    - <name_reconstruction>_original_reconstruction.csv: Original finale reconstruction
    - <name_reconstruction>_renormalized_reconstruction.csv: Final reconstructions rescaled to the mean and the standard error of the input climate index
    - <name_reconstruction>_testing_samples.csv: An array that contains the R randomly drawn testing samples
    - <name_reconstruction>_training_samples.csv: An array that contains the R randomly drawn training samples
