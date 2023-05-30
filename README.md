# Heano-Lab-Codes

These are the C++ Codes for the paper - Abuabakar et al., Computational Modeling of Locoregional Recurrence with Spatial Structure identifies Tissue-specific Carcinogenic Profiles. Frontiers in Oncology, 2023

Authors - Sharafudeen Dahiru Abubakar and Hiroshi Haeno

How to use

For Cancer Initiation Patterns,

1. Assign a value for all parameters - r1, r2, r3, u1, u2 and u3

2. Compile in terminal

g++ 'file_name'.cpp - 'output_of_file_name'.out

3. Run the 'output_of_file_name'.out

./ 'output_of_file_name'.out

If you want to know time taken for each simulation, run

time ./ 'output_of_file_name'.out

4. Run different combination of parameters to get different patterns

5. Best to run each combination in triplicate.


For Cancer Recurrence dependency,

1. Choose a parameter you want to find its dependency on recurrence time.

2. Select range of values and number of simulations you want to run

3. Compile in terminal

g++ 'file_name'.cpp - 'output_of_file_name'.out

4. Run the 'output_of_file_name'.out

./ 'output_of_file_name'.out

If you want to know time taken for each simulation, run

time ./ 'output_of_file_name'.out



For Fitting to clinical data,

1. For each cancer type, run a Kapler-Meier curve.

2. Extract percentage of patients with recurrence at each time point

3. Choose the value of SRUN and put the appropritate percentages as the values in the dataset

4. Choose the number of simulation trials you want to run

4. Compile in terminal

g++ 'file_name'.cpp - 'output_of_file_name'.out

5. Run the 'output_of_file_name'.out

./ 'output_of_file_name'.out

If you want to know time taken for each simulation, run

time ./ 'output_of_file_name'.out

6. Might take some time to complete depending on number of trials.
