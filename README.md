# Online learning of the transfer matrix of dynamic scattering media
Codes to simulate wavefront shaping experiments and optimally recover the transfer matrix of dynamic scattering media in an online and recursive fashion.


Its use in a conventional wavefront shaping experiment is illustrated in the figure below.


<img src="https://github.com/laboGigan/online_learning_TM/blob/main/pics/method_summary.png" width="1000"/>


The core of the process lies in the so-called Recursive Least-Squares (RLS) algorithm, described in pseudo-code below.


<img src="https://github.com/laboGigan/online_learning_TM/blob/main/pics/algo.png" width="600"/>


### Publication
Please find more information in our accompanying paper, "Online learning of the transfer matrix of dynamic scattering media: wavefront shaping meets multidimensional time series".


### Functions
* main_function.m : simulates wavefront shaping experiments at varying stability times of the scattering medium and noise levels,                       and produces figures comparing a conventional transfer-matrix approach with the one based on the RLS algorithm;
* TM.m : reconstructs the transfer matrix using a conventional algorithm employing N measurements, where N is the number of inputs degrees of freedom;
* RLS_TM.m : reconstructs the transfer matrix based on the RLS algorithm;
* RLS_update.m : the updating routine of the RLS algorithm, as described in the figure above;
* error_compute.m : computes the root-mean-square normalized error between the ground truth transfer matrix and its estimate.
