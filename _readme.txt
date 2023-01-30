The code in this folder implements the various GPS/IMU integration schemes investigated in the thesis
"GPS/IMU Integration for Space Applications" by Jordan Mitchell. Each file has a brief description on
what it does and a reference to the report for a deeper explanation. Some of the code relies on
functions from the Navigation Toolbox.

To test the algorithms simply open "main.m" and go to line 23 and replace "EKF" with either "EKFBias",
"EKF", "UKF", "SS" and "simple" depending on which integration algorithm you wish to run. Clicking 
run will then generate a trajectory, simulate gps/imu readings, pass it to the algorithm and plot the
results. To change the trajectory go to "trajectory.m" and uncomment the desired trajectory.

To characterise an IMU sensor open "allan.m" and click run. By default this will use the data in 
"ADIS16507 Gyro 8 minutes.csv". The program will figure out the Angular Random Walk, Bias Instability,
and Rate Random Walk parameters and compares them against the values in the datasheet. It also
will simulate some IMU measurements and compare the time domain and allan plots against the real data.
Note that by default the program assumes the units are in degrees to match with the datasheet.