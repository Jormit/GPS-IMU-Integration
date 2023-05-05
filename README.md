## GPS-IMU-Integration 🛰️
An exploration into various GPS/IMU integration techniques conducted as part of undergraduate thesis. Various algorithms are included and tested using simulated IMU/GPS data. MATLAB Navigation Toolbox is required to run code. 

Included Algorithms:
- IMU Only Positioning
- IMU Only Positioning with Periodic GPS Reset
- IMU/GPS Extended Kalman Filter
- IMU/GPS Unscented Kalman Filter
- IMU/GPS Extended Kalman Filter with Bias Estimation

Also included is a method for determining IMU random walk and bias instability using Allan Deviation method.

## Usage
To test the algorithms simply open "main.m" and go to line 23 and replace "EKF" with either "EKFBias",
"EKF", "UKF", "SS" and "simple" depending on which integration algorithm you wish to run. Clicking 
run will then generate a trajectory, simulate gps/imu readings, pass it to the algorithm and plot the
results. To change the trajectory go to "trajectory.m" and uncomment the desired trajectory.

To characterise an IMU sensor open "allan.m" and click run. By default this will use the data in 
"ADIS16507 Gyro 8 minutes.csv". The program will figure out the Angular Random Walk, Bias Instability,
and Rate Random Walk parameters and compares them against the values in the datasheet. It also
will simulate some IMU measurements and compare the time domain and allan plots against the real data.
Note that by default the program assumes the units are in degrees to match with the datasheet.

## Example Outputs

<img src="https://user-images.githubusercontent.com/15094591/236585556-1477a875-5ae6-4bf7-bf48-7f3547e1fe91.png" width="700">
