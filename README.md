# Unscented Kalman Filter Project
In this project I implement an Unscented Kalman Filter to make Sensor Fusion with Radar and Lidar measurements.
The code is written in C++, because this language provides great performance in terms of memory and speed (if we do things the right way).
This work is part of the Self-Driving Car Engineer Nanodegree Program.

---

[//]: # (Image References)
[simulation]: ./img/simulation.png
[sample1]: ./img/px_vs_py.png
[nis]: ./img/nis.png

## Results
### Simulator
In order to test the program I used a visualization tool provided by Udacity that plots noisy lidar and radar measurements as a simulated car drives in a figure eight pattern. 
However, we should note that in the real world application the radar will be placed fixed to the car, so the origin of the sensor axis will always be at the same position as the car. In the real case we are tracking pedestrians, not the car itself.
In the next figure we can see the lidar and radar position measurements (as red and blue circles respectively) as well as the filter estimations (as green triangles).

![Simulator][simulation]

We can also see the Root Mean Squared Error (RMSE) of the estimations for each step in the UI. 
This numbers may be considered as standard deviation errors for the estimations. 
That means that our position estimates have an error of approximately 0.1m while for the velocity we have 0.4 m/s.

### Input measurements
The program was also tested with a set of input measurements from the files in the folder `data\`.

#### Data
- For the input `obj_pose-laser-radar-synthetic-input.txt` the RMSE values obtained were:
  ```
  Accuracy - RMSE:
  px: 0.0726829
  py: 0.0841594
  vx: 0.346685
  vz: 0.236226
  ```
  
Here it's a plot of the x and y position measured, estimated and ground truth values

![Data 1][sample1]

In order to tune the process noise parameters I plotted the Normalized Innovation Squared (NIS):

![NIS][nis]


## Source files hierarchy and description
- main.cpp : Main program. Reads in measurements and apply Sensor Fusion algorithm.
  - ukf.cpp/.h : Implements Unscented Kalman Filter for both Lidar and Radar measurements.
    - tools.cpp/.h : Common math tools.
    - measurement_package.h : Class for measurements values.
  - ground_truth_package.h : Class for ground truth values.

## Application dependencies
To build the program there are some applications you have to install on your computer

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make` 
   * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt output.txt`

