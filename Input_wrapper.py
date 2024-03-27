import scipy.io
import numpy as np

class Parameters:
    
    Flag_strategy = 'Scenario 2'    # 'Strategy 1' or 'Strategy 2'
    sim_length = 24 # duration of the simulation in hours
    sim_size = 5 # duration of simulation time step in minutes
    N_step = int(sim_length*60/sim_size) # number of simulation steps
    N_max_cars = 100000 # number of cars in orignal origin-destination matrix  
    EVpen = 0.1  # share of EV penetration   
    N_colonnine_each_station = 5 # number of car slots per charging station
    
    EVpen_all = np.arange(0.01, 0.4, 0.04)  # different shares of EV penetration for MC simulations
    N_colonnine_each_station_all = [1,2,3,4,5,6,7,8,9,10]  # different number of car slots per charging statio for MC simulations
    
    OD = np.genfromtxt('OD.csv', delimiter=',')
    
    A1 = scipy.io.loadmat('path_A1_complete.mat', appendmat = False)
    path_charging_station = list(A1["path_charging_station"].flatten())
    path_caselli = A1['path_caselli'].flatten()
    path_caselli_inversed = (max(path_caselli) - path_caselli)
    path_charging_station_inversed = np.flip(max(path_caselli) - path_charging_station)
    
    # vehicles data
    SOC_min = 30 # minimun SOC to stop on the closest charging station
    SOC_min2 = 15 # minimun SOC to stop on the second closest charging station
    SOC_min3 = 15 # minimun SOC to stop on the second closest charging station for Strategy 1
    SOC_min_exit = 30 # minimun SOC to exit highway
    SOC_max = 80
    SOC_initial_max = 100 # max SOC initial of a car entering a toll booth
    SOC_initial_min = 50  # min SOC initial of a car entering a toll booth
    N_car_type = 10 # number of vehicle types
    
    
    # coefficients for consumption as a function of speed (2nd degree polynomial)
    a=0.0101
    b=-0.0557
    c=81.569
    
    # parameters for double Gaussian distribution
    a1=1
    a2=1
    mu1=8
    sigma1=2
    mu2=18
    sigma2=2
    
    EV_data = {}
    EV_data['EV1'] = {'Size': 'vehicle S', 'Battery': 50, 'Rate':30}
    EV_data['EV2'] = {'Size': 'vehicle M', 'Battery': 60, 'Rate':80}
    EV_data['EV3'] = {'Size': 'vehicle M', 'Battery': 60, 'Rate':80}
    EV_data['EV4'] = {'Size': 'vehicle M', 'Battery': 60, 'Rate':80}
    EV_data['EV5'] = {'Size': 'vehicle B', 'Battery': 80, 'Rate':120}
    EV_data['EV6'] = {'Size': 'vehicle B', 'Battery': 80, 'Rate':120}
    EV_data['EV7'] = {'Size': 'vehicle B', 'Battery': 80, 'Rate':120}
    EV_data['EV8'] = {'Size': 'vehicle B', 'Battery': 80, 'Rate':120}
    EV_data['EV9'] = {'Size': 'vehicle B', 'Battery': 80, 'Rate':120}
    EV_data['EV10'] = {'Size': 'vehicle P', 'Battery': 100, 'Rate':140}
    
    
    

    
    
