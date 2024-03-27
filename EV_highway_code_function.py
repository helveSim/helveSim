import numpy as np
import math
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from station_update import station_update
import state_updates as su
import data_pr
from Input_wrapper import Parameters
import matplotlib.ticker as ticker

np.random.seed(10)

Params = Parameters()
Params.EV_data

EVpen = Params.EVpen
N_colonnine_each_station = Params.N_colonnine_each_station

 
ODbase = Params.OD
OD=Params.OD/(sum(sum(Params.OD)))
OD=OD*Params.N_max_cars
Num_EV=math.ceil(sum(sum(OD))*EVpen)

All_th=np.zeros(1+len(OD)*len(OD))
count1=0
for d in range(OD.shape[0]):
    for e in range(OD.shape[1]):
        All_th[count1+1]=All_th[count1]+OD[d,e]/sum(sum(OD))
        count1=count1+1
   
OD_EV=np.zeros_like(OD)

for i in range(Num_EV):
    place=np.random.rand()
    for f in range(len(All_th)):      
        if (All_th[f-1] <= place < All_th[f]):
            iii = math.ceil(f/OD.shape[0])
            jjj = f % OD.shape[0]
            if jjj == 0:
                jjj = OD.shape[1]
            OD_EV[iii-1,jjj-1]=OD_EV[iii-1,jjj-1]+1
 
M_in_out=OD_EV
All_charging_station_info = data_pr.charging_stations(Params, EVpen, N_colonnine_each_station)

# models of EVs
EV_data = Params.EV_data

# double Gaussian distribution for  distribution of EV entrances in time
X=np.arange(0,60*24+Params.sim_size, Params.sim_size)/60

pdf1=norm.pdf(X, Params.mu1, Params.sigma1)
pdf2=norm.pdf(X, Params.mu2, Params.sigma2)

pdf=Params.a1*pdf1+Params.a2*pdf2 
In_time_th=np.cumsum(pdf)/sum(pdf)
In_time_step=np.zeros(Params.N_step)
  
for h in range(int(M_in_out.sum())):
    inde=np.random.rand()
    for g in range(len(In_time_th)-1):
        if In_time_th[g] <= inde < In_time_th[g+1]:
            In_time_step[g]=In_time_step[g]+1


# allocation of vehicles in the vector creation space with OD info
vector_position = np.zeros([len(M_in_out[np.nonzero(M_in_out)]), 3])
counterx=0
for i in range(M_in_out.shape[0]):
    for j in range(M_in_out.shape[1]):
        if M_in_out[i,j]!=0:
            vector_position[counterx,0]=M_in_out[i,j] 
            vector_position[counterx,1]=i+1
            vector_position[counterx,2]=j+1
            counterx=counterx+1

All_car_info = data_pr.EV_dynamics(Params, In_time_step, vector_position, All_charging_station_info) 

#%% simulation of the EV circulation on a highway
count5= np.zeros(len(All_car_info)).astype(int)
for ii in range(Params.N_step):  
    for kkk in range(len(All_car_info)):
        if All_car_info[f"dict {kkk+1}"]['time stamp'][count5[kkk]] == (ii+1)*Params.sim_size:      
            if All_car_info[f"dict {kkk+1}"]['state flag'][-1] == 1: #%se guido
                All_car_info[f"dict {kkk+1}"]['time stamp'].append((ii+2)*Params.sim_size) #%update the time stamp
                count5[kkk] += 1
                All_car_info = su.driving(Params, All_car_info, kkk)

                # updating information on charging stations EV passes during the trip
                All_car_info = station_update(Params,All_car_info, All_charging_station_info, kkk)  

                # Update the state (form driving (1), EV can either exit highway (3), stop on the charing station (0) or continue driving (1)
                cond1 = True if any(All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1]) else False
                cond2 = (All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1][0][3] < Params.SOC_min) if any(All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1]) else False
                cond3 = (All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][count5[kkk]-1]-(All_car_info[f"dict {kkk+1}"]['km to do'][-2]*All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)']/All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)'])<Params.SOC_min_exit) if (len(All_car_info[f"dict {kkk+1}"]['charging station before exit'][-2]) == 1) else False   
                cond4 = (All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1][1][3]< Params.SOC_min2) if len(All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1])>=2 else False                
                cond5 = All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1][0][2] <= round(All_car_info[f"dict {kkk+1}"]['speed'][-1]*Params.sim_size/60,1) if any(All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1]) else False                 

                if any(All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1]):
                    cond6_idx = All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1][0][1]
                    cond6 = any(All_charging_station_info[f"dict {cond6_idx}"]["colonnine info"][:,1]==0)
                else:
                    False
                cond7 = All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1][0][3] < Params.SOC_min3 if len(All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1]) >= 1 else False 
                cond8=(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][count5[kkk]-1]-(All_car_info[f"dict {kkk+1}"]['km to do'][-2]*All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)']/All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']))<Params.SOC_min3 if len(All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1])==1 else False 
                cond9=All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1][1][3] < Params.SOC_min3 if len(All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1])>=2 else False 

                if Params.Flag_strategy == 'Scenario 1':
                    # if EV arrived to the destination (state 3)
                    if All_car_info[f"dict {kkk+1}"]['km to do'][-1] <= All_car_info[f"dict {kkk+1}"]['speed'][-1]*Params.sim_size/60:
                        All_car_info = su.exiting(Params, All_car_info, kkk, ii)       
                    elif (cond1 and (cond2 or cond3 or cond4) and cond5):  
                            All_car_info[f"dict {kkk+1}"]['state flag'].append(0)
                            All_car_info[f"dict {kkk+1}"]['speed'][-1] = 0;
                            All_car_info[f"dict {kkk+1}"]['charging station index'].append([0,0])
                            All_car_info[f"dict {kkk+1}"]['charging station index'][-1][0] = All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1][0][1]                               
                    else:   # otherwise continue driving (state 1)
                            All_car_info[f"dict {kkk+1}"]['state flag'].append(1)
                            All_car_info[f"dict {kkk+1}"]['speed'][-1]=All_car_info[f"dict {kkk+1}"]['speed'][-1]
                            All_car_info[f"dict {kkk+1}"]['charging station index'].append([])
                
                if Params.Flag_strategy == 'Scenario 2':    
                    if All_car_info[f"dict {kkk+1}"]['km to do'][-1] <= All_car_info[f"dict {kkk+1}"]['speed'][-1]*Params.sim_size/60:
                        All_car_info = su.exiting(Params, All_car_info, kkk, ii)
                    elif (cond1 and (cond2 or  cond3 or cond4) and cond5 and cond6) or (cond1 and (cond8 or cond7 or cond9) and cond5 and (not cond6)):    
                            All_car_info[f"dict {kkk+1}"]['state flag'].append(0)
                            All_car_info[f"dict {kkk+1}"]['speed'][-1] = 0;
                            All_car_info[f"dict {kkk+1}"]['charging station index'].append([0,0])
                            All_car_info[f"dict {kkk+1}"]['charging station index'][-1][0] = All_car_info[f"dict {kkk+1}"]['charging station before exit'][count5[kkk]-1][0][1]       
                    else: 
                            All_car_info[f"dict {kkk+1}"]['state flag'].append(1)
                            All_car_info[f"dict {kkk+1}"]['speed'][-1]=All_car_info[f"dict {kkk+1}"]['speed'][-1]
                            All_car_info[f"dict {kkk+1}"]['charging station index'].append([])
                         
    
    for kkk in range(len(All_car_info)):
        
        if All_car_info[f"dict {kkk+1}"]['time stamp'][count5[kkk]] == (ii+1)*Params.sim_size:
 
            # if EV stopped (state 0), we need to check if it can charge or should wait for it's turn
            if All_car_info[f"dict {kkk+1}"]['state flag'][-1] == 0:
                All_car_info[f"dict {kkk+1}"]['time stamp'].append((ii+2)*Params.sim_size)    
                count5[kkk] += 1
                
                All_car_info[f"dict {kkk+1}"]['actual position (km)'].append(All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1]) # update actual position
                All_car_info[f"dict {kkk+1}"]['speed'].append(All_car_info[f"dict {kkk+1}"]['speed'][-1])  # update actual speed 
                All_car_info[f"dict {kkk+1}"]['km to do'].append(All_car_info[f"dict {kkk+1}"]['km to do'][-1]) # update the "km to do" 
                station_index = All_car_info[f"dict {kkk+1}"]['charging station index'][-1][0] # index of the station EV stopped
                direction = All_car_info[f"dict {kkk+1}"]['direction']
                
                colonnina_index = 0
                check_queue = 0 # 0 - not in the line; 1 - in the line
                check_ready = 0 # 0 - not the first in the line; 1 - first in the line
                check_charge = 0 # avaliable charging slot on the station
                
                if direction == 1:
                    if All_car_info[f"dict {kkk+1}"]['car number'] in All_charging_station_info[f"dict {station_index}"]["queue list"][0]:
                            check_queue = 1

                    if check_queue !=1: # enter the line, if not there yet 
                        All_charging_station_info[f"dict {station_index}"]["counter queue"][0] += 1
                        All_charging_station_info[f"dict {station_index}"]["queue list"][0].append(All_car_info[f"dict {kkk+1}"]['car number'])
                        check_queue = 1
     
                    # if in the line, check if first 
                    if check_queue == 1:
                        if All_charging_station_info[f"dict {station_index}"]["queue list"][0][0] == 0 and len(All_charging_station_info[f"dict {station_index}"]["queue list"][0])>1:
                            All_charging_station_info[f"dict {station_index}"]["queue list"][0] = All_charging_station_info[f"dict {station_index}"]["queue list"][0][1:]
                        if All_charging_station_info[f"dict {station_index}"]["queue list"][0][0] == All_car_info[f"dict {kkk+1}"]['car number']: 
                            check_ready = 1  
                            for hhh in range(len(All_charging_station_info[f"dict {station_index}"]["colonnine info"])):
                                if All_charging_station_info[f"dict {station_index}"]["colonnine info"][hhh,1] == 0:
                                    check_charge = 1 # there is a free charging slot
                                    colonnina_index = hhh # index of free column
                                    break
                
                    if check_charge == 1: # EV starts charging, if possible
    
                        All_car_info[f"dict {kkk+1}"]['state flag'].append(2) 
                        
                        All_car_info[f"dict {kkk+1}"]['actual SOC (%)'].append(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1] + 100*All_car_info[f"dict {kkk+1}"]['charge rate (kW)']*(Params.sim_size/60)/All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']) 
                        All_car_info[f"dict {kkk+1}"]['range to SOC min (km)'].append(round(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1]*All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']/All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)'],1))                     
                        All_car_info[f"dict {kkk+1}"]['charging station index'].append([All_car_info[f"dict {kkk+1}"]['charging station index'][-1][0],colonnina_index+1]) 
                        N_dict = All_car_info[f"dict {kkk+1}"]['charging station index'][-1][0]          
                        All_charging_station_info[f"dict {N_dict}"]['number of cars recharged'][0] += 1 
                        All_charging_station_info[f"dict {N_dict}"]['colonnine info'][colonnina_index,1] = 1 
                        
                        if All_charging_station_info[f"dict {station_index}"]["counter queue"][0] == 1:
                            All_charging_station_info[f"dict {station_index}"]["counter queue"][0] -= 1
                            All_charging_station_info[f"dict {station_index}"]["queue list"][0] = [0]
                             
                        if All_charging_station_info[f"dict {station_index}"]["counter queue"][0]>1:
                            support = [0]*(All_charging_station_info[f"dict {station_index}"]["counter queue"][0]-1)
                            for i in range(len(support)):
                                support[i] = All_charging_station_info[f"dict {station_index}"]["queue list"][0][i+1]
                            All_charging_station_info[f"dict {station_index}"]["counter queue"][0] = All_charging_station_info[f"dict {station_index}"]["counter queue"][0] - 1
                            All_charging_station_info[f"dict {station_index}"]["queue list"][0] = support                                     
                 
                        All_charging_station_info[f"dict {N_dict}"]['station charging power vs time each vehicle'][ii,All_car_info[f"dict {kkk+1}"]['car number']] = All_car_info[f"dict {kkk+1}"]['charge rate (kW)']
                        All_car_info = station_update(Params,All_car_info, All_charging_station_info, kkk)    
            
                    if check_charge != 1: 
                        All_car_info = su.waiting(Params, All_car_info, kkk)
                        All_car_info = station_update(Params, All_car_info, All_charging_station_info, kkk)
                 
                    
                if direction == -1:       
                    if All_car_info[f"dict {kkk+1}"]['car number'] in All_charging_station_info[f"dict {station_index}"]["queue list"][1]:
                        check_queue = 1
      
                    if check_queue !=1:  
                        All_charging_station_info[f"dict {station_index}"]["counter queue"][1] += 1
                        All_charging_station_info[f"dict {station_index}"]["queue list"][1].append((All_car_info[f"dict {kkk+1}"]['car number']))
                        check_queue = 1
                     
                    if check_queue == 1:

                        if All_charging_station_info[f"dict {station_index}"]["queue list"][1][0] == 0 and len(All_charging_station_info[f"dict {station_index}"]["queue list"][1])>1:
                            All_charging_station_info[f"dict {station_index}"]["queue list"][1] = All_charging_station_info[f"dict {station_index}"]["queue list"][1][1:]
                        if All_charging_station_info[f"dict {station_index}"]["queue list"][1][0] == All_car_info[f"dict {kkk+1}"]['car number']: 
                            check_ready=1  
                            for hhh in range(len(All_charging_station_info[f"dict {station_index}"]["colonnine info"])):
                                if All_charging_station_info[f"dict {station_index}"]["colonnine info"][hhh,2] == 0:
                                    check_charge = 1 
                                    colonnina_index = hhh 
                                    break
                                
                    if check_charge == 1: 
                        All_car_info[f"dict {kkk+1}"]['state flag'].append(2) 
                            
                        All_car_info[f"dict {kkk+1}"]['actual SOC (%)'].append(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1] + 100*All_car_info[f"dict {kkk+1}"]['charge rate (kW)']*(Params.sim_size/60)/All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']) 
                        All_car_info[f"dict {kkk+1}"]['range to SOC min (km)'].append(round(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1]*All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']/All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)'],1))                    
                        All_car_info[f"dict {kkk+1}"]['charging station index'].append([All_car_info[f"dict {kkk+1}"]['charging station index'][-1][0],colonnina_index+1])                            
                        N_dict = All_car_info[f"dict {kkk+1}"]['charging station index'][-1][0]
                        All_charging_station_info[f"dict {N_dict}"]['number of cars recharged'][1] += 1 
                        All_charging_station_info[f"dict {N_dict}"]['colonnine info'][colonnina_index,2] = 1 
    
                      
                        if All_charging_station_info[f"dict {station_index}"]["counter queue"][1] == 1:
                            All_charging_station_info[f"dict {station_index}"]["counter queue"][1] -= 1
                            All_charging_station_info[f"dict {station_index}"]["queue list"][1] = [0]
         
                        if All_charging_station_info[f"dict {station_index}"]["counter queue"][1] > 1:
        
                            support = [0]*(All_charging_station_info[f"dict {station_index}"]["counter queue"][1]-1)
                            for i in range(len(support)):
                                support[i] = All_charging_station_info[f"dict {station_index}"]["queue list"][1][i+1]   
                            All_charging_station_info[f"dict {station_index}"]["counter queue"][1] = All_charging_station_info[f"dict {station_index}"]["counter queue"][1] - 1
                            All_charging_station_info[f"dict {station_index}"]["queue list"][1] = support

                        
                        All_charging_station_info[f"dict {N_dict}"]['station charging power vs time each vehicle'][ii,All_car_info[f"dict {kkk+1}"]['car number']] = All_car_info[f"dict {kkk+1}"]['charge rate (kW)']
                        All_car_info = station_update(Params,All_car_info, All_charging_station_info, kkk)    
            
                    if check_charge != 1: 
                        All_car_info = su.waiting(Params, All_car_info, kkk)
                        All_car_info = station_update(Params, All_car_info, All_charging_station_info, kkk)    

  
    for kkk in range(len(All_car_info)):
        
        if All_car_info[f"dict {kkk+1}"]['time stamp'][count5[kkk]] == (ii+1)*Params.sim_size:  

            # If EV is charging, it can either continue charging, or start moving 
            if All_car_info[f"dict {kkk+1}"]['state flag'][-1] == 2:
                All_car_info[f"dict {kkk+1}"]['time stamp'].append((ii+2)*Params.sim_size)    
                count5[kkk] += 1
                if All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1] < Params.SOC_max: # check if SOC max is reached
                    All_car_info = su.charging(Params, All_car_info, kkk)               # update the EV data
                    temp_idx = All_car_info[f"dict {kkk+1}"]['charging station index'][-1][0]
                    All_charging_station_info[f"dict {temp_idx}"]['station charging power vs time each vehicle'][ii,All_car_info[f"dict {kkk+1}"]['car number']] = All_car_info[f"dict {kkk+1}"]['charge rate (kW)']
                    All_car_info = station_update(Params,All_car_info, All_charging_station_info, kkk)

                else:
                    # update to state 1
                    All_car_info[f"dict {kkk+1}"]['state flag'].append(1)
                    All_car_info = su.driving(Params, All_car_info, kkk  )
      
                    index_s = All_car_info[f"dict {kkk+1}"]['charging station index'][-1][0]
                    index_c = All_car_info[f"dict {kkk+1}"]['charging station index'][-1][1]
                    if direction == 1:
                        All_charging_station_info[f"dict {index_s}"]['colonnine info'][index_c-1,1] = 0   
                    else:
                        All_charging_station_info[f"dict {index_s}"]['colonnine info'][index_c-1,2] = 0                        
                    All_car_info[f"dict {kkk+1}"]['charging station index'].append([])
                    All_car_info = station_update(Params,All_car_info, All_charging_station_info, kkk)


p_mean = []
for l in range(len(All_charging_station_info)):
    All_charging_station_info[f"dict {l+1}"]['station charging power vs time totale'][:,1] = All_charging_station_info[f"dict {l+1}"]['station charging power vs time each vehicle'][:,1:math.ceil(EVpen*Params.N_max_cars)].sum(axis = 1)
    temp = np.mean(All_charging_station_info[f"dict {l+1}"]['station charging power vs time totale'][:,1])
    p_mean.append(temp)

# check how many vehicles are present at the charging stations at the same time (charging or queuing)
for ii in range(Params.N_step):
    for kkk in range(len(All_car_info)):
        for yyy in range(len(All_car_info[f"dict {kkk+1}"]['time stamp'])):
            if (All_car_info[f"dict {kkk+1}"]['time stamp'][yyy] == (ii+1)*Params.sim_size) and (any(All_car_info[f"dict {kkk+1}"]['charging station index'])):    
                if any(All_car_info[f"dict {kkk+1}"]['charging station index'][yyy-1]):
                    ind=All_car_info[f"dict {kkk+1}"]['charging station index'][yyy-1][0]
                    if All_car_info[f"dict {kkk+1}"]['direction'] == 1:
                        All_charging_station_info[f"dict {ind}"]['number of cars vs time'][ii,1] += 1
                    else:
                        All_charging_station_info[f"dict {ind}"]['number of cars vs time'][ii,2] += 1
  
CAR_output_info = np.zeros((len(All_car_info),10))
counter_car_recharge_need = 0

for kkk in range(len(All_car_info)): 
    counter_stateflag2 = 0
    counter_stateflag0 = 0
    CAR_output_info[kkk,0] = All_car_info[f"dict {kkk+1}"]['km to do'][0]
    CAR_output_info[kkk,7] = All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']
    CAR_output_info[kkk,8] = All_car_info[f"dict {kkk+1}"]['arrival']   
    CAR_output_info[kkk,9] = All_car_info[f"dict {kkk+1}"]['destination']  
    
    CAR_output_info[kkk,1] = round(abs(All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1]-All_car_info[f"dict {kkk+1}"]['actual position (km)'][0])) 
    CAR_output_info[kkk,2] = round(All_car_info[f"dict {kkk+1}"]['time stamp'][-1]-All_car_info[f"dict {kkk+1}"]['time stamp'][0])/60 

    for yyy in range(len(All_car_info[f"dict {kkk+1}"]['time stamp'])-1):
        first = All_car_info[f"dict {kkk+1}"]['state flag'][yyy]
        second = All_car_info[f"dict {kkk+1}"]['state flag'][yyy+1]        
        if first == 0 and second == 2:
            CAR_output_info[kkk,3] += 1                    
        elif All_car_info[f"dict {kkk+1}"]['state flag'][yyy] == 2:
                counter_stateflag2 += 1 
                CAR_output_info[kkk,4] = counter_stateflag2 * Params.sim_size 

        if CAR_output_info[kkk,3] != 0:
           CAR_output_info[kkk,6] = CAR_output_info[kkk,4]/(CAR_output_info[kkk,4] + CAR_output_info[kkk,5])

        if first == 0 and second == 0:
                    counter_stateflag0=counter_stateflag0 + 1
                    CAR_output_info[kkk,5] = counter_stateflag0*Params.sim_size

counter_car_waiting = 0
for mmm in range(len(All_car_info)):
    if CAR_output_info[mmm,5] != 0:
        counter_car_waiting += 1

only_recharged_CAR_info = np.zeros((len(All_car_info),10))
for ooo in range(len(All_car_info)):
    if CAR_output_info[ooo,3] != 0:
        only_recharged_CAR_info[ooo,:] = CAR_output_info[ooo]

only_recharged_CAR_info = only_recharged_CAR_info[~np.all(only_recharged_CAR_info == 0, axis=1)]        


only_recharged_CAR_infoPD = pd.DataFrame(only_recharged_CAR_info)
for kkk in range(len(only_recharged_CAR_info)):
    only_recharged_CAR_infoPD.loc[kkk,10] = only_recharged_CAR_info[kkk,4] / only_recharged_CAR_info[kkk,3]


# output information
CAR_output_infoPD = pd.DataFrame(CAR_output_info)
N_EV_cars = len(All_car_info)
perc_car_recharged = round((100*(len(only_recharged_CAR_info)/len(All_car_info))),1) # percentage of total vehicles that stopped to recharge at least once
perc_car_waiting = round((100*(counter_car_waiting/len(only_recharged_CAR_info))),1) # percentage of vehicles that stopped to charge that had to wait 
t_wait_mean = round(only_recharged_CAR_infoPD[5].mean()) # average waiting time in minutes
t_wait_max = round(only_recharged_CAR_infoPD[5].max())  # maximum waiting time in minutes
t_recharge_mean = round(only_recharged_CAR_infoPD[10].mean(),1) # average charging time in minutes
tot_km_EV = CAR_output_infoPD[0].sum() # total distance traveled by all EVs
VTMG_EV = tot_km_EV/Params.path_caselli[-1] 
N_tot_recharges = only_recharged_CAR_infoPD[3].sum() # number of recharges for all vehicles
t_tot_recharges = only_recharged_CAR_infoPD[4].sum()/62
p_mean = np.round(p_mean)
      
    