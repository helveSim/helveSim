import numpy as np
import math

# def charging_stations(path_charging_station, path_charging_station_inversed, N_colonnine_each_station, N_step, sim_size):
def charging_stations(Param, EVpen, N_colonnine_each_station):
    All_charging_station_info={} 

    for kk in range(len(Param.path_charging_station)):
        dict_keys = ('number', 'position', 'number of colonnine','colonnine info','number of cars recharged',
                     'number of cars vs time', 'counter queue', 'queue list', 'station charging power vs time each vehicle', 'station charging power vs time totale')
        charge_station_sample = dict.fromkeys(dict_keys)
        
        All_charging_station_info[f"dict {kk+1}"]=charge_station_sample  # use the same structure in every cell
        All_charging_station_info[f"dict {kk+1}"]['number']=kk+1
        All_charging_station_info[f"dict {kk+1}"]['position'] = [Param.path_charging_station[kk]]
        All_charging_station_info[f"dict {kk+1}"]['position'].append(np.flip(Param.path_charging_station_inversed)[kk])
        All_charging_station_info[f"dict {kk+1}"]['number of colonnine'] = [N_colonnine_each_station]
        All_charging_station_info[f"dict {kk+1}"]['number of colonnine'].append(N_colonnine_each_station)
        

        All_charging_station_info[f"dict {kk+1}"]['colonnine info']=np.zeros([All_charging_station_info[f"dict {kk+1}"]['number of colonnine'][0],3])
        All_charging_station_info[f"dict {kk+1}"]['number of cars recharged'] = [0,0]
        All_charging_station_info[f"dict {kk+1}"]['number of cars vs time']=np.zeros([Param.N_step,3])
        All_charging_station_info[f"dict {kk+1}"]['number of cars vs time'][:,0]=np.arange(Param.sim_size,Param.N_step*Param.sim_size+Param.sim_size,Param.sim_size, dtype=int)
        All_charging_station_info[f"dict {kk+1}"]['counter queue'] = [0,0]
        All_charging_station_info[f"dict {kk+1}"]['queue list'] = [[0]]
        All_charging_station_info[f"dict {kk+1}"]['queue list'].append([0])
        All_charging_station_info[f"dict {kk+1}"]['station charging power vs time each vehicle'] = np.zeros([Param.N_step, math.ceil(EVpen*Param.N_max_cars)+1])  
        All_charging_station_info[f"dict {kk+1}"]['station charging power vs time each vehicle'][:,0]=np.arange(Param.sim_size,Param.N_step*Param.sim_size+Param.sim_size,Param.sim_size, dtype=int)
        All_charging_station_info[f"dict {kk+1}"]['station charging power vs time totale']=np.zeros([Param.N_step,2])
        All_charging_station_info[f"dict {kk+1}"]['station charging power vs time totale'][:,0]=np.arange(Param.sim_size,Param.N_step*Param.sim_size + Param.sim_size, Param.sim_size, dtype=int)
        
        for ll in range(All_charging_station_info[f"dict {kk+1}"]['number of colonnine'][0]):
            All_charging_station_info[f"dict {kk+1}"]['colonnine info'][ll,0]=150 # maximum charging power of a charging slot
            All_charging_station_info[f"dict {kk+1}"]['colonnine info'][ll,1]=0 # status of sharging slot (0=available, 1=occupied)
    return All_charging_station_info


    
def EV_dynamics(Params, In_time_step, vector_position, All_charging_station_info):    
    All_car_count=0
    All_car_info={}

    # create the allocation matrix of the entering vehicles each time step
    M_in_step = np.zeros((len(Params.path_caselli), len(Params.path_caselli)))

    # temp = 0    
    for ii in range(Params.N_step):
        for jj in range(int(In_time_step[ii])):
            if In_time_step[ii] > 0:           
                if len(vector_position) != 1:        
                    casella = np.random.randint(1, len(vector_position))  
                else:
                    casella = 1
                iii = int((vector_position[casella-1][1]))
                jjj = int((vector_position[casella-1][2]))
                M_in_step[iii-1,jjj-1] += 1  
                vector_position[casella-1][0] -= 1  
                
                del_index = []
                for pp in range(len(vector_position)):
                    if vector_position[pp][0]==0:
                        del_index.append(pp)
                        
                for ele in sorted(del_index, reverse = True): 
                    vector_position = np.delete(vector_position,ele,0)
                

                            
                car_type_index = np.random.randint(1, Params.N_car_type)  
                All_car_count=All_car_count + 1
                Acc = All_car_count
                car_sample_keys = ('car number','enter time (min)','arrival','destination','actual position (km)','SOC initial (%)','battery energy size (kWh)',
                                   'consumption (Wh/km)','charge rate (kW)','range to SOC min (km)','state flag','actual SOC (%)','km to do','charging station before exit',
                                   'time stamp','direction','speed','charging station index')
                car_sample = dict.fromkeys(car_sample_keys)
                
                All_car_info[f"dict {Acc}"] = car_sample
                
                All_car_info[f"dict {Acc}"]['car number'] = All_car_count # car number
                All_car_info[f"dict {Acc}"]['enter time (min)'] = (ii+1)*Params.sim_size # enter time
                All_car_info[f"dict {Acc}"]['arrival'] = iii # arrival
                All_car_info[f"dict {Acc}"]['destination'] = jjj # destination   
                
                if iii>jjj:
                    All_car_info[f"dict {Acc}"]['direction'] = -1    
                    All_car_info[f"dict {Acc}"]['actual position (km)'] = [Params.path_caselli_inversed[iii-1]] # actual position: to be updated during the simulation
                    All_car_info[f"dict {Acc}"]['km to do'] = [Params.path_caselli_inversed[jjj-1] - Params.path_caselli_inversed[iii-1]]  # distance to drive until exit: to be updated during the simulation

                if iii<jjj:
                    All_car_info[f"dict {Acc}"]['direction'] = 1 
                    All_car_info[f"dict {Acc}"]['actual position (km)'] = [Params.path_caselli[iii-1]] # actual position: to be updated during the simulation
                    All_car_info[f"dict {Acc}"]['km to do'] = [abs(Params.path_caselli[iii-1] - Params.path_caselli[jjj-1])] 
                
                All_car_info[f"dict {Acc}"]['speed'] = [np.random.randint(100, 130)] #[int(Rand_speed[temp])] # average speed
                All_car_info[f"dict {Acc}"]['SOC initial (%)'] = [np.random.randint(Params.SOC_initial_min, Params.SOC_initial_max)] # intial SOC 
                All_car_info[f"dict {Acc}"]['battery energy size (kWh)'] = Params.EV_data[f'EV{car_type_index}']['Battery']    # battery energy
                All_car_info[f"dict {Acc}"]['consumption (Wh/km)'] = ((Params.a*(All_car_info[f"dict {Acc}"]['speed'][-1])**2)+Params.b*(All_car_info[f"dict {Acc}"]['speed'][-1])+Params.c)/10 # consumption 
                All_car_info[f"dict {Acc}"]['charge rate (kW)'] = Params.EV_data[f'EV{car_type_index}']['Rate'] #%charge rate
                All_car_info[f"dict {Acc}"]['range to SOC min (km)'] = [round((((All_car_info[f"dict {Acc}"]['SOC initial (%)'][-1]-0)/100)*All_car_info[f"dict {Acc}"]['battery energy size (kWh)'])/(All_car_info[f"dict {Acc}"]['consumption (Wh/km)'] /100))] #%range to SOC zero
                All_car_info[f"dict {Acc}"]['state flag'] = [1] # when a car is created, it is always in the condition of driving
                All_car_info[f"dict {Acc}"]['actual SOC (%)'] = All_car_info[f"dict {Acc}"]['SOC initial (%)'] # to be updated during the simulation

                All_car_info[f"dict {Acc}"]['charging station before exit'] = [0]
                All_car_info[f"dict {Acc}"]['charging station index'] = [[]]
                # temp = temp + 1  
                

                if All_car_info[f"dict {Acc}"]['direction'] == -1:
                    counter4=1
                    CSSI = []
                    for mm in list(reversed(range(1,len(Params.path_charging_station)+1))):
                        if All_car_info[f"dict {Acc}"]['actual position (km)'][-1] < All_charging_station_info[f"dict {mm}"]['position'][1] < Params.path_caselli_inversed[jjj-1]:

                            CSSI0 = counter4
                            CSSI1 = All_charging_station_info[f"dict {mm}"]['number'] # charging station index  
                            CSSI2 = All_charging_station_info[f"dict {mm}"]['position'][1] - All_car_info[f"dict {Acc}"]['actual position (km)'][-1] #%km to that charging station
                            CSSI3 = np.round((All_car_info[f"dict {Acc}"]['SOC initial (%)'][-1]-(CSSI2 * All_car_info[f"dict {Acc}"]['consumption (Wh/km)'])/All_car_info[f"dict {Acc}"]['battery energy size (kWh)']),1) #%arrival SOC to that charging station
        
                            CSSI.append([CSSI0,CSSI1,CSSI2,CSSI3])
                            counter4=counter4+1
                    All_car_info[f"dict {Acc}"]['charging station before exit'].append(CSSI) # charging station seen info
                    
                    
                
                if All_car_info[f"dict {Acc}"]['direction'] == 1: # check the direction (if EV moves in progressive order from zero km)
                    counter4=1
                    CSSI = []
                    for mm in range(1,len(Params.path_charging_station)+1):
                        if All_charging_station_info[f"dict {mm}"]['position'][0]>All_car_info[f"dict {Acc}"]['actual position (km)'][-1] and All_charging_station_info[f"dict {mm}"]['position'][0] < Params.path_caselli[jjj-1]:
                            CSSI0 = counter4
                            CSSI1 = All_charging_station_info[f"dict {mm}"]['number'] #%charging station index  
                            CSSI2 = abs(All_charging_station_info[f"dict {mm}"]['position'][0] - All_car_info[f"dict {Acc}"]['actual position (km)'][-1]) # km to that charging station
                            CSSI3 = np.round((All_car_info[f"dict {Acc}"]['SOC initial (%)'][-1]-(CSSI2 * All_car_info[f"dict {Acc}"]['consumption (Wh/km)'])/All_car_info[f"dict {Acc}"]['battery energy size (kWh)']),1) # arrival SOC to the charging station
        
                            CSSI.append([CSSI0,CSSI1,CSSI2,CSSI3])
                            counter4=counter4+1

                    All_car_info[f"dict {Acc}"]['charging station before exit'].append(CSSI) # charging station seen info
                                    
                All_car_info[f"dict {Acc}"][ 'time stamp']=[(ii+1)*Params.sim_size] # vehicle progressive time      
                
    for nc in range(len(All_car_info)):
        if All_car_info[f"dict {nc+1}"]['charging station before exit'][0] == 0:
            All_car_info[f"dict {nc+1}"]['charging station before exit'].pop(0)
    return All_car_info