# def station_update(All_car_info, All_charging_station_info, path_charging_station, path_caselli, path_caselli_inversed, kkk):  
def station_update(Params,All_car_info, All_charging_station_info, kkk):                        
    if All_car_info[f"dict {kkk+1}"]['direction'] == 1: # check the direction of EV
        counter4=1
        CSSI = []
        for mm in range(1,len(Params.path_charging_station)+1):
            if All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1] < All_charging_station_info[f"dict {mm}"]['position'][0] < Params.path_caselli[(All_car_info[f"dict {kkk+1}"]['destination'])-1]:
                CSSI0 = counter4
                CSSI1 = All_charging_station_info[f'dict {mm}']['number']  #%charging station index
                CSSI2 = abs(All_charging_station_info[f'dict {mm}']['position'][0]-All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1]) #%km to that charging station
                CSSI3 = round(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1]-(CSSI2*All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)'])/All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)'],1)
                
                CSSI.append([CSSI0,CSSI1,CSSI2,CSSI3])
                counter4 += 1

    
        if any(CSSI):
            All_car_info[f"dict {kkk+1}"]['charging station before exit'].append(CSSI)
        else:
            All_car_info[f"dict {kkk+1}"]['charging station before exit'].append([]) 
        CSSI = []
        
        return All_car_info


    if All_car_info[f"dict {kkk+1}"]['direction'] == -1: # check the direction of EV
        counter4=1
        CSSI = []
        for mm in list(reversed(range(1,len(Params.path_charging_station)+1))):
            if All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1] < All_charging_station_info[f"dict {mm}"]['position'][1] < Params.path_caselli_inversed[(All_car_info[f"dict {kkk+1}"]['destination'])-1]:
                CSSI0 = counter4
                CSSI1 = All_charging_station_info[f'dict {mm}']['number']  #%charging station index
                CSSI2 = All_charging_station_info[f'dict {mm}']['position'][1] - All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1] #%km to that charging station
                CSSI3 = round(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1]-(CSSI2*All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)'])/All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)'],1)
                  
                CSSI.append([CSSI0,CSSI1,CSSI2,CSSI3])
                counter4 += 1

    
        if any(CSSI):
            All_car_info[f"dict {kkk+1}"]['charging station before exit'].append(CSSI)
        else:
            All_car_info[f"dict {kkk+1}"]['charging station before exit'].append([]) 
        CSSI = []
        
        return All_car_info