def driving(Params, All_car_info, kkk):
    All_car_info[f"dict {kkk+1}"]['actual position (km)'].append(round((All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1] + All_car_info[f"dict {kkk+1}"]['speed'][0]*Params.sim_size/60),1)) #%update the actual position
    if any(All_car_info[f"dict {kkk+1}"]['actual SOC (%)']):
        All_car_info[f"dict {kkk+1}"]['actual SOC (%)'].append(round(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1] - 100*((All_car_info[f"dict {kkk+1}"]['speed'][0]*Params.sim_size/60)*(1/100)*All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)'])/All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)'],1)) #%update actual SOC
    else:    
        All_car_info[f"dict {kkk+1}"]['actual SOC (%)'].append(round((All_car_info[f"dict {kkk+1}"]['SOC initial (%)'][-1] - 100*((All_car_info[f"dict {kkk+1}"]['speed'][0]*Params.sim_size/60)*(1/100)*All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)'])/All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']),1)) #%update actual SOC
    All_car_info[f"dict {kkk+1}"]['km to do'].append(round(All_car_info[f"dict {kkk+1}"]['km to do'][-1]-abs(All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1] - All_car_info[f"dict {kkk+1}"]['actual position (km)'][-2]),1)) #%update the "km to do"
    All_car_info[f"dict {kkk+1}"]['range to SOC min (km)'].append(round(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1]*All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']/All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)'],1)) #%update the "range to SOC zero"  
    All_car_info[f"dict {kkk+1}"]['speed'].append(All_car_info[f"dict {kkk+1}"]['speed'][0])
    return All_car_info


def charging(Params, All_car_info, kkk):
    All_car_info[f"dict {kkk+1}"]['actual position (km)'].append(All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1]) #%update actual position (uguale a quella che c'era prima)
    All_car_info[f"dict {kkk+1}"]['speed'].append(All_car_info[f"dict {kkk+1}"]['speed'][-1])  #%update actual speed (stessa di prima)
    All_car_info[f"dict {kkk+1}"]['km to do'].append(All_car_info[f"dict {kkk+1}"]['km to do'][-1]) #%update the "km to do" (stessi di prima)
    All_car_info[f"dict {kkk+1}"]['state flag'].append(2) # continuo a caricare
    
    All_car_info[f"dict {kkk+1}"]['actual SOC (%)'].append(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1] + 100*All_car_info[f"dict {kkk+1}"]['charge rate (kW)']*(Params.sim_size/60)/All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']) #%update actual SOC
    All_car_info[f"dict {kkk+1}"]['range to SOC min (km)'].append(round(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1]*All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']/All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)'],1)) #%update the "range to SOC zero"
    All_car_info[f"dict {kkk+1}"]['charging station index'].append(All_car_info[f"dict {kkk+1}"]['charging station index'][-1]) # in quale stazione sono
    return All_car_info
    

def waiting(Params, All_car_info, kkk):
    All_car_info[f"dict {kkk+1}"]['state flag'].append(0) #%update state flag 
    All_car_info[f"dict {kkk+1}"]['actual SOC (%)'].append(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1]) # update actual SOC (stesso di prima)
    All_car_info[f"dict {kkk+1}"]['range to SOC min (km)'].append(round(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1]*All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']/All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)'],1)) #%update the "range to SOC zero"                                                
    All_car_info[f"dict {kkk+1}"]['charging station index'].append(All_car_info[f"dict {kkk+1}"]['charging station index'][-1]) # update in che stazione e colonnina sono
    return All_car_info
    

def exiting(Params, All_car_info, kkk, ii):
    All_car_info[f"dict {kkk+1}"]['time stamp'].append((ii+3)*Params.sim_size)
    All_car_info[f"dict {kkk+1}"]['state flag'].append(1)
    All_car_info[f"dict {kkk+1}"]['actual position (km)'].append(round((All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1] + All_car_info[f"dict {kkk+1}"]['speed'][-1]*Params.sim_size/60),1)) #%update the actual position
    All_car_info[f"dict {kkk+1}"]['actual SOC (%)'].append(round((All_car_info[f"dict {kkk+1}"]['SOC initial (%)'][-1] - abs(100*(( All_car_info[f"dict {kkk+1}"]['speed'][-1]*Params.sim_size/60)*(1/100)*All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)']))/All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']),1)) #%update actual SOC
    All_car_info[f"dict {kkk+1}"]['km to do'].append(round(All_car_info[f"dict {kkk+1}"]['km to do'][-1]-abs(All_car_info[f"dict {kkk+1}"]['actual position (km)'][-1] - All_car_info[f"dict {kkk+1}"]['actual position (km)'][-2]),1)) #%update the "km to do"
    All_car_info[f"dict {kkk+1}"]['range to SOC min (km)'].append(round(All_car_info[f"dict {kkk+1}"]['actual SOC (%)'][-1]*All_car_info[f"dict {kkk+1}"]['battery energy size (kWh)']/All_car_info[f"dict {kkk+1}"]['consumption (Wh/km)'],1)) #%update the "range to SOC zero"
    All_car_info[f"dict {kkk+1}"]['state flag'].append(3)
    All_car_info[f"dict {kkk+1}"]['speed'].append(All_car_info[f"dict {kkk+1}"]['speed'][-1])
    All_car_info[f"dict {kkk+1}"]['charging station index'].append([])  
    return All_car_info