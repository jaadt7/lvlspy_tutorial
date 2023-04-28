import re
import math

import numpy as np
import lvlspy.level as lv
import lvlspy.spcoll as lc
import scipy.special as spc
import lvlspy.species as ls 
import lvlspy.transition as lt

from gslconsts.consts import *

def rate_mag(E_i, E_f, j, A):
    '''
    Calculate the magnetic transition rate between two states

    Args:
    - E_i: energy of the initial state (in keV)
    - E_f: energy of the final state (in keV)
    - j  : the angular momentum of the state
    - A  : mass number


    Returns:
    - The magnetic contribution to the Weisskopf estimate between the two states
    '''
    dE = E_i - E_f

    S = 2 * (j + 1)/(j*np.power(spc.factorial2(2*j+1),2)) * np.power(3/(j + 3),2)

    return(0.55 * S  * np.power(A,-2/3)    
        * np.power(dE/197000.0,2*j + 1)
        * np.power(1.4 * np.power(A, 1.0/3.0), 2*j)
        * GSL_CONST_NUM_ZETTA
    )

def rate_elec(E_i, E_f, j, A):
    '''
    Calculate the electric transition rate between two states

    Args:
    - E_i: energy of the initial state (in keV)
    - E_f: energy of the final state (in keV)
    - j  : the angular momentum of the state
    - A  : mass number


    Returns:
    - The electric contribution to the Weisskopf estimate between the two states
    '''
    
    dE = E_i - E_f

    S = 2 * (j + 1)/(j*np.power(spc.factorial2(2*j+1),2)) * np.power(3/(j + 3),2)

    return ( 2.4 * S 
        * np.power(dE / 197000.0, 2 * j + 1)
        * np.power(1.4 * np.power(A, 1.0 / 3.0), 2 * j)
        * GSL_CONST_NUM_ZETTA
    )



def weisskopf_estimate(E_i, E_f, J_i, J_f,P_i,P_f, A,mod_B):
    """
    Calculates the Weisskopf estimate for a transition between two states.

    Args:
    - E_i: energy of the initial state (in keV)
    - E_f: energy of the final state (in keV)
    - J_i: total angular momentum of the initial state
    - J_f: total angular momentum of the final state
    - P_i: parity of initial state
    - P_f: parity of final state
    - A  : mass number
    - mod_B: a string indicating whether a reduced matrix for a transition exists


    Returns:
    - The total Weisskopf estimate for the transition rate (in per second)
    """
    
    # avoid type difference in range function, cast the sum and different to integers
    diff = int(J_i - J_f)
    sm = int(J_i + J_f)
    j = range(max(1,abs(diff)),sm + 1) #range of photon angular momentum

    W = 0
    
    for jj in j:
        if np.power(-1,jj)*P_i == P_f:
            if mod_B[1] == 'E' and mod_B[2] == jj:
                B = float(mod_B[5:-1])
                r = rate_elec(E_i,E_f,jj,A)
                W += r/B
            else:
                r = rate_elec(E_i,E_f,jj,A)
                W += r/10 #Weisskopf tends to over estimate by a factor of 10 on average

        else:
            if mod_B[1] == 'M' and mod_B[2] == jj:
                B = float(mod_B[5:-1])
                r = rate_mag(E_i,E_f,jj,A)
                W += r/B
            else:
                r = rate_mag(E_i,E_f,jj,A)
                W += r/10 #Weisskopf tends to over estimate by a factor of 10 on average

    return W

def get_mass_num_and_file_name(species):
    '''
    Takes the species input, extracts the mass number and appends 0 to the start of the number to determine which file to load

    Args:
    - species (str):specifies specified by the user

    Returns:
    - file_name (str) : ENSDF file to be loaded
    - A (int): the mass number
    - l_identifier (str): the key used to help signify which lines to be read while reading in levels data
    - g_identifier (str): the key used to help signify which lines to be read while reading in gamma data
    '''

    match = re.search(r"\d+",species)
    
    if match:
        A = int(match.group())
    else:
        raise ValueError("Invalid Atomic Species Symbol")

    if len(match.group()) == 1:

        file_name  = 'ensdf.00'+match.group()
        identifier = '  ' + species
            
    elif len(match.group()) == 2:

        file_name  = 'ensdf.0'+match.group()
        identifier = ' ' + species
        
    else:

        file_name  = 'ensdf.'+match.group()
        identifier = species
    
    sym_len = len(species.replace(match.group(),''))

    if sym_len == 1:
        
        l_identifier = identifier + '   L'
        g_identifier = identifier + '   G'
    

    else:
        
        l_identifier = identifier + '  L'
        g_identifier = identifier + '  G'

    b_identifier = identifier + 'B '
    lib_sp = species.replace(match.group(),'').lower() + str(A)
    return file_name, A, lib_sp,[l_identifier, g_identifier,b_identifier]



def extract_multi_parity(jpi):
    '''
    Takes jpi as the input and extracts the j and the parity and calculates the multiplicity

    Args:
    
    - jpi (str): specifies the j and parity of the level

    Returns:
    - multi (int) : the multiplicity of the level
    - parity (str): the parity of the level
    '''

    #first strip any available parentheses
    jpi = jpi.replace('(','')
    jpi = jpi.replace(')','')
    
    if '+' not in jpi and '-' not in jpi: 
        parity = '+'
        multi = int(2*eval(jpi) + 1)
    else:
        parity = jpi[-1]       
        multi = int(2*eval(jpi[0:-1]) + 1)

    return multi, parity

def get_level_and_tran_data(file_loc,identifiers):
    '''
    Takes the file and species' name and extracts the level data from the file

    Args:
    - file_name: the ENSDF formatted file
    - speces: species specified by the user

    Returns
    - lvls: a 2D array containing the level data with the following structure 'Energy','Multiplicity','Parity' 
    - trans: a  2D array containing the indices of the level transitions and any experimental reducted matrix coefficient    
    '''
    #lvls format is (energy, multiplicity, parity)
    lvls  = []
    #trans format is (top level, bottom level, reduced matrix)
    trans = []
    zero_counter = 0 #zero counter required as to only read in the adopted values
    
    with open(file_loc,'r') as f:

        for line in f:
            #reading in level
            if line.startswith(identifiers[0]):

               energy = line[9:18]
               energy = energy.replace('X+','')
               energy = float(energy.strip()) #strip the spaces and cast to float
               if energy == 0.0: zero_counter += 1
               if zero_counter == 2: break # do not continue reading forward
               
               jpi    = line[21:38].strip() #strip the spaces
               if jpi == '':
                   print(str(energy) + r' not added to list since no $J^{\pi}$ specified')
    
               elif 'TO' in jpi or ',' in jpi or ':' in jpi:
                   print(str(energy) + ' not added since multiplet state')
               else:
                   multi,parity = extract_multi_parity(jpi)
                   lvls.append([energy, multi, parity])
                   print(energy,multi,parity)
            
            #reading in gamma info
            if line.startswith(identifiers[1]):
                E_g = line[9:18] #gamma ray energy
                E_g = float(E_g.strip())

                dE = abs(lvls[-1][0] - E_g)
               
                index = -1
                print(E_g,dE)
                B = 'N/A'
                for i in range(len(lvls)):
                    if math.isclose(abs(dE - lvls[i][0]),0.0,abs_tol= 1.0): index = i

                        
                        
                #reading in reduced matrix element in weiskopf units
                if index == -1:
                    print('transition not found from level ' + str(lvls[-1][0]) + ' with E_gamma = ' + str(E_g))
                else:
                    trans.append([len(lvls)-1,index,B])
            
            if line.startswith(identifiers[2]): 
                B = line.split()[2]
                trans[-1][2] = B        
            
          

    return lvls,trans


species = ['26AL','85KR','180TA']
directory = '../ENSDF_Database/'
s_coll = []

for i in range(len(species)):
    print('Extracting Species: '+species[i])

    file_name, A, lib_sp,identifiers = get_mass_num_and_file_name(species[i])
    file_loc = directory + file_name
    levels,transitions = get_level_and_tran_data(file_loc,identifiers)

    #generating the levels
    levs = []
    for j in range(len(levels)):
        levs.append(lv.Level(levels[j][0],levels[j][1]))
        levs[j].update_properties({'parity' : levels[j][2]})
    
    sp = ls.Species(lib_sp,levels = levs)
    
    #using the lvlspy api to help make Einstein A calculations and transition building easier
    lvs = sp.get_levels()

    for j in range(len(transitions)):
        E_i = lvs[transitions[j][0]].get_energy()
        E_f = lvs[transitions[j][1]].get_energy()

        J_i = (lvs[transitions[j][0]].get_multiplicity() - 1)/2
        J_f = (lvs[transitions[j][1]].get_multiplicity() - 1)/2

        p_i = lvs[transitions[j][0]].get_properties()['parity']
        p_f = lvs[transitions[j][1]].get_properties()['parity']

        if p_i == '+':
            p_i = 1
        else:
            p_i = -1
        
        if p_f == '+':
            p_f = 1
        else:
            p_f = -1

        ein_A = weisskopf_estimate(E_i,E_f,J_i,J_f,p_i,p_f,A,transitions[j][2])
        t = lt.Transition(lvs[transitions[j][0]],lvs[transitions[j][1]],ein_A)
        sp.add_transition(t)        


    s_coll.append(sp)

species_collection = lc.SpColl(s_coll)
species_collection.write_to_xml('nuc_collection.xml')
