import numpy as np
import pandas as pd
import math
import scipy


def get_aly_restr_lig_borecsh(res_info_csv, temperature, K_r, K_thA, K_thB, K_phiA, K_phiB, K_phiC):
    K = 8.314472*0.001  # Gas constant in kJ/mol/K
    V = 1.66            # standard volume in nm^3
    T = temperature #300.0           # Temperature in Kelvin
    res_df = pd.read_csv(res_info_csv, )
    r0 = float(res_df.iloc[0, 4])/10 # distance in A
    # K_r = 4184.0 # force constant for distance (kJ/mol/nm^2)

    thA = float(res_df.iloc[0, 5]) # Angle in rad
    # K_thA = 41.84 # force constant for angle (kJ/mol/rad^2)

    thB = float(res_df.iloc[0, 6]) # Angle in rad
    # K_thB = 41.84 # force constant for angle (kJ/mol/rad^2)

    # K_phiA = 41.84 # force constant for angle (kJ/mol/rad^2)
    # K_phiB = 41.84 # force constant for angle (kJ/mol/rad^2)
    # K_phiC = 41.84 # force constant for angle (kJ/mol/rad^2)

    arg =(
        (8.0 * math.pi**2.0 * V) / (r0**2.0 * math.sin(thA) * math.sin(thB)) 
        * 
        (
            ( (K_r * K_thA * K_thB * K_phiA * K_phiB * K_phiC)**0.5 ) / ( (2.0 * math.pi * K * T)**(3.0) )
        )
    )
    dG = - K * T * math.log(arg)
    return abs(dG)/4.184

def numerical_distance_integrand(r, r0, kr, R, T):
    r_eff = abs(r-r0)
    if r_eff <0 :
        r_eff =0
    return (r**2)*np.exp(-(kr*r_eff**2)/(2*R*T))

def numerical_angle_integrand(theta, theta0, spring_constant, R, T):
    return np.sin(theta) * np.exp(-spring_constant / (2 * R*T) * (theta - theta0) ** 2)

def numerical_torsion_integrand(phi, phi0, spring_constant, R, T):
    d_tor = phi-phi0
#     dphi = d_tor
    dphi = d_tor - np.floor(d_tor / (2 * np.pi) + 0.5) * (2 * np.pi)
    return np.exp(-spring_constant/(2 * R*T)*dphi**2)

def get_aly_restr_lig_single_atom(res_info_csv, temperature, K_r, K_ang, K_tor):
    R = 8.314472*0.001  # Gas constant in kJ/mol/K
    v0 = 1.66            # standard volume in nm^3
    T = temperature          # Temperature in Kelvin
    # K_r = 4184.0*4 # force constant for distance (kJ/mol/nm^2)
    res_df = pd.read_csv(res_info_csv, )
    r0 = float(res_df.iloc[0, 4])/10 # distance in A to distance in nm
    thetaA = float(res_df.iloc[0, 5]) # Angle in rad
    # K_ang = 41.84 # force constant for angle (kJ/mol/rad^2)
    phiA = float(res_df.iloc[0, 7]) # Torsion in rad
    # K_tor = 41.84 # force constant for torsion (kJ/mol/rad^2)
#     One distance
    rmin = max(0, r0-4*np.sqrt(R*T/K_r))
    rmax = r0+4*np.sqrt(R*T/K_r) # # Dist. which gives restraint energy = 8 RT
    I = lambda r: numerical_distance_integrand(r, r0, K_r, R, T)
    z_r = scipy.integrate.quad(I, rmin, rmax)[0]
    dg = z_r
#    One angle
    I = lambda theta: numerical_angle_integrand(theta, thetaA, K_ang, R, T)
    z_ang = scipy.integrate.quad(I, 0, np.pi)[0]
    dg *= z_ang
#    One torsion
    I = lambda phi: numerical_torsion_integrand(phi, phiA, K_tor, R, T)
    z_tor = scipy.integrate.quad(I, -np.pi, np.pi)[0]
    dg *= z_tor

    dg = -R*T*np.log(v0/dg)
    # print(dg/4.184)
    return abs(dg/4.184)

if __name__ == '__main__':
    print(get_aly_restr_lig_borecsh('res_databystd.csv'))