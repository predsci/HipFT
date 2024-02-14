#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import os

# Version 1.0.1

def argParsing():
  parser = argparse.ArgumentParser(description='HipFt History Plots.')

  parser.add_argument('-histfiles',
    help='A comma separated list of history files',
    type=str,
    required=False,
    default=' ')

  return parser.parse_args()


def stats(data):

  data_stats = np.zeros(4)

  data_stats[0] = np.min(data)
  data_stats[1] = np.max(data)
  data_stats[2] = np.mean(data)
  data_stats[3] = np.std(data)

  return data_stats


def run(args):  


  arg_dict = vars(args)
  hist_list = arg_dict['histfiles'].split(',')

  if arg_dict['histfiles'] == ' ':
    hist_list=[]
    wDir=os.getcwd()
    for file in os.listdir(wDir):
      if "hipft_history_sol_r" in file and file.endswith(".out"):
        hist_list.append(wDir+'/'+file)
    hist_list = sorted(hist_list)

  ######################

  for dire in hist_list:
    h_file_name = dire
    hist_sol = pd.read_table(h_file_name,header=0,sep='\s+')

    time  = np.array(hist_sol['TIME'])
    fluxp = np.array(hist_sol['FLUX_POSITIVE'])
    fluxm = np.array(hist_sol['FLUX_NEGATIVE'])
    
    fluxp_pn = np.array(hist_sol['NPOLE_FLUX_POSITIVE'])
    fluxm_pn = np.array(hist_sol['NPOLE_FLUX_NEGATIVE'])
    area_pn = np.array(hist_sol['NPOLE_AREA'])
    
    fluxp_ps = np.array(hist_sol['SPOLE_FLUX_POSITIVE'])
    fluxm_ps = np.array(hist_sol['SPOLE_FLUX_NEGATIVE'])
    area_ps = np.array(hist_sol['SPOLE_AREA'])
    
    ax_dipole = np.array(hist_sol['AX_DIPOLE'])
    eq_dipole = np.array(hist_sol['EQ_DIPOLE'])

    brmin = np.array(hist_sol['BR_MIN'])
    brmax = np.array(hist_sol['BR_MAX'])
    brabsmin = np.array(hist_sol['BR_ABS_MIN'])
    valerr = np.array(hist_sol['VALIDATION_ERR_CVRMSD'])


    #Compute derived quantities:
    flux_tot_un = np.abs(fluxp) + np.abs(fluxm)
    flux_tot_s = fluxp + fluxm
    flux_imb = flux_tot_un*0.0
    flux_imb[flux_tot_un < 1e-15] = 0.0
    flux_imb[flux_tot_un >= 1e-15] = 100.0*flux_tot_s[flux_tot_un >= 1e-15]/flux_tot_un[flux_tot_un >= 1e-15]
    pole_n_avg_field = (fluxp_pn+fluxm_pn)/area_pn
    pole_s_avg_field = (fluxp_ps+fluxm_ps)/area_ps

    ###### INFO ##########

    brmax_stats = stats(brmax)
    brminabs_stats = stats(np.abs(brmin))
    brabsmin_stats = stats(brabsmin)

    print('File : '+h_file_name)
    print('Quantity  Min  Max  Mean  StdDev')
    print('max(Br)    ',brmax_stats)
    print('|min(Br)|  ',brminabs_stats)
    print('min(|Br|)  ',brabsmin_stats)
    print('')


def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()

