#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import os

# Version 1.1.0

def argParsing():
  parser = argparse.ArgumentParser(description='HipFt History Plots.')

  parser.add_argument('-histfiles',
    help='A comma separated list of history files',
    type=str,
    required=False,
    default=' ')

  parser.add_argument('-rlist',
    help='A comma separated list of history file numbers to include',
    type=str,
    required=False,
    default=' ')

  parser.add_argument('-rexclude',
    help='A comma separated list of history file numbers to ignore (Overrides rlist)',
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

def GetRealization(ll,x,rl):
  for idx in range(len(ll)+1):
    if x <= ll[idx]:
        return rl[idx]
  return -1


def stats2(data,ll,rl):

  stat0 = np.min(data)
  stat1 = GetRealization(ll,np.argmin(data)+1,rl)
  stat2 = np.max(data)
  stat3 = GetRealization(ll,np.argmax(data)+1,rl)
  stat4 = np.mean(data)
  stat5 = np.std(data)
  return [stat0,stat1,stat2,stat3,stat4,stat5]


def run(args):  

  fluxp_list = []
  fluxm_list = []
  fluxp_pn_list = []
  fluxm_pn_list = []
  area_pn_list = []
  fluxp_ps_list = []
  fluxm_ps_list = []
  area_ps_list = []
  ax_dipole_list = []
  eq_dipole_list = []
  brmax_list = []
  brminabs_list = []
  brabsmin_list = []
  valerr_list = []
  lengthArrays = []
  rList = []
  length = 0

  arg_dict = vars(args)
  hist_list = arg_dict['histfiles'].split(',')
  rexclude_list = arg_dict['rexclude'].split(',')
  rlist_list = arg_dict['rlist'].split(',')

  if arg_dict['histfiles'] == ' ':
    hist_list=[]
    wDir=os.getcwd()
    for file in os.listdir(wDir):
      if "hipft_history_sol_r" in file and file.endswith(".out"):
        hist_list.append(wDir+'/'+file)
    hist_list = sorted(hist_list)

  if rlist_list[0] == ' ':
    rlist_list='all'


  ######################

  for dire in hist_list:
    r=int((dire.split('/')[-1]).replace('hipft_history_sol_r','').replace('.out',''))
    if str(r) in rexclude_list:
      continue
    elif str(r) in rlist_list or 'all' in rlist_list:
      rList.append(r)

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

      fluxp_list.append(fluxp)
      fluxm_list.append(fluxm)
      fluxp_pn_list.append(fluxp_pn)
      fluxm_pn_list.append(fluxm_pn)
      area_pn_list.append(area_pn)
      fluxp_ps_list.append(fluxp_ps)
      fluxm_ps_list.append(fluxm_ps)
      area_ps_list.append(area_ps)
      ax_dipole_list.append(ax_dipole)
      eq_dipole_list.append(eq_dipole)
      brmax_list.append(brmax)
      brminabs_list.append(np.abs(brmin))
      brabsmin_list.append(brabsmin)
      valerr_list.append(valerr)
      length+=len(time)
      lengthArrays.append(length)

      #Compute derived quantities:
      flux_tot_un = np.abs(fluxp) + np.abs(fluxm)
      flux_tot_s = fluxp + fluxm
      flux_imb = flux_tot_un*0.0
      flux_imb[flux_tot_un < 1e-15] = 0.0
      flux_imb[flux_tot_un >= 1e-15] = 100.0*flux_tot_s[flux_tot_un >= 1e-15]/flux_tot_un[flux_tot_un >= 1e-15]
      pole_n_avg_field = (fluxp_pn+fluxm_pn)/area_pn
      pole_s_avg_field = (fluxp_ps+fluxm_ps)/area_ps

      ###### INFO ##########

      fluxp_stats = stats(fluxp)
      fluxm_stats = stats(fluxm)
      
      fluxp_pn_stats = stats(fluxp_pn)
      fluxm_pn_stats = stats(fluxm_pn)
      area_pn_stats = stats(area_pn)

      fluxp_ps_stats = stats(fluxp_ps)
      fluxm_ps_stats = stats(fluxm_ps)
      area_ps_stats = stats(area_ps)

      ax_dipole_stats = stats(ax_dipole)
      eq_dipole_stats = stats(eq_dipole)

      brmax_stats = stats(brmax)
      brminabs_stats = stats(np.abs(brmin))
      brabsmin_stats = stats(brabsmin)
      valerr_stats = stats(valerr)


      print('File : '+h_file_name)
      print('Quantity  Min  Max  Mean  StdDev')
      print('flux +     ',fluxp_stats)
      print('flux -     ',fluxm_stats)
      print('flux + NP  ',fluxp_pn_stats)
      print('flux - NP  ',fluxm_pn_stats)
      print('Area NP    ',area_pn_stats)     
      print('flux + SP  ',fluxp_ps_stats)
      print('flux - SP  ',fluxm_ps_stats)
      print('Area SP    ',area_ps_stats) 
      print('Ax Dipole  ',ax_dipole_stats )
      print('Eq Dipole  ',eq_dipole_stats)                 
      print('max(Br)    ',brmax_stats)
      print('|min(Br)|  ',brminabs_stats)
      print('min(|Br|)  ',brabsmin_stats)
      print('Valerr     ',valerr_stats)
      print('')

  print('Global Statistics')
  print('Quantity  Min Min_Realization Max Max_Realization Mean  StdDev')
  print('flux +     ',stats2(fluxp_list,lengthArrays,rList))
  print('flux -     ',stats2(fluxm_list,lengthArrays,rList))
  print('flux + NP  ',stats2(fluxp_pn_list,lengthArrays,rList))
  print('flux - NP  ',stats2(fluxm_pn_list,lengthArrays,rList))
  print('Area NP    ',stats2(area_pn_list,lengthArrays,rList))
  print('flux + SP  ',stats2(fluxp_ps_list,lengthArrays,rList))
  print('flux - SP  ',stats2(fluxm_ps_list,lengthArrays,rList))
  print('Area SP    ',stats2(area_ps_list,lengthArrays,rList))
  print('Ax Dipole  ',stats2(ax_dipole_list,lengthArrays,rList))
  print('Eq Dipole  ',stats2(eq_dipole_list,lengthArrays,rList))
  print('max(Br)    ',stats2(brmax_list,lengthArrays,rList))
  print('|min(Br)|  ',stats2(brminabs_list,lengthArrays,rList))
  print('min(|Br|)  ',stats2(brabsmin_list,lengthArrays,rList))
  print('Valerr     ',stats2(valerr_list,lengthArrays,rList))
  print('')

def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()

