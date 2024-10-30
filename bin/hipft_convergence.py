#!/usr/bin/env python3
import argparse
import os

def argParsing():
  parser = argparse.ArgumentParser(description='HipFt Convergence Tests.')

  parser.add_argument('-test',
    help='Choose test: AdvPhG, AdvPhSc, AdvThG, AdvPhThG, DiffSc, DiffAdvPhSc',
    dest='test',
    type=str,
    required=True)

  parser.add_argument('-ofile',
    help='Output Convergence file name',
    dest='ofile',
    type=str,
    required=False)

  parser.add_argument('-strang',
    help='Turn on Strang splitting',
    dest='strang',
    action='store_true',
    required=False,
    default=False)

  parser.add_argument('-fnm',
    help='Flow num method : 1: FE+Upwind, 2: RK3TVD+Upwind, 3: RK3TVD+WENO3, 4: SSPRK(4,3)+WENO3, 5: Central Differencing (Default: 4)',
    dest='fnm',
    type=int,
    required=False,
    default=4)

  parser.add_argument('-dnm',
    help='Diffusion num method : 1: Explicit Euler, 2: Explicit RKL2 Super Time-stepping, 3: Explicit RKG2 Super Time-stepping (Default: 3)',
    dest='dnm',
    type=int,
    required=False,
    default=3)

  parser.add_argument('-weno_eps',
    help='Weno epsilon value. Zero sets epsilon to grid spacing. (Default: 0)',
    dest='weno_eps',
    type=float,
    required=False,
    default=0)

  parser.add_argument('-np',
    help='Starting resolution in phi (Default: 64 for phi tests; 1024 for theta tests only; 2x theta for theta/phi tests)',
    dest='np',
    type=int,
    required=False)

  parser.add_argument('-nt',
    help='Starting resolution in theta (Default: 64 for theta tests; 512 for phi tests only)',
    dest='nt',
    type=int,
    required=False)

  parser.add_argument('-dt',
    help='dt to use (Defualt is dt_cfl)',
    dest='dt',
    type=float,
    required=False)

  parser.add_argument('-tend',
    help='End time to use (Defualt is 672.0)',
    dest='tend',
    type=float,
    required=False,
    default=672.0)

  parser.add_argument('-pts',
    help='Number of resolutions to try (doubling the resolution each time) (Default: 7)',
    dest='pts',
    type=int,
    required=False,
    default=7)

  parser.add_argument('-hipft',
    help='HipFT location',
    dest='hipft',
    type=str,
    required=False)

  parser.add_argument('-htest',
    help='HipFT testsuite location',
    dest='htest',
    type=str,
    required=False)

  parser.add_argument('-errtype',
    help='Error to use for plotting: 1: RMSD, 2: MAXABS, 3: CVRMSD, 4: MAPE, 5: MAXAPE, 6: HHABS (Default: 6)',
    dest='errtype',
    type=int,
    required=False,
    default=6)

  return parser.parse_args()


def run(args):

  hipft=args.hipft if args.hipft else os.popen('which hipft').read().replace('\n','')
  htest=args.htest if args.htest else hipft.replace('bin/hipft','testsuite')
  hconv=args.htest if args.htest else hipft.replace('bin/hipft','bin/hipft_convergence_plot.py')
  hconv+=' -errtype '+str(args.errtype)

  if hipft == "":
    print("No hipft found")
    return

  testtype=args.test

  if testtype == "AdvPhG":
    maketestdir(testtype)
    text=testtype+"-fnm"+str(args.fnm)
    if args.dt:
      text=text+"-dt"+str(args.dt)
    if args.weno_eps != 0:
      text=text+"-weno_eps"+str(args.weno_eps)
    convergencefile = args.ofile if args.ofile else text+'.txt'
    createConFile(text,convergencefile)
    np_s=args.np if args.np else 64
    nt=args.nt if args.nt else 512
    np_list=[np_s * 2 ** (n - 1) for n in range(1, args.pts + 1)]
    for np in np_list:
      testname=testtype+"-nt"+str(nt)+"-np"+str(np)+"-fnm"+str(args.fnm)
      if args.dt:
        testname=testname+"-dt"+str(args.dt)
      if args.weno_eps != 0:
        testname=testname+"-weno_eps"+str(args.weno_eps)
      preparedir(testname,testtype,htest,hipft)
      print("=> Modifying hipft.in input file with chosen parameters")
      sed("res_nt",str(nt))
      sed("res_np",str(np))
      sed("time_end",str(args.tend))
      if args.dt:
        sed("dt_max",str(args.dt))
      if args.fnm == 5:
        sedi("flow_num_method",str(1))
        sedi("upwind",str(0.))
      else:
        sedi("flow_num_method",str(args.fnm))
      if args.weno_eps != 0:
        sedi("weno_eps",str(args.weno_eps))
      run_hipft()
      out=os.popen("cat diffh.log").read().replace('\n','')
      print("=> Exiting run directory "+testname)
      os.chdir("..")
      with open(convergencefile, 'a') as f:
        f.write('[ '+str(np)+','+str(nt)+','+out+' ], \\ \n')
    with open(convergencefile, 'a') as f:
      f.write(']) \n')
    print("=> Plotting Convergence")
    os.system(hconv+' -cfile '+convergencefile)
    print("=> Exiting test directory "+testtype)
    os.chdir("..")

  elif testtype == "AdvPhSc":
    maketestdir(testtype)
    text=testtype+"-fnm"+str(args.fnm)
    if args.dt:
      text=text+"-dt"+str(args.dt)
    if args.weno_eps != 0:
      text=text+"-weno_eps"+str(args.weno_eps)
    convergencefile = args.ofile if args.ofile else text+'.txt'
    createConFile(text,convergencefile)
    nt_s=args.nt if args.nt else 64
    np_s=args.np if args.np else nt_s*2
    np_list=[np_s * 2 ** (n - 1) for n in range(1, args.pts + 1)]
    nt_list=[nt_s * 2 ** (n - 1) for n in range(1, args.pts + 1)]
    for nt,np in zip(nt_list,np_list):
      testname=testtype+"-nt"+str(nt)+"-np"+str(np)+"-fnm"+str(args.fnm)
      if args.dt:
        testname=testname+"-dt"+str(args.dt)
      if args.weno_eps != 0:
        testname=testname+"-weno_eps"+str(args.weno_eps)
      preparedir(testname,testtype,htest,hipft)
      print("=> Modifying hipft.in input file with chosen parameters")
      sed("res_nt",str(nt))
      sed("res_np",str(np))
      sed("time_end",str(args.tend))
      if args.dt:
        sed("dt_max",str(args.dt))
      sed("n_realizations",'1')
      sed("advance_diffusion",".true.")
      sed("diffusion_coef_constant","0")
      if args.fnm == 5:
        sedi("flow_num_method",str(1))
        sedi("upwind",str(0.))
      else:
        sedi("flow_num_method",str(args.fnm))
      if args.weno_eps != 0:
        sedi("weno_eps",str(args.weno_eps))
      run_hipft()
      out=os.popen("cat diffh.log").read().replace('\n','')
      print("=> Exiting run directory "+testname)
      os.chdir("..")
      with open(convergencefile, 'a') as f:
        f.write('[ '+str(np)+','+str(nt)+','+out+' ], \\ \n')
    with open(convergencefile, 'a') as f:
      f.write(']) \n')
    print("=> Plotting Convergence")
    os.system(hconv+' -cfile '+convergencefile)
    print("=> Exiting test directory "+testtype)
    os.chdir("..")

  elif testtype == "AdvThG":
    maketestdir(testtype)
    text=testtype+"-fnm"+str(args.fnm)
    if args.dt:
      text=text+"-dt"+str(args.dt)
    if args.weno_eps != 0:
      text=text+"-weno_eps"+str(args.weno_eps)
    convergencefile = args.ofile if args.ofile else text+'.txt'
    createConFile(text,convergencefile)
    np=args.np if args.np else 1024
    nt_s=args.nt if args.nt else 64
    nt_list=[nt_s * 2 ** (n - 1) for n in range(1, args.pts + 1)]
    for nt in nt_list:
      testname=testtype+"-nt"+str(nt)+"-np"+str(np)+"-fnm"+str(args.fnm)
      if args.dt:
        testname=testname+"-dt"+str(args.dt)
      if args.weno_eps != 0:
        testname=testname+"-weno_eps"+str(args.weno_eps)
      preparedir(testname,testtype,htest,hipft)
      print("=> Modifying hipft.in input file with chosen parameters")
      sed("res_nt",str(nt))
      sed("res_np",str(np))
      sed("time_end",str(args.tend))
      if args.dt:
        sed("dt_max",str(args.dt))
      if args.fnm == 5:
        sedi("flow_num_method",str(1))
        sedi("upwind",str(0.))
      else:
        sedi("flow_num_method",str(args.fnm))
      if args.weno_eps != 0:
        sedi("weno_eps",str(args.weno_eps))
      run_hipft()
      out=os.popen("cat diffh.log").read().replace('\n','')
      print("=> Exiting run directory "+testname)
      os.chdir("..")
      with open(convergencefile, 'a') as f:
        f.write('[ '+str(np)+','+str(nt)+','+out+' ], \\ \n')
    with open(convergencefile, 'a') as f:
      f.write(']) \n')
    print("=> Plotting Convergence")
    os.system(hconv+' -cfile '+convergencefile)
    print("=> Exiting test directory "+testtype)
    os.chdir("..")

  elif testtype == "AdvPhThG":
    maketestdir(testtype)
    text=testtype+"-fnm"+str(args.fnm)
    if args.dt:
      text=text+"-dt"+str(args.dt)
    if args.weno_eps != 0:
      text=text+"-weno_eps"+str(args.weno_eps)
    convergencefile = args.ofile if args.ofile else text+'.txt'
    createConFile(text,convergencefile)
    nt_s=args.nt if args.nt else 64
    np_s=args.np if args.np else nt_s*2
    np_list=[np_s * 2 ** (n - 1) for n in range(1, args.pts + 1)]
    nt_list=[nt_s * 2 ** (n - 1) for n in range(1, args.pts + 1)]
    for nt,np in zip(nt_list,np_list):
      testname=testtype+"-nt"+str(nt)+"-np"+str(np)+"-fnm"+str(args.fnm)
      if args.dt:
        testname=testname+"-dt"+str(args.dt)
      if args.weno_eps != 0:
        testname=testname+"-weno_eps"+str(args.weno_eps)
      preparedir(testname,testtype,htest,hipft)
      print("=> Modifying hipft.in input file with chosen parameters")
      sed("res_nt",str(nt))
      sed("res_np",str(np))
      sed("time_end",str(args.tend))
      if args.dt:
        sed("dt_max",str(args.dt))
      if args.fnm == 5:
        sedi("flow_num_method",str(1))
        sedi("upwind",str(0.))
      else:
        sedi("flow_num_method",str(args.fnm))
      if args.weno_eps != 0:
        sedi("weno_eps",str(args.weno_eps))
      run_hipft()
      out=os.popen("cat diffh.log").read().replace('\n','')
      print("=> Exiting run directory "+testname)
      os.chdir("..")
      with open(convergencefile, 'a') as f:
        f.write('[ '+str(np)+','+str(nt)+','+out+' ], \\ \n')
    with open(convergencefile, 'a') as f:
      f.write(']) \n')
    print("=> Plotting Convergence")
    os.system(hconv+' -cfile '+convergencefile)
    print("=> Exiting test directory "+testtype)
    os.chdir("..")

  elif testtype == "DiffSc":
    maketestdir(testtype)
    text=testtype+"-dnm"+str(args.dnm)
    if args.dt:
      text=text+"-dt"+str(args.dt)
    if args.strang:
      text=text+"-strang"
    convergencefile = args.ofile if args.ofile else text+'.txt'
    createConFile(text,convergencefile)
    nt_s=args.nt if args.nt else 64
    np_s=args.np if args.np else nt_s*2
    np_list=[np_s * 2 ** (n - 1) for n in range(1, args.pts + 1)]
    nt_list=[nt_s * 2 ** (n - 1) for n in range(1, args.pts + 1)]
    for nt,np in zip(nt_list,np_list):
      testname=testtype+"-nt"+str(nt)+"-np"+str(np)+"-dnm"+str(args.dnm)
      if args.dt:
        testname=testname+"-dt"+str(args.dt)
      if args.strang:
        testname=testname+"-strang"
      preparedir(testname,testtype,htest,hipft)
      print("=> Modifying hipft.in input file with chosen parameters")
      sed("res_nt",str(nt))
      sed("res_np",str(np))
      sed("time_end",str(args.tend))
      if args.dt:
        sed("dt_max",str(args.dt))
      sed("n_realizations",'1')
      sed("diffusion_coef_constant","500")
      sedi("diffusion_num_method",str(args.dnm))
      run_hipft()
      out=os.popen("cat diffh.log").read().replace('\n','')
      print("=> Exiting run directory "+testname)
      os.chdir("..")
      with open(convergencefile, 'a') as f:
        f.write('[ '+str(np)+','+str(nt)+','+out+' ], \\ \n')
    with open(convergencefile, 'a') as f:
      f.write(']) \n')
    print("=> Plotting Convergence")
    os.system(hconv+' -cfile '+convergencefile)
    print("=> Exiting test directory "+testtype)
    os.chdir("..")

  elif testtype =="DiffAdvPhSc":
    maketestdir(testtype)
    text=testtype+"-fnm"+str(args.fnm)+"-dnm"+str(args.dnm)
    if args.dt:
      text=text+"-dt"+str(args.dt)
    if args.strang:
      text=text+"-strang"
    if args.weno_eps != 0:
      text=text+"-weno_eps"+str(args.weno_eps)
    convergencefile = args.ofile if args.ofile else text+'.txt'
    createConFile(text,convergencefile)
    nt_s=args.nt if args.nt else 64
    np_s=args.np if args.np else nt_s*2
    np_list=[np_s * 2 ** (n - 1) for n in range(1, args.pts + 1)]
    nt_list=[nt_s * 2 ** (n - 1) for n in range(1, args.pts + 1)]
    for nt,np in zip(nt_list,np_list):
      testname=testtype+"-nt"+str(nt)+"-np"+str(np)+"-fnm"+str(args.fnm)+"-dnm"+str(args.dnm)
      if args.dt:
        testname=testname+"-dt"+str(args.dt)
      if args.strang:
        testname=testname+"-strang"
      if args.weno_eps != 0:
        testname=testname+"-weno_eps"+str(args.weno_eps)
      preparedir(testname,testtype,htest,hipft)
      print("=> Modifying hipft.in input file with chosen parameters")
      sed("res_nt",str(nt))
      sed("res_np",str(np))
      sed("time_end",str(args.tend))
      if args.dt:
        sed("dt_max",str(args.dt))
      if args.strang:
        sedi("strang_splitting",".true.")
      if args.weno_eps != 0:
        sedi("weno_eps",str(args.weno_eps))
      sed("n_realizations",'1')
      sed("diffusion_coef_constant","500")
      sedi("diffusion_num_method",str(args.dnm))
      if args.fnm == 5:
        sedi("flow_num_method",str(1))
        sedi("upwind",str(0.))
      else:
        sedi("flow_num_method",str(args.fnm))
      run_hipft()
      out=os.popen("cat diffh.log").read().replace('\n','')
      print("=> Exiting run directory "+testname)
      os.chdir("..")
      with open(convergencefile, 'a') as f:
        f.write('[ '+str(np)+','+str(nt)+','+out+' ], \\ \n')
    with open(convergencefile, 'a') as f:
      f.write(']) \n')
    print("=> Plotting Convergence")
    os.system(hconv+' -cfile '+convergencefile)
    print("=> Exiting test directory "+testtype)
    os.chdir("..")



def maketestdir(testdir):
  print("=> Making test directory "+testdir)
  os.makedirs(testdir, exist_ok=True)
  print("=> Entering test directory "+testdir)
  os.chdir(testdir)


def createConFile(text,convergencefile):
  with open(convergencefile, 'w') as f:
    f.write(text+'\n')
    f.write('NP, NT, RMSD, MAXABS, CVRMSD, MAPE, MAXAPE \n')
    f.write('\n')
    f.write(text.replace('-','_')+' = np.array([ \\ \n')


def preparedir(testname,testtype,htest,hipft):
  print("=> Making run directory "+testname)
  os.makedirs(testname, exist_ok=True)
  print("Removing old run (if exists)")
  os.system('rm '+testname+'/*.out 2>/dev/null')
  os.system('rm '+testname+'/*.log 2>/dev/null')
  os.system('rm '+testname+'/*.err 2>/dev/null')
  os.system('rm '+testname+'/*.png 2>/dev/null')
  os.system('rm '+testname+'/*.in 2>/dev/null')
  os.system('rm '+testname+'/hipft 2>/dev/null')
  print("Copying input template")
  if testtype == "AdvPhG":
    os.system('cp '+htest+'/advect_gaussians_phi/input/hipft.in '+testname+'/')
  elif testtype =="AdvPhSc":
    os.system('cp '+htest+'/diffuse_advect_soccer/input/hipft.in '+testname+'/')
  elif testtype == "AdvThG":
    os.system('cp '+htest+'/advect_gaussians_theta/input/hipft.in '+testname+'/')
  elif testtype == "AdvPhThG":
    os.system('cp '+htest+'/advect_gaussians_phi_theta/input/hipft.in '+testname+'/')
  elif testtype == "DiffSc":
    os.system('cp '+htest+'/diffuse_soccer/input/hipft.in '+testname+'/')
  elif testtype =="DiffAdvPhSc":
    os.system('cp '+htest+'/diffuse_advect_soccer/input/hipft.in '+testname+'/')
  print("=> Linking to hipft binary located at '"+hipft+"'")
  os.system('ln -s '+hipft+' '+testname+'/hipft')
  print("=> Entering run directory "+testname)
  os.chdir(testname)


def sed(match,value):
  os.system('sed -i "s/.*'+match+'.*/  '+match+' = '+value+'/" "hipft.in"')


def sedi(match,value):
  os.system("sed -i '16i\ \ "+match+" = "+value+"' hipft.in")


def run_hipft():
  print("=> Running HipFT with command:")
  print('   mpiexec -np 1 ./hipft 1>hipft.log 2>hipft.err ')
  os.system('mpiexec -np 1 ./hipft 1>hipft.log 2>hipft.err ')
  print("=> Run complete!")
  print("=> Computing solution errors")
  os.system('hipft_map_diff.py hipft_brmap_final.h5 hipft_brmap_final_analytic.h5 > diffh.log')
  print("=> Plotting")
  dpi=str(72)
  os.system('plot2d -tp -cmin -1 -cmax 1 hipft_brmap_initial.h5 -dpi '+dpi)
  os.system('plot2d -tp -cmin -1 -cmax 1 hipft_brmap_final.h5 -dpi '+dpi)
  os.system('plot2d -tp -cmin -1 -cmax 1 hipft_brmap_final_analytic.h5 -dpi '+dpi)
  print("=> Solution errors:")
  print("")
  tmp = open("diffh.log", "r").read()
  print(tmp)
  print("")


def main():
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()
