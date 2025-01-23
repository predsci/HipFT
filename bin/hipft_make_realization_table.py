#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import pandas as pd
import argparse

def argParsing():
  parser = argparse.ArgumentParser(description='Create table of variales that change with each run.')

  parser.add_argument('-rpfile',
    help='Realization parameters file (default is hipft_realization_parameters.out in current directory)',
    type=str,
    dest='rpfile',
    default='hipft_realization_parameters.out',
    required=False)

  parser.add_argument('-o',
    help='Name of output png',
    dest='oFile',
    required=False)

  parser.add_argument('-clist',
    help='A comma separated list of indices of columns to plot from the hipft_realization_parameters.out file (Default plots columns that change)',
    dest='clist',
    required=False)

  return parser.parse_args()


def save_dataframe_as_table(df, column_nums, filename):

    column_list = [df.columns[cl] for cl in column_nums]

    column_names=["r", "DA: $\mu_{latlim}$ (deg)", "DA: $\mu_{pow}$", "DA: $\mu_{lim}$", "DR: ${p0}$ (m/s)", \
                  "DR: ${p2}$ (m/s)", "DR: ${p4}$ (m/s)", "MF: ${p1}$ (m/s)", "MF: ${p3}$ (m/s)", "MF: ${p5}$ (m/s)", \
                    "FA (Gauss)","D: $\\nu$ (km$^2$/s)","RFE: $\Phi$ ($10^{21}$ Mx/hr)"]
    selected_columns = [column_names[i] for i in column_nums]

    if column_list is not None:
        df = df[column_list]

    column_widths = [len(col) for col in selected_columns]
    fig_width = sum(column_widths) * 0.1

    fig, ax = plt.subplots(figsize=(fig_width, len(df) * 0.1))
    ax.axis("tight")
    ax.axis("off")
    
    cell_text = df.astype(str).values

    table = ax.table(cellText=cell_text, colLabels=selected_columns, cellLoc='center', loc='center')
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)
    
    for i in range(len(df.columns) + 1):
        table.auto_set_column_width([i])
    
    for i in range(len(df) + 1): 
        for j in range(len(df.columns)):
            cell = table[(i, j)]
            cell.set_edgecolor("black")
            cell.set_linewidth(0.5)
    
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)


def run(args):
  df = pd.read_csv(args.rpfile, sep='\\s+')
  arg_dict = vars(args)

  if args.clist is not None:
    clist_list = list(map(int, arg_dict['clist'].split(',')))
  else:
    clist_list = [df.columns.get_loc(col) for col in df.columns if df[col].nunique() > 1]
       
  if not args.oFile:
     args.oFile = "hipft_realization_parameters_table.pdf"

  save_dataframe_as_table(df, clist_list, args.oFile)


def main():
  args = argParsing()
  run(args)


if __name__ == '__main__':
  main()
