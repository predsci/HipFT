#!/usr/bin/env python3
########################################################################
#
#  Version 0.0.1
#
########################################################################
#
#   'crossall' expects in format [x,y,...]
#   'sweep1d' expects in format [x] or [[x],[x,y,...]] where x is default
#   'manual' expects in format [x,y,...]
#
########################################################################

import signal
import sys
import argparse
import yaml
import os
import numpy as np

keys = [
        'diffusion_coefs',
        'flow_attenuation_values',
        'data_assimilation_mu_limits',
        'data_assimilation_mu_powers',
        'data_assimilation_lat_limits',
        'random_flux_amounts',
        'flow_meridional_p1',
        'flow_meridional_p3',
        'flow_meridional_p5',
        'flow_differential_p0',
        'flow_differential_p2',
        'flow_differential_p4'
    ]

hipft_keys = [
        'diffusion_coef_constant',
        'flow_attenuate_value',
        'assimilate_data_mu_limit',
        'assimilate_data_mu_power',
        'assimilate_data_lat_limit',
        'source_rfe_total_unsigned_flux_per_hour',
        'flow_mf_coef_p1_value',
        'flow_mf_coef_p1_value',
        'flow_mf_coef_p1_value',
        'flow_dr_coef_p0_value',
        'flow_dr_coef_p2_value',
        'flow_dr_coef_p4_value'
    ]

def signal_handler(signal, frame):
                print('You pressed Ctrl+C! Stopping!')
                sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)


def argParsing():
    parser = argparse.ArgumentParser(description='Create variable ranges for hipft.in.')

    parser.add_argument('input_yaml_file',
        help='Input yaml file.')

    parser.add_argument('hipft_file',
        help='Input hipft file.')

    return parser.parse_args()


def read_yaml(args):
    input_data = None
    with open(args.input_yaml_file, 'r') as stream:
        input_data = yaml.safe_load(stream)
    if not input_data:
        check_err(1,'Yaml file empty!')
    return input_data


def expand_cross(repeat_name, to_repeat, to_tile, repeats):
    output_list = {}
    loc_repeats = len(to_repeat)
    to_repeat = np.repeat(to_repeat, repeats) if repeats > 0 else to_repeat
    output_list[repeat_name] = to_repeat
    for item in to_tile:
        item_arr = to_tile[item]
        item_arr = np.tile(item_arr, loc_repeats) if item_arr.size > 0 else item_arr
        output_list[item] = item_arr
    repeats *= loc_repeats
    return repeats, output_list


def expand_sweep1d(repeat, arr, oarr):
    len_repeats = len(oarr[repeat]) - 1
    for item in oarr:
        if item == repeat:
            arr[item] = np.hstack([arr[item], oarr[item]])
        else:
            arr[item] = np.hstack([arr[item], np.repeat(oarr[item][0], len_repeats)])
    return arr


def check_lengths(output_dict):
    lengths = [len(value) for value in output_dict.values()]
    return all(length == lengths[0] for length in lengths)


def expand(real_params, type, keys):
    non_empty = {}
    for key in keys:
        value = real_params.get(key, [])
        arr = np.array(value, dtype=object)
        if arr.size > 0:
            if arr.ndim == 1 and isinstance(arr[0], (int, float)):
                non_empty[key] = arr
            elif arr.ndim == 1:
                if len(arr) == 2:
                    lst1, lst2 = arr
                    if len(lst1) > 1:
                        check_err(1,'invalid array format.')
                    non_empty[key] = np.concatenate([np.array(lst1), np.array(list(filter(lambda x: x != lst1[0], lst2)))])
                else:
                    check_err(1,'invalid array format.')
            elif  arr.ndim == 2 and len(arr) == 2:
                lst1, lst2 = arr
                if len(lst1) > 1:
                    check_err(1,'invalid array format.')
                non_empty[key] = np.concatenate([np.array(lst1), np.array(list(filter(lambda x: x != lst1[0], lst2)))])
            else:
                check_err(1,'invalid array format.')
    output_list = {}

    if non_empty:
        if type == 'crossall':
            for i, item in enumerate(non_empty):
                if i == 0:
                    repeats = len(non_empty[item])
                    output_list[item] = non_empty[item]
                else:
                    repeats, output_list = expand_cross(item, non_empty[item], output_list, repeats)
            return output_list
        elif type == 'sweep1d':
            output_list = {key: np.array([]) for key in non_empty.keys()}
            for item in non_empty:
                output_list = expand_sweep1d(item, output_list, non_empty)
        elif type == 'manual':
            equal_lengths = check_lengths(non_empty)
            if not equal_lengths:
                check_err(1,'The inputs for manual mode are not equal lengths.')
            output_list = non_empty
        else:
            check_err(1,'realization_parameters must be "crossall", "sweep1d", or "manual"!')
    return output_list


def write_hipft_file(output_list, hipft_key_dict, hipft_file):
    if output_list:
        n_real = len(next(iter(output_list.values())))
        try:
            with open(hipft_file, 'r') as f:
                lines = f.readlines()
            for key, value in output_list.items():
                match = hipft_key_dict.get(key)
                match_found = False
                for i, line in enumerate(lines):
                    if match in line:
                        lines[i] = f"  {match}s = {', '.join(map(str, value))}\n"
                        match_found = True
                        break
                if not match_found:
                    lines.append(f"  {match}s = {', '.join(map(str, value))}\n")
            with open(hipft_file, 'w') as f:
                f.writelines(lines)
        except FileNotFoundError:
            print(f"Error: The file '{hipft_file}' does not exist.")
        except KeyError as e:
            print(f"Error: The key '{e}' was not found in the keys list.")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
    else:
        n_real = 1
    sed('n_realizations', n_real, hipft_file)


def sed(match, value, file):
    with open(file, 'r') as f:
        lines = f.readlines()
    updated = False
    for i, line in enumerate(lines):
        if match in line:
            lines[i] = f"  {match} = {value}\n"
            updated = True
            break
        if line.strip() == '/':
            del lines[i]
    if not updated:
        lines.append(f"  {match} = {value}\n")
        lines.append('/\n')
        lines = ['!\n' if line == '\n' else line for line in lines]
        lines.append('!\n')
    with open(file, 'w') as f:
        f.writelines(lines)


def check_err(ierr,message):
    if ierr > 0:
        print(' ')
        print(message)
        print('Value of error code: '+str(ierr))
        sys.exit(1)


def main():
    hipft_key_dict = {key: hipft_keys[i] for i, key in enumerate(keys)}
    args = argParsing()
    yaml_data = read_yaml(args)
    hipft_params = yaml_data.get('hipft', None)
    if not hipft_params:
        check_err(1,'hipft section of yaml empty!')
    real_params = hipft_params.get('realization_parameters', None)
    if not real_params:
        check_err(1,'realization_parameters section of yaml empty!')
    realization_combination_mode = hipft_params.get('realization_combination_mode', None)
    if not realization_combination_mode:
        check_err(1,'realization_combination_mode empty!')
    output_list =    expand(real_params, realization_combination_mode, keys)
    write_hipft_file(output_list, hipft_key_dict, args.hipft_file)


if __name__ == '__main__':
    main()
