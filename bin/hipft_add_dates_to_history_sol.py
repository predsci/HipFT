#!/usr/bin/env python3
import argparse
import signal
import sys
import numpy as np
from astropy.time import Time, TimeDelta
from pathlib import Path

def signal_handler(signal, frame):
  print('You pressed Ctrl+C! Stopping!')
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

# Version 1.0.0

def argParsing():
  parser = argparse.ArgumentParser(description='hipft_add_dates_to_history_sol:  This tool converts a hipft history sol output text file to one with dates based on a chosen start datetime (YYYY-MM-DDTHH:MM:SS).  It outputs both a file using UTC and one using TAI.')

  parser.add_argument('utstart',
                      help='Start date of first HipFT dump (regardless of HipFT''s start_time) in UTC: YYYY-MM-DDTHH:MM:SS')

  parser.add_argument('-histfiles',
                      help='A comma separated list of history files',
                      type=str,
                      required=False,
                      default='')

  return parser.parse_args()


def find_history_files():
    """Finds HipFT history files in the current directory, excluding generated _utc.out and _tai.out files."""
    return [
        file for file in Path.cwd().glob("hipft_history_sol*.out")
        if not (file.stem.endswith("_utc") or file.stem.endswith("_tai"))
    ]


def process_file(file: Path, start_date_utc):
    """Processes a single history file and generates output files with UTC and TAI times."""
    try:
        lines = file.read_text().splitlines()
        
        if len(lines) < 2:
            print(f"Skipping {file}: Not enough data.")
            return

        # Extract time values from file (skipping header)
        times = np.array([float(line.split()[1]) for line in lines[1:]])

        # Normalize times to start at 0
        times -= times[0]

        # Convert times to UTC and TAI
        start_date_tai = start_date_utc.tai
        out_times_utc = start_date_utc + TimeDelta(times * 3600.0, format="sec")
        out_times_tai = start_date_tai + TimeDelta(times * 3600.0, format="sec")

        # Define output filenames
        utc_file = file.with_name(f"{file.stem}_utc.out")
        tai_file = file.with_name(f"{file.stem}_tai.out")

        # Write to output files
        with utc_file.open("w") as utc_out, tai_file.open("w") as tai_out:
            header = lines[0].strip()
            utc_out.write(f"{header} UTC(sec) UTC\n")
            tai_out.write(f"{header} TAI(sec) TAI\n")

            for i, line in enumerate(lines[1:]):
                utc_out.write(f"{line.strip()} {out_times_utc[i].unix:.3f} {out_times_utc[i]}\n")
                tai_out.write(f"{line.strip()} {out_times_tai[i].unix_tai:.3f} {out_times_tai[i]}\n")

    except Exception as e:
        print(f"Error processing {file}: {e}")


def main():
    args = argParsing()

    try:
        # Validate and parse the start date
        start_date_utc = Time(args.utstart, format="isot", scale="utc", out_subfmt="date_hms")
    except ValueError:
        print("Error: Invalid date format. Use YYYY-MM-DDTHH:MM:SS")
        sys.exit(1)

    # Get list of history files
    hist_list = [Path(file) for file in args.histfiles.split(",") if file] if args.histfiles else find_history_files()

    if not hist_list:
        print("No history files found. Exiting.")
        sys.exit(1)

    total_files = len(hist_list)

    # Process each file
    print(f"\rProcessed 0 / {total_files} history files", end="", flush=True)
    for i, file in enumerate(hist_list, start=1):
        process_file(file, start_date_utc)
        print(f"\rProcessed {i} / {total_files} history files", end="", flush=True)

    print()


if __name__ == '__main__':
  main()
