#!/bin/bash -l

cd ..
python code/1_0_generate_hmm_init_signal_peaks.py
python code/1_1_init_border_cpu.py
python code/1_2_generate_hmm_input.py
