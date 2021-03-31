# Time-frequency estimation of SCG signals

_This repository is needs cleaning. Please re-visit in a few weeks!_

This repository estimates time-frequency distribution of SCG signals using Short-time Fourier transform, chirplet transform, and smoothed-pseudo Wigner-Ville distribution.

## Example
Clone this repository. Set the MATLAB directory to this clone. Adjust the values of the input parameters of the _featureExtr.m_, and run the following command:
```ruby
[all_dom_freq2,T_spectro1] = featureExtr(Sig, Fs, t_IF, PCT_input_1, PCT_input_2, STFT_input_1, STFT_input_2, STFT_input_3,SPWV_input_1, SPWV_input_2, SPWV_input_3,Rem_out_peak)
```

## Reference
Please refer to the following manuscript:

A. Taebi and H. A. Mansy, "Analysis of seismocardiographic signals using polynomial chirplet transform and smoothed pseudo Wigner-Ville distribution," 2017 IEEE Signal Processing in Medicine and Biology Symposium (SPMB), Philadelphia, PA, 2017, pp. 1-6, doi: [10.1109/SPMB.2017.8257022](https://doi.org/10.1109/SPMB.2017.8257022).
