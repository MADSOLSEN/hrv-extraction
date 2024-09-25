# Heart Rate Variability (HRV) Extraction from ECG Signals

This repository contains MATLAB code for extracting Heart Rate Variability (HRV) from Electrocardiogram (ECG) signals. The code has been applied in research on autonomic arounsals[^1] and the detection of sleep-disordered breathing[^2].

![HRV extraction algorithm on raw ECG signals](/Resources/images/RRtachogram_new-1.png)

## Features

- R peak detection algorithm from raw ECG data
- Ectopic beat detection
- Extract HRV signals from raw ECG data
- Suitable for use in sleep disorder research

## Requirements

- **MATLAB**: Version R2018b or later
- **Toolboxes**:
  - Signal Processing Toolbox
  - Bioinformatics Toolbox

## Installation

1. Clone this repository to your local machine:

   ```bash
   git clone https://github.com/MADSOLSEN/hrv-extraction.git
   ```

2. Navigate to the project directory:

   ```bash
   cd hrv-extraction
   ```
3. Open MATLAB and run the main script test.m to process a sample ECG:

   ```bash
   test.m
   ```

## Publications

This code was used in the following journal publications:

[^1]: Olsen, M., Schneider, L. D., Cheung, J., Peppard, P. E., Jennum, P. J., Mignot, E., & Sørensen, H. B. D. (2018). Automatic, electrocardiographic-based detection of autonomic arousals and their association with cortical arousals, leg movements, and respiratory events in sleep. *Sleep, 41*(3), zsy006. [https://doi.org/10.1093/sleep/zsy006](https://doi.org/10.1093/sleep/zsy006)

[^2]: Olsen, M., Mignot, E., Jennum, P. J., & Sørensen, H. B. D. (2020). Robust, ECG-based algorithm for Sleep Disordered Breathing detection in large population-based cohorts using an automatic, data-driven approach. *Sleep, 43*(5), Article zsz276. [https://doi.org/10.1093/sleep/zsz276](https://doi.org/10.1093/sleep/zsz276)

## Contact

For any inquiries or issues, please contact [somnio.ai@gmail.com](mailto:somnio.ai@gmail.com).



