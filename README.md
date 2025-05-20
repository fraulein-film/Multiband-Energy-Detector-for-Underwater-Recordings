# Multiband-Energy-Detector-for-Underwater-Recordings

This MATLAB project provides an implementation of a **Multiband Energy Detector (ED)** for detecting recreational boat passages in **shallow-water underwater recordings**. The detector uses **STA/LTA (Short-Term Average / Long-Term Average)** analysis with frequency band partitioning to improve spectral sensitivity and robustness in noisy environments.

## Features

- Multiband STA/LTA-based energy detection
- Designed for shallow-water hydrophone recordings
- Configurable bandpass filtering and windowing
- Supports `.wav` files as input
- No third-party libraries needed beyond MATLAB’s Signal Processing Toolbox

---

## Requirements

- MATLAB (tested with R2023 or later)
- Signal Processing Toolbox

---

## File Overview

| File                             | Description                                         |
|----------------------------------|-----------------------------------------------------|
| `Multiband_Energy_Detector.m`| Main script for running the energy detector         |
| `sec_to_ddhhmmss.m`             | Utility function to format seconds into hh:mm:ss   |
| `README.md`                     | This readme file                                    |
| `LICENSE`                       | MIT license for reuse                               |

> **Note**: Example audio recordings are **not included** due to data sharing restrictions. Please use your own underwater `.wav` files to test the script.

---

## Usage

1. **Clone or download** this repository.

2. **Open `Multiband_Energy_Detector.m` in MATLAB**.

3. **Save `sec_to_ddhhmmss.m`** in the **same folder** as `Multiband_Energy_Detector.m`.

4. Add your `.wav` files to the folder and define them in the recording array inside the script:

```matlab
recording = {"your_file.wav"};
% OR
recording = {"your_file1.wav", "your_file2.wav"};
```

5. Click the **Run** button in the MATLAB interface or run from the Command Window.

---

## Input Format

- Input files must be **.wav** recordings. 
- Recommended: Underwater recordings with detectable boat passages. 

---

## Citation 

If you use this code in your work, please cite: 

T. Muhr, “Automatic detection of boats in underwater recordings using multiband energy detection,” Master’s Thesis, Chalmers University of Technology, 2025.

---

### `LICENSE` (MIT)

MIT License

Copyright (c) 2025 Theresia Muhr

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## Acknowledgments

Thanks to IVL Swedish Environmental Research Institute and Björn Wrede for support and field recordings used in the original thesis.