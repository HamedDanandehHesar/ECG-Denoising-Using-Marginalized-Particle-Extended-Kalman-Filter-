<img width="560" height="420" alt="untitled" src="https://github.com/user-attachments/assets/efdfa7db-6839-44fb-b417-15ef7a7306f4" />
 
# ECG Denoising using Marginalized Particle Extended Kalman Filter (MP‑EKF)

MATLAB implementation of the **Marginalized Particle Extended Kalman Filter (MP‑EKF)** for model‑based ECG denoising.

This implementation corresponds to the method proposed in:

**H. Hesar, M. Mohebbi**  
*ECG Denoising Using Marginalized Particle Extended Kalman Filter with an Automatic Particle Weighting Strategy*  
IEEE Journal of Biomedical and Health Informatics, Vol. 21, No. 3, pp. 635–644, 2016  
https://ieeexplore.ieee.org/document/7494628/

---

# Overview

Electrocardiogram (ECG) signals are frequently corrupted by different types of noise such as:

• baseline wander  
• muscle artifacts  
• electrode motion artifacts  
• measurement noise  

Accurate ECG denoising is essential for reliable cardiac analysis such as:

• arrhythmia detection  
• heart rate variability analysis  
• QRS detection  
• morphological analysis  

This project implements a **Bayesian nonlinear state‑space model of the ECG signal** combined with a **Marginalized Particle Extended Kalman Filter (MP‑EKF)** to reconstruct the clean ECG signal from noisy observations.

The method exploits a **hybrid filtering strategy**:

• nonlinear states are estimated using **particle filtering**  
• linear states are estimated using **Kalman filtering**

This marginalization significantly **reduces computational complexity while improving estimation accuracy**.

---

# Key Contributions of the Method

The algorithm integrates several important ideas:

• Model‑based ECG signal representation  
• Phase‑domain modeling of the cardiac cycle  
• Gaussian mixture approximation of ECG morphology  
• Marginalized particle filtering  
• Extended Kalman filtering for linear sub‑states  
• Automatic particle weighting strategy using a synthetic ECG model  

The combination of these components produces **robust ECG denoising even at low SNR levels**.

---

# Algorithm Workflow

The complete denoising pipeline consists of the following stages.

---

# 1. Loading ECG Data

The script loads ECG data from a MATLAB `.mat` file containing:

```
x   → ECG signal matrix
fs  → sampling frequency
```

The first channel of `x` is used as the ECG signal:

```
ecg = x(1,:);
```

---

# 2. Noise Injection (for Evaluation)

To evaluate denoising performance, white Gaussian noise is added:

```
SNR = 6
x_noisy = awgn(ecg,SNR,'measured')
```

This allows quantitative performance analysis using SNR improvement.

---

# 3. R‑Peak Detection

R‑peaks are detected using the **Pan‑Tompkins algorithm**:

```
[qrs_positions] = pantompkins_qrs(x_noisy,fs)
```

R‑peaks are essential because they define the **cardiac cycle boundaries**.

---

# 4. Phase Domain Representation

ECG is modeled in the **phase domain** instead of the time domain.

The function

```
calculate_linear_phase_ver2
```

computes a **linear phase signal** based on the RR intervals.

Phase is mapped to:

```
[-π , π]
```

Each cardiac cycle therefore corresponds to a complete phase rotation.

---

# 5. Mean ECG Morphology Estimation

To obtain a representative ECG waveform:

1. ECG samples are grouped according to phase bins.
2. Mean ECG amplitude is computed in each bin.

This is performed using:

```
ECGsd_extractor_ver1
```

Outputs:

• `ECGmean` → average ECG waveform  
• `ECGsd` → standard deviation of ECG amplitude  

The mean ECG is then further smoothed using **wavelet denoising**:

```
ECGmean = wdenoise(ECGmean)
```

---

# 6. Gaussian Mixture Modeling of ECG

The ECG morphology is approximated using a **sum of Gaussian kernels**.

The ECG model is defined as:

\[
ECG(\theta) =
\sum_{i=1}^{L} a_i \exp\left(
-\frac{(\theta-\theta_i)^2}{2b_i^2}
\right)
\]

where

• \(a_i\) = amplitude  
• \(b_i\) = width  
• \(θ_i\) = center phase  

These parameters are estimated using **Particle Swarm Optimization (PSO)**.

Maximum number of Gaussian kernels:

```
L = 50
```

After estimation, only the **strongest Gaussian components** are retained.

---

# 7. Instantaneous Angular Frequency Estimation

The instantaneous angular frequency is derived from the RR intervals:

\[
\omega = \frac{2\pi fs}{RR}
\]

This variable becomes part of the **state vector**.

---

# 8. State‑Space Model

The ECG system is modeled using a nonlinear state‑space representation.

State vector:

```
x = [ECG amplitude
     phase
     angular frequency]
```

Measurement vector:

```
y = [noisy ECG
     phase
     angular frequency]
```

The nonlinear dynamics describe how ECG amplitude evolves along the phase trajectory.

---

# 9. Marginalized Particle Extended Kalman Filter

The MP‑EKF separates the state variables into:

Nonlinear states

```
ECG amplitude
ECG phase
```

Linear state

```
angular frequency
```

The algorithm operates as follows:

Particle Filter  
estimates nonlinear states

Kalman Filter  
estimates linear states conditioned on particles

This hybrid strategy significantly reduces particle dimensionality.

---

# 10. Particle Weighting Strategy

A special **automatic particle weighting strategy** is used.

Instead of only comparing particles with noisy measurements:

```
p(y | particle)
```

the algorithm also compares them with a **synthetic ECG signal** generated from the Gaussian mixture model:

```
x_synth
```

Particle weight:

```
p(y | particle) + p(x_synth | particle)
```

This improves particle selection and reduces degeneracy.

---

# 11. Resampling

When the **effective number of particles** decreases:

```
Neff = 1 / sum(w^2)
```

systematic resampling is performed.

A **roughening step** is applied to avoid particle impoverishment.

---

# 12. State Estimation

The final ECG estimate is computed as the **mean of the particles**:

```
x_estimate = mean(particles)
```

The denoised ECG corresponds to:

```
MPEKF_output(1,:)
```

---

# Output

The algorithm generates:

Denoised ECG signal

```
MPEKF_output
```

Signal‑to‑Noise Ratio improvement

```
MPEKF_SNR
```

Visualization plots:

• original ECG  
• noisy ECG  
• MP‑EKF reconstructed ECG  

---

# Figures Generated

The script produces figures including:

Figure 1  
Noisy ECG signal with detected R‑peaks.

Figure 2  
Comparison of:

Original ECG  
Noisy ECG  
MP‑EKF reconstructed ECG.

---

# MATLAB Requirements

The following toolboxes may be required:

Signal Processing Toolbox  
Wavelet Toolbox  
Global Optimization Toolbox  
Statistics Toolbox  

External helper function:

```
pantompkins_qrs.m
```

---

# Applications

This algorithm is suitable for:

ECG denoising  
Biomedical signal processing research  
Wearable ECG monitoring  
Preprocessing for arrhythmia detection  
Model‑based cardiac signal analysis  

---

# Reference

Hamed Hesar, Maryam Mohebbi  
ECG Denoising Using Marginalized Particle Extended Kalman Filter with an Automatic Particle Weighting Strategy  
IEEE Journal of Biomedical and Health Informatics  
Vol. 21, No. 3, pp. 635–644  
2016  

https://ieeexplore.ieee.org/document/7494628/

---
