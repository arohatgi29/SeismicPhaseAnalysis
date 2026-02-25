# Seismic Phase Analysis Toolkit

Tools for analyzing, quantifying, and correcting seismic phase distortions in land seismic data, with emphasis on near-surface scattering effects, phase variability, and circular-statistics–based methods.

## Overview

Seismic phase is highly sensitive to small-scale near-surface heterogeneities that scatter the wavefield and introduce trace-dependent distortions. These distortions degrade waveform coherence, stacking power, resolution, and effective bandwidth.

This repository provides algorithms and workflows for:

- Phase extraction in time and frequency domains  
- Circular statistics for wrapped phase analysis  
- Phase variability quantification  
- Phase masking and stabilization  
- Acquisition and processing quality control (QC)  
- Visualization of phase behavior in prestack data  

The toolkit is intended for research, acquisition QC, advanced processing, and development of new phase-based seismic attributes.

---

## Key Features

### Phase Analysis

- Instantaneous phase computation  
- Spectral phase estimation via FFT  
- Wrapped phase handling (−π to π)  
- Optional phase unwrapping  
- Trace-dependent phase characterization  

### Circular Statistics

Phase is treated as a circular variable:

- Circular mean  
- Mean resultant length  
- Circular variance  
- von Mises concentration parameter (κ) estimation  
- Phase clustering and dispersion analysis  

### Phase QC Metrics

- Phase variability across traces  
- Frequency-dependent stability  
- Offset-dependent behavior  
- Noise-cone analysis  
- Comparison of raw vs processed data  

### Phase Masking & Stabilization

- Local reference phase estimation  
- Neighborhood circular averaging  
- Frequency-domain masking  
- Amplitude-preserving corrections  
- Coherence improvement without kinematic distortion  

### Visualization Tools

- Phase maps on acquisition geometry  
- Offset-domain displays  
- Frequency-dependent plots  
- Wiggle / variable-area seismic sections  
- Circular histograms and rose diagrams  

---

## Papers & Publications

### Journal Article

- Rohatgi, A., Bakulin, A., and Fomel, S., 2025,  
  **Data-driven analysis of seismic phase using circular statistics**:  
  *The Leading Edge*, 44(9), 683–691.

---

### Conference Papers & Expanded Abstracts

- Rohatgi, A., Bakulin, A., and Fomel, S., 2024,  
  **Analyzing the impact of additive and multiplicative noise on seismic data analysis**:  
  SEG International Exposition and Annual Meeting, Expanded Abstracts, SEG-2024-4086176.

- Rohatgi, A., Bakulin, A., and Fomel, S., 2024,  
  **Phase pilot recovery: A foundation for mitigating speckle scattering noise in seismic data**:  
  International Meeting for Applied Geoscience & Energy (IMAGE), Expanded Abstracts, 43.

- Rohatgi, A., Bakulin, A., and Fomel, S., 2025,  
  **Seismic phase spectral analysis: Field-data insights from circular statistics**:  
  International Meeting for Applied Geoscience & Energy (IMAGE), Expanded Abstracts, 44.

- Rohatgi, A., Bakulin, A., and Fomel, S., 2026,  
  **Theory of seismic phase analysis using circular statistics**:  
  International Meeting for Applied Geoscience & Energy (IMAGE), Expanded Abstracts, 44.

- Rohatgi, A., Bakulin, A., Fomel, S., and Badger, J., 2026,  
  **Impact of near-surface topography on reflection distortions: From diffractions to speckle noise**:  
  International Meeting for Applied Geoscience & Energy (IMAGE), Expanded Abstracts, 44.

---

### Other Conference Contributions

- Bakulin, A., Rohatgi, A., and Fomel, S., 2025,  
  **Statistical analysis of seismic phase variability in dense data**:  
  86th EAGE Annual Conference & Exhibition, 1–5.

- Rohatgi, A., Bakulin, A., and Fomel, S., 2024,  
  **Seismic speckle noise: Recognizing scattering noise caused by near-surface heterogeneities**:  
  American Geophysical Union (AGU) Fall Meeting.

- Bakulin, A., Shuster, M., Bhattacharya, S., Delshad, M., Alhotan, M., Li, C., et al., 2024,  
  **Field-test design for geophysical monitoring of hydrogen injection in water-bearing layer**:  
  American Geophysical Union (AGU) Fall Meeting.

