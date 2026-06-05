# FRB_pdf_LSS

This repository implements the computation of the one-point probability distribution function (PDF) of dispersion measure (DM) within the baryonification framework(https://arxiv.org/pdf/2507.07892). The relevant formulas and theoretical background can be found in (https://arxiv.org/pdf/2601.18784).

---

## 📁 Repository Structure

The repository contains two main components:

### 1. `pdf`
A full implementation of the code, optimized for execution on the Marvin HPC system(https://wiki.hpc.uni-bonn.de/).

This version is designed for high-performance computation and large-scale runs, but can be adapted to other computing environments if needed.

#### Running the code

1. Configure the desired BFC parameters in `pipeline.py`.
2. Submit the job to the scheduler:

```bash
sbatch start.sh
```

---

### 2. `Test_notebook`
A lightweight Jupyter notebook version of the `pdf` code, intended for execution on a personal computer.

- Reduced numerical resolution for faster runtime  
- Not intended for precision scientific analysis  
- Useful for understanding the structure of the code, functions, and numerical methods  

---

---

### 3. `Emulator`
This emulator was used to run the MCMC analysis and fit the PDFs.

Users should be aware that the emulator may be less accurate at very low DM values.

---

## ⚙️ Requirements

To run this code, the **BFC code** is required as a prerequisite.

---

## 📌 Notes

- The HPC version (`pdf`) is intended for production-level computations.  
- The notebook version is intended for learning, testing, and exploration.  
- Results from the notebook should be treated as qualitative rather than quantitative.

---

## 📚 References

- Baryonification framework: https://arxiv.org/pdf/2507.07892  
- Theoretical background: https://arxiv.org/pdf/2601.18784  

---
