# Adaptive Design & Sequential Estimation for Heterogeneous Data
*(Bradley–Terry–Davidson Model)*

This repository contains code associated with the paper:

> **"Adaptive Design and Sequential Estimation for Heterogeneous Data with an Application to Bradley–Terry–Davidson Models"**  

**Key points:**
- *Adaptive sampling* strategies (optimal design, uniform, uncertainty sampling).
- *MLE solvers* for BTD.
- *Soft ranking* loss and \(\ell_2\)-optimal design.


## How to Use
1. **Install** the required R packages (e.g., `igraph`, `nloptr`, `MASS`, `Kendall`, etc.).
2. **Clone** or download this repository.
3. **Open** an R session in the project’s directory.
4. **Run**:
   ```r
   source("R/basic_functions.R")  # core functions
   source("R/simulation_1.R")     # or simulation_2.R
