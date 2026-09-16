# Low-Rank Matrix Approximation for Image Compression using SVD

**Undergraduate Thesis Project** — Galatasaray University  
**Major**: Mathematics  

**Author**: Berrïn Karakaya  
**Advisor**: Oğuzhan Kaya  

---

## 📌 Project Overview

This project investigates how tools from the **Spectral Theorem** allow us to analyze high-dimensional matrices and obtain a **low-rank approximation** of a digital image (represented as a matrix) while minimizing information loss relative to the original representation.

The primary objectives are:
1. Review the foundational concepts of the Spectral Theorem.
2. Extend these concepts to rectangular matrices via **Singular Value Decomposition (SVD)**.
3. Apply SVD to digital images to construct approximations using a reduced amount of data, illustrating **dimensionality reduction**.

---

## 📖 Thesis Content

| Section | Content |
|---------|---------|
| **Preliminary Results** | Matrix representation of linear maps, rank, transpose, symmetric and Hermitian matrices, eigenvalues/eigenvectors, characteristic polynomial, algebraic/geometric multiplicities, diagonal matrices, inner product, Euclidean norm, orthogonal and unitary matrices, diagonalizability. |
| **Spectral Theorem** | Orthogonal diagonalization of real symmetric matrices: $A = QDQ^\top$. Proofs of lemmas (real eigenvalues, orthogonality of eigenvectors, linear independence). |
| **Singular Value Decomposition** | Geometry of SVD, Rayleigh Quotient, Low-rank matrix approximation (Eckart-Young theorem). |
| **Application to Image Compression** | Representing an image as a matrix, reconstruction via singular value truncation, empirical examples, and quality analysis. |
| **Conclusion & Appendices** | Summary and supplementary documents. |

---

## 🔬 Core Mathematical Concept

A grayscale image can be viewed as a matrix $A \in \mathbb{R}^{m \times n}$.  
The **SVD** allows us to factorize the matrix as:

$$A = U \Sigma V^\top = \sum_{i=1}^{r} \sigma_i \, u_i v_i^\top$$

where $\sigma_1 \geq \sigma_2 \geq \dots \geq \sigma_r > 0$ are the singular values.

The **rank-$k$ approximation** (where $k \ll r$) is given by:

$$A_k = \sum_{i=1}^{k} \sigma_i \, u_i v_i^\top$$

According to the Eckart-Young theorem, this approximation is optimal with respect to both the Frobenius norm and the spectral $L_2$ norm. It enables storing the image using significantly fewer coefficients while preserving essential visual information.

---

## 📂 Repository Structure

* `code/`: Python scripts for SVD analysis and image visualization (`visualizationdesvd.py`).
* `rapports/`: Complete LaTeX source files (`main.tex`), compiled PDF report (`main_current.pdf`), bibliography reference files (`.bib`), and performance metrics charts (PSNR / SSIM).
* `Journal-de-Bord.md`: Project development log and tracking notes.

---

## 🛠️ Tech Stack & Tools

* **Language:** Python (NumPy, Matplotlib, OpenCV)
* **Typesetting:** LaTeX
* **Domain:** Linear Algebra, Dimensionality Reduction, Image Processing

