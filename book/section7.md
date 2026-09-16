# 7. Annexes

## Visualisation de la DVS avec Python

```python
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

A = np.array([[2, 1], [1, 2]])
U, S, Vt = np.linalg.svd(A)
Sigma = np.diag(S)

lim = 4
theta = np.linspace(0, 2*np.pi, 100)
circle = np.array([np.cos(theta), np.sin(theta)])

def setup_ax(ax, title):
    ax.set_aspect('equal')
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axhline(0, color='black', lw=0.5, alpha=0.3, ls='--')
    ax.axvline(0, color='black', lw=0.5, alpha=0.3, ls='--')
    ax.set_title(title, fontsize=10)
    ax.tick_params(labelsize=7)

fig, axs = plt.subplots(2, 2, figsize=(8, 8))
plt.subplots_adjust(wspace=0.5, hspace=0.5)

setup_ax(axs[0, 0], "Cercle unite et directions V")
axs[0, 0].plot(circle[0], circle[1], lw=1.5)
axs[0, 0].quiver(0, 0, Vt[0,0], Vt[0,1], color='red', angles='xy', scale_units='xy', scale=1)
axs[0, 0].quiver(0, 0, Vt[1,0], Vt[1,1], color='green', angles='xy', scale_units='xy', scale=1)
axs[0, 0].text(Vt[0,0]-0.6, Vt[0,1]-0.6, '$v_1$', color='red', weight='bold')
axs[0, 0].text(Vt[1,0]-0.6, Vt[1,1]+0.3, '$v_2$', color='green', weight='bold')

setup_ax(axs[0, 1], "Transformation directe par A")
trans_A = A @ circle
axs[0, 1].plot(trans_A[0], trans_A[1], color='red', lw=1.5)
v1f, v2f = A @ Vt[0], A @ Vt[1]
axs[0, 1].quiver(0, 0, v1f[0], v1f[1], color='red', angles='xy', scale_units='xy', scale=1)
axs[0, 1].quiver(0, 0, v2f[0], v2f[1], color='green', angles='xy', scale_units='xy', scale=1)
axs[0, 1].text(v1f[0]-1.0, v1f[1]-0.2, r'$\sigma_1 u_1$', color='red')
axs[0, 1].text(v2f[0]-1.0, v2f[1]+0.4, r'$\sigma_2 u_2$', color='green')

setup_ax(axs[1, 0], "Rotation par $V^T$")
trans_Vt = Vt @ circle
axs[1, 0].plot(trans_Vt[0], trans_Vt[1], lw=1.5)
axs[1, 0].quiver(0, 0, 1, 0, color='red', angles='xy', scale_units='xy', scale=1)
axs[1, 0].quiver(0, 0, 0, 1, color='green', angles='xy', scale_units='xy', scale=1)
axs[1, 0].text(1.2, 0.2, '$e_1$', color='red')
axs[1, 0].text(0.2, 1.2, '$e_2$', color='green')

setup_ax(axs[1, 1], "Etirement par $\Sigma$")
trans_Sigma = Sigma @ trans_Vt
axs[1, 1].plot(trans_Sigma[0], trans_Sigma[1], lw=1.5)
axs[1, 1].quiver(0, 0, S[0], 0, color='red', angles='xy', scale_units='xy', scale=1)
axs[1, 1].quiver(0, 0, 0, S[1], color='green', angles='xy', scale_units='xy', scale=1)
axs[1, 1].text(S[0]+0.2, 0.2, r'$\sigma_1 e_1$', color='red')
axs[1, 1].text(0.2, S[1]+0.2, r'$\sigma_2 e_2$', color='green')

plt.show()
```

## Chargement des images et calcul du rang matriciel

```python
from PIL import Image
import numpy as np
img1 = Image.open("chat_simple.jpg").convert('L')
img2 = Image.open("chat_detaille.jpg").convert('L')
print("Tailles originales :")
print(img1.size, img2.size)
target_size = (
    min(img1.size[0], img2.size[0]),
    min(img1.size[1], img2.size[1])
)
img1 = img1.resize(target_size, Image.LANCZOS)
img2 = img2.resize(target_size, Image.LANCZOS)
print("Taille commune :", target_size)
mat1 = np.array(img1, dtype=np.float64)
mat2 = np.array(img2, dtype=np.float64)
print("Shapes :", mat1.shape, mat2.shape)
rank1 = np.linalg.matrix_rank(mat1)
rank2 = np.linalg.matrix_rank(mat2)
print("Rank chat simple :", rank1)
print("Rank chat detaille :", rank2)
print("Rank max theorique :", min(mat1.shape))
```

## Calcul de l'energie conservee via DVS

```python
from PIL import Image
import numpy as np
image1_path = "chat_simple.jpg"
image2_path = "chat_detaille.jpg"
img1 = Image.open(image1_path).convert('L')
img2 = Image.open(image2_path).convert('L')
target_size = (
    min(img1.size[0], img2.size[0]),
    min(img1.size[1], img2.size[1])
)
img1 = img1.resize(target_size, Image.LANCZOS)
img2 = img2.resize(target_size, Image.LANCZOS)
mat1 = np.array(img1, dtype=np.float64)
mat2 = np.array(img2, dtype=np.float64)
print("Shape :", mat1.shape)
def apply_svd_and_print(matrix, image_name):
    U, S, Vt = np.linalg.svd(matrix, full_matrices=False)
    total_energy = np.sum(S**2)
    print(f"\nImage : {image_name}")
    print(f"Singular values : {len(S)}")
    print("\nTop 20 singular values")
    for i in range(20):
        print(f"sigma{i+1} = {S[i]:.2f}")
    print("\nEnergy retention (%)")
    for k in [5, 10, 20, 50, 100]:
        energy = np.sum(S[:k]**2) / total_energy * 100
        print(f"k={k:3d} -> {energy:6.2f}%")
    return S
print("DVS analysis")
S1 = apply_svd_and_print(mat1, "chat_simple.jpg")
S2 = apply_svd_and_print(mat2, "chat_detaille.jpg")
```

## Analyse de l'energie, erreur de Frobenius et PSNR via DVS

```python
from PIL import Image
import numpy as np
image1_path = "chat_simple.jpg"
image2_path = "chat_detaille.jpg"
img1 = Image.open(image1_path).convert('L')
img2 = Image.open(image2_path).convert('L')
target_size = (
    min(img1.size[0], img2.size[0]),
    min(img1.size[1], img2.size[1])
)
img1 = img1.resize(target_size, Image.LANCZOS)
img2 = img2.resize(target_size, Image.LANCZOS)
mat1 = np.array(img1, dtype=np.float64)
mat2 = np.array(img2, dtype=np.float64)
print("Shape :", mat1.shape)
def analyze_svd(matrix, image_name):
    U, S, Vt = np.linalg.svd(matrix, full_matrices=False)
    total_energy = np.sum(S**2)
    total_fro_norm = np.sqrt(total_energy)
    print(f"\nImage : {image_name}")
    print("k | Energy (%) | Frobenius error | PSNR (dB)")
    print("-" * 70)
    for k in [5, 10, 20, 50, 100]:
        energy = np.sum(S[:k]**2) / total_energy * 100
        fro_error = np.sqrt(np.sum(S[k:]**2)) / total_fro_norm
        reconstructed = (
            U[:, :k]
            @ np.diag(S[:k])
            @ Vt[:k, :]
        )
        mse = np.mean((matrix - reconstructed)**2)
        if mse > 0:
            psnr = 20 * np.log10(255 / np.sqrt(mse))
        else:
            psnr = float("inf")
        print(
            f"{k:3d} | "
            f"{energy:10.2f}% | "
            f"{fro_error:14.4f} | "
            f"{psnr:9.2f}"
        )
    print("-" * 70)
print("DVS analysis start")
analyze_svd(mat1, "chat_simple.jpg")
analyze_svd(mat2, "chat_detaille.jpg")
```

## Approximation de rang k et compression DVS

```python
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
image_path = "chat_simple.jpg"
original_img = Image.open(image_path).convert('L')
img = np.array(original_img, dtype=np.float64)
print("Matrix shape :", img.shape)
print("=" * 70)
print("APPROXIMATION DE RANG k (DVS)")
print("=" * 70)
print("A\approx A_k = U_k Sigma_k V_k^T")
print("U_k : k singular vectors (left)")
print("Sigma_k : diagonal matrix of k singular values")
print("V_k^T : k singular vectors (right)")
print("Eckart-Young theorem: best rank-k approximation")
print("=" * 70)
U, S, Vt = np.linalg.svd(img, full_matrices=False)
ks = [5, 10, 20, 50, 100]
fig, axes = plt.subplots(2, 3, figsize=(10, 7))
fig.suptitle("SVD rank-k approximations", fontsize=14)
axes[0, 0].imshow(original_img, cmap='gray')
axes[0, 0].set_title('Original')
axes[0, 0].axis('off')
for i, k in enumerate(ks):
    approx = U[:, :k] @ np.diag(S[:k]) @ Vt[:k, :]
    approx = np.clip(approx, 0, 255).astype(np.uint8)
    row = (i + 1) // 3
    col = (i + 1) % 3
    axes[row, col].imshow(approx, cmap='gray')
    axes[row, col].set_title(f'k = {k}')
    axes[row, col].axis('off')
plt.tight_layout()
plt.savefig("chat_simple_svd_compression.png", dpi=300, bbox_inches='tight')
plt.show()
print("Saved: chat_simple_svd_compression.png")
```

## Approximation de rang k par DVS (image detaillee)

```python
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
image_path = "chat_detaille.jpg"
original_img = Image.open(image_path).convert('L')
img = np.array(original_img, dtype=np.float64)
print("Matrix shape :", img.shape)
print("=" * 70)
print("RANK-k MATRIX APPROXIMATION (DVS)")
print("=" * 70)
print("A approx A_k = U_k Sigma_k V_k^T")
print("U_k : first k left singular vectors")
print("Sigma_k : diagonal matrix of k singular values")
print("V_k^T : first k right singular vectors")
print("Eckart-Young theorem: optimal rank-k approximation")
print("=" * 70)
U, S, Vt = np.linalg.svd(img, full_matrices=False)
ks = [5, 10, 20, 50, 100]
fig, axes = plt.subplots(2, 3, figsize=(12, 7))
fig.suptitle("SVD rank-k approximations (chat detaille)", fontsize=14)
axes[0, 0].imshow(original_img, cmap='gray')
axes[0, 0].set_title('Original image')
axes[0, 0].axis('off')
for i, k in enumerate(ks):
    approx = U[:, :k] @ np.diag(S[:k]) @ Vt[:k, :]
    approx = np.clip(approx, 0, 255).astype(np.uint8)
    row = (i + 1) // 3
    col = (i + 1) % 3
    axes[row, col].imshow(approx, cmap='gray')
    axes[row, col].set_title(f'k = {k}')
    axes[row, col].axis('off')
plt.tight_layout()
plt.savefig("chat_detaille_svd_compression.png", dpi=300, bbox_inches='tight')
plt.show()
print("Saved: chat_detaille_svd_compression.png")
```

## Analyse DVS : energie, erreur Frobenius et PSNR

```python
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
image_path = "chat_detaille.jpg"
img = Image.open(image_path).convert('L')
mat = np.array(img, dtype=np.float64)
print("Matrix shape :", mat.shape)
U, S, Vt = np.linalg.svd(mat, full_matrices=False)
ks = np.arange(0, 201, 5)
energy = np.zeros(len(ks))
fro_error = np.zeros(len(ks))
psnr_values = []
total_energy = np.sum(S**2)
for i, k in enumerate(ks):
    if k == 0:
        energy[i] = 0.0
        fro_error[i] = 1.0
        psnr_values.append(0.0)
    else:
        cum_energy = np.sum(S[:k]**2)
        energy[i] = cum_energy / total_energy
        fro_error[i] = np.sqrt(1 - energy[i])
        recon = U[:, :k] @ np.diag(S[:k]) @ Vt[:k, :]
        mse = np.mean((mat - recon)**2)
        psnr = 20 * np.log10(255 / np.sqrt(mse)) if mse > 0 else 100
        psnr_values.append(psnr)
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
fig.suptitle("SVD analysis - chat detaille")
axes[0].plot(S[:200], linewidth=1.5)
axes[0].set_title("Singular values (sigma)")
axes[0].set_xlabel("k")
axes[0].set_ylabel("sigma")
axes[0].grid(True)
axes[1].plot(ks, fro_error * 100, linewidth=2)
axes[1].set_title("Frobenius error (%)")
axes[1].set_xlabel("k")
axes[1].set_ylabel("error (%)")
axes[1].grid(True)
plt.tight_layout()
plt.savefig("chat_detaille_svd_analysis_corrected.png", dpi=300, bbox_inches='tight')
plt.show()
print("Saved: chat_detaille_svd_analysis_corrected.png")
```

## Analyse DVS : image simple

```python
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
image_path = "chat_simple.jpg"
img = Image.open(image_path).convert('L')
mat = np.array(img, dtype=np.float64)
print("Matrix shape :", mat.shape)
U, S, Vt = np.linalg.svd(mat, full_matrices=False)
ks = np.arange(0, 201, 5)
energy = np.zeros(len(ks))
fro_error = np.zeros(len(ks))
psnr_values = []
total_energy = np.sum(S**2)
for i, k in enumerate(ks):
    if k == 0:
        energy[i] = 0.0
        fro_error[i] = 1.0
        psnr_values.append(0.0)
    else:
        cum_energy = np.sum(S[:k]**2)
        energy[i] = cum_energy / total_energy
        fro_error[i] = np.sqrt(1 - energy[i])
        recon = U[:, :k] @ np.diag(S[:k]) @ Vt[:k, :]
        mse = np.mean((mat - recon)**2)
        psnr = 20 * np.log10(255 / np.sqrt(mse)) if mse > 0 else 100
        psnr_values.append(psnr)
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
fig.suptitle("SVD analysis - chat simple")
axes[0].plot(S[:200], linewidth=1.5)
axes[0].set_title("Singular values (sigma)")
axes[0].set_xlabel("k")
axes[0].set_ylabel("sigma")
axes[0].grid(True)
axes[1].plot(ks, fro_error * 100, linewidth=2)
axes[1].set_title("Frobenius error (%)")
axes[1].set_xlabel("k")
axes[1].set_ylabel("error (%)")
axes[1].grid(True)
plt.tight_layout()
plt.savefig("chat_simple_svd_analysis_corrected.png", dpi=300, bbox_inches='tight')
plt.show()
print("Saved: chat_simple_svd_analysis_corrected.png")
```

## Comparaison PSNR et SSIM via DVS

```python
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from skimage.metrics import structural_similarity as ssim
print("=" * 80)
print("PSNR AND SSIM DEFINITIONS")
print("=" * 80)
print("PSNR: pixel-based error metric")
print("PSNR = 20 * log10(255 / sqrt(MSE))")
print("Higher PSNR = better reconstruction\n")
print("SSIM: structural similarity metric")
print("Measures luminance, contrast, structure")
print("Range: 0 to 1 (1 = identical image)")
print("=" * 80)
def compute_metrics(image_path, name):
    original = Image.open(image_path).convert('L')
    orig_array = np.array(original, dtype=np.float64)
    U, S, Vt = np.linalg.svd(orig_array, full_matrices=False)
    ks = [5, 10, 20, 50, 100]
    psnr_list = []
    ssim_list = []
    orig_uint8 = np.clip(orig_array, 0, 255).astype(np.uint8)
    print(f"\nResults for {name}")
    print("k | PSNR | SSIM")
    print("-" * 40)
    for k in ks:
        recon = U[:, :k] @ np.diag(S[:k]) @ Vt[:k, :]
        recon = np.clip(recon, 0, 255).astype(np.uint8)
        mse = np.mean((orig_array - recon)**2)
        psnr = 20 * np.log10(255 / np.sqrt(mse)) if mse > 0 else 100
        ssim_val = ssim(orig_uint8, recon, data_range=255)
        psnr_list.append(psnr)
        ssim_list.append(ssim_val)
        print(f"{k:3d} | {psnr:8.2f} | {ssim_val:.4f}")
    return ks, psnr_list, ssim_list
ks, psnr_detail, ssim_detail = compute_metrics("chat_detaille.jpg", "detaille")
ks, psnr_simple, ssim_simple = compute_metrics("chat_simple.jpg", "simple")
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
fig.suptitle("PSNR and SSIM comparison")
axes[0].plot(ks, psnr_detail, 'r-o', label='detaille')
axes[0].plot(ks, psnr_simple, 'b-o', label='simple')
axes[0].set_title('PSNR')
axes[0].set_xlabel('k')
axes[0].set_ylabel('PSNR')
axes[0].legend()
axes[0].grid(True)
axes[1].plot(ks, ssim_detail, 'r-o', label='detaille')
axes[1].plot(ks, ssim_simple, 'b-o', label='simple')
axes[1].set_title('SSIM')
axes[1].set_xlabel('k')
axes[1].set_ylabel('SSIM')
axes[1].legend()
axes[1].grid(True)
plt.tight_layout()
plt.savefig("psnr_ssim_comparison.png", dpi=300, bbox_inches='tight')
plt.show()
print("Saved: psnr_ssim_comparison.png")
```
