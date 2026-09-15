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
