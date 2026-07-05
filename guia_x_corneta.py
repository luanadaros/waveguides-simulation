from pathlib import Path
import matplotlib.pyplot as plt
from SimulacaoFDTD import SimulacaoFDTD

Path("resultados").mkdir(exist_ok=True)

guia = SimulacaoFDTD(
    f0=10e9,
    a=0.020,
    A=0.020,
    L_guia=0.12,
    L_corneta=0,
    L_livre=0.12,
    largura_dominio=0.05,
    usar_corneta=False
)
guia.main_loop(3000)

corneta = SimulacaoFDTD(
    f0=10e9,
    a=0.020,
    A=0.060,
    L_guia=0.12,
    L_corneta=0.10,
    L_livre=0.15,
    largura_dominio=0.12,
    usar_corneta=True
)
corneta.main_loop(3000)

vmax = 1.0  # mesma escala de cor pros dois -> comparação justa
vmin = -vmax

fig, ax = plt.subplots(1, 2, figsize=(13, 5))

im0 = ax[0].imshow(
    guia.Ey,
    origin="lower",
    cmap="RdBu",
    aspect="auto",
    extent=[0, guia.z[-1], 0, guia.x[-1]],
    vmin=vmin, vmax=vmax
)
ax[0].set_title("Guia (terminação abrupta)")
ax[0].set_xlabel("z (m)")
ax[0].set_ylabel("x (m)")
fig.colorbar(im0, ax=ax[0], label="Ey")

im1 = ax[1].imshow(
    corneta.Ey,
    origin="lower",
    cmap="RdBu",
    aspect="auto",
    extent=[0, corneta.z[-1], 0, corneta.x[-1]],
    vmin=vmin, vmax=vmax
)
ax[1].set_title("Corneta (transição gradual)")
ax[1].set_xlabel("z (m)")
ax[1].set_ylabel("x (m)")
fig.colorbar(im1, ax=ax[1], label="Ey")

plt.tight_layout()
plt.savefig(
    "resultados/simulacao4_comparacao.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()