from pathlib import Path
import matplotlib.pyplot as plt

from SimulacaoFDTD import SimulacaoFDTD

Path("resultados").mkdir(exist_ok=True)

sim = SimulacaoFDTD(
    f0=10e9,
    a=0.020,
    A=0.020,
    L_guia=0.12,
    L_corneta=0,
    L_livre=0,
    largura_dominio=0.05,
    usar_corneta=False
)

frames = sim.main_loop(
    n_passos=2500,
    capturar_em=[2499]
)

Ey = frames[2499]

plt.figure(figsize=(10,4))

plt.imshow(
    Ey.T,
    origin="lower",
    cmap="RdBu",
    aspect="auto",
    extent=[0, sim.L_guia, 0, sim.Nx*sim.dx]
)

plt.colorbar(label="Ey")
plt.xlabel("z (m)")
plt.ylabel("x (m)")
plt.title("Modo TE10 propagando no guia")

plt.tight_layout()

plt.savefig(
    "resultados/simulacao1_prop_guia.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()