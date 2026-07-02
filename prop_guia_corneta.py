from pathlib import Path
import matplotlib.pyplot as plt

from SimulacaoFDTD import SimulacaoFDTD

Path("resultados").mkdir(exist_ok=True)

sim = SimulacaoFDTD(
    f0=10e9,
    a=0.020,
    A=0.060,
    L_guia=0.12,
    L_corneta=0.10,
    L_livre=0.15,
    largura_dominio=0.12,
    usar_corneta=True
)

frames = sim.main_loop(
    n_passos=4000,
    capturar_em=[3999]
)

Ey = frames[3999]

plt.figure(figsize=(12,5))

plt.imshow(
    Ey.T,
    origin="lower",
    cmap="RdBu",
    aspect="auto",
    extent=[0, sim.z[-1], 0, sim.x[-1]]
)

plt.colorbar(label="Ey")
plt.xlabel("z (m)")
plt.ylabel("x (m)")
plt.title("Transição Guia → Corneta → Espaço Livre")

plt.tight_layout()

plt.savefig(
    "resultados/simulacao2_prop_guia_corneta.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()