from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from SimulacaoFDTD import SimulacaoFDTD

Path("resultados").mkdir(exist_ok=True)

sim = SimulacaoFDTD(
    f0=5e9,
    a=0.020,
    A=0.020,
    L_guia=0.20,
    L_corneta=0,
    L_livre=0,
    largura_dominio=0.05,
    usar_corneta=False
)

sim.main_loop(3500)

ix = sim.Nx//2

Ey = np.abs(sim.Ey[ix,:])

plt.figure(figsize=(8,4))

plt.plot(sim.z,Ey)

plt.grid(True)

plt.xlabel("z (m)")
plt.ylabel("|Ey|")
plt.title("Modo Evanescente")

plt.tight_layout()

plt.savefig(
    "resultados/simulacao3_prop_evanescente.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()