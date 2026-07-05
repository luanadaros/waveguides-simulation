# Modo evanescente em guia de onda retangular
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

ix = sim.Nx // 2
T = 1 / sim.f0
n_por_periodo = int(round(T / sim.dt))

n_pre = 3000  # deixa o regime estacionário se estabelecer
for n in range(n_pre):
    sim.passo(n)

# envelope = maior |Ey| ao longo de 1 período completo, em cada z
envelope = np.zeros(sim.Nz)
for k in range(n_por_periodo):
    sim.passo(n_pre + k)
    envelope = np.maximum(envelope, np.abs(sim.Ey[ix, :]))

plt.figure(figsize=(8, 4))
plt.semilogy(sim.z * 100, envelope)  # log no eixo y evidencia a exponencial
plt.grid(True, which="both")
plt.xlabel("z (cm)")
plt.ylabel("|Ey| (envelope, escala log)")
plt.title(f"Modo evanescente — f0={sim.f0/1e9:.1f} GHz < fc10={sim.fc10/1e9:.1f} GHz")
plt.tight_layout()
plt.savefig(
    "resultados/simulacao3_prop_evanescente.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()