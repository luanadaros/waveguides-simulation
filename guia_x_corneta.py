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

fig,ax = plt.subplots(1,2,figsize=(13,5))

ax[0].imshow(
    guia.Ey.T,
    origin="lower",
    cmap="RdBu",
    aspect="auto"
)

ax[0].set_title("Guia")

ax[1].imshow(
    corneta.Ey.T,
    origin="lower",
    cmap="RdBu",
    aspect="auto"
)

ax[1].set_title("Corneta")

plt.tight_layout()

plt.savefig(
    "resultados/simulacao4_comparacao.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()