import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

from SimulacaoFDTD import SimulacaoFDTD, C0

OUT = "output"
os.makedirs(OUT, exist_ok=True)

# modo propagante
f0 = 10e9          # 10 GHz
a = 0.02           # largura do guia: 2 cm  -> fc10 = 7.5 GHz  (f0 > fc10 -> propaga)
A = 0.06           # boca da corneta: 6 cm
L_guia = 0.04
L_corneta = 0.06
L_livre = 0.06
largura_dominio = 0.10

sim = SimulacaoFDTD(
    f0=f0, a=a, A=A, L_guia=L_guia, L_corneta=L_corneta, L_livre=L_livre,
    largura_dominio=largura_dominio, pontos_por_lambda=20,
)
print(sim.info())

n_passos = 1600

passos = list(range(0, n_passos, 8))
frames = sim.main_loop(n_passos, capturar_em=passos)

vmax = 0.9 * sim.E0
extent = [0, sim.L_guia + sim.L_corneta + sim.L_livre, 0, largura_dominio] # z (m), x (m)

def plot_frame(Ey, titulo, caminho, marcar_regioes=True):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    im = ax.imshow(
        Ey.T, origin='lower', extent=extent, aspect='auto',
        cmap='RdBu_r', norm=mcolors.Normalize(vmin=-vmax, vmax=vmax),
    )
    # Contorno das paredes PEC
    parede = np.ma.masked_where(~sim.mask_pec, np.ones_like(sim.mask_pec, dtype=float))
    ax.imshow(parede.T, origin='lower', extent=extent, aspect='auto', cmap='gray_r', vmin=0, vmax=1, alpha=0.9)

    if marcar_regioes:
        ax.axvline(L_guia, color='k', ls='--', lw=0.8)
        ax.axvline(L_guia + L_corneta, color='k', ls='--', lw=0.8)
        ax.text(L_guia/2, largura_dominio*0.97, 'GUIA', ha='center', va='top', fontsize=9)
        ax.text(L_guia + L_corneta/2, largura_dominio*0.97, 'CORNETA', ha='center', va='top', fontsize=9)
        ax.text(L_guia + L_corneta + L_livre/2, largura_dominio*0.97, 'ESPAÇO LIVRE', ha='center', va='top', fontsize=9)

    ax.set_xlabel('z (m) - direção de propagação')
    ax.set_ylabel('x (m) - direção transversal')
    ax.set_title(titulo)
    fig.colorbar(im, ax=ax, label='Ey (u.a.)')
    fig.tight_layout()
    fig.savefig(caminho, dpi=150)
    plt.close(fig)

passos_ordenados = sorted(frames.keys())

# passo em que a frente de onda está no guia (início da simulação)
n_guia = passos_ordenados[len(passos_ordenados)//6]
# passo em que a onda já atravessou a corneta e entra no espaço livre
n_transicao = passos_ordenados[len(passos_ordenados)//2]
# passo final: campo já irradiado no espaço livre (regime quase estacionário)
n_livre = passos_ordenados[-1]

plot_frame(frames[n_guia], f'Campo Ey — onda ainda no guia (passo {n_guia})',
           f'{OUT}/1_campo_no_guia.png')
plot_frame(frames[n_transicao], f'Campo Ey — transição na corneta (passo {n_transicao})',
           f'{OUT}/2_campo_na_corneta.png')
plot_frame(frames[n_livre], f'Campo Ey — irradiação no espaço livre (passo {n_livre})',
           f'{OUT}/3_campo_espaco_livre.png')

print("Figuras estáticas salvas.")
