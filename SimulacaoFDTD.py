import numpy as np
import matplotlib.pyplot as plt
import imageio
from pathlib import Path

# parametros fisicos vacuo
EPS0 = 8.8541878128e-12       # permissividade
MU0 = 4 * np.pi * 1e-7        # permeabilidade
C0 = 1.0 / np.sqrt(EPS0 * MU0)  # velocidade da luz

class SimulacaoFDTD:
    def __init__(self, f0, a, A, L_guia, L_corneta, L_livre,
                 largura_dominio, espessura_parede_cel=2,
                 pontos_por_lambda=20, margem_cfl=0.99, usar_corneta=True):
        """
        f0 : Frequência da fonte [Hz]
        a : Largura do guia de onda [m]
        A : Largura da boca da corneta [m] (A > a)
        L_guia : Comprimento da seção reta do guia [m].
        L_corneta : Comprimento da seção da corneta [m].
        L_livre : Comprimento da região de espaço livre após a corneta [m].
        largura_dominio : Extensão total do domínio na direção x [m].
        espessura_parede_cel : Espessura da parede metálica, em número de células.
        pontos_por_lambda : Quantidade de células por comprimento de onda mínimo simulado.
        margem_cfl : Fator de segurança (<1) aplicado ao critério CFL.
        usar_corneta : Se False, simula apenas o guia reto (usado para o estudo do modo evanescente).
        """
        self.f0 = f0
        self.omega = 2 * np.pi * f0
        self.a = a
        self.A = A if usar_corneta else a
        self.L_guia = L_guia
        self.L_corneta = L_corneta if usar_corneta else 0.0
        self.L_livre = L_livre
        self.usar_corneta = usar_corneta

        # frequência de corte do modo dominante TE10
        self.fc10 = C0 / (2 * a)

        # dx, dz, dt, Nx, Nz
        lambda_min = C0 / f0
        self.dx = lambda_min / pontos_por_lambda # precisão da simulação
        self.dz = self.dx  # dx = dz = delta (grade quadrada)
        self.Nx = int(round(largura_dominio / self.dx)) + 1
        self.Nz = int(round((L_guia + self.L_corneta + L_livre) / self.dz)) + 1

        # critério de estabilidade CFL 
        self.dt = margem_cfl * self.dx / (C0 * np.sqrt(2))

        # posição central do guia no eixo x
        self.x = np.arange(self.Nx) * self.dx
        self.z = np.arange(self.Nz) * self.dz
        self.xc = largura_dominio / 2.0


        # Ey definido nos nós inteiros (Nx, Nz)
        # Hx definido com deslocamento em z (Nx, Nz-1)
        # Hz definido com deslocamento em x (Nx-1, Nz)
        self.Ey = np.zeros((self.Nx, self.Nz))
        self.Hx = np.zeros((self.Nx, self.Nz - 1))
        self.Hz = np.zeros((self.Nx - 1, self.Nz))

        # máscara PEC
        self.espessura_parede_cel = espessura_parede_cel
        self.mask_pec, self.iz_boca = self._construir_mascara_pec()

        # posição da fonte em z
        self.iz_fonte = max(3, int(round(0.15 * L_guia / self.dz))) 
        ix_low, ix_up = self._paredes_em_z(self.iz_fonte)
        self.ix_fonte = np.arange(ix_low + 1, ix_up)  # células internas ao guia
        self.x_fonte = self.x[self.ix_fonte]
        self.x_low_fonte = self.x[ix_low]

        # Amplitude da fonte e rampa suave de ligamento (evita transiente forte)
        self.E0 = 1.0
        self.n_rampa = int(round(3 * (1 / f0) / self.dt))  # ~3 períodos

    def _semi_largura(self, z_val):
        if z_val <= self.L_guia:
            return self.a / 2.0
        elif z_val <= self.L_guia + self.L_corneta and self.usar_corneta:
            t = (z_val - self.L_guia) / self.L_corneta
            return self.a / 2.0 + (self.A / 2.0 - self.a / 2.0) * t
        else:
            return None  # fora do guia/corneta -> sem parede (espaço livre)

    def _paredes_em_z(self, iz):
        semi = self._semi_largura(self.z[iz])
        if semi is None:
            return None, None
        ix_inf = int(round((self.xc - semi) / self.dx))
        ix_sup = int(round((self.xc + semi) / self.dx))
        return ix_inf, ix_sup

    def _construir_mascara_pec(self):
        mask = np.zeros((self.Nx, self.Nz), dtype=bool)
        iz_boca = self.Nz - 1
        for iz in range(self.Nz):
            ix_inf, ix_sup = self._paredes_em_z(iz)
            if ix_inf is None:
                if iz_boca == self.Nz - 1:
                    iz_boca = iz  # primeiro ponto após a boca da corneta
                continue
            t = self.espessura_parede_cel
            for d in range(t):
                if 0 <= ix_inf - d < self.Nx:
                    mask[ix_inf - d, iz] = True
                if 0 <= ix_sup + d < self.Nx:
                    mask[ix_sup + d, iz] = True
        return mask, iz_boca

    def _injetar_fonte(self, n):
        t = n * self.dt
        rampa = min(1.0, t / (self.n_rampa * self.dt)) if self.n_rampa > 0 else 1.0
        perfil_x = np.sin(np.pi * (self.x_fonte - self.x_low_fonte) / self.a)
        self.Ey[self.ix_fonte, self.iz_fonte] = (
            self.E0 * perfil_x * np.sin(self.omega * t) * rampa
        )

    def _mur_abc(self, Ey_ant):
        coef = (C0 * self.dt - self.dx) / (C0 * self.dt + self.dx)

        # bordas em x (esquerda / direita)
        self.Ey[0, :] = Ey_ant[1, :] + coef * (self.Ey[1, :] - Ey_ant[0, :])
        self.Ey[-1, :] = Ey_ant[-2, :] + coef * (self.Ey[-2, :] - Ey_ant[-1, :])

        # bordas em z (início / fim do domínio)
        self.Ey[:, 0] = Ey_ant[:, 1] + coef * (self.Ey[:, 1] - Ey_ant[:, 0])
        self.Ey[:, -1] = Ey_ant[:, -2] + coef * (self.Ey[:, -2] - Ey_ant[:, -1])

    def passo(self, n):
        # atualizar Hx e Hz
        self.Hx -= (self.dt / (MU0 * self.dz)) * (self.Ey[:, 1:] - self.Ey[:, :-1])
        self.Hz += (self.dt / (MU0 * self.dx)) * (self.Ey[1:, :] - self.Ey[:-1, :])

        # atualizar Ey (pontos internos)
        Ey_ant = self.Ey.copy()
        self.Ey[1:-1, 1:-1] += (self.dt / EPS0) * (
            (self.Hz[1:, 1:-1] - self.Hz[:-1, 1:-1]) / self.dx
            - (self.Hx[1:-1, 1:] - self.Hx[1:-1, :-1]) / self.dz
        )

        # (fronteiras absorventes, calculadas antes da injeção de fonte/PEC)
        self._mur_abc(Ey_ant)

        # injetar fonte TE10
        self._injetar_fonte(n)

        # aplicar máscara PEC (zerar Ey nas paredes)
        self.Ey[self.mask_pec] = 0.0

    def main_loop(self, n_passos, capturar_em=None):
        frames = {}
        capturar_em = set(capturar_em or [])
        for n in range(n_passos):
            self.passo(n)
            if n in capturar_em:
                frames[n] = self.Ey.copy()
        return frames

    def info(self):
        return (
            f"f0 = {self.f0/1e9:.2f} GHz | fc10 = {self.fc10/1e9:.2f} GHz "
            f"({'propagante' if self.f0 > self.fc10 else 'EVANESCENTE'})\n"
            f"dx = dz = {self.dx*1e3:.3f} mm | dt = {self.dt*1e12:.3f} ps\n"
            f"Nx = {self.Nx} | Nz = {self.Nz}"
        )

    def salvar_gif(self, frames, path):
        pasta = Path("resultados/tmp")
        pasta.mkdir(parents=True, exist_ok=True)

        imagens = []

        for passo, Ey in frames.items():

            fig, ax = plt.subplots(figsize=(10,5))

            im = ax.imshow(
                Ey.T,
                origin="lower",
                aspect="auto",
                cmap="RdBu",
                extent=[0, self.z[-1], 0, self.x[-1]]
            )

            plt.colorbar(im, ax=ax, label="Campo elétrico Ey")

            ax.set_xlabel("z (m)")
            ax.set_ylabel("x (m)")
            ax.set_title(f"t = {passo*self.dt*1e9:.2f} ns")

            arquivo = pasta / f"{passo:05d}.png"

            plt.savefig(
                arquivo,
                dpi=120,
            )

            plt.close()

            imagens.append(imageio.imread(arquivo))

        imageio.mimsave(
            f"resultados/{path}.gif",
            imagens,
            duration=0.08
        )