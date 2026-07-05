# Relatório — Simulação Computacional FDTD 2D de Guia de Onda e Antena Corneta

Este README documenta o que compõe o relatório exigido na **Seção 6** do enunciado
(*Trabalho Prático: Simulação Computacional de Guias de Onda e Antenas Corneta*),
descrevendo os quatro itens obrigatórios e como cada um é atendido pelos scripts
deste projeto.

## Sumário

- [1. Apresentação Visual](#1-apresentação-visual)
- [2. Análise de Transição](#2-análise-de-transição)
- [3. Estudo de Caso — Modo Evanescente](#3-estudo-de-caso--modo-evanescente)
- [4. Código Fonte](#4-código-fonte)
- [Estrutura do repositório](#estrutura-do-repositório)
- [Como reproduzir os resultados](#como-reproduzir-os-resultados)
- [Parâmetros físicos utilizados](#parâmetros-físicos-utilizados)

---

## 1. Apresentação Visual

**Exigência do enunciado:** mapas de calor (amplitude de `Ey`) no guia, na corneta e no espaço livre.

Gerados por `main.py`, que produz três figuras estáticas em `output/`:

| Arquivo | Conteúdo |
|---|---|
| `1_campo_no_guia.png` | Campo `Ey` capturado quando a frente de onda ainda está dentro da seção reta do guia |
| `2_campo_na_corneta.png` | Campo capturado durante a transição, com a onda já dentro da corneta |
| `3_campo_espaco_livre.png` | Campo em regime já estabelecido, com irradiação visível no espaço livre |

**Detalhe importante de implementação:** os três instantes de tempo (`n_guia`,
`n_transicao`, `n_livre`) **não são escolhidos arbitrariamente** — são calculados a
partir do tempo físico real de trânsito da onda até cada marco geométrico, usando a
**velocidade de grupo** do modo TE10:

```
vg = c · √(1 − (fc10/f0)²)
```

em vez da velocidade da luz `c`. Isso é essencial porque, dentro do guia, a energia
se propaga a `vg ≈ 66%` de `c` para os parâmetros usados aqui — se o tempo de trânsito
for estimado com `c`, os frames escolhidos ficam adiantados demais e mostram o campo
já em regime estabelecido em todos os três "estágios", em vez de capturar de fato o
guia, a transição e o espaço livre separadamente. O cálculo também soma o tempo de
rampa de acionamento da fonte (~3 períodos) e considera a posição real da fonte
(que não fica em `z=0`).

## 2. Análise de Transição

**Exigência do enunciado:** explicar a mudança da frente de onda (plana → cilíndrica)
e sua relação com o casamento de impedância promovido pela corneta.

Esta parte é **textual**, redigida com base na comparação visual entre as três figuras
acima (e opcionalmente `simulacao4_comparacao.png`, que compara lado a lado um guia
com terminação abrupta e um guia com corneta). Pontos a desenvolver no texto do
relatório:

- Dentro do guia, a frente de onda é aproximadamente **plana**, com o perfil
  transversal senoidal do modo TE10 mantido constante ao longo de `z`.
- Ao entrar na corneta, a seção transversal cresce gradualmente (afunilamento
  linear de `a` até `A`), e a frente de onda passa a se curvar, aproximando-se de
  uma **frente cilíndrica/esférica** conforme se aproxima da boca.
- Essa transição gradual — em vez de uma abertura abrupta — funciona como um
  **casamento de impedância**: reduz reflexões na transição guia→espaço livre,
  permitindo que mais energia seja efetivamente irradiada em vez de refletida de
  volta para o guia. A comparação entre `1_campo_no_guia`/`simulacao1` (terminação
  abrupta) e `simulacao2`/`2_campo_na_corneta` (transição suave) evidencia esse
  efeito.
- Pode-se citar a impedância característica do modo TE10 no guia,
  `Zte = η₀ / √(1 − (fc10/f0)²)` (onde `η₀ ≈ 377 Ω`), que difere da impedância do
  espaço livre; a corneta faz essa impedância variar suavemente entre os dois
  valores ao longo do seu comprimento.

## 3. Estudo de Caso — Modo Evanescente

**Exigência do enunciado:** simular com `f < fc10` e demonstrar visualmente e
matematicamente o decaimento exponencial sem propagação.

Gerado por `simulacao_evanescente.py`, com `f0 = 5 GHz < fc10 = 7,5 GHz`
(`a = 0,02 m`), guia reto sem corneta (`usar_corneta=False`). Saída:
`simulacao3_prop_evanescente.png`.

**Metodologia:** em vez de capturar `|Ey(z)|` num único instante de tempo arbitrário
(o que produz uma curva "ruidosa" com ondulações de fase), o script:

1. deixa o regime se estabelecer (~3000 passos de "aquecimento");
2. varre um período completo de oscilação, guardando o **envelope**
   (`np.maximum` do módulo do campo) em cada posição `z`;
3. plota o envelope em **escala logarítmica** (`plt.semilogy`).

Numa onda evanescente pura, o log da amplitude decai **linearmente** com `z` — por
isso a escala log é a evidência visual mais forte de decaimento exponencial
(um trecho reto no gráfico log-linear).

**Validação quantitativa** (não obrigatória, mas fortalece o relatório): medindo a
inclinação do trecho reto e comparando com a constante de atenuação teórica

```
α = (2π/c) · √(fc10² − f0²)
```

o valor medido no FDTD (~103 Np/m) ficou em torno de **88% do valor teórico**
(~117 Np/m), diferença esperada de dispersão numérica do método FDTD, especialmente
próximo à frequência de corte. Após um certo `z`, o gráfico atinge um patamar de
**ruído numérico de fundo** (não um artefato físico) — vale mencionar isso
explicitamente no texto para deixar claro que não é um erro de implementação.

## 4. Código Fonte

**Exigência do enunciado:** código Python organizado e comentado.

| Arquivo | Função |
|---|---|
| `SimulacaoFDTD.py` | Classe com o núcleo do método FDTD 2D (grade de Yee, atualização de `Hx`/`Hz`/`Ey`, máscara PEC, fonte TE10, ABC de Mur 1ª ordem) |
| `main.py` | Gera as 3 figuras principais (guia / corneta / espaço livre) com escolha física dos instantes de captura |
| `simulacao_guia.py` | Caso isolado: propagação em guia reto sem corneta |
| `simulacao_corneta.py` | Caso isolado: guia + corneta + espaço livre |
| `simulacao_evanescente.py` | Estudo de caso do modo evanescente (`f0 < fc10`) |
| `simulacao_comparacao.py` | Comparação lado a lado: terminação abrupta vs. corneta |

## Estrutura do repositório

```
.
├── SimulacaoFDTD.py
├── main.py
├── simulacao_guia.py
├── simulacao_corneta.py
├── simulacao_evanescente.py
├── simulacao_comparacao.py
├── output/
│   ├── 1_campo_no_guia.png
│   ├── 2_campo_na_corneta.png
│   └── 3_campo_espaco_livre.png
└── resultados/
    ├── simulacao1_prop_guia.png
    ├── simulacao2_prop_guia_corneta.png
    ├── simulacao3_prop_evanescente.png
    └── simulacao4_comparacao.png
```

## Como reproduzir os resultados

```bash
pip install numpy matplotlib imageio

python main.py                     # figuras principais (Seção 6.1)
python simulacao_guia.py           # caso isolado: guia reto
python simulacao_corneta.py        # caso isolado: guia + corneta
python simulacao_evanescente.py    # estudo do modo evanescente (Seção 6.3)
python simulacao_comparacao.py     # comparação guia vs. corneta
```

## Parâmetros físicos utilizados

| Parâmetro | Caso propagante (guia+corneta) | Caso evanescente |
|---|---|---|
| `f0` | 10 GHz | 5 GHz |
| `a` (largura do guia) | 0,02 m | 0,02 m |
| `fc10 = c/(2a)` | 7,5 GHz | 7,5 GHz |
| Condição | `f0 > fc10` → propaga | `f0 < fc10` → evanescente |
| `A` (boca da corneta) | 0,06 m | — |
| Resolução espacial | `λmin/20` (mais fino que o mínimo exigido `λmin/10`) | idem |
| Critério CFL | `Δt = 0,99·Δ/(c√2)` | idem |
| Condição de contorno | PEC nas paredes + ABC de Mur 1ª ordem nas bordas externas | idem |