# Progetto 1 – Trasmissione e ricezione testo (MATLAB-only)

Sistema Tx–Rx per trasmissione/ricezione di un testo:
testo→bit → (codifica a ripetizione N=3) → BPSK → canale → demod → (decod) → testo

## Ambienti simulati
- **AWGN**
- **Rayleigh flat fading + AWGN** (con equalizzazione ideale)
- **Ostacolo**: attenuazione addizionale modellata come riduzione dello SNR (loss = 6 dB)
- **AWGN_SYNC**: AWGN con imperfezioni di sincronizzazione (offset di fase/frequenza, senza compensazione)

## Output
- Curva **BER vs SNR**: `BER_vs_SNR.png`

## Parametri principali
- SNR = 0:2:14 dB
- Nframes = 50
- Codifica: ripetizione N=3
