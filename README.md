# Progetto 1 – Trasmissione e ricezione testo (MATLAB-only)

Sistema Tx–Rx per trasmissione/ricezione di un testo:
testo → bit → (codifica) → BPSK → canale → demod → (decod) → testo

## Funzionalità implementate
- Conversione testo↔bit (senza toolbox)
- Modulazione/Demodulazione **BPSK**
- Misura prestazioni tramite **BER vs SNR** (media su più frame)

## Codifiche di canale
- **NONE** (nessuna codifica)
- **REP3**: codifica a ripetizione N=3 (majority voting in ricezione)
- **HAMMING(7,4)**: FEC con correzione di 1 errore per codeword

## Ambienti simulati
- **AWGN**
- **Rayleigh flat fading + AWGN** (con equalizzazione ideale)
- **Ostacolo**: attenuazione addizionale modellata come riduzione dello SNR (loss = 6 dB)
- **AWGN_SYNC**: errori tipici hardware (timing offset + offset di fase/frequenza), senza compensazione
- **AWGN_SYNC_CORR**: stima e correzione della sincronizzazione tramite **preamble** (timing + CFO + fase)

## Output
- `BER_sync_comparison.png` → confronto: AWGN vs AWGN_SYNC vs AWGN_SYNC_CORR
- `BER_coding_comparison.png` → confronto codifiche: NONE vs REP3 vs HAMMING(7,4)

## Parametri principali
- SNR = 0:2:14 dB
- Nframes = 50
- Preamble = 2×64 simboli (due metà uguali per stima CFO)
- Impairments (AWGN_SYNC): timingDelay = 2, freqOffsetNorm = 0.001, phaseOffset = 10°

## File principali nel repository
- `main_progetto1.m` (codice completo commentato)
- `Sistemi di Comunicazione – Prova finale.pdf` (relazione)
- `results.mat` (risultati salvati)

## Esecuzione
1. Aprire `main_progetto1.m` in MATLAB
2. Premere **Run**
3. Il codice genera automaticamente:
   - `BER_sync_comparison.png`
   - `BER_coding_comparison.png`
   - `results.mat`
   e stampa a console un esempio di testo ricostruito lato ricevitore.
