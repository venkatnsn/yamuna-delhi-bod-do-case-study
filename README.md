# Yamuna at Delhi — a water-quality modelling case study

This is a short technical case study on the BOD–DO problem in the 22 km Wazirabad → Okhla reach of the Yamuna in Delhi. I put it together as a portfolio piece for water-quality modelling roles, and as a counterpart to my earlier Elk Valley selenium case study: a different region (India under CPCB), a different pollutant family (organic load + DO sag rather than metalloid speciation), and a different physical core (coupled Streeter–Phelps ODEs along a multi-source 1D reach), but the same end-to-end workflow — calibrate, quantify uncertainty, translate to compliance probability, frame the trade-off.

Everything here is built from public data (CPCB Yamuna Action Plan reports, CWC Old Railway Bridge gauge records, CPCB station medians 2018–2022). The case is a teaching example; it is not a site-specific regulatory model.

**Author.** Dr Venkat Nara — Research Fellow, Water Quality Modelling, Geochemistry & Hydrology, University of Melbourne.
**Contact.** [venkat.nsn@gmail.com](mailto:venkat.nsn@gmail.com) · +61 451 589 633 · [LinkedIn](https://linkedin.com/in/sivanagavenkatnara)

---

## Why Yamuna at Delhi

The Wazirabad → Okhla reach is a highly challenging urban-river BOD–DO case under dry-season conditions, and one of the most-studied reaches of organic-load surface water in South Asia. It sits at the intersection of several things that a water-quality modeller actually has to deal with: an order-of-magnitude imbalance between dry-season headwater flow (~4 m³/s after Wazirabad diversions) and combined municipal drain inflow (~52 m³/s through Najafgarh, Shahdara, and Sahibabad), a CPCB Class C compliance framework with two simultaneous limits (DO ≥ 4 mg/L *and* BOD ≤ 3 mg/L), a regulator that publishes monthly station-level data, and a planning horizon dominated by capital-allocation decisions on STP capacity and managed environmental flow. That combination is exactly what I wanted to demonstrate I can model end-to-end.

The framing throughout is consulting-grade rather than purely academic: a screening-level river water-quality case study that combines calibration, uncertainty quantification, and intervention prioritisation under a CPCB compliance framework — the kind of model used to rapidly test whether investment should prioritise treatment, flow release, interception, or a combined pathway before moving to higher-complexity tools (QUAL2K, WASP).

The same workflow — coupled-ODE process model, calibration against routine monitoring, LHS Monte Carlo, scenario-as-compliance-probability — transfers straightforwardly to BOD/DO modelling on the Hindon, Sabarmati, Musi, or any Indian municipal-impacted reach under a CPCB Class designation; or to ammonia and fecal coliform compliance under the same framework with a swap of the kinetics block.

---

## What's in the repository

| File | What it is |
|---|---|
| `Yamuna_Delhi_Case_Study.pdf` | The primary deliverable. Five-page technical report: executive summary, management-insight box, context, study area, methodology (process model, calibration, uncertainty, compliance metric), key findings, scenario table, limitations, transferability, references. This is what I attach to applications. |
| `Yamuna_Delhi_Methods_Appendix.pdf` | Two-page auditable appendix: 11 governing equations, all five Monte Carlo input distributions, station-by-station calibration statistics, scenario summary with marginal and joint compliance probabilities, sensitivity ranking. |
| `notebook/Yamuna_BODDO_Analysis.ipynb` | The Python notebook behind the report (process model, drain-mixing, calibration grid, 2000-run LHS Monte Carlo, scenarios). The notebook also contains an *illustrative* NH₃ screening (Step 13) — labelled in the notebook as not a model output. End-to-end runtime ≈ 60 s. |
| `notebook/Yamuna_Fig1…Fig7.png` | Source-term, longitudinal profiles, calibration scatter, parameter SSR surface, MC scatter + Spearman ranking, scenario compliance, illustrative NH₃ screening (notebook only — not used in the case study). |
| `build/` | Builders that regenerate the notebook (`make_yamuna_nb.py`) and the two DOCX files (`make_case_study.js`, `make_methods_appendix.js`). |

---

## Headline results

Every number here appears in the report. I'm listing them again so anyone reading the repo first doesn't have to open the PDF.

| Metric | Value |
|---|---|
| Monte Carlo runs | 2000 (Latin Hypercube) |
| Uncertain parameters | 5 (Q_head, k_d,20, k_a multiplier, drain BOD, drain Q) |
| Calibration stations | 4 (Palla, Nizamuddin, ITO bridge, Okhla) |
| Calibration data | CPCB monthly medians 2018–2022 (pre- and post-monsoon used; monsoon held back) |
| Pooled NSE_BOD | 0.78 ("Very Good" under Moriasi et al. 2007) |
| Pooled NSE_DO  | 0.95 |
| Pooled PBIAS_BOD | +1.6% |
| Pooled PBIAS_DO  | −0.9% |
| Best k_d,20 | 0.46 1/d |
| Best k_a multiplier (× O'Connor–Dobbins) | 1.10 |
| **S0** baseline — P(DO ≥ 4) | 4.8% |
| **S1** drain interception + STP retrofit (η ≈ 50%) | 46.4% |
| **S2** aggressive STP upgrade (η ≈ 70%) | 81.5% ← clears my 75% planning threshold on drain treatment alone |
| **S3** e-flow release at Wazirabad (×3) | 7.8% |
| **S4** S2 + S3 combined | 86.9% |
| P(BOD ≤ 3) under S4 | 0% — BOD compliance is effectively unreachable under the tested scenarios |
| Dominant Spearman driver of Okhla BOD | drain BOD multiplier, ρ = +0.93 |
| Dominant Spearman driver of Okhla DO | drain BOD ρ = −0.56; k_a ρ = +0.53; k_d ρ = −0.54 |
| Headwater Q sensitivity | |ρ| < 0.1 — essentially decoupled at dry-season Q_head/Q_drain ratio |

---

## Methods, briefly

**Governing equations** (full list with units in the Methods Appendix):

- BOD: ∂L/∂x = −(k_d / u_km)·L
- DO: ∂O/∂x = (k_a / u_km)·(O_s − O) − (k_d / u_km)·L
- Temperature: k_d(T) = k_{d,20}·1.047^(T−20); k_a(T) = k_{a,20}·1.024^(T−20)
- Reaeration: k_{a,20} = 3.93·u^0.5 / H^1.5 (O'Connor–Dobbins)
- DO saturation: O_s(T) = 14.652 − 0.41022·T + 0.00799·T² − 7.78×10⁻⁵·T³ (Benson–Krause, sea level)
- Tributary mixing: C_mix = (Q_u·C_u + Q_t·C_t) / (Q_u + Q_t)
- Treatment: C_out = C_in·(1 − η)
- Compliance: P(comply) = (1/N)·Σ 𝟙{ DO_i ≥ 4 ∩ BOD_i ≤ 3 }, N = 2000
- NH₃ receptor screening (illustrative only, notebook Step 13 — not part of the calibrated model): f_NH3 = 1 / [1 + 10^(pKa − pH)], pKa = 0.0902 + 2729.92/(T+273.15) (Emerson et al. 1975)

**Uncertainty.** Lognormal on headwater Q multiplier (σ=0.35); Uniform on k_{d,20} (0.20–0.55 1/d, widened from the Chapra 1997 urban-stream range to bracket the calibrated optimum at ~0.46 1/d) and k_a multiplier (0.7–1.4); Lognormal on drain BOD multiplier (σ=0.25); Triangular on drain Q multiplier (0.6–1.0–1.3). I deliberately keep the priors bounded rather than over-specifying distributions the available data don't support, and DO is floored at zero in the MC outputs because negative DO is unphysical. Spearman rank correlations identify drain BOD as the dominant uncertainty driver (ρ ≈ +0.93 vs Okhla BOD): source-strength uncertainty dominates over kinetic uncertainty in the dry-season regime.

**Scenarios.**

| Tag | Description | Median BOD (mg/L) | Median DO (mg/L) | P(DO ≥ 4) | P(BOD ≤ 3) | P(joint) |
|---|---|---|---|---|---|---|
| S0 | Baseline (current operations) | 27.4 | 1.82 | 4.8% | 0% | 0% |
| S1 | Drain interception + STP retrofit (η ≈ 50%) | 13.8 | 3.91 | 46.4% | 0% | 0% |
| S2 | Aggressive STP upgrade (η ≈ 70%) | 8.3 | 4.77 | 81.5% | 0% | 0% |
| S3 | E-flow release at Wazirabad (×3) | 24.3 | 2.29 | 7.8% | 0% | 0% |
| S4 | S2 + S3 combined | 7.7 | 4.91 | 86.9% | 0% | 0% |

**The binding constraint is BOD, not DO.** Joint Class C is effectively unreachable on this reach under the tested scenarios in the pre-monsoon: even S4 leaves median Okhla BOD at ~7.7 mg/L, and the BOD-marginal compliance probability is essentially zero across all five scenarios. Reaching BOD ≤ 3 mg/L would require either >90% drain treatment efficiency or a dry-season e-flow that doesn't exist as a credible policy lever. **DO ≥ 4 mg/L *is* achievable** — under aggressive treatment alone (S2), P(DO ≥ 4) at Okhla rises from 4.8% to 81.5%, and S4 (S2 + ×3 e-flow) lifts it to 86.9%. Aggressive drain treatment is the high-leverage dry-season lever; e-flow augmentation alone (S3) moves DO compliance only to ~8%.

---

## Scope and limitations

This is a screening-level, public-data case study built on a steady-state 1D representation. It does not resolve diurnal DO swings, sediment oxygen demand, nitrification kinetics, or the slack-water reservoirs above each barrage where eutrophic algal cycles dominate the DO budget. Drain BOD and Q are treated as time-invariant within each season; in reality both have a diurnal and a stormwater signature. Headwater Q is the most contestable input — CWC and CPCB report different dry-season Wazirabad release figures — but not a high-leverage input on this reach: the Spearman ranking shows |ρ| < 0.1 for headwater Q against Okhla BOD, because dry-season Q_head is an order of magnitude smaller than the combined drain inflow. I have erred on the conservative side wherever the data are ambiguous. The model is intended for dry-season planning and should not be used to infer monsoon BOD compliance quantitatively without a storage-aware dynamic model.

---

## Reproducibility

**Requirements** (Python 3.10+): `numpy`, `pandas`, `scipy`, `matplotlib`, `jupyter`, `nbformat`.

```bash
git clone <this-repo>
cd yamuna-delhi-bod-do-case-study
pip install numpy pandas scipy matplotlib jupyter nbformat
jupyter notebook notebook/Yamuna_BODDO_Analysis.ipynb
```

Then `Kernel → Restart & Run All`. End-to-end runtime on a laptop is about 60 seconds; the slowest steps are the 25 × 25 grid sweep (~5 s) and each 2000-run LHS Monte Carlo block (~10 s).

To regenerate the case study and methods appendix:

```bash
# Notebook (rebuild structure, then execute):
python build/make_yamuna_nb.py
python -m jupyter nbconvert --to notebook --execute --inplace notebook/Yamuna_BODDO_Analysis.ipynb


```

---

## Key references

1. Chapra, S.C. (1997). *Surface Water-Quality Modeling.* McGraw-Hill, New York.
2. CPCB (2021). *Yamuna Water Quality Status Report 2018–2020.* Central Pollution Control Board, New Delhi.
3. Emerson, K., Russo, R.C., Lund, R.E., Thurston, R.V. (1975). Aqueous ammonia equilibrium calculations: effect of pH and temperature. *J. Fish. Res. Bd. Canada* 32: 2379–2383.
4. McKay, M.D., Beckman, R.J., Conover, W.J. (1979). A comparison of three methods for selecting values of input variables in the analysis of output from a computer code. *Technometrics* 21(2): 239–245.
5. Moriasi, D.N., Arnold, J.G., Van Liew, M.W., et al. (2007). Model evaluation guidelines for systematic quantification of accuracy in watershed simulations. *Trans. ASABE* 50(3): 885–900.
6. O'Connor, D.J., Dobbins, W.E. (1958). Mechanism of reaeration in natural streams. *Trans. ASCE* 123: 641–684.
7. Streeter, H.W., Phelps, E.B. (1925). *A Study of the Pollution and Natural Purification of the Ohio River.* US Public Health Service Bulletin 146.
8. USEPA (2013). *Aquatic Life Ambient Water Quality Criteria for Ammonia — Freshwater 2013.* EPA-822-R-13-001.

---

## License

The written case study and methods appendix are released under CC BY 4.0 — feel free to share and adapt with attribution. The Python / Jupyter code is released under the MIT License (see `LICENSE`).

---

## Contact

For questions, collaboration, or regenerated model outputs / figures:

**Dr Venkat Nara**
Research Fellow — Water Quality Modelling, Geochemistry & Hydrology
University of Melbourne, Australia
[venkat.nsn@gmail.com](mailto:venkat.nsn@gmail.com)  |  +61 451 589 633
[LinkedIn](https://linkedin.com/in/sivanagavenkatnara)
