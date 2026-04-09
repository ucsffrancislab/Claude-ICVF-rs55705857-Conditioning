# Analysis Framework: Testing Whether ICVF Polygenic Architecture Emerges for IDHmt Glioma After Conditioning on rs55705857

## Background & Hypothesis

In our PGS catalog-wide analysis, 18 ICVF PGS were significantly protective for IDHwt glioma (OR 0.85–0.93) but **none reached significance for IDHmt**. We hypothesize that the ICVF polygenic signal is present for IDHmt but is masked by the overwhelming effect of rs55705857 (CCDC26/MYC, OR ~6.5), which acts as a monogenic switch for OPC self-renewal — the same biological axis that ICVF PGS capture polygenically for IDHwt.

**Prediction:** If we remove or condition on rs55705857, an ICVF polygenic signal should emerge for IDHmt glioma risk among non-carriers or within genotype strata.

---

## Analysis Design

### Required Data

- Individual-level genotype data for IDHmt cases and controls
- rs55705857 genotype (chr8:128,748,020 on GRCh38; A>G; risk allele = G)
- Pre-computed ICVF PGS for each individual (the 18 FDR-significant scores from our primary analysis)
- IDHmt case/control status (and if possible, IDHmt-intact vs. IDHmt-codel)
- Standard covariates: age, sex, principal components (PCs 1–10)

### Analysis 1: Stratified Association Testing

**Approach:** Split the sample by rs55705857 genotype and test ICVF PGS–IDHmt associations within each stratum.

```
Strata:
  - AA homozygotes (non-carriers of risk allele) — expected to be ~85–90% of sample
  - AG heterozygotes + GG homozygotes (carriers) — expected to be ~10–15%

Model within each stratum:
  logit(IDHmt case) ~ ICVF_PGS + age + sex + PC1–PC10

Test each of the 18 ICVF PGS separately within the AA stratum.
Apply FDR correction across 18 tests.
```

**Expected result if hypothesis is correct:** ICVF PGS associations become nominally or FDR-significant in the AA (non-carrier) stratum, with protective direction (OR < 1) consistent with IDHwt results.

**Expected result if hypothesis is wrong:** No signal in either stratum — the absence of ICVF signal for IDHmt is real, not masked.

### Analysis 2: Conditional Association Testing (Full Sample)

**Approach:** Include rs55705857 genotype as a covariate and test whether ICVF PGS are associated with IDHmt risk after adjustment.

```
Model:
  logit(IDHmt case) ~ ICVF_PGS + rs55705857_dosage + age + sex + PC1–PC10

Test each of the 18 ICVF PGS.
Apply FDR correction across 18 tests.
Compare effect sizes (OR) and p-values to:
  (a) the unconditional IDHmt results (our current analysis)
  (b) the IDHwt results (our current analysis)
```

**Key comparison:** If conditioning on rs55705857 makes the ICVF effect estimates for IDHmt converge toward the IDHwt estimates (in magnitude and significance), this supports the hypothesis that the two architectures are functional analogs.

### Analysis 3: Interaction Test

**Approach:** Formally test whether the ICVF PGS effect on IDHmt risk differs by rs55705857 carrier status.

```
Model:
  logit(IDHmt case) ~ ICVF_PGS + rs55705857_dosage + ICVF_PGS × rs55705857_dosage + age + sex + PC1–PC10

Test the interaction term (ICVF_PGS × rs55705857_dosage).
```

**Interpretation:**
- Significant negative interaction: rs55705857 carriers are less sensitive to ICVF polygenic effects (consistent with the monogenic switch dominating)
- Non-significant interaction but main effect of ICVF becomes significant: the signals are additive/independent (still supports the hypothesis, but with a different architecture)

### Analysis 4: ICVF PGS Excluding rs55705857

**Approach:** Reconstruct ICVF PGS after removing rs55705857 (and any variants in LD, r² > 0.1 within ±1Mb) to ensure the signal is not driven by the shared variant itself.

```
For each of the 18 ICVF PGS:
  1. Identify rs55705857 and LD proxies (r² > 0.1, ±1Mb window) in the PGS variant list
  2. Remove these variants and recompute the PGS
  3. Re-test: logit(IDHmt case) ~ ICVF_PGS_pruned + age + sex + PC1–PC10
```

**Purpose:** This distinguishes between two scenarios:
- (a) The ICVF–IDHmt signal (if it emerges) is driven entirely by the shared CCDC26/MYC variant → signal disappears after pruning
- (b) The ICVF–IDHmt signal reflects a genuine distributed polygenic architecture beyond the 8q24 locus → signal persists after pruning

**Scenario (b) is the stronger evidence** for our hypothesis. It would mean the OPC differentiation program captured by ICVF PGS influences IDHmt risk through many variants, not just the one shared locus.

---

## Power Considerations

- rs55705857 risk allele frequency is ~5–7% in European populations
- The non-carrier stratum (AA) will retain ~85–90% of the sample → minimal power loss for Analysis 1
- The key limiting factor is total IDHmt case count — estimate power for detecting ORs in the 0.85–0.95 range (the effect sizes we see for IDHwt) at the available sample size
- If individual-level data are limited, Analysis 2 (conditional on rs55705857 dosage) is preferred over Analysis 1 (stratified) because it uses the full sample
- Consider meta-analyzing across available IDHmt cohorts if single-cohort power is insufficient

## Interpretation Guide

| Result Pattern | Interpretation |
|---------------|---------------|
| ICVF signal emerges in AA stratum (Analysis 1) AND after conditioning (Analysis 2) AND survives LD pruning (Analysis 4) | **Strong support** — ICVF polygenic architecture is a true functional analog of rs55705857 for IDHmt, operating through distributed OPC differentiation variants |
| ICVF signal emerges after conditioning but disappears after LD pruning | **Partial support** — the signal is real but driven by the shared 8q24 locus, not a distributed architecture |
| ICVF signal emerges in AA stratum but interaction is non-significant (Analysis 3) | **Additive model** — ICVF and rs55705857 contribute independently to IDHmt risk (still supports the differentiation–transformation axis, but they're not redundant) |
| No ICVF signal in any analysis | **Hypothesis not supported** — the absence of ICVF signal for IDHmt reflects genuine biological differences between IDHwt and IDHmt risk architecture, not masking by rs55705857 |

## Priority Order

If time/resources are limited, run in this order:
1. **Analysis 2** (conditional, full sample) — most efficient, answers the primary question
2. **Analysis 4** (LD-pruned PGS) — critical for distinguishing shared-locus vs. distributed architecture
3. **Analysis 1** (stratified) — provides intuitive visualization of effect within non-carriers
4. **Analysis 3** (interaction) — formal test but requires largest sample size for adequate power
