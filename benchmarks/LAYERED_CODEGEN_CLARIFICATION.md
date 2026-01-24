# LayeredCodegen Representation in Benchmark Plots - Clarification

**Date**: 2026-01-16
**Status**: Analysis Complete - Action Plan Ready

---

## Executive Summary

After thorough investigation of all benchmark plots and data, I've identified **how LayeredCodegen appears across different benchmarks**. The confusion arises because **not all benchmarks include LayeredCodegen** - some benchmarks only compare earlier implementations (TMP, Layered, Symbolic).

### Key Finding: LayeredCodegen is ONLY in Hermite Coefficients Benchmarks

**LayeredCodegen (impl=3) is present ONLY in:**
- Hermite expansion coefficients benchmarks (`hermite_coefficients.json`)
- Appears in 3 plots explicitly labeled as "LayeredCodegen"

**LayeredCodegen is NOT present in:**
- Coulomb Hermite integrals benchmarks (only TMP vs Layered)
- J/K matrix benchmarks (uses LayeredCodegen internally but doesn't benchmark it separately)

---

## Complete Implementation Index Mapping

Based on analysis of all benchmark files and plotting scripts:

### Hermite Coefficients Benchmarks (hermite_coefficients.json)
- **impl = 0**: TMP (Original template metaprogramming)
- **impl = 1**: Layered (Hand-written layered implementation)
- **impl = 2**: Symbolic (SymPy-generated polynomial expressions)
- **impl = 3**: LayeredCodegen (NEW - automatically generated layered code)

### Coulomb Hermite Benchmarks (coulomb_hermite.json)
- **impl = 0**: TMP
- **impl = 1**: Layered
- **No impl=3**: LayeredCodegen NOT benchmarked for Coulomb integrals

### J/K Matrix Benchmarks (jk_alkanes_results.csv)
- Only measures end-to-end J/K matrix construction time
- Uses LayeredCodegen internally for Hermite coefficients
- Does NOT have separate implementation comparisons
- Shows "baseline vs RECURSUM" (where RECURSUM = LayeredCodegen approach)

---

## Plot-by-Plot Breakdown

### Plots WITH LayeredCodegen Label:

#### 1. hermite_coefficients_comparison.png/pdf
- **Location**: `benchmarks/results/figures/`
- **Shows**: 4-way comparison including LayeredCodegen
- **Labels**: "TMP (Original)", "Layered (Hand-written)", "Symbolic (SymPy)", "LayeredCodegen (NEW)"
- **Status**: ✅ Correct - LayeredCodegen properly labeled

#### 2. hermite_coefficients_vs_L.png/pdf
- **Location**: `benchmarks/results/figures/`
- **Shows**: Performance scaling with angular momentum L
- **Labels**: "TMP", "Layered", "Symbolic", "LayeredCodegen"
- **Status**: ✅ Correct - LayeredCodegen properly labeled

#### 3. hermite_layered_codegen_speedup.png/pdf
- **Location**: `benchmarks/results/figures/`
- **Shows**: LayeredCodegen speedup bar chart
- **Labels**: "vs Layered", "vs TMP", "vs Symbolic"
- **Status**: ✅ Correct - explicitly about LayeredCodegen

### Plots WITHOUT LayeredCodegen (Correctly):

#### 4. coulomb_hermite_comparison.png/pdf
- **Location**: `benchmarks/results/figures/`
- **Shows**: TMP vs Layered for Coulomb R integrals
- **Labels**: "TMP", "Layered"
- **Why no LayeredCodegen**: Coulomb integrals benchmark doesn't include impl=3
- **Status**: ✅ Correct - LayeredCodegen not benchmarked for this recurrence type

#### 5. coulomb_hermite_scaling.png/pdf
- **Location**: `benchmarks/results/figures/`
- **Shows**: Scaling analysis for Coulomb integrals
- **Labels**: "TMP", "Layered"
- **Why no LayeredCodegen**: Same as above
- **Status**: ✅ Correct

### Plots with Indirect LayeredCodegen Reference:

#### 6-9. J/K Matrix Plots (jk_*.pdf)
- **Location**: `benchmarks/results/figures/`
- **Plots**: jk_comparison, jk_scaling_analysis, jk_recursum_impact, jk_combined_overview
- **Labels**: "Hand-Written (Baseline)" vs "RECURSUM LayeredCodegen" or just "LayeredCodegen"
- **Meaning**: These compare:
  - **Baseline**: What performance would be WITHOUT LayeredCodegen (9.8× slower Hermite coefficients)
  - **RECURSUM**: What performance IS with LayeredCodegen (fast Hermite coefficients)
- **Status**: ✅ Correct - shows downstream impact of LayeredCodegen on application performance

---

## Manuscript Consistency Check

### Manuscript Claims:
1. "LayeredCodegen achieves 9.8× speedup over hand-written Layered" ✅ Supported by Plot 3
2. "LayeredCodegen achieves 1.9× speedup over TMP" ✅ Supported by Table and Plots 1-2
3. "LayeredCodegen systematically outperforms all implementations" ✅ Supported by Plots 1-3
4. J/K matrix performance enabled by LayeredCodegen ✅ Supported by Plots 6-9

### Plot References in Manuscript:
- Section 5.2 (LayeredCodegen benchmarks): References Figures 1-3 (hermite plots) ✅
- Section 5.3 (Coulomb benchmarks): References Figures 4-5 (coulomb plots, no LayeredCodegen) ✅
- Section 5.4 (J/K matrices): References Figure 6+ (jk plots show LayeredCodegen impact) ✅

**Conclusion**: Manuscript is consistent with plot labeling.

---

## Why Coulomb Benchmarks Don't Include LayeredCodegen

After examining the code structure:

1. **LayeredCodegen Implementation Status**:
   - ✅ Fully implemented for Hermite E coefficients (3-index linear recurrence)
   - ⚠️ Under development for Coulomb R integrals (4-index tetrahedral recurrence)

2. **Manuscript Acknowledgment**:
   - Program Summary section states: "LayeredCodegen backend: currently supports 3-index Hermite coefficients; 4-index tetrahedral Coulomb integrals under development"
   - This explains why Coulomb plots show only TMP vs Layered

3. **No Correction Needed**: The absence of LayeredCodegen from Coulomb plots is intentional and documented.

---

## Validation: Cross-Check with Data Files

### hermite_coefficients.json:
```bash
$ grep '"impl":' hermite_coefficients.json | sort -u
"impl": 0.0  # TMP
"impl": 1.0  # Layered
"impl": 2.0  # Symbolic
"impl": 3.0  # LayeredCodegen ✅
```

### coulomb_hermite.json:
```bash
$ grep '"impl":' coulomb_hermite.json | sort -u
"impl": 0.0  # TMP
"impl": 1.0  # Layered
# No impl=3 ✅ Confirmed
```

---

## Assessment: Are Plots Correct?

### ✅ All Plots Are Correctly Labeled

After comprehensive analysis:

1. **Hermite plots (1-3)**: Correctly show all 4 implementations including LayeredCodegen
2. **Coulomb plots (4-5)**: Correctly show only TMP vs Layered (LayeredCodegen not implemented for Coulomb)
3. **J/K plots (6-9)**: Correctly show LayeredCodegen's downstream performance impact
4. **Manuscript text**: Correctly describes what each plot shows
5. **Figure captions**: Correctly specify which implementations are compared

### No Corrections Needed

The confusion likely arose from:
- Not all benchmarks include LayeredCodegen (by design)
- J/K plots use "RECURSUM LayeredCodegen" label (correct but different wording)
- Coulomb benchmarks intentionally exclude LayeredCodegen (not yet implemented)

---

## Recommendation: Documentation Enhancement

While plots are correct, we can improve clarity:

### Option A: Add Clarifying Text to Manuscript (Minimal Change)

Add one sentence to Section 5 introduction:

> "Note: LayeredCodegen benchmarks are presented for Hermite coefficients only (Section 5.2, Figures 1-3), as the 4-index Coulomb integral backend is under development (Section 5.3 uses TMP and hand-written Layered implementations). Section 5.4 demonstrates LayeredCodegen's downstream impact on J/K matrix construction."

### Option B: Add Figure Caption Notes (More Explicit)

Update figure captions:

- **Figure 4 (Coulomb comparison)**: Add note "(LayeredCodegen backend for Coulomb integrals under development)"
- **Figure 6+ (J/K matrices)**: Add note "(uses LayeredCodegen for Hermite coefficient computation)"

### Option C: Create Summary Table (Comprehensive)

Add table showing which implementations are benchmarked for each recurrence type:

| Recurrence Type | TMP | Layered | Symbolic | LayeredCodegen |
|-----------------|-----|---------|----------|----------------|
| Hermite E coeffs | ✓ | ✓ | ✓ | ✓ |
| Coulomb R integrals | ✓ | ✓ | ✗ | ✗ |
| J/K matrices (application) | - | - | - | ✓ (internal) |

---

## Action Items

### Required: None
All plots are correct and consistent with manuscript.

### Optional (Documentation Enhancement):
- [ ] Add clarifying sentence to Section 5 introduction (Option A)
- [ ] Update figure captions for Coulomb and J/K plots (Option B)
- [ ] Add implementation coverage table (Option C)

### If You Want to Regenerate Plots:
All plotting scripts are available and documented:
- `benchmarks/results/figures/regenerate_all_plots.py` - Main 5 plots
- `benchmarks/scripts/plot_jk_alkanes.py` - J/K matrix plots
- Uses data from `benchmarks/results/raw/*.json`

**No regeneration needed** - current plots are correct.

---

## Summary for Manuscript Integration

**Key Points to Remember:**

1. ✅ **Hermite benchmarks show LayeredCodegen** (impl=3) - 3 plots
2. ✅ **Coulomb benchmarks don't show LayeredCodegen** - by design (not implemented yet)
3. ✅ **J/K benchmarks show LayeredCodegen impact** - downstream performance benefit
4. ✅ **All manuscript claims are supported** by correctly labeled plots
5. ✅ **No plot corrections needed** - labeling is accurate

**Bottom Line**: The "confusion" stems from not all benchmarks including LayeredCodegen. This is correct and intentional - LayeredCodegen is only benchmarked where it has been implemented (Hermite coefficients). The manuscript properly documents this in the Program Summary section.

---

**Analysis completed**: 2026-01-16
**Conclusion**: All plots correctly represent benchmark data. Optional documentation enhancements suggested above.
