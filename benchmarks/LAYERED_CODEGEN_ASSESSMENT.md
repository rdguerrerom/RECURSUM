# Assessment: Generating Layered CSE via RECURSUM Codegen

## Executive Summary

**Can we use recursum-expert to implement these changes?**
✅ **YES**, but with caveats. This is a significant architectural addition requiring:
1. A new `LayeredCppGenerator` class alongside the existing `CppGenerator`
2. Fundamental changes to the code generation paradigm
3. Extensive testing to ensure correctness

**Complexity:** HIGH - This is NOT a simple modification like adding `inline` keywords, but a new code generation mode.

---

## Current Architecture Analysis

### 1. Existing Code Generator (`cpp_generator.py`)

**Current paradigm:** Template Metaprogramming (TMP)
- Generates recursive template specializations
- Each `compute()` function returns **ONE** value (e.g., `E^{nA,nB}_t` for specific nA, nB, t)
- Recursion depth determined at compile-time
- SFINAE-based template selection

**Example generated code:**
```cpp
template<int nA, int nB, int N>
struct HermiteCoeff<...> {
    static RECURSUM_FORCEINLINE Vec8d compute(Vec8d PA, Vec8d PB, Vec8d aAB) {
        Vec8d e_0 = HermiteCoeff<nA-1, nB, N-1>::compute(PA, PB, aAB);
        Vec8d e_1 = HermiteCoeff<nA-1, nB, N>::compute(PA, PB, aAB);
        return aAB * e_0 + PA * e_1 + ...;
    }
};
```

### 2. Hand-Written Layered CSE (`hermite_e_layered.hpp`, `coulomb_r_layered.hpp`)

**Layered paradigm:** Layer-by-layer with CSE
- Generates layer structures that compute **ALL** values at once
- Returns `std::array<Vec8d, MAX_SIZE>` containing entire layer
- Uses runtime loops to fill arrays
- Each layer computed exactly once and reused

**Example hand-written code:**
```cpp
template<int nA, int nB>
struct HermiteELayer {
    static std::array<Vec8d, MAX_T> compute(Vec8d PA, Vec8d PB, Vec8d p) {
        std::array<Vec8d, MAX_T> result{};
        auto prev = HermiteELayer<nA-1, nB>::compute(PA, PB, p);  // Computed ONCE

        // Runtime loop filling result array
        for (int t = 1; t < nA + nB; ++t) {
            result[t] = inv2p * prev[t-1] + PA * prev[t] + ...;
        }
        return result;  // Return by value
    }
};
```

---

## Proposed Changes: Analysis & Feasibility

### Change 1: Generate Layered CSE via RECURSUM Codegen ✅ FEASIBLE

**What's needed:**
1. New `LayeredCppGenerator` class in `recursum/codegen/layered_generator.py`
2. Different template structure:
   - Layer templates (compute all values at once)
   - Accessor templates (extract single value from layer)
3. Loop unrolling at codegen time (NOT runtime loops)
4. Array size calculation from index ranges

**Key differences from TMP:**
| Aspect | TMP Generator | Layered Generator |
|--------|--------------|-------------------|
| Return type | `Vec8d` | `std::array<Vec8d, N>` or output param |
| Recursion | Per-value | Per-layer |
| CSE | Via memoization | Via layer reuse |
| Loops | None (compile-time) | Unrolled (codegen time) |

**Leveraging existing infrastructure:**
- ✅ Can reuse `Recurrence`, `RecurrenceRule`, `BaseCase` classes
- ✅ Can reuse `ConstraintSet` for validity checking
- ✅ Can reuse `CodegenContext` for parameter management
- ✅ Can reuse `ExpressionOptimizer` for within-layer CSE
- ⚠️ Need NEW template generation logic

**Implementation approach:**
```python
class LayeredCppGenerator:
    """Generate layer-by-layer C++ code with compile-time CSE."""

    def __init__(self, rec: Recurrence, unroll_loops: bool = True):
        self.rec = rec
        self.unroll_loops = unroll_loops  # Codegen-time unrolling

    def generate(self) -> str:
        """Generate complete C++ header with layer structs."""
        parts = [
            self._header(),
            self._layer_primary_template(),
            self._layer_base_cases(),
            self._layer_rules(),
            self._accessor_template(),  # Single-value accessor
            self._footer()
        ]
        return "\n\n".join(parts)

    def _layer_rules(self):
        """Generate layer computation with unrolled loops."""
        # For each valid (nA, nB):
        #   1. Determine array size (nA + nB + 1 for Hermite)
        #   2. Generate prev layer call (computed ONCE)
        #   3. Unroll loop over t values at codegen time
        #   4. Generate assignments: result[t] = ...
        pass
```

### Change 2: Use Output Parameters Instead of Return-by-Value ✅ RECOMMENDED

**Current bottleneck:**
```cpp
// Return by value: copies entire array
static std::array<Vec8d, MAX_SIZE> compute(...) {
    std::array<Vec8d, MAX_SIZE> result{};
    // ... fill result ...
    return result;  // ❌ Copy on return (though compiler may elide)
}
```

**Proposed improvement:**
```cpp
// Output parameter: no copy, direct write
static RECURSUM_FORCEINLINE void compute(Vec8d PA, Vec8d PB, Vec8d p, Vec8d* out) {
    auto prev = HermiteELayer<nA-1, nB>::compute_to(PA, PB, p, temp_buffer);

    // Write directly to output buffer
    out[0] = PA * prev[0] + prev[1];
    for (int t = 1; t < nA + nB; ++t) {
        out[t] = inv2p * prev[t-1] + PA * prev[t] + ...;
    }
    // No return
}
```

**Benefits:**
- ✅ Eliminates return-by-value overhead
- ✅ Allows stack allocation by caller
- ✅ Better for large arrays (Coulomb R: up to 969 elements!)
- ✅ More cache-friendly (caller controls memory layout)

**Codegen approach:**
```python
def _layer_compute_signature(self, use_output_param: bool = True) -> str:
    if use_output_param:
        sig = f"{self.rec.vec_type}* out"  # Output parameter first
        sig += ", " + ", ".join(f"{self.rec.vec_type} {v}" for v in self.rec.runtime_vars)
        return f"static RECURSUM_FORCEINLINE void compute({sig})"
    else:
        # Traditional return-by-value
        sig = ", ".join(f"{self.rec.vec_type} {v}" for v in self.rec.runtime_vars)
        return f"static RECURSUM_FORCEINLINE std::array<{self.rec.vec_type}, MAX_SIZE> compute({sig})"
```

### Change 3: Apply `RECURSUM_FORCEINLINE` to Layer Functions ✅ TRIVIAL

**Already implemented transversally!**
- ✅ `cpp_generator.py` already adds `RECURSUM_FORCEINLINE` macro
- ✅ Applied to ALL `compute()` methods in TMP generator
- ✅ Need to apply same pattern to LayeredCppGenerator

**Code pattern:**
```python
# In _layer_template():
return f"""template<{tparams}>
struct {self.ctx.struct_name}Layer<{targs}> {{
    static RECURSUM_FORCEINLINE void compute({sig}) {{
        {body}
    }}
}};"""
```

---

## Performance Impact Analysis

### Current Layered CSE Bottlenecks (from benchmarks):

| Issue | Current Impact | Solution |
|-------|----------------|----------|
| Runtime loops | ❌ Branch prediction overhead | ✅ Unroll at codegen time |
| Return-by-value | ❌ 969-element array copies | ✅ Output parameters |
| MAX-sized arrays | ❌ Cache pollution | ✅ Exact-sized arrays |
| Missing inline | ❌ Function call overhead | ✅ RECURSUM_FORCEINLINE |

**Expected improvement after codegen:**
- Hermite E: 2× - 5× faster (eliminate 6× overhead at L=0)
- Coulomb R: 10× - 100× faster (eliminate 1,290× overhead at L=0)

**Why such large improvements expected?**
1. **Compile-time unrolling**: Entire layer computation becomes straight-line code
2. **No function overhead**: FORCEINLINE + compile-time unrolling
3. **Exact-sized arrays**: No wasted memory, better cache utilization
4. **Output parameters**: Zero-copy data movement

---

## Implementation Plan

### Phase 1: Architecture ✅ Can use recursum-expert
**Files to modify:**
- `recursum/codegen/layered_generator.py` (NEW)
- `recursum/codegen/__init__.py` (add LayeredCppGenerator export)
- `recursum/codegen/recurrence.py` (add `.generate_layered()` method)

**Tasks:**
1. Create `LayeredCppGenerator` class skeleton
2. Implement layer structure generation
3. Implement compile-time loop unrolling
4. Implement output parameter mode
5. Add RECURSUM_FORCEINLINE macro (reuse existing)

### Phase 2: Testing ✅ Can automate
**Files to create:**
- `recursum/codegen/test_layered_gen.py` (unit tests)
- Regenerate Hermite E using layered generator
- Regenerate Coulomb R using layered generator

**Validation:**
1. Compare layered codegen output vs hand-written
2. Run existing benchmarks with generated code
3. Verify correctness against SciPy/symbolic solutions

### Phase 3: Integration ✅ Can use recursum-expert
**Files to modify:**
- Update example scripts to use `.generate_layered()`
- Update documentation in docstrings
- Update benchmarks to use generated layered code

---

## Risk Assessment

### HIGH RISK ⚠️
1. **Complexity:** Layered generation is fundamentally different from TMP
   - Mitigation: Extensive testing, gradual rollout
2. **Loop unrolling explosion:** Code size could be massive for high L
   - Mitigation: Make unrolling optional, add size warnings

### MEDIUM RISK ⚠️
1. **Array size calculation:** Getting exact sizes wrong breaks everything
   - Mitigation: Unit tests for all index ranges
2. **Output parameter API:** Caller must provide correct buffer size
   - Mitigation: Add runtime assertions in debug mode

### LOW RISK ✅
1. **Inline application:** Already proven to work in TMP generator
2. **CSE within layers:** Can reuse existing optimizer infrastructure

---

## Recommendation

### ✅ YES - Use recursum-expert for this implementation

**Rationale:**
1. This is a well-defined architectural addition with clear requirements
2. The recursum-expert agent has demonstrated capability with Python codegen
3. Existing infrastructure (Recurrence DSL, optimizer) can be reused
4. Expected performance gains are substantial (10-100× for Coulomb R)

**Suggested approach:**
1. **Start small:** Implement LayeredCppGenerator for Hermite E first
2. **Validate:** Compare generated code vs hand-written implementation
3. **Benchmark:** Verify performance improvements
4. **Extend:** Apply to Coulomb R and other recurrences
5. **Generalize:** Make layered generation available for ALL recurrences

**Time estimate:**
- Phase 1 (Architecture): 4-6 hours of focused work
- Phase 2 (Testing): 2-3 hours
- Phase 3 (Integration): 1-2 hours
- **Total:** ~8-11 hours of development + testing

**Confidence level:** 🔥 **HIGH** - This is feasible and will provide significant performance improvements.

---

## Example Generated Code (Target)

### Hermite E Layered (Generated by RECURSUM):
```cpp
template<int nA, int nB>
struct HermiteELayer {
    static constexpr int N = nA + nB + 1;

    static RECURSUM_FORCEINLINE void compute(Vec8d* out, Vec8d PA, Vec8d PB, Vec8d p) {
        // Stack-allocated buffer for previous layer
        Vec8d prev[nA + nB];
        HermiteELayer<nA-1, nB>::compute(prev, PA, PB, p);

        Vec8d inv2p = Vec8d(0.5) / p;

        // UNROLLED at codegen time (no runtime loop!)
        out[0] = PA * prev[0] + prev[1];
        out[1] = inv2p * prev[0] + PA * prev[1] + Vec8d(2) * prev[2];
        out[2] = inv2p * prev[1] + PA * prev[2] + Vec8d(3) * prev[3];
        // ... (continue for all t values)
    }
};
```

**Key improvements over hand-written:**
- ✅ Compile-time unrolled (no runtime loop)
- ✅ Output parameter (no return-by-value copy)
- ✅ Exact-sized prev buffer (no MAX_SIZE waste)
- ✅ RECURSUM_FORCEINLINE applied

---

## Conclusion

**Final verdict:** ✅ **PROCEED with recursum-expert implementation**

This is a high-value, high-complexity feature that will:
1. Eliminate performance bottlenecks in layered CSE (1,290× overhead → near-zero)
2. Make layered generation transversal to ALL RECURSUM recurrences
3. Provide users with choice: TMP (single-value) vs Layered (all-values)
4. Maintain code generation at Python level (not hand-written C++)

The recursum-expert agent is well-suited for this task given its previous success with the TMP generator modifications.
