#pragma once
// Per-quartet HGP-OS recurrence scalars (the prerequisites the generated
// kernel consumes). Field names match osx emit.py:240-249. The caller fills
// these once per primitive quartet (outside the timed kernel), mirroring how
// libint2's low-level API receives its Libint_t prerequisites.
struct ScalarPack {
    double PAx, PAy, PAz;
    double QCx, QCy, QCz;
    double WPx, WPy, WPz;
    double WQx, WQy, WQz;
    double ABx, ABy, ABz;
    double CDx, CDy, CDz;
    double inv_2zp, inv_2zq, inv_2zpq;
    double frac_q_over_pq, frac_p_over_pq;
};
