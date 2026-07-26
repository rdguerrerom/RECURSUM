"""Turn results_eri.json into a RECURSUM-vs-libint2 comparison table + liveness
stats. Reports ns/quartet, ns/integral, and speedup per Cartesian class."""
import json
import sys

from recursum_dag_emit import emit_kernel, class_name, CANON_LADDER

CLASSES = CANON_LADDER


def load(path):
    d = json.load(open(path))
    t = {}  # (impl, cls) -> {mean, stddev, nout}
    for b in d["benchmarks"]:
        name = b["name"]
        if not name.endswith("_mean"):
            continue
        parts = name.split("/")            # e.g. RECURSUM / ssps / min_time:0.500_mean
        impl, cls = parts[0], parts[1]
        t[(impl, cls)] = {"time": b["real_time"]}
    return t


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "results_eri.json"
    t = load(path)

    print(f"{'class':>6} {'nout':>5} {'nodes':>6} {'peak':>5} "
          f"{'RECURSUM':>10} {'libint2':>10} {'speedup':>8}   ns/integral")
    print("-" * 78)
    for cls in CLASSES:
        name = class_name(*cls)
        _, no, nn, pk = emit_kernel(*cls)
        r = t.get(("RECURSUM", name)); l = t.get(("libint2", name))
        if not r or not l:
            print(f"{name:>6}  (missing timing)")
            continue
        rt, lt = r["time"], l["time"]
        sp = lt / rt
        marker = "RECURSUM" if sp > 1 else "libint2 "
        print(f"{name:>6} {no:>5} {nn:>6} {pk:>5} "
              f"{rt:>9.2f}n {lt:>9.2f}n {sp:>7.2f}x   "
              f"R={rt/no:5.2f} L={lt/no:5.2f}  [{marker} faster]")
    print("-" * 78)
    print("Times are ns per full Cartesian class (all nout integrals), "
          "mean of 20 reps.\nBoth built g++ -O3 -march=native -ffast-math, "
          "scalar double, same quartet.")


def load_stats(path):
    import json
    d=json.load(open(path)); t={}
    for b in d["benchmarks"]:
        nm=b["name"]; 
        for suf in ("_mean","_median","_stddev"):
            if nm.endswith(suf):
                impl,cls=nm.split("/")[0],nm.split("/")[1]
                t.setdefault((impl,cls),{})[suf[1:]]=b["real_time"]
    return t

if __name__=="__main__" and __import__("sys").argv[-1]=="stats":
    import sys
    st=load_stats("results_eri.json")
    from recursum_dag_emit import class_name, CANON_LADDER
    print(f"{'class':>6} {'R_med':>9} {'R_sd':>7} {'L_med':>9} {'L_sd':>7} {'speedup':>8}")
    for cls in CANON_LADDER:
        n=class_name(*cls); r=st.get(("RECURSUM",n)); l=st.get(("libint2",n))
        if not r or not l: continue
        rm,rs=r.get("median",r.get("mean")),r.get("stddev",0)
        lm,ls=l.get("median",l.get("mean")),l.get("stddev",0)
        print(f"{n:>6} {rm:9.2f} {rs:7.3f} {lm:9.2f} {ls:7.3f} {lm/rm:7.2f}x")
