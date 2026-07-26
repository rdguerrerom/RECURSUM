"""Direct measurement of each kernel's STACK FRAME from the compiled object --
the concrete realization of the peak-liveness slot count. Reads the
`sub $0xNNN,%rsp` prologue of recursum_eri_<class>[_naive]."""
import re
import subprocess
import sys

from recursum_dag_emit import emit_kernel, class_name, CANON_LADDER


def frame_bytes(obj, fnname):
    try:
        d = subprocess.run(["objdump", "-C", "-d", obj],
                           capture_output=True, text=True).stdout
    except FileNotFoundError:
        return None
    lines = d.splitlines()
    inside = False
    label = re.compile(r"<.*\brecursum_eri_%s\b" % re.escape(fnname))
    biggest = 0
    for ln in lines:
        if label.search(ln) and ln.rstrip().endswith(">:"):
            inside = True
            continue
        if inside:
            if re.match(r"^[0-9a-f]+ <", ln):   # next function label
                break
            m = re.search(r"sub\s+\$0x([0-9a-f]+),%rsp", ln)
            if m:
                biggest = max(biggest, int(m.group(1), 16))
    return biggest


if __name__ == "__main__":
    print(f"{'class':>6} {'peakON':>7} {'peakNv':>7} {'frameON':>9} {'frameNv':>9} "
          f"{'frame_reduce':>12}")
    for cls in CANON_LADDER:
        n = class_name(*cls)
        pk = emit_kernel(*cls, reuse=True)[3]
        nv = emit_kernel(*cls, reuse=False)[3]
        fon = frame_bytes(f"k_{n}.o", n)
        fnv = frame_bytes(f"k_{n}_naive.o", n + "_naive")
        red = (f"{100*(1-fon/fnv):.0f}%" if (fon and fnv) else "-")
        fnv_s = "-" if fnv is None else str(fnv)
        print(f"{n:>6} {pk:>7} {nv:>7} {str(fon):>9} {fnv_s:>9} {red:>12}")
    print("\nframe = bytes subtracted from %rsp in the kernel prologue "
          "(stack working set). ON = peak-liveness slot reuse; Nv = naive "
          "(one slot per intermediate). ffff has no naive variant (ON-only).")
