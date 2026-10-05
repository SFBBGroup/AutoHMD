import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- config
# Replicates are auto-detected: any number of rep*/analysis/iRMSD*.dat
REP_GLOB   = 'rep*/analysis/iRMSD*.dat'
MEAN_CUTOFF = 2.5    # mean iRMSD <= this -> predicted high-quality
MAX_CUTOFF  = 6.0    # max  iRMSD <= this -> predicted high-quality
MEAN_ACC    = 90.0   # expected accuracy of the mean cutoff (%) — reminder only
MAX_ACC     = 83.3   # expected accuracy of the max cutoff (%)  — reminder only
TOTAL_TIME_NS = 70.0
LOG_FILE   = 'iRMSD_analysis.log'
CITATION   = ("AutoHMD: Scalable and Accurate Validation of Antibody-Antigen "
              "Complexes by Heated Molecular Dynamics")

BASE_COLORS = ['darkorange', 'darkviolet', 'steelblue', 'seagreen', 'crimson',
               'goldenrod', 'teal', 'slateblue', 'sienna', 'deeppink']

# ---------------------------------------------------------------- find replicates
def rep_num(path):
    m = re.search(r'rep(\d+)', path)
    return int(m.group(1)) if m else 0

files = sorted(glob.glob(REP_GLOB), key=rep_num)
if not files:
    raise SystemExit(f"No replicate files matched '{REP_GLOB}' in {os.getcwd()}")

irmsd_list, labels, colors = [], [], []
for i, f in enumerate(files):
    data = np.loadtxt(f, skiprows=1)
    y = data if data.ndim == 1 else data[:, 1]
    if len(y) == 0:
        print(f"skipping empty file: {f}")
        continue
    irmsd_list.append(y)
    labels.append(f"Replica {rep_num(f)}")
    colors.append(BASE_COLORS[i % len(BASE_COLORS)])
print(f"{len(irmsd_list)} replicate(s): {', '.join(labels)}")

# ---------------------------------------------------------------- plot
plt.rcParams['legend.fontsize'] = 13
plt.figure(figsize=(10, 6))
for irmsd, lab, col in zip(irmsd_list, labels, colors):
    t = np.linspace(0, TOTAL_TIME_NS, len(irmsd))
    plt.plot(t, irmsd, label=lab, color=col)

# temperature-phase vertical lines
for x in (30, 42.5, 55):
    plt.axvline(x=x, color='black', linestyle='--')

# classification cutoffs (grey): dashed = mean, dotted = max
plt.axhline(y=MEAN_CUTOFF, color='0.45', linestyle='--', linewidth=1.6)
plt.axhline(y=MAX_CUTOFF,  color='0.65', linestyle=':',  linewidth=1.8)

for pos, temp in zip([15, 36, 48, 62], ['310 K', '330 K', '360 K', '390 K']):
    plt.text(pos, 11.5, temp, color='black', va='top', ha='center', fontsize=22)

plt.xlabel('Time (ns)', fontsize=24, fontweight='bold')
plt.ylabel('iRMSD (Å)', fontsize=24, fontweight='bold')
plt.xticks(np.arange(0, TOTAL_TIME_NS + 5.0, 5.0), fontsize=22, rotation=45)
plt.yticks(fontsize=22)
plt.ylim(0, 12)
plt.xlim(0, TOTAL_TIME_NS)
plt.legend(loc='upper right', bbox_to_anchor=(0.82, 1.20),
           ncol=min(len(labels), 5), fontsize=16)
plt.tight_layout()
plt.savefig('iRMSD_replicas.png', transparent=True, dpi=400, format='png', pad_inches=0.5)
plt.grid(False)

# ---------------------------------------------------------------- log
per_rep = []
for lab, irmsd in zip(labels, irmsd_list):
    per_rep.append(dict(
        label=lab, rmax=float(np.max(irmsd)), rmean=float(np.mean(irmsd)),
        rmin=float(np.min(irmsd)), rstd=float(np.std(irmsd)),
        viol_mean=100.0 * np.mean(irmsd > MEAN_CUTOFF),
        viol_max=100.0 * np.mean(irmsd > MAX_CUTOFF)))

pose_mean = np.mean([r['rmean'] for r in per_rep])
pose_max  = np.mean([r['rmax'] for r in per_rep])
pred_mean = 'high-quality' if pose_mean <= MEAN_CUTOFF else 'incorrect'
pred_max  = 'high-quality' if pose_max  <= MAX_CUTOFF  else 'incorrect'

with open(LOG_FILE, 'w') as log:
    log.write(f"iRMSD analysis — {os.path.basename(os.getcwd())}\n")
    log.write(f"replicates ({len(per_rep)}): {', '.join(labels)}\n")
    log.write(f"cutoffs: mean <= {MEAN_CUTOFF:.2f} A, max <= {MAX_CUTOFF:.2f} A "
              f"(statistics over the full trajectory, all heating phases)\n\n")
    log.write("per-replicate statistics (A) and cutoff violation (% of frames)\n")
    log.write(f"{'replica':10}{'max':>8}{'mean':>8}{'min':>8}{'std':>8}"
              f"{'%>'+str(MEAN_CUTOFF):>9}{'%>'+str(MAX_CUTOFF):>9}\n")
    for r in per_rep:
        log.write(f"{r['label']:10}{r['rmax']:8.2f}{r['rmean']:8.2f}{r['rmin']:8.2f}"
                  f"{r['rstd']:8.2f}{r['viol_mean']:9.1f}{r['viol_max']:9.1f}\n")
    log.write("\npose-level prediction (statistic averaged over replicates)\n")
    log.write(f"  mean iRMSD = {pose_mean:.2f} A  ->  {pred_mean}  "
              f"(cutoff {MEAN_CUTOFF:.1f} A)\n")
    log.write(f"  max  iRMSD = {pose_max:.2f} A  ->  {pred_max}  "
              f"(cutoff {MAX_CUTOFF:.1f} A)\n")
    log.write(f"\nreminder: at these pre-established cutoffs the expected "
              f"classification accuracy is {MEAN_ACC:.1f}% (mean) and "
              f"{MAX_ACC:.1f}% (max).\n")
    log.write(f"\nif you use this analysis, please cite:\n  {CITATION}\n")

print(f"saved iRMSD_replicas.png and {LOG_FILE}")
print(f"\nIf you use this analysis, please cite:\n  {CITATION}")
# plt.show()  # uncomment for interactive use
