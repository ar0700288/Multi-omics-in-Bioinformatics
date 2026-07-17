"""
Random Survival Forest (RSF) on the same top 50 Factor 6/12 genes and the
same time/event data used by the Cox model in 04_survival_analysis.Rmd.
Unlike Cox, RSF doesn't assume each gene has a linear effect on risk.
Fits the model, reports its concordance index, and saves a one-bar image
of that score - no survival curves or feature importance.
"""

import os
import pandas as pd
import pyreadr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv

DATA_DIR = "../data"
OUT_DIR = "../Outputs/06_random_survival_forest"
RANDOM_STATE = 42

os.makedirs(OUT_DIR, exist_ok=True)

# 386 patients x 50 genes, patient ID is the row index
genes = pd.read_csv(f"{DATA_DIR}/top_extracted50_genes.csv", index_col=0)

# clin_final.rds is an R data frame - pyreadr reads it straight into pandas
clinical = pyreadr.read_r(f"{DATA_DIR}/clin_final.rds")[None]
clinical = clinical.set_index("submitter_id")

# same time/event definition as the Cox model: time to death if the
# patient died, otherwise time followed up for; event = 1 means died
clinical["time"] = clinical["days_to_death"].where(
    clinical["vital_status"] == "Dead", clinical["days_to_last_follow_up"]
)
clinical["event"] = (clinical["vital_status"] == "Dead").astype(int)

# combine genes with time/event, drop patients with no usable follow-up time
data = genes.join(clinical[["time", "event"]], how="inner").dropna(subset=["time"])
data = data[data["time"] > 0]
print(f"patients used for RSF: {len(data)}")

X = data.drop(columns=["time", "event"])
# sksurv wants (event, time) packed into one structured array, not two columns
y = Surv.from_arrays(event=data["event"].astype(bool), time=data["time"])

# stratify on event so both splits have a similar death/censoring ratio
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, stratify=data["event"], random_state=RANDOM_STATE
)
print(f"train: {len(y_train)} patients | test: {len(y_test)} patients")

rsf = RandomSurvivalForest(n_estimators=200, min_samples_leaf=10, random_state=RANDOM_STATE, n_jobs=-1)
rsf.fit(X_train, y_train)

# concordance index - the survival-analysis version of accuracy, i.e. how
# often the model ranks two patients' risk in the correct order
c_index = rsf.score(X_test, y_test)
print(f"concordance index (test set): {c_index:.3f}")

plt.figure(figsize=(3, 4))
plt.bar(["Concordance\nIndex"], [c_index], color="#2E4A7A", width=0.5)
plt.axhline(0.5, linestyle="--", color="#C0392B", linewidth=1)   # 0.5 = random guessing
plt.ylim(0, 1)
plt.text(0, c_index + 0.02, f"{c_index:.3f}", ha="center", fontweight="bold")
plt.title("Random Survival Forest\nTest-Set Accuracy")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "concordance_index.png"), dpi=150)
print("saved: concordance_index.png")
