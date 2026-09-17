"""
Predict 3-year survival (died within 3 years vs. survived past 3 years)
from the top 50 Factor 6/12 genes, using Random Forest and XGBoost.

The label is imbalanced (more survivors than deaths), so this script
tries three standard ways of dealing with that and compares them on
accuracy, precision and recall, using 5-fold stratified cross-validation:
  1. class_weight  - tell the model to care more about the minority class
  2. oversampling  - duplicate/synthesize minority-class rows (SMOTE)
  3. undersampling - drop some majority-class rows
"""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyreadr

from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score
from xgboost import XGBClassifier
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler

DATA_DIR = "../data"
OUT_DIR = "../Outputs/05_ml_classification"
RANDOM_STATE = 42
STRATEGIES = ["class_weight", "oversample", "undersample"]

os.makedirs(OUT_DIR, exist_ok=True)

# ── Step 1: load the gene expression features and the survival label ──

# 386 patients x 50 genes (the Factor 6/12 genes from 02_mofa_integration.Rmd)
genes = pd.read_csv(os.path.join(DATA_DIR, "top_extracted50_genes.csv"), index_col=0)

# clin_final.rds is an R data frame - pyreadr reads it straight into pandas
clinical = pyreadr.read_r(os.path.join(DATA_DIR, "clin_final.rds"))[None]
clinical = clinical.set_index("submitter_id")

# label: 0 = died within 3 years, 1 = survived past 3 years (see 01_QC_EDA.Rmd)
# some patients have NaN (not enough follow-up to call it either way) - drop those
data = genes.join(clinical["label"], how="inner").dropna(subset=["label"])
X = data.drop(columns="label")
y = data["label"].astype(int)

# ── Step 2: EDA - how imbalanced is the label, really? ──

counts = y.value_counts().sort_index()
percentages = (counts / len(y) * 100).round(1)

print(f"patients with a usable label: {len(y)}")
print("class balance (0 = died <3yr, 1 = survived >3yr):")
for label_value in counts.index:
    print(f"  class {label_value}: {counts[label_value]} patients ({percentages[label_value]}%)")

plt.figure(figsize=(4, 4))
counts.plot.bar(color=["#C0392B", "#2E4A7A"])
plt.title("Class Balance - 3-Year Survival Label")
plt.xlabel("Class (0 = died <3yr, 1 = survived >3yr)")
plt.ylabel("Number of Patients")
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "class_balance.png"), dpi=150)
print("saved: class_balance.png\n")


# ── Helper 1: apply one imbalance strategy to a training fold ──
# "class_weight" leaves the rows alone - the model itself upweights the
# minority class instead (see build_model below). The other two actually
# add/remove rows, and are only ever applied to the training fold.
def resample(strategy, X_train, y_train):
    if strategy == "oversample":
        return SMOTE(random_state=RANDOM_STATE).fit_resample(X_train, y_train)
    if strategy == "undersample":
        return RandomUnderSampler(random_state=RANDOM_STATE).fit_resample(X_train, y_train)
    return X_train, y_train


# ── Helper 2: build a fresh, untrained model for one (model, strategy) pair ──
def build_model(model_name, strategy, y_train):
    use_class_weight = strategy == "class_weight"
    if model_name == "RandomForest":
        return RandomForestClassifier(
            class_weight="balanced" if use_class_weight else None,
            random_state=RANDOM_STATE,
        )
    n_died, n_survived = y_train.value_counts().sort_index()
    return XGBClassifier(
        scale_pos_weight=(n_died / n_survived) if use_class_weight else 1,
        eval_metric="logloss",
        random_state=RANDOM_STATE,
    )


# ── Step 3: 5-fold stratified cross-validation ──
# each fold keeps the same death/survival ratio as the full data, and
# resampling is rebuilt inside every fold using only that fold's training
# rows - the validation rows are always left untouched.

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
fold_results = []

for fold_num, (train_idx, val_idx) in enumerate(cv.split(X, y), start=1):
    X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
    y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

    for strategy in STRATEGIES:
        X_tr, y_tr = resample(strategy, X_train, y_train)

        for model_name in ["RandomForest", "XGBoost"]:
            model = build_model(model_name, strategy, y_tr)
            model.fit(X_tr, y_tr)
            pred = model.predict(X_val)

            fold_results.append({
                "fold":      fold_num,
                "strategy":  strategy,
                "model":     model_name,
                "accuracy":  accuracy_score(y_val, pred),
                "precision": precision_score(y_val, pred),
                "recall":    recall_score(y_val, pred),
            })

print(f"ran {cv.get_n_splits()} folds x {len(STRATEGIES)} strategies x 2 models\n")

# ── Step 4: average performance across the 5 folds ──

fold_results_df = pd.DataFrame(fold_results)
results_df = (
    fold_results_df.groupby(["strategy", "model"])[["accuracy", "precision", "recall"]]
    .mean()
    .reset_index()
    .sort_values("accuracy", ascending=False)
)
print("average performance across 5 folds:")
print(results_df.to_string(index=False))

# ── Step 5: plot the averaged results ──

results_df["combo"] = results_df["model"] + " + " + results_df["strategy"]
results_df.set_index("combo")[["accuracy", "precision", "recall"]].plot.bar(
    figsize=(9, 5), color=["#2E4A7A", "#C0392B", "#27AE60"]
)
plt.title("Classification Performance by Model and Imbalance Strategy")
plt.ylabel("Score")
plt.xlabel("")
plt.xticks(rotation=30, ha="right")
plt.ylim(0, 1)
plt.legend(title="Metric")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "classification_results.png"), dpi=150)
print("\nsaved: classification_results.png")
