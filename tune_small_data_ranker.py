#!/usr/bin/env python3
"""Round-two model-selection diagnostic; never used as discovery evidence."""

import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor, HistGradientBoostingRegressor

from pipeline.catalyst_spaces import encode_population
from pipeline.discovery import candidate_id
from run_divide_conquer_pilot import _training_frame


from run_divide_conquer_pilot import ROUND
manifest = json.loads(Path(f"results/pilot/divide_conquer_{ROUND}/manifest.json").read_text())
paths = {"turquoise_hydrogen": Path(f"results/screening/pilot/divide_conquer_{ROUND}_pyrolysis.csv"),
         "fuel_cell_orr": Path(f"results/fuel_cell/pilot/divide_conquer_{ROUND}_orr.csv")}
results = []
for record in manifest["records"]:
    application = record["application"]
    outcome = "E_act" if application == "turquoise_hydrogen" else "orr_overpotential_V"
    train = _training_frame(application)
    train = train[np.isfinite(pd.to_numeric(train[outcome], errors="coerce"))]
    train_genomes = [ast.literal_eval(g) if isinstance(g, str) else tuple(g) for g in train.genome]
    x_train = encode_population(train_genomes)
    y_train = train[outcome].to_numpy(float)
    test = pd.read_csv(paths[application])
    test = test[test.valid.eq(True)].copy()
    test_genomes = [ast.literal_eval(g) for g in test.genome]
    x_test = encode_population(test_genomes)
    y_test = test[outcome].to_numpy(float)
    cutoff = float(np.quantile(y_test, .20))
    budget = record["budget"]
    models = []
    for leaf in (1, 2, 3, 5):
        for features in ("sqrt", 0.5, 1.0):
            models.append((f"extra_leaf{leaf}_feat{features}", ExtraTreesRegressor(
                n_estimators=512, min_samples_leaf=leaf, max_features=features,
                random_state=20260720, n_jobs=-1)))
    models.extend([
        ("random_forest", RandomForestRegressor(n_estimators=512, min_samples_leaf=2,
                                                random_state=20260720, n_jobs=-1)),
        ("hist_gradient", HistGradientBoostingRegressor(max_iter=200, l2_regularization=2.0,
                                                        random_state=20260720)),
    ])
    for name, model in models:
        model.fit(x_train, y_train)
        prediction = model.predict(x_test)
        selected = np.argsort(prediction)[:budget]
        results.append({"application": application, "model": name,
                        "spearman": float(spearmanr(prediction, y_test).statistic),
                        "hits": int(np.sum(y_test[selected] <= cutoff)),
                        "mean_selected": float(y_test[selected].mean())})
    if application == "fuel_cell_orr":
        adsorption = ["dG_OH_eV", "dG_O_eV", "dG_OOH_eV"]
        finite = train[adsorption].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
        x_ads = x_train[finite.to_numpy()]
        y_ads = train.loc[finite, adsorption].to_numpy(float)
        for leaf in (1, 2, 3, 5):
            model = ExtraTreesRegressor(n_estimators=1024, min_samples_leaf=leaf,
                                        max_features=1.0, random_state=20260721, n_jobs=-1)
            model.fit(x_ads, y_ads)
            d_oh, d_o, d_ooh = model.predict(x_test).T
            prediction = 1.23 + np.maximum.reduce([d_ooh - 4.92, d_o - d_ooh,
                                                   d_oh - d_o, -d_oh])
            selected = np.argsort(prediction)[:budget]
            results.append({"application": application, "model": f"adsorption_extra_leaf{leaf}",
                            "spearman": float(spearmanr(prediction, y_test).statistic),
                            "hits": int(np.sum(y_test[selected] <= cutoff)),
                            "mean_selected": float(y_test[selected].mean())})

Path(f"results/pilot/divide_conquer_{ROUND}/model_selection.json").write_text(
    json.dumps(results, indent=2) + "\n")
for application in ("turquoise_hydrogen", "fuel_cell_orr"):
    subset = [r for r in results if r["application"] == application]
    print(application)
    for row in sorted(subset, key=lambda r: (-r["hits"], -r["spearman"]))[:6]:
        print(row)
