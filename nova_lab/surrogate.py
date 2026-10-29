"""Surrogate bundle: main LightGBM model plus optional auxiliary-target models.

Auxiliary models are trained on the large score sets of past targets (scores
recomputed with today's formula). Boltz-2 affinity has a strong ligand-only
component, so their predictions transfer as extra features.
"""
from dataclasses import dataclass, field

import joblib
import numpy as np

from . import DATA_DIR


@dataclass
class Surrogate:
    main: object
    aux: dict = field(default_factory=dict)  # target -> model

    def augment(self, X: np.ndarray) -> np.ndarray:
        if not self.aux:
            return X
        extra = np.column_stack([m.predict(X) for m in self.aux.values()])
        return np.hstack([X, extra.astype(X.dtype)])

    def predict(self, X: np.ndarray) -> np.ndarray:
        return self.main.predict(self.augment(X))

    def save(self, target: str):
        joblib.dump(self, DATA_DIR / "dataset" / f"{target}_surrogate.joblib")

    @staticmethod
    def load(target: str) -> "Surrogate":
        return joblib.load(DATA_DIR / "dataset" / f"{target}_surrogate.joblib")
