"""Feature-map kernels: raw, increments, Hermite projection, DCT high-pass, oracle fractional difference.

Each class keeps exactly the arithmetic of the EXP7 switch it replaces (exp7_core.transform_features and
gram_full), so saved EXP7 rows reproduce.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kmem import memory_kernels as mk
from kmem.legacy import mbcore as mb


def _gauss(x: np.ndarray, sigma: float) -> np.ndarray:
    return mb.gauss_gram(np.abs(x[:, None] - x[None, :]), sigma)


@dataclass(frozen=True)
class Raw:
    name: str = "raw"
    kind: str = "feature"
    lag: int = 0

    def features(self, f: np.ndarray, x: np.ndarray, H: float | None = None) -> np.ndarray:
        return f

    def gram(self, x: np.ndarray, sigma: float, H: float | None = None) -> np.ndarray:
        return mk.double_center(_gauss(x, sigma))


@dataclass(frozen=True)
class Increment:
    """D K D' with (D f)_t = f_{t+lag} - f_t, double centred."""
    lag: int
    name: str = ""
    kind: str = "feature"

    def __post_init__(self):
        if not self.name:
            object.__setattr__(self, "name", f"inc{self.lag}")

    def features(self, f: np.ndarray, x: np.ndarray, H: float | None = None) -> np.ndarray:
        return f[self.lag:] - f[:-self.lag]

    def gram(self, x: np.ndarray, sigma: float, H: float | None = None) -> np.ndarray:
        return mk.increment(_gauss(x, sigma), self.lag)


@dataclass(frozen=True)
class Hermite:
    """Q_r K Q_r: removes 1, He_1(z), ..., He_r(z) of the normal scores z of x."""
    r: int
    name: str = ""
    kind: str = "feature"
    lag: int = 0

    def __post_init__(self):
        if not self.name:
            object.__setattr__(self, "name", f"herm{self.r}")

    def features(self, f: np.ndarray, x: np.ndarray, H: float | None = None) -> np.ndarray:
        u = mk.hermite_basis(x, self.r)
        return f - u @ (u.T @ f)

    def gram(self, x: np.ndarray, sigma: float, H: float | None = None) -> np.ndarray:
        return mk.hermite(_gauss(x, sigma), x, self.r)


@dataclass(frozen=True)
class HighPass:
    """F K F: removes the lowest frac of DCT-II frequencies, the constant included."""
    frac: float
    name: str
    kind: str = "feature"
    lag: int = 0

    def features(self, f: np.ndarray, x: np.ndarray, H: float | None = None) -> np.ndarray:
        u = mk.dct_basis(x.size, mk.n_low(x.size, self.frac))
        return f - u @ (u.T @ f)

    def gram(self, x: np.ndarray, sigma: float, H: float | None = None) -> np.ndarray:
        return mk.highpass(_gauss(x, sigma), self.frac)


@dataclass(frozen=True)
class FracOracle:
    """Fractional difference of order d = H - 1/2 with the true H (an oracle: H must be passed)."""
    name: str = "frac_oracle"
    kind: str = "feature"
    lag: int = 0

    def features(self, f: np.ndarray, x: np.ndarray, H: float | None = None) -> np.ndarray:
        if H is None:
            raise ValueError("frac_oracle needs the true H")
        d = H - 0.5
        return f - f.mean(axis=0, keepdims=True) if d == 0.0 else mk.fracdiff_features(f - f.mean(axis=0, keepdims=True), d)

    def gram(self, x: np.ndarray, sigma: float, H: float | None = None) -> np.ndarray:
        return mk.fracdiff(_gauss(x, sigma), H - 0.5)
