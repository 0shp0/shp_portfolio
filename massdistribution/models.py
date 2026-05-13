"""
models.py
=========
GWTC BBH 질량 분포 모델 구현 및 KS 검정.

두 모델을 공통 기반 클래스(KSTestBase) 위에 구현한다.
    - KSTruncated : 단순 멱함수 (절단 구간 적용)
    - KSPP        : 멱함수 + 가우시안 혼합 (POWER + PEAK, smoothing window 포함)

적분 방법
    - KSTruncated._model() : Gauss-Legendre 수치적분으로 CDF 직접 계산
    - KSPP._model()        : 복합 분포 형태가 복잡하여 np.trapz로 수치 정규화

데이터
    GWTC BBH 카탈로그 (GWTC-1/2/3)의 m1 중앙값.
    데이터 파일은 저작권 문제로 레포에 포함하지 않는다. data/README.md 참고.
"""

from __future__ import annotations
from abc import ABC, abstractmethod

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# ---------------------------------------------------------------------------
# 공통 기반 클래스
# ---------------------------------------------------------------------------

class KSTestBase(ABC):
    """KS 검정 기반 클래스.

    데이터 히스토그램, Gauss-Legendre 수치적분, KS 통계량 계산 등
    두 모델에 공통되는 로직을 제공한다.
    subclass는 _model()을 구현하여 모델 PDF/CDF를 정의한다.

    Parameters
    ----------
    data   : KS 검정에 사용할 관측 데이터 (1D array)
    m_min  : 질량 하한 [M☉]
    m_max  : 질량 상한 [M☉]
    binwidth : 히스토그램 bin 너비 [M☉]
    """

    def __init__(
        self,
        data:     np.ndarray,
        m_min:    float,
        m_max:    float,
        binwidth: float,
    ):
        self.data     = data
        self.m_min    = m_min
        self.m_max    = m_max
        self.binwidth = binwidth

        self.bins = np.arange(0, 110 + self.binwidth, self.binwidth, dtype=float)

        self.data_pdf, self.data_cdf, self.midpoints = self._data()
        self.model_m, self.model_pdf, self.model_cdf = self._model()
        self.df, self.max_diff_row = self._compare()

    # --- 데이터 히스토그램 ---

    def _data(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """관측 데이터를 히스토그램으로 변환하여 PDF, CDF, 구간 중점을 반환한다."""
        data_pdf, bin_edges = np.histogram(self.data, bins=self.bins, density=True)
        data_cdf, _         = np.histogram(self.data, bins=self.bins, density=True)
        data_cdf = np.cumsum(data_cdf * self.binwidth)
        midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2
        return data_pdf, data_cdf, midpoints

    # --- Gauss-Legendre 수치적분 ---

    def _integrate_gauss(self, f, m_min: float, m_max: float, N: int = 50) -> float:
        """Gauss-Legendre 구적법으로 f를 [m_min, m_max] 구간에서 적분한다.

        Parameters
        ----------
        f     : 피적분 함수 (callable)
        m_min : 적분 하한
        m_max : 적분 상한
        N     : Gauss-Legendre 노드 수 (기본값 50)
        """
        t, w = np.polynomial.legendre.leggauss(N)
        x        = 0.5 * (m_max - m_min) * t + 0.5 * (m_max + m_min)
        integral = 0.5 * (m_max - m_min) * np.sum(w * f(x))
        return integral

    # --- 모델 정의 (subclass 구현) ---

    @abstractmethod
    def _model(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """모델 PDF와 CDF를 계산하여 (model_m, model_pdf, model_cdf)를 반환한다."""

    # --- KS 통계량 ---

    def _compare(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """모델 CDF와 데이터 CDF의 차이를 계산한 DataFrame을 반환한다."""
        df = pd.DataFrame({
            "mass":      self.midpoints,
            "data_cdf":  self.data_cdf,
            "model_cdf": self.model_cdf,
        })
        df["difference"] = np.abs(df["model_cdf"] - df["data_cdf"])
        max_diff_row = df[df["difference"] == df["difference"].max()]
        return df, max_diff_row

    def ks_statistic(self) -> float:
        """KS 통계량 D_n (모델-데이터 CDF 최대 차이)을 반환한다."""
        return self.df["difference"].max()

    # --- 공통 CDF 플롯 ---

    def plot_cdf(self, model_label: str = "Model") -> None:
        """데이터와 모델의 CDF를 비교 플롯한다."""
        data_m_ext  = np.hstack([[0], self.midpoints, [110]])
        data_cdf_ext  = np.hstack([[0], self.data_cdf,  [1]])
        model_m_ext = np.hstack([[0], self.model_m,   [110]])
        model_cdf_ext = np.hstack([[0], self.model_cdf, [1]])

        plt.figure(figsize=(8, 5))
        plt.step(data_m_ext,  data_cdf_ext,  where='mid', label='GWTC BBH', color='green')
        plt.step(model_m_ext, model_cdf_ext, where='mid', label=model_label, color='blue')
        plt.xlabel(r"$m_1\ [M_\odot]$", fontsize=20)
        plt.ylabel("CDF", fontsize=20)
        plt.title(self._plot_title())
        plt.xlim(0, 110)
        plt.ylim(0, 1.1)
        plt.legend(loc="lower right", fontsize=20)
        plt.grid()
        plt.show()

    @abstractmethod
    def _plot_title(self) -> str:
        """플롯 제목 문자열 (subclass에서 파라미터 포함하여 반환)."""


# ---------------------------------------------------------------------------
# Truncated Power Law 모델
# ---------------------------------------------------------------------------

class KSTruncated(KSTestBase):
    """절단 멱함수(Truncated Power Law) 분포와 GWTC BBH m1의 KS 검정.

    [m_min, m_max] 구간에서 p(m) ∝ m^(-alpha) 이고 구간 밖은 0.
    CDF는 Gauss-Legendre 수치적분으로 계산한다.

    Parameters
    ----------
    data     : 관측 데이터
    m_min    : 멱함수 적용 하한 [M☉]
    m_max    : 멱함수 적용 상한 [M☉]
    binwidth : bin 너비 [M☉]
    alpha    : 멱함수 지수
    """

    def __init__(
        self,
        data:     np.ndarray,
        m_min:    float,
        m_max:    float,
        binwidth: float,
        alpha:    float,
    ):
        self.alpha = alpha
        super().__init__(data, m_min, m_max, binwidth)

    def _power_law(self, x: np.ndarray) -> np.ndarray:
        return x ** (-self.alpha)

    def _model(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        model_m    = self.midpoints
        m_min_idx  = np.searchsorted(model_m, self.m_min)
        m_max_idx  = np.searchsorted(model_m, self.m_max)

        # 해석적 정규화 상수
        norm = (1 - self.alpha) / (
            model_m[m_max_idx - 1] ** (1 - self.alpha)
            - model_m[m_min_idx]   ** (1 - self.alpha)
        )

        # PDF
        model_pdf = np.where(
            (model_m > self.m_min) & (model_m < self.m_max),
            norm * model_m ** (-self.alpha),
            0.0,
        )

        # CDF — Gauss-Legendre 수치적분
        model_cdf = np.zeros_like(model_m)
        model_cdf[m_min_idx:m_max_idx] = np.array([
            self._integrate_gauss(self._power_law, model_m[m_min_idx], m) * norm
            for m in model_m[m_min_idx:m_max_idx]
        ])
        model_cdf[m_max_idx:] = model_cdf[m_max_idx - 1]

        return model_m, model_pdf, model_cdf

    def _plot_title(self) -> str:
        return (
            rf"Truncated Power Law  "
            rf"$\alpha$={self.alpha},  "
            rf"$m_{{\rm min}}$={self.m_min},  "
            rf"$m_{{\rm max}}$={self.m_max},  "
            rf"bin={self.binwidth}"
        )


# ---------------------------------------------------------------------------
# Power + Peak (PP) 모델
# ---------------------------------------------------------------------------

class KSPP(KSTestBase):
    """멱함수 + 가우시안 혼합 (POWER + PEAK) 분포와 GWTC BBH m1의 KS 검정.

    p(m) ∝ [(1-peak)·m^(-alpha) + peak·N(mu, sig²)] · S(m)

    여기서 S(m)은 m_min 근방의 급격한 절단을 부드럽게 처리하는
    smooth window 함수이다.

    PDF 정규화는 np.trapz를 사용한다 (복합 분포 형태로 인해 해석적 정규화 불가).

    Parameters
    ----------
    data     : 관측 데이터
    m_min    : 분포 하한 [M☉]
    m_max    : 분포 상한 [M☉]
    dx       : bin 너비 [M☉]
    peak     : 가우시안 성분 혼합 비율 (0~1)
    alpha    : 멱함수 지수
    delta    : smooth window 폭 [M☉]
    mu       : 가우시안 평균 [M☉]
    sig      : 가우시안 표준편차 [M☉]
    """

    def __init__(
        self,
        data:  np.ndarray,
        m_min: float,
        m_max: float,
        dx:    float,
        peak:  float,
        alpha: float,
        delta: float,
        mu:    float,
        sig:   float,
    ):
        self.dx    = dx
        self.peak  = peak
        self.alpha = alpha
        self.delta = delta
        self.mu    = mu
        self.sig   = sig
        super().__init__(data, m_min, m_max, binwidth=dx)

    def _smooth_window(self, x: np.ndarray) -> np.ndarray:
        """m_min 근방을 부드럽게 절단하는 window 함수 S(m).

        x < m_min          : S = 0
        m_min ≤ x < m_min+delta : S = [exp(...) + 1]^-1  (logistic 형태)
        x ≥ m_min+delta    : S = 1
        """
        S     = np.ones_like(x)
        below = x <= self.m_min
        trans = (x > self.m_min) & (x < self.m_min + self.delta)

        S[below] = 0.0
        S[trans] = (
            np.exp(
                self.delta / (x[trans] - self.m_min)
                + self.delta / (x[trans] - self.m_min - self.delta)
            ) + 1
        ) ** -1
        return S

    def pp_pdf(self, x: np.ndarray) -> np.ndarray:
        """PP 모델 PDF를 계산한다.

        각 성분(멱함수, 가우시안)을 smooth window를 적용한 x 범위에서
        개별 정규화한 뒤 peak 비율로 혼합하고, 전체를 다시 정규화한다.
        """
        S  = self._smooth_window(x)

        gauss_raw    = np.exp(-0.5 * ((x - self.mu) / self.sig) ** 2)
        powerlaw_raw = x ** (-self.alpha)

        # 각 성분을 smooth window 적용 범위에서 정규화
        gauss    = gauss_raw    / np.trapz(gauss_raw    * S, x)
        powerlaw = powerlaw_raw / np.trapz(powerlaw_raw * S, x)

        fx  = ((1 - self.peak) * powerlaw + self.peak * gauss) * S
        fx /= np.trapz(fx, x)
        return fx

    def _model(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        model_m   = self.midpoints
        model_pdf = self.pp_pdf(model_m)

        # trapz 정규화
        area      = np.trapz(model_pdf, model_m)
        model_pdf = model_pdf / area
        model_cdf = np.cumsum(model_pdf * self.dx)
        model_cdf = model_cdf / model_cdf[-1]

        return model_m, model_pdf, model_cdf

    def plot_pdf(self) -> None:
        """모델 PDF와 데이터 히스토그램을 비교 플롯한다."""
        plt.figure(figsize=(8, 5))
        plt.plot(self.midpoints, self.model_pdf, label='PP model', color='blue')
        plt.hist(
            self.data, bins=self.bins, density=True,
            alpha=0.4, label='GWTC BBH', color='green',
        )
        plt.xlabel(r"$m_1\ [M_\odot]$", fontsize=20)
        plt.ylabel("PDF", fontsize=20)
        plt.title(self._plot_title())
        plt.legend()
        plt.grid()
        plt.show()

    def _plot_title(self) -> str:
        return (
            rf"Power + Peak  "
            rf"$\alpha$={self.alpha},  peak={self.peak},  "
            rf"$\mu$={self.mu},  $\sigma$={self.sig},  "
            rf"$m_{{\rm min}}$={self.m_min},  "
            rf"$m_{{\rm max}}$={self.m_max},  "
            rf"bin={self.dx}"
        )