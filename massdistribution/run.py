"""
run.py
======
GWTC BBH m1 분포 피팅 및 KS 검정 실행 스크립트.

데이터 파일 경로를 수정하여 사용한다. (data/README.md 참고)
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from models import KSTruncated, KSPP

plt.rcParams['font.family']      = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

# ---------------------------------------------------------------------------
# 데이터 로드
# ---------------------------------------------------------------------------

GWTC_BBH = pd.read_csv("./data/GWTC_BBH.txt", sep=r'\s+', encoding='utf-8')
BBH_m1   = GWTC_BBH['m1_median'].round(1).to_numpy()

# ---------------------------------------------------------------------------
# Truncated Power Law
# ---------------------------------------------------------------------------

trunc = KSTruncated(
    data     = BBH_m1,
    m_min    = 2,
    m_max    = 90,
    binwidth = 0.1,
    alpha    = 1.5,
)

print("=== Truncated Power Law ===")
print(trunc.df)
print("KS statistic D_n =", trunc.ks_statistic())
trunc.plot_cdf(model_label="Truncated")

# ---------------------------------------------------------------------------
# Power + Peak (PP)
# ---------------------------------------------------------------------------

pp = KSPP(
    data  = BBH_m1,
    m_min = 5,
    m_max = 40,
    dx    = 1.0,
    peak  = 0.1,
    alpha = 2.3,
    delta = 2.0,
    mu    = 40.0,
    sig   = 2.0,
)

print("\n=== Power + Peak ===")
print(pp.df)
print("KS statistic D_n =", pp.ks_statistic())
pp.plot_pdf()
pp.plot_cdf(model_label="PP")