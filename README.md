# shp_portfolio

**한국어** | [English](#english)

물리학 석사 과정에서 수행한 수치 천체물리 및 중력파 데이터 분석 관련 코드 모음입니다.

---

## 한국어

### 프로젝트 목록

#### [2bodymerger](./2bodymerger)
**4차 Hermite predictor-corrector를 이용한 쌍성 궤도 시뮬레이션 및 중력파 방출에 의한 합체 계산**

2024 GW-NR 겨울학교 수치상대론 부문 문제풀이 코드.
뉴턴 중력 궤도적분과 2.5PN 복사 반력 항을 구현하여 블랙홀 쌍성의 나선형 합체 과정을 시뮬레이션한다.
Hermite 적분기가 요구하는 PN jerk 항을 직접 유도하여 구현.

`Python` · `Numerical Integration` · `Post-Newtonian`

---

#### [massdistribution](./massdistribution)
**GWTC BBH 질량 분포 모델 구현 및 KS 검정**

동기의 중력파 데이터 분석 연구 지원 과정에서 두 질량 분포 모델 구현을 전담.
Abbott et al. (2021) 의 Truncated Power Law 및 Power+Peak 모델을 구현하고
GWTC 카탈로그 데이터와 KS 검정으로 적합도를 평가한다.

`Python` · `Statistical Modeling` · `Gravitational Wave Astronomy`

---

#### [ET_matched_filter](./ET_matched_filter)
**Einstein Telescope MDC 데이터에 대한 matched filtering 및 network SNR 계산**

Einstein Telescope International R&D 활동의 일환으로 수행한 Mock Data Challenge(MDC) 분석.
PyCBC를 이용해 BBH/BHNS/NSNS 이벤트에 대한 CBC 템플릿을 생성하고,
3개 검출기(E1, E2, E3)의 antenna pattern을 적용한 matched filtering으로 network SNR을 계산한다.

`Python` · `PyCBC` · `Matched Filtering` · `Einstein Telescope`

---

## English

Graduate-level numerical astrophysics and gravitational wave data analysis projects.

### Projects

#### [2bodymerger](./2bodymerger)
**Binary star orbit simulation and GW-driven merger via 4th-order Hermite predictor-corrector**

Competition solution for the 2024.01 GW-NR Winter School. (Numerical Relativity part)
Implements Newtonian orbit integration and 2.5PN radiation-reaction to simulate
the inspiral merger of black hole binaries.
The PN jerk term required by the Hermite integrator was analytically derived and implemented.

`Python` · `Numerical Integration` · `Post-Newtonian`

---

#### [massdistribution](./massdistribution)
**GWTC BBH mass distribution modeling and KS test**

Implemented two BBH mass distribution models in support of a collaborator's GW data analysis research.
Based on Abbott et al. (2021): Truncated Power Law and Power+Peak models are fitted
to the GWTC catalog and evaluated via the Kolmogorov-Smirnov test.

`Python` · `Statistical Modeling` · `Gravitational Wave Astronomy`

---

#### [ET_matched_filter](./ET_matched_filter)
**Matched filtering and network SNR calculation on Einstein Telescope MDC data**

Analysis performed as part of the Einstein Telescope International R&D activity (Mock Data Challenge).
Generates CBC waveform templates (IMRPhenomXPHM / IMRPhenomPv2_NRTidal) using PyCBC,
applies antenna patterns for the 3-detector ET network (E1, E2, E3),
and computes network SNR via matched filtering.

`Python` · `PyCBC` · `Matched Filtering` · `Einstein Telescope`