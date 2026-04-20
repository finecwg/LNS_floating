# Floating Object Solver v2 — 발표 노트 (한국어)

> **대상 독자**: 기본 미적분학과 선형대수만 아는 청중
> **목적**: 발표 중 참고용 + 교육 자료
> **스타일**: 각 슬라이드마다 "핵심 메시지 → 배경지식 → 자세한 설명 → 발표 멘트" 구조

---

## Slide 1 — Title

**핵심 메시지**: v2 솔버는 유체 위의 부유체(floating body)의 운동을 "단일(monolithic) 제약조건" 방식으로 일관되게 푸는 방법이다.

**발표 멘트**:
> "안녕하세요. 오늘은 Floating Object Solver v2에 대해 발표합니다. 이 v2는 기존 v1과 달리, 유체 위 강체의 운동을 스스로 결정하게 하는 방식으로 — 즉 외부에서 운동을 미리 지정하지 않고 뉴턴 제2법칙을 풀어서 — 계산합니다. 핵심 키워드는 'self-consistent body dynamics'와 'monolithic constraint solve' 입니다."

**배경**: v1은 물체의 위치를 $z_b(t) = -\text{draft} + A\sin(\omega t)$처럼 *미리 지정*했지만, v2에서는 물체의 무게와 유체 압력이 스스로 물체 운동을 결정한다.

---

## Slide 2 — Outline

**발표 멘트**:
> "오늘 발표는 여섯 부분으로 구성됩니다. 먼저 왜 v2가 필요한지 동기를 설명하고, 물리 모델을 보여드린 뒤, 수치적인 방식, 코드 구현, 구체적 예시, 그리고 한계와 향후 과제 순입니다. 강조하고 싶은 점은 — 발표하는 모든 수식이 실제로 `floating_object_v2.py` 코드 안에 들어있다는 것입니다. 논리와 코드가 1:1로 대응됩니다."

**시간 배분 권장**: 총 30~45분이면 각 파트 5~8분씩.

---

## Slide 3 — Problem setup (문제 설정)

**문제 상황**:
- 2차원 유체(위: 공기, 아래: 물), 수평 방향은 주기적(periodic in $x$)
- 밑바닥은 고정벽 ($z=-D$), 비미끄럼(no-slip) 조건
- 유체 표면 위에 **평평한 바닥을 가진 강체** 하나가 떠 있음
- 반폭 $a$, 단위 길이당 질량 $m$ (2D 문제이므로 "단위 out-of-plane 길이당" 질량)
- $t=0$에서 물체가 유체 표면에 막 접촉한 상태

**우리가 구하고 싶은 것**:
- $z_b(t)$: 물체 바닥의 높이 (시간에 따라)
- $v_b(t)$: 물체 연직 속도
- $\eta(x,t)$: 유체 표면 변위
- $p_s(x,t)$: 표면 압력 — 특히 물체 아래 영역에서

**청중을 위한 설명**:
> "2D로 단순화한 이유는, 구형 물체나 3D 강체 대신 '연직 방향 운동만' 집중하기 위해서입니다. 물이 양 옆으로 무한히 뻗어있다고 가정하기 위해 주기적 경계조건을 씁니다. 실험으로 치면 긴 수조 중앙에 나무토막을 올려놓고 가라앉는 것을 관찰하는 셈이죠."

---

## Slide 4 — v1 vs v2 at a glance (한눈에 비교)

**핵심 차이**:

| 항목 | v1 | v2 |
|---|---|---|
| 물체 운동 | 사인파로 **강제 진동** | 중력 + 유체 압력으로 **스스로 결정** |
| 물체 역할 | 파도 생성기 (wavemaker) | 진짜 강체 역학 |
| 표면 압력 $p_s$ | 계산 안 함 | **라그랑주 승수** (핵심!) |
| 물체 아래 BC | $\eta \leftarrow z_b(t)$로 덮어쓰기 | 기하·운동·수평 no-slip 3가지 제약 |
| 불안정성 위험 | 강제 진동이라 안정 | Fully implicit → added-mass 문제 없음 |

**"Key insight" 박스 설명**:
> "여기가 이 발표의 핵심 아이디어입니다: 물체 운동을 풀어넣어도, 새로 생긴 식들이 전부 **선형이고 계수가 상수**여서 — v1에서 쓰던 'LU 분해 한 번만' 트릭이 그대로 작동합니다. 이 때문에 implicit 풀이가 비싸지 않습니다."

**Added-mass 용어 설명 (질문 나올 수 있음)**: 물체가 유체 속에서 움직일 때, 물체만이 아니라 주변 유체도 함께 가속되어야 해서 "실제 질량보다 큰 관성"처럼 느껴지는 효과. 이 added mass를 암묵적(implicit)으로 처리하지 않으면 가벼운 물체에서는 수치적 불안정이 발생한다.

---

## Slide 5 — Helmholtz decomposition of the velocity field

**핵심 아이디어**:
- 유체의 속도장 $\mathbf{u}$를 두 부분으로 나눈다:
$$\mathbf{u} = \nabla\phi + \mathbf{w}$$
- $\nabla\phi$: 포텐셜(irrotational) 부분 → 라플라스 방정식 만족
- $\mathbf{w}$: 와류(vortical) 부분 → 발산이 0 ($\nabla\cdot\mathbf{w}=0$)

**왜 이렇게 쪼개나?**:
- 전체 Navier-Stokes 방정식은 비선형이라 매우 어려움
- 대부분의 영역(bulk)에서는 유체 운동이 거의 비회전성(potential flow)에 가까움
- 점성(viscous) 효과는 경계(표면, 바닥) 근처에만 집중됨
- 이 분해는 "Bulk는 단순한 선형 PDE로 + 경계에 점성 효과 응축"이라는 장점

**교육 포인트** (선형대수 관점):
> "벡터장을 두 성분으로 나누는 건 선형대수에서 벡터를 '기저 두 개로 분해' 하는 것과 비슷합니다. 여기서는 '회전 없는 부분'과 '발산 없는 부분'이 각각의 기저 역할을 합니다."

**$\phi$의 의미**: 속도 포텐셜. 전자기학의 전기 포텐셜처럼, **그 기울기(gradient)가 속도의 일부**를 준다.

**Reference**: Galeano-Rios, Milewski, Vanden-Broeck (JFM 2017).

---

## Slide 6 — Bulk equations (유체 내부 방정식)

**지배 방정식 (내부, $-D \le z \le 0$)**:

1. **라플라스 방정식** (포텐셜):
$$\Delta\phi = 0$$

2. **열 방정식** (와류 성분들):
$$\partial_t w^i = \frac{1}{\text{Re}}\Delta w^i, \quad i=1,3$$

**청중을 위한 설명**:

- $\Delta = \partial_{xx} + \partial_{zz}$ = 2D 라플라시안. "함수의 곡률의 합"이라고 생각하면 됨.
- $\Delta_H = \partial_{xx}$ = 수평 라플라시안. 표면에서만 사용.
- **라플라스 방정식의 직관**: "공간 어디서나 주변 평균값과 같은 함수" → 정상상태 전도, 포텐셜 흐름의 특성.
- **열 방정식의 직관**: 점성(viscosity) $\nu = 1/\text{Re}$에 비례해 와류가 시간에 따라 퍼진다(diffusion).
- **Reynolds 수 $\text{Re}$**: 관성/점성의 비율. 크면 점성이 약함(물), 작으면 점성이 강함(꿀).

**발표 멘트**:
> "Bulk에서는 이 두 가지 간단한 선형 PDE만 풀면 됩니다. 포텐셜은 시간에 독립인 라플라스 방정식, 와류는 확산(열)방정식. 2D라서 $w^1$(수평), $w^3$(연직) 두 성분만 있고, $w^2=0$은 대칭성에서 자동으로 0."

---

## Slide 7 — Free-surface boundary conditions (자유표면 BC, 물체 바깥)

**네 가지 표면 조건** ($z=0$에서, 물체 바깥 영역):

1. **운동학적 BC (2.15c)**: $\partial_t \eta = \partial_z \phi + w^3$
   - 표면이 유체 흐름을 따라 움직인다 (물 표면이 수직방향 속도로 올라감)

2. **동역학적 BC (2.15d)** with $p_s = 0$:
$$\partial_t \phi_s = -\frac{\eta}{\text{Fr}} + \frac{\kappa[\eta]}{\text{We}} + \frac{2}{\text{Re}}\Delta_H\phi_s - \frac{2}{\text{Re}}\partial_z w^3$$
   - 압력 균형: 중력(Fr), 표면장력(We), 점성(Re) 효과의 합
   - $\kappa[\eta] \approx \partial_{xx}\eta$ = 선형화된 곡률

3. **와류 BC (2.15e)**: $\partial_t w^3 = \frac{2}{\text{Re}}\Delta_H(\partial_z\phi + w^3)$
   - 표면 근처 와류의 생성/확산 방정식

4. **응력 없음 BC (2.11a)**: $\partial_z w^1 = -(\partial_{xz}\phi + \partial_x w^3)$
   - 표면에서 접선 방향 응력이 0

**"Remark" 박스 설명**:
> "많은 구현에서 간단화를 위해 (2.17) $w^3 = \frac{2}{\text{Re}}\Delta_H\eta$ 라는 **근사식**을 씁니다. 이는 초기 과도상태가 끝난 후에만 맞고, 바닥의 no-slip 조건을 제대로 다룰 수 없습니다. 이 v2는 **(2.15e)를 직접** 풀어서 더 정확하게 합니다."

**교육 포인트**: 경계조건이 많아 보여도, 물리적으로는 각각 명확한 의미가 있음:
- 운동학 = "기하학적으로 어디 있는가"
- 동역학 = "힘/압력이 어떻게 균형잡히는가"
- 와류 = "소용돌이가 어떻게 진화하는가"

---

## Slide 8 — Under the body: three new constraints (물체 밑 세 가지 제약)

**핵심**: 물체 바닥(wetted region $\mathcal{B} = \{x : |x-x_c| < a\}$)에서는 자유표면 조건이 **세 가지 제약으로 교체**된다.

**제약 (1): 기하학적 일치**
$$\eta(x,t) = z_b(t), \quad \forall x \in \mathcal{B}$$
- **의미**: 물체 바닥과 유체 표면이 **같은 위치**여야 함. 물체가 떠 있다면 유체 표면이 물체 바닥 높이에 와있다.

**제약 (2): 운동학적 일치 (no-penetration)**
$$\partial_z\phi + w^3 = v_b(t), \quad \forall x \in \mathcal{B},\ z=0$$
- **의미**: 유체가 물체를 뚫고 지나갈 수는 없으므로, 물체와 맞닿은 유체의 수직 속도는 물체의 수직 속도 $v_b$와 같아야 함.
- (1)과 시간미분으로 일관성: $\dot\eta = \dot z_b = v_b$.

**제약 (3): 수평 no-slip**
$$\partial_x\phi + w^1 = 0, \quad \forall x \in \mathcal{B},\ z=0$$
- **의미**: 물체는 강체이므로 바닥에서 유체가 수평으로 미끄러질 수 없음. 2.11a(응력 없음)는 공기-물 계면에서 맞지만, 물체와 접하는 곳에서는 성립 안함.

**발표 멘트**:
> "자유표면 조건 네 개 중에서, 운동학적 BC는 기하+운동학 제약 두 개로 쪼개지고, 응력 없음 조건은 수평 no-slip으로 대체됩니다. 동역학적 BC(2.15d)는 **여전히 모든 곳에서** 성립하지만, 물체 아래에서는 $p_s$가 미지수가 됩니다. 이 점이 v2의 핵심입니다."

---

## Slide 9 — Surface pressure $p_s$: the Lagrange multiplier

> **이 슬라이드는 v2의 가장 중요한 개념입니다. 다음 다섯 부분으로 차근차근 설명합니다.**
> (A) 라그랑주 승수의 기초 복습 (미적분)
> (B) 역학에서의 구속력 해석
> (C) PDE와 안장점 구조 (Stokes 방정식 비유)
> (D) 우리 문제에서 왜 $p_s$가 라그랑주 승수인가 (상세 유도)
> (E) 이 접근의 장점

---

### (A) 라그랑주 승수의 기초 복습 (미적분 수준)

**문제**: 제약조건 $g(x,y) = 0$ 하에서 $f(x,y)$를 최소화하라.

**간단한 예시**: $f(x,y) = x^2 + y^2$ 를 $g(x,y) = x + y - 1 = 0$ 제약 하에서 최소화.

**해법 (라그랑주)**: 보조함수
$$\mathcal{L}(x,y,\lambda) = f(x,y) - \lambda\, g(x,y) = x^2 + y^2 - \lambda(x + y - 1)$$
을 정의하고, 모든 변수에 대해 $\nabla\mathcal{L} = 0$:
$$\frac{\partial\mathcal{L}}{\partial x} = 2x - \lambda = 0, \quad
  \frac{\partial\mathcal{L}}{\partial y} = 2y - \lambda = 0, \quad
  \frac{\partial\mathcal{L}}{\partial \lambda} = -(x+y-1) = 0$$

풀면 $x = y = 1/2$, $\lambda = 1$. $\lambda$가 바로 **라그랑주 승수**.

**기하학적 해석**:
- $f$의 등고선과 제약면 $g=0$이 서로 **접할 때** 극값이 나옴.
- 접점에서 두 곡면의 법선 벡터가 평행: $\nabla f = \lambda \nabla g$.
- $\lambda$는 두 법선의 비율 = "제약을 조금 완화하면 목적함수가 얼마나 변하는가"의 비율.

**$\lambda$의 의미 (경제학적 해석, 섀도 프라이스)**:
> "제약을 $\epsilon$만큼 완화했을 때 (즉 $g = \epsilon$로 바꿨을 때), 최적값이 $\lambda\epsilon$만큼 변합니다. 이 '민감도'가 라그랑주 승수입니다."

---

### (B) 역학에서의 라그랑주 승수: 구속력 (constraint force)

**예시**: 질량 $m$인 구슬이 곡선 $y = f(x)$ 위를 미끄러짐 (중력만 작용, 마찰 없음).

**제약**: $g(x,y) = y - f(x) = 0$ (구슬이 곡선 위에 있어야 함).

**운동 방정식 (뉴턴)**:
$$m\ddot{\mathbf{r}} \;=\; \mathbf{F}_{\text{중력}} + \lambda \nabla g$$

여기서:
- 첫 항: 외력 (중력 $-mg\hat{y}$)
- 둘째 항: **구속력**. $\lambda$는 라그랑주 승수, $\nabla g$는 곡선에 수직 방향.

**물리적 의미**: $\lambda$는 곡선이 구슬을 밀어주는 **수직항력(normal force)의 크기**. 즉 구슬이 곡선을 뚫고 내려가지 않도록 "필요한 만큼" 반작용으로 밀어주는 힘.

**핵심 통찰**:
> "라그랑주 승수는 단지 수학적 편법이 아니라, 제약조건을 지키기 위해 실제로 존재해야 하는 **물리적 힘**을 대변합니다."

---

### (C) PDE에서의 라그랑주 승수: Stokes 방정식과 안장점 구조

이 개념은 유체역학에서 이미 익숙합니다. **비압축성 Stokes 방정식**을 봅시다:
$$-\nu\Delta\mathbf{u} + \nabla p = \mathbf{f}, \qquad \nabla\cdot\mathbf{u} = 0$$

**해석**:
- 우리는 "모멘텀 방정식"을 풀면서 동시에 "비압축성 제약" $\nabla\cdot\mathbf{u} = 0$을 만족시켜야 함.
- 이때 **압력 $p$는 비압축성 제약에 대응하는 라그랑주 승수**로 등장.
- 즉, $p$는 유체가 압축되지 않도록 "밀어주는 필요한 압력"입니다.

**행렬 형태로 쓰면 (유한차원 이산화 후)**:
$$
\begin{pmatrix} A & B^T \\ B & 0 \end{pmatrix}
\begin{pmatrix} \mathbf{u} \\ p \end{pmatrix}
= \begin{pmatrix} \mathbf{f} \\ 0 \end{pmatrix}
$$
이것이 바로 **안장점 문제(saddle-point problem)**. 위쪽 블록 $A\mathbf{u} + B^T p = \mathbf{f}$ 는 모멘텀식, 아래쪽 $B\mathbf{u} = 0$ 은 제약식, $p$가 라그랑주 승수.

**우리 v2 솔버의 행렬 구조도 정확히 이 패턴**을 따릅니다.

---

### (D) 우리 문제에서의 구체적 적용 (상세 유도)

#### 단계 1: 제약 식별

물체 아래 영역 $\mathcal{B}$에서 강제해야 할 **제약**:
$$C(x,t) \;=\; \partial_z \phi + w^3 - v_b(t) \;=\; 0, \quad \forall x \in \mathcal{B}$$

즉 "유체의 수직 속도가 물체 속도와 같아야 함" (운동학적 일치, no-penetration).

#### 단계 2: 동역학 BC를 보조 방정식(Lagrangian)처럼 보기

자유표면에서 동역학 BC (2.15d)는:
$$\partial_t \phi_s = -\frac{\eta}{\text{Fr}} + \frac{\kappa[\eta]}{\text{We}} + \frac{2}{\text{Re}}\Delta_H\phi_s - \frac{2}{\text{Re}}\partial_z w^3 - p_s$$

이 방정식이 **물체 아래에서도 유지**되는데, $p_s$ 항이 있습니다. 이 $p_s$를 라그랑주 승수로 간주합니다.

#### 단계 3: 제약과 승수의 "쌍(pair)"

- **자유 영역** ($x \notin \mathcal{B}$): 제약 없음 $\Rightarrow$ $p_s = 0$ (자유 $\phi_s$ 발전).
- **물체 아래** ($x \in \mathcal{B}$): 제약 $C = 0$ 있음 $\Rightarrow$ $p_s$가 미지수(승수)로 조정됨.

이것이 전형적인 **라그랑주 쌍**:
$$(\text{제약 }C = 0) \quad \longleftrightarrow \quad (\text{승수 }p_s \text{ 미지수})$$

#### 단계 4: 왜 $p_s$가 정확히 "피드백"을 담당하나?

- 동역학 BC에서 $p_s$ 항은 $\phi_s$의 시간 발전에 직접 영향 → $\phi$ (bulk)을 통해 $\partial_z\phi|_{z=0}$도 영향받음.
- $p_s$를 조정하면 $\partial_z\phi|_{z=0}$이 바뀌고, 결국 $C = \partial_z\phi + w^3 - v_b$ 값이 바뀜.
- 시스템이 $p_s$를 **정확히 $C = 0$이 되도록 스스로 결정**합니다.

**수식으로 말하면**:
$$p_s \;=\; \text{"제약을 지키기 위해 fluid가 body bottom에 가해야 하는 수직 응력"}$$

이는 역학 예시 (B)에서 "구슬이 곡선을 뚫고 내려가지 않도록 곡선이 미는 수직항력"과 **완전히 동일한 논리**.

#### 단계 5: 행렬 구조로 다시 보기

v2 monolithic matrix의 핵심 블록을 봅시다:

| 행 | $\phi_s$ 열 | $p_s$ 열 | 기타 |
|---|---|---|---|
| **$\phi_s$ 동역학 BC** (모든 $j$) | 대각 블록 | $+\Delta t$ | $w^3, \eta^n$ 결합 |
| **$p_s$ 제약식** (body $j$) | $+1/\Delta z$, $-1/\Delta z$ | $0$ | $w^3, v_b$ 결합 |
| **$p_s$ pin** (free $j$) | $0$ | $+1$ | — |
| **Newton의 $v_b$ 행** | $0$ | $-\Delta t \Delta x$ (body $j$만) | $M$ |

**읽는 법**:
1. 동역학 BC 행은 모든 표면 $j$에서 $\phi_s - \Delta t\, p_s = \text{RHS}$ 형태. $p_s$가 "힘"처럼 들어감.
2. $p_s$ 행은 물체 아래서는 $C = 0$ 제약식을 직접 강제.
3. $p_s$ 행은 자유 영역에서는 $p_s = 0$을 강제 (pin).
4. Newton의 $v_b$ 행은 $p_s$를 물체 위에서 적분해 힘으로 사용.

이 구조가 **Stokes 방정식의 saddle-point 구조와 정확히 일치**하며, 일반적 $A\mathbf{x} + B^T\boldsymbol{\lambda} = \mathbf{f}$, $B\mathbf{x} = \mathbf{g}$ 형식의 확장판입니다.

---

### (E) 이 접근의 장점

**1. 제약의 기계 정밀도 강제**:
- v1처럼 $\eta \leftarrow z_b$로 덮어쓰기 하면 제약은 지켜지지만 **압력 정보가 소실**됨.
- v2에서는 $p_s$가 unknown이 되어 제약이 모든 행렬 방정식과 **함께** 풀림 → 제약 오차는 SPARSE SOLVER의 수치 정밀도 수준.

**2. 물체에 가해지는 힘을 직접 측정**:
- 시뮬레이션 결과로 $p_s(x,t)$ 분포가 바로 출력됨 → $\int p_s\,dx$ = 물체에 작용하는 총 힘.
- 이 힘으로 stress distribution, impact force, resonance 등 분석 가능.

**3. Added-mass 불안정 회피 (핵심 수치적 장점)**:
- Staggered 방식 (fluid를 먼저 풀고 → body를 나중에 풂)에서는 **가벼운 물체**일수록 수치 불안정.
- 이유: body 가속도가 fluid 압력에 의존하는데, 반복이 엇갈려 잘못된 타이밍으로 업데이트됨.
- Monolithic implicit 방식에서는 $v_b$와 $p_s$가 **동시에** 결정 → 물리적으로 일관.

**4. 확장 용이**:
- 곡면 물체: 제약 $C$만 바꿔서 같은 구조 유지.
- 다중 물체: Lagrange 쌍을 여러 개 추가.
- 탄성 막: 다른 BVP지만 라그랑주 쌍 구조 그대로.

---

### 발표 시 멘트 (요약)

> "라그랑주 승수는 세 단계로 이해하면 쉽습니다.
> **첫째**, 미적분에서 제약 있는 최적화의 '섀도 프라이스'.
> **둘째**, 역학에서 제약(예: 막대 위의 구슬)을 지키기 위한 **수직항력**.
> **셋째**, 유체 PDE에서 비압축성 제약을 지키기 위한 **압력** (Stokes).
>
> 이 세 개념이 전부 같은 수학적 구조 — **제약 + 승수 쌍** — 을 가집니다.
>
> v2에서 $p_s$는 '유체가 물체 바닥에 가해야 하는 압력'으로, 운동학적 일치 $\phi_z + w^3 = v_b$ 를 강제하는 역할을 합니다. 자유 영역에서는 $p_s = 0$으로 pinned, 물체 아래에서는 unknown으로 남아 자동 결정. 이 구조 덕분에 implicit 결합이 실용적으로 돌아갑니다."

---

### 비유 (청중이 물리/수학에 익숙하지 않을 때)

> "자동차가 평평한 도로를 달릴 때, 땅이 위로 바퀴를 밀어주는 힘 (수직항력)이 있잖아요. 그 크기는 자동차 무게와 운동에 따라 자동으로 결정됩니다.
>
> v2에서 물체가 물 위에 떠 있을 때, 물이 물체 바닥을 위로 밀어주는 '수직 압력' $p_s$도 마찬가지로 자동 결정돼요. 우리는 '물이 물체를 뚫지 않는다'는 조건만 강제하고, 그 조건을 만족시키기 위해 필요한 압력은 수학이 알아서 찾아줍니다. 그게 바로 라그랑주 승수의 정확한 역할입니다."

---

## Slide 10 — Body dynamics: Newton's 2nd law

**운동 방정식 (무차원)**:

1. 운동학: $\dfrac{dz_b}{dt} = v_b$
2. 뉴턴 제2법칙: $M\dfrac{dv_b}{dt} = -\dfrac{M}{\text{Fr}} + \int_{\mathcal{B}} p_s\,dx$

**각 항 해석**:
- $M = m/(\rho L_c^2)$: 무차원 질량 (2D → per unit length)
- $-M/\text{Fr}$: 중력(무게). Fr이 작을수록 중력이 강함.
- $+\int_\mathcal{B} p_s\,dx$: 물체 아래에서 유체가 위로 밀어올리는 총 힘.
- **평형**에서: $\int_\mathcal{B} p_s\,dx = M/\text{Fr}$ ⇔ 유체 압력이 물체 무게를 정확히 지탱.

**교육 포인트 (대학 물리)**:
> "고등학교/대학 물리에서 배운 뉴턴 제2법칙 $F = ma$ 그대로입니다. $F$는 중력과 유체 반작용의 합이고, $a = dv_b/dt$. 특이한 점은 '유체 반작용'이 미리 알려지지 않고, 유체 방정식과 함께 풀어야 한다는 것."

---

## Slide 11 — Non-dimensionalisation (무차원화)

**왜 무차원화?**:
- SI 단위 대신 "특성 길이 $L_c$, 특성 속도 $U$" 기준의 비율로 모든 양을 나타내면, **물리 문제가 무차원 수 몇 개(Fr, We, Re)로 완전히 특징지어짐**.
- 같은 Fr/We/Re를 가진 실험은 (크기가 달라도) 같은 유체역학적 거동을 보임 (동역학적 유사성).

**네 가지 무차원 수**:

| 기호 | 정의 | 의미 |
|---|---|---|
| Fr | $U^2/(gL_c)$ | 관성 vs 중력 |
| We | $\rho U^2 L_c/\sigma$ | 관성 vs 표면장력 |
| Re | $UL_c/\nu$ | 관성 vs 점성 |
| $M$ | $m/(\rho L_c^2)$ | 물체 질량 (무차원) |

**코드에서 쓰는 기본값**: $L_c = 2.5$ cm, $T_c = 0.05$ s → $U = 50$ cm/s
- Fr = 1.02 → 관성과 중력이 거의 같은 정도
- We = 35.7 → 표면장력은 상대적으로 약함
- Re = 250,000 → 사실상 거의 점성 없는 상태(고 Re)

**발표 멘트**:
> "이 값들이 실제 물리적으로 무엇을 의미하는지 감을 잡아둡시다. Re가 25만이니까 점성은 매우 약하고, 거의 포텐셜 흐름에 가깝습니다. Fr이 1 근처니까 중력과 관성이 비슷한 스케일이고, We가 큰 편이니 표면장력은 부차적입니다."

---

## Slide 12 — Spatial discretisation (공간 이산화)

**격자 설정**:
- 균일 직사각형 격자: $n_x \times n_z$ 노드
- $x$방향 주기적, $z$방향 유한 (바닥부터 표면까지)
- $i=0$이 바닥 ($z=-D$), $i=n_z-1$이 표면 ($z=0$)

**유한차분 연산자들**:
- 수평 라플라시안: $(\Delta_H u)_j = (u_{j+1} - 2u_j + u_{j-1})/\Delta x^2$
  - 2차 테일러 전개에서 유도 — 3-점 stencil
- 1차 미분: $(D_x u)_j = (u_{j+1} - u_{j-1})/(2\Delta x)$ — centered
- 2D 라플라시안: 5-점 stencil (중심 + 상하좌우)

**body mask**:
```python
body_mask = np.abs(x - xc_obj) < a_obj
```
- **Boolean 배열** $\chi_j$: 컬럼 $j$가 물체 아래에 있는지 여부 (1 또는 0).
- 물체 기하가 바뀌지 않는다면 시간독립.

**교육 포인트**:
> "PDE를 수치로 풀 때는 연속적인 미분을 '근처 값들의 차이'로 근사합니다. 2차 테일러 전개에서 $u(x+\Delta x) - 2u(x) + u(x-\Delta x) \approx u''(x)\cdot\Delta x^2$. 이걸 정렬해 $u''$의 근사식이 됩니다."

---

## Slide 13 — The augmented state vector (확장된 상태 벡터)

**v1의 상태 벡터**: $[\vec{w}^3, \vec{w}^1, \vec{\phi}_s, \vec{\phi}_{\text{bulk}}, \vec{\eta}]$

**v2의 상태 벡터** (3개 추가):
$$\mathbf{x} = \begin{bmatrix} \vec{w}^3 \\ \vec{w}^1 \\ \vec{\phi} \\ \vec{\eta} \\ \vec{p}_s \\ z_b \\ v_b \end{bmatrix}, \quad N_{\text{total}} = 3 n_x n_z + 2 n_x + 2$$

**새로 추가된 것들**:
- $\vec{p}_s$ (크기 $n_x$): 표면 압력 — 물체 밖에서는 0이고, 물체 아래에서는 미지수
- $z_b, v_b$ (크기 1씩): 물체 상태 스칼라

**크기 추정**: $n_x=400, n_z=80$이면 $N_{\text{total}} \approx 96{,}800$.

**교육 포인트 (선형대수)**:
> "연립 선형방정식을 푼다는 건 $A\mathbf{x} = \mathbf{b}$ 를 푸는 것. 여기서 $A$는 약 97,000 × 97,000 행렬이지만, 대부분이 0(sparse)이라서 실제로는 다룰 수 있습니다. 보통 행당 5~6개 정도만 0이 아닌 값을 가져요."

---

## Slide 14 — Monolithic matrix: block structure (블록 구조)

**"monolithic" 의미**: 유체·표면·물체 방정식을 **한 번에**(쪼개지 않고) 하나의 큰 행렬로 합쳐서 푸는 방식.

**블록 구조 읽는 법**:
- 행 = 어느 방정식, 열 = 어느 미지수
- 대각 블록 (파란 ✓): "각 필드가 자기 자신에 대해 암묵적으로(implicit) 발전"
- 비대각 블록 (빨간 $\star$): "다른 필드와의 결합"

**주요 결합들**:
- $w^3 \leftrightarrow \phi$: 바닥 non-slip 조건
- $w^1 \leftrightarrow \phi, w^3$: 바닥 non-slip + 표면 BC
- $\phi \leftrightarrow w^3, p_s$: 표면 동역학 BC
- $\eta \leftrightarrow \phi, w^3, z_b$: 자유표면/물체 제약
- $p_s \leftrightarrow \phi, w^3, v_b$: 자유(0) / 물체(kinematic match)
- $z_b \leftrightarrow v_b$: 운동학
- $v_b \leftrightarrow p_s$: 뉴턴 법칙의 압력적분

**핵심 속성 (빨간 박스)**:
> "모든 행렬 계수가 격자·$\Delta t$·무차원 수·body mask 만으로 결정됨. **해(solution)에 의존하지 않음** → 한 번만 LU 분해하고 영원히 재사용."

---

## Slide 15 — $w^3$ block: three row types (세 종류 행)

**$w^3$ 방정식은 위치에 따라 세 가지 형태**:

**(1) 내부** ($1 \le i \le n_z - 2$): 열방정식 암묵적 이산화
$$\left(I - \frac{\Delta t}{\text{Re}}\Delta_{2D}\right) w^{3,n+1} = w^{3,n}$$
- 왼쪽: $n+1$ 시점 값 (implicit)
- 오른쪽: $n$ 시점 값 (known)

**(2) 바닥** ($i=0$): 비미끄럼 (implicit, $\phi$와 결합)
$$w^{3,n+1}_{0,j} + \frac{\phi^{n+1}_{1,j} - \phi^{n+1}_{0,j}}{\Delta z} = 0$$
- 중심차분으로 $\partial_z\phi$를 근사, $w^3 = -\partial_z\phi$ 즉 $u_z = \phi_z + w^3 = 0$ (바닥에서 수직속도 0).

**(3) 표면** ($i=n_z-1$): (2.15e)에서 유도
$$\left(I - \frac{2\Delta t}{\text{Re}}\Delta_H\right) w^{3,n+1}_s = w^{3,n}_s + \frac{2\Delta t}{\text{Re}}\Delta_H\left(\frac{\phi^n[-1] - \phi^n[-2]}{\Delta z}\right)$$
- RHS는 $\phi^n$ (explicit 사용) — LHS에 여분의 off-diagonal 결합을 피하기 위함.

**교육 포인트**:
> "`implicit vs explicit` 의 차이를 짚자면: explicit는 새 값을 '이미 아는 값들'로 직접 계산 (간단하지만 $\Delta t$ 제약 심함). Implicit는 '새 값들이 서로 엮인 연립방정식'을 풀어야 함 (한 스텝이 비싸지만 안정성 좋음). 우리는 공간도함수가 들어간 점성 항에 implicit를 쓰므로 $\Delta t$ 제약이 훨씬 느슨해집니다."

---

## Slide 16 — $w^1$ block: piecewise surface BC

**내부·바닥 행**: $w^3$과 동일한 구조 (열 방정식 + 바닥 non-slip).

**표면 행은 위치($\chi_j$)에 따라 두 가지로 분기**:

**자유 표면 ($\chi_j = 0$, 응력 없음 2.11a)**:
$$w^{1,n+1}_{-1,j} - w^{1,n+1}_{-2,j} + \Delta z (D_x w^{3,n+1}_{-1})_j = -\Delta z (\partial_{xz}\phi^n)_j$$
- $w^1$ 표면과 $w^3$ 표면이 (implicit) 결합, $\phi^n$은 explicit(RHS).

**물체 아래 ($\chi_j = 1$, 수평 no-slip)**:
$$w^{1,n+1}_{-1,j} + (D_x \phi^{n+1}_{-1})_j = 0$$
- $w^1$ 표면과 $\phi$ 표면이 (implicit) 결합.

**핵심 포인트**:
> "행렬에서 **한 줄**에 불과한데, `body_mask[j]` 값에 따라 다른 행렬 계수로 채워집니다. 코드로 보면 `if body_mask[j]: ... else: ...` 구조. 이 분기가 있어도 행렬 자체는 여전히 상수고 LU 재사용 가능."

---

## Slide 17 — $\phi$ block: Laplace + dynamic BC

**$\phi$ 블록은 세 종류의 행**:

**바닥 ($i=0$)**: 고스트 포인트로 no-slip 포함한 라플라스
$$(-2\alpha - 2)\phi_{0,j} + \alpha(\phi_{0,j\pm 1}) + 2\phi_{1,j} - 2\Delta z\, w^3_{0,j} = 0$$
- $\alpha = (\Delta z/\Delta x)^2$
- 고스트 포인트 trick: 바닥 바깥 점을 가상으로 두고, non-slip ($\partial_z\phi|_{z=-D} = -w^3_0$)로부터 표현.

**내부** ($1 \le i \le n_z-2$): 5점 스텐실 라플라시안
$$(-2\alpha-2)\phi_{i,j} + \alpha\phi_{i,j\pm 1} + \phi_{i\pm 1,j} = 0$$

**표면** ($i=n_z-1$): 동역학적 BC (2.15d), **모든 노드에서 동일**
$$\left(I - \frac{2\Delta t}{\text{Re}}\Delta_H\right)\phi^{n+1}_{-1,j} + \frac{2\Delta t}{\text{Re}}\frac{w^3_{-1,j}-w^3_{-2,j}}{\Delta z} + \Delta t\, p^{n+1}_{s,j} = \phi^n_{-1,j} - \frac{\Delta t}{\text{Fr}}\eta^n_j + \frac{\Delta t}{\text{We}}\kappa[\eta^n]_j$$

**핵심 아이디어**:
> "$\phi$ 표면 행은 모든 $j$에서 **같은 형식** — $p_s$ 항이 항상 들어있어요. 그럼 자유 영역에서 $p_s$를 어떻게 0으로 만드나? → 별도의 $p_s$ 행이 '자유에선 $p_s = 0$, 물체 아래서는 운동학 제약'을 담당합니다."

---

## Slide 18 — $\eta$ and $p_s$ blocks: the constraint switch

**두 블록 모두 '자유/물체' 로 분기**:

**$\eta$ 행**:
- 자유 ($\chi_j = 0$): 운동학 BC 이산화 (2.15c)
  $$\eta^{n+1}_j - \frac{\Delta t}{\Delta z}(\phi^{n+1}_{-1,j} - \phi^{n+1}_{-2,j}) - \Delta t\, w^{3,n+1}_{-1,j} = \eta^n_j$$
- 물체 ($\chi_j = 1$): 기하 제약
  $$\eta^{n+1}_j - z^{n+1}_b = 0$$

**$p_s$ 행**:
- 자유 ($\chi_j = 0$): 압력 0으로 고정
  $$p^{n+1}_{s,j} = 0$$
- 물체 ($\chi_j = 1$): 운동학 일치 (Lagrange 승수 역할)
  $$\frac{\phi^{n+1}_{-1,j} - \phi^{n+1}_{-2,j}}{\Delta z} + w^{3,n+1}_{-1,j} - v^{n+1}_b = 0$$

**발표 멘트**:
> "여기가 v2의 진짜 아름다운 지점입니다. 물체 아래서 $\eta$ 행이 'η = z_b' 제약으로 바뀌고, $p_s$ 행은 '운동학 일치'로 바뀝니다. 이 두 변화가 $p_s$를 라그랑주 승수로 자동 조정시키는 구조를 만듭니다. 자유 영역에서는 $p_s$를 그냥 '0으로 pin'해놓고."

**↪ Slide 9와 연결**:
> "이 행렬 구조가 정확히 Slide 9의 라그랑주 쌍(제약 + 승수)을 구현한 것입니다.
> - **$p_s$ 행, body node = 제약 방정식** (유체 수직속도 = 물체 속도)
> - **$\phi_s$ 동역학 BC의 $p_s$ 열 = 승수 결합** (제약을 지키기 위해 필요한 압력이 $\phi_s$ 발전에 피드백)
> - **$p_s$ 행, free node = 승수 자유도 pinning** ($p_s = 0$ 이므로 자유 영역에서는 승수 없음)
>
> Stokes 방정식의 saddle-point 구조와 동일한 **$(A, B^T; B, 0)$** 패턴이 여기서 반복됩니다."

---

## Slide 19 — Body rows: discretised Newton's 2nd law (물체 방정식)

**$z_b$ 방정식** (Implicit Euler of $\dot z_b = v_b$):
$$z^{n+1}_b - \Delta t\, v^{n+1}_b = z^n_b$$

**$v_b$ 방정식** (뉴턴 제2법칙):
$$M v^{n+1}_b - \Delta t \Delta x \sum_{j\in\mathcal{B}} p^{n+1}_{s,j} = M v^n_b - \frac{\Delta t M}{\text{Fr}}$$

**유도**:
1. 시작: $M\dot v_b = -M/\text{Fr} + \int_\mathcal{B} p_s\,dx$
2. Implicit Euler: $M(v^{n+1}_b - v^n_b)/\Delta t = -M/\text{Fr} + \Delta x \sum_j p^{n+1}_{s,j}$
3. $\Delta t$ 곱하고 미지수는 좌변으로, 아는 값은 우변으로.

**midpoint 합** ($\Delta x \sum_{j\in\mathcal{B}} p_{s,j}$): 적분 $\int_\mathcal{B} p_s\,dx$의 이산 근사 (사각형 근사).

**포인트**:
> "이 두 행만 추가하면 물체 역학이 유체 방정식과 완전히 결합됩니다. 중요한 건 모든 계수가 여전히 상수라는 것 — $M$, $\Delta t$, $\Delta x$, $\text{Fr}$ 전부 시간에 변하지 않죠."

---

## Slide 20 — Why it all works: constant coefficient matrix

**결정적 관찰**:
모든 행렬 항이 다음 중 하나에만 의존:
- 격자 파라미터 ($n_x, n_z, \Delta x, \Delta z$)
- 시간 스텝 $\Delta t$
- 무차원 수 (Fr, We, Re, M)
- body mask $\chi$

**해(solution)에 의존하지 않음!** → 행렬 $A$는 시간에 상수.

**결과: Sparse LU 분해 한 번**:
$$A = LU \quad (\text{setup에서 }\textbf{한 번만}\text{ 계산})$$

매 시간 스텝:
$$\mathbf{x}^{n+1} = U^{-1}(L^{-1}\mathbf{b}(\mathbf{x}^n)), \quad \text{cost: } \mathcal{O}(N) \text{ back-substitution}$$

**LU 분해란? (선형대수 복습)**:
> "$Ax = b$를 풀 때, $A = LU$ (하삼각 × 상삼각)로 분해해놓으면, $Ly = b$와 $Ux = y$를 차례로 풀 수 있어요. 두 번의 삼각계 풀이는 각각 $O(N)$이어서 총 비용이 선형입니다. Dense 행렬이면 분해 자체가 $O(N^3)$이지만, sparse 행렬이면 훨씬 싸고 — 더 중요한 건 **한 번만 분해하면 됨**."

**성능 비교** (v1에서 관찰):
- v1 sequential solver: 34 ms/step
- v1 monolithic: 1.9 ms/step (~18× 빠름)
- v2 monolithic: 유사한 개념, 약간 큰 시스템이지만 여전히 ms/step 단위.

---

## Slide 21 — Algorithm: one full time step (알고리즘)

**Setup (시뮬 시작 시 1회만)**:
```python
A, offsets = build_monolithic_matrix_v2(nx, nz, dx, dz, dt,
                                        Fr, We, Re, M, body_mask)
LU = splu(A)   # O(N^(3/2)) 정도, 한 번만
```

**매 시간 스텝 (back-substitution만)**:
```python
def _step(w3, w1, phi, eta, ps, zb, vb):
    # 1. 시간 n의 상태로 RHS 구성 (explicit 항들)
    b = build_monolithic_rhs_v2(w3, w1, phi, eta, zb, vb, ...)

    # 2. A x = b 풀기 (back-sub, 저렴)
    x = LU.solve(b)    # O(N)

    # 3. 시간 n+1의 상태로 unpack
    return unpack_state(x, offsets)
```

**흐름 요약**:
1. 현재 상태 → RHS 조립
2. LU 사용해 새 상태 벡터 1회 solve
3. 벡터를 (w³, w¹, φ, η, p_s, z_b, v_b)로 분해

**발표 멘트**:
> "실제 코드가 이 세 줄로 요약됩니다. 여기서 LU.solve가 핵심이고, 이건 몇 밀리초 수준입니다. 오히려 RHS 조립이 더 비쌀 수도 있어요."

---

## Slide 22 — Configuration dataclasses (설정 데이터클래스)

**Python dataclass 사용**: 파라미터를 깔끔하게 조직.

**PhysicalParams**: 차원 있는 물리량 (cm, g, s).
- `L_cm`, `D_cm`: 영역 크기
- `g`, `sigma`, `rho`, `nu`: 물리 상수 (중력, 표면장력, 밀도, 점성)

**BodyParams**: 물체 파라미터.
- `half_width_cm`, `center_x_cm`: 기하
- `mass_per_length_g_per_cm`: 단위 길이당 질량 (무차원 $M$으로 변환됨)
- `initial_z_cm`, `initial_vz_cm_s`: 초기 조건

**NumericalParams**: 수치 파라미터.
- `nx`, `nz`: 격자
- `nt`, `total_time_s`: 시간 스텝 수, 총 시간
- `unit_length`, `unit_time`: 무차원화 스케일

**발표 멘트**:
> "사용자는 이 세 dataclass만 건드리면 시뮬을 구성할 수 있어요. 실험 파라미터 바꾸고 싶으면 `BodyParams(mass_per_length_g_per_cm=5.0)` 이렇게만 쓰면 됩니다."

---

## Slide 23 — Matrix builder: body dynamics rows (행렬 조립, 물체 행)

**Row block 6: $z_b$ 방정식**:
```python
A[zb, zb] = 1.0    # z_b 자신
A[zb, vb] = -dt    # - Δt · v_b
```

**Row block 7: $v_b$ 방정식**:
```python
A[vb, vb] = M                      # inertia
for j in range(nx):
    if body_mask[j]:
        A[vb, ps(j)] = -dt * dx    # 압력 적분 결합
```

**Row block 5: $p_s$ 행 (라그랑주 승수)**:
```python
for j in range(nx):
    row = ps(j)
    if body_mask[j]:
        # 운동학 일치: φ_z + w³ - v_b = 0
        A[row, phi(nz-1, j)] =  1.0/dz
        A[row, phi(nz-2, j)] = -1.0/dz
        A[row, w3(nz-1, j)]  =  1.0
        A[row, vb]           = -1.0
    else:
        A[row, ps(j)] = 1.0    # 자유 영역은 0으로 pin
```

**발표 멘트**:
> "이게 실제 코드입니다. 수식과 직접 1:1 대응하죠. `A[vb, ps(j)] = -dt*dx`는 뉴턴 법칙의 $-\Delta t \Delta x \sum p_s$ 항을 정확히 표현. `body_mask[j]`가 분기 지점입니다."

---

## Slide 24 — Solver class: FloatingObjectSolverV2

**클래스 구조**:

```python
class FloatingObjectSolverV2:
    def __init__(self, config):
        self.cfg = config
        self._setup()       # 격자, 마스크, 연산자

        # A를 한 번만 빌드, 한 번만 factor
        A, self.off = build_monolithic_matrix_v2(...)
        self.LU = splu(A)

    def _step(self, w3, w1, phi, eta, ps, zb, vb):
        b = build_monolithic_rhs_v2(
            w3, w1, phi, eta, zb, vb,
            self.off, self.dx, self.dz, self.dt,
            self.Fr, self.We, self.Re, self.M,
            self.body_mask, self.Delta_H, self.Dx,
        )
        return unpack_state(self.LU.solve(b), self.off)

    def run(self):
        # 초기화 → nt 스텝 루프 → history dict 반환
        ...
```

**사용 예**:
```python
solver = FloatingObjectSolverV2(config)
result = solver.run()
# result["z_b"], result["v_b"], result["eta_history"], ...
```

**발표 멘트**:
> "객체 하나 만들고 `run()` 호출하면 끝. 내부에서 한 번만 LU 분해하고, 매 스텝마다 back-substitution. 깔끔한 API입니다."

---

## Slide 25 — Theoretical expectation: hydrostatic equilibrium

**정지 상태에서의 힘 균형** (뉴턴 법칙 → 0 = ...):
$$0 = -\frac{M}{\text{Fr}} + \int_\mathcal{B} p_s\,dx \implies \int_\mathcal{B} p_s\,dx = \frac{M}{\text{Fr}}$$

**$p_s$가 물체 아래 대체로 균일하다고 가정**:
- 동역학 BC (2.15d)에서 정지 상태 + 평탄 표면 조건 → $p_s \approx -\eta/\text{Fr}$
- 물체 아래 $\eta = z_b$이므로 $p_s \approx -z_b/\text{Fr}$
- 적분: $p_s \cdot 2a = -z_b \cdot 2a / \text{Fr} = M/\text{Fr}$
- 정리: $\boxed{|z_b^{\text{eq}}| = \dfrac{M}{2a} = \dfrac{m}{\rho \cdot 2a}}$

**아르키메데스 원리의 재발견!**:
- 유체에 잠긴 부피당 물리적 무게(= 밀도 × 부피 × 중력) = 물체 무게
- 2D에서 "잠긴 부피" = (폭 $2a$) × (잠긴 깊이 $|z_b|$)

**예시**:
- $m = 0.5$ g/cm, $2a = 10$ cm, 물 ($\rho = 1$ g/cm³)
- $|z_b^{\text{eq}}| = 0.5/10 = 0.05$ cm $=$ **0.5 mm**

**발표 멘트**:
> "이 예측을 스모크 테스트로 확인할 겁니다. 0.5mm는 매우 작아서, 시뮬에서 빠르게 수렴하는 것을 볼 수 있어야 합니다."

---

## Slide 26 — Smoke test: release from rest (스모크 테스트)

**파라미터** (빠른 iteration용 작은 격자):
- $n_x \times n_z = 80 \times 20$, $\Delta t = 5$ ms, 30 스텝 ($t_{\text{end}} = 0.15$ s)
- $m = 0.5$ g/cm, $2a = 10$ cm, Fr = 1.02, We = 35.7, Re = $2.5\times 10^5$
- 초기 조건: $z_b(0) = 0$, $v_b(0) = 0$ (정지 상태에서 놓음)

**솔버 출력** (실제 결과):
```
FloatingObjectSolverV2 initialised
  Grid: 80x20, Body: 2a=10.0 cm, M_dim=0.0800
  Body nodes: 15/80
  Monolithic size N = 4962, nnz = 23938
  Matrix build: 0.02s,  LU factor: 0.01s
  ...
  Step  30 | t=0.150s | z_b=-0.593mm | v_b=-0.317cm/s | F=+7.88e-02
```

**해석**:
- $N = 4962$개 미지수, 그중 약 24k개 nonzero → 행당 평균 5개
- LU 분해 0.01s — 매우 빠름
- 0.15초 후 $z_b = -0.593$ mm까지 가라앉음 (평형 $\sim 0.5$ mm에 근접)
- 여전히 $v_b < 0$이므로 약간 더 깊이 가라앉다가 반발할 예정

**발표 멘트**:
> "이론 예측이 $|z_b^{\text{eq}}| = 0.5$ mm 였고, 0.15초 만에 $z_b = -0.59$ mm 에 도달했으니 아직 평형을 약간 지나친 상태입니다. 앞으로 더 진행되면 진동하면서 평형으로 수렴할 것으로 예상."

---

## Slide 27 — Smoke test: sanity checks (검증)

**체크할 네 가지**:

| 양 | 예상 | 관찰 ($t=0.15$ s) |
|---|---|---|
| $p_s$ (물체 바깥) | 기계 정밀도로 0 | $\max|p_s^{\text{free}}| = 0$ ✓ |
| $p_s$ (물체 아래) | 가장자리에 집중, 양수 | $[0.033, 0.021, \ldots, 0.018, \ldots, 0.033]$ ✓ |
| $\int_\mathcal{B} p_s\,dx$ | $M/\text{Fr} = 0.0784$ | $0.0788$ (+0.5%) ✓ |
| 물체 운동 | 중력으로 가라앉고 감속 | $z_b: 0 \to -0.59$ mm, $v_b$ 방향 전환 중 |
| 성능 | 한 번 factor, 빠른 풀이 | LU 0.01 s, 스텝당 0.37 ms |

**판정 (Verdict 박스)**:
> "네 가지 정성적 서명 모두 관찰: 제약조건 정확히 강제됨, 힘이 아르키메데스 예측으로 수렴, 물체 역학이 물리적, 계산 비용 무시할 만함."

**+0.5% 편차 해설**:
> "0.0788 vs 0.0784의 차이는 $M\dot v_b$ 가속도 항 때문입니다. 즉 물체가 아직 감속 중이어서 완전한 평형은 아니라는 의미. 더 시간이 지나면 편차가 줄어들 것."

**$p_s$ 프로파일 ([0.033, 0.021, ..., 0.018, ..., 0.033])**:
- 가장자리 값(0.033)이 중앙(0.018)보다 큼
- 물리적 이유: 가장자리에서 곡률 점프 + 유체 흐름이 집중되어 압력이 높음
- 정확한 Hertzian 접촉과 비교할 수 있는 점

---

## Slide 28 — Another scenario: drop impact

**파라미터**: $m = 5$ g/cm (10배 더 무거움), $z_b(0) = 0$, $v_b(0) = -20$ cm/s (하강 중 충격).

**예상 동역학**:
1. 물체가 표면 아래로 가라앉음 (중력으로 가속, 유체 반작용으로 감속)
2. 대략 7 mm 깊이까지 내려감 (운동에너지 → 위치에너지 대략 추정)
3. 반발: 물체 상승, 평형을 지나침, 진동 시작
4. 진동 감쇠: 파 방출 + 점성 소산
5. 최종 평형: $m/(\rho \cdot 2a) = 5/10 = 5$ mm

**진동 주기 (added-mass 추정)**:
$$T \approx 2\pi\sqrt{\frac{M_{\text{eff}}}{2a \cdot (\rho g / L_c)}}, \quad M_{\text{eff}} = m + m_{\text{added}}$$

**Added mass 설명 (물리)**:
> "평평한 2D strip이 유체 표면에서 움직이면, 주변 유체도 함께 움직여야 하므로 **추가 질량** $m_{\text{added}} \sim \rho a^2$가 관성에 보태집니다. 무거운 물체는 천천히 진동하고, 가벼운 물체는 빨리 반응해요."

**발표 멘트**:
> "이 시나리오가 실험에서 흔한 '물체를 떨어뜨려 튀는 것을 관찰' 이랑 직결됩니다. 무게와 초기 속도를 바꿔가며 진동 주기, 감쇠율 등을 측정할 수 있습니다."

---

## Slide 29 — Caveat: periodic boundary ⇒ wave build-up

**문제 현상**:
- 물체가 방출한 파가 주기 경계로 감싸져 **다시 물체로 돌아옴**
- 물체 - 복귀파 상호작용이 인위적으로 발생

**증상**:
- $x=0/L$ 경계 근처 파 진폭이 시간에 따라 증가
- 물체가 완전한 정수압 평형에 도달하지 못함
- 에너지가 증가하는 것처럼 보임 (닫힌 영역에서 빠져나가지 못해서)

**진단법**: $\eta(x,t)$ space-time 다이어그램 보기 — 외향 파가 경계에서 반대로 재진입하면 wrap-around 모드.

**해결책**:
- **짧은 run**: $t_{\text{end}} < L/(2c_{\text{wave}})$ — 파가 돌아오기 전에 끝냄.
- **큰 영역**: $L \gg c_{\text{wave}} t_{\text{end}}$.
- **흡수층 (absorbing/sponge layer)**: 경계 근처에 감쇠 추가. 아직 미구현.

**발표 멘트**:
> "이건 솔버 버그가 아니라 주기 경계의 **물리적 특성** 때문입니다. 실제로 긴 수조를 시뮬하려면 흡수층을 넣거나 충분히 큰 $L$을 써야 해요. 기존에 유저가 2초에서 10초로 늘렸더니 이상해 보였던 GIF 현상이 바로 이거였습니다."

---

## Slide 30 — Limitations of v2 (한계)

**7가지 주요 한계**:

1. **선형화된 자유표면**: $|\eta| \ll D$ 가정. 강한 충격이나 매우 무거운 물체는 한계.
2. **Implicit Euler 시간 적분**: 1차 정확도. 안정성은 좋지만 고정밀 장시간 적분에는 부족.
3. **평평한 바닥 강체만**: 곡면/변형 가능/비강체 접촉 지원 안 함.
4. **2D 직교 좌표만**: axisymmetric, 3D 불가.
5. **주기 수평 BC**: 파 wrap-around, 열린 경계 방사 불가.
6. **접촉선 물리 없음**: 물체 가장자리의 급격한 불연속 — 실제 물리는 더 복잡.
7. **분리/재접촉 없음**: 접촉 후 평생 접촉 상태 유지.

**발표 멘트**:
> "이 한계들이 '이 모델이 틀렸다'는 뜻은 아닙니다. 적용 범위가 명확하다는 의미예요. 잔잔한 진동·작은 진폭·강체 접촉 시나리오에서는 아주 잘 작동합니다. 더 넓은 범위로 확장하려면 아래 슬라이드의 '확장' 방법들이 필요합니다."

---

## Slide 31 — Natural extensions (자연스러운 확장)

**시간 적분**:
- **VS-BDF2**: 가변 스텝 후진차분공식, 2차 정확도. 행렬은 여전히 상수, Euler 계수만 달라짐. [Hairer & Wanner 1996]

**도메인·기하**:
- **흡수/sponge 층**: 주기 경계 반사 억제
- **Axisymmetric 좌표**: 구 충돌 문제 (Galeano-Rios 2017, Agüero 2022)
- **다중 물체**: 라그랑주 승수 블록 여러 개 추가. 행렬은 여전히 sparse & 상수.

**물리**:
- **곡면 물체**: 기하 제약을 $\eta = z_b + s(x-x_c)$로 일반화 ($s$ = 물체 프로파일).
- **Faraday 강제**: 연직 진동 항 추가 → pilot-wave 유사체 (Milewski 2015).
- **탄성 막(elastic membrane)**: 강체 대신 막. 다른 BVP지만 제약-승수 구조는 동일.

**발표 멘트**:
> "이 확장들의 공통점은 **모두 v2의 기본 골조를 유지**한다는 것입니다. 새 물리를 넣어도 '선형, 상수 계수, monolithic, LU 한 번' 전략이 계속 살아있죠. 그게 이 접근의 강력함."

---

## Slide 32 — References

**주요 참고 논문**:

1. **Galeano-Rios, Milewski, Vanden-Broeck (2017)**. Non-wetting impact of a sphere onto a bath. *J. Fluid Mech.* 826, 97-127. [LNS 프레임워크, 출발점]
2. **Milewski, Galeano-Rios, Nachbin, Bush (2015)**. Faraday pilot-wave dynamics. *J. Fluid Mech.* 778, 361-388.
3. **Galeano-Rios (2016)**. Hydrodynamic Pilot-waves. PhD thesis, IMPA. [세부 유도]
4. **Agüero, Alventosa, Harris, Galeano-Rios (2022)**. Impact of a rigid sphere onto an elastic membrane. *Proc. R. Soc. A* 478, 20220340.
5. **Lamb (1932)**. *Hydrodynamics*, 6th ed., CUP. (고전)
6. **Batchelor (1967)**. *An Introduction to Fluid Dynamics*, CUP. (교과서)
7. **Hairer & Wanner (1996)**. *Solving ODEs II: Stiff and DAE Problems*, Springer. (BDF2 등)

**코드**: `03_Floating_Object/floating_object_v2.py`
**원고**: `99_overleaf_Floating_LNS_Wongyung/main.tex` (§2.5, §3.4)

**발표 멘트**:
> "이 프로젝트는 완전히 새로운 것이 아니라, Galeano-Rios 팀의 LNS 프레임워크를 '부유체 문제'로 확장한 작업입니다. 중요한 건 — 원 논문 저자들도 공개 코드를 내놓지 않았기 때문에, 여기 구현이 이 분야의 첫 공개 구현 중 하나라는 점이에요."

---

## Slide 33 — Summary & Questions

**세 박스로 요약**:

**What v2 does**:
> 강체가 점성 유체 위에서 떠 있음. 중력과 자기일관적 유체 압력으로 운동 결정 — 외부에서 지정 안 함.

**How it works**:
> 미지수에 $(p_s, z_b, v_b)$ 추가. 자유표면 BC를 세 가지 제약으로 대체. 표면 압력을 라그랑주 승수로 사용. 하나의 monolithic 행렬 조립, 한 번 분해, 매 스텝 back-solve.

**Why it works**:
> 모든 방정식이 선형이고 계수가 해에 독립. 스텝당 $\mathcal{O}(N)$ back-substitution으로 **fully implicit** coupling이 실용적임.

**마지막 멘트**:
> "세 문장으로 오늘 발표를 정리합니다. 무엇을 — 유체 위 물체의 자기일관 역학. 어떻게 — 확장된 monolithic 시스템. 왜 효율적인가 — 상수 계수 덕분에 LU 한 번만. 질문 받겠습니다."

---

## 발표 팁 (종합)

### 시간 배분 (45분 발표 기준)
- Part 1 (Context, 슬라이드 1-4): 5분
- Part 2 (Physical model, 5-11): 10분
- Part 3 (Numerical method, 12-21): 15분 — **여기가 가장 밀도 높음**
- Part 4 (Implementation, 22-24): 5분
- Part 5 (Examples, 25-29): 7분
- Part 6 (Outlook, 30-33): 3분

### 예상 질문과 대응

**Q1. "v1에서도 됐는데 왜 v2가 필요한가?"**
> A. v1은 물체를 **강제 진동시키는** 장치처럼 다뤄서 '진짜 물리'가 아니었음. 물체가 유체와 어떻게 상호작용하는지(떨어지고, 반발하고, 평형 찾는) 궁금하면 v2 필요.

**Q2. "라그랑주 승수 접근이 왜 필요한가? 그냥 $\eta = z_b$ 대입하면 되지 않나?"**
> A. 그게 v1의 '$\eta \leftarrow z_b$ 덮어쓰기' 방식. 문제는 압력 분포를 얻을 수 없다는 것. v2에서는 $p_s$가 자동으로 결정되어서 물체에 작용하는 힘을 알 수 있음.

**Q3. "Implicit Euler는 1차 정확도라 부정확하지 않나?"**
> A. 맞음. 하지만 안정성이 매우 중요한 coupled 문제에서는 1차라도 괜찮음. 정확도가 필요하면 VS-BDF2로 확장 가능 (행렬 형태 유지됨).

**Q4. "주기 경계가 문제라면서 왜 그걸 썼나?"**
> A. 수치적으로 가장 단순하기 때문. 흡수층을 넣는 것은 향후 과제. 단기 시뮬(파가 wrap around 하기 전)에서는 주기 경계도 문제없음.

**Q5. "왜 2D로 제한했나? 실제 물체는 3D 아닌가?"**
> A. 2D로 이미 핵심 물리 다 담김. 3D 확장은 사소하게는 아니지만(격자 크기 × 10-100배), 수식 구조는 동일. 우선 2D에서 잘 검증한 뒤 3D로.

---

## 용어 참고 (한-영)

| 한국어 | English |
|---|---|
| 자유표면 | free surface |
| 포텐셜 흐름 | potential flow |
| 와류 성분 | vortical component |
| 라플라시안 | Laplacian |
| 열 방정식 | heat equation |
| 운동학적 BC | kinematic BC |
| 동역학적 BC | dynamic BC |
| 비미끄럼 | non-slip |
| 응력 없음 | stress-free |
| 라그랑주 승수 | Lagrange multiplier |
| 제약조건 | constraint |
| 구속력 | constraint force |
| 수직항력 | normal force |
| 안장점 문제 | saddle-point problem |
| 섀도 프라이스 | shadow price (Lagrange multiplier의 경제학적 해석) |
| 비압축성 | incompressibility |
| KKT 조건 | Karush--Kuhn--Tucker conditions |
| 운동학적 일치 | kinematic match |
| 강제 (강제 방정식) | enforcement / constraint enforcement |
| 평형 | equilibrium |
| 아르키메데스 원리 | Archimedes' principle |
| 무차원화 | non-dimensionalisation |
| 유한차분 | finite difference |
| 암함수/양함수 | implicit / explicit |
| 희소행렬 | sparse matrix |
| LU 분해 | LU decomposition |
| Back-substitution | 후진 대입 |
| Reynolds 수 | Reynolds number |
| Froude 수 | Froude number |
| Weber 수 | Weber number |
| Added mass | 유체 부가질량 |
| 주기 경계 | periodic boundary |
| 흡수층 | absorbing / sponge layer |

---

*이 노트는 `slides_v2.pdf` (33 슬라이드)에 1:1 대응합니다. 발표 중 슬라이드 번호를 따라 참조하세요.*
