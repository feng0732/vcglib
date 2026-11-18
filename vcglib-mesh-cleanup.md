# 网格清理对实体拓扑与体积解释的影响

**仓库**：vcglib（分支 `task/29`，HEAD `5cae2ab`）
**分析对象**：`vcg/complex/algorithms/clean.h`、`vcg/complex/algorithms/inertia.h`、`vcg/complex/algorithms/hole.h`、`vcg/complex/algorithms/inside.h`、`vcg/complex/algorithms/update/topology.h`、`vcg/simplex/face/topology.h`、`vcg/space/triangle3.h`
**方法**：只使用仓库已有语义与已实现的函数。为取得可观察数字，在仓库**之外**（`/tmp`）用 `g++ 13.3 -std=c++17 -O1` 直接引用仓库头文件编写探针程序，构造 `tri::TriMesh` 并调用 `Clean`/`Hole`/`Inertia`。除特别标注 "NDEBUG" 的条目外，所有数字均来自带 `assert` 的默认构建（失败会 abort，因此表中未 abort 的行即"函数正常返回"）。探针源码与完整输出未进入仓库。

---

## 0. 结论速览

1. `tri::Inertia::Mass()` 在代码层面**不是**"实体体积"，而是边界曲面的 **x 方向通量** ∮ x·n_x dA。闭合流形上它才等于体积；开放网格上它等于"体积减去缺失片通量"，而且**可能凑巧等于正确答案**（§2、§6）。
2. 仓库把"清理"和"拓扑重建"**显式解耦**：`RemoveDuplicateVertex` / `RemoveDuplicateFace` 都声明不改拓扑（`clean.h:202-206`）。清理后不重建 `FF`，查询会读到过期拓扑（实测 `CountHoles` 仍是 12，重建后才变 0，§3.4）。
3. 四类缺陷对实体判断的作用点完全不同：
   * **重复顶点**：改变**连通性与边界环**（一个闭合盒子裂成 12 个片、12 条边界环），但**不改变** `Mass` 的数值；
   * **零面积面**：分"拓扑退化"（同一顶点指针出现两次）与"几何零面积"（共线）两种；前者会**打断拓扑构建**，后者只制造共面/非流形边并污染欧拉数；
   * **反向法线**：完全不改变 `IsWaterTight`（它只看边界/非流形边），只改变 `Mass` 的**符号和大小**，且影响量 ∝ 该面 `n_x·A`；`n_x=0` 的面翻转后 `Mass` 不变（易被误判为"无问题"）；
   * **边界孔洞**：`IsWaterTight` 正确地给出 `false`，但 `Mass` 变成"缺了一块的曲面通量"，可能恰好等于正确答案。
4. 仓库中最强的"封闭"判定 `IsWaterTight` 有明确的**盲区**：不看可定向性、不看非流形顶点（`clean.h:1191-1198`）；实测"两张三角形背对背粘成的零厚度片"、"两个仅在一点相接的实心锥体"、"两个完全重合的闭合盒子"都被判为 `WT=1`（§5.4、§7）。
5. 因此，**"清理成功"≠"实体可解释"**。§8 给出一套只依赖仓库函数的判定流程，以及必须标记为"不确定"的完整清单。

---

## 1. 仓库给出的网格语义（作为讨论边界）

| 概念 | 仓库实现 | 语义要点 |
|---|---|---|
| 面朝向 | 顶点绕序（`face::SwapEdge` 交换 V0/V1，`clean.h:1542-1551`） | 存储的 `N()` **不参与**体积/朝向判定 |
| 面法向 | `NormalizedTriangleNormal` = 叉积单位化（`triangle3.h:88-92`） | 由绕序即时重算 |
| 面积 | `DoubleArea` = 叉积模长（`triangle3.h:241-245`），即**两倍面积** | 命名上易混淆 |
| 边界 | `face::IsBorder` = `cFFp(j)==&f`（`face/topology.h:55-63`） | 依赖 FF 邻接 |
| 非流形边 | `face::IsManifold` = FF 对合检查（`face/topology.h:41-48`） | 依赖 FF 邻接 |
| 朝向一致性 | `CheckOrientation`：相邻面共享边两端顶点是否"反向"（`face/topology.h:380-394`） | 只在**分量内部**比较，无全局符号 |
| "封闭" | `IsWaterTight`：无边界边且无非流形边（`clean.h:1199-1204`） | 不需要 FF；不查可定向性 |
| 体积 | `Inertia::Mass()`（`inertia.h:233`）；`Stat::ComputeMeshVolume` 只是它的包装（`stat.h:285-289`） | 文档说"Meaningful only if the mesh is watertight"（`inertia.h:231`） |
| 内部判定 | `Inside::Is_Inside` 用**存储的** `f->N()` 与最近点方向（`inside.h:68`） | 与体积路径的语义来源不同 |

**关键顺序约束**（仓库自己写明）：`RemoveDuplicateVertex`"does not update any topology relation … usually performed BEFORE building any topology information"（`clean.h:202-206`）。同理 `RemoveDuplicateFace`（`clean.h:329-334` 注释）。

---

## 2. `Mass()` 的真实语义：x 方向通量，而非"体积积分"

### 2.1 代码事实

`Inertia::Compute` 对每个面（跳过 `IsD()` 与 `DoubleArea <= std::numeric_limits<float>::min()` 的面，`inertia.h:196`）做 Mirtich 投影积分，然后

```cpp
// inertia.h:211
T0 += fn[X] * ((A == X) ? Fa : ((B == X) ? Fb : Fc));
```

`Mass()` 直接返回 `T0`（`inertia.h:233-236`）。`CenterOfMass()` 是 `T1[k]/T0`（`inertia.h:242-249`），**分母为 0 时得到 NaN**。

注意 `A/B/C` 只是为数值稳定挑选的投影轴（`inertia.h:204-207`），`T0` 累加的物理量是**穿过曲面的 x 方向通量**。

### 2.2 可观察验证：`Mass == Σ A·n_x·c_x`

用独立的逐面求和 `Σ (DoubleArea/2)·n_x·centroid_x` 对照：

| 网格 | `Mass()` | 独立 x-flux | 解析体积 |
|---|---|---|---|
| 闭合盒子 `[-1,1]^3` | 8.0000 | 8.0000 | 8 |
| Torus R=2 r=1（64×32） | 39.1623 | 39.1623 | 39.4784（−0.80%） |
| Sphere subdiv 4 | 4.1797 | 4.1797 | 4.1888（−0.22%） |
| 盒子缺一个 +z 顶面三角形 | 8.0000 | 8.0000 | 8（**碰巧正确**） |
| 盒子缺一个 +x 墙面三角形 | 6.0000 | 6.0000 | 8（错误） |

闭合时二者相等正是散度定理：∮ x·n_x dA = ∫ ∂x/∂x dV = V。所以"体积"是由闭合性**推出的**，不是算法保证的。

### 2.3 开放网格的解析式

设网格是闭合曲面挖掉补片 P，则

```
Mass = V − ∮_P x·n_x dA
```

实测（删除单个三角形，然后沿 x 平移 +10）：

| 删除的面 | n_x·A | Mass(t=0) | Mass(t_x=+10) | Δ | 预测 −10·n_x·A |
|---|---|---|---|---|---|
| face 0 (n=+z) | 0 | 8.0000 | 8.0000 | 0 | 0 |
| face 2 (n=+x) | +2 | 6.0000 | **−14.0000** | −20 | −20 |
| face 10 (n=−x) | −2 | 6.0000 | **26.0000** | +20 | +20 |
| face 4 (n=+y) | 0 | 8.0000 | 8.0000 | 0 | 0 |

同样一个"缺一个三角形"的网格，绕 z 旋转 45° 时 `Mass` 为 7.3333，旋转 90° 又变回 8.0000：

```
open box (one +x wall tri missing), rotZ  0.0 -> Mass = 6.0000
                                 rotZ 30.0 -> Mass = 6.7887
                                 rotZ 45.0 -> Mass = 7.3333
                                 rotZ 90.0 -> Mass = 8.0000
```

**这条现象是报告后面所有"必须标记不确定"的根源**：开放网格的 `Mass` 依赖坐标系，并且可能恰好等于正确答案。任何只检查"体积是否等于期望值"的验收都会被它骗过。

---

## 3. 重复顶点（duplicate vertices）

### 3.1 实现语义

`RemoveDuplicateVertex`（`clean.h:206-276`）：

* 按 `cP()` 排序（`clean.h:222`），比较用的是 `Point3::operator==` 的**逐位精确相等**（`point3.h:369-372`）；
* 命中的顶点被 `DeleteVertex`（`clean.h:237`），并把**面/边/四面体**里指向它的指针改写为保留者（`clean.h:247-268`）；
* 末尾默认 `RemoveDegenerateFace`（`clean.h:270`），并顺带删退化边/重复边（`clean.h:271-274`）。
* 它**不重建 FF/VE/TT**。

`MergeCloseVertex(radius)`（`clean.h:1904-1910`）= `ClusterVertex`（`clean.h:1912-1945`，把半径内顶点坐标**吸附**到代表点）+ `RemoveDuplicateVertex`。

### 3.2 全缝合盒子（每个面 3 个独立顶点）：拓扑崩解但体积不变

构造：12 个三角面、36 个顶点，同一几何位置上有 4 份拷贝。

| 状态 | V | F | E | 边界边 | 孔洞环 | 连通分量 | χ | `Mass` |
|---|---|---|---|---|---|---|---|---|
| 缝合原始 | 36 | 12 | 36 | **36** | **12** | **12** | 12 | 8.0000 |
| 精确合并（删 28 个顶点） | 8 | 12 | 18 | 0 | 0 | 1 | 2 | 8.0000 |

结论：
* 边界环 36/12、分量 12 —— 一个闭合实体在拓扑上被读成 **12 张互不相连的三角片**，`IsWaterTight=0`；
* 但 `Mass` 前后都是 8.0000：**重复顶点不进入 `Mass` 的计算路径**（`Inertia::Compute` 只按每个面的 3 个顶点坐标逐面积分，`inertia.h:196-222`）。因此"体积正确"完全不能证明顶点已合并。
* 这正是仓库强调"先合并顶点、再建拓扑"的原因（`clean.h:202-206`）。

### 3.3 精确相等是脆弱的：1e-4 抖动使函数完全失效

同一构造，给每份拷贝加 ~1e-4 的各向异性抖动：

| 操作 | 删/合并数 | 结果 | `Mass` |
|---|---|---|---|
| 原始（抖动 1e-4） | — | V=36, bE=36, 孔洞 12, cc=12 | 8.0001 |
| `RemoveDuplicateVertex` | **0** | 完全无变化 | 8.0001 |
| `MergeCloseVertex(1e-3)` | 28 | V=8, bE=0, 孔洞 0, cc=1 | 8.0001 |

因此：扫描件的"接缝重复点"需要 `MergeCloseVertex`，而该函数**必须给定半径**；半径过大会成实体破坏（§3.5）。

### 3.4 合并后不重建拓扑 → 查询读到过期结果

```
seamed box, before cleanup:              CountHoles = 12
after weld, WITHOUT rebuilding FF:       CountHoles = 12   IsWaterTight = 1
after rebuilding FF:                     CountHoles = 0    IsWaterTight = 1
```

`IsWaterTight` 走 `FillEdgeVector`（`update/topology.h:200-208`）只读顶点指针，故立刻正确；`CountHoles`（`clean.h:1236-1259`）与 `CountConnectedComponents`（`clean.h:1265-1269`，分量遍历在 `1271-1306`）读 FF，于是给出**自相矛盾的两个答案**。任何"清理后顺手查询"的流程都会踩到这一点。

### 3.5 合并会**静默删面**（拓扑被改动，且不报告）

构造：加一个三角面，三个顶点是**不同对象**的拷贝——V0 拷贝盒子顶点 0、V1 拷贝盒子顶点 1、V2 再拷贝一次顶点 0（即 V0 与 V2 同位置、但指针不同），因此 `RemoveDegenerateFace` 的指针判据看不到它，`DoubleArea` 为 0。

```
before:  DoubleArea = 0，F = 13
RemoveDuplicateVertex: deleted verts = 3，F 13 -> 12（因 `clean.h:270` 的 `RemoveDegenerateFace` 静默丢弃）
after:   V=8 F=12 bE=0 WT=1 Mass=8.0000
```

返回值只报告"删了几个顶点"，**不报告"顺手删了一个面"**。这是一个真实的语义变化：一个几何退化的补片在合并后消失了，而调用方只看到顶点计数。

`MergeCloseVertex` 半径过大的破坏性更强（薄板 2×2×0.02）：

| 半径 | 合并数 | 结果 | `Mass` |
|---|---|---|---|
| 0.001 | 0 | V=8 F=12 WT=1 | 0.0800 |
| **0.02** | 0 | 不变（正好等于板厚，未触发 `dist < radius`） | 0.0800 |
| **0.05** | 4 | V=4 **F=4**，非流形边 1，`WT=0` | **0.0000** |

`ClusterVertex` 用 `dist < radius`（`clean.h:1935`）严格小于，所以阈值等于板厚时"侥幸"没塌；一旦超过，盒子被压成 4 个顶点、4 个面，体积归零，且没有警告。

### 3.6 悬空顶点：污染 `MeshGenus` 的奇偶性

`MeshGenus(m)` 直接用 `m.vn`（`clean.h:1402-1412`），而 `E` 只由面边统计（`FillEdgeVector` 只遍历面）。悬空顶点不出现在 E 中，却进入 V：

| 网格 | V | E | F | χ=V−E+F | `MeshGenus` |
|---|---|---|---|---|---|
| 闭合盒子 | 8 | 18 | 12 | 2 | 0 |
| +1 个悬空顶点 | 9 | 18 | 12 | 3 | **0**（整数除法 1/2 = 0） |
| +2 个悬空顶点 | 10 | 18 | 12 | 4 | **−1**（2/2 = 1） |
| 调用 `RemoveUnreferencedVertex` 后 | 8 | 18 | 12 | 2 | 0 |

即：**未清理的悬空顶点会让亏格按"每两个顶点 −1"的奇偶规律漂移**，且 1 个顶点时"看起来正常"。这类结果必须通过 `Clean::RemoveUnreferencedVertex`（`clean.h:396-433`）消除后才能解释。

### 3.7 分类

| 情形 | 清理后能否解释为实体 |
|---|---|
| 精确重复顶点，合并后重建 FF，`WT=1, COH=1, cc=1, Mass>0` | **可解释** |
| 抖动重复顶点，未做 `MergeCloseVertex`（仍 bE>0, cc>1） | **不可解释**（`Mass` 正确也不能采信） |
| 合并时半径过大导致 F 骤减 / nmE>0 | **必须标记不确定** |
| 合并静默删面（§3.5） | **必须标记不确定**（实体边界被改过，且无报告） |
| 存在悬空顶点未清理 | **不确定**（亏格/OE 不可解释） |

---

## 4. 零面积面与退化面

仓库明确区分了两个不同概念，这一点直接影响清理顺序。

### 4.1 两种"退化"

`RemoveDegenerateFace`（`clean.h:477-493`）的注释（`clean.h:469-476`）：

> Degenerate faces are faces that are **Topologically** degenerate, i.e. have two or more vertex reference that link the same vertex (and not only two vertexes with the same coordinates). All Degenerate faces are zero area faces BUT not all zero area faces are degenerate. We do not take care of topology because **when we have degenerate faces the topology calculation functions crash**.

判据是**指针相等**（`clean.h:484-486`）。

`RemoveFaceOutOfRangeArea(m, MinAreaThr=0)`（`clean.h:838-856`）与 `RemoveZeroAreaFace`（`clean.h:858`）用 `DoubleArea<=2*MinAreaThr`，即**面积为零**（`clean.h:847-848`）。注意 `MinAreaThr*=2`（`clean.h:841`）与 `DoubleArea` 已是两倍面积，两者相互抵消。

实测两者互补：

| 构造 | `DoubleArea` | `RemoveDegenerateFace` | `RemoveZeroAreaFace` |
|---|---|---|---|
| 同一顶点指针出现两次（V0==V1） | 0 | **删 1** | 删 1 |
| 两个不同顶点对象、同一坐标 + 一个真实第三点 | — | **0**（不删） | **删 1** |
| 三点共线（不同对象） | 0 | 0 | 删 1 |

### 4.2 拓扑退化面会让**所有**拓扑查询崩溃（debug）

```
calling IsWaterTight with a topologically degenerate face...
update/topology.h:167: PEdge::Set: Assertion `v[0] != v[1]' failed. -> abort
```

因为 `CountEdgeNum` → `FillEdgeVector` 会构造 `PEdge`，其中 `assert(v[0]!=v[1])`（`update/topology.h:167`）。**`IsWaterTight` 也逃不掉**。NDEBUG 下不崩溃，但结果不可信：

```
NDEBUG: E=19 bE=1 nmE=1 WT=0 holes=1   （真实值：闭合盒子，E=18 bE=0 nmE=0）
```

所以顺序上**必须先 `RemoveDegenerateFace`**，再建拓扑、再谈封闭性。

### 4.3 几何退化面：共面、非流形边、欧拉数

**情形 A：悬空共线面（不与实体相连）**

| 状态 | V | F | E | bE | nmE | 孔洞 | cc | χ | `Mass` |
|---|---|---|---|---|---|---|---|---|---|
| 盒子 + 共线面 | 11 | 13 | 21 | 3 | 0 | 1 | 2 | 3 | 8.0000 |
| `RemoveZeroAreaFace` | 11 | 12 | 18 | 0 | 0 | 0 | 1 | **5** | 8.0000 |
| 再 `RemoveUnreferencedVertex` | 8 | 12 | 18 | 0 | 0 | 0 | 1 | 2 | 8.0000 |

中间态 χ=5、`MeshGenus` 因整数除法给 −1（0.5 截断）。**删完面必须删悬空顶点**，否则"封闭且亏格为负"的自相矛盾结果会被输出。

**情形 B：零面积"鳍"贴到实体边上**

| 状态 | V | F | E | bE | nmE | `WT` | `Mass` |
|---|---|---|---|---|---|---|---|
| 盒子 + 鳍 | 9 | 13 | 20 | 2 | **1** | 0 | 8.0000 |
| `RemoveZeroAreaFace` | 9 | 12 | 18 | 0 | 0 | 1 | 8.0000 |

零面积面把一条实体边变成 3 面共边，从而**整个实体被判为不封闭**。清理后恢复。注意此时 `CountHoles`/`MeshGenus` 会因 nmE>0 进入未定义区域（debug 下 `Pos::FlipF` 断言失败，见 `vcg/simplex/face/pos.h:196-206`）；使用前必须先用 `CountNonManifoldEdgeFF`（`clean.h:1089`）门控。

### 4.4 数值上"零面积"是相对 `ScalarType` 的

`DoubleArea` 用 `ScalarType` 计算（float 网格即 float）。近共线三角形的高度 h 从 1e-7 降到 1e-12，`DoubleArea` 输出都是 0.0，于是 `RemoveZeroAreaFace` 一并删除：

```
h=0        DoubleArea=0.000000e+00  removed by RemoveZeroAreaFace=1
h=1e-7     DoubleArea=0.000000e+00  removed by RemoveZeroAreaFace=1
h=1e-9     DoubleArea=0.000000e+00  removed by RemoveZeroAreaFace=1
h=1e-12    DoubleArea=0.000000e+00  removed by RemoveZeroAreaFace=1
```

这是"清理的判定精度 = 网格标量精度"，不是几何精度。清理会把"实际不是零面积、但在 float 下算出来是零"的面删掉，属于**有损且不可逆**的拓扑改动。

### 4.5 真正危险的一类：**非零**面积的小三角形（不会被清理）

`Inertia::Compute` 的守卫是 `DoubleArea > std::numeric_limits<float>::min()`（`inertia.h:196`，≈1.18e−38），与 `ScalarType` 无关。因此面积 1e−3、但位于 x=1e6 的小三角形**通过守卫**并主导体积：

| 小三角形所在平面 x | `DoubleArea` | n_x | 预测 ΔV = A·n_x·X | 实测 `Mass` | `RemoveZeroAreaFace` |
|---|---|---|---|---|---|
| 1 | 1.0e−3 | +1.00 | 0.0005 | 8.0005 | 0 |
| 1e3 | 1.0e−3 | +1.00 | 0.5000 | 8.5000 | 0 |
| 1e6 | 1.0e−3 | +1.00 | **500.0** | **508.0000** | **0** |

一个闭合盒子 + 一张 0.0005 面积、远离原点的游离三角形：`E=21, bE=3, 孔洞=1, cc=2, χ=3, IsWaterTight=0, COH=1`，`Mass` 被放大 63 倍（8 → 508）。这里的 `WT=0` 来自游离片自身的 3 条边界边。

把这块游离面换成**小闭合壳**（1e−2 缩放的四面体，挂在 x=X、y=+5 处）后，所有判据都通过，但放大效应消失了：

| 小壳所在 x | F | bE | nmE | 孔洞 | cc | χ | `WT` | `COH` | `Mass` |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 16 | 0 | 0 | 0 | 2 | 4 | 1 | 1 | 8.000003 |
| 1e3 | 16 | 0 | 0 | 0 | 2 | 4 | 1 | 1 | 8.000003 |
| 1e6 | 16 | 0 | 0 | 0 | 2 | 4 | 1 | 1 | 8.000000（float 下壳已塌成零面积） |

对比同为 x=1e3、面积量级相当的**开放**三角形：`Mass=8.5000`（+0.5，随位置线性增长），而闭合小壳只增加 +3e−6（≈它自身体积，与位置无关）。这是通量分解的直接后果：

* 额外的**闭合分量**贡献 = 它自身的体积，**与位置无关**（散度定理）；
* 额外的**开放补片**贡献 = −∮_P x·n_x dA，**与位置成正比**，可以任意大。

所以"小面积"不是安全信号，但它**也不是**自动危险：这两种补片恰好被 `IsWaterTight` 区分开（开放补片有 3 条边界边，`WT=0`）。真正的教训是：

* `Inertia` 的守卫只看**绝对面积** `DoubleArea > float::min`（`inertia.h:196`），与"这个面是否属于该实体"无关；
* 一旦跳过了 `IsWaterTight` / `CountHoles` 的前置检查（例如只做 `RemoveZeroAreaFace` 就调 `Mass`），x=1e6 处的 5e−4 面积补片会让 8 变成 508；
* 反过来，`IsWaterTight` 能挡住这种补片，**却不能**挡住"两个完全重合的闭合壳"（§7.4）——二者 `Mass` 都是错误的，但只有前者被判据捕获。

因此"只清理零面积面"的策略对体积解释**既不够也不必要**：它漏掉零面积的拓扑效应（§4.2），也对非零面积的非实体补片无效。

### 4.6 分类

| 情形 | 清理后能否解释为实体 |
|---|---|
| 拓扑退化面（重复指针），先 `RemoveDegenerateFace` | **可解释** |
| 共线面已删、悬空顶点已删、FF 已重建 | **可解释** |
| 近共线面被 float 判零而删除 | **不确定**（几何上非零，被有损删除） |
| 面已删但悬空顶点未删（χ≠2、g<0） | **不确定** |
| 几何退化的定义只按"面积恰好为 0"（`clean.h:847-848`），而 `Mass` 的守卫只按 `DoubleArea > float::min`（`inertia.h:196`） | 两张判据都不覆盖"面积很小但位置很远"的补片；开放补片会被 `WT` 挡住，闭合小壳则改变体积/分量数（§4.5） |

---

## 5. 反向法线

### 5.1 唯一真值是绕序

`FlipMesh` 通过 `face::SwapEdge<FaceType,false>(f,0)` 交换 V0/V1（`clean.h:1542-1551`）。`Inertia::Compute` 只用 `NormalizedTriangleNormal`（`inertia.h:199`），**从不读取存储的 `f.N()`**。用无 `Normal3f` 组件的网格调用 `Inertia` 仍得到 `Mass=8`，`FlipMesh` 后得到 `−8`，可以确认这一点。

### 5.2 `IsWaterTight` 与朝向完全无关

下表全部是"闭合、非流形边 0、孔洞 0"，但朝向不同：

| 网格 | bE | nmE | 孔洞 | `IsWaterTight` | `IsCoherentlyOrientedMesh` | `Mass` |
|---|---|---|---|---|---|---|
| 闭合盒子 | 0 | 0 | 0 | 1 | 1 | 8.0000 |
| 翻转 1 个 +z 顶面三角 | 0 | 0 | 0 | **1** | **0** | 8.0000 |
| 翻转 1 个 +x 墙面三角 | 0 | 0 | 0 | **1** | **0** | **4.0000** |
| 翻转 1 个 +y 侧面三角 | 0 | 0 | 0 | **1** | **0** | 8.0000 |
| 全部面翻转 | 0 | 0 | 0 | **1** | **1** | **−8.0000** |

要点：
* `IsWaterTight` **必须**与 `IsCoherentlyOrientedMesh`（`clean.h:1472-1483`，逐边 `CheckOrientation`）联用，否则"封闭但朝向不一致"会被当成实体；
* 全部翻转时 `COH=1`（一致性只管相邻面，不管全局符号），`Mass` 变负 —— 符号是唯一可观察的全局朝向信号；
* 单面翻转对 `Mass` 的影响 ∝ 该面 `n_x·A`：+x 墙面 −4，+z/+y 面 0。这与 §2 的通量解释一致，也说明**"体积没变"不能证明没有反向面**。

### 5.3 修复函数的实际行为与局限

`OrientCoherentlyMesh`（`clean.h:1485-1538`）在分量内做 BFS，把相邻面翻到一致；**全局符号由起始面决定**——而起始面是 `m.face` 里第一个未访问的面（`clean.h:1494-1500`），因此修复结果依赖"被翻的面是不是第一个面"：

| 初始缺陷 | 修复后 `Mass` | `IsOriented` / `IsOrientable` |
|---|---|---|
| 翻转 face 0（即第一个面） | **−8.0000** | 0 / 1 |
| 翻转 face 1 … face 11 中任一个 | **+8.0000** | 0 / 1 |
| 全部面翻转（本身一致） | −8.0000 | **1 / 1**（不动） |

（`IsOriented=0` 表示"原本不是一致的、修复时确实翻过面"；`IsOrientable=1` 表示可定向。最后一行 `IsOriented=1` 表示"本来就是一致的"，因此函数整体不翻。）

也就是说：只修朝向一致性**不会**把体积变成正值；`Mass` 的符号仍可能为负，必须另调 `FlipNormalOutside`（且它要求单分量）。

`FlipNormalOutside`（`clean.h:1557-1604`）是**全局投票**：在 7 个方向上取极值点的法向做投票，超过一半才 `FlipMesh` 整个网格（`clean.h:1600`）。文档前置条件写得很清楚（`clean.h:1552-1556`）：

> It assumes that mesh has already has coherent normals. mesh is **watertight and signle component**.

实测其边界：

| 网格 | `Mass` 前 | `FlipNormalOutside` | `Mass` 后 |
|---|---|---|---|
| 全部翻转的盒子 | −8.0000 | flipped=1 | +8.0000 |
| 两个分离盒子，第二个反向（cc=2） | 0.0000 | flipped=1（全局翻） | **0.0000** |
| 外盒外向 + 内盒外向（cc=2） | 72.0000 | flipped=0 | 72.0000 |
| 外盒外向 + 内盒反向（空腔，cc=2） | 56.0000 | flipped=0 | 56.0000 |

关键结论：
* **多分量网格的总符号无法用仓库函数修复**。两个盒子一正一反时 `Mass=0`，`FlipNormalOutside` 只会整体翻转（0 仍是 0），`OrientCoherentlyMesh` 也不跨分量（`COH=1`，因为每分量内部都一致）。
* **空腔（nested shells）的"正确"体积依赖内外壳的相反朝向**：56 = 64 − 8 才是带空腔实体；72 = 64 + 8 是"两块材料"。仓库里**没有** API 区分"material side"，`FlipNormalOutside` 明确假设单分量，所以 56 / 72 哪一个"是实体"必须由用户语义给出，不能由代码判定。
* `trimeshinfo` 的 `if (mi.Volume < 0.0) mi.Volume = -mi.Volume;`（`trimeshinfo.cpp:819-823`）只抹掉符号，因此把 −8 报成 8，也把"一正一反"的 0 报成 0，仍无法恢复多分量语义。

### 5.4 零厚度双层片：全部"封闭性"检查都通过

用 3 个顶点、2 个三角形、绕序相反拼成一张平面片：

| V | F | E | bE | nmE | 孔洞 | cc | χ | `IsWaterTight` | `COH` | `Mass` | `CenterOfMass` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 2 | 3 | 0 | 0 | 0 | 1 | **2** | **1** | **1** | 0.0000 | **(NaN,NaN,NaN)** |

这是一个拓扑上"闭合"（每条边恰好两面、χ=2）但几何上零厚度的伪流形。`Mass=0` 可以解释为"空实体"，但 `CenterOfMass()` 因 `T0=0` 变成 NaN（`inertia.h:245-247`），`SelfIntersections` 也返回 0（§7.4）。实体判断在此必须显式标注"退化/体积为零，质心未定义"。

---

## 6. 边界孔洞

### 6.1 "封闭"的仓库定义与其盲区

```cpp
/// Very simple test of water tightness. No boundary and no non manifold edges.
/// Assume that it is orientable.
/// The rationale of not testing orientability here is that
/// it requires FFAdj while this test do not require any adjacency.
static bool IsWaterTight(MeshType & m)   // clean.h:1191-1204
{
  int edgeNum=0,edgeBorderNum=0,edgeNonManifNum=0;
  CountEdgeNum(m, edgeNum, edgeBorderNum, edgeNonManifNum);
  return (edgeBorderNum==0) && (edgeNonManifNum==0);
}
```

`CountEdgeNum` 用 `FillEdgeVector` + 排序统计每条边被几个面引用（`clean.h:1206-1232`）：恰 1 次 = 边界，>2 次 = 非流形。它**不查可定向性、不查非流形顶点**。`CountNonManifoldVertexFF`（`clean.h:1141-1191`）是独立的第二个判据。

`CountHoles`（`clean.h:1236-1259`）沿 `Pos::NextB` 走边界环，每个环计 1，**依赖 FF**，且对非流形边不安全。

### 6.2 孔洞对体积的影响：可能"碰巧正确"

§2.3 已给出公式 `Mass = V − ∮_P x·n_x dA`。所以在**一个**三角形大小的孔上：

* 缺 +z 顶面 → `Mass=8.0000`（正确）；
* 缺 +x 墙面 → `Mass=6.0000`（错误，且平移 10 后变 −14）。

两者都满足 `孔洞=1, WT=0`。也就是说：**`Mass` 在带孔网格上既可能对、也可能错，且无任何迹象区分**。这一点必须在报告/UI 中显式声明。

### 6.3 孔填充：能恢复封闭，但几何是"猜"的

`Hole::EarCuttingFill<Ear>(m, sizeHole, selected, cb)`（`hole.h:530-552`）先用 `GetInfo`（`hole.h:608-653`）找边界环，对 `size < sizeHole` 的环调用 `FillHoleEar`（`hole.h:474-527`）。`FillHoleEar` 预分配 `holeSize-2` 个面（`hole.h:482`），用优先队列选"耳"（`hole.h:489-520`），最后把没用掉的预分配面删除（`hole.h:521-526`）。

| 场景 | 填充前 | 填充后 | 体积 |
|---|---|---|---|
| 盒子缺 1 个面（3 边环） | F=11, 孔洞 1, `WT=0`, `Mass=8.0000` | F=12, 孔洞 0, `WT=1` | **8.0000（精确）** |
| 球去赤道带（102 边环） | F=4336, 孔洞 1, `Mass=3.6173` | F=4436, 孔洞 0, `WT=1` | 3.8816（真值 4.1797，**−7.1%**） |
| 球去顶盖 | `Mass=3.5266` | 孔洞 0, `WT=1` | 3.5320（真值 4.1797，**−15.5%**） |

**填充把"不可解释"变成"看起来可解释"**：`WT` 从 0 变 1，`Mass` 变成一个稳定数值，但它只是"最小权重三角化"的耳切结果，不是原始实体的体积。球顶盖的例子误差 −0.65（15.5%），而修复前后 `WT`/`孔洞`/`cc`/`χ` 全部"正确"。

阈值语义也值得一提（`hole.h:546`：`if((*ith).size < sizeHole)`，**严格小于**）：

```
3-edge boundary loop, sizeHole=3  -> returned 0  (F=11, 孔洞 1, WT=0)
3-edge boundary loop, sizeHole=30 -> returned 1  (F=12, 孔洞 0, WT=1)
```

`sizeHole` 的单位是边界**边数**（`GetInfo` 中 `++holesize`，`hole.h:636-642`），不是面积/长度。

返回值也是"处理了几个孔"，不是"成功封了几个孔"——在未合并顶点的缝合网格上，它返回 11 而 `F` 完全没变：

```
open box, every face on its own (seamed): V=33 F=11 bE=33 孔洞=11 cc=11
  EarCuttingFill returned 11                -> V=33 F=11 bE=33 孔洞=11 cc=11  (完全未改变)
  after weld                                -> V=8  F=11 bE=3  孔洞=1  cc=1
  fill after weld returned 1                -> V=8  F=12 bE=0  孔洞=0  WT=1
```

原因可以直接从代码读出：每个"3 边环"就是一个孤立三角形，它的三个耳都落在**同一张面**上，`TrivialEar::Close` 的第一条检查 `if(e0.f==e1.f) return false;`（`hole.h:215-219`；`MinimumWeightEar` 未重写它，`SelfIntersectionEar` 最后也回落到它，`hole.h:389`）必然失败；`FillHoleEar` 预分配的 `holeSize-2 = 1` 个面从未被写入，循环结束后被 `while(f!=m.face.end()) DeleteFace(...)` 删掉（`hole.h:521-526`）。而 `EarCuttingFill` 在调用前已经 `++holeCnt`（`hole.h:544-548`），所以返回值统计的是"尝试填充的孔数"，不是"成功封住的孔数"。**必须先合并顶点再填孔**，否则填充是空操作却返回非零计数。

### 6.4 填充不检查自交

`EarCuttingFill` 只有 `EarCuttingIntersectionFill`（`hole.h:559`）带自交检查，普通耳切没有。对非平面/非凸边界，"补片"可能穿出实体或与自身相交。实测的两例（球带、球盖）都未触发 `SelfIntersections`，但体积误差说明补片形状远离真实表面。因此：**填充后的 `WT=1` 只说明"闭合"，不说明"正确"。**

### 6.5 分类

| 情形 | 清理后能否解释为实体 |
|---|---|
| 无边、无 nm 边、`COH=1`、`Mass>0`、无 nm 顶点 | **可解释** |
| 有孔洞（`CountHoles>0` / `IsWaterTight=0`） | **不可解释**（即使 `Mass` 数值恰好正确） |
| 孔洞已填、但填的是非平面/大孔 | **必须标记不确定**（体积误差可达 15%+） |
| 填充前未合并顶点（`EarCuttingFill` 空操作） | **必须标记不确定**（返回值不可信） |
| 孔洞边界含非流形顶点/边（`CheckValidity` 只检查全为边界边，`hole.h:451-462`） | **必须标记不确定**（`CountHoles` 本身不安全） |

---

## 7. 方向不一致与"封闭"的其它假阳性

### 7.1 非流形顶点：`IsWaterTight` 看不见

两个实心四面体（锥体）在**一个顶点**相接，合并共享顶点后：

| 状态 | V | F | E | bE | nmE | **nmV** | 孔洞 | cc | χ | `WT` | `COH` | `Mass` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 未合并（两个独立分量） | 8 | 8 | 12 | 0 | 0 | 0 | 0 | 2 | 4 | 1 | 1 | 0.3333 |
| **合并后（共享顶点）** | 7 | 8 | 12 | 0 | 0 | **1** | 0 | 2 | 3 | **1** | 1 | 0.3333 |
| `SplitNonManifoldVertex` | 8 | 8 | 12 | 0 | 0 | 0 | 0 | 2 | 4 | 1 | 1 | 0.3333 |

`IsWaterTight=1` 但 `CountNonManifoldVertexFF=1`：这就是 `clean.h:1193-1198` 注释里承认的"不测可定向性/不管顶点"的代价。**判实体必须把 `CountNonManifoldVertexFF==0` 作为独立条件**，这正是 `trimeshinfo` 的 `VManifold`/`FManifold` 设计意图（`trimeshinfo.cpp:786-787`）。

`SplitNonManifoldVertex`（`clean.h:648-707`）能把顶点复制开、恢复 χ=2C，且 `Mass` 不变（0.3333）—— 说明这类缺陷是**可修复的**，属于"清理后仍可解释"，前提是修复后重新判定。

`RemoveNonManifoldVertex`（`clean.h:511-525`）则是**破坏性**的：对同一对锥体，它删掉 1 个面并留下 F=2、bE=6 的残骸（`WT=0`），且 `Mass` 仍是 0.3333（因为剩下的面仍构成 x 通量）。这不是"清理"，是"切掉"。

### 7.2 共享一面的两个盒子：合并顶点造成 nmE

| 步骤 | V | F | E | bE | **nmE** | `WT` | `COH` | `Mass` |
|---|---|---|---|---|---|---|---|---|
| 未合并（两分量，共享面重合） | 16 | 24 | 36 | 0 | 0 | 1 | 1 | 16.0000 |
| `RemoveDuplicateVertex` | 12 | 24 | 31 | 0 | **5** | 0 | 0 | 16.0000 |
| 再 `RemoveDuplicateFace` | 12 | 22 | 31 | 0 | **4** | 0 | 0 | **20.0000** |
| 手工删掉**两份**内壁 | 12 | 20 | 30 | 0 | 0 | 1 | 1 | 16.0000 |

* 只合并顶点 → 内壁变成 4 面共边（nmE=5），实体判定彻底失效；
* `RemoveDuplicateFace`（`clean.h:334-356`，按**顶点索引**排序三元组去重）只删一份内壁，剩下**单层内壁**，`Mass` 从 16 变 20；
* 只有"两份都删"才是真正的并集（16.0000，χ=2）。

**并且保留哪一份取决于 append 顺序**（`RemoveDuplicateFace` 删 `fvec[i]`，即排序后靠前的一份）：

```
append order 0, dedup faces(2), 存活内壁 n_x=+1 -> Mass = 20.0000
append order 1, dedup faces(2), 存活内壁 n_x=-1 -> Mass = 12.0000
```

同一个几何输入，两个仅顺序不同的构造，得到 20 与 12 两个不同体积。**这是"必须标记不确定"的典型案例：清理函数的输出依赖输入元素顺序。**

### 7.3 反向重复面：去重"删一个"，删哪个由索引顺序决定

`RemoveDuplicateFace` 对排序后的三元组只删 `fvec[i]`（`clean.h:349-355`），即**保留索引较大的一份**。因此当两份朝向相反时，去重结果取决于哪一份是后加入的：

| 构造（盒子 + face 2 的副本） | 附加顺序 | 去重前 `Mass` | 去重后 `Mass` | 说明 |
|---|---|---|---|---|
| 同向副本 | 盒子先 / 副本先 | 10.0000 | 8.0000 / 8.0000 | 结果与顺序无关 |
| 反向副本（V0,V2,V1） | 盒子先（副本索引 12） | 6.0000 | **4.0000** | 原 +x 墙面被删，**反向副本存活** |
| 反向副本（V0,V2,V1） | 副本先（盒子面索引 12） | 6.0000 | **8.0000** | 反向副本被删，原墙面存活 |

逐面核对（`order=0`，附加副本在后）：

```
after weld, before dedup: face2 n_x=+1.00, face3 n_x=+1.00, face12 n_x=-1.00, Mass=6.0000
RemoveDuplicateFace -> 1
after dedup:              face3 n_x=+1.00, face12 n_x=-1.00                Mass=4.0000
```

注意这里 `RemoveDuplicateFace` 把实体原本的 +x 墙面删掉、把**反向的重复面**留下，`Mass` 从 6 掉到 4 —— 它**制造**了新的朝向缺陷，而不是修正。同样的调用在 §7.2 的共享面场景里又制造了 20。函数只按"顶点索引三元组相同"去重，完全不看朝向，也不区分"重复面"与"实体间内壁"这两种语义。

### 7.4 `SelfIntersections` 的漏检（不确定性的放大器）

`SelfIntersections`（`clean.h:1750-1785`）用均匀网格 + `GetInBoxFace` 取候选，再 `TestFaceFaceIntersection`（`clean.h:1869-1897`）；后者在 `sv==3`（完全共顶点）时**直接返回 true**（`clean.h:1872`）。

| 场景 | `SelfIntersections` | 说明 |
|---|---|---|
| 两个**完全重合**的闭合盒子（独立顶点） | **0** | 网格包围盒查询退化：对 face 2 的精确包围盒，`GetInBoxFace` 返回 0 个候选（把 bbox 外扩 1e-4 后返回 10 个）。该网格 `WT=1, COH=1, cc=2, χ=4, Mass=16.0000` —— 体积被**重复计算**两次 |
| 两个共面重叠盒子（平移 (1,0,0)） | **0** | 共面重叠未被 `NoDivTriTriIsect` 检出 |
| 两个横向穿插盒子（平移 (1.5,0.3,0.2)） | 1（28 个面） | 正常检出 |
| 零厚度双层片 | 0 | — |

因此 `SelfIntersections=0` **不能**作为"无自交"的证明；而自交本身既破坏"实体"语义，又破坏 `Mass` 的体积解释（两个穿插盒子的 `Mass=16` 是两个体积之和，不是并集体积）。凡是 `IsWaterTight=1` 但 `cc>1` 或包围盒重叠的输入，都必须标记"可能自交，体积语义未定义"。

---

## 8. 判定流程：清理后仍可解释 vs 必须标记不确定

### 8.1 建议顺序（全部使用仓库函数）

```
1. Clean::RemoveDegenerateFace(m)            // 必须先做，否则拓扑函数 debug 下断言
2. Clean::RemoveDuplicateVertex(m)            // 或 Clean::MergeCloseVertex(m, r)，r 需人工给定
3. Clean::RemoveDuplicateFace(m)              // 必须在 2 之后（按索引判定，clean.h:334-369 注释）
4. Clean::RemoveZeroAreaFace(m)               // 几何零面积
5. Clean::RemoveUnreferencedVertex(m)         // 必须，否则 χ/亏格被悬空点污染
6. UpdateTopology::FaceFace(m)                // 必须，所有查询都依赖
7. Clean::CountNonManifoldEdgeFF(m)
   Clean::CountNonManifoldVertexFF(m)
   Clean::IsWaterTight(m)
   Clean::CountHoles(m)
   Clean::CountConnectedComponents(m)
   Clean::IsCoherentlyOrientedMesh(m)
   Clean::MeshGenus(m)
8. Inertia::Mass() / CenterOfMass()
```

如果只做"精确清理"（不做 `MergeCloseVertex`），第 2 步只能消除逐位相等的重复点；坐标有抖动时它会返回 0（§3.3）。

### 8.2 "清理后仍可解释的实体"（全部条件同时成立）

| # | 条件 | 依据 |
|---|---|---|
| E1 | 拓扑退化面计数 = 0（否则 query 不可用） | `clean.h:469-476`；`update/topology.h:167` |
| E2 | `CountNonManifoldEdgeFF(m) == 0` | `clean.h:1089` |
| E3 | `CountNonManifoldVertexFF(m) == 0` | `clean.h:1141`；§7.1 |
| E4 | `IsWaterTight(m) == true`（无边、无 nm 边） | `clean.h:1199` |
| E5 | `CountHoles(m) == 0`（与 E4 冗余地交叉验证） | `clean.h:1236` |
| E6 | `IsCoherentlyOrientedMesh(m) == true` | `clean.h:1472`；§5.2 |
| E7 | 无未清理悬空顶点（`CountUnreferencedVertex(m)==0`） | `clean.h:385`；§3.6 |
| E8 | `Mass()` 有限且 > 0 | `inertia.h:233`；§2 |
| E9 | 无零厚度/退化面（可用 `RemoveZeroAreaFace` 返回 0 + 包围盒非零厚度代理） | §5.4 |
| E10 | 单分量（或用户显式接受"材料 + 空腔"语义，见 E11） | §5.3 |
| E11 | 多分量时，`Mass` 是**带符号求和**；"空腔"（56）与"两块材料"（72）无法由仓库判定 | §5.3 |
| E12 | 包围盒重叠的两分量需人工排除自交（`SelfIntersections` 有漏检） | §7.4 |

满足 E1–E9 且单分量时，`Mass` 是**有符号体积**；配合 E8（>0）即"实体体积"。这是仓库能给出的最强保证。

### 8.3 必须标记为"不确定"的情况（清单）

| 情况 | 观察到的现象（实测） | 为什么不能采信 |
|---|---|---|
| 存在孔洞（`CountHoles>0` 或 `WT=0`） | 缺 +z 面 `Mass=8`（对）；缺 +x 面 `Mass=6`；平移后 `−14` | `Mass` 是通量，退化为 `V−∮_P x n_x`；可能**碰巧正确**（§2.3） |
| 孔洞已填 | 球顶盖填充后 `WT=1`、误差 −15.5% | 补片是耳切猜测，非真实表面（§6.3） |
| 填充前未合并顶点 | `EarCuttingFill` 返回 11 但 F/bE/孔洞全不变 | 返回值不能作为"已封孔"证据（§6.3） |
| `WT=1` 但 `nmV>0` | 两锥共享顶点：`WT=1, nmV=1` | `IsWaterTight` 不检查非流形顶点（§7.1） |
| `WT=1` 但 `cc>1` 且体积相互抵消 | 两盒一正一反：`WT=1, COH=1, Mass=0` | 符号无法由 `OrientCoherentlyMesh`/`FlipNormalOutside` 跨分量修复（§5.3） |
| 多分量嵌套（空腔） | 56 vs 72 两者都 `WT=1, COH=1, χ=4` | 无"材料侧"API；`FlipNormalOutside` 假设单分量（§5.3） |
| `Mass` 对坐标敏感（开放网格） | 同一网格 rotZ 0/30/45/90 → 6.0/6.79/7.33/8.0 | 坐标系依赖即非体积（§2.3） |
| 零厚度双层片 | `WT=1, COH=1, χ=2, Mass=0`，`CoM=(NaN,NaN,NaN)` | 拓扑"封闭"但无实体；质心未定义（§5.4） |
| 小面积、远离原点的补片 | 面积 5e−4 的三角形在 x=1e6 → `Mass=508`（真值 8） | `Inertia` 的守卫只按 `DoubleArea > float::min`，与位置无关（§4.5） |
| 近共线面被 float 判零 | `h=1e−7` 与 `1e−12` 都 `DoubleArea=0`，被 `RemoveZeroAreaFace` 删除 | 清理精度 = 标量精度，有损删除（§4.4） |
| 合并顶点时静默删面 | `deleted verts=3` 同时 F 13→12 | 返回值不报告面被删（§3.5） |
| `MergeCloseVertex` 半径过大 | 薄板 r=0.05 → F=4, `WT=0`, `Mass=0` | 拓扑被压塌且无告警（§3.5） |
| 清理后未重建 FF | 同一时刻 `IsWaterTight=1` 与 `CountHoles=12` | 两个答案自相矛盾（§3.4） |
| 悬空顶点未删（尤其 1 个） | +1 顶点 `g=0`（假正常）；+2 顶点 `g=−1` | `MeshGenus` 用 `m.vn`，E 只含面边（§3.6） |
| 存在未清理的非流形边 | 零面积鳍贴合实体 → 整个盒子 `WT=0`；`CountHoles`/`MeshGenus` 进入未定义区 | debug 下 `Pos::FlipF` 断言失败（§4.3） |
| 包围盒重叠 / 共面重叠 | `SelfIntersections` 对完全重合与共面重叠返回 0 | 漏检；体积语义未定义（§7.4） |
| NaN / Inf 坐标 | NaN 顶点：`WT=1` 但 `Mass=6.0`；Inf 顶点：`WT=1` 但 `Mass=NaN` | `Inertia` 静默跳过 `DoubleArea>min` 为假的面；`RemoveDegenerateVertex` 只处理 NaN（`clean.h:439-467`），Inf 不被清理 |
| 相同的几何、不同的构造顺序 | 共享面盒子：append 顺序 0/1 → `Mass=20` / `12` | `RemoveDuplicateFace` 保留排序靠前者（§7.2） |

### 8.4 关于仓库中已有的"门控"写法（不应作为规范）

`apps/unsupported/trimeshinfo/trimeshinfo.cpp:815` 用

```cpp
if ((mi.VManifold && mi.FManifold) && (mi.Oriented) && (!mi.numholes))
```

作为"体积可解释"的闸门，并在 `:821-823` 取绝对值。但该文件在 `apps/CMakeLists.txt` 中未被构建（只加入 `sample, metro, tridecimator, embree`），且其调用与当前 `Clean` API 已不一致（实测编译错误）：

* `RemoveFaceOutOfRangeArea<false>(m,0)`（`:776`）—— 模板参数已不存在，按类型比较报错；
* `CountEdges(m, …)`（`:790`, `:796`）—— `Clean` 中无此成员（现为 `CountEdgeNum`）；
* `IsOrientedMesh(m, …)`（`:804`）—— 无此成员（现为 `OrientCoherentlyMesh` 的 out 参数形式）。

同时 `:786-787` 把 `CountNonManifoldVertexFF(m) > 0` 赋给 `VManifold`、`CountNonManifoldEdgeFF(m) > 0` 赋给 `FManifold`，而 `:391-394` 在 `VManifold && FManifold` 时打印 `Volume: <值>`、`:418-421` 打印 `Manifold: YES` —— 即"检出非流形"被当成"是流形"，布尔含义与该名字相反（XML 输出 `:456-459` 同样反向）。这些可以说明"用单一布尔量做实体判定的脆弱性"，但不构成可用参考实现。

---

## 9. 一张总表（全部为实测，`Mass` 由 `Inertia::Mass()`）

| 场景 | V | F | E | bE | nmE | nmV | 孔洞 | cc | χ | `WT` | `COH` | `Mass` | 实体可解释？ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 闭合盒子 | 8 | 12 | 18 | 0 | 0 | 0 | 0 | 1 | 2 | 1 | 1 | 8.0000 | ✅ |
| 闭合 torus（亏格 1） | 2048 | 4096 | 6144 | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 1 | 39.1623 | ✅ |
| 闭合 sphere | 2562 | 5120 | 7680 | 0 | 0 | 0 | 0 | 1 | 2 | 1 | 1 | 4.1797 | ✅ |
| 缺 +z 顶面 1 个三角 | 8 | 11 | 18 | 3 | 0 | 0 | 1 | 1 | 1 | 0 | 1 | 8.0000 | ❌（碰巧对） |
| 缺 +x 墙面 1 个三角 | 8 | 11 | 18 | 3 | 0 | 0 | 1 | 1 | 1 | 0 | 1 | 6.0000 | ❌ |
| 上者平移 x+10 | 8 | 11 | 18 | 3 | 0 | 0 | 1 | 1 | 1 | 0 | 1 | −14.0000 | ❌ |
| 全缝合盒子（未合并） | 36 | 12 | 36 | 36 | 0 | 0 | 12 | 12 | 12 | 0 | 1 | 8.0000 | ❌ |
| 全缝合盒子（合并后） | 8 | 12 | 18 | 0 | 0 | 0 | 0 | 1 | 2 | 1 | 1 | 8.0000 | ✅ |
| 缝合 + 1e−4 抖动，未 `MergeCloseVertex` | 36 | 12 | 36 | 36 | 0 | 0 | 12 | 12 | 12 | 0 | 1 | 8.0001 | ❌ |
| 闭合盒 + 2 悬空点 | 10 | 12 | 18 | 0 | 0 | 0 | 0 | 1 | 4 | 1 | 1 | 8.0000 | ❌（g=−1） |
| 盒 + 悬空共线面（删面、未删点） | 11 | 12 | 18 | 0 | 0 | 0 | 0 | 1 | 5 | 1 | 1 | 8.0000 | ❌（χ=5） |
| 盒 + 零面积鳍（贴边） | 9 | 13 | 20 | 2 | 1 | 0 | — | 1 | 2 | 0 | 0 | 8.0000 | ❌ |
| 同上，`RemoveZeroAreaFace` 后 | 9 | 12 | 18 | 0 | 0 | 0 | 0 | 1 | 3 | 1 | 1 | 8.0000 | ❌（χ=3） |
| 盒 + 重复 +x 墙面，反向 | 8 | 13 | 18 | 0 | 3 | 0 | — | 1 | 3 | 0 | 0 | 6.0000 | ❌ |
| 闭合盒，翻转 1 个 +x 面 | 8 | 12 | 18 | 0 | 0 | 0 | 0 | 1 | 2 | **1** | **0** | 4.0000 | ❌ |
| 闭合盒，全部翻转 | 8 | 12 | 18 | 0 | 0 | 0 | 0 | 1 | 2 | 1 | **1** | −8.0000 | ⚠️（符号） |
| 两盒一正一反 | 16 | 24 | 36 | 0 | 0 | 0 | 0 | 2 | 4 | **1** | **1** | **0.0000** | ❌ |
| 外盒 + 内盒同向 | 16 | 24 | 36 | 0 | 0 | 0 | 0 | 2 | 4 | 1 | 1 | 72.0000 | ⚠️（语义） |
| 外盒 + 内盒反向（空腔） | 16 | 24 | 36 | 0 | 0 | 0 | 0 | 2 | 4 | 1 | 1 | 56.0000 | ⚠️（语义） |
| 零厚度双层片 | 3 | 2 | 3 | 0 | 0 | 0 | 0 | 1 | 2 | **1** | **1** | 0.0000（CoM=NaN） | ❌ |
| 两面锥共享顶点（合并后） | 7 | 8 | 12 | 0 | 0 | **1** | 0 | 2 | 3 | **1** | 1 | 0.3333 | ❌（nmV=1） |
| 共享一面的两盒（只合并顶点） | 12 | 24 | 31 | 0 | 5 | 0 | — | 1 | 5 | 0 | 0 | 16.0000 | ❌ |
| 同上 + `RemoveDuplicateFace` | 12 | 22 | 31 | 0 | 4 | 0 | — | 1 | 3 | 0 | 0 | **20.0000** | ❌ |
| 同上，内壁两份都删 | 12 | 20 | 30 | 0 | 0 | 0 | 0 | 1 | 2 | 1 | 1 | 16.0000 | ✅ |
| 盒 + 远场小三角形（x=1e6） | 11 | 13 | 21 | 3 | 0 | 0 | 1 | 2 | 3 | 0 | 1 | **508.0000** | ❌ |
| 盒 + 远处小闭合壳（x=1e3） | 12 | 16 | 24 | 0 | 0 | 0 | 0 | 2 | 4 | 1 | 1 | 8.000003 | ⚠️（多分量） |
| 两盒**完全重合**（未合并，cc=2） | 16 | 24 | 36 | 0 | 0 | 0 | 0 | 2 | 4 | **1** | **1** | **16.0000** | ❌（体积重复计） |
| 同两个盒合并顶点后（重复面） | 8 | 24 | 18 | 0 | **18** | 0 | — | 1 | **14** | 0 | 0 | 16.0000 | ❌ |

（"—"表示因 nmE>0 跳过 `CountHoles`/`MeshGenus`；debug 构建下对这类输入遍历会触发 `Pos::FlipF` 断言。）

---

## 10. 给使用者的最小规则集

1. **不要单独使用 `Mass()`**。至少在调用前断言 `IsWaterTight && CountNonManifoldEdgeFF==0 && CountNonManifoldVertexFF==0 && IsCoherentlyOrientedMesh && CountHoles==0 && CountUnreferencedVertex==0`。
2. **`IsWaterTight` 不是"可解释实体"的充分条件**（§7.1、§5.4、§7.4），必须叠加 `COH` 与 `nmV`。
3. **`Mass` 正确不能反证清理正确**（§3.2、§6.2）。
4. **清理后必须重建 FF**（§3.4），且**先拓扑退化面、后拓扑查询**（§4.2）。
5. **带孔网格的体积必须报"未定义"**，不能因为数值合理就采信（§2.3）。
6. **孔填充后必须报"体积不确定"**，除非另有原始实体的对照（§6.3）。
7. **多分量/嵌套/自交的输入必须由用户给出语义**，仓库函数无法判定（§5.3、§7.2、§7.4）。
8. **任何依赖元素顺序的清理结果（§7.2、§7.3）都要在报告中声明**：同一份共享面几何，`RemoveDuplicateFace` 后可以是 20 或 12；同一份反向重复面，去重后可以是 4 或 8。
