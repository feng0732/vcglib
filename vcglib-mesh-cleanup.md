# 网格清理对实体拓扑与体积解释的影响（vcglib 源码分析）

本报告基于仓库中 `vcg/complex/algorithms/clean.h`、`inertia.h`、`stat.h`、`hole.h`、`inside.h`、
`update/topology.h`、`update/flag.h` 与 `vcg/simplex/face/topology.h` 的实现，分析四类典型缺陷
（重复顶点、零面积面、反向法线、边界孔洞）如何影响"该网格是否表示一个实体"的判断以及体积的
可解释性，并区分**清理后仍可解释为实体**与**必须标记为不确定**的情形。讨论仅限于仓库已有的
网格语义。

---

## 1. 仓库中"实体判断"的语义基础

vcglib 没有显式的 solid 类型；"实体性"由一组可观察的拓扑/几何量共同构成：

| 可观察量 | 实现 | 判定含义 |
|---|---|---|
| 水密性 | `Clean::IsWaterTight`（clean.h:1199-1204） | 无边界边且无非流形边（`edgeBorderNum==0 && edgeNonManifNum==0`） |
| 边计数 | `Clean::CountEdgeNum`（clean.h:1206-1232） | 按顶点对排序统计每条边被面引用的次数：1 次=边界边，>2 次=非流形边 |
| 孔洞数 | `Clean::CountHoles`（clean.h:1236-1259） | 沿边界 Pos 环（`NextB`）计边界环数 |
| 连通分量 | `Clean::ConnectedComponents`（clean.h:1271-1306） | 基于 FF 邻接的洪泛 |
| 亏格 | `Clean::MeshGenus`（clean.h:1397-1412） | `V-E+F = 2C-2G-B` |
| 方向一致性 | `Clean::IsCoherentlyOrientedMesh` / `OrientCoherentlyMesh`（clean.h:1472-1538） | 共享边被两面以相反方向遍历 |
| 有符号体积 | `Inertia::Mass`（inertia.h:233-236）、`Stat::ComputeMeshVolume`（stat.h:285-288） | 散度定理面积分，注释明确 "Meaningful only if the mesh is watertight"（inertia.h:231） |
| 内外测试 | `Inside::Is_Inside`（inside.h:43-80） | 用最近面法线符号判断内外，依赖法线朝外 |

关键约定：

- **体积是有符号积分**。`Inertia::Compute`（inertia.h:189-227）对每个面累加
  `T0 += fn[X] * Fa`（inertia.h:211），其中 `fn` 是由顶点绕序决定的法线
  （`vcg::NormalizedTriangleNormal`，inertia.h:199）。翻转一个面，其贡献变号。
- **方向一致性约定**在 `face::CheckOrientation`（simplex/face/topology.h:380-392）：
  共享边两侧的面必须以相反方向经过该边（`f.V0(z) == g->V1(gi)`）。
- **水密性测试不检查可定向性**，注释说明："Assume that it is orientable. It could be
  debated if a closed non orientable surface is watertight or not... it requires FFAdj
  while this test do not require any adjacency"（clean.h:1192-1198）。因此
  `IsWaterTight==true` 是实体性的必要而非充分条件。
- **FF 邻接对非流形边的处理**：`UpdateTopology::FaceFace`（update/topology.h:395-435）
  把共享同一顶点对的所有面按排序顺序连成一个环，>2 个面共享边时邻接关系成为环状，
  后续基于 Pos 的遍历（孔洞计数、定向传播）在此类边上行为退化。

---

## 2. 四类缺陷的影响与清理语义

### 2.1 重复顶点 —— `Clean::RemoveDuplicateVertex`（clean.h:206-276）

**缺陷影响（可观察结果）**：

- 顶点按**指针/索引**而非坐标参与拓扑匹配（`FillEdgeVector`，update/topology.h:200-209；
  `PEdge` 按顶点索引排序比较）。两个坐标相同但指针不同的顶点使本应配对的边无法匹配：
  每条边只被引用 1 次 → `CountEdgeNum` 的 `boundary_e` 增大 → `IsWaterTight` 返回 false。
- 裂缝被 `CountHoles` 计为边界环（孔洞数虚增），`MeshGenus` 公式中 `B` 增大、`G` 被污染。
- FF 邻接沿裂缝断开，`ConnectedComponents` 可能把一个壳计成多个分量。
- 体积数值本身**不变**（散度积分只看面几何），但因为水密性失败，按 inertia.h:231 的
  注释该体积**不可解释**为实体体积。

**清理语义**：

- 按坐标排序合并顶点（`RemoveDuplicateVert_Compare`，clean.h:193-199），并重写所有面、
  边、四面体的顶点指针（clean.h:247-268，四面体见 264-268）。
- 注释明确顺序约束："it does not update any topology relation... usually performed
  BEFORE building any topology information"（clean.h:202-205）——必须在建拓扑前调用。
- 级联副作用：默认 `RemoveDegenerateFlag=true`，合并后产生两个顶点引用相同的面会被
  `RemoveDegenerateFace` 删除（clean.h:270）。**可观察结果**：`fn` 减少，Euler 示性数
  改变；若这些退化面原本是"粘合"两瓣曲面的零厚度连接，删除后拓扑结构发生实质变化。

**结论**：重复顶点清理是把"几何闭合但拓扑断开"的网格恢复为可判断实体的**必要前提**；
清理后若其余条件满足，`IsWaterTight` 恢复 true，体积重新可解释。

### 2.2 零面积面 / 退化面 —— `RemoveDegenerateFace`（clean.h:477-493）与 `RemoveZeroAreaFace`（clean.h:858）

**仓库的区分**（clean.h:469-476 注释）：

> "Degenerate faces are faces that are Topologically degenerate... All Degenerate faces
> are zero area faces BUT not all zero area faces are degenerate. We do not take care of
> topology because when we have degenerate faces the topology calculation functions crash."

- **拓扑退化面**（两个顶点引用相同）：必须删除，否则拓扑计算直接崩溃——这是硬性前提。
- **几何零面积面**（三顶点共线但引用不同）：`RemoveFaceOutOfRangeArea(m,0)`
  （clean.h:838-858）按 `DoubleArea` 阈值删除。

**缺陷影响（可观察结果）**：

- 零面积面对体积贡献恒为零（`DoubleArea` 为叉积模长，space/triangle3.h:241；法线为零
  向量），因此 `Mass()` 数值不受其影响。
- 但零面积面的归一化法线是 0/0（`NormalizedTriangleNormal`，triangle3.h:88），会污染
  顶点法线平均、`FlipNormalOutside` 的投票（clean.h:1591-1603）和 `Inside::Is_Inside`
  的符号判断（inside.h:60-66）。

**清理的双刃剑效应**：

- 删除零面积面**不改变体积数值**，但如果该面在拓扑上承担着闭合曲面的连接角色
  （例如细缝两侧的极限三角形），删除后 `CountEdgeNum` 的 `boundary_e` 增加、
  `CountHoles` 增加、`IsWaterTight` 由 true 变 false。
- **可观察的判定翻转**：清理前网格"水密但含病态面"，体积形式上可解释；清理后网格
  "健康但带孔"，体积反而不可解释。是否删除应结合 `CountEdgeNum` 前后对比判断。

### 2.3 反向法线 / 方向不一致 —— `IsCoherentlyOrientedMesh`、`OrientCoherentlyMesh`、`FlipNormalOutside`

**缺陷影响（可观察结果）**：

- 体积是有符号量：部分面翻转时，`Mass()` 等于"外向部分体积减去内翻部分扫过的体积"，
  不对应任何实体；全部面翻转时 `Mass()` 为负，但 `CenterOfMass()`（`T1/T0`，
  inertia.h:242-249）仍然正确——**质心对整体翻转不变、体积变号**是一个可观察的诊断特征。
- `IsCoherentlyOrientedMesh`（clean.h:1472-1483）逐边调用 `CheckOrientation`，
  任一边方向不相反即返回 false。
- `Inside::Is_Inside` 用最近面法线与 `(test-nearest)` 的点积符号判内外
  （inside.h:60-66）；方向不一致时内外判定逐面翻转，结果无意义。

**清理语义**：

- `OrientCoherentlyMesh`（clean.h:1485-1538）沿 FF 邻接洪泛，跨**流形边**传播方向
  （`face::IsManifold` 检查，clean.h:1509），发现方向冲突且邻面未访问时执行
  `face::SwapEdge` 翻转（clean.h:1518）。输出两个布尔：
  - `_IsOriented==false`：原网格方向不一致，但已被修复（可定向）；
  - `_IsOrientable==false`：发现已访问面需要再次翻转（clean.h:1519-1523），即不可定向
    曲面（如 Möbius 带）——**必须标记为不确定**，仓库内无任何后续例程能赋予其实体语义。
- `FlipNormalOutside`（clean.h:1557-1604）用 7 个方向上极值顶点的法线投票决定是否
  整体翻转，但其前置假设写在注释里："mesh has already has coherent normals. mesh is
  watertight and single component"（clean.h:1552-1556）。在不满足假设的网格上调用，
  投票结果无定义。
- 注意 `RemoveDuplicateFace`（clean.h:334-356）的 `SortedTriple` 对三顶点索引排序
  （clean.h:310），**反向重复面（同三顶点、绕序相反）也被视为重复**并删除其一。
  对双层叠合的"薄壳"网格，该清理会把双层删成单层，边界边数爆炸——可观察为清理后
  `IsWaterTight` 反而失败。

### 2.4 边界孔洞 —— `CountHoles`（clean.h:1236）与 `Hole::EarCuttingFill`（hole.h:529-556）

**缺陷影响（可观察结果）**：

- 每个边界环使 `IsWaterTight` 失败；`MeshGenus` 中 `B` 增大。
- 体积积分缺少"孔洞盖面"的贡献：`Mass()` 的偏差等于缺失盖片所扫过的有符号体积，
  大小和符号都无先验界限 → **未补洞的带孔网格体积必须标记为不确定**。
- 对照：`Stat::ComputeShellBarycenter`（stat.h:259-272）注释明确 "Works for any
  triangulated model (no problem with open, nonmanifold selfintersecting models)"——
  薄壳质心是带孔网格上**唯一仍然良定义**的整体几何量；`ComputeMeshArea`
  （stat.h:290-299）同样始终可算但不说明实体性。

**清理（补洞）语义**：

- `Hole::EarCuttingFill<EAR>` 按边界边数（`Info::size`，hole.h:417-432）填充小于阈值
  的洞；`EAR` 可选 `TrivialEar` / `MinimumWeightEar` / `SelfIntersectionEar`，
  `MinimumWeightFill`（hole.h:850）用最小权三角剖分选择盖面几何。
- `TrivialEar` 的前置条件写在注释里："Require FF-adajcncy and edge-manifoldness around
  the mesh (at most two triangles per edge)"（hole.h:52-54）——非流形边未清理时不能补洞。
- 补洞后 `IsWaterTight` 可恢复 true，体积重新可解释，**但体积数值包含补洞算法对缺失
  几何的假设**（不同 EAR 策略给出不同盖面、不同体积）。报告体积时应注明该依赖性。

---

## 3. 三种网格情形的判定流程

### 3.1 闭合流形网格

清理链：`RemoveDuplicateVertex` → `RemoveDegenerateFace`（级联）→ `RemoveDuplicateFace`
→ 建 FF 邻接 → `OrientCoherentlyMesh` → `FlipNormalOutside`。

- 可观察终态：`IsWaterTight==true`、`IsCoherentlyOrientedMesh==true`、
  `CountNonManifoldEdgeFF/CountNonManifoldVertexFF==0`、`MeshGenus` 为有限整数。
- `Mass()` 的绝对值即实体体积，符号指示法线整体朝向；参考实现用解析解校验：
  `apps/sample/trimesh_inertia/trimesh_inertia.cpp` 中 2×4×10 盒子的 `Mass()==80`、
  环面 `Mass()≈2π²a²c≈39.478`。
- **结论：清理后可解释为实体。** 自相交是例外：`SelfIntersections`
  （clean.h:1750-1785）非空时，水密性与方向一致性都可能通过，但体积按环绕数加权，
  不等于任何实体体积——必须标记为不确定。

### 3.2 带孔网格

- 去重、去退化、定向一致化均可正常执行，且执行后 `CountHoles` 不变（这些操作不补洞）。
- 可观察终态：`IsWaterTight==false`（`edgeBorderNum>0`）、`CountHoles==B>0`、
  `MeshGenus` 公式中 `B` 项非零。注释提醒："A closed sphere, a sphere with one hole
  (e.g. a disk) and a sphere with two holes (e.g. a tube) all of them have Genus == 0"
  （clean.h:1392-1393）——亏格不能区分闭合与带孔。
- **结论：未补洞时体积必须标记为不确定**（缺项无界）；薄壳质心、面积仍可报告。
  显式调用 `EarCuttingFill` 后转入 3.1 情形，但体积依赖补洞假设。

### 3.3 方向不一致网格

- 可观察：`IsCoherentlyOrientedMesh==false`；`Mass()` 为带正负抵消的伪体积；
  `Inside::Is_Inside` 结果逐区域翻转。
- `OrientCoherentlyMesh` 修复后：若 `_IsOrientable==true`，网格转为 3.1 或 3.2 情形，
  体积可解释（注意此时 `Mass()` 符号可能整体为负，由 `FlipNormalOutside` 归一化）。
- 若 `_IsOrientable==false`（不可定向）或网格含非流形边导致定向传播无法覆盖
  （`OrientCoherentlyMesh` 只跨 `IsManifold` 边传播，clean.h:1509），**必须标记为
  不确定**。非流形情形可先尝试 `SplitNonManifoldVertex`（clean.h:648-707，复制顶点并
  按 `moveThreshold` 分离）或 `RemoveNonManifoldFace`（clean.h:799-835，按面积升序删面）
  再重新判定；后者会制造孔洞，转入 3.2 情形。

---

## 4. 判定矩阵：可解释为实体 vs 必须标记为不确定

| 状态（清理后） | 可观察证据 | 实体判定 | 体积判定 |
|---|---|---|---|
| 水密 + 方向一致 + 无自交 | `IsWaterTight && IsCoherentlyOrientedMesh && !SelfIntersections` | 是 | `|Mass()|` 有效 |
| 水密 + 方向一致 + 自交 | `SelfIntersections` 非空（clean.h:1750） | 否 | 不确定（环绕数加权） |
| 水密 + 不可定向 | `OrientCoherentlyMesh` 给出 `_IsOrientable==false`（clean.h:1521） | 否 | 不确定（`IsWaterTight` 注释声明不覆盖此情形，clean.h:1192-1198） |
| 带孔（未补洞） | `CountEdgeNum` 的 `boundary_e>0`、`CountHoles>0` | 否 | 不确定（缺项无界）；`ComputeShellBarycenter`、面积仍有效 |
| 带孔（已 `EarCuttingFill`） | 补洞后 `IsWaterTight==true` | 是（在补洞假设下） | 有效但依赖 EAR 策略 |
| 非流形边/顶点 | `CountNonManifoldEdgeFF`/`CountNonManifoldVertexFF>0`（clean.h:1089、1141） | 否 | 不确定；`RemoveNonManifoldFace`/`SplitNonManifoldVertex` 后重新判定 |
| 双层叠合薄壳被 `RemoveDuplicateFace` 删成单层 | 清理后 `boundary_e` 激增 | 否（清理引入的破坏） | 不确定；应跳过 `RemoveDuplicateFace` 或改用选择性交手 |

---

## 5. 操作建议的清理顺序（依据源码中的前置条件注释）

1. `RemoveDegenerateVertex`（NaN 坐标，clean.h:439-467）；
2. `RemoveDuplicateVertex`（必须在建拓扑前，clean.h:202-205）——级联删除拓扑退化面；
3. `RemoveDuplicateFace`（注意 2.3 节所述反向重复面风险）；
4. `RemoveUnreferencedVertex`（clean.h:396-433）；
5. 建 FF 邻接（`UpdateTopology::FaceFace`），统计 `CountEdgeNum`/`CountHoles`/
   `CountNonManifoldEdgeFF`/`CountNonManifoldVertexFF` 作为基线；
6. 非流形处理（`SplitNonManifoldVertex` 或 `RemoveNonManifoldFace`）；
7. `OrientCoherentlyMesh` → 检查 `_IsOrientable`；
8. （可选）`Hole::EarCuttingFill` 补洞；
9. `FlipNormalOutside`（要求水密、单分量、方向已一致，clean.h:1552-1556）；
10. 最终判定：`IsWaterTight && IsCoherentlyOrientedMesh && !SelfIntersections`
    三者同时成立时，`Inertia::Mass()` 才可报告为实体体积，否则标记为不确定。

---

## 附：主要源码证据索引

- 水密性与边计数：`vcg/complex/algorithms/clean.h:1192-1232`
- 孔洞计数：`vcg/complex/algorithms/clean.h:1236-1259`
- 亏格公式：`vcg/complex/algorithms/clean.h:1362-1412`
- 重复顶点/面清理：`vcg/complex/algorithms/clean.h:202-276、329-356`
- 退化面语义注释：`vcg/complex/algorithms/clean.h:469-476`
- 零面积面删除：`vcg/complex/algorithms/clean.h:837-858`
- 非流形计数与清理：`vcg/complex/algorithms/clean.h:1089-1191、511-525、648-707、799-835`
- 定向：`vcg/complex/algorithms/clean.h:1472-1604`；`vcg/simplex/face/topology.h:380-392`
- 体积/惯量：`vcg/complex/algorithms/inertia.h:184-249`；`vcg/complex/algorithms/stat.h:274-299`
- 补洞：`vcg/complex/algorithms/hole.h:52-54、398-556、850-901`
- 内外测试：`vcg/complex/algorithms/inside.h:43-80`
- FF 邻接构建（非流形边成环）：`vcg/complex/algorithms/update/topology.h:395-435`
- 参考数值（盒子体积 80、环面 ≈39.478）：`apps/sample/trimesh_inertia/trimesh_inertia.cpp`
