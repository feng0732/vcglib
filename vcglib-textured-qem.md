# 带纹理三角网格的 QEM 边折叠简化：机制、约束与现状分析

> 分析对象：本仓库（VCGLib）中带每楔（wedge）纹理坐标的三角网格二次误差测度（Quadric Error Metric, QEM）边折叠实现。
> 所有结论均来自仓库内源码静态阅读；未修改仓库任何文件，未使用外部服务，未运行构建（本报告为只读分析）。
> 文中 `文件:行号` 均相对仓库根目录。

---

## 1. 涉及的源码与调用链

| 角色 | 文件 |
|---|---|
| 带纹理 5D QEM 折叠算子（核心） | `vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric_tex.h`（802 行） |
| 5 维二次型（xyz+uv 联合误差） | `vcg/math/quadric5.h`（738 行） |
| 纯几何 3D QEM 折叠算子（对照） | `vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h`（655 行） |
| 折叠基类 / 可行性接口 / 堆调度 | `vcg/complex/algorithms/local_optimization/tri_edge_collapse.h`、`vcg/complex/algorithms/local_optimization.h` |
| 实际折叠执行与 Link 条件（拓扑有效性） | `vcg/complex/algorithms/edge_collapse.h`（331 行） |
| 退化/重复元素清理 | `vcg/complex/algorithms/clean.h` |
| 三角形质量（形状）度量 | `vcg/space/triangle3.h:188`（`QualityFace`） |
| 唯一使用该算子的示例程序 | `apps/unsupported/texdeci/texdeci.cpp`（**未纳入 CMake 构建**，仅有 `.pro`） |
| 纯几何简化器（对照，已纳入构建） | `apps/tridecimator/tridecimator.cpp` |

调用链（以 `texdeci.cpp:94-110` 为典型用法）：

```
LocalOptimization::Init<MyTriEdgeCollapseQTex>()
  -> TriEdgeCollapseQuadricTex::Init()           // VF 拓扑、边界标志、InitQuadric、建堆
LocalOptimization::DoOptimization()              // local_optimization.h:201-232
  -> pop_heap -> IsUpToDate() -> IsFeasible() -> Execute() -> UpdateHeap()
```

---

## 2. 边折叠候选如何处理几何误差

### 2.1 数据结构：每个顶点持有"几何二次型 + 按 UV 分组的 5D 二次型"

`QuadricTexHelper`（`tri_edge_collapse_quadric_tex.h:90-172`）为每个顶点维护两份临时数据：

- 一个经典 3D 二次型 `math::Quadric<double>`（`Qd3`，`...tex.h:161`）；
- 一个 `vector<pair<TexCoord2f, Quadric5<double>>>`（`Vd`，`...tex.h:96,164`）：**顶点上每出现一个不同的 (u,v)，就有一个独立的 5D 二次型**。这是处理纹理接缝的关键设计——接缝处同一 3D 顶点在不同面片上携带不同 UV，它们互不合并。

### 2.2 初始化：5D 二次型按面片构建

`InitQuadric`（`...tex.h:563-590`）对每个非删除、可读且三顶点均可读的面调用 `Quadric5::byFace`（`quadric5.h:218-271`）：

1. 先以 `onlygeo=true` 构建仅含 xyz 的几何二次型，并通过 `AddtoQ3` 累加到三个顶点的 3D 二次型（`quadric5.h:196-213, 227-230`）。
2. 再以 `onlygeo=false` 构建完整 5D 二次型，即把每个顶点视为 5 维点 `(x,y,z,u,v)`，误差是该点到三角形在 5D 空间所张平面的平方距离（`ComputeQuadricFromE1E2`，`quadric5.h:407-451`，对应 Garland–Heckbert 98 论文记号）。**几何误差与纹理误差因此被耦合进同一个二次型**，而非两个独立加权项。
3. 对边界边（或启用 `QualityQuadric` 时的所有边）追加额外 5D 二次型以约束边不漂移：边界边按 `BoundaryWeight`（默认 0.5）缩放，质量约束按 0.05 缩放（`quadric5.h:232-257`）。
4. 对顶点上每个可写的不同 UV 调用 `Alloc`/`SumAll`（`...tex.h:105-131, 576-586`）。`Sum3`（`quadric5.h:665-690`）把 3D 几何二次型"提升"进 5D：在 uv 对角块上各加 1（`a[12]+=1; a[14]+=1`），即对该 (u,v) 施加保持原位的二次约束。源码注释明确承认：这是在模拟本不存在的 UV 连续性，**接缝顶点的优先级因此估算不准，需要用户用 `ExtraTCoordWeight` 补偿**（`quadric5.h:655-664`）。

### 2.3 候选优先级：模拟折叠 → 求最优 5D 落点 → 质量惩罚

`ComputePriority` → `ComputeMinimalsAndPriority`（`...tex.h:384-466`）：

1. `GetTexCoords`（`...tex.h:239-276`）扫描折叠边 v0–v1 上**至多两个**共面邻面，收集每侧 v0/v1 的楔 UV。两组 UV 完全相同返回 1（连续区域），否则返回 2（接缝边）。
2. 对每组 UV，把两端对应二次型相加（`qsum_1`/`qsum_2`），用 `Quadric5::Minimum`（5×5 高斯消元，`quadric5.h:554-591, 453-515`）求使联合误差最小的 5D 落点；奇异时退回端点/中点择优（`ComputeMinimal`，`...tex.h:468-516`）。
3. `ComputeTexPriority`（`...tex.h:305-382`）**临时把两端移动到候选落点**，扫描折叠后存活的面：
   - `QuadErr = qsum.Apply(vv)`：5D 联合误差（`...tex.h:326`），下限钳到 1e-15（`:367`）；
   - 计算最坏面形状质量 `MinQual = min QualityFace(...)`（`QualityFace` 为内切/外接圆半径比，退化面→0，正三角形→1，见 `triangle3.h:186-191`）；
   - 最终 `priority = QuadErr / MinQual`（`...tex.h:369`）——会产生瘦长/退化三角形的折叠被**软性惩罚**（优先级升高、排到堆后），但**不是硬性拒绝**；
   - 可选 `NormalCheck`：法向变化超过 `CosineThr` 时优先级再 ×1000（`...tex.h:372-374`），默认关闭。
4. 接缝边（ncoords==2）：两组分别求优先级，取 `max(priority1, priority2)`；并对较差的一组调用 `ComputeMinimalWithGeoContraints`（`quadric5.h:519-551`），**锁定 xyz 等于另一组的几何落点、只优化 uv**，保证两套 UV 解共享同一 3D 位置（`...tex.h:451-463`，`Execute` 注释 `:713-714`）。
5. `ExtraTCoordWeight`（默认 **0.0**）按两端 UV 槽数量对多 UV 顶点的折叠再加惩罚（`...tex.h:433, 463`）；默认值意味着该补偿机制出厂时不生效。

### 2.4 执行：一次折叠保留两套新 UV

`Execute`（`...tex.h:687-794`）：

- 几何上交给 `EdgeCollapser::Do`（`edge_collapse.h:212-325`）：删除含该边的面、把 v0 的邻面重指到 v1、删除 v0，v1 取新 3D 位置；
- 连续区域只写一套新 UV；接缝边写两套新 UV（`newtcoord1/2`，`:725-742`）；
- 遍历存活顶点 v1 的所有面，按**逐分量精确相等的 float 比较**把旧 UV 重写到对应新 UV（`:745-768`）；匹配不上既不属于两套 UV 的楔 UV 走 `else` 分支，按其 UV 槽合并二次型后追加（`:769-787`）。

---

## 3. 纹理接缝与拓扑有效性如何限制候选选择

### 3.1 拓扑有效性：Link 条件始终生效（注意与几何版的差异）

- 带纹理算子**没有重写** `IsFeasible`（对该头文件 `grep IsFeasible` 无命中），因此继承基类 `TriEdgeCollapse::IsFeasible`（`tri_edge_collapse.h:194-196`），**无条件**调用 `EdgeCollapser::LinkConditions`。
- 对照：纯几何版 `TriEdgeCollapseQuadric::IsFeasible`（`tri_edge_collapse_quadric.h:149-156`）只有在 `pp->PreserveTopology == true` 时才检查 Link 条件。
- 后果：带纹理算子的参数类里虽然也有 `PreserveTopology`（`...tex.h:50,74`，默认 false），但该标志**对可行性判断完全不起作用**——拓扑保持是硬性、不可关闭的。会改变连通分支数、产生非流形或自交的折叠在 `DoOptimization` 中直接不执行（`local_optimization.h:219-227`）。
- Link 条件本身（`edge_collapse.h:133-207`，依据 Dey–Edelsbrunner 1999）通过对两端邻域的顶点/边计数实现：共享边计数必须为 0、共享顶点计数必须等于边的 link 顶点数；并对边界端点注入"虚拟顶点/虚拟边"以正确处理开网格边界。

### 3.2 接缝的限制方式：不阻止几何合并，但强制 UV 不被焊接

- Link 条件只看拓扑，不看 UV，因此**接缝边本身可以折叠**。
- 接缝边两侧 UV 不同 → `GetTexCoords` 返回 2 → 两侧各自维护独立 5D 二次型与独立新 UV，优先级取两侧最大值（§2.3 第 4 点）。几何约束保证折叠后仍是一个 3D 顶点，但携带两套 (u,v)，楔坐标按来源分别重写（§2.4）。**即接缝在几何上可以被压平，但纹理坐标不会被平均成一个值**，避免把两张纹理图集的采样位置焊死在一起。
- `ExtraTCoordWeight > 0` 时，UV 槽越多的顶点（多重接缝交汇处）折叠优先级越高、越晚被折叠，从而间接保护接缝（`...tex.h:433,463`）；默认 0 时此保护关闭。
- `PreserveBoundary=true` 时，边界边端点被清除可写标志而不进堆（`...tex.h:601-612`），结束后在 `Finalize` 恢复（`:213-225`）。这同时保护几何边界与 UV 图集裁切边（二者在数据中通常重合）。

### 3.3 其他限制

- 过期候选：端点被删除或 mark 落后于 `GlobalMark` 的堆元素被 `IsUpToDate` 丢弃（`tri_edge_collapse.h:198-210`），失败计数进 `FailStat::OutOfDate`。
- 可写/可读标志：初始化建堆只收 `IsRW()` 顶点（`...tex.h:617-637`）。
- **不存在**纯几何版中的硬性 `AreaCheck`、`HardQualityCheck`、`HardNormalCheck`（这些是 `TriEdgeCollapseQuadricParameter` 的成员，`tri_edge_collapse_quadric.h:87-91, 413-421`）；带纹理参数类（`...tex.h:40-84`）没有对应字段。纹理版对面形状只有 `QuadErr/MinQual` 软惩罚。

---

## 4. 三类场景下应出现的结果

### 4.1 正常网格（流形、无接缝、无退化）

**预期与当前行为**

- 每个顶点只有一个 UV 槽，`GetTexCoords` 返回 1，走 `...tex.h:432-433` 的单组路径。
- 5D 二次型同时度量几何位置与 UV 漂移；堆按误差/质量比从低到高执行折叠，等价于在 5D 空间做标准 Garland–Heckbert QEM，并带面形状软惩罚。
- Link 条件始终通过（流形内部边的 link 计数天然满足），折叠只受误差排序控制。
- 边界边通过 `byFace` 的边界二次型（权重 0.5）被隐式保留；显式保留需开 `PreserveBoundary`。
- 目标面数由 `LocalOptimization::SetTargetSimplices` 控制，堆空或达标即停（`texdeci.cpp:98-108`）。

### 4.2 接缝邻域（同一 3D 顶点携带多个不同 UV）

**预期与当前行为**

- 顶点的 `Vd()` 含多个 `(uv, Quadric5)` 槽（`...tex.h:96, 164`）；接缝边 `ncoords==2`，产生两套候选解，优先级取 `max`（`...tex.h:463`），即**接缝处的折叠代价由误差更大的那一侧决定**，接缝边相对连续区域被推迟。
- 两套解通过 `MinimumWithGeoContraints` 共享 xyz、各自优化 uv（`...tex.h:453-459`）；`Execute` 后存活顶点保留两套新 UV 并按楔来源精确匹配重写（`...tex.h:734-768`）。接缝的"一个几何顶点、多个纹理坐标"语义在简化后保持。
- 交汇于同一顶点但不参与当前边的第三套 UV：在 `Execute` 的 `else` 分支中原样保留并合并二次型（`...tex.h:769-787`）。

**当前实现的风险点（静态阅读结论）**

1. **精确 float 相等匹配**：UV 重写（`:750-764`）、`Contains`/`Qd`（`:133-160`）、`GetTexCoords` 同一性判定（`:262-267`）全部用 `==` 比较 `float`。数值上应相同但经不同浮点路径产生的 UV 会被当作不同槽，可能累积出冗余 UV 槽或落入 `Qd` 找不到坐标时的 `assert(0)`（`:158`）/ `Execute` 中的 `assert(0)`（`:784`，仅 debug 生效；release 下行为未定义）。
2. **非流形边隐患**：`GetTexCoords` 只记录前两个共面邻面。流形网格一条边最多两面，没有问题；但若输入含 >2 面的非流形边，第三面的 UV 不参与优先级计算，且 Link 条件通常先把该边判为不可行而掩盖此问题——两个机制的交互没有测试覆盖。
3. `Sum3` 对 UV 连续性是"伪造"的（`quadric5.h:655-664` 注释自陈），默认 `ExtraTCoordWeight=0` 时多 UV 顶点的优先级偏差没有任何补偿。

### 4.3 退化三角形（零面积 / 重复顶点引用）

区分两种"退化"（源码在 `clean.h:470-475` 明确区分）：

**(a) 拓扑退化面：两个顶点引用指向同一顶点**

- `Clean::RemoveDegenerateFace` 直接删除（`clean.h:477-493`）；`RemoveDuplicateVertex(..., RemoveDegenerateFlag=true)` 也会级联清理（`clean.h:206, 270-272`）。
- `texdeci.cpp:73-77` 在简化前默认执行重复顶点/未引用顶点清理。**建议把清理作为简化前置条件**——注释同时指出拓扑退化面会让拓扑计算函数崩溃（`clean.h:474-476`），而 Link 条件与 VF 迭代都依赖完整拓扑。

**(b) 几何零面积但三顶点引用不同（三点共线/重合坐标）**

- `QualityFace` 返回接近 0（`triangle3.h:196-210` 中 `area2<=0 → 0`）。
- 5D 二次型构建分两支（`quadric5.h:222-270`）：
  - 几何质量 `q==0` 但 **UV 面积非零**（2D 叉积判别，`:260-263`）：仍构建完整 5D 二次型，即退化几何面仍可贡献纹理误差；
  - UV 面积也为零：整个二次型置零（`:265-269`），该面对误差无贡献，折叠删除它时几乎无 QEM 代价。
- 折叠模拟中若产生零面积存活面，`MinQual→0`，使 `QuadErr/MinQual` 变得很大（`...tex.h:366-369`；`QuadErr` 被钳到 1e-15 下限，避免了除零但优先级数值巨大），该折叠被排到堆尾——**是软抑制而非硬拒绝**；当堆中没有其他候选时仍可能被执行，算法本身不保证输出零退化。
- 对照：纯几何版可通过 `HardQualityCheck`/`AreaCheck` 把这类折叠优先级直接置为 `max()` 硬拒绝（`tri_edge_collapse_quadric.h:413-418`），纹理版没有等价开关。

---

## 5. 测试与构建证据

**当前测试覆盖（基于仓库文件系统静态核查）**

- 仓库没有针对边折叠/QEM 的单元测试：`apps/test/` 下只有 `camerashot`、`quaternion`、`trackball`、`tetramesh` 四个无关测试程序；全仓 `grep TriEdgeCollapse` 的命中仅 `tridecimator`（纯几何）与 `texdeci`（纹理）两个 app。
- 唯一的纹理简化程序 `apps/unsupported/texdeci/texdeci.cpp` 位于 `unsupported/`，该目录**不在 `apps/CMakeLists.txt:5-8` 的 `add_subdirectory` 列表中**，只有 qmake 的 `.pro` 文件——默认 CMake 构建既不编译也不运行它，CI 层面没有任何可执行证据覆盖 `tri_edge_collapse_quadric_tex.h`。
- 可引用的"证据"因此只有：(i) 源码内部断言（如 `...tex.h:423-424, 727-728` 的 `IsValid()`/非 NaN 检查，`quadric5.h:510-511`）；(ii) `FailStat` 失败统计框架（`tri_edge_collapse.h:50-67`，但纹理版因继承基类 `IsFeasible` 不会写入其中的 Link 条件计数器——计数只在几何版 `tri_edge_collapse_quadric.h:153-154` 自增）；(iii) 纯几何 `tridecimator` 对同一 `EdgeCollapser::Do`/`LinkConditions` 底层的间接使用。

---

## 6. 当前行为 vs. 建议

| # | 主题 | 当前行为（源码证据） | 建议 |
|---|---|---|---|
| 1 | 拓扑保持开关 | 纹理版继承基类 `IsFeasible`，Link 条件**始终**生效，参数 `PreserveTopology` 被忽略（`tri_edge_collapse.h:194-196`；`...tex.h` 无 `IsFeasible` 重写） | 若要与几何版语义对齐，应像 `tri_edge_collapse_quadric.h:149-156` 一样按 `PreserveTopology` 条件化；或在文档/参数注释中明确"纹理版恒保拓扑"。当前属于易误用点 |
| 2 | 退化/翻转面 | 仅 `QuadErr/MinQual` 软惩罚，无硬拒绝（`...tex.h:366-369`）；无 `AreaCheck`/`HardQualityCheck` 等价物 | 把几何版的硬性面积/质量/翻转检查移植到纹理参数类；简化前强制 `Clean::RemoveDegenerateFace` + `RemoveDuplicateVertex` |
| 3 | UV 匹配 | 全部用 float `==` 精确比较（`...tex.h:133-145, 262-267, 750-764`），失配触发 `assert(0)`（`:158, :784`） | 改为带容差的 UV 比较（或在导入时对 UV 做量化/焊接预处理），并以可恢复的失败分支取代裸 `assert(0)` |
| 4 | 接缝权重 | `ExtraTCoordWeight` 默认 0.0（`...tex.h:69`），`Sum3` 自陈优先级估算偏差无补偿（`quadric5.h:655-664`） | 在接缝敏感场景默认给一个小的正权重，或在文档中给出取值指导 |
| 5 | 法向检查 | `NormalCheck` 默认关；开启时使用面**存储**法向 `x.F()->N()` 而非现算法向（`...tex.h:344-345, 357`），未初始化法向时会得到 NaN/垃圾值（几何版改用 `NormalizedTriangleNormal` 现算，`tri_edge_collapse_quadric.h:351-358`） | 与几何版统一为现算法向，或在 `Init` 中强制 `UpdateNormal` |
| 6 | 非流形边 | `GetTexCoords` 只取前两面（`...tex.h:249-273`），>2 面边的第三面 UV 不进优先级 | 显式检测非流形边并跳过/报统计，不要依赖 Link 条件隐式掩盖 |
| 7 | 测试与构建 | `texdeci` 不进 CMake、无任何折叠单元测试（`apps/CMakeLists.txt:5-8`；`apps/test/` 无相关用例） | 为三类场景（正常/接缝/退化）添加最小合成网格的单元测试：断言折叠后面数、UV 槽数、接缝两侧 UV 不被焊接、Link 拒绝计数；把 `texdeci` 纳入构建或迁移出 unsupported |
| 8 | 失败可观测性 | 纹理版不更新 `FailStat::LinkCondition*` 计数（计数点仅存在于几何版） | 在纹理版 `IsFeasible` 中接入同一统计框架，便于量化"多少候选因拓扑被拒" |

### 三类场景的结论速览

- **正常网格**：行为即 5D 空间标准 QEM + 面形状软惩罚 + 恒开 Link 条件，预期可稳定达到目标面数。
- **接缝邻域**：几何可折叠、UV 不焊接，双解共享 xyz、优先级取两侧最大；但 float 精确匹配与默认关闭的接缝权重是主要脆弱点。
- **退化三角形**：拓扑退化依赖简化前 `Clean` 清理（否则拓扑计算可能崩溃）；零面积面在误差侧被软抑制，算法不硬性保证输出无退化，需要前置清理或移植硬检查。

---

## 附：关键代码位置索引

- 每顶点多 UV 二次型容器：`tri_edge_collapse_quadric_tex.h:90-172`
- 候选 UV 收集（至多两面）：`tri_edge_collapse_quadric_tex.h:239-276`
- 双解优先级与几何约束：`tri_edge_collapse_quadric_tex.h:384-466`
- 质量软惩罚与可选法向惩罚：`tri_edge_collapse_quadric_tex.h:305-382`
- 折叠执行与 UV 重写：`tri_edge_collapse_quadric_tex.h:687-794`
- Link 条件（拓扑有效性）：`edge_collapse.h:133-207`
- 基类恒开可行性检查：`tri_edge_collapse.h:194-196`
- 几何版条件化可行性（对照）：`tri_edge_collapse_quadric.h:149-156`
- 5D 二次型构建/退化分支：`quadric5.h:218-271`
- `Sum3` 对接缝连续性的近似及注释：`quadric5.h:655-690`
- 几何约束最小化：`quadric5.h:519-551`
- 退化面/重复顶点清理：`clean.h:206-272, 470-493`
- 面质量度量：`triangle3.h:186-210`
- 唯一纹理简化示例（未纳入 CMake）：`apps/unsupported/texdeci/texdeci.cpp:68-115`
