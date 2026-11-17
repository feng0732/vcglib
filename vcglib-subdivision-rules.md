# VCGLib 细分的内部面 / 边界 / 锐边规则与属性传递分析

本文只基于仓库现有实现（分支 `task/21`，HEAD `5cae2ab`）进行静态阅读与外部实测，
**未修改仓库内任何代码**。所有实测均在仓库外的临时目录 `/tmp/vcgprobe/` 中用独立小程序
`g++ -std=c++17 -I<repo> -I<repo>/eigenlib` 编译运行，用于验证源码读出的权重与标志传播。
行号均指当前 HEAD。

---

## 0. 涉及的两套细分设施

VCGLib 里"混合三角形与多边形"的细分实际由**两套互不相同的设施**承担：

| 家族 | 入口 | 输入 | 细分粒度 | 明确支持边界规则 | 认识"锐边" |
|---|---|---|---|---|---|
| 三角网格 | `tri::RefineE` / `tri::Refine` / `tri::RefineOddEvenE` | 纯三角形 + `FFAdjacency` | 单条边 | 是（谓词 + 各 `MIDPOINT` 仿函数） | 只传标志，不改几何 |
| 多边形网格 | `tri::CatmullClark<>::Refine` / `tri::DooSabin<>::Refine` | 带 `face::PolyInfo` 的多边形 | 整张面 | **否**（无任何边界分支） | **否**（全文件无 crease 相关代码） |
| 桥接 | `PolygonSupport::MergeFlatFaces` + `ImportFromTriMesh`；`BitQuadCreation::MakePureByRefine` | 三角形 → faux 边编码的多边形 / bit-quad | — | 按 faux 边 | 部分（有 TODO，见 §4.4） |

多边形在三角网格里的表示方式是"把多边形内部的对角边标成 **faux**"，见
`vcg/complex/algorithms/refine_catmullclark.h:34-38`：

> In a trimesh, the generic polygons with n sides are codified represented by
> tagging the internal edge of the face as 'faux' with the SetF().

因此下文中的"内部面 / 边界 / 锐边"在三角家族里分别对应
**内部边、`BORDER0..2` 位、`FACEEDGESEL0..2` 位**（`vcg/simplex/face/base.h:114-123`），
在多边形家族里只有"内部 / 边界"之分（锐边概念缺失）。

---

## 1. 三角家族的结构：模板选择与三种拓扑情形的关系

### 1.1 模板只按"哪几条边被切"选取

`RefineE` 的细分模板是 8 行查表 `SplitTab`（`refine.h:95-106`），索引由三条边是否
被判据命中决定：

```cpp
int ind = ((vv[3] != NULL) ? 1 : 0) + ((vv[4] != NULL) ? 2 : 0) + ((vv[5] != NULL) ? 4 : 0);
```

（`refine.h:457`）。**模板本身完全不知道边界和锐边**——边界/锐边只通过下面两条路径
影响结果：

1. **是否切这条边**：由 `EDGEPRED` 决定（`refine.h:364`, `refine.h:402`）；
2. **切了以后新面拿到什么标志**：`SplitTab[ind].TE` 把"新三角形的第 j 条边"映射回
   "原三角形的第几条边"（`refine.h:487-504`），据此复制 `B`（边界）与
   `FACEEDGESEL`（锐边）位。`TE == 3` 表示新边是内部生成的对角线，强制清零。

两种"两刀"情形还需要选对角线，规则是**比较两条候选对角线的长度，取短者**
（`refine.h:506-511`）：

```cpp
if(SplitTab[ind].TriNum==3 &&
   SquaredDistance(vv[...swap[0][0]]->P(), vv[...swap[0][1]]->P()) <
   SquaredDistance(vv[...swap[1][0]]->P(), vv[...swap[1][1]]->P()) )
{ (*nf[2]).V(1)=(*nf[1]).V(0); (*nf[1]).V(1)=(*nf[2]).V(0); ... }
```

注意对角线选择是**几何的**（依赖顶点位置），所以即便拓扑完全相同的两个三角形，
顶点位置不同会得到不同的三角剖分。实测（三个顶点 A=(0,0,0)、B=(1,0,0)、C=(x,y,0)，
只切 AB、BC 两条边）：

| C | 输出三张面 | 说明 |
|---|---|---|
| (0,1,0) 或 (-0.9,0.05,0) | `[3,1,4] [0,4,2] [4,0,3]` | `swap[0]` 更短，交换发生 |
| (0.5,0.4,0) | `[3,1,4] [0,3,2] [4,2,3]` | `swap[1]` 更短，不交换 |

（`/tmp/vcgprobe/p8.cpp`，顶点 3=AB 中点、4=BC 中点。）

### 1.2 边界边在 `RefineE` 中只被"标记一次"

```cpp
RD[edgeCur.F()].ep[edgeCur.E()]=true;  ++NewFaceNum; ++NewVertNum;
assert(edgeCur.IsManifold());
if(!edgeCur.IsBorder()) { edgeCur.FlipF(); edgeCur.F()->SetV(); RD[...].ep[...]=true; ++NewFaceNum; }
```

（`refine.h:366-376`）。即边界边只增加 1 张面、1 个顶点；内部边增加 2 张面。
第二遍填充 `vp[]` 时同样只在非边界时把中点指针写到对侧（`refine.h:402-413`）。
**这是 `RefineE` 里唯一一处真正的边界分支**，它只影响计数与指针共享，
不影响中点的位置公式——位置公式由具体的 `MIDPOINT` 仿函数决定。
`refine.h:281-283` 特别警告该仿函数必须关于 `Pos` 翻转对称，否则会崩。

---

## 2. 内部面

### 2.1 位置规则

| 仿函数 | 内部规则 | 源码 |
|---|---|---|
| `MidPoint` | `P = (P0+P1)/2`（先按指针排序保证对称） | `refine.h:143-145` |
| `MidPointButterfly` | `(vl+vr)/2 + (vu+vd)/8 - (vul+vur+vdl+vdr)/16` | `refine.h:627` |
| `OddPointLoop` | `3/8(l+r) + 1/8(u+d)`（两端点中若有价 6 顶点则改用 `weight` 表；注释承认"其余情形未处理"） | `refine_loop.h:386-399`、`:393` |
| `EvenPointLoop`（k 价内部点） | `(1 - k·β(k))·P + β(k)·Σ 邻点`，β 见 `LoopWeight::beta` | `refine_loop.h:491-497`、`:67-69` |

实测（3×3 三角网格，中心顶点 z 抬到 1，切全部边）：

* 内部边的 Loop 奇点 `(1,0.5,0.375)`：`3/8·1 = 0.375`，与式一致；
  若为 `MidPoint` 则为 `(1,0.5,0)`。
* 6 价内部顶点：`k=6, β=1/16, 1-kβ=0.625`，实测 `(1,1,0.625)`（`p28`）。
  若误用 8 环则应为 `0.679458`，实测值证明 `EvenPointLoopGeneric` 的
  `NextE()` 环计数走的就是价数（`refine_loop.h:452-457`）。

`MidPointButterfly` 的内部权重在 `refine.h:650-674` 还有一份按"邻点是否重合"分情况的
`#if 0` 旧实现，当前**不参与编译**；实际生效的是 `:627` 那条。

### 2.2 两种"两刀"与"三刀"情形

`SplitTab` 覆盖 8 种组合（`refine.h:97-105`）：1→1、1→2（三种）、1→3（三种）、1→4。
1→4 时中心三角形的顶点顺序是 `(3,4,5)`（三个中点），且 `TE` 全为 3
（`refine.h:105`），即**中心面的三条边都是"内部边"，`B` 与 `FACEEDGESEL` 一律清零**
（`refine.h:499-503`）。这是"内部面"最直观的拓扑判据：中心面永远不会被继承任何边界/锐边标记。

---

## 3. 边界

### 3.1 `MidPoint`：边界上没有任何特殊规则

`MidPoint::operator()`（`refine.h:139-160`）只做端点平均，因此边界边的中点就是几何中点。
实测边界边 (0,0,0)-(1,0,0) 的新点 `(0.5,0,0)`，与内部边使用同一公式。

### 3.2 Loop：边界边 = 1/2，边界点 = 3/4 + 1/8 + 1/8

奇点（新点）：

```cpp
if (he.IsBorder()) {
    proj.addVertex(*l, 0.5);
    proj.addVertex(*r, 0.5);
    ...
} else { ... 3/8,3/8,1/8,1/8 ... }
```

（`refine_loop.h:370-377`）。偶点（原顶点）：

```cpp
if (he.IsBorder()) {          // Border rule
    if(valence) { ... (*valence)[he.V()] = std::max(2*(k-1), 3); }
    he.FlipV(); r = he.v; he.FlipV(); he.NextB(); l = he.v;
    proj.addVertex(*curr, 3.0/4.0);
    proj.addVertex(*l, 1.0/8.0);
    proj.addVertex(*r, 1.0/8.0);
    proj.project(nv);
}
```

（`refine_loop.h:459-485`）。注意 `l`、`r` 由 `NextB()` 沿**边界**取得
（`pos.h:269-282`），所以偶点只受两个边界邻居影响，**完全忽略内部邻居**。

可观察差异（6×2 三角带，抬升边界点 x=2）：

| 顶点 | 输入 | 输出 | 说明 |
|---|---|---|---|
| 边界顶点 x=2（z=1，两边界邻居 z=0） | (2,0,1) | (2,0,0.75) | `3/4·1 + 1/8·0 + 1/8·0` |
| 同点，若内部邻居也抬到 1 | — | 仍为 0.75 | 内部邻居被边界规则完全忽略 |
| 边界边 (2,0)-(3,0) 奇点 | — | (2.5,0,0.5) | 纯 1/2 |
| 内部边 (2,1)-(3,1) 奇点（邻点全 0） | — | (2.5,1,0) | 1/8 项参与 |

（`p24.cpp`、`p11.cpp`。）

另外，边界顶点会被写入一个**伪造价数** `max(2*(k-1),3)`（`refine_loop.h:467-468`，
注释自认 "not perfect, but better than nothing"），这个价数随后被 `OddPointLoop`
用来选择 `weight`（`refine_loop.h:386-392`），从而间接影响与边界相邻的**内部边**的奇点位置。

`EvenPointLoopGeneric` 没有"角点(comer)规则"，`3/4+1/8+1/8` 对折角同样生效。
实测一个 90° 折线端点（同时是边界角点）`(0,0,0) → (0.125,0,0.125)`，
被拉离了折线（`p26.cpp`）。

### 3.3 Butterfly：9/16 − 1/16

```cpp
if(he.IsBorder()) {
    he.NextB(); vr0=&he.v->P(); he.FlipV(); he.NextB();
    assert(&he.v->P()==vl); he.NextB(); vl0=&he.v->P();
    nv.P()=((*vl)+(*vr))*(9.0/16.0)-((*vl0)+(*vr0))/16.0;
}
```

（`refine.h:595-605`）。实测抬升一个边界顶点后，相邻两条边界边上的新点均为
`z = 9/16·1 = 0.5625`（而非 0.5 的中点），内部边仍为 `1/2, 1/8, -1/16` 组合
（`p11.cpp`, `p27.cpp`）。`MidPointButterfly2` 在两端任一为边界时直接回退到
`MidPointButterfly`（`refine.h:742`）。

### 3.4 边界标志（`B`）如何随拓扑传播

只在 `TE != 3` 时按 `TE` 回填，其它一律清零：

```cpp
if(SplitTab[ind].TE[i][j]!=3) {
    if(orgflag & (MESH_TYPE::FaceType::BORDER0<<(SplitTab[ind].TE[i][j])))
        (*nf[i]).SetB(j); else (*nf[i]).ClearB(j);
    if(orgflag & (MESH_TYPE::FaceType::FACEEDGESEL0<<(SplitTab[ind].TE[i][j])))
        (*nf[i]).SetFaceEdgeS(j); else (*nf[i]).ClearFaceEdgeS(j);
} else { (*nf[i]).ClearB(j); (*nf[i]).ClearFaceEdgeS(j); }
```

（`refine.h:487-504`）。`orgflag` 在改写面之前捕获（`refine.h:477`）。
两刀交换对角线时，还需要把 `B`/`FACEEDGESEL` 从第 0 号楔搬到第 1 号楔
（`refine.h:517-525`），否则标记会错位。

实测：8 张面的 3×3 网格切全部边后 `B` 位与重算的 FF 完全一致
（`B-vs-FF mismatch = 0`），边界半边的总数 8 → 16，锐边标记 2 → 4
（`p3.cpp` / `p4.cpp`），说明边界/锐边标记在 1→4 细分下是**按边数线性扩张**的。

但有两处副作用需要注意：

1. **新顶点不带 `B`（顶点级边界位）**。`RefineE` 只维护面级 `B`；
   实测新建的 16 个顶点没有一个 `IsB()`，必须再调用
   `UpdateFlags::VertexBorderFromFaceBorder` 才能得到 16 个边界顶点（`p29.cpp`）。
2. **`V`（visited）位泄漏**。为去重，`RefineE` 在对侧面上 `SetV()`
   （`refine.h:373`），而从不清理；`ImportData` 又把整份 `Flags()` 拷贝给子面
   （`vcg/simplex/face/component.h:367-371`）。实测 8 张面细化后仍有 4 张 `IsV()`
   （`p20.cpp`）。若下游代码（如 `PolygonSupport::ExtractPolygon`）依赖 `V` 位，
   必须先 `FaceClearV`。

### 3.5 `RefineSelected` 造成的"选择边界"

`RefineSelected` 会跳过"对侧未选中"的边（`refine.h:363`, `refine.h:400`），
边界边因 `FFlip()` 指向自身而总被视为选中。实测：

* 只选中一个**内部**三角形 → **完全不切**（VN/FN 不变），
  因为它的三条边对侧都未选中（`p29.cpp`）；
* 只选中一个**含一条边界边**的三角形 → 仅切那条边界边，1→2，
  其余边在对侧留下 T 顶点（`p35.cpp`）。

`TrivialMidPointRefine` 的文档明确写了这个后果（`refine.h:1000-1002`：
"if only a subset of the faces are selected T-vertices will be created"），
而 `RefineE` 的结果亦然。

### 3.6 多边形家族：Catmull-Clark 的"边界规则"实际上是把边界点拉向内部

`refine_catmullclark.h` 全文没有任何边界判断。它只依靠一个**计数归一化**：

```cpp
// the only normalization step we do is for edge vertices to handle the fact
// that on boundary faces we do not have two contributions so the weight should be doubled.
for(auto &v : refinedOut.vert)
    if(v_cnH[v]>0) v.P() /= v_cnH[&v];
```

（`refine_catmullclark.h:136-139`）。但代码里写入的**不是半个和**，而是
`0.5·面重心 + 0.25·V0 + 0.25·V1`（`:129-131`，权重和已经是 1）。
于是：

* 内部边：两面各贡献一次，除以 2 → `0.25(F_i+F_j) + 0.25(V0+V1)`，
  恰好等于标准 CC 边点 `(V0+V1+F_i+F_j)/4`；
* 边界边：只贡献一次，除以 1 → `0.5F + 0.25V0 + 0.25V1`，**被面重心拉向内部**。

顶点位置更新（`:141-163`）用

```cpp
float val1 = (val -2.0)/(val*val);
float val2 = 1.0/(val*val);
...
vert[...].P() += f.V(i)->P() * val1;
vert[...].P() += center * val2;
vert[...].P() += edgep  * val2;   // 每个面只加"出边"的那一个边点
```

于是整体为

```
P' = (n-2)/n · P + ( ΣF_i + ΣE_i ) / n²
```

（`n` 为入射面数，`F_i` 为面重心，`E_i` 为对应该面的边点。）而教科书 CC 为

```
P' = [ F̄ + 2R̄ + (n-3)·P ] / n         (F̄ = ΣF/n, R̄ = Σ中点/n)
```

两者有两处实质差异：本实现中 `P` 的系数是 `(n-2)/n`（而非 `(n-3)/n`，多出 `1/n`），
且边点用的是 `0.5·F + 0.25·V0 + 0.25·V1`（`refine_catmullclark.h:129-131`）
而不是 `(V0+V1)/2`，所以边项不是标准式子的简单缩放。
用上式代回实测数据可精确复现两种情形：内部 4 价中心 `0.65625`、
边界 2 价顶点 `(1.125, 0.4375, 0)`（见下）。

可观察差异（3×3 四边形网格，中心顶点 z=1，其余 z=0，一步）：

| 量 | 本实现实测 | 标准 CC 应为 |
|---|---|---|
| 内部中心点 | `(1,1,0.65625)` | `(1,1,0.5625)` |
| 边(0,0)-(1,0) 的边点 | `(0.5,0.25,0)` | `(0.5,0,0)` |
| 边界顶点 (1,0) | `(1.125,0.4375,0)` | 保持在 `(1,0,0)`（或按三次 B 样条边界规则） |
| 边界行 max y | `0.4375`（输入全 0） | `0` |

（`p9.cpp`、`p23.cpp`、`p25.cpp`。）也就是说：**Catmull-Clark 在开边界上会
把整条边界向内收缩并改变其形状**，不能用于需要保持边界曲线的情形；
`RefineSingleStep` 的注释（`:136-138`）描述的意图与代码行为不符。

### 3.7 DooSabin：边界直接产生退化多边形

`DooSabin::Refine`（`refine_doosabin.h:48-160`）分三批建面：

1. 每个原面 → 一个面（角点取 `(V+重心)/2`，`:102`）；
2. 每个原顶点 → 一个面，面数 = 该顶点的**关联面数**（`degreeVH`，`:78-87`），
   用 `VFOrderedStarFF` 取环（`:118`），并断言环长等于度数（`:119`）；
3. 每条边 → 一个四边形，条件是 `fi->FFp(j)->IsV()`（`:136`）。

对边界：

* 边界顶点只被 1~2 个面共享 → 对应"顶点面"退化为 1 边形 / 2 边形；
* 边界边满足 `FFp(j) == fi` 且 `fi` 已被 `SetV()`，于是**仍然建四边形**，
  但对侧 `FlipF()` 原地反弹，四个角里只有两个不同顶点（重引用）。

实测 3×3 四边形网格一步 DooSabin 后（`p18.cpp`、`p32.cpp`）：

```
out FN=25: 17 个四边形 + 4 个 1 边形 + 4 个 2 边形
其中 8 个"四边形"含重复顶点引用，例如
  [13] (0.25,0.25) (0.75,0.25) (0.75,0.25) (0.25,0.25)
```

后果：`UpdateTopology::FaceFace` 会在这张输出网格上直接断言失败
（`update/topology.h:167` 的 `assert(v[0] != v[1])`），即**输出网格的拓扑
无法用标准流程重算**。且 `DooSabin::Refine` 不写 FF 邻接，返回时
全部 80 个半边的 `FFp` 是野指针（`p32.cpp` 计数 80/80 越界），
必须自行重建。

---

## 4. 锐边（crease）

### 4.1 仓库中的锐边表示

锐边 = `face::FACEEDGESEL0..2`（`vcg/simplex/face/base.h:120-123`），
由 `UpdateFlags::FaceEdgeSelSignedCrease` / `FaceEdgeSelCrease` /
`FaceEdgeSelBorder` 计算（`update/flag.h:472-524`），
按**有符号二面角**判定，并要求两侧标记互反
（`MeshAssert::ConsistentFaceEdgeSelection`，`mesh_assert.h:179-189`）。
`face::Pos::IsEdgeS()` 是读取入口（`pos.h:343`）；
`Pos::NextEdgeS()`（`pos.h:325-340`）可沿锐边链行走，但**在细分代码中没有被使用**。

### 4.2 三角家族：锐边只影响标志，不影响几何

`RefineE` 中与锐边有关的只有两处：初始化时刷新 `B` 位（`refine.h:337`），
以及 `TE` 回填时的 `FACEEDGESEL` 拷贝（`refine.h:494-497`、`:522-525`）。
**没有任何 `MIDPOINT` 仿函数读取 `IsEdgeS()`**——
`MidPoint`、`MidPointButterfly`、`OddPointLoop`、`EvenPointLoop` 的代码里
都没有 crease 分支（`refine.h:139-160`、`:595-628`、`refine_loop.h:358-407`、
`:444-501`）。

实测（把网格中一条内部边标成锐边，其余完全相同）：

* `MidPoint`：锐边中点与不标锐边时**逐位相同**，`(0.5,1,0.5)`（`p4.cpp` §1）；
* Loop：25 个顶点上 `max|P_crease − P_nocrease| = 0`（`p5.cpp` L2）——
  **锐边完全不改变 Loop 的位置**。

也就是说，**Loop/Butterfly 会把锐边当作普通内部边去平滑**，锐边被圆化：

* 90° 折痕上的顶点（价 6）一步 Loop 后偏离折痕线 `0.1768`
  （`(1,0,0) → (1.125,0.125,0.125)`），折痕上的奇点也全部离开折痕线，
  落在折痕线上的顶点数 = 0（`p26.cpp`）。

标志层面的传播是"一条锐边 → 两条锐边（每条子边）"：实测标 1 条内部锐边后，
新网格中 `FACEEDGESEL` 半边数 2 → 4，且两侧互反（`p3.cpp`、`p20.cpp`），
`MeshAssert::ConsistentFaceEdgeSelection` 通过。但**几何上已经不再锐**。

### 4.3 唯一能"保留"锐边的路径：先切断（`CreaseCut`）

```cpp
template<class MESH_TYPE>
void CreaseCut(MESH_TYPE &m, float angleRad)
{
  tri::UpdateFlags<MESH_TYPE>::FaceEdgeSelSignedCrease(m, -angleRad, angleRad);
  CutMeshAlongSelectedFaceEdges(m);
}
```

（`crease_cut.h:36-40`）。`CutMeshAlongSelectedFaceEdges`（`:51-127`）沿标记边
**复制顶点**（`importData` 后改指针，`:120-125`），把锐边变成真正的边界边。

于是细分时走的是 §3.2/§3.3 的**边界规则**，曲线被保留：

| 方案 | 折痕顶点 `(1,0,0)` 一步后 | 在折痕线上的顶点数 | 最大偏离 |
|---|---|---|---|
| Loop，不切 | `(1.125,0.125,0.125)` | 0 | 0.1768 |
| Loop + `CreaseCut(30°)` | `(1,0,0)` | 7 / 31 | 0（仅折痕内部顶点） |

（`p25.cpp`、`p26.cpp`。）注意这不是"sharp subdivision 的 crease 规则"
（如 DeRose 的半锐/全锐规则），而是**用拓扑切割换取边界规则**：

* 锐边上的顶点被复制，网格不再水密（`p17.cpp`：VN 9 → 13）；
* 锐边**端点**（"...line" 的尽头）在切割后成为边界**角点**，而 Loop 没有角点规则，
  实测 `(0,0,0) → (0.125,0,0.125)`，端点仍被拉离折痕；
* 切割依赖 `FFAdjacency` 与流形性（`:57-60` 内部断言）。

### 4.4 多边形家族：没有锐边概念；桥接层会吞掉锐边

`refine_catmullclark.h` 与 `refine_doosabin.h` 全文不含 `crease` / `IsEdgeS` /
`IsF` 的判定（普查确认）。三角形转多边形时：

```cpp
// polygon_support.h:61-75
for (int w=0; w<3; w++) {
    FaceType *fb = fa->FFp(w);
    if ( (fb>fa) && (fa->N()*fb->N() > minDist) ) { fa->SetF(w); fb->SetF(fa->FFi(w)); }
}
```

`MergeFlatFaces` **只看二面角/法线**，没有检查 `IsFaceEdgeS`；
`ImportFromTriMesh`（`:149-176`）只搬运顶点属性与面的颜色/质量（`:172-174`），
**完全不搬运 `FACEEDGESEL`**。实测：标了 6 条锐边（其中含一条对角线）的网格
经 `MergeFlatFaces` + `ImportFromTriMesh` 后，多边形网格里
`FACEEDGESEL` 半边数 = 0（`p33.cpp`）；而锐边所在的面还被合并成一个 6 边形。
即 **锐边在 tri → poly 转换中被静默丢弃，同时多边形的边界因此改变**。

bit-quad 路线（`BitQuadCreation`）略有保留，但也有隐患：

* `MakePureByRefine` 通过"把边标成 faux"把三角形拼成四边形
  （`bitquad_creation.h:407-693`），面级标志（含 `FACEEDGESEL`）不做改动，
  实测锐边标记被保留（半边数 2 → 4，`p22.cpp`）；
* 但 `selectBestDiag` 里只有一句 TODO：

  ```cpp
  // todo: check creases? (continue if edge k is a crease)
  ```

  （`bitquad_creation.h:121`）。实测把一条内部锐边标出后跑一遍
  `MakePureByRefine`，该边**被标成 faux**（`f0 E1 faux=1 crease=1`），
  即"锐边"与"被合并掉的多边形内部边"重合，后续任何按 faux 还原多边形的
  工序都会把锐边当成内部边丢掉（`p22.cpp`）。
* `MakePureByCatmullClark` 就是 `MakePureByRefine` 调两次（`:689-693`）。

`Clean::HasConsistentPerFaceFauxFlag`（`clean.h:1031-1049`）只检查
"faux 互反且不在边界上"，**不检查 faux 与 crease 是否冲突**，所以上述
"faux 且 crease"的状态能通过一致性检查。

### 4.5 三条边界/锐边路径对比（可观察差异一览）

| 情形 | 三角家族（`RefineE`） | Catmull-Clark | DooSabin |
|---|---|---|---|
| 内部边/面 | 1→4（或按谓词 1→2/1→3），模板见 `SplitTab` | 每面一个面点，每边一个边点，全部四边形 | 每面/每顶点/每边各建面 |
| 边界 | `B` 位按 `TE` 继承；中点公式各有专门规则（Loop 1/2 与 3/4+1/8；Butterfly 9/16−1/16；MidPoint 无） | 无规则；边界点被拉向内部（y: 0 → 0.25/0.4375） | 边界顶点 → 1/2 边形；边界边 → 重引用四边形 |
| 锐边 | 只传标志，几何被平滑（偏离 0.1768）；需 `CreaseCut` 才保留 | 不存在 | 不存在 |

---

## 5. 属性如何随拓扑变化传递

### 5.1 `RefineE`（三角家族）— 逐项

新顶点由 `MIDPOINT` 自己写（`refine.h:139-160`），或由用户提供的
`InterpolatorFunctorType` 追加（`refine.h:158-159`，
模板参数默认 `BaseInterpolator`，`refine.h:109-115`）：

| 属性 | 规则 | 源码 |
|---|---|---|
| `P` | 各方案自己的位置规则 | 见 §2.1 |
| `N` | `(N0+N1).normalized()` | `refine.h:147-148` |
| `C`（顶点色） | `lerp(C0,C1,0.5)` | `refine.h:150-151` |
| `Q`（顶点质量） | `(Q0+Q1)/2` | `refine.h:153-154` |
| `T`（顶点纹理坐标） | `(T0+T1)/2` | `refine.h:156-157` |
| 用户自定义 | `(*intFunc)(nv, ep)` 可任意覆写 | `refine.h:158-159` |

实测（v0: C=(10,0,0) Q=0；v1: C=(0,20,0) Q=1）：
新点 `v9` 为 `Q=0.5`、`C=(10,10,50)`（`p4.cpp`），与上式一致。

新面：

```cpp
nf[0]=&*fi;                            // 原面本身复用
for(i=1;i<SplitTab[ind].TriNum;++i){
    nf[i]=&*lastf; ++lastf;
    if(RefineSelected || (*fi).IsS()) (*nf[i]).SetS();
    nf[i]->ImportData(*fi);            // 逐字复制全部面属性
}
```

（`refine.h:459-467`）。`ImportData` 会复制面法线、面颜色、面质量、
**以及整份 `Flags()`**（`vcg/simplex/face/component.h:253-257`、`:367-371`、
`:390-392`）。实测子面的面质量/颜色与父面完全相同（`Q=100/101`、`C=(9,8,7)`），
即使几何已被细分——**这是"复制"而非"重算"**。

楔属性（per-wedge）只对纹理坐标做了重映射：

```cpp
wtt[i]=(*fi).WT(i);
wtt[3+i]=mid.WedgeInterp((*fi).WT(i),(*fi).WT((i+1)%3));
...
(*nf[i]).WT(j) = wtt[SplitTab[ind].TV[i][j]];
```

（`refine.h:470-475`、`:483-484`）。而 `face::WedgeColor` / `WedgeNormal`
的 `ImportData` 是**按原楔序号 1:1 拷贝**（`component.h:411-419`、`:275`），
不做 `TV` 重映射。实测：子面 `WT` 已被重排（`u=0.5,1.5,1`），
但 `WC(0,1,2)` 仍是父面的 `(0,10,20)`，且新建的非首个子面直接是默认值
（`p21.cpp`）——**楔色/楔法线与顶点的对应关系在细分后错位**。

动态属性（`AddPerVertexAttribute` / `AddPerFaceAttribute`）**不会被插值**：
`AddVertices`/`AddFaces` 只把 `SimpleTempData` 扩容
（`allocate.h:189-208`），新元素值保持默认构造值。实测新顶点的
`float` 属性全为 0（`p11.cpp`、`p12.cpp`）。

还有一个容易被忽略的点：`ImportData` 复制的是**整份标志位**，其中包含 `F`（faux）。
`RefineE` 之后并没有按 `TE` 重映射 faux 位，于是 1→4 细分会产生
**faux 互反性破坏**：

```
pre  Clean::IsFaceFauxConsistent = 1
post Clean::IsFaceFauxConsistent = 0，8 个半边与其邻居的 faux 位不一致
```

（`p37.cpp`；判定函数见 `clean.h:877-894`。）如果输入是 bit-quad 表示的多边形，
细分后必须重新整理 faux 位，否则后续 `PolygonSupport` 的合并会出错。

### 5.2 Loop（`RefineOddEvenE`）：偶点/奇点的属性来源不同

`RefineOddEvenE`（`refine_loop.h:553-621`）先用原位置的**偶点规则**把结果存进
`newEven`，再调用 `RefineE` 生成奇点，最后写回偶点：

```cpp
std::vector<std::pair<CoordType,CoordType> > newEven(m.vn);
...
even(newEven[index], aux);
...
RefineE< MESH_TYPE, ODD_VERT > (m, odd, edgePred, RefineSelected, cbOdd);
for(i...) { m.vert[i].P()=newEven[i].first; m.vert[i].N()=newEven[i].second; }
```

（`refine_loop.h:577`、`:599`、`:607`、`:608-613`。）

因此：

* **偶点只更新 `P` 和 `N`**，`C`/`Q`/`T` 保持原值（实测 v0..v8 的 `Q`、`T` 未变，`p14.cpp`）；
* **奇点只由 `OddPointLoopGeneric` 写**，其中仅 `P`、`N`、`C` 被赋值
  （`refine_loop.h:367-368`、`:375-376`），`Q`/`T` 保持新增元素的默认值
  （实测奇点 `Q=0`、`T=(0,0)`，`p14.cpp`）；
* **偶点的颜色规则是"与第一个入射面的下一个顶点取平均"**：

  ```cpp
  (*fi).V(i)->C().lerp((*fi).V0(i)->C() , (*fi).V1(i)->C(),0.5f);
  ```

  （`refine_loop.h:589-591`）。其中 `V0(i)==V(i)`（`component.h` 的多边形访问器
  `component_polygon.h:82-93`），所以这是 `(C_self + C_next)/2`，
  而 `V(i+1)` 取决于遍历到该顶点时**第一个遇到的面**。实测同一网格仅改变
  面在数组中的顺序，v0 的颜色从 `(5,10,0)` 变成 `(25,25,0)`（`p14.cpp` / `p14b.cpp`）。
  这个改动发生在 `RefineE` 之前，所以**奇点的颜色平均读到的是已经被改过的偶点颜色**
  （`p15.cpp` 打印的端点颜色是 `(5,10,0)` 而不是原始的 `(10,0,0)`）。
* **法线在默认 `Centroid` 投影下是未定义值**。`Centroid::project` 只写 `nv.first`：

  ```cpp
  inline void project(std::pair<CoordType,CoordType> &nv) const {
      LVector position = sumP / sumW;
      nv.first = CoordType(...);        // nv.second 从未赋值
  }
  ```

  （`refine_loop.h:224-227`），而调用方会执行 `nv.N()=pp.second`
  （`:375-376`）或 `proj.project(nv)`（`:481`）。`Point3f()` 的默认构造
  不初始化成员（`space/point3.h:59`），所以 `N()` 得到垃圾值。
  实测 Loop 细化后法线为 `1e-41` 级或 `2.19e4` 级的无意义向量（`p13.cpp`）。
  只有改用 `LS3Projection`（`:240-340`，它会写 `nv.second`）才有合法法线。
* 价数属性是临时的：`RefineOddEvenE` 申请 `PerVertexAttribute<int>`、
  传给奇/偶两个仿函数、结束后删除（`refine_loop.h:570-572`、`:615-618`）。
  用户自定义的 per-vertex 属性同样不被插值（同 §5.1）。

### 5.3 Catmull-Clark / DooSabin（多边形家族）：只算位置，属性全部丢弃

两者都是**新建输出网格**（`refinedOut.Clear()`，`refine_catmullclark.h:67`；
`refine_doosabin.h:60`），只写坐标，从不 `ImportData`。实测：

| 量 | 输入 | CatmullClark 输出 | DooSabin 输出 |
|---|---|---|---|
| 顶点 `Q` | `1+x` | `0` | `0` |
| 顶点 `C` | `(1,2,3)` | `(255,255,255)`（白） | `(255,255,255)` |
| 顶点 `N` | `(0,0,1)` | `(0,0,0)` | `(0,0,0)` |
| 面 `Q`/`C` | `9` / `(4,5,6)` | `0` / 白 | `0` |
| 自定义 `PerVertexAttribute` | 有 | **不存在** | **不存在** |

（`p12.cpp`。）另外两者都**不写 FF 邻接**：`RefineSingleStep` 返回后
64 个半边的 `FFp` 全是越界值，`DooSabin` 80/80 越界；必须再调用
`UpdateTopology::FaceFace`（CC 可行，DooSabin 因退化面会断言失败，见 §3.7）。

CC 里也确实完全不使用属性：`PolyBarycenter` 只读坐标（`refine_catmullclark.h:97`），
边点只读坐标与计数（`:129-133`），顶点更新只读 `val`/坐标（`:148-162`）。

`DooSabin` 同样只写坐标：`(fi->V(j)->P()+b)/2.0f`（`refine_doosabin.h:102`），
其余全部是索引搬运。

### 5.4 桥接层对属性的处理

* `PolygonSupport::ImportFromTriMesh`：顶点走 `ImportData`（全属性），
  面只搬颜色与质量（`polygon_support.h:160-162`、`:172-174`），
  **不搬面法线、不搬 `FACEEDGESEL`、不搬结楔属性**（实测锐边丢失，`p33.cpp`）；
* `BitQuadCreation`：新顶点用 `ImportData(*(fi->V0(0)))` 整份复制（"lazy"，
  `bitquad_creation.h:470`、`:567`、`:648`），再用 `GeometricInterpolator::Apply`
  覆写位置（`bitquad_support.h:92-96`）；`GeometricInterpolator` **只插值 `P`**，
  另外在端点为边界时给新点 `SetB()`（`:95`）。自定义属性同样是"整份复制"
  而不是加权平均。

---

## 6. 结论与落到代码上的注意点

1. **模板与位置规则是正交的**。`SplitTab` 只描述"切了哪几条边"的三角剖分
   （`refine.h:95-106`、`:457`）；边界/锐边对几何的影响只来自
   `MIDPOINT` 仿函数里显式的 `IsBorder()` 分支——`MidPoint` 没有，
   `OddPointLoop`/`EvenPointLoop` 与 `MidPointButterfly` 有。
2. **锐边在三角家族里是"标志位"而非"几何约束"**。`RefineE` 只传播
   `FACEEDGESEL`（`refine.h:494`），Loop/Butterfly/MidPoint 都不读它，
   所以锐边会被平滑（实测偏离 0.1768）。要保留锐边只能 `CreaseCut`
   （`crease_cut.h:36-40`）把它变成边界，代价是顶点复制、网格不水密、
   且折痕端点（角点）仍无规则可用。
3. **多边形家族没有边界规则，也没有锐边概念**。Catmull-Clark 的
   "边界归一化"实际把边界点拉向面重心：边点 `(0.5,0.25,0)` 而非 `(0.5,0,0)`；
   边界顶点 y 从 0 变到 0.4375（`(1,0,0)→(1.125,0.4375,0)`）；
   顶点更新式为 `(n-2)/n·P + (ΣF+ΣE)/n²`，与教科书 CC 的
   `(F̄+2R̄+(n-3)P)/n` 不同（内部中心实测 `0.65625`，教科书式 `0.5625`）。
   DooSabin 在边界上产生 1/2 边形与重引用四边形，输出拓扑无法用
   `UpdateTopology::FaceFace` 重建。
4. **属性传递是"按组件类型分档"的，不是一致的**：
   * 三角家族：顶点 `P/N/C/Q/T` 有定义明确的插值；面属性逐字复制；
     楔纹理坐标重映射而楔色/楔法线不重映射；动态属性归零；
   * Loop：偶点只更新 `P/N`，奇点只更新 `P/N/C`；偶点颜色依赖面遍历顺序；
     默认 `Centroid` 投影下法线未定义；
   * 多边形家族：**全部丢失**；
   * faux 标志在 `RefineE` 后不再互反，需要重算。
5. **`RefineE` 的"隐式状态"**：新顶点不带 `B`；`V` 位残留在原面上；
   `RefineSelected` 会在未选中一侧留下 T 顶点。这些都不是 bug 报告，
   而是使用本库时必须自行补齐的后处理。

---

## 附录：仓库外验证程序

以下程序位于 `/tmp/vcgprobe/`（不属于仓库、未提交），均以
`g++ -std=c++17 -O0 -I<repo> -I<repo>/eigenlib` 编译。

| 文件 | 验证内容 | 关键输出 |
|---|---|---|
| `p3.cpp` | `B`/`FACEEDGESEL` 在 1→4 后的计数与互反 | 边界半边 8→16，锐边 2→4，mismatch=0 |
| `p4.cpp` | 锐边是否影响 `MidPoint` 位置 | 标/不标锐边，v12 均为 `(0.5,1,0.5)` |
| `p5.cpp` | 锐边是否影响 Loop 位置 | `max\|P_crease−P_nocrease\| = 0` |
| `p8.cpp` | 两刀情形的对角线选择 | 位置相关，两种面序 |
| `p9/p23/p25.cpp` | CC 权重与边界收缩 | 中心 `0.65625`；边界边点 `(0.5,0.25,0)` |
| `p11/p12.cpp` | 属性传递（三角/Loop/CC/DS） | 见 §5 |
| `p13.cpp` | Loop 默认投影下的法线 | 未初始化值（如 `2.19e4`） |
| `p14/p14b.cpp` | Loop 偶点颜色的面序依赖 | `(5,10,0)` vs `(25,25,0)` |
| `p15.cpp` | 奇点颜色读到已改写的偶点色 | 端点 `C=(5,10,0)` |
| `p17/p26.cpp` | `CreaseCut` 前后折痕保留 | 偏离 `0.1768` → 0 |
| `p18/p32.cpp` | DooSabin 退化面与 FF 野指针 | 8 个退化面；80/80 越界 |
| `p20/p29.cpp` | `V` 位泄漏、新点无 `B` | 4/8 面 `IsV`；新点 `B=0` |
| `p21.cpp` | 楔属性重映射不对称 | `WT` 重排、`WC` 原样 |
| `p22.cpp` | `MakePureByRefine` 把锐边变 faux | `faux=1 crease=1` |
| `p33.cpp` | tri→poly 丢失锐边 | `FACEEDGESEL` 半边 0 |
| `p37.cpp` | `RefineE` 后 faux 不再互反 | `IsFaceFauxConsistent` 1→0 |
