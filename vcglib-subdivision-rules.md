# VCGLib 混合三角形/多边形细分中的边界与锐边规则

分析对象：`/srv/fb2/feng0732/0021-vcglib-3-22d2a3b8`（分支 `task/21`，HEAD `5cae2ab`）。
**只读分析，未修改仓库任何文件。** 实验代码在 `/tmp/vcgtest/`（仓库外），下文所有数字均为实测输出。

---

## 0. 结论速览

| 问题 | 仓库实际行为（均已实测） |
|---|---|
| 边界边改变细分模板吗？ | **拓扑模板不变**（`SplitTab` 只按"哪几条边被分割"索引）；边界只改变**分割向邻面的传播次数**，并让 Loop/Butterfly 走**专门的取点公式**。 |
| 锐边（`FACEEDGESEL`）改变细分模板吗？ | **不。完全不改变。** 三条细分路径（`RefineE`/Loop/Catmull-Clark/Doo-Sabin）都不读该位。逐顶点比对：置位与不置位的 Loop/MidPoint/Butterfly 输出**位移恒为 0**。 |
| 锐边怎么才能生效？ | 必须先用 `tri::CreaseCut` / `CutMeshAlongSelectedFaceEdges` 把 crease **切开成边界**。但切开后跑 Loop 会**在折痕两端裂开**（实测缝隙 0.177–0.25），因为 Loop 只有"边界规则"、没有"折痕规则"。 |
| 多边形方案的边界规则？ | **不存在。** `CatmullClark` 对边界顶点会给出**跑到原网格包围盒之外**的位置；`DooSabin` 在边界上生成 1 边 / 2 边 / 重复顶点的**退化面**。两者输出都**不建 FF 邻接**。 |
| 属性传递？ | 三条互不相同的路径：① `RefineE` 走"插值 + `ImportData` + 位重映射"（较完整）；② Loop 依赖 `RefineE` 但**破坏原顶点颜色**并写入**未初始化的法向**（UB）；③ Catmull-Clark / Doo-Sabin **丢弃全部顶点与面属性**。 |

---

## 1. 实现地图与三种位的语义

| 文件 | 角色 |
|---|---|
| `vcg/complex/algorithms/refine.h` | 通用三角细分 `RefineE`/`Refine`/`RefineMidpoint`；`SplitTab[8]` 拓扑模板；`MidPoint`、`MidPointButterfly`、`MidPointButterfly2` |
| `vcg/complex/algorithms/refine_loop.h` | Loop：`OddPointLoopGeneric`（边中点）、`EvenPointLoopGeneric`（原顶点）、`RefineOddEvenE` 驱动 |
| `vcg/complex/algorithms/refine_catmullclark.h` | 多边形 Catmull-Clark |
| `vcg/complex/algorithms/refine_doosabin.h` | 多边形 Doo-Sabin |
| `vcg/complex/algorithms/polygon_support.h` | 三角↔多边形转换；`FAUX` 位编码多边形面 |
| `vcg/complex/algorithms/crease_cut.h` | `CreaseCut` / `CutMeshAlongSelectedFaceEdges`：把 crease 打开成边界 |
| `vcg/complex/algorithms/update/flag.h` | `FaceBorderFromFF`、`FaceEdgeSelCrease`、`FaceEdgeSelBorder` |
| `vcg/simplex/face/base.h` | 位定义；`pos.h`/`topology.h` 提供 `IsBorder()`、`NextB()` |

```cpp
// vcg/simplex/face/base.h:111,120,125
BORDER0      = 0x00000040,   // 该半边为流形边界
FACEEDGESEL0 = 0x00008000,   // 用户/角度标出的特征边（锐边）
FAUX0        = 0x00040000,   // 多边形面内部被三角化出来的假边
```

三者可同时置在同一条边上，互不推导。两个容易误解的点：

- **`IsBorder()` 不读 BORDER 位**，而是实时查 FF：`return f.cFFp(j)==&f;`（`vcg/simplex/face/topology.h:55-59`）。BORDER 位只是缓存。
- **`FaceEdgeSelSignedCrease` 默认不标边界**（`MarkBorderFlag=false`，`flag.h:472`）；要"边界也算特征边"须显式传 `true` 或调用 `FaceEdgeSelBorder`（`flag.h:499`）。

---

## 2. 拓扑模板：只看"哪些边被分割"

`RefineE` 用 `SplitTab[8]`（`refine.h:95-106`）按三位掩码选 1/2/3/4 个子三角形：

```cpp
// refine.h:457
int ind = ((vv[3] != NULL) ? 1 : 0) + ((vv[4] != NULL) ? 2 : 0) + ((vv[5] != NULL) ? 4 : 0);
```

- 三边全分 → 1→4。
- 两边分割 → 1→3，并按两条候选对角线长度的比较切换：

```cpp
// refine.h:506-511
if(SplitTab[ind].TriNum==3 &&
   SquaredDistance(vv[SplitTab[ind].swap[0][0]]->P(), vv[SplitTab[ind].swap[0][1]]->P()) <
   SquaredDistance(vv[SplitTab[ind].swap[1][0]]->P(), vv[SplitTab[ind].swap[1][1]]->P()) )
{ (*nf[2]).V(1)=(*nf[1]).V(0); (*nf[1]).V(1)=(*nf[2]).V(0); ... }
```

实测（`t_split.cpp`，三角形 (0,0),(3,0),(0,1)，只分割两条长边）：输出 3 三角形，实际选中 `v0–mp12`（长 1.5811），而非 `mp01–v2`（长 1.8028）——即**取较短对角线**。

**边界与锐边都不能改变 `ind`。** 它们只影响下面两件事。

### 2.1 边界 → 改变"分割传播的面数"

```cpp
// refine.h:366-376
RD[edgeCur.F()].ep[edgeCur.E()]=true;
++NewFaceNum; ++NewVertNum;
assert(edgeCur.IsManifold());
if(!edgeCur.IsBorder())
{
    edgeCur.FlipF();
    edgeCur.F()->SetV();
    RD[edgeCur.F()].ep[edgeCur.E()]=true;   // 只有非边界边才传播到邻面
    ++NewFaceNum;                            // 内部边 +2 面，边界边 +1 面
}
```

- 内部边被分割：1 新点、**2** 新面。
- 边界边被分割：1 新点、**1** 新面。
- 单三角形三边全分（全是边界边）仍得 1→4（每条边各自只分裂一次）：实测 `t_checks.cpp` → `ok=1 FN=4 VN=6`。

`RefineMidpoint`（`refine.h:1077`）用 `NextF()` 绕一圈标记同一几何边的所有面（`refine.h:1122-1132`），边界只绕一次，并用 `--NewFaceNum` 修正重复计数。

### 2.2 锐边 → 对拓扑**零影响**（源码 + 实测双证据）

`refine.h` 中 `IsFaceEdgeS`/`IsEdgeS` 只出现在"位搬运"处（`refine.h:487-503`、`522-525`、`1246-1261`、`1281-1284`），**没有一处用于选择权重或决定分割**。`refine_catmullclark.h`、`refine_doosabin.h` 中 `FaceEdgeS|IsF(|SetF(` **零命中**。

逐顶点实测（`t_ignore.cpp`，闭合八面体，一条内部边置 crease）：

```
Loop      VN 都不变=18；置 crease 后位移顶点数 = 0，最大位移 = 0.000000000000
MidPoint  位移顶点数 = 0，最大位移 = 0.000000000000
Butterfly 位移顶点数 = 0，最大位移 = 0.000000000000
```

### 2.3 faux 位：多边形在三角网格中的编码，且不会被重映射

`polygon_support.h` 用 `FAUX*` 编码多边形面（`MergeFlatFaces` 双向置位：`polygon_support.h:70-71`；`ImportFromPolyMesh`：`polygon_support.h:129-131`）。

`refine.h` 完全不看 faux 位。**但 `face::BitFlags::ImportData` 会整体复制 `Flags()`**（`vcg/simplex/face/component.h:367-371`），而 TE 重映射块**只覆盖 BORDER 与 FACEEDGESEL，不覆盖 FAUX**。结果是 faux 位被复制到所有子面、且边编号已变，语义错位。实测（`t_faux6.cpp`，共面对角线三角 + `MergeFlatFaces`，再 `Refine`）：

```
faux edges total: on-diagonal=4 off-diagonal=4
  faux edge NOT on diagonal: f0 e2 (0.50,0.50)-(0.50,0.00)
  faux edge NOT on diagonal: f1 e0 (0.50,0.50)-(0.50,1.00)
  ...
```

**后果**：把该网格送回 `PolygonSupport::ImportFromTriMesh` → `ExtractPolygon`（`polygon_support.h:192`），其 faux 游走（`polygon_support.h:222-238`）会**死循环**。实测 `t_faux4.cpp 1` 30s 超时（exit=124）；手工重放循环（`t_faux5.cpp`）3 步内进入无限 `FlipF`：

```
step 0: f2 e0 v0
  !! inner loop runaway
```

> **实践结论：不要对 faux 编码的三角网格直接调用 `tri::Refine`。**

---

## 3. 情形一：内部面

| 方案 | 新边中点 | 原顶点 |
|---|---|---|
| `MidPoint` | `(l+r)/2` | 不动 |
| Loop `OddPointLoopGeneric` | `3/8(l+r) + 1/8(u+d)` | `(1-kβ)·v + β·Σneighbors` |
| Butterfly | `(l+r)/2 + (u+d)/8 − (ul+ur+dl+dr)/16` | 不动 |
| Catmull-Clark | `0.25(a+b+C₁+C₂)` | 见 §5.1（与教科书权重不同） |
| Doo-Sabin | 角落中点 + 面心，建 1 面/角 + 1 面/边 | — |

Loop 奇点内部模板：

```cpp
// refine_loop.h:394-399（"unhandled case ... use Loop's weights" 分支）
proj.addVertex(*l, 3.0/8.0);
proj.addVertex(*r, 3.0/8.0);
proj.addVertex(*u, 1.0/8.0);
proj.addVertex(*d, 1.0/8.0);
```

偶点内部模板：

```cpp
// refine_loop.h:491-497
ScalarType beta = weight.beta(k);
proj.addVertex(*curr, 1.0 - ScalarType(k) * beta);
do { proj.addVertex(*he.VFlip(), beta); he.NextE(); } while(he != heStart);
```

`LoopWeight::beta(k)`：`k>3` 时用经典 `(5/8 − (3/8+cos(2π/k)/4)²)/k`，`k==3` 时 `3/16`（`refine_loop.h:67-69`）。
实测闭合八面体（全 valence 4）：`VN 6→18`，`v0 (1,0,0) → (0.5156,0,0)`；四面体 Loop：`VN 4→10`，`v0 (0,0,0) → (0.1875,0.1875,0.1875)`。

---

## 4. 情形二：边界

### 4.1 Loop 奇点（边中点）：二点平均

```cpp
// refine_loop.h:370-377
if (he.IsBorder()) {
    proj.addVertex(*l, 0.5);
    proj.addVertex(*r, 0.5);
    proj.project(pp);
    nv.P()=pp.first;
    nv.N()=pp.second;
}
```

实测（`t_mid.cpp`，开放正方形 (0,0),(1,0),(1,1,0.2),(0,1)）：

```
边界边中点: v4 (0.5000,0.0000,0.0000)   v5 (1.0000,0.5000,0.1000)
           v7 (0.5000,1.0000,0.1000)   v8 (0.0000,0.5000,0.0000)
内部对角线中点: v6 (0.5000,0.5000,0.0750)
```

`v4` 的 z 严格是 `(0+0)/2=0`，**对角顶点 z=0.2 没有参与** —— 边界规则阻断了跨边界影响（这是正确行为）。

### 4.2 Loop 偶点（原顶点）：`3/4, 1/8, 1/8`

```cpp
// refine_loop.h:459-485
if (he.IsBorder()) {   // Border rule
    if(valence) {
        k = 0; do { he.NextE(); k++; } while(!he.IsBorder());
        (*valence)[he.V()] = std::max(2*(k-1), 3);   // ← 覆写共享 valence 属性
    }
    he.FlipV(); r = he.v;
    he.FlipV(); he.NextB(); l = he.v;
    proj.addVertex(*curr, 3.0/4.0);
    proj.addVertex(*l, 1.0/8.0);
    proj.addVertex(*r, 1.0/8.0);
    proj.project(nv);
}
```

实测：`v0 (0,0,0) → (0.125,0.125,0)`、`v3 (0,1,0) → (0.125,0.875,0.025)`，与 `0.75·v + 0.125·l + 0.125·r` 逐位一致。

**两个必须注意的边界事实：**

1. **角点没有特殊处理。** valence=2 的角点也套 `3/4,1/8,1/8`，于是 (0,0,0) 内缩到 (0.125,0.125,0)。标准做法应固定角点，代码没有。
2. **边界分支会覆写 valence 属性**（`refine_loop.h:461-468`）。该值随后被**奇点规则**读取（`refine_loop.h:386-392` 的 `==6` 判据），于是边界顶点会改变邻近边中点走哪个权重分支。这是"边界改写锐边附近权重"的唯一通道，也是 §4.4 讨论的机制。

### 4.3 Butterfly 边界：4 点三点式，角点处会外溢

```cpp
// refine.h:595-605
if(he.IsBorder())
{
    he.NextB();              vr0 = &he.v->P();
    he.FlipV(); he.NextB();  assert(&he.v->P()==vl);
    he.NextB();              vl0 = &he.v->P();
    nv.P() = ((*vl)+(*vr))*(9.0/16.0) - ((*vl0)+(*vr0))/16.0;
}
```

实测（`t_bfly2.cpp`，边界 x 链 0,1,3,6）：
```
formula 9/16*(1+3) - 1/16*(0+6) = 1.8750
out new v13 (1.87500,0.00000,0.00000)          ← 精确吻合
```

**角点缺陷**：`NextB()` 会绕过角点取到另一侧的顶点。重放路径（`t_bfly4.cpp`）：

```
BORDER edge (3,0)-(6,0):
  r=(6,0)
  after NextB (get vr0): v=(6,1)     ← 绕过角点，取到 y=1 的顶点
  after NextB again (get vl0): v=(3,0)
```

代入公式得 `y = 9/16·(0+0) − 1/16·(0+1) = −0.0625`，即**新点跑到边界矩形之外**。实测 `t_bfly3.cpp`：4 个新顶点的 y 为 `−0.0625 / 1.0625`。这是"边界 4 点规则在角点处需替换为角点规则"的缺失，不是数值噪声。

### 4.4 边界 → 锐边的间接耦合（唯一通道）

Loop 的边界偶点分支把 valence 改写成 `max(2(k-1),3)`（`refine_loop.h:467`），而奇点规则用 `valence==6` 选择分支（`refine_loop.h:386-392`，否则落到 `3/8,3/8,1/8,1/8`）。因此"边界"会影响"邻近边中点"的权重选择；这只是副作用，代码里没有任何"折痕"概念。

---

## 5. 情形三：锐边

### 5.1 直接置 crease 位 → 完全无效

`t_ignore.cpp` 逐顶点比对证明位移恒为 0（见 §2.2）。源码层面：`refine.h` 只在位搬运处出现 `IsFaceEdgeS`，从不用于分支。

### 5.2 正确用法：先切开成边界

```cpp
// crease_cut.h:36-39
void CreaseCut(MESH_TYPE &m, float angleRad)
{
  tri::UpdateFlags<MESH_TYPE>::FaceEdgeSelSignedCrease(m, -angleRad, angleRad);
  CutMeshAlongSelectedFaceEdges(m);
}
```

`CutMeshAlongSelectedFaceEdges`（`crease_cut.h:51`）绕顶点星走一圈，每次跨过 `IsEdgeS()` 就换一个新的顶点索引（`crease_cut.h:96-112`），从而把内部锐边拆成两条重合的边界边。

实测（`t_final2.cpp`，3×3 网格沿 x=1 开 crease）：`FACEEDGESEL` 半边数 `before=4 after=4` —— **CreaseCut 不清 crease 位**，可迭代使用。

### 5.3 但"切开 + Loop"会裂开

Loop 只有边界规则、没有折痕规则，于是折痕两侧各自作为独立边界向内收缩，彼此不再贴合。

实测（`t_seamfocus.cpp`，开放正方形沿对角线切开后 Loop）：

```
after cut:  v0(0,0,0) v1(1,0,0) v2(1,1,0.4) v3(0,1,0) v4(0,0,0) v5(1,1,0.4)
after Loop: seam 端点 v0=(0.12500,0.25000,0.05000)  v4=(0.25000,0.12500,0.05000)
            GAP = 0.17678  （不闭合）
            seam 内部中点 v8=(0.5,0.5,0.2)  v9=(0.5,0.5,0.2)  GAP = 0
```

另一组（`t_crease_gap.cpp`，网格沿 x=1 切开，端点落在原网格边界上）：

```
pre-refine duplicate pairs: (0,11) (1,9) (4,10) (7,12)
after Loop:  v0-v11 gap=0.177218   v1-v9 gap=0.250000
             v4-v10 gap=0.000000   v7-v12 gap=0.250000
```

闭合 crease 环（`t_ring3.cpp`，方管底部环切开）：18 对重合点中 10 对保持重合、8 对裂开（最大 0.1768）。

**结论**：折痕"看起来是折痕"依赖两侧局部邻域对称；在端点、角点、以及两侧邻域不对称处会产生可见裂缝。要求几何严密的折痕，必须在此方案外固定顶点或改用真正的 crease-aware 细分。

### 5.4 CreaseCut 自身的两个限制

**(a) 悬空锐边会留下无引用顶点。** 实测（`t_verify_all.cpp` A：闭合八面体，仅一条内部 crease 边）：

```
before cut VN=6 unused=0
after  cut VN=8 unused=2      ← 2 个孤立顶点
```

B（crease 度数为 3 的环+对角线）：`after cut VN=18 unused=4`。
`CutMeshAlongSelectedFaceEdges` 只用"跨过 crease"计数来分配新索引，不检查该顶点上 crease 的度数，因此 `crease-degree==1` 的**端点**会白拿一个没人引用的副本。使用时须配 `Clean::RemoveUnreferencedVertex`。

**(b) 边界上不需要切的地方也会切。** `t_seampos.cpp` 中 `v11=(0,0,0)` 与 `v0` 重合，但 `v0` 并不在 crease 上——因为该顶点是边界顶点，`borderVertexFlag` 分支让它从边界边起走，跨界时仍split。多出来的副本需靠 `RemoveDuplicateVertex` 收拾。

---

## 6. 多边形方案：Catmull-Clark 与 Doo-Sabin

### 6.1 Catmull-Clark：权重不是教科书 CC，且**没有边界规则**

```cpp
// refine_catmullclark.h:80-82  valence = 该顶点被多少个面引用一次
for(auto &f : baseIn.face) for(int i=0;i<f.VN();i++) v_valH[f.V(i)]++;

// refine_catmullclark.h:150-161
float val1 = (val -2.0)/(val*val);      // 自身权重
float val2 = 1.0/(val*val);             // 每个面心 / 每个邻边的权重
refinedOut.vert[...].P() += f.V(i)->P() * val1;
refinedOut.vert[...].P() += center     * val2;
refinedOut.vert[...].P() += edgep      * val2;
```

按面累加后总权重：自身 `n·val1=(n−2)/n`、面心 `n·val2=1/n`、邻边 `n·val2=1/n`（合计 1）。教科书 CC 是自身 `(n−3)/n`、面心 `1/n`、邻边 **`2/n`**：

```
n=3: code self=0.3333 face=0.3333 edge=0.3333 | textbook self=0.0000 face=0.3333 edge=0.6667
n=4: code self=0.5000 face=0.2500 edge=0.2500 | textbook self=0.2500 face=0.2500 edge=0.5000
```

**边权重只有标准的一半**，因此对不规则 valence 与闭合体角点结果不同。实测立方体角点：`(0.194444,0.194444,0.194444)`，教科书值 `2/9≈0.222222`；规则平面网格则恰好正确（`t_grid.cpp`：内部点 `(1.00000,1.00000,0.00000)` 不动）。

**边界完全没处理**。`refine_catmullclark.h` 里 `IsBorder` / `BORDER0` 零命中。实测：

```
# t_cc1.cpp — 单个平面四边形 (0,0) (1,0) (1,1) (0,1)
slot0 ← (0,0) : (1.0000,0.7500)     slot1 ← (1,0) : (0.2500,1.0000)
slot2 ← (1,1) : (0.0000,0.2500)     slot3 ← (0,1) : (0.7500,0.0000)

# t_cc_boundary.cpp — 双四边形条带
k=1 角点 (0,0) → (1.0000,0.7500)    k=1 角点 (2,0) → (1.2500,1.0000)
k=2 边界中点 (1,0) → (1.1250,0.4375)  k=2 边界中点 (1,1) → (0.8750,0.5625)
```

`(2,0)→(1.25,1.0)` 已**离开原网格包围盒**。我用手算复现了 `(0,0)→(1.0,0.75)`：`val=1 ⇒ val1=−1, val2=1`，即 `P' = −P + C + E`（`t_ccverify.cpp` 预测值与实测完全一致）。边界顶点没有被固定、也没有折线保形。

边中点的边界"补偿"也不对。代码注释（`refine_catmullclark.h:136-139`）说"边界面只有一份贡献，权重要加倍"，但实际归一化只除以计数：

```cpp
// refine_catmullclark.h:138-139
for(auto &v : refinedOut.vert)
    if(v_cnH[v]>0) v.P() /= v_cnH[&v];
```

内部边（计数 2）= `0.25(C₁+C₂)+0.25(a+b)`，等价于标准 CC 边点 `(a+b+F₁+F₂)/4` ✔。边界边（计数 1）= `0.5·C + 0.25(a+b)` ✘ —— 实测 (0,1) 边的中点为 `(0.5000,0.2500)`，而正确的二点边界规则应给 `(0.5,0)`（`t_ccverify.cpp`）。**注释描述的"加倍"并没有发生。**

**输出没有 FF 邻接**（`refine_catmullclark.h` 中零 `UpdateTopology`）。实测 `t_checks.cpp`：

```
CC out: halfedges with NULL FFp = 16 / 16
```

### 6.2 Doo-Sabin：边界上产生退化面

```cpp
// refine_doosabin.h:130-153 — 每条边一个"四边形"
tri::UpdateFlags<PolyMeshType>::FaceClearV(base);
for(auto fi=base.face.begin();fi!=base.face.end();++fi) {
    fi->SetV();
    for(int j=0;j<fi->VN();++j)
        if(fi->FFp(j)->IsV()) {          // 边界上 FFp(j)==fi，IsV 为真
            auto newf = tri::Allocator<PolyMeshType>::AddFaces(refined,1);
            newf->Alloc(4);
            face::Pos<FaceType> startPos(&*fi,j);
            ... startPos.FlipV(); ... startPos.FlipF(); ... startPos.FlipV();
```

边界上 `FlipF()` 是"原地翻转"（`FFp(j)==&f`），四次 `FlipV/FlipF` 会重复取到同一角点。实测（`t_poly.cpp`，单个四边形）：

```
ds: VN=4 FN=9
  f0 VN=4: 0 1 2 3
  f1 VN=1: 0            ← 1 边面
  f2 VN=1: 1
  f3 VN=1: 2
  f4 VN=1: 3
  f5 VN=4: 0 1 1 0      ← 只有 2 个不同顶点的"四边形"
  f6 VN=4: 1 2 2 1
  f7 VN=4: 2 3 3 2
  f8 VN=4: 3 0 0 3
```

开放双四边形（`t_ds2.cpp`）：`f2 VN=1: 0`、`f3 VN=2: 4 1`、`f8 VN=4: 0 1 1 0` 等。闭合立方体则干净（6 个 4 边面 + 8 个 3 边面 + 12 个 4 边面 = 26）。

**Doo-Sabin 只适用于闭合多边形网格。** 且输出同样不建 FF：`DooSabin out: halfedges with NULL FFp = 24 / 24`。

### 6.3 两者都不读 crease / faux

`grep "FaceEdgeS\|IsF(\|SetF(\|IsAnyF" refine_catmullclark.h refine_doosabin.h` → 零命中。多边形方案无法表达折痕；`PolyMeshType` 也没有 faux 概念。

---

## 7. 属性如何随拓扑变化传递

### 7.1 `RefineE` 路径（最完整）

**(a) 新顶点**：由 `MIDPOINT` 函子写，基类行为：

```cpp
// refine.h:139-160
nv.P()=   (V0->P()+V1->P())/2.0;
if( tri::HasPerVertexNormal(*mp))     nv.N()= (V0->N()+V1->N()).normalized();
if( tri::HasPerVertexColor(*mp))      nv.C().lerp(V0->C(),V1->C(),.5f);
if( tri::HasPerVertexQuality(*mp))    nv.Q() = (V0->Q()+V1->Q()) / 2.0;
if( tri::HasPerVertexTexCoord(*mp))   nv.T().P() = (V0->T().P()+V1->T().P()) / 2.0;
if(intFunc) (*intFunc)(nv,ep);        // 用户插值器钩子
```

- 组件属性会被平均；**自定义动态 per-vertex 属性不插值，保持 0**。实测 `t_attr.cpp`：`new v4 Q=0.5 C=(5,10,15) custom=0`。
- 边界边与内部边**在属性上无区别**（都取两端点），只有几何取点才有边界特例。

**(b) 新面**：

```cpp
// refine.h:461-467
for(i=1;i<SplitTab[ind].TriNum;++i){
    nf[i]=&*lastf; ++lastf;
    if(RefineSelected || (*fi).IsS()) (*nf[i]).SetS();
    nf[i]->ImportData(*fi);      // 只复制编译期组件
}
```

实测 `t_attr2.cpp`（`face::Color4b=9` + 动态 `fattr=42`）：`face component color: 9 9 9 9; dynamic fattr: 42 0 0 0` —— 只有复用原槽位的 `nf[0]` 保留动态属性。

**(c) 位重映射**（锐边/边界随拓扑传递的唯一机制）：

```cpp
// refine.h:487-503
if(SplitTab[ind].TE[i][j]!=3) {   // TE!=3：新边对应原三角形某条边
    if(orgflag & (MESH_TYPE::FaceType::BORDER0<<(SplitTab[ind].TE[i][j])))
        (*nf[i]).SetB(j);  else (*nf[i]).ClearB(j);
    if(orgflag & (MESH_TYPE::FaceType::FACEEDGESEL0<<(SplitTab[ind].TE[i][j])))
        (*nf[i]).SetFaceEdgeS(j); else (*nf[i]).ClearFaceEdgeS(j);
}
else { (*nf[i]).ClearB(j); (*nf[i]).ClearFaceEdgeS(j); }   // TE==3 内部新边
```

`TE`（`refine.h:91-92,95-106`）把"新三角形第 j 边"映射回"原三角形哪条边"，`3` 表示内部新边。实测（`t_te.cpp`，对角线置 crease、f0 的 0–1 置 faux）：

```
input flagged
  f0: 0 1 2  B=(110) ES=(001) F=(100)
  f1: 0 2 3  B=(011) ES=(100) F=(000)
after full refine (B/ES/F)
  f0: 4 5 6  B=(000) ES=(000) F=(100)
  f2: 0 4 6  B=(100) ES=(001) F=(100)   ← 边界 + crease 正确落到对角线两段
  f3: 4 1 5  B=(110) ES=(000) F=(100)   ← 原边界边 B 正确
  f4: 6 5 2  B=(010) ES=(001) F=(100)
  f5: 0 6 8  B=(001) ES=(100) F=(000)
  f6: 6 2 7  B=(010) ES=(100) F=(000)
  f7: 8 7 3  B=(011) ES=(000) F=(000)
```

B 与 ES 都正确按 TE 传播（对角线被分成两段，每段的两个三角形各带一份）。**但 F（faux）保持 `(100)` 未重映射**——与 §2.3 一致。

3 分模板切换对角线时还手动交换了位（`refine.h:517-525`），但**只交换 B 与 ES，不交换 F**：

```cpp
if((*nf[1]).IsB(0)) (*nf[2]).SetB(1); else (*nf[2]).ClearB(1);
...
if((*nf[1]).IsFaceEdgeS(0)) (*nf[2]).SetFaceEdgeS(1); ...
```

### 7.2 Loop 路径：破坏原顶点颜色，并写入未初始化法向

`RefineOddEvenE` 先算偶点、再调 `RefineE`、最后只回写 `P`/`N`：

```cpp
// refine_loop.h:606-613
RefineE< MESH_TYPE, ODD_VERT > (m, odd, edgePred, RefineSelected, cbOdd);
for(size_t i=0;i<newEven.size();++i) {
    if(updatedList[i]) {
        m.vert[i].P() = newEven[i].first;
        m.vert[i].N() = newEven[i].second;   // 只回写 P 和 N
    }
}
```

**(a) 原顶点颜色被就地覆盖**（even 预扫描时 `lerp` 写回原顶点，随后奇点又从这些被改过的颜色插值）：

```cpp
// refine_loop.h:589-592
if( tri::HasPerVertexColor(m) ) {
    (*fi).V(i)->C().lerp((*fi).V0(i)->C() , (*fi).V1(i)->C(),0.5f);
}
```

实测 `t_color.cpp`（输入 R 通道 `0 40 80 120`）：

```
in : R = 0 40 80 120
out: R = 20 60 50 70 40 55 35 60 45      ← 原顶点已不是 0/40/80/120
```

**(b) 法向是未初始化内存（UB）**。`Centroid::project` 只写 `nv.first`：

```cpp
// refine_loop.h:224-228
inline void project(std::pair<CoordType,CoordType> &nv) const {
    LVector position = sumP / sumW;
    nv.first = CoordType(position.X(), position.Y(), position.Z());   // second 未写
}
```

而 `Point3()` 默认构造**不做任何初始化**（`vcg/space/point3.h:57`: `inline Point3 () { }`），调用方却直接读 `pp.second`：

```cpp
// refine_loop.h:400-403
std::pair<CoordType,CoordType>pp;
proj.project(pp);
nv.P()=pp.first;
nv.N()=pp.second;      // ← 读未初始化值
```

实测 `t_norm.cpp`（输入法向全 `(0,0,1)`，`-O0`/`-O2`）输出法向全为 0；用污染栈（`t_norm_poison.cpp`）则读到 `3.33565e-41, -3.05187e-35` 之类的垃圾位模式——**确认是未初始化读取**。细分后必须 `UpdateNormal`。

自定义 per-vertex 动态属性同样为 0：`t_attr.cpp` → `Loop: v4 Q=0 C=(10) custom=0 | v0 Q=0 custom=100`。

### 7.3 Catmull-Clark / Doo-Sabin：属性全部丢弃

```cpp
// refine_catmullclark.h:85-99
for(auto &v : baseIn.vert)
    tri::Allocator<PolyMeshType>::AddVertex(refinedOut,Point3f(0,0,0));   // 默认构造，不 ImportData
for(auto &f : baseIn.face)
    auto vp = tri::Allocator<PolyMeshType>::AddVertex(refinedOut,PolyBarycenter(f));

// refine_catmullclark.h:192
tri::Allocator<PolyMeshType>::AddQuadFace(refinedOut, v0p, v1p, v2p, v3p);   // 不复制面属性
```

```cpp
// refine_doosabin.h:102
auto newv = tri::Allocator<PolyMeshType>::AddVertex(refined,(fi->V(j)->P()+b)/2.0f);
// refine_doosabin.h:97-100
auto newf = tri::Allocator<PolyMeshType>::AddFaces(refined,1);
newf->Alloc(fi->VN());
```

实测 `t_attr2.cpp`（输入 `C=(10,20,30) Q=5 N=(0,0,1)`、面 `C=(200,100,50)`）：

```
CC  out v0: C=(255,255,255) Q=0 N=(0,0,0) ; new face C=(255,255,255)
DS  out v0: C=(255,255,255) Q=0 ;          face C=(255,255,255)
```

### 7.4 一个正面案例

`PolygonSupport::ImportFromTriMesh`（`polygon_support.h:149`）会传属性：顶点 `ImportData`（`:160`）、面颜色/质量（`:180-182`）：

```cpp
if(tri::HasPerFaceColor(tm) && tri::HasPerFaceColor(pm)) pfi->C()=tfi->C();
if(tri::HasPerFaceQuality(tm) && tri::HasPerFaceQuality(pm)) pfi->Q()=tfi->Q();
```

这是"三角 → 多边形"三步流程中唯一传属性的地方。

---

## 8. 三种情形的可观察差异汇总

| | 拓扑输出 | 几何（新点） | 位 / 属性 |
|---|---|---|---|
| **内部面** | 三边全分 1→4；内部边 1 新点 **2** 新面；两边分时按较短对角线 1→3 | Loop 奇点 `3/8,3/8,1/8,1/8`，偶点 `(1−kβ,β…)`；Butterfly 8 点；CC 边点 = `0.25(C₁+C₂)+0.25(a+b)` | B/ES 按 TE 正确传播；FAUX 整体复制而错位；顶点组件平均，动态属性为 0（`nf[0]` 除外） |
| **边界** | 边界边只 **+1** 新面（不跨面传播） | Loop 奇点 `(l+r)/2`、偶点 `3/4·v+1/8·l+1/8·r`（角点不固定，内缩 0.125）；Butterfly `9/16(l+r)−1/16(l0+r0)`，角点处**越界**（y=−0.0625）；**CC/DS 无规则**：CC 角点 `(0,0)→(1.0,0.75)`、边界边点 `(0.5,0.25)`≠`(0.5,0)`，DS 生成 1/2 边退化面 | B 按 TE 传播；CC/DS 属性全丢 |
| **锐边** | **拓扑零变化**（实测位移恒 0） | 直接置位 → 几何逐位不变；`CreaseCut` 后 → 两侧各自按边界规则收缩，端点**裂开 0.177–0.25** | ES 按 TE 传播到两段新边；`CreaseCut` 不清 ES 位，但会留孤立顶点 / 多余副本 |

---

## 9. 实战建议

| 场景 | 推荐 |
|---|---|
| 三角网格 + 普通边界 | `Loop` 或 `MidPoint`。Loop 自带边界规则（但角点会内缩，需要固定时自己钉住）。 |
| 三角网格 + 开放边界 + 不希望外溢 | **不要用 Butterfly**：角点处 4 点规则会取到绕角顶点，实测产生包裹盒外的点。 |
| 三角网格 + 锐边 | `UpdateFlags::FaceEdgeSelCrease(m, angle)` → `CreaseCut` → 细分。预期折痕处会有 0.1–0.25 量级的缝隙，需要后处理（`RemoveDuplicateVertex` + 投影）或固定缝上顶点。 |
| faux 编码的混合网格 | 先 `PolygonSupport::ImportFromTriMesh` 转 `PolyMeshType`，**不要**在三角网格上直接 `Refine`（faux 错位 → 下游死循环）。 |
| 多边形 + 闭合 | `CatmullClark` 或 `DooSabin`；注意输出**无 FF**，必须 `UpdateTopology::FaceFace`；属性需自行传；CC 权重与教科书不同。 |
| 多边形 + 开放边界 | **当前实现不可用**（CC 角点出包围盒，DS 退化面）。建议先补洞/裁剪，或三角化后走 Loop。 |
| 依赖颜色/quality/自定义属性 | 只有 `RefineE` + `MidPoint`/`Butterfly` 能传组件属性；Loop 会破坏颜色；CC/DS 全丢。 |
| 依赖法向 | **任何细分后都必须 `UpdateNormal`**：Loop 写未初始化法向（UB），CC/DS 输出 0。 |

---

## 10. 复现清单

全部在 `/tmp/vcgtest/`（仓库外）：

| 文件 | 验证内容 |
|---|---|
| `t_mid.cpp` | Loop 边界/内部模板差异；MidPoint 边界无特例；对角线 z 不进入边界点 |
| `t_ignore.cpp` | **逐顶点证明** crease 位对 Loop/MidPoint/Butterfly 位移恒 0 |
| `t_te.cpp` | B / ES / FAUX 在 `RefineE` 后的重映射（B、ES 正确；FAUX 错位） |
| `t_split.cpp` | 两分边选较短对角线 |
| `t_bfly2.cpp` | Butterfly 边界公式 `9/16(l+r)−1/16(l0+r0)` 精确吻合 |
| `t_bfly3/4.cpp` | Butterfly 角点绕行 → 新点 y=−0.0625 越界，并给出路径重放 |
| `t_final2.cpp` | `CreaseCut` 不清 ES 位（4→4） |
| `t_seamfocus.cpp` | 切开后 Loop：端点裂开 0.17678，内部中点贴合 |
| `t_crease_gap.cpp` | 有效内部折痕切开后：0.177 / 0.25 / 0 / 0.25 四对端点缝隙 |
| `t_ring3.cpp` | 闭合 crease 环：18 对重合点中 8 对裂开 |
| `t_verify_all.cpp` | CreaseCut 留孤立顶点（`unused=2` / `unused=4`） |
| `t_faux6.cpp` / `t_faux4.cpp` / `t_faux5.cpp` | faux 编排放错位；细分后 `ExtractPolygon` 死循环 |
| `t_cc1.cpp` / `t_cc_boundary.cpp` / `t_ccverify.cpp` | CC 边界塌陷；`(2,0)→(1.25,1.0)` 出包围盒；手算复现 `−P+C+E` |
| `t_cc_check2.cpp` | CC 权重与教科书对比（边权重仅一半） |
| `t_grid.cpp` / `t_ccc.cpp` | CC 规则网格正确、立方体角点 0.1944≠0.2222 |
| `t_poly.cpp` / `t_ds2.cpp` | DS 边界退化面（VN=1/2 与 `0 1 1 0`）；闭合立方体正常 |
| `t_checks.cpp` | CC/DS 输出 FF 全 NULL（16/16、24/24） |
| `t_attr.cpp` / `t_attr2.cpp` | 组件属性平均、动态属性 0（`nf[0]` 例外）、CC/DS 全丢 |
| `t_color.cpp` | Loop 破坏原顶点颜色（0,40,80,120 → 20,60,50,70） |
| `t_norm.cpp` / `t_norm_poison.cpp` | Loop 法向为 0；污染栈后读到 `3.3e-41` 垃圾位模式（UB） |

编译（任一实验）：

```bash
g++ -std=c++17 -O0 \
  -I/srv/fb2/feng0732/0021-vcglib-3-22d2a3b8 \
  -I/srv/fb2/feng0732/0021-vcglib-3-22d2a3b8/eigenlib \
  /tmp/vcgtest/t_mid.cpp -o /tmp/vcgtest/t_mid && /tmp/vcgtest/t_mid
```

---

## 11. 源码级定论（可直接 grep 核对）

```
1) 锐边不参与细分决策
   grep -n "IsFaceEdgeS\|IsEdgeS" vcg/complex/algorithms/refine.h
     → 仅 487-503 / 522-525 / 1246-1261 / 1281-1284（位搬运），无分支用途
   grep -n "FaceEdgeS\|IsF(\|SetF(" vcg/complex/algorithms/refine_catmullclark.h refine_doosabin.h
     → 零命中

2) 边界只影响传播与取点，不影响 SplitTab
   refine.h:457            ind 仅由三条边是否被分割组成
   refine.h:370-376        边界边不把分割传给邻面（NewFaceNum 只 +1）
   refine.h:595-605、320-322、370-377、459-485   各方案的边界取点公式

3) 多边形方案没有边界规则
   grep -n "IsBorder\|BORDER0" vcg/complex/algorithms/refine_catmullclark.h
     → 零命中；边界顶点落在同一加权公式里（refine_catmullclark.h:150-161）

4) CC/DS 输出不建 FF
   grep -n "UpdateTopology" vcg/complex/algorithms/refine_catmullclark.h → 零命中
   refine_doosabin.h:62 只在输入上建一次 FaceFace(base)
```
