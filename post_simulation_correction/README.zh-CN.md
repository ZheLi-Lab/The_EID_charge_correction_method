# 模拟后的解析电荷校正

这是对原仓库的**独立补充方法**，不替换原有 EID 拟合程序，也不替换
`eid_additional_correction.py`。不能默认把两套校正同时加到同一个自由能结果上。

本目录的 `charge_correction.py` 和 `apbs_tmp.in` 与整理前的原文件逐字节一致；
此次仅增加组织、说明与基础测试，没有修改公式或数值流程。

## 安装

建议 Python 3.9+；另外需要单独安装 AMBER/AmberTools 的 `ambpdb` 和 APBS，
并确保两者可从 `PATH` 调用。

在仓库根目录执行：

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r post_simulation_correction/requirements.txt
command -v ambpdb
command -v apbs
```

默认最多并行运行 5 个 APBS 任务，每个网格为 257³，需要注意内存。

## 准备输入并运行

需要复合物和配体溶液体系各一套匹配的、含溶剂的 AMBER 拓扑与坐标文件，
两侧盒子都应近似立方体。复合物坐标还必须是可通过标准输入交给 `ambpdb`
读取的文本格式，不能直接使用二进制轨迹。

每个体系在独立空目录运行，避免中间文件覆盖。以下命令从仓库根目录开始：

```sh
mkdir -p runs/system_01
cd runs/system_01
python ../../post_simulation_correction/charge_correction.py \
  /absolute/path/to/inputs/complex.rst7 \
  /absolute/path/to/inputs/complex.prmtop \
  /absolute/path/to/inputs/ligand.rst7 \
  /absolute/path/to/inputs/ligand.prmtop \
  -1 3 -r MOL -o results.csv
```

请把示例绝对路径替换为实际输入。参数顺序为：

```text
复合物坐标 复合物拓扑 配体溶液坐标 配体溶液拓扑 配体电荷 复合物溶质电荷
```

- `-1`、`3` 仅为示例，单位为元电荷。这里的复合物电荷不包含溶剂和中和离子；
  程序中蛋白电荷 = 复合物溶质电荷 − 配体电荷。
- `-r MOL` 是配体残基名，默认 `MOL`。
- 当前有效半径计算不支持配体电荷为零。
- `-o results.csv` 会额外保存汇总文件，不会取消默认输出。

## 输出及符号

- `charge_corrections.csv`：complex、ligand、binding 三组结果，
  包含 NET、USV、RIP、EMP、DSC 及总和，单位同时提供 kJ/mol、kcal/mol。
- `charge_corrections_details.csv`：两侧盒长、水分子数、电荷、RIP、有效半径与分项结果。
- `complex/`、`ligand/`：APBS 输入、DX 电势网格、PQR 和 RIP 报告。
- 当前运行目录另有 PDB、PQR 等中间文件；重复运行会覆盖同名文件。

**代码规定 binding = ligand − complex。** 应先核对你所采用的热力学循环和
湮灭/去耦自由能符号，再决定如何应用校正。

## 发布前必须注意

完整记录见 [VALIDATION.md](VALIDATION.md)。尤其是：

1. 代码虽然检查 AMBER 周期盒，但真正用于计算的是 PDB 坐标跨度，
   单侧公式还会把 x 方向跨度截断为整数 Å，并非直接使用周期盒长。
2. 配体溶液侧的 PQR 来自复合物中的配体，而不是独立配体坐标。
3. 蛋白筛选使用所有 `ATOM` 行，配体使用残基名字串搜索，
   可能重复选择配体或带入水/离子，必须人工检查生成的 PQR。
4. 水分子数假定三原子水且残基名为 `WAT`；默认溶剂参数按 TIP3P 设置。
5. 详细输出的无 pandas 回退存在变量未初始化问题，因此依赖清单包含 pandas。
6. 缺少真实分子输入与参考结果，目前不能声称完整工作流已验证。

此次刻意不静默修复上述问题，避免改变文章已有结果。建议后续以单独版本
修复，并对比修复前后的数值。

## 基础测试

在仓库根目录运行：

```sh
python -m unittest discover -s post_simulation_correction/tests -v
```

仅验证合成数据下的公式内部一致性、binding 差分和 CSV 格式，
不等于分子体系上的科学验证。

更多输入输出细节和参考文献见 [English README](README.md)。
