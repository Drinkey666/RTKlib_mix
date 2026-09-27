# Android TXT 与 RINEX 的 PPP 全链路等价性检验

## 结论

同一可执行文件、产品和 `configure_ppp()` 下，两条链路均得到 **628/639 个 Q6**，但这**不能**证明坐标等价。以 RINEX 解为对照，`TXT_FULL` 在 628 个共同 Q6 历元的高程差 RMS 为 **0.8075 m**，最终高程差 **+0.7370 m**；ZTD 差 RMS 为 **0.1656 m**。`TXT_COMMON_ONLY` 的高程差 RMS 为 **0.000045 m**，ZTD 差 RMS 为 **0.000014 m**。因此，原始数据到 `obsd_t` 的共同观测转换基本正确，主要差异来自**送入 PPP 的信号集合**，不是 Q6 判定、坐标输出或必须修改 PPP 模型。

`TXT_RINEX_COMPATIBLE` 的高程差 RMS 为 **0.00749 m**，适合作为当前捕获数据的临时对照策略；它仍需读取该次 RINEX 来学习频槽代码，**尚不能原样用于 Android 实时版**。`TXT_FULL` 对当前产品/偏差模型不宜直接作为实时默认模式。以下所有“误差”均是 **TXT－RINEX 差值**，不是对真实天线坐标的绝对精度评价。

## 数据、方法与输出

- TXT：`E:\RTKLIB_Data\OBS\RawData_20260922_141724.txt`。
- RINEX：`E:\RTKLIB_Data\OBS\GEOL00CHN_R_20262650617_00U_01S_MO.rnx`。
- 639 个历元；四种输入都用同一 Release `gnss_replay.exe`、同一 NAV/SP3/CLK/BIA/IONEX/VMF3/ATX 和 PPP 配置，`--adr-unc-max 1.0`、快速回放、持续同一 `rtk_t`。没有修改 `ppp.c`、`rtklib.h` 或参数、产品与随机模型。
- `--result-file` 逐历元存储 Q、ns、ECEF、经纬高和形式标准差；`--state-dump` 保存 `$CLK/$RCB/$TROP/$ION` 等滤波状态；`--use-dump` 保存 sat/code/slot、P/L、C/N0、高度角、ADR uncertainty、OSB 有效性、当前接收机码偏差槽标签、载波采用标志及后验残差缓存。`--trace-level 3` 仅对 `TXT_FULL` 另跑一遍，用首/末次模型残差核对 prefit/postfit；性能统计均来自 trace 2。码偏差槽标签是对 `ppp.c::ppp_ifb_index()` 的只读诊断映射，若以后更改 PPP 模型须同步检查这个标签。
- 四组逐历元结果与状态流位于 `C:\Users\Drinkey\Documents\Codex\2026-09-21\ppp-ppp\work\rtklib_audit\equiv_{rinex,full,compatible,common}.{csv,stat,trace}`，对应 `_use.csv`。分析脚本：[analyze_equiv.py](analyze_equiv.py)；其输入为该目录。状态流的 TOW 是三位小数，匹配采用 0.01 s 容差；坐标直接按同历元周/TOW 配对。ENU 坐标系取对应 RINEX 解位置。

四组模式：`RINEX` 直接解码 RINEX OBS；`TXT_FULL` 使用适配器当前全部选择结果；`TXT_RINEX_COMPATIBLE` 先从 **当前 `rinex.c` 实际解出的频槽代码**学习系统/频槽白名单，再在原始信号选择时限制代码；`TXT_COMMON_ONLY` 在选择之前用同历元、同卫星、同 code 的 RINEX 记录过滤原始测量，末尾再做一次安全核对，仅为诊断。后两者没有在代码中硬编码“北斗永远选 2I”或“永远选 1P”。

## 四组完整 PPP 结果

表中的 ENU RMS、最终 ENU、ZTD RMS 均为 **TXT－RINEX**；RINEX 一行是对照本身，不代表绝对零误差。`phase` 是送入 `rtkpos()` 的载波观测数；`used phase` 是 trace 中已进入 PPP 残差统计的载波数。时间为本次 PC 快速回放的每历元处理时间，非 Android 性能。

| 路径 | Q5/Q6 | code / phase / Doppler 输入 | Q6 平均有效 ns | used phase | PPP_REJECT 行数 | E/N/U RMS (m) | 最终 E/N/U (m) | ZTD 差 RMS (m) | mean / P95 (ms) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| RINEX | 11/628 | 29,387 / 17,650 / 29,387 | 17.705 | 15,659 | 205 | 对照 | 对照 | 对照 | 4.266 / 5.825 |
| TXT_FULL | 11/628 | 33,278 / 21,599 / 34,221 | 20.115 | 19,259 | 235 | 0.125 / 0.105 / **0.807** | −0.016 / −0.099 / **+0.737** | **0.166** | 4.755 / 6.893 |
| TXT_RINEX_COMPATIBLE | 11/628 | 29,387 / 17,651 / 30,310 | 17.705 | 15,659 | 205 | 0.00445 / 0.00059 / **0.00749** | +0.00015 / +0.00007 / −0.00061 | 0.000069 | 4.253 / 5.773 |
| TXT_COMMON_ONLY | 11/628 | 29,387 / 17,650 / 29,387 | 17.705 | 15,659 | 205 | 0.000016 / 0.000063 / **0.000045** | +0.000004 / −0.000038 / +0.000009 | 0.000014 | 4.225 / 5.780 |

628 个共同 Q6 历元没有 Q 状态不一致。逐历元 TXT－RINEX ENU 均值及最大绝对差：

| 路径 | E/N/U 均值 (m) | E/N/U 最大绝对差 (m) | 最终 X/Y/Z 差 (m) |
| --- | ---: | ---: | ---: |
| TXT_FULL | +0.0548 / −0.0606 / +0.7638 | 0.4383 / 0.3163 / 1.7952 | −0.1507 / +0.6677 / +0.2911 |
| TXT_RINEX_COMPATIBLE | +0.00158 / −0.00017 / +0.00163 | 0.0456 / 0.00324 / 0.0826 | −0.000010 / −0.000580 / −0.000248 |
| TXT_COMMON_ONLY | +0.000011 / −0.000059 / +0.000033 | 0.000036 / 0.000110 / 0.000143 | −0.000010 / +0.000025 / −0.000028 |

经纬高逐历元和 ECEF 坐标一同记录在上述 CSV；ENU 由同一 RINEX 位置计算，最终大地高差分别为 +0.737010、−0.000609、+0.000009 m。形式 XYZ 标准差差的 RMS：FULL 为 0.0033/0.0088/0.0036 m，COMPATIBLE 为 0.00079/0.00033/0.00021 m，COMMON_ONLY 小于 0.000001 m。形式标准差相近并不意味着系统偏差不存在。

### PPP 状态也不等价

下表为 628 个共同 Q6 历元逐历元 TXT－RINEX 状态差 RMS；时钟/ISB 单位 ns，其余单位 m。ION 是共有卫星-历元状态的 RMS（16,482 对）。

| 模式 | GPS clock | GAL ISB | BDS ISB | BDS C7I RCB | GAL C5Q RCB | GPS C5Q RCB | ION |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| FULL | 2.584 | 0.385 | **1.379** | **0.850** | 0.0394 | 0.0654 | **0.1082** |
| RINEX_COMPATIBLE | 0.0262 | 0.0040 | 0.0078 | 0.00747 | 0.00134 | 0.00238 | 0.00167 |
| COMMON_ONLY | 0.00038 | 0.00043 | 0.00036 | 0.000074 | 0.000075 | 0.000047 | 0.000023 |

FULL 的 ZTD 差均值为 −0.1350 m、RMS 0.1656 m；COMPATIBLE 为 −0.000006/0.000069 m；COMMON_ONLY 为 +0.000002/0.000014 m。当前 BDS C5P RCB 在这些运行中未激发，差为 0，**不能**据此称该偏差已受验证。FULL 的有效卫星数比 RINEX 平均多 2.409，表明额外信号确实改变了解算几何和状态，而非仅被输入层记录。

## 多出的 3,955 条相位到底是什么

以同历元、同卫星、同 code 的 RINEX phase 为基准，TXT_FULL 的 **3,955 条额外相位**按星座/code/频槽分组如下。输入净数仅多 3,949，是因为 FULL 同时挤掉了 6 条 RINEX 的北斗 2I 相位：17,650 − 6 + 3,955 = 21,599。

| 星座/code/槽 | 额外 phase | 平均 C/N0 (dB-Hz) | 平均高度角 | ADR 不确定度均值 (m) | OSB valid / missing | PPP 实际 used phase | phase PPP_REJECT | used phase 后验 RMS (m) | 对应 code 后验 RMS (m) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BDS 5Q / 3 | **3,677** | 31.74 | 47.00° | 0.5016 | **0 / 3,677** | **3,491** | 6 | 0.1538 | 无已采用 code |
| BDS 1P / 1 | 277 | 28.25 | 22.01° | 0.00309 | 277 / 0 | 277 | 0 | 0.0463 | 14.32（这 277 条对应缓存） |
| GAL 1C / 1 | 1 | 31.60 | 0° | 0.5009 | 1 / 0 | 0 | 0 | 不适用 | 不适用 |

`used phase` 由现有 `$PPP_DIAG_SIG` 与频槽 `vsat`、后验残差交叉核对，并非仅统计 adapter 输出；trace 3 首/末次模型残差进一步复核了 BDS 5Q 的 used phase 后验 RMS 0.1538 m、BDS 1P 的 0.0463 m。对应初次（prefit）模型残差 RMS 分别为 34.50、32.79 m，包含尚未初始化的模糊度，**不能**拿它单独断言信号质量差。BDS 1P 这 277 条的 code 初次/最终模型残差 RMS 为 28.76/14.32 m；FULL 全部 BDS 1P 中真正参与 PPP 的 code 共 380 条、后验 RMS 17.56 m。BDS 5Q 在 FULL 中没有被当前 PPP 模型采用的 code（`nP=0`），尽管有 3,491 条 phase 参与；它的 code residual RMS 因而不可计算。`PPP_REJECT` 是显式拒绝事件计数，不等于输入减 used；其余未用项还可能因高度角、有效性或模型条件被跳过。

分类判断：

- **BDS 5Q**：原始 Android 中存在，而本次 RINEX reader 未输出，是 A 类；同时该 BIA 在 `nav->osb_valid[sat-1][code]` 上 **3,677/3,677 缺 OSB**，是 C 类。轨道/钟差覆盖在运行 trace 中显示 `IN_RANGE`，没有证据把它归为 D 类。phase 后验 RMS 高于 BDS 2I（约 0.0436 m），但不能仅凭后验残差把它定为异常 E；强耦合状态可吸收偏差。当前 PPP 接收机码偏差槽不覆盖 BDS 5Q code。
- **BDS 1P**：原始有、RINEX 本次未输出，是 A 类；OSB **277/277 有效**，不是 C 类。与 BDS 2I 共享槽 1，是 F 类频槽竞争；当前模型没有单独的 BDS 1P 接收机码偏差状态。code 后验残差较大，状态/坐标差异值得警惕，但不能把所有偏差唯一归因于 1P code。
- **GAL 1C**：仅 1 条，0° 高度角且未用于 PPP，属于本次 RINEX 未出现的少量输入差异，不足以判定映射 B 错误。

在共同信号上，先前严格观测对照有 **17,644 条 phase 匹配，米级相位差 RMS 0.000066 m、最大 0.000127 m**，LLI 全部一致；剩余 6 条属于共享频槽选择丢失，不是 ADR 数值换算错误。COMMON_ONLY 预先过滤后，17,650 条 RINEX phase 全部匹配，code 也无不一致。由此可把主要原因定位在信号选择及其产品/偏差支持，而不是 Android P/L/D/LLI 转换。

### 两类信号的消融验证

`--exclude-raw-signal` 在信号选择前屏蔽一个 code，其余全部保持不变，仍以 RINEX 解为对照。该实验**不**证明哪个解更接近独立真值，但可以说明哪些额外信号驱动轨迹偏离。

| TXT_FULL 的消融 | Q6 | 输入 phase | Q6 高程差均值/RMS (m) | 最终高程差 (m) | ZTD 差 RMS (m) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 无消融 | 628 | 21,599 | +0.764 / 0.807 | +0.737 | 0.166 |
| 排除 BDS 5Q | 628 | 17,922 | +0.540 / 0.608 | +0.380 | 0.093 |
| 排除 BDS 1P | 628 | 21,328 | +0.208 / 0.235 | +0.326 | 0.059 |
| RINEX-compatible（两类均未进入） | 628 | 17,651 | +0.00163 / 0.00749 | −0.00061 | 0.000069 |

所以，额外 phase **确实进入 PPP 并改变估计**。对“与 RINEX 等价”的目标而言，FULL 的额外 BDS 信号造成显著偏差；对“绝对定位是否更准”的目标，还需要同一手机天线相位中心的独立真值、多次采集及完整产品支持，不能仅按 RINEX 一致性判断。

## 6 条 BDS 2I/1P 频槽冲突与策略

原 FULL 对 C23 的 5 个历元、C33 的 1 个历元在槽 1 选了 1P，RINEX reader 实际给出 2I。`rinex.c` 根据观测头代码列表和 `getcodepri()/set_index()` 选择同频槽代码；本次头文件/reader 实际输出 G:1C/5Q、E:1C/7Q/5Q、C:2I/7I。新的 `rinex-compatible` 是从解出的 RINEX 观测**学习本次系统/槽代码**后在 adapter 选取前限制，因而恢复这 6 条 2I，没有把 2I 写成普适优先规则。结果是 17,650 条 RINEX phase 全部匹配，另有 1 条未被 PPP 采用的 GAL 1C。`common-only` 再按历元精确匹配，消除这最后 1 条。

当前 Android 实时默认策略应是**产品/模型支持感知的、可配置的 RINEX 式信号优先级**：先按 RTKLIB 槽与 `getcodepri` 规则排序，再检查该 code 的 OSB、接收机偏差模型与实际残差表现；同槽冲突要有确定规则和日志。不能使用需要旁边 RINEX 文件的 `common-only`；也不应把当前 `rinex-compatible` 的学习机制误称为可离线独立运行的生产实现。更不能简单地永远选 BDS 2I 或永远选 1P。`full` 可保留为研究开关，不宜在这套产品/模型下直接上 Android 默认。

## 是否已可进入 Android JNI

**可进入 JNI 接口原型与速度测试，但尚不满足把它作为可靠 Android 实时 PPP 默认链路的条件。**理由：共同 P/L/D/code/LLI 和状态/坐标在 COMMON_ONLY 下已达到紧密一致，说明转换边界可复用；但 FULL 的状态、高程轨迹显著不同，且最大的额外信号 BDS 5Q 缺 OSB，当前高等价模式仍依赖本次 RINEX 文件。进入生产前应完成无 RINEX 的可配置信号优先级/产品支持检查，针对其它手机与采集时段复验，再用同一天线的独立真值评估**绝对**定位精度。此轮不调整 PPP 滤波、权重、IONEX、VMF3 或偏差模型。

## 复现命令

在项目根目录 PowerShell 中，以下命令的公共产品采用程序目前的默认路径；换日期需显式提供相应产品。四次必须使用同一 exe 和产品，并分别指定不同输出文件：

```powershell
$exe = 'E:\GNSS\rtklib_app\rtklib_Project2\gnss_replay\x64\Release\gnss_replay.exe'
$txt = 'E:\RTKLIB_Data\OBS\RawData_20260922_141724.txt'
$rnx = 'E:\RTKLIB_Data\OBS\GEOL00CHN_R_20262650617_00U_01S_MO.rnx'
$out = 'C:\Users\Drinkey\Documents\Codex\2026-09-21\ppp-ppp\work\rtklib_audit'
& $exe --source rinex --rinex-obs $rnx --mode fast --adr-unc-max 1.0 --result-file "$out\equiv_rinex.csv" --state-dump "$out\equiv_rinex.stat" --use-dump "$out\equiv_rinex_use.csv" --trace "$out\equiv_rinex.trace"
& $exe --source txt --txt $txt --rinex-obs $rnx --signal-mode full --mode fast --adr-unc-max 1.0 --result-file "$out\equiv_full.csv" --state-dump "$out\equiv_full.stat" --use-dump "$out\equiv_full_use.csv" --trace "$out\equiv_full.trace"
& $exe --source txt --txt $txt --rinex-obs $rnx --signal-mode rinex-compatible --mode fast --adr-unc-max 1.0 --result-file "$out\equiv_compatible.csv" --state-dump "$out\equiv_compatible.stat" --use-dump "$out\equiv_compatible_use.csv" --trace "$out\equiv_compatible.trace"
& $exe --source txt --txt $txt --rinex-obs $rnx --signal-mode common-only --mode fast --adr-unc-max 1.0 --result-file "$out\equiv_common.csv" --state-dump "$out\equiv_common.stat" --use-dump "$out\equiv_common_use.csv" --trace "$out\equiv_common.trace"
python 'E:\GNSS\rtklib_app\rtklib_Project2\gnss_replay\analyze_equiv.py' $out
```

`--trace-level 3` 能重跑 FULL 并生成卫星级 prefit/postfit 诊断，但其大量日志写入会影响耗时比较。`--adr-unc-max off` 继续可用；本数据有效 ADR 不确定度最大约 0.51 m，1.0 与 off 的 phase/Q6 计数相同。此不确定度目前只用于门限与统计，不加入 PPP 权重。
