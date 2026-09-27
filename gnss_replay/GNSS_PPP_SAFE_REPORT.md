# 无 RINEX OBS 依赖的 Android PPP 信号策略验证

## 结论

已实现独立的 `--signal-mode ppp-safe`，并设为 TXT 回放默认模式。639 历元测试中，`TXT_PPP_SAFE` 与 RINEX 均为 **Q5=11、Q6=628**；共同 Q6 历元的 TXT－RINEX E/N/U 轨迹差 RMS 为 **0.00564/0.00727/0.02356 m**，最终 E/N/U 差为 **+0.00250/−0.00250/+0.01366 m**。这达到本阶段的厘米量级坐标差目标，远小于 `TXT_FULL` 的高程差 RMS **0.80748 m**。但这是**相对 RINEX**的一致性，不是对独立真值的绝对定位精度证明。

`PPP_SAFE` 在提供和不提供 `--rinex-obs` 两次回放中，639 历元的解算字段、逐历元滤波状态、观测导出和 PPP trace 完全相同；仅 `processing_ms` 随运行时负载变化。因此，可以把 `gnss_adapter` + `gnss_signal_policy` **冻结为本次数据验证过的 Android JNI 输入层 v1**，继续做 JNI 接线和多设备/多日期验证。正式实时运行**不需要 RINEX OBS**。这不表示所有未来手机与产品组合已通过精度验证，也不表示 PPP 模型本身不需独立审计。

## 实现边界与信号表

新增 [gnss_signal_policy.h](../src/gnss_signal_policy.h) 和 [gnss_signal_policy.c](../src/gnss_signal_policy.c)：`gnss_ppp_signal_priority(sys,slot,code)` 返回策略优先级，`gnss_ppp_signal_allowed()` 返回是否可进入当前 PPP/RCB 模型。`gnss_adapter` 保留 Android 频率/code 到 RTKLIB code/slot 的完整识别；`ppp-safe` 先以 signal identity 过滤和排序，再在**同一信号**的重复 Raw 记录间按码时间、ADR、多普勒和 C/N0 选择测量。`full` 继续保留原有研究用的宽信号选择，未删除 BDS 1P/5Q 的映射。

| 星座 | RTKLIB 频槽 | `ppp-safe` 允许 | 本数据送入的 P/L 数 | PPP 实际采用的 code/phase 数 |
| --- | ---: | --- | ---: | ---: |
| GPS | 1 | 1C | 4,430 / 1,950 | 3,620 / 1,922 |
| GPS | 3 | 5Q | 3,602 / 3,532 | 3,194 / 3,140 |
| Galileo | 1 | 1C | 4,323 / 2,001 | 4,113 / 1,735 |
| Galileo | 3 | 5Q | 3,798 / 3,333 | 3,797 / 3,311 |
| BeiDou | 1 | 2I | 8,822 / 3,547 | 8,796 / 3,505 |
| BeiDou | 2 | 7I | 650 / 0 | 532 / 0 |
| BeiDou | 3 | 5P | 0 / 0 | 0 / 0 |

“PPP 实际采用”来自已有 `$PPP_DIAG_SIG` 聚合，不只是 adapter 的输入数。BDS 5P 是正式策略允许的 code，但**这批原始数据没有该信号**，故不能声称已验证其解算效果。Galileo E5b（本次主要为 7Q）仍由 adapter 正确识别，但当前 `ppp.c` 对它已有隔离，`ppp-safe` 不送入；保留在 `full` 供研究。BDS 1P 与 2I 共用槽 1，当前 PPP 主码按 BDS 2I 建模，故正式策略先选 **2I**，不让 1P 靠 C/N0/ADR 分数抢槽。BDS 5Q 在本次 BIA **无 OSB**，而当前模型没有 BDS C5Q 的 RCB 槽；它在 `full` 中曾有 3,491 条 phase 真正进入 PPP 并显著改变坐标/ZTD，因此从 `ppp-safe` 排除，但并未从原始转换能力中删除。

## 产品支持扫描

程序加载 BIA 后扫描 TXT 中实际出现的卫星/code，并输出 `SIGNAL_PRODUCT`；下表是**不同卫星数**，不是历元观测条数。`OSB valid/missing` 读取 `nav->osb_valid[sat-1][code]`，代表当前加载文件的覆盖标志；它本身不代替各历元的 Bias 时间有效性检验。

| 星座/code | 观测到的卫星 | OSB valid | OSB missing | safe 处理 |
| --- | ---: | ---: | ---: | --- |
| GPS 1C | 12 | 11 | 1 | 保留 |
| GPS 5Q | 6 | 6 | 0 | 保留 |
| GAL 1C | 10 | 10 | 0 | 保留 |
| GAL 7Q | 6 | 6 | 0 | 隔离，不送当前 PPP |
| GAL 5Q | 6 | 6 | 0 | 保留 |
| BDS 2I | 14 | 8 | 6 | 保留 |
| BDS 1P | 7 | 7 | 0 | 让位于 BDS 2I 主频槽 |
| BDS 7I | 2 | 0 | 2 | 保留，当前模型有 BDS C7I RCB/legacy DCB 路径 |
| BDS 5P | 0 | 0 | 0 | 策略支持，本次无观测 |
| BDS 5Q | 6 | 0 | 6 | 排除：当前无 RCB，且本次 BIA 无 OSB |

因此**没有**实行“OSB missing 就删掉所有观测”的规则。尤其 BDS 2I、7I 仍可走当前 PPP/legacy DCB 逻辑；本轮未改偏差模型，也未宣称 fallback 的绝对精度已验证。

## 四组 639 历元回放

样本 TXT：`E:\RTKLIB_Data\OBS\RawData_20260922_141724.txt`；对照 RINEX：`E:\RTKLIB_Data\OBS\GEOL00CHN_R_20262650617_00U_01S_MO.rnx`。四组使用同一 Release 可执行文件、同一 NAV/SP3/CLK/BIA/IONEX/VMF3/ATX 和 `configure_ppp()`，ADR uncertainty 默认/显式 1.0 m。没有修改 `ppp.c`、`rtklib.h`、PPP 参数或产品模型。表中 E/N/U 与状态指标均为 **TXT－RINEX**；RINEX 行是基准，不是绝对真值。

| 路径 | Q5/Q6 | 输入 code/phase | PPP used phase | PPP_REJECT | Q6 E/N/U RMS (m) | 最终 E/N/U (m) | ZTD 差 RMS (m) | mean/P95 (ms) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| RINEX | 11/628 | 29,387 / 17,650 | 15,659 | 205 | 对照 | 对照 | 对照 | 4.323 / 5.960 |
| TXT_FULL | 11/628 | 33,278 / 21,599 | 19,259 | 235 | 0.12491 / 0.10489 / **0.80748** | −0.01595 / −0.09908 / **+0.73701** | 0.16560 | 4.765 / 6.846 |
| TXT_PPP_SAFE，有 RINEX 对照 | 11/628 | 25,625 / 14,363 | 13,613 | 164 | **0.00564 / 0.00727 / 0.02356** | **+0.00250 / −0.00250 / +0.01366** | **0.02114** | 3.998 / 5.444 |
| TXT_PPP_SAFE，无 RINEX | 11/628 | 25,625 / 14,363 | 13,613 | 164 | 同上 | 同上 | 同上 | 4.075 / 5.448 |

`PPP_SAFE` 与 RINEX 的 E/N/U 差均值为 −0.00077/−0.00460/−0.01011 m，最大绝对差为 0.02150/0.01841/0.04412 m；628 个共同 Q6 历元的状态均为 Q6。平均有效卫星数比 RINEX 少约 0.022。相比旧 `rinex-compatible` 的 U RMS 0.00749 m，正式策略的 0.02356 m 较大，主要设计差异是它按要求**不送 Galileo E5b**，不是偷偷依赖 RINEX 调参。不要再为追求逐字节 RINEX 一致而重新放宽正式信号集合。

### 滤波状态差异

628 个共同 Q6 历元 TXT－RINEX 的状态差 RMS：

| 模式 | GPS clock (ns) | GAL ISB (ns) | BDS ISB (ns) | BDS C7I RCB (m) | GAL C5Q RCB (m) | GPS C5Q RCB (m) | ZTD (m) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| TXT_FULL | 2.584 | 0.385 | 1.379 | 0.84969 | 0.03943 | 0.06542 | 0.16560 |
| TXT_PPP_SAFE | 0.1017 | **0.7230** | 0.0177 | 0.00446 | 0.04157 | 0.00280 | 0.02114 |

`PPP_SAFE` 的 ZTD 差均值是 −0.00756 m，GAL ISB 差均值 −0.578 ns。坐标轨迹已经是厘米量级，但 **GAL ISB 并非与 RINEX 完全等价**；这是应保留在跨数据集验证清单里的状态差异。BDS C5P RCB 的差为 0 只因本次没有 BDS 5P，不能视为该状态通过测试。`full` 与 `safe` 都有 Q6=628，说明 Q6 数量本身仍不足以判断观测集合是否合适。

## `--rinex-obs` 独立性与默认模式测试

`PPP_SAFE` 有/无 RINEX 两次均为 639 历元；对照 RINEX 读取使用**独立 `nav_t`**，不改变解算器已经加载的产品。比较结果：

- 解算 CSV 去掉运行耗时 `processing_ms` 后，**所有行、所有字段完全相同**；耗时列按运行负载自然不同。
- `.stat` 逐历元滤波状态文件、`_obs.csv` 观测文件和 PPP `.trace` 均**逐字节相同**。
- 不写 `--signal-mode`、也不提供 RINEX 的默认模式另跑 639 历元，解算字段、状态、观测和 trace 与显式 `ppp-safe` 无 RINEX 组同样完全相同。

因此正式信号选择与 PPP 解算已摆脱 RINEX OBS 运行依赖。`rinex-compatible`、`common-only` 仍保留诊断用途；`full` 是研究开关。RINEX 对照不用于 `ppp-safe` 的观测选择。

## 构建、复现和冻结范围

Visual Studio 2022 的 x64 Release 与 Debug 均应构建 `gnss_replay/gnss_replay.vcxproj`；639 历元结果保存在 `C:\Users\Drinkey\Documents\Codex\2026-09-21\ppp-ppp\work\rtklib_audit\ppp_safe`。用 [analyze_ppp_safe.py](analyze_ppp_safe.py) 可重算 ENU/ZTD/ISB/RCB 差和有/无 RINEX 的一致性检查：

```powershell
python 'E:\GNSS\rtklib_app\rtklib_Project2\gnss_replay\analyze_ppp_safe.py' 'C:\Users\Drinkey\Documents\Codex\2026-09-21\ppp-ppp\work\rtklib_audit\ppp_safe'
```

最小的正式入口不需要 `--rinex-obs`，也无需显式指定 signal mode 或 ADR 门限：

```powershell
& 'E:\GNSS\rtklib_app\rtklib_Project2\gnss_replay\x64\Release\gnss_replay.exe' --source txt --txt 'E:\RTKLIB_Data\OBS\RawData_20260922_141724.txt' --mode fast
```

本次冻结的是**TXT/Android Raw → `obsd_t` 的 v1 信号选择与数值转换边界**，不是冻结所有产品配置或宣称厘米级绝对定位。JNI 接入后应在不同手机、不同日期与不同卫星/OSB 覆盖下重跑同样的 P/L/D/LLI、状态和独立真值检查；BDS 5P 与更多缺 OSB 情形尚未由本样本覆盖。ADR reset/slip 仍保留当前历元 `L=ADR/λ` 并设置 `LLI_SLIP`；默认 uncertainty 门限为 1.0 m，`off` 保留作研究，未把 uncertainty 写入 PPP 权重。
