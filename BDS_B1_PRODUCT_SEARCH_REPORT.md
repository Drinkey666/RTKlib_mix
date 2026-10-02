# BDS B1I/B1C 本地产品完整性搜索

## 结论

本次只检查本地产品和已有观测共存记录，未运行或修改 PPP。扫描范围为 `E:\GNSS`、`E:\RTKLIB_Data`、`E:\collect data`、`E:\BaiduNetdiskDownload` 及 `C:\Users\Drinkey\Documents\Codex` 的全部子目录，共 95 个目标文件：14 个 BIA、11 个 GNSS ANTEX、29 个 CLK、41 个 SP3；没有发现 BSX 或 `.bias` 文件。运行基点为 commit `7055987`。逐条观测使用实际 `C2I/L2I` 与 `C1P/L1P` 精确匹配，OSB 有效期按左闭右开区间判断，不外推、不用 `L1X` 或 `L1D` 代替 `L1P`。

| Session（GPST） | RAW coexist | BIAS complete | ATX complete | SP3 complete | CLK complete | FULL_MODEL_COMPLETE |
|---|---:|---:|---:|---:|---:|---:|
| 2026-09-22 06:17:45.471–06:28:23.999，GPS week 2437/day 2 | 4304 | 0 | 4304 | 4304 | 4304 | **0** |
| 2026-03-21 02:11:12.000–02:21:10.999，GPS week 2410/day 6 | 4800 | 4800 | 4800 | 4800 | 4800 | **4800** |

计数单位是同历元同 PRN 的 B1I/B1C 共存配对。`BIAS complete` 要求同一份 BIA 对该配对同时给出精确的 `C2I/L2I/C1P/L1P`；`ATX complete` 要求该 PRN 在历元有效的 `C02/C01`；SP3/CLK 要求该 PRN 的有效记录夹住观测历元。完整配对必须同时满足全部五项。CSV 中保留逐历元和逐星证据。

## 9 月产品：更新了，但目标时段仍有相位缺口

更新后的 [WUM RTS BIA](E:/RTKLIB_Data/BIA/WUM0MGXRTS_20262650000_01D_05M_OSB.BIA) SHA-256 为 `A44AE401ED51DAC7B1893EED3C597B348428DE6992A5359FFF3349AD7D165652`；旧副本 SHA-256 为 `192E8C5ECFA4C5A4EA46AFCE42433FAE160FA2CA1400317CDA85B80451BBFA73`。旧副本的 `C25 L2I/L1P` 相位 OSB 只有 47 段，最晚到 04:00；更新版有 264 段，并非只覆盖至 04:00。不过更新版对全部 7 颗共存卫星（C06/C23/C24/C25/C33/C41/C42）的 `L2I`、`L1P` 都在 **06:05–06:45 GPST** 缺档，而观测处于 06:17:45–06:28:24。不能根据文件整体首尾为 00:00–24:00 就判断中间连续覆盖，也不能保持 06:05 的常数至 06:45。

| 精确信号 | 更新 RTS BIA 在 4304 个目标配对中的覆盖 |
|---|---:|
| `C2I` code | 4304 |
| `L2I` phase | 0 |
| `C1P` code | 4304 |
| `L1P` phase | 0 |
| `C1D/L1D/C1X/L1X` | 0/0/4304/0；不用于 `1P` 替代 |

原配置所用 [IGS20 ATX](E:/RTKLIB_Data/Tables/igs20.atx) 的 `C06` 只有 `C02/C06/C07`，没有 `C01`。本机另一份 [更新的 IGS20 ATX](<E:/collect data/2026.5.19N/igs20.atx>) 中，`C06` 新天线块自文件所列 2026-04-10 16:00 起包含 **`C01/C02/C05/C07`**，且 9 月所有 7 颗参与卫星都具备这四个频点。因此当前本机 **不缺 C06/C01 校准**；此前以旧 ATX 得出的缺失结论需要更新。详见 [逐星 ATX 表](gnss_replay/bds_b1_product_search/atx_frequency_matrix.csv)。

9 月 [WUM RTS SP3](E:/RTKLIB_Data/SP3/WUM0MGXRTS_20262650000_01D_05M_ORB.SP3) 仅到 03:55，[WUM RTS CLK](E:/RTKLIB_Data/CLK/WUM0MGXRTS_20262650000_01D_05S_CLK.CLK) 仅到 03:59:30。可覆盖目标历元的 [WUM NRT SP3](E:/RTKLIB_Data/SP3/WUM0MGXNRT_20262641000_02D_05M_ORB.SP3) 与 [WUM NRT CLK](E:/RTKLIB_Data/CLK/WUM0MGXNRT_20262641000_02D_05M_CLK.CLK) 是另一产品系列；RTS BIA + NRT SP3/CLK 只能列为 **CANDIDATE_UNVERIFIED**。即使暂时忽略系列兼容问题，目标时段相位 OSB 仍为 0/4304，无法启用 B1I+B1C PPP。

## 3 月产品：本机已有覆盖完整的同系列组合

3 月原始 [RINEX OBS](E:/RTKLIB_Data/OBS/GNSS00GEO_R_20260800210_10M_01S_MO.rnx) 的 BDS 观测码确为 `C1P/L1P` 和 `C2I/L2I`。原配置的 3 月 **RTS** BIA 对 `L2I` 覆盖 4800/4800，但仅有 `L1X`，没有 `L1P`；因此它仍不满足 B1C 相位需求。本机另有同一天 **WUM FIN** 组合：

| 类型 | 最多完整配对的文件 |
|---|---|
| BIA | [WUM FIN OSB](E:/RTKLIB_Data/PPP_Result/diag_20260321/WUM0MGXFIN_20260800000_01D_01D_OSB.BIA) |
| SP3 | [WUM FIN orbit](E:/RTKLIB_Data/PPP_Result/diag_20260321/WUM0MGXFIN_20260800000_01D_05M_ORB.SP3) |
| CLK | [WUM FIN clock](E:/RTKLIB_Data/PPP_Result/diag_20260321/WUM0MGXFIN_20260800000_01D_30S_CLK.CLK) |
| ATX | [本机 IGS20](E:/RTKLIB_Data/Tables/igs20.atx) |

该 FIN BIA 对 C24/C25/C26/C34/C35/C39/C40/C44 的 `C2I/L2I/C1P/L1P` 均覆盖 4800/4800；FIN SP3、FIN CLK 逐星覆盖全部历元，IGS20 ATX 对这 8 颗卫星的 `C01/C02/C05/C07` 均有效。BIA、SP3、CLK 同属 `WUM0MGXFIN` / WHU 系列，按本次字段完整性门槛得到 **4800/4800 FULL_MODEL_COMPLETE**。CLK 头部标注使用 `IGS20_2408.ATX`；本机文件名为 `igs20.atx`，确有目标频点，但无法从本地文件名证明它与 CLK 生产时所用的 ANTEX 版本逐字相同。这个版本溯源问题应在真正运行实验前核实；本次不因此改变 PPP 配置。

## 产品组合与界限

所有候选组合列在 [product_combinations.csv](gnss_replay/bds_b1_product_search/product_combinations.csv)，含 SP3/CLK/BIAS 中心、产品系列、ATX 文件和可完整覆盖的配对数。9 月全部 22 个时段相交的本地组合最多为 **0**；3 月 WUM FIN 同系列组合可到 **4800**。不同系列的组合保留 `CANDIDATE_UNVERIFIED` 标记；没有把混用当成可直接运行的 PPP 产品链。

本次只是**产品完整性检查**。它不证明偏差基准、天线版本和 PPP 观测模型已完全兼容，更不代表 B1C filter 已经启用或定位性能已获验证。特别是 9 月仍缺目标时段的相位 OSB。需要的外部产品和确切时段列在 [MISSING_BDS_B1_PRODUCTS.md](MISSING_BDS_B1_PRODUCTS.md)。

## 可复核产物

- [product_inventory.csv](gnss_replay/bds_b1_product_search/product_inventory.csv)：95 个文件的路径、类型、中心、内容时段、采样及系统。
- [bias_signal_catalog.csv](gnss_replay/bds_b1_product_search/bias_signal_catalog.csv)：每份 BIA 的所有 BDS 卫星对八种指定信号的记录数、首末时刻和最大内部缺档。
- [bias_exact_signal_coverage.csv](gnss_replay/bds_b1_product_search/bias_exact_signal_coverage.csv)：每份候选 BIA、每颗星、每个精确信号的历元覆盖。
- [bias_candidate_coverage.csv](gnss_replay/bds_b1_product_search/bias_candidate_coverage.csv)：各 BIA 的目标配对覆盖汇总。
- [atx_frequency_matrix.csv](gnss_replay/bds_b1_product_search/atx_frequency_matrix.csv)：逐 PRN 的 `C01/C02/C05/C07`。
- [BDS_B1_PRODUCT_MATRIX.csv](gnss_replay/bds_b1_product_search/BDS_B1_PRODUCT_MATRIX.csv)：逐 session、PRN、signal 的产品矩阵。
- [pair_completeness.csv](gnss_replay/bds_b1_product_search/pair_completeness.csv)：9104 个配对的各项布尔门槛。
- [pair_counts.csv](gnss_replay/bds_b1_product_search/pair_counts.csv)：4304/4800 汇总。

复现：从仓库根目录运行 `python tests/search_bds_b1_products.py`。脚本只读本机产品和既有共存 CSV，写入独立诊断目录；不修改 PPP 源码、滤波开关或产品配置。
