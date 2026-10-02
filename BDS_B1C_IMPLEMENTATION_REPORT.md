# 北斗 B1I + B1C 同时参与 PPP：实现与验证

## 使用方式与边界

默认 `ppp-safe` / 原 VS 后处理仍用原来的三频配置，逐历元坐标保持不变。新功能必须显式启用：

- `gnss_replay.exe`：`--signal-mode ppp-b1c`，支持 `--source rinex` 和 `--source txt-fast`。无需 RINEX OBS 才能处理 TXT。
- VS 主程序：环境变量 `RTK_PPP_B1C=1`。可用 `RTK_PPP_OBS`、`RTK_PPP_NAV`、`RTK_PPP_SP3`、`RTK_PPP_CLK`、`RTK_PPP_BIA`、`RTK_PPP_ATX`、`RTK_PPP_OUTPUT` 指向同一会话的输入、产品和独立结果文件。未设置开关时仍是旧模式。

新模式中北斗 `2I/B1I` 位于 slot 0、`5P/B2a` 位于 slot 2、`1P/B1C` 位于 slot 3；`7I/B2b` 不进入 PPP。B1I 和 B1C 具有不同物理频率、独立模糊度、独立周跳/质量状态。B1C 伪距另估一个相对 `C2I` 的接收机 `BDS_C1P` 码偏差状态。首版只开放精确的 `1P`，RINEX 同槽同时存在 `1P/1D` 时优先 `1P`；不能把 `1D/1X` 当作 `1P` 的偏差替代。

产品门控在每颗星、每个历元执行：`C2I+C1P` 码 OSB 均在有效期内，且卫星 ANTEX `C02+C01` 频点和天线块有效，B1C 码才进入 PPP；`L2I+L1P` 相位 OSB 也均在有效期内且相位观测存在，B1C 相位才进入 PPP。相位 OSB 按精确观测码读取、转换为米后从相位观测中扣除。缺相位产品时允许已校正的 B1C 码单独使用，但不伪称双相位 PPP。trace 中 `$B1C_GATE` 给出逐历元数量和跳过原因。

`NFREQ` 与相关结构由 3 扩为 4，因此**必须将所有使用本项目 `rtklib.h` 的程序一起重新编译**；不能把新旧对象文件、DLL/JNI 库混用。本次未修改独立的 Android App 工程，Android 实时端尚未启用这一新模式。

## 固定输入测试

2026-03-21 RINEX 使用同一天 WUM FIN SP3/CLK/BIA、`igs20_2408.atx`、同一 IONEX/VMF3。600 历元：

| 模式 | Q0/Q5/Q6 | B1C 门控通过码 | B1C 门控通过相位 | B1C 最终 `phase_vsat=1` | PPP_REJECT | 平均/P95耗时 |
|---|---:|---:|---:|---:|---:|---:|
| 旧 `ppp-safe` | 0/11/589 | 0 | 0 | 0 | 76 | 7.25/11.13 ms |
| 新 `ppp-b1c` | 0/11/589 | 4800 | 4768 | 4689 | 76 | 8.75/12.78 ms |

两模式共同轨迹三维差 RMS 0.390 m、最大 0.634 m，末历元椭球高 477.353 m → 477.448 m。这证明新增观测确实影响估计，不证明在所有数据上都提升精度。按项目既有独立参考坐标计算，该**单次**三维 RMS 约 5.378 m → 5.056 m，仍不是厘米级。VS 主程序独立运行得到 600 个 `.pos` 历元、11 个 Q5、589 个 Q6，并输出 `BDS_C1P` 接收机偏差状态。

2026-09-22 TXT 使用当前本地 NRT 轨道/钟差、RTS BIA 组合，仅用于入口与门控测试，**不作为完整产品兼容或定位精度证明**。639 历元中，B1C 原始候选 4407、码门控通过 3581、相位门控通过 0；Q0/Q5/Q6 为 0/11/628。该时段 `L2I/L1P` 相位 OSB 缺档，代码正确禁止 B1C 相位入滤波，未拿 `L1X` 代替 `L1P`。

旧模式回归：与修改前已保存的 3 月 600 历元、9 月 639 历元结果比较，排除耗时字段后，两组逐历元坐标和状态均逐行一致。Release 主程序、Debug 主程序、Release `gnss_replay` 均编译通过。新编译的 ANTEX 测试覆盖 8 颗卫星 × 5 个频点，读取差异为 0；`tests/test_bds_b1c_model_readiness.py` 通过。

测试输出：`gnss_replay/bds_b1c_implementation_smoke/`，重点查看 `final_rinex_b1c_solution.csv`、`final_rinex_b1c_use.csv`、`final_rinex_b1c.trace`、`final_txt_b1c_solution.csv`、`final_txt_b1c.trace` 及 `vs_b1c.pos.stat`。

## 注意

开关只保证在产品与观测满足门槛时按上述模型参加解算；不能保证换日期、手机或 BIA 后精度必然改善。B1C 相位缺失时，trace 会给出原因；此时 Q6 可由其他信号维持，不能把 Q6 不变当成 B1C 相位已使用。
