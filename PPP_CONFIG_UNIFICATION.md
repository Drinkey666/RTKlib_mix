# 电脑后处理、replay、Android PPP 配置统一报告

日期：2026-09-30。配置版本：`pc-baseline-20260930-v1`。

## 1. 结果与范围

四个解算入口统一调用 `smartphone_ppp_configure()`：

- 电脑 `rnx2rtkp.c` 后处理；
- 电脑 `gnss_replay.c` 的 TXT 和 RINEX 回放；
- Android `rtk_jni.c` RINEX 后处理；
- Android `gnss_jni.cpp` 实时原始观测 PPP。

以修改前的**当前电脑端默认配置**为准，没有重新优化 PPP 参数。电脑端 `ppp.c`、`rtklib.h`、`rtkcmn.c`、gnss_adapter 和 ppp-safe 均未因本次任务改动。没有将参考坐标注入解算。

## 2. 统一后的设置

| 项目 | 所有入口的默认值 |
| --- | --- |
| 解算 | 静态、前向、浮点 PPP；dynamics=0，AR 关闭 |
| 星座 | GPS + Galileo + BeiDou；GLONASS 不默认输入 |
| 槽位容量/使用数量 | NFREQ=4 / opt.nf=4 |
| 卫星布局 | ENAGLO 已编译；不编译 ENAQZS；MAXSAT=197 |
| 产品模型 | 精密轨道钟差；IONOOPT_EST；TROPOPT_EST |
| 截止高度角 | 15° |
| 相位噪声常数/高度角项 | 0.008 / 0.012 m |
| code/phase ratio | F1=100，F2=70，F3=60，F4=100 |
| ambiguity / ionosphere / trop 过程噪声 | 1e-4 / 1e-3 / 1e-4，沿用电脑端原定义和单位 |
| maxout / GF 周跳阈值 | 30 / 0.20 m |
| 天线与改正 | 卫星 PCV、相位缠绕、地球潮汐、日食处理启用；接收机 PCV 关闭 |
| 强制 code-only 预热 | DOPPWARM=0 |
| Galileo 7Q | GALE5BPHASE=0 |
| BeiDou 码方差倍率 | BDSCODEVAR=1 |
| IONEX 软约束 | IONCONS=1.5，IONCONSINT=1 |
| VMF3 初始化/湿延迟软约束 | VMF3SIG=0.15，VMF3ZWDSIG=0.30，VMF3ZWDINT=300 |
| 其他预处理与抗差 | 完全沿用当前电脑端，不调整阈值 |
| TXT ADR 门限/信号策略 | 1.0 m / ppp-safe，未改 adapter 与 policy |
| 后处理输出 | 椭球高，sstat=2，maxsolstd=0 |

`pppopt`：

```text
-GAP_RESION=120 -IONCONS=1.5 -IONCONSINT=1 -VMF3SIG=0.15 -VMF3ZWDSIG=0.30 -VMF3ZWDINT=300 -DOPPSM=0.90 -DOPPWARM=0 -PREPROC=1 -DOPPSLIP=0.50 -CODEJUMP=30 -MWTHRES=5 -BDSCODEVAR=1 -BDSCODEWARM=120 -WGTELCN=0 -PPPDIAG=1 -PPPQUAR=3,120 -GALE5BPHASE=0
```

这是一套**静止手机测试配置**，不能把它直接称为已经验证的步行/车载动态配置。

## 3. 修正了哪些不一致

原 replay 和 Android 两个入口仍为 nf=3、DOPPWARM=10，缺少第四槽 eratio 以及显式 GALE5BPHASE=0。Android CMake 同时带有 NFREQ=3 和 ENAQZS，造成观测结构容量、卫星编号布局与当前 PC 不同。

Android 活跃核心位于 `app/src/main/cpp/ppp_live/rtklib`，两个 JNI 库实际都从这里编译；旧的 `app/src/main/cpp/rtklib` 并不参与当前这两个库的构建。本次同步活跃核心中的：

- rtklib.h：四槽观测布局；
- rtkcmn.c：当前 PC 的 BDS B1I/B1C 槽位、priority 和 ANTEX 映射；
- ppp.c：当前 PC 的独立 BDS C1P 接收机码偏差状态和 Galileo 7Q 等已有处理；
- geoid.c、pntpos.c：仅同步一处空白差异，无算法改动。

原三个有实质差异的 Android 文件已保存到 `E:/android/test/FTPGet/tools/config_backup_20260930`。这不是另行设计 Android 滤波模型，而是移除版本差异。编译时全部 C/C++ 单元均使用相同的 NFREQ=4，避免 JNI 与 C 核心结构不一致。

**TXT/实时 B1C 仍按原 ppp-safe 策略不输入。** 四槽容量统一不会自动新增 TXT B1C；本次没有改信号映射或策略。RINEX 后处理能够沿用当前 PC 的 B1C 支持。因此，配置相同不等于两条入口拥有完全相同的观测集合。

## 4. 回归证据

### 4.1 K80 2026-03-17，同一 RINEX 和产品

| 测试 | 历元 | Q5 | Q6 |
| --- | ---: | ---: | ---: |
| 修改前电脑默认结果 | 1200 | 3 | 1197 |
| 本次电脑后处理 | 1200 | 3 | 1197 |
| 本次 RINEX replay | 1200 | 3 | 1197 |

- 修改前后电脑 `.pos` 的全部数据行逐字相同。
- 修改前后电脑 `.stat` 逐字相同。
- PC 与 replay 的全部 **138605 行状态文本**一致，包括 POS、CLK、RCB、TROP、ION、SAT；文件换行格式不同不算数值差异。
- PC 与 replay 纬经高按相同输出精度一致。原始导出的小数位不同导致纬/经度最多约 5e-10°、高程最多 0.00005 m 的舍入差，不代表滤波差异。
- CSV 的 time 标签为原始输入时间，状态/.pos 标签为 sol.time；SPP 可修正接收机时钟。本数据二者最多约 3 ms，比较需辨别时间标签含义。两条解算的状态时间也完全一致。
- replay 平均/P95/max 处理耗时：7.780 / 8.719 / 11.949 ms；>1000 ms 处理超时 0 次。耗时与运行机器和负载有关。

证据：`gnss_replay/config_unification/k80_pc.pos`、`k80_pc.pos.stat`、`k80_replay.csv`、`k80_replay.stat`。
可重跑 `python tests/test_config_unification.py`，其检查的不仅是 Q 数量。

### 4.2 已有 TXT/RINEX 样本：完整 639 历元快速回放

| 入口 | Q5 | Q6 | adapter/reader 相位数 | 平均/P95/max ms |
| --- | ---: | ---: | ---: | --- |
| TXT ppp-safe | 2 | 637 | 14363 | 6.944 / 9.187 / 13.387 |
| RINEX | 2 | 637 | 17650 | 6.901 / 9.181 / 12.369 |

- 两组 639 个输入历元按 week + 毫秒 tow 匹配一致，无丢历元；原始浮点时间最多差约 3.6e-8 s。
- 637 个共同 Q6 历元，TXT − RINEX 的 E/N/U RMS：**0.000992 / 0.001292 / 0.003038 m**。
- 末历元 dE/dN/dU：**+0.000405 / +0.000440 / -0.000903 m**。
- 两组 PPP_REJECT 均为 635；这是当前版本/数据下的诊断计数，不作为单独正确性证据。
- 两组不同相位数来自保留的不同输入信号策略；没有为了匹配而删改策略或调整滤波参数。
- 本实验只检验入口差异，不能把上述 RMS 当作对独立真值的定位精度。

证据：`txt-fast_full.csv/.stat/.log`、`rinex_full.csv/.stat/.log`。

### 4.3 TXT 实时调度冒烟测试

同一 20 历元在 TXT_FAST、TXT_REALTIME、RINEX 路径均为 Q5=1、Q6=19。

- TXT_FAST 与 TXT_REALTIME 的 ENU 差为 0，状态及 Q 一致，无丢历元。
- TXT_REALTIME 平均/P95/max：6.652 / 8.790 / 9.324 ms，deadline miss=0。
- 20 历元 TXT_REALTIME − RINEX 的共同 Q6 ENU RMS：0.004478 / 0.005793 / 0.014567 m。
- 20 历元与完整数据的 RMS 不同，统计区间不同，不应混用。

详细结果位于 `gnss_replay/config_unification/txt_smoke/REALTIME_VS_RINEX_REPORT.md`。
只运行了 20 历元 1× 调度测试，**没有声称完整 639 历元都以 1× 回放完毕**。

## 5. 编译和 Android 验证边界

- PC 后处理：x64 Debug、Release 成功。
- PC replay：x64 Debug、Release 成功。
- Android assembleDebug 成功；两个 JNI 库在 arm64-v8a、armeabi-v7a、x86、x86_64 全部编译成功。
- 检查每种 ABI 的 45 个编译单元，均为 NFREQ=4、无 ENAQZS。
- 强制重跑 8 项 JVM 单元测试，失败/错误均为 0。它们不是手机实时 PPP 精度测试。
- PC/Android 活跃核心的 29 个源码/头文件规范化 SHA256 全部相同。
- 正常 CMake 源码校验通过，故意损坏的独立测试快照被拒绝。
- **尚未在手机上跑同一 RINEX 的数值回归，也未采集新的 live 数据。** 目前证明配置/源码/结构布局一致和 PC 入口回归通过，不能宣称已实测 PC/Android 浮点轨迹完全一致。

APK：`E:/android/test/FTPGet/app/build/outputs/apk/debug/app-debug.apk`。

## 6. 后续如何保持一致

PC 唯一默认参数源：`src/smartphone_ppp_config.c/.h`。不要再分别修改四处参数。

Android 独立项目保留这份文件及核心的版本化快照。构建时 `verify_ppp_profile.cmake` 校验 `pc_profile.sha256`；快照被无意修改会阻止构建。检查当前 PC 是否又更新，运行：

```powershell
# 在 PC 项目根目录执行；默认只检查，不复制。
.\tools\sync_smartphone_ppp.ps1 -AndroidRoot E:\android\test\FTPGet

# 仅在新的 PC 基线验证后，显式同步；旧 Android 文件自动备份。
.\tools\sync_smartphone_ppp.ps1 -AndroidRoot E:\android\test\FTPGet -UpdateAndroid
```

Android CMake 可选参数 `-DPPP_PC_SOURCE_DIR=<电脑src路径>` 可以连同 PC 当前源码一起校验；不提供时仍检查 Android 快照自身，保证 Android 项目独立可构建。它**不会在没有指定 PC 路径时自动发现另一个仓库的变化**。

PC 的 RTK_PPP_NF / RTK_PPP_OPTS / RTK_PPP_SYSTEMS / RTK_PPP_GLO / RTK_PPP_EXSATS 显式实验覆盖保留；replay/Android 不自动继承这些环境变量。统一保证的是默认配置，设置实验覆盖后不再属于相同默认条件。

控制台、trace 或 Android PPP_POST/PPP_LIVE logcat 的启动日志会给出 `PPP_CONFIG`。当前 PC 默认签名为 `240c4d93`；签名覆盖日志中列出的配置，不包含产品内容、exsats 等所有 prcopt_t 字段。比较产品还应核对文件哈希和时间覆盖，不能只看配置签名。
