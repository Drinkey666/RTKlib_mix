# PC 文件选择窗口：RINEX / Android 原始 TXT 二选一

## 使用方法

1. 双击 `E:/GNSS/rtklib_app/rtklib_Project2/x64/Release/rtklib_Project2.exe`。正常启动打开中文文件选择窗口。
2. 先选择输入类型：**RINEX OBS 后处理**，或 **Android 原始 TXT（1×模拟实时）**。只出现一行观测输入；切换类型清空观测与输出路径，保留产品选择。
3. 观测行点“浏览”，选择对应的 `.rnx/.obs` 或原始 `RawData_*.txt`；再选择NAV、SP3、CLK、BIA、IONEX、VMF3前后两期、格网高程及ATX。结果行点“浏览”，选择保存目录和文件名：RINEX保存POS，TXT保存CSV。
4. 点“开始解算”。RINEX仍在控制台显示后处理进度；TXT会自动打开“模拟实时 PPP 坐标”窗口，逐历元刷新坐标。回放完成后点击“关闭窗口”。

也可以在VS打开解决方案，选择 **x64 / Debug或Release**，按 `Ctrl+F5`。主项目的x64构建输出已统一到根目录 `x64/<Configuration>`。

若以前在VS调试环境中设置了 `RTK_PPP_OBS`，程序会保留无界面的自动测试行为；在“项目属性 → 调试 → 命令参数”填 `--select-files` 即可强制打开窗口。

## 选择时注意

- RINEX模式保持原来的 `postpos()` 后处理链路。TXT模式复用已有 `gnss_replay` → `gnss_adapter` → `ppp-safe` → `rtkpos()`，完全不生成或读取RINEX OBS，也不需要提供对应RINEX。两条链互斥，一次只运行一条。
- TXT须是当前采集App/GNSS Logger风格的Android原始Raw记录，不是任意TXT坐标文件。原始解析、adapter及信号规则没有在本次改动中扩展。
- TXT以原始GPST间隔1×回放，用单调时钟计算目标时间，而非每历元固定等待1秒；通常10分钟数据约需10分钟。保留原始时刻间隔，不重采样、不填补采集中断，不将首次产品加载耗时计入GNSS时间轴。全程只初始化一次滤波器，每个有效历元调用一次 `rtkpos()`。
- 产品应对应同一次采集的日期和时间。窗口检查文件存在、可读取和输出路径，不会自动匹配产品日期，也不能仅凭文件存在证明产品有效。
- VMF3前、后两期应夹住观测时间；同时选择独立的 `orography_ell_5x5` 格网高程文件。
- ATX、格网高程可以一直使用本地文件；没有下载步骤。
- 初次窗口不预填旧日期产品。之后记住“开始解算”时的选择，保存到 `%LOCALAPPDATA%/RTKLIB_PPP/file_selection.ini`。
- 选择观测后，若结果路径为空，会建议带时间戳的 `.pos`（RINEX）或 `.csv`（TXT）文件名；可用保存对话框改到 `PPP_Result` 文件夹。
- RINEX结果旁生成 `.pos.stat`、`.pos.trace`，以及原有events文件。TXT结果旁生成 `.csv.stat`、`.csv.trace`、`.csv.obs.csv`（P/L/D/code/SNR/LLI），控制台显示逐历元位置、Q值和耗时，并汇总平均/P95/最大处理耗时和deadline miss。已有结果/诊断文件会提示确认覆盖；取消不启动解算。
- TXT结果CSV字段为 `week,tow,Q,ns,nobs,x_m,y_m,z_m,lat_deg,lon_deg,h_m,std_x_m,std_y_m,std_z_m,processing_ms`，不是RTKLIB POS格式；湿延迟/钟差/RCB等状态在同名 `.stat` 中。
- Unicode文件对话框支持中文显示；RTKLIB核心仍使用Windows当前代码页的文件接口。无法无损转换或过长的路径会被拒绝，不会用乱码路径继续解算。

## 改动范围

文件窗口：`src/pc_file_picker.c/.h`；TXT桥接新增 `src/pc_txt_replay.c/.h`、`src/gnss_replay.h`；已有 `gnss_replay.c` 仅提取可复用入口并允许接收PC入口已生成的同一份 `prcopt_t`。主VS工程以 `GNSS_REPLAY_EMBEDDED` 编译它，独立 `gnss_replay.exe` 仍保留自己的main。

新增 `tests/test_pc_observation_source.py`；既有 `tests/test_pc_file_picker.py` 继续检查RINEX。没有启动子进程拼接命令，也没有将TXT转成临时RINEX。

PPP算法、默认参数、`ppp.c`、`rtklib.h`、`rtkpos.c`、adapter、signal policy及Android均未修改，哈希已核对。

Release x64字符集设置与Debug统一为MultiByte，以匹配现有RTKLIB的窄字符Windows文件接口；新窗口明确调用Unicode API。原生窗口编译依赖Windows自带库，不需要额外GUI框架。

保留 `--no-gui`、`--check-files`、`--help`。新增 `--source rinex|txt-realtime|txt-fast`，其中txt-fast仅用于快速诊断。自动测试可用 `RTK_PPP_OBS` 或 `RTK_PPP_TXT`，二者同时非空立即返回输入冲突；显式source与环境输入冲突也拒绝。

TXT自动模式必须指定 `RTK_PPP_TXT`、所有对应产品路径和新的 `RTK_PPP_OUTPUT` CSV输出路径，不回退旧RINEX观测或旧结果。可用 `--max-epochs 20` 做短实时测试。TXT默认 `ppp-safe`、ADR uncertainty上限1.0 m，与现有输入层一致；已有PPP配置及环境参数原样传给replay，不重设随机模型。

正常人工使用不需要手填路径或改源码。若VS调试环境设有RTK_PPP_OBS/TXT，请填 `--select-files` 强制显示选择窗口。历史产品路径仅用于原有无界面模式，不自动填入首次GUI。

选定VMF3无法解析时停止，不静默回退。主入口还检查输出是否本次更新且包含有效解算行，避免读取失败后误报“结果保存成功”；并未在本任务重写postpos内部所有错误分支。

## 验证结果

- Debug/Release x64均构建通过，正式根目录EXE已更新。已有其他模块的编译警告未在本任务全面清理。
- 已实际检查中文窗口布局与原生文件浏览对话框。用户在使用测试窗口时停止了自动界面操作，没有代替用户提交其文件选择。
- 8项入口/路径检查通过：help、未知参数、正常检查、缺OBS、缺VMF3、输出选成目录、输出父目录不存在、输出覆盖OBS。失败返回2且不启动PPP。
- 实际完整运行新程序03-17、03-21，以及旧程序03-17，共3组，全部输出到隔离目录。
- 03-17：1200行，Q6=1197，Q5=3；与修改前EXE的坐标行相同，`.stat`逐字节相同；与已验证replay的138470行状态逐行相同。
- 03-21：600行，Q6=599，Q5=1；与已验证的同日replay基线全部状态逐行相同。没有用03-17产品处理它。
- 对replay的状态比较仅忽略CRLF/LF换行差异，没有舍弃坐标或状态字段。
- 验证目录：`gnss_replay/file_picker_audit/run_20261002_191403_097736`，包含日志、pos/stat/trace及summary.json。用户原始解算结果未覆盖。

复现：

```text
python tests/test_pc_file_picker.py --exe x64/Release/rtklib_Project2.exe
```

## 2026-10-02 新增TXT入口验证

- 主程序Release/Debug x64、独立replay Release均编译通过；已有其他模块警告仍存在。
- RINEX回归：03-17的1200历元（Q6=1197、Q5=3）、03-21的600历元（Q6=599、Q5=1），完整状态仍与既有基线逐行相同。结果位于 `gnss_replay/file_picker_audit/run_20261002_194049_617869`。
- TXT 09-22完整快速回放639历元：Q6=637、Q5=2、Q0=0；平均/P95/max处理耗时7.507/10.130/32.502 ms。这个数量是当前已有算法和配置的实测，不是早期版本628/639的承诺，也不是精度判断依据。
- TXT前20历元1×回放与快速回放、独立replay比较：全部坐标和Q/ns/std一致（不比较耗时），完整stat及观测dump逐字节一致，ENU差为0，没有少历元；Q6=19、Q5=1。
- 20历元实时测试实际总耗时19.997 s（包含启动加载）；平均/P95/max处理耗时5.704/6.623/6.680 ms，>1000 ms处理和调度deadline miss均为0。本次没有执行完整639历元的1×长时回放，不将短时验证表述为长时验证。
- 8项输入冲突/路径失败测试及TXT文件检查通过；另用GC诊断参数验证原PC配置确实传入replay，而非被恢复为默认星座。
- `ppp.c`、`rtklib.h`、`rtkpos.c`、默认共享配置、adapter、signal policy哈希未变。Android工程未改动。
- 验证目录 `gnss_replay/observation_source_audit/20261002_194141_594504` 中保存本次所有CSV/stat/trace/log/summary.json，未覆盖原始文件。
- computer-use只读检查了新增单选按钮和中文窗口布局；发现用户在操作文件窗口后停止界面自动操作，没有代替用户开始解算或保存选择。

复现TXT测试：

```text
python tests/test_pc_observation_source.py --exe x64/Release/rtklib_Project2.exe --replay-exe gnss_replay/x64/Release/gnss_replay.exe
```

这是观测输入界面改进，不是精度改进。TXT对比RINEX本身可能仍有信号差异；本次证明入口复用和实时等待不会额外改变已有TXT解算。原先03-21高程系统偏差仍需按系统审查计划另外处理。

## TXT实时坐标显示窗口（2026-10-02）

选择Android原始TXT并开始解算后，自动显示独立窗口：

- 当前GNSS历元GPST日期/时间、week/TOW（不是电脑当前时间）。
- Q值及说明：Q=5为单点解，Q=6为PPP浮点解，Q=0为无有效解。Q=6不保证已收敛到厘米级。
- 使用卫星数 `sol.ns` 和adapter输入卫星数 `nobs`，二者分开展示。
- 纬度、经度、高程（椭球高）；无有效解时清空坐标显示，不显示上一历元旧坐标。
- adapter与滤波处理耗时、已处理历元数，以及最近的滚动记录。快速诊断模式下界面可合并刷新，完整数据仍逐历元保存在CSV。
- 结果保存路径。回放结束后保留最后历元，便于查看。

“停止回放”会请求后台在当前加载/历元结束后停止，关闭产品及结果文件，保存已经完成的历元；窗口仍保留最后结果。运行中点击窗口关闭按钮也会先安全停止再关闭。加载产品期间停止可能需要等待当前产品读取结束。停止后不继续同一个会话，重新开始会创建新的回放会话。

CSV及stat在逐历元通知窗口前刷新缓冲；trace、obs dump继续沿用原日志机制，正常结束/停止时关闭并保存。此刷新不是断电保护，也不保证强制结束进程后所有诊断日志完整。

### 实现边界

新增 `src/pc_live_view.c/.h`。主线程处理窗口，单独后台线程执行已有TXT回放；回放通过 `gnss_replay_epoch_t` 只读快照通知界面。窗口没有 `rtk_t/nav_t` 指针，不能改滤波状态、产品或参数。等待仍按单调时钟的原始GNSS时间目标进行，仅增加可取消检查；没有每秒重初始化滤波器。

RINEX入口未增加此窗口，原后处理行为不变；TXT与RINEX依然二选一。该窗口模拟采集节奏，不是手机在线连接，也未接入NTRIP。

无人值守输入路径也可加 `--live-view` 显示窗口，例如 `--source txt-realtime --live-view --max-epochs 20`；该选项只用于TXT。窗口模式完成后需人工关闭；普通不带该选项的命令行回放仍自动结束。

### 新窗口验证结果

- 主项目Release/Debug x64、独立replay Release和回调测试程序均构建通过。核心6个文件（ppp、rtklib.h、rtkpos、共享默认配置、adapter、signal policy）SHA256保持原值。未修改Android。
- 639历元只读回调测试：每个回调的时间、Q、ns/nobs、ECEF、纬经高均与同时保存的CSV字段相同；与增加窗口前的639历元相比，除处理耗时外全部CSV字段相同，stat及obs dump逐字节相同。
- 20历元1×回调回放耗时20.003 s（含启动加载），结果与原20历元回放相同。
- 使用computer-use技能检查真实中文窗口布局、状态、坐标和滚动记录；实际窗口20历元回放的坐标及完整stat与无窗口回放相同，完成后关闭测试窗口，未操作用户自己的采集或选择设置。
- 第3历元后请求停止：返回明确停止码3，恰好保存3历元及3次坐标回调，约3.020 s；没有强制终止工作线程。
- RINEX回归：03-17共1200历元，Q6=1197、Q5=3；03-21共600历元，Q6=599、Q5=1；完整状态仍与基线逐行一致。
- 本次实际窗口只测试20历元；不将其表述为完整639历元的长时窗口实时验证。

测试数据和摘要：`gnss_replay/live_view_audit/20261002_204019_212712`；本次RINEX回归：`gnss_replay/file_picker_audit/run_20261002_203848_529911`。原始结果未覆盖。

复现只读接口及停止测试：

```text
python tests/test_pc_live_view.py
```

加 `--ui` 会运行实际20历元窗口测试，完成后需点击“关闭窗口”让测试继续检查结果。
