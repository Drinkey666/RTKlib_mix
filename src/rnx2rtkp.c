/*------------------------------------------------------------------------------
* rnx2rtkp.c : read rinex obs/nav files and compute receiver positions
*
*          Copyright (C) 2007-2016 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:55:16 $
* history : 2007/01/16  1.0 new
*           2007/03/15  1.1 add library mode
*           2007/05/08  1.2 separate from postpos.c
*           2009/01/20  1.3 support rtklib 2.2.0 api
*           2009/12/12  1.4 support glonass
*                           add option -h, -a, -l, -x
*           2010/01/28  1.5 add option -k
*           2010/08/12  1.6 add option -y implementation (2.4.0_p1)
*           2014/01/27  1.7 fix bug on default output time format
*           2015/05/15  1.8 -r or -l options for fixed or ppp-fixed mode
*           2015/06/12  1.9 output patch level in header
*           2016/09/07  1.10 add option -sys
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "./rtklib.h"

#define PROGNAME    "rnx2rtkp"          /* program name */
#define MAXFILE     16                  /* max number of input files */

/* help text -----------------------------------------------------------------*/
static const char *help[]={
"",
" usage: rnx2rtkp [option]... file file [...]",
"",
" Read RINEX OBS/NAV/GNAV/HNAV/CLK, SP3, SBAS message log files and compute ",
" receiver (rover) positions and output position solutions.",
" The first RINEX OBS file shall contain receiver (rover) observations. For the",
" relative mode, the second RINEX OBS file shall contain reference",
" (base station) receiver observations. At least one RINEX NAV/GNAV/HNAV",
" file shall be included in input files. To use SP3 precise ephemeris, specify",
" the path in the files. The extension of the SP3 file shall be .sp3 or .eph.",
" All of the input file paths can include wild-cards (*). To avoid command",
" line deployment of wild-cards, use \"...\" for paths with wild-cards.",
" Command line options are as follows ([]:default). A maximum number of",
" input files is currently set to 16. With -k option, the",
" processing options are input from the configuration file. In this case,",
" command line options precede options in the configuration file.",
"",
" -?        print help",
" -k file   input options from configuration file [off]",
" -o file   set output file [stdout]",
" -ts ds ts start day/time (ds=y/m/d ts=h:m:s) [obs start time]",
" -te de te end day/time   (de=y/m/d te=h:m:s) [obs end time]",
" -ti tint  time interval (sec) [all]",
" -p mode   mode (0:single,1:dgps,2:kinematic,3:static,4:static-start,",
"                 5:moving-base,6:fixed,7:ppp-kinematic,8:ppp-static,9:ppp-fixed) [2]",
" -m mask   elevation mask angle (deg) [15]",
" -sys s[,s...] nav system(s) (s=G:GPS,R:GLO,E:GAL,J:QZS,C:BDS,I:IRN) [G|R]",
" -f freq   number of frequencies for relative mode (1:L1,2:L1+L2,3:L1+L2+L5) [2]",
" -v thres  validation threshold for integer ambiguity (0.0:no AR) [3.0]",
" -b        backward solutions [off]",
" -c        forward/backward combined solutions [off]",
" -i        instantaneous integer ambiguity resolution [off]",
" -h        fix and hold for integer ambiguity resolution [off]",
" -bl bl,std     baseline distance and stdev",
" -e        output x/y/z-ecef position [latitude/longitude/height]",
" -a        output e/n/u-baseline [latitude/longitude/height]",
" -n        output NMEA-0183 GGA sentence [off]",
" -g        output latitude/longitude in the form of ddd mm ss.ss' [ddd.ddd]",
" -t        output time in the form of yyyy/mm/dd hh:mm:ss.ss [sssss.ss]",
" -u        output time in utc [gpst]",
" -d col    number of decimals in time [3]",
" -s sep    field separator [' ']",
" -r x y z  reference (base) receiver ecef pos (m) [average of single pos]",
"           rover receiver ecef pos (m) for fixed or ppp-fixed mode",
" -l lat lon hgt reference (base) receiver latitude/longitude/height (deg/m)",
"           rover latitude/longitude/height for fixed or ppp-fixed mode",
" -y level  output solution status (0:off,1:states,2:residuals) [0]",
" -x level  debug trace level (0:off) [0]",
" --version display release version",
};
/* show message --------------------------------------------------------------*/
extern int showmsg(const char *format, ...)
{
    va_list arg;
    va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
    fprintf(stderr,"\r");
    return 0;
}
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}

/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i=0;i<(int)(sizeof(help)/sizeof(*help));i++) fprintf(stderr,"%s\n",help[i]);
    exit(0);
}
/* rnx2rtkp main -------------------------------------------------------------*/
int main() {
    int i, n, ret;
    double tint = 0.0;              /* 求解时间间隔 (0:使用观测数据默认间隔) */
    gtime_t ts = { 0 }, te = { 0 }; /* 历元时段始末控制变量 */
    char* infile[MAXFILE], outfile[MAXSTRPATH] = { '\0' };

    // ========================================================================
    // 1. 【核心修改区：文件路径配置】 (注意 Windows 下路径用双斜杠 \\)
    // ========================================================================

    /* 设置结果输出文件的绝对路径 */
    char result_file[] = "E:\\RTKLIB_Data\\PPP_Result\\smartphone_ppp_solution.pos";
    strcpy(outfile, result_file);

    /* * 输入文件列表 (PPP 核心文件)
     */
    char infile_[MAXFILE][MAXSTRPATH] = {
        /*"E:\\RTKLIB_Data\\OBS\\GNSS00GEO_R_20261421459_11M_01S_MO.rnx", */
        "E:\\RTKLIB_Data\\OBS\\smoothed_output.rnx",
        "E:\\RTKLIB_Data\\NAV\\BRDC00WRD_S_20261420000_01D_MN.rnx",       /* 2. 广播星历文件 */
        "E:\\RTKLIB_Data\\SP3\\WHU0MGXRTS_20261420000_01D_30S_ORB.SP3",   /* 3. 精密星历文件 (SP3) */
        "E:\\RTKLIB_Data\\CLK\\WHU0MGXRTS_20261420000_01D_05S_CLK.CLK",   /* 4. 精密钟差文件 (CLK) */
         "E:\\RTKLIB_Data\\IONEX\\whug1420.26i",
        ""
    };

    /* 模型参数文件路径设置 (PPP 强依赖这些模型文件) */
    filopt_t filopt = { "", "", "", "", "", "", "", "" };
    // 手机本身没有天线相位中心(PCV)模型，但必须提供卫星的 .atx 文件，否则依然会有系统误差
    strcpy(filopt.satantp, "E:\\RTKLIB_Data\\Tables\\igs20.atx"); /* 卫星天线参数文件 */
    strcpy(filopt.iono, "E:\\RTKLIB_Data\\IONEX\\whug1420.26i");
    strcpy(filopt.dcb, "E:\\RTKLIB_Data\\BIA\\WUM0MGXRTS_20261420000_01D_05M_OSB.BIA");
    // ========================================================================
    // 2. 【时间段设置】
    // ========================================================================
  // i 是遍历 infile_ (源数组) 的读取索引
// n 是填入 infile (目标指针数组) 的写入索引兼计数器
    for (i = 0, n = 0; i < MAXFILE; i++) {

        // strcmp 是 C 语言比较字符串的函数。
        // 如果 infile_[i] 不是空字符串 "" (返回值不等于 0)
        if (strcmp(infile_[i], "") != 0) {

            // 核心逻辑：获取字符串的内存地址，存入指针数组
            // n++ 表示存入当前位置后，有效文件数量 n 自动加 1
            infile[n++] = &infile_[i][0];
        }
    }

    /* ========================================================================
        【解算参数配置：移动端全参数估计 (EST) 极限约束策略】
        ======================================================================== */
    prcopt_t prcopt = prcopt_default;
    solopt_t solopt = solopt_default;

    /* --- 核心估计模式 --- */
    prcopt.mode = PMODE_PPP_STATIC;
    prcopt.navsys = SYS_GPS | SYS_GLO | SYS_GAL | SYS_CMP;
    prcopt.nf = 2;

    prcopt.ionoopt = IONOOPT_EST; /* 坚守电离层估计 */
    prcopt.tropopt = TROPOPT_EST; /* 坚守对流层估计 */

    /* 🌟 1. 物理遮挡：切除低仰角，防止非线性大气误差 🌟 */
    prcopt.elmin = 20.0 * D2R; /* 必须 20 度！20度以下的对流层投影函数极其不准，会直接撕裂滤波器 */

    /* 🌟 2. 状态过程噪声约束 (Process Noise - 极其核心) 🌟 */
    // prn[1] 控制电离层每秒允许的变化量。默认可能偏大，这里压到 1E-4。
    prcopt.prn[1] = 1E-3; /* 电离层过程噪声：从 1E-4 放宽到 1E-3，加速收敛 */

    // prn[2] 控制对流层每秒允许的变化量。对流层极其稳定，压死到 1E-5。
    // 这样滤波器即使拿到几十米的伪距噪声，也绝对不敢把它塞进对流层里！
    prcopt.prn[2] = 1E-5;

    /* 🌟 3. 初始方差约束 (Initial Variance) 🌟 */
    // 既然你传入了 IONEX (.26i) 文件，电离层初值是有一定准度的，限制其初始搜索范围
    prcopt.std[1] = 2.0; /* 电离层初始标准差 (米) */
    prcopt.std[2] = 0.2; /* 对流层初始标准差 (米) - 天顶延迟一般在 2.3m 左右，盲猜也不会偏太多 */

    /* 🌟 4. 彻底抛弃伪距信任 (Variance Mapping) 🌟 */
    // 因为你在估计全参数，伪距的毒性会被成倍放大。必须极致降权！
    prcopt.err[1] = 500.0; /* 伪距方差是相位的 500 倍 */
    prcopt.err[2] = 0.015; /* 稍微放宽一点相位的包容度，防假周跳 */
    prcopt.err[3] = 10.0;  /* 仰角惩罚指数 */

    prcopt.maxinno[0] = 0.0; /* 依然关闭剔除，靠方差权重硬抗 */
    prcopt.maxinno[1] = 0.0;

    /* --- 常规设置 --- */
    prcopt.tidecorr = 1;
    prcopt.posopt[0] = 0;
    prcopt.posopt[1] = 0;
    prcopt.posopt[2] = 1;
    // ========================================================================
    // 4. 【执行解算】
    // ========================================================================
    printf("Starting Smartphone Un-combined PPP Processing...\n");
    traceopen("E:\\RTKLIB_Data\\ppp_debug.trace");
    tracelevel(3);
    long t1 = clock();
    ret = postpos(ts, te, tint, 0.0, &prcopt, &solopt, &filopt, infile, n, outfile, "", "");
    long t2 = clock();

    if (!ret) fprintf(stderr, "%40s\r", "Processing completed.");

    printf("\n* Result saved to: %s\n", outfile);
    printf("* The total time for running the program: %6.3f seconds\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
    printf("Press 'Enter' key to exit...\n");
    getchar();

    return ret;
}
