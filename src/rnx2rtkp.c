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
#include "smartphone_ppp_config.h"

#define PROGNAME    "rnx2rtkp"          /* program name */
#define MAXFILE     16                  /* max number of input files */

/* Optional overrides keep controlled A/B runs from overwriting the normal
 * observation, solution and trace files. Empty variables retain defaults. */
static int env_copy(char *dst, size_t capacity, const char *name)
{
    const char *value = getenv(name);

    if (!value || !*value) return 1;
    if (strlen(value) >= capacity) {
        fprintf(stderr, "%s is too long\n", name);
        return 0;
    }
    strcpy(dst, value);
    return 1;
}

/* help text -----------------------------------------------------------------*/
static const char* help[] = {
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
extern int showmsg(const char* format, ...)
{
    va_list arg;
    va_start(arg, format); vfprintf(stderr, format, arg); va_end(arg);
    fprintf(stderr, "\r");
    return 0;
}
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}

/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i = 0; i < (int)(sizeof(help) / sizeof(*help)); i++) fprintf(stderr, "%s\n", help[i]);
    exit(0);
}
/* rnx2rtkp main -------------------------------------------------------------*/
int main() {
    int i, n, ret;
    double tint = 0.0;              /* 求解时间间隔 (0:使用观测数据默认间隔) */
    gtime_t ts = { 0 }, te = { 0 }; /* 历元时段始末控制变量 */
    const char* infile[MAXFILE];
    char outfile[MAXSTRPATH] = { '\0' };

    // ========================================================================
    // 1. 【核心修改区：文件路径配置】 (注意 Windows 下路径用双斜杠 \\)
    // ========================================================================

    /* 设置结果输出文件的绝对路径 */
    char result_file[] = "E:\\RTKLIB_Data\\PPP_Result\\K80-GNSS00GEO_R_20260760533_20M_01S_MO.pos";
    strcpy(outfile, result_file);

    /* * 输入文件列表 (PPP 核心文件)
     */
    char infile_[MAXFILE][MAXSTRPATH] = {
        /*"E:\\RTKLIB_Data\\OBS\\GNSS00GEO_R_20261421459_11M_01S_MO.rnx", */
        "E:\\RTKLIB_Data\\OBS\\K80-GNSS00GEO_R_20260760533_20M_01S_MO.rnx",
        "E:\\RTKLIB_Data\\NAV\\BRDC00IGS_R_20260760000_01D_MN.rnx",       /* 2. 广播星历文件 */
        "E:\\RTKLIB_Data\\SP3\\WUM0MGXFIN_20260760000_01D_05M_ORB.SP3",   /* 3. 精密星历文件 (SP3) */
        "E:\\RTKLIB_Data\\CLK\\WUM0MGXFIN_20260760000_01D_30S_CLK.CLK",   /* 4. 精密钟差文件 (CLK) */
         "E:\\RTKLIB_Data\\IONEX\\COD0OPSFIN_20260760000_01D_01H_GIM.INX",
        ""
    };

    if (!env_copy(outfile, sizeof(outfile), "RTK_PPP_OUTPUT") ||
        !env_copy(infile_[0], sizeof(infile_[0]), "RTK_PPP_OBS")) return -1;

    /* 模型参数文件路径设置 (PPP 强依赖这些模型文件) */
    filopt_t filopt = { 0 };
    // 手机本身没有天线相位中心(PCV)模型，但必须提供卫星的 .atx 文件，否则依然会有系统误差
    strcpy(filopt.satantp, "E:\\RTKLIB_Data\\Tables\\igs20.atx"); /* 卫星天线参数文件 */
    strcpy(filopt.iono, "E:\\RTKLIB_Data\\IONEX\\COD0OPSFIN_20260760000_01D_01H_GIM.INX");
    strcpy(filopt.dcb, "E:\\RTKLIB_Data\\BIA\\WUM0MGXFIN_20260760000_01D_01D_OSB.BIA");
    if (!env_copy(filopt.dcb, sizeof(filopt.dcb), "RTK_PPP_BIA")) return -1;
    // 开启 WGS84 椭球高 到 正常高(海拔高) 的转换

    // 指定高程异常模型文件的绝对路径 (注意双斜杠)

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
    /* All entry points share the current PC baseline. */
    smartphone_ppp_configure(&prcopt, &solopt);
    {
        /* Reproducible B1C ablation without changing any other PPP option. */
        const char *nf_env = getenv("RTK_PPP_NF");
        if (nf_env && *nf_env) {
            if (!strcmp(nf_env, "3")) prcopt.nf = 3;
            else if (strcmp(nf_env, "4")) {
                fprintf(stderr, "RTK_PPP_NF must be 3 or 4\n");
                return -1;
            }
        }
    }
    if (!env_copy(prcopt.pppopt, sizeof(prcopt.pppopt), "RTK_PPP_OPTS")) return -1;
    {
        const char *systems = getenv("RTK_PPP_SYSTEMS");
        if (systems && *systems) {
            int mask = 0;
            const char *c;
            for (c = systems; *c; c++) {
                switch (*c) {
                case 'G': mask |= SYS_GPS; break;
                case 'R': mask |= SYS_GLO; break;
                case 'E': mask |= SYS_GAL; break;
                case 'C': mask |= SYS_CMP; break;
                case 'J': mask |= SYS_QZS; break;
                default:
                    fprintf(stderr, "Invalid RTK_PPP_SYSTEMS character: %c\n", *c);
                    return -1;
                }
            }
            prcopt.navsys = mask;
        }
    }
    {
        /* GLONASS FDMA code needs receiver-bias validation before it can be
         * a safe default for phone PPP. This explicit opt-in wins over the
         * general system mask for reproducible A/B processing. */
        const char *glo = getenv("RTK_PPP_GLO");
        if (glo && *glo) {
            if (!strcmp(glo, "1")) prcopt.navsys |= SYS_GLO;
            else if (!strcmp(glo, "0")) prcopt.navsys &= ~SYS_GLO;
            else {
                fprintf(stderr, "RTK_PPP_GLO must be 0 or 1\n");
                return -1;
            }
        }
    }
    if ((prcopt.navsys & SYS_GLO) && NSATGLO == 0) {
        fprintf(stderr, "GLONASS requested but ENAGLO is not enabled in this build\n");
        return -1;
    }
    {
        /* Controlled satellite leave-one-out tests; no default exclusion. */
        const char *excluded = getenv("RTK_PPP_EXSATS");
        if (excluded && *excluded) {
            char ids[256], *id;
            if (strlen(excluded) >= sizeof(ids)) {
                fprintf(stderr, "RTK_PPP_EXSATS is too long\n");
                return -1;
            }
            strcpy(ids, excluded);
            for (id = strtok(ids, ", "); id; id = strtok(NULL, ", ")) {
                int sat = satid2no(id);
                if (sat < 1 || sat > MAXSAT) {
                    fprintf(stderr, "Invalid RTK_PPP_EXSATS satellite: %s\n", id);
                    return -1;
                }
                prcopt.exsats[sat - 1] = 1;
            }
        }
    }

    {
        char config_error[256];
        if (!smartphone_ppp_validate(&prcopt, config_error, sizeof(config_error))) {
            fprintf(stderr, "PPP configuration error: %s\n", config_error);
            return -1;
        }
    }

    // ========================================================================
    // 4. 【执行解算】
    // ========================================================================
    printf("Starting Smartphone Un-combined PPP Processing...\n");
    {
        char trace_file[MAXSTRPATH] = "E:\\RTKLIB_Data\\ppp_debug.trace";
        const char *level = getenv("RTK_PPP_TRACE_LEVEL");
        if (!env_copy(trace_file, sizeof(trace_file), "RTK_PPP_TRACE")) return -1;
        traceopen(trace_file);
        tracelevel(level ? atoi(level) : 3);
    }
    smartphone_ppp_log("pc-postprocess", &prcopt);
    if (!pppvmf3load("E:\\RTKLIB_Data\\tro\\VMF3_20260317.H00",
                     "E:\\RTKLIB_Data\\tro\\VMF3_20260317.H06",
                      "E:\\RTKLIB_Data\\tro\\orography_ell_5x5")) {
        fprintf(stderr,
            "WARNING: VMF3 files were not loaded; using the internal troposphere model.\n");
    }
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
