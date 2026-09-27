/*------------------------------------------------------------------------------
* preceph.c : precise ephemeris and clock functions
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* references :
*     [1] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-c),
*         12 February, 2007
*     [2] J.Ray, W.Gurtner, RINEX Extensions to Handle Clock Information,
*         27 August, 1998
*     [3] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [4] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*     [5] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-d),
*         February 21, 2016
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2009/01/18 1.0  new
*           2009/01/31 1.1  fix bug on numerical error to read sp3a ephemeris
*           2009/05/15 1.2  support glonass,galileo,qzs
*           2009/12/11 1.3  support wild-card expansion of file path
*           2010/07/21 1.4  added api:
*                               eci2ecef(),sunmoonpos(),peph2pos(),satantoff(),
*                               readdcb()
*                           changed api:
*                               readsp3()
*                           deleted api:
*                               eph2posp()
*           2010/09/09 1.5  fix problem when precise clock outage
*           2011/01/23 1.6  support qzss satellite code
*           2011/09/12 1.7  fix problem on precise clock outage
*                           move sunmmonpos() to rtkcmn.c
*           2011/12/01 1.8  modify api readsp3()
*                           precede later ephemeris if ephemeris is NULL
*                           move eci2ecef() to rtkcmn.c
*           2013/05/08 1.9  fix bug on computing std-dev of precise clocks
*           2013/11/20 1.10 modify option for api readsp3()
*           2014/04/03 1.11 accept extension including sp3,eph,SP3,EPH
*           2014/05/23 1.12 add function to read sp3 velocity records
*                           change api: satantoff()
*           2014/08/31 1.13 add member cov and vco in peph_t sturct
*           2014/10/13 1.14 fix bug on clock error variance in peph2pos()
*           2015/05/10 1.15 add api readfcb()
*                           modify api readdcb()
*           2017/04/11 1.16 fix bug on antenna offset correction in peph2pos()
*           2020/11/30 1.17 support SP3-d [5] to accept more than 85 satellites
*                           support NavIC/IRNSS in API peph2pos()
*                           LC defined GPS/QZS L1-L2, GLO G1-G2, GAL E1-E5b,
*                            BDS B1I-B2I and IRN L5-S for API satantoff()
*                           fix bug on reading SP3 file extension
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))

#define NMAX        10              /* order of polynomial interpolation */
#define MAXDTE      900.0           /* max time difference to ephem time (s) */
#define EXTERR_CLK  1E-3            /* extrapolation error for clock (m/s) */
#define EXTERR_EPH  5E-7            /* extrapolation error for ephem (m/s^2) */
#define MAX_BIAS_SYS 4              /* # of constellations supported */

/* table to translate code to code bias table index  */
/* =============================================================================
 * 变量与函数：code_bias_ix & init_bias_ix
 * 功能：初始化卫星码偏差（Code Bias）查询索引表。
 * 解释：GNSS 卫星会发射多种信号（如 GPS 的 L1C, L1W, L2W 等），由于卫星内部硬件电路不同，
 * 不同信号到达天线发射口的时间有微小差异（DCB/ISC）。这个表将字符串形式的“信号代号”
 * 映射为简单的数字索引（0, 1, 2...），以便在内存中极速查找对应的硬件延迟改正值。
 * -----------------------------------------------------------------------------
 * -1 = 不支持的码类型
 * 0 = 选定的基准码（作为基准，自身偏差记为 0）
 * 1-3 = 其他码在偏差表中的对应索引
 * ============================================================================= */
static int8_t code_bias_ix[MAX_BIAS_SYS][MAXCODE+1];

static void init_bias_ix(void)
{
    int i, j;

    for (i = 0; i < MAX_BIAS_SYS; i++) {
        for (j = 0; j <= MAXCODE; j++) code_bias_ix[i][j] = -1;
    }

    /* GPS: slot 0=L1, slot 1=L2, slot 2=L5. */
    code_bias_ix[0][CODE_L1W] = 0;
    code_bias_ix[0][CODE_L1C] = 1;
    code_bias_ix[0][CODE_L1L] = 2;
    code_bias_ix[0][CODE_L1X] = 3;

    code_bias_ix[0][CODE_L2W] = 0;
    code_bias_ix[0][CODE_L2L] = 1;
    code_bias_ix[0][CODE_L2S] = 2;
    code_bias_ix[0][CODE_L2X] = 3;

    code_bias_ix[0][CODE_L5X] = 0;
    code_bias_ix[0][CODE_L5I] = 1;
    code_bias_ix[0][CODE_L5Q] = 2;

    /* GLONASS. */
    code_bias_ix[1][CODE_L1P] = 0;
    code_bias_ix[1][CODE_L1C] = 1;
    code_bias_ix[1][CODE_L2P] = 0;
    code_bias_ix[1][CODE_L2C] = 1;

    /* Galileo: E1 / E5b / E5a. */
    code_bias_ix[2][CODE_L1X] = 0;
    code_bias_ix[2][CODE_L1C] = 1;
    code_bias_ix[2][CODE_L1B] = 2;

    code_bias_ix[2][CODE_L7X] = 0;
    code_bias_ix[2][CODE_L7I] = 1;
    code_bias_ix[2][CODE_L7Q] = 2;

    code_bias_ix[2][CODE_L5X] = 0;
    code_bias_ix[2][CODE_L5I] = 1;
    code_bias_ix[2][CODE_L5Q] = 2;

    /* BeiDou: B1I / B2b(B2I) / B2a. */
    code_bias_ix[3][CODE_L2I] = 0;
    code_bias_ix[3][CODE_L2Q] = 1;
    code_bias_ix[3][CODE_L2X] = 2;

    code_bias_ix[3][CODE_L7I] = 0;
    code_bias_ix[3][CODE_L7Q] = 1;
    code_bias_ix[3][CODE_L7X] = 2;

    code_bias_ix[3][CODE_L5X] = 0;
    code_bias_ix[3][CODE_L5D] = 1;
    code_bias_ix[3][CODE_L5P] = 2;
}

/* =============================================================================
 * 函数：code2sys
 * 功能：将 SP3 文件中的卫星系统字符标识转换为 RTKLIB 内部的系统宏定义。
 * 解释：SP3 文件中卫星通常写为 "G01", "R05", "C08"。这里提取首字母进行归类。
 * ============================================================================= */
static int code2sys(char code)
{
    if (code == 'G' || code == ' ') return SYS_GPS; /* 美国 GPS */
    if (code == 'R') return SYS_GLO;            /* 俄罗斯 GLONASS */
    if (code == 'E') return SYS_GAL;            /* 欧洲 Galileo (SP3-d 标准) */
    if (code == 'J') return SYS_QZS;            /* 日本 QZSS (SP3-d 标准) */
    if (code == 'C') return SYS_CMP;            /* 中国 北斗 (SP3-d 标准) */
    if (code == 'I') return SYS_IRN;            /* 印度 NavIC (SP3-d 标准) */
    if (code == 'L') return SYS_LEO;            /* 低轨卫星 LEO */
    return SYS_NONE;
}

/* =============================================================================
 * 函数：readsp3h
 * 功能：读取并解析 SP3 精密星历文件的文件头（Header）。
 * 解释：提取 SP3 文件的起止时间、卫星总数 (ns)、时间坐标系以及用于推算
 * 坐标/钟差标准差的基准因子 (bfact)。它是正确读取后续数据块的前提。
 * ============================================================================= */
static int readsp3h(FILE* fp, gtime_t* time, char* type, int* sats,
    double* bfact, char* tsys)
{
    int i = 0, j, k = 0, ns = 0, nl = 5, sys, prn;
    char buff[1024];

    trace(3, "readsp3h:\n");

    /* 逐行读取文件头 */
    while (fgets(buff, sizeof(buff), fp)) {

        /* 第一行：解析 SP3 版本 (c 或 d) 以及参考起始历元时间 */
        if (buff[0] == '#' && (buff[1] == 'c' || buff[1] == 'd')) {
            *type = buff[2];
            if (str2time(buff, 3, 28, time)) return 0;
        }
        /* '+' 号开头的行：解析该 SP3 文件中包含的所有卫星的 PRN 号 */
        else if (buff[0] == '+' && buff[1] == ' ') {
            if (i == 2) {
                ns = (int)str2num(buff, 3, 3); /* 提取卫星总数 */
                if (ns > 85) nl = ns / 17 + (ns % 17 != 0); /* 根据总数计算需要读取几行 '+' */
            }
            for (j = 0; j < 17 && k < ns; j++) {
                sys = code2sys(buff[9 + 3 * j]);
                prn = (int)str2num(buff, 10 + 3 * j, 2);
                if (k < MAXSAT) sats[k++] = satno(sys, prn); /* 将系统+编号映射为内部卫星 ID */
            }
        }
        /* '%c' 行：解析时间系统 (通常为 GPS 或 UTC) */
        else if (i == 2 * nl + 2) {
            memcpy(tsys, buff + 9, 3); tsys[3] = '\0';
        }
        /* '%f' 行：解析位置和钟差精度的 Base Factor (基底因子) */
        else if (i == 2 * nl + 4) {
            bfact[0] = str2num(buff, 3, 10);
            bfact[1] = str2num(buff, 14, 12);
        }
        /* '%i' 行：读取到最后一行注释头，跳出循环 */
        else if (i == 2 * nl + 11) {
            break; /* at end of header */
        }
        i = i + 1; /* line counter */
    }
    return ns; /* 返回解析出的卫星总数 */
}

/* =============================================================================
 * 函数：addpeph
 * 功能：向导航数据结构（nav_t）中动态添加一条精密星历记录。
 * 解释：由于不确定 SP3 文件到底包含多少个历元，这里采用动态扩容池的策略。
 * 每次内存不足时，将数组长度增加 256，防止内存溢出 (OOM)。
 * ============================================================================= */
static int addpeph(nav_t* nav, peph_t* peph)
{
    peph_t* nav_peph;

    /* 如果当前记录数达到了最大容量，触发扩容 */
    if (nav->ne >= nav->nemax) {
        nav->nemax += 256;
        if (!(nav_peph = (peph_t*)realloc(nav->peph, sizeof(peph_t) * nav->nemax))) {
            trace(1, "readsp3b malloc error n=%d\n", nav->nemax);
            free(nav->peph); nav->peph = NULL; nav->ne = nav->nemax = 0;
            return 0;
        }
        nav->peph = nav_peph;
    }
    /* 将传入的新星历记录存入数组末尾 */
    nav->peph[nav->ne++] = *peph;
    return 1;
}

/* =============================================================================
 * 函数：readsp3b
 * 功能：读取 SP3 文件的核心数据主体（Body），获取精密坐标、速度和钟差。
 * 核心：★ 无 fseek 的前瞻读取状态机补丁 ★
 * 解释：通过一次读取并缓存下一行，解决了传统使用 fseek 回退文件指针导致的兼容性
 * 和死循环问题。极大地提高了对于缺失卫星行、空行等劣质格式 SP3 文件的容错率。
 * ============================================================================= */
static void readsp3b(FILE* fp, char type, int* sats, int ns, double* bfact,
    char* tsys, int index, int opt, nav_t* nav)
{
    peph_t peph;
    gtime_t time;
    double val, std, base;
    int i, j, sat, sys, prn, pred_o, pred_c, v;
    char buff[1024];

    trace(3, "readsp3b: type=%c ns=%d index=%d opt=%d\n", type, ns, index, opt);

    /* 先瞻读取文件的第一行主体数据 */
    if (!fgets(buff, sizeof(buff), fp)) return;

    while (!feof(fp)) {

        if (!strncmp(buff, "EOF", 3)) break; /* 到达文件真实结尾 */

        /* 历元时间行以 '*' 开头 */
        if (buff[0] != '*' || str2time(buff, 3, 28, &time)) {
            trace(2, "sp3 invalid epoch %31.31s\n", buff);
            if (!fgets(buff, sizeof(buff), fp)) break;
            continue;
        }

        /* 如果时间是 UTC，转换到 GPS 时间 */
        if (!strcmp(tsys, "UTC")) time = utc2gpst(time);
        peph.time = time;
        peph.index = index;

        /* 初始化当前历元所有卫星的坐标/速度矩阵为 0 */
        for (i = 0; i < MAXSAT; i++) {
            for (j = 0; j < 4; j++) {
                peph.pos[i][j] = 0.0;
                peph.std[i][j] = 0.0f;
                peph.vel[i][j] = 0.0;
                peph.vst[i][j] = 0.0f;
            }
            for (j = 0; j < 3; j++) {
                peph.cov[i][j] = 0.0f;
                peph.vco[i][j] = 0.0f;
            }
        }

        /* 初始状态重置 */
        v = 0;
        pred_o = pred_c = 0;

        /* ===================================================================
           🌟 终极修复：废弃依据 ns(文件头声明卫星数) 的内层 for 循环计数限制。
           改为“边界符驱动”的 while(1) 动态读取，完美兼容实际卫星数超标的残缺星历。
           =================================================================== */
        while (1) {
            if (!fgets(buff, sizeof(buff), fp)) break;
            if (!strncmp(buff, "EOF", 3)) break;

            /* 【绝杀补丁】：碰到下一个历元的 '*'，代表当前历元的数据确确实实读完了！
               立刻退出内层读取，并把装有 '*' 的 buff 平滑移交给外层 while 处理。 */
            if (buff[0] == '*') {
                break;
            }

            /* 跳过无效空行或错误格式 */
            if (strlen(buff) < 4 || (buff[0] != 'P' && buff[0] != 'V')) {
                continue;
            }

            /* 解析卫星系统和 PRN 号 */
            sys = buff[1] == ' ' ? SYS_GPS : code2sys(buff[1]);
            prn = (int)str2num(buff, 2, 2);
            if (sys == SYS_SBS) prn += 100;
            else if (sys == SYS_QZS) prn += 192; /* extension to sp3-c */

            /* 如果是未开启的星座，跳过该行，继续读下一行 */
            if (!(sat = satno(sys, prn))) {
                continue;
            }

            /* 判定预测标志位 */
            if (buff[0] == 'P') {
                pred_c = strlen(buff) >= 76 && buff[75] == 'P';
                pred_o = strlen(buff) >= 80 && buff[79] == 'P';
            }

            /* 解析坐标/速度的 X, Y, Z 以及 钟差 (共 4 维) */
            for (j = 0; j < 4; j++) {
                if (j < 3 && (opt & 1) && pred_o) continue;
                if (j < 3 && (opt & 2) && !pred_o) continue;
                if (j == 3 && (opt & 1) && pred_c) continue;
                if (j == 3 && (opt & 2) && !pred_c) continue;

                val = str2num(buff, 4 + j * 14, 14);
                std = str2num(buff, 61 + j * 3, j < 3 ? 2 : 3);

                if (buff[0] == 'P') { /* 读取位置 (Position) 和 钟差 */
                    if (val != 0.0 && fabs(val - 999999.999999) >= 1E-6) {
                        peph.pos[sat - 1][j] = val * (j < 3 ? 1000.0 : 1E-6);
                        v = 1; /* 标记该历元存在有效数据 */
                    }
                    if ((base = bfact[j < 3 ? 0 : 1]) > 0.0 && std > 0.0) {
                        peph.std[sat - 1][j] = (float)(pow(base, std) * (j < 3 ? 1E-3 : 1E-12));
                    }
                }
                else if (v) { /* 读取速度 (Velocity) 和 钟速 */
                    if (val != 0.0 && fabs(val - 999999.999999) >= 1E-6) {
                        peph.vel[sat - 1][j] = val * (j < 3 ? 0.1 : 1E-10);
                    }
                    if ((base = bfact[j < 3 ? 0 : 1]) > 0.0 && std > 0.0) {
                        peph.vst[sat - 1][j] = (float)(pow(base, std) * (j < 3 ? 1E-7 : 1E-16));
                    }
                }
            }
        }

        /* 如果该历元包含有效卫星，存入内部结构体 */
        if (v) {
            if (!addpeph(nav, &peph)) return;
        }
    }
}

/* =============================================================================
 * 函数：cmppeph
 * 功能：供 qsort (快速排序) 调用的比较函数，按时间顺序排列精密星历历元。
 * ============================================================================= */
static int cmppeph(const void* p1, const void* p2)
{
    peph_t* q1 = (peph_t*)p1, * q2 = (peph_t*)p2;
    double tt = timediff(q1->time, q2->time);
    return tt < -1E-9 ? -1 : (tt > 1E-9 ? 1 : q1->index - q2->index);
}

/* =============================================================================
 * 函数：combpeph
 * 功能：合并内存中的精密星历。
 * 解释：当用户提供多个 SP3 文件时，往往存在时间重叠或重复历元。此函数会将
 * nav->peph 数组按照时间排序，并将时刻相同、但信息互补的卫星数据拼接到同一个历元中。
 * ============================================================================= */
static void combpeph(nav_t* nav, int opt)
{
    int i, j, k, m;

    trace(3, "combpeph: ne=%d\n", nav->ne);

    /* 1. 按时间升序快速排序 */
    qsort(nav->peph, nav->ne, sizeof(peph_t), cmppeph);

    if (opt & 4) return;

    /* 2. 遍历合并相同时间（时间差小于 1 纳秒）的重复历元 */
    for (i = 0, j = 1; j < nav->ne; j++) {

        if (fabs(timediff(nav->peph[i].time, nav->peph[j].time)) < 1E-9) {

            /* 相同时间下，将记录 j 中有数据的卫星，填补到记录 i 的空缺中 */
            for (k = 0; k < MAXSAT; k++) {
                if (norm(nav->peph[j].pos[k], 4) <= 0.0) continue;
                for (m = 0; m < 4; m++) nav->peph[i].pos[k][m] = nav->peph[j].pos[k][m];
                for (m = 0; m < 4; m++) nav->peph[i].std[k][m] = nav->peph[j].std[k][m];
                for (m = 0; m < 4; m++) nav->peph[i].vel[k][m] = nav->peph[j].vel[k][m];
                for (m = 0; m < 4; m++) nav->peph[i].vst[k][m] = nav->peph[j].vst[k][m];
            }
        }
        /* 如果时间不同，作为独立新历元保留 */
        else if (++i < j) nav->peph[i] = nav->peph[j];
    }
    nav->ne = i + 1; /* 更新去重后的实际历元总数 */

    trace(4, "combpeph: ne=%d\n", nav->ne);
}
/* =============================================================================
 * 函数：readsp3
 * 功能：精密星历 (SP3) 文件读取总调度器。
 * 解释：支持带有通配符（如 IGS000*.sp3）的文件路径，自动展开并逐个调用
 * readsp3h (读头) 和 readsp3b (读体)，最后调用 combpeph 进行多文件合并去重。
 * ============================================================================= */
extern void readsp3(const char* file, nav_t* nav, int opt)
{
    FILE* fp;
    gtime_t time = { 0 };
    double bfact[2] = { 0 };
    int i, j, n, ns, sats[MAXSAT] = { 0 };
    char* efiles[MAXEXFILE], * ext, type = ' ', tsys[4] = "";

    trace(3, "readpephs: file=%s\n", file);

    /* 为展开通配符后的多个文件路径分配内存 */
    for (i = 0; i < MAXEXFILE; i++) {
        if (!(efiles[i] = (char*)malloc(1024))) {
            for (i--; i >= 0; i--) free(efiles[i]);
            return;
        }
    }
    /* 解析通配符，获取匹配到的实际文件列表和数量 n */
    n = expath(file, efiles, MAXEXFILE);

    for (i = j = 0; i < n; i++) {
        if (!(ext = strrchr(efiles[i], '.'))) continue;

        /* 仅处理扩展名为 .sp3 或 .eph 的文件 */
        if (!strstr(ext, ".sp3") && !strstr(ext, ".SP3") &&
            !strstr(ext, ".eph") && !strstr(ext, ".EPH")) continue;

        if (!(fp = fopen(efiles[i], "r"))) {
            trace(2, "sp3 file open error %s\n", efiles[i]);
            continue;
        }
        /* 依次调用 readsp3h(读文件头) 和 readsp3b(读主体) */
        ns = readsp3h(fp, &time, &type, sats, bfact, tsys);
        readsp3b(fp, type, sats, ns, bfact, tsys, j++, opt, nav);

        fclose(fp);
    }
    for (i = 0; i < MAXEXFILE; i++) free(efiles[i]);

    /* 如果读取到了有效数据，执行排序和合并 */
    if (nav->ne > 0) combpeph(nav, opt);
}

/* =============================================================================
 * 函数：readsap (Read Satellite Antenna Parameters)
 * 功能：读取卫星天线相位中心参数文件 (通常为 .atx 格式)。
 * 解释：精密星历给出的卫星坐标通常是卫星质心 (CoM)，而 GPS 信号是从天线发射出去的。
 * 我们必须根据高度角和方位角，利用 ANTEX 文件中的修正参数将质心改正到天线相位中心 (APC)。
 * ============================================================================= */
extern int readsap(const char* file, gtime_t time, nav_t* nav)
{
    pcvs_t pcvs = { 0 };
    pcv_t pcv0 = { 0 }, * pcv;
    int i;

    char tstr[40];
    trace(3, "readsap : file=%s time=%s\n", file, time2str(time, tstr, 0));

    /* 解析 ANTEX 文件内容到临时结构体 pcvs */
    if (!readpcv(file, &pcvs)) return 0;

    /* 遍历所有可能的卫星，查找在当前观测时间 time 下有效的 APC 参数 */
    for (i = 0; i < MAXSAT; i++) {
        pcv = searchpcv(i + 1, "", time, &pcvs);
        nav->pcvs[i] = pcv ? *pcv : pcv0; /* 存入导航数据核心结构 nav->pcvs 中 */
    }
    free(pcvs.pcv);
    return 1;
}

/* =============================================================================
 * 函数：readdcbf
 * 功能：读取老格式的差分码偏差 (DCB) 文件。
 * 解释：主要处理类似于 "P1-C1" (P码与C码的偏差) 或 "P2-C2" 的传统 DCB 文件。
 * ============================================================================= */
static int readdcbf(const char* file, nav_t* nav, const sta_t* sta)
{
    FILE* fp;
    double cbias;
    char buff[256], str1[32], str2[32] = "";
    int i, j, sat, type = 0;

    trace(3, "readdcbf: file=%s\n", file);

    if (!(fp = fopen(file, "r"))) {
        trace(2, "dcb parameters file open error: %s\n", file);
        return 0;
    }
    while (fgets(buff, sizeof(buff), fp)) {

        /* 识别当前读取的行属于哪种频间/码间偏差 */
        if (strstr(buff, "DIFFERENTIAL (P1-C1) CODE BIASES")) type = 1;
        else if (strstr(buff, "DIFFERENTIAL (P2-C2) CODE BIASES")) type = 2;

        str2[0] = '\0';
        if (!type || sscanf(buff, "%31s %31s", str1, str2) < 1) continue;

        /* 解析出具体的纳秒级偏差值 cbias */
        if ((cbias = str2num(buff, 26, 9)) == 0.0) continue;

        /* 判断该项是接收机(Receiver)端偏差还是卫星(Satellite)端偏差 */
        if (sta && (!strcmp(str1, "G") || !strcmp(str1, "R"))) { /* receiver DCB */
            for (i = 0; i < MAXRCV; i++) {
                if (!strcmp(sta[i].name, str2)) break; /* 匹配测站接收机名 */
            }
            if (i < MAXRCV) {
                j = !strcmp(str1, "G") ? 0 : 1;
                /* 核心转换：将文件中的纳秒 (ns) 乘以光速 (CLIGHT) 转化为米 (m) */
                nav->rbias[i][j][type - 1] = cbias * 1E-9 * CLIGHT;
            }
        }
        else if ((sat = satid2no(str1))) { /* satellite dcb (通过类似 G01 提取卫星号) */
            nav->cbias[sat - 1][type - 1][0] = cbias * 1E-9 * CLIGHT; /* 同样转为距离 米 */
        }
    }
    fclose(fp);

    return 1;
}

/* =============================================================================
 * 函数：sys2ix
 * 功能：将系统宏转换为索引。
 * ============================================================================= */
 /* satellite system to index */
static int sys2ix(int sys)
{
    switch (sys) {
    case SYS_GPS: return 0;
    case SYS_SBS: return 0;
    case SYS_GLO: return 1;
    case SYS_GAL: return 2;
    case SYS_CMP: return 3;
    case SYS_QZS: return 4;
    case SYS_IRN: return 5;
    }
    return -1;
}

/* =============================================================================
 * 函数：code2bias_ix
 * 功能：安全地从全局表格 code_bias_ix 中查询指定系统和观测码类型的偏差表索引号。
 * ============================================================================= */
extern int code2bias_ix(int sys, int code)
{
    int sys_ix = sys2ix(sys);

    if (sys_ix < 0 || sys_ix >= MAX_BIAS_SYS) return -1;
    if (code <= CODE_NONE || code > MAXCODE) return -1;
    return code_bias_ix[sys_ix][code];
}

/* =============================================================================
 * 函数：readbiaf
 * 功能：读取现代格式的 SINEX Bias (BIA 或 BSX) 文件。
 * 解释：现代多频多系统定位采用更严密的绝对偏差(OSB)和差分偏差(DSB)表示法。
 * 这个函数将文件中每种信号的纳秒级偏差提取出来，并转换为统一参考基准下的距离改正(米)。
 * ============================================================================= */
/* Bias-SINEX times in the supported GPS time system: YYYY:DOY:SSSSS. */
static int bia_time(const char *s, gtime_t *time)
{
    double ep[6] = {0};
    int year, doy, second;
    char extra;

    if (sscanf(s, "%d:%d:%d%c", &year, &doy, &second, &extra) != 3 ||
        year < 1980 || year > 2199 || doy < 1 || doy > 366 ||
        second < 0 || second > 86400) return 0;
    ep[0] = year; ep[1] = 1; ep[2] = 1;
    *time = timeadd(epoch2time(ep), (doy - 1) * 86400.0 + second);
    return 1;
}

static int add_codeosb(nav_t *nav, int sat, int code, gtime_t start,
                       gtime_t end, double bias)
{
    codeosb_t *records;
    int nmax;

    if (nav->nosb >= nav->nosbmax) {
        nmax = nav->nosbmax ? nav->nosbmax * 2 : 1024;
        if (!(records = (codeosb_t *)realloc(nav->osbs,
                                              sizeof(codeosb_t) * nmax))) return 0;
        nav->osbs = records;
        nav->nosbmax = nmax;
    }
    records = nav->osbs + nav->nosb;
    records->start = start;
    records->end = end;
    records->bias = bias;
    records->next = nav->osb_head[sat - 1][code];
    nav->osb_head[sat - 1][code] = ++nav->nosb;
    return 1;
}

/* Return an OSB only inside its declared half-open validity interval. */
extern int codeosb_at(const nav_t *nav, gtime_t time, int sat, int code,
                      double *bias)
{
    int index, found = 0;
    const codeosb_t *record;

    if (!nav || sat < 1 || sat > MAXSAT || code < 1 || code > MAXCODE) return 0;
    for (index = nav->osb_head[sat - 1][code]; index;
         index = record->next) {
        record = nav->osbs + index - 1;
        found = 1;
        if (timediff(time, record->start) >= 0.0 &&
            timediff(time, record->end) < 0.0) {
            if (bias) *bias = record->bias;
            return 1;
        }
    }
    return found ? -1 : 0;
}

static int readbiaf(const char* file, nav_t* nav)
{
    FILE* fp;
    char buff[512], work[512], *tok[24], *p;
    int nt, t, sat, sys, code1, code2, freq1, freq2;
    int bias_ix1, bias_ix2, nread = 0;
    double cbias;
    gtime_t start, end;

    trace(3, "readbiaf: file=%s\n", file);

    if (!(fp = fopen(file, "r"))) {
        trace(2, "bias file open error: %s\n", file);
        return 0;
    }

    while (fgets(buff, sizeof(buff), fp)) {
        char *kind = NULL, *obs1 = NULL, *obs2 = NULL, *prn = NULL;
        char *unit = NULL, *estimate = NULL;
        int have_start = 0, have_end = 0;

        if (strstr(buff, "TIME_SYSTEM") && !strstr(buff, "G")) {
            trace(1, "readbiaf: unsupported time system in %s\n", file);
            fclose(fp);
            return 0;
        }

        strncpy(work, buff, sizeof(work) - 1);
        work[sizeof(work) - 1] = '\0';

        nt = 0;
        for (p = strtok(work, " \t\r\n"); p && nt < (int)(sizeof(tok) / sizeof(tok[0]));
             p = strtok(NULL, " \t\r\n")) {
            tok[nt++] = p;
        }
        if (nt < 2) continue;

        kind = tok[0];
        if (strcmp(kind, "OSB") && strcmp(kind, "DSB")) continue;

        /*
         * Satellite SINEX-BIAS OSB/DSB records are "TYPE SVN PRN ...".
         * Do NOT scan from the SVN field: after enabling BDS C59-C63 an SVN
         * such as C061 is itself a syntactically valid PRN and would send C02
         * biases to C61.  Token 2 is the actual PRN for satellite records.
         */
        if (nt < 3) continue;
        prn = tok[2];
        sat = satid2no(prn);
        if (!sat) continue;

        sys = satsys(sat, NULL);
        if (sys2ix(sys) < 0) continue;

        /* Find one/two code observables such as C1C, C5Q, C7X.
         * Be strict: BeiDou SVN/PRN tokens also start with 'C' (C080/C01).
         */
        for (t = 1; t < nt; t++) {
            size_t len = strlen(tok[t]);
            if (len != 3 || tok[t][0] != 'C' ||
                tok[t][1] < '0' || tok[t][1] > '9' ||
                !((tok[t][2] >= 'A' && tok[t][2] <= 'Z') ||
                  (tok[t][2] >= 'a' && tok[t][2] <= 'z'))) continue;
            if (!obs1) obs1 = tok[t];
            else if (!obs2) {
                obs2 = tok[t];
                break;
            }
        }
        if (!obs1) continue;

        for (t = 0; t < nt; t++) {
            gtime_t epoch;
            if (!bia_time(tok[t], &epoch)) continue;
            if (!have_start) { start = epoch; have_start = 1; }
            else { end = epoch; have_end = 1; break; }
        }

        /* Find the estimate after the ns unit token. */
        for (t = 0; t + 1 < nt; t++) {
            if (!strcmp(tok[t], "ns") || !strcmp(tok[t], "NS")) {
                unit = tok[t];
                estimate = tok[t + 1];
                break;
            }
        }
        if (!unit || !estimate) {
            /* Fallback for products matching the traditional fixed columns. */
            cbias = str2num(buff, 70, 21);
        }
        else {
            cbias = atof(estimate);
        }
        if (!isfinite(cbias)) continue;

        code1 = obs2code(obs1 + 1);
        if (code1 == CODE_NONE) continue;

        /*
         * Bias-SINEX OSB sign convention is
         *     unbiased observable = observed observable - OSB.
         * Store OSB by exact observation code before native-frequency mapping:
         * this also preserves modern signals (e.g. BDS B1C) that are outside
         * the classic NFREQ=3 slots but may be used by a future Android adapter.
         */
        if (!strcmp(kind, "OSB")) {
            /* Never turn a time-bounded product into a timeless correction. */
            if (!have_start || !have_end || timediff(end, start) <= 0.0) continue;
            if (!add_codeosb(nav, sat, code1, start, end,
                             cbias * 1E-9 * CLIGHT)) {
                trace(1, "readbiaf: out of memory for code OSB\n");
                fclose(fp);
                return 0;
            }
            nav->osb[sat - 1][code1] = cbias * 1E-9 * CLIGHT;
            nav->osb_valid[sat - 1][code1] = 1;
            nread++;
            continue;
        }

        /* Legacy DSB/DCB representation is frequency-slot based. */
        freq1 = code2idx(sys, (uint8_t)code1);
        if (freq1 < 0 || freq1 >= MAX_CODE_BIAS_FREQS) continue;

        bias_ix1 = code2bias_ix(sys, code1);
        if (bias_ix1 < 0 || bias_ix1 > MAX_CODE_BIASES) continue;

        /* DSB: only same-frequency code pairs fit this cbias structure. */
        if (!obs2) continue;
        code2 = obs2code(obs2 + 1);
        if (code2 == CODE_NONE) continue;
        freq2 = code2idx(sys, (uint8_t)code2);
        if (freq2 != freq1 || freq2 < 0 || freq2 >= MAX_CODE_BIAS_FREQS) continue;

        bias_ix2 = code2bias_ix(sys, code2);
        if (bias_ix2 < 0 || bias_ix2 > MAX_CODE_BIASES) continue;

        if (bias_ix1 == 0 && bias_ix2 > 0) {
            nav->cbias[sat - 1][freq1][bias_ix2 - 1] =
                cbias * 1E-9 * CLIGHT;
            nread++;
        }
        else if (bias_ix2 == 0 && bias_ix1 > 0) {
            nav->cbias[sat - 1][freq1][bias_ix1 - 1] =
                -cbias * 1E-9 * CLIGHT;
            nread++;
        }
    }
    fclose(fp);

    trace(2, "readbiaf: loaded %d code-bias records from %s\n", nread, file);
    return nread > 0;
}

/* =============================================================================
 * 函数：readdcb
 * 功能：对外暴露的 硬件延迟偏差(DCB) 文件加载主入口！
 * 解释：支持读取旧版 .DCB 格式以及新版 .BIA / .BSX 格式。它会自动调用 init_bias_ix
 * 建立查询表，展开文件通配符，并自动根据后缀名分发给 readdcbf 或 readbiaf 去处理。
 * ============================================================================= */
extern int readdcb(const char* file, nav_t* nav, const sta_t* sta)
{
    int i, j, k, n, dcb_ok = 0;
    char* efiles[MAXEXFILE] = { 0 };

    trace(3, "readdcb : file=%s\n", file);

    /* 调用上面看到的函数，初始化 测距码-表索引 查询表 */
    init_bias_ix();

    /* 安全起见，清空导航结构体中原有的所有码偏差数据 */
    for (i = 0; i < MAXSAT; i++) {
        int code;
        for (j = 0; j < MAX_CODE_BIAS_FREQS; j++) for (k = 0; k < MAX_CODE_BIASES; k++) {
            nav->cbias[i][j][k] = 0.0;
        }
        for (code = 0; code <= MAXCODE; code++) {
            nav->osb[i][code] = 0.0;
            nav->osb_valid[i][code] = 0;
            nav->osb_head[i][code] = 0;
        }
    }
    free(nav->osbs); nav->osbs = NULL; nav->nosb = nav->nosbmax = 0;

    /* 为展开多个文件分配内存 */
    for (i = 0; i < MAXEXFILE; i++) {
        if (!(efiles[i] = (char*)malloc(1024))) {
            for (i--; i >= 0; i--) free(efiles[i]);
            return 0;
        }
    }
    /* 展开文件通配符 (例如 /data/CODE*.BIA) */
    n = expath(file, efiles, MAXEXFILE);

    /* 遍历匹配到的文件列表 */
    for (i = 0; i < n; i++) {
        /* 根据不同的文件扩展名路由给相应的解析器 */
        if (strstr(efiles[i], ".BIA") || strstr(efiles[i], ".bia") ||
            strstr(efiles[i], ".BSX") || strstr(efiles[i], ".bsx"))
            dcb_ok = readbiaf(efiles[i], nav); /* 新版 SINEX 格式 */
        else if (strstr(efiles[i], ".DCB") || strstr(efiles[i], ".dcb"))
            dcb_ok = readdcbf(efiles[i], nav, sta); /* 旧版 DCB 格式 */
    }

    /* 释放路径内存 */
    for (i = 0; i < MAXEXFILE; i++) free(efiles[i]);

    return dcb_ok; /* 如果任意一个文件被成功读取，返回 1 */
}
/* =============================================================================
 * 功能：基于 Neville 算法的多项式插值
 * 说明：用于从离散的精密星历数据点中插值计算出任意时刻的卫星位置。
 * ============================================================================= */
static double interppol(const double* x, double* y, int n)
{
    int i, j;

    for (j = 1; j < n; j++) {
        for (i = 0; i < n - j; i++) {
            y[i] = (x[i + j] * y[i] - x[i] * y[i + 1]) / (x[i + j] - x[i]);
        }
    }
    return y[0];
}

/* =============================================================================
 * 功能：利用精密星历计算卫星坐标 (XYZ) 和轨道误差方差
 * 说明：通过二分查找确定时间窗口，利用 Neville 多项式插值计算卫星位置，
 * 并加入地球自转修正和外推误差补偿。
 * ============================================================================= */
static int pephpos(gtime_t time, int sat, const nav_t* nav, double* rs,
    double* dts, double* vare, double* varc)
{
    double t[NMAX + 1], p[3][NMAX + 1], c[2], * pos, std = 0.0, s[3], sinl, cosl;
    int i, j, k, index;

    char tstr[40];
    trace(4, "pephpos : time=%s sat=%2d\n", time2str(time, tstr, 3), sat);

    rs[0] = rs[1] = rs[2] = dts[0] = 0.0;

    /* 边界检查：确保当前时间在精密星历的覆盖范围内 */
    if (nav->ne < NMAX + 1 ||
        timediff(time, nav->peph[0].time) < -MAXDTE ||
        timediff(time, nav->peph[nav->ne - 1].time) > MAXDTE) {
        trace(3, "no prec ephem %s sat=%2d\n", time2str(time, tstr, 0), sat);
        return 0;
    }
    /* 二分查找找到时间对应的最近邻数据点索引 */
    for (i = 0, j = nav->ne - 1; i < j;) {
        k = (i + j) / 2;
        if (timediff(nav->peph[k].time, time) < 0.0) i = k + 1; else j = k;
    }
    index = i <= 0 ? 0 : i - 1;

    /* 确定插值区间，确保区间内有足够的点进行多项式插值 */
    i = index - (NMAX + 1) / 2;
    if (i < 0) i = 0; else if (i + NMAX >= nav->ne) i = nav->ne - NMAX - 1;

    for (j = 0; j <= NMAX; j++) {
        t[j] = timediff(nav->peph[i + j].time, time);
        if (norm(nav->peph[i + j].pos[sat - 1], 3) <= 0.0) {
            trace(3, "prec ephem outage %s sat=%2d\n", time2str(time, tstr, 0), sat);
            return 0;
        }
    }
    /* 坐标系转换：在插值前对轨道点应用地球自转补偿 */
    for (j = 0; j <= NMAX; j++) {
        pos = nav->peph[i + j].pos[sat - 1];
        sinl = sin(OMGE * t[j]);
        cosl = cos(OMGE * t[j]);
        p[0][j] = cosl * pos[0] - sinl * pos[1];
        p[1][j] = sinl * pos[0] + cosl * pos[1];
        p[2][j] = pos[2];
    }
    /* 执行 Neville 插值得到三维坐标 */
    for (i = 0; i < 3; i++) {
        rs[i] = interppol(t, p[i], NMAX + 1);
    }
    /* 轨道误差处理：计算位置方差并加入外推引起的误差增长 */
    if (vare) {
        for (i = 0; i < 3; i++) s[i] = nav->peph[index].std[sat - 1][i];
        std = norm(s, 3);

        if (t[0] > 0.0) std += EXTERR_EPH * SQR(t[0]) / 2.0;
        else if (t[NMAX] < 0.0) std += EXTERR_EPH * SQR(t[NMAX]) / 2.0;
        *vare = SQR(std);
    }
    /* 钟差处理：对精密钟差进行线性插值 */
    t[0] = timediff(time, nav->peph[index].time);
    t[1] = timediff(time, nav->peph[index + 1].time);
    c[0] = nav->peph[index].pos[sat - 1][3];
    c[1] = nav->peph[index + 1].pos[sat - 1][3];

    if (t[0] <= 0.0) {
        if ((dts[0] = c[0]) != 0.0) {
            std = nav->peph[index].std[sat - 1][3] * CLIGHT - EXTERR_CLK * t[0];
        }
    }
    else if (t[1] >= 0.0) {
        if ((dts[0] = c[1]) != 0.0) {
            std = nav->peph[index + 1].std[sat - 1][3] * CLIGHT + EXTERR_CLK * t[1];
        }
    }
    else if (c[0] != 0.0 && c[1] != 0.0) {
        dts[0] = (c[1] * t[0] - c[0] * t[1]) / (t[0] - t[1]);
        i = t[0] < -t[1] ? 0 : 1;
        std = nav->peph[index + i].std[sat - 1][3] * CLIGHT + EXTERR_CLK * fabs(t[i]);
    }
    else {
        dts[0] = 0.0;
    }
    if (varc) *varc = SQR(std);
    return 1;
}

/* =============================================================================
 * 功能：利用精密钟差 (PCLK) 文件计算卫星钟差
 * 说明：当存在比 SP3 文件更高精度的 CLK 文件时，优先使用此逻辑进行线性插值。
 * ============================================================================= */
static int pephclk(gtime_t time, int sat, const nav_t* nav, double* dts,
    double* varc)
{
    double t[2], c[2], std;
    int i, j, k, index;

    char tstr[40];
    trace(4, "pephclk : time=%s sat=%2d\n", time2str(time, tstr, 3), sat);

    if (nav->nc < 2 ||
        timediff(time, nav->pclk[0].time) < -MAXDTE ||
        timediff(time, nav->pclk[nav->nc - 1].time) > MAXDTE) {
        trace(3, "no prec clock %s sat=%2d\n", time2str(time, tstr, 0), sat);
        return 1;
    }
    /* 二分查找 */
    for (i = 0, j = nav->nc - 1; i < j;) {
        k = (i + j) / 2;
        if (timediff(nav->pclk[k].time, time) < 0.0) i = k + 1; else j = k;
    }
    index = i <= 0 ? 0 : i - 1;

    /* 精密钟差的线性插值 */
    t[0] = timediff(time, nav->pclk[index].time);
    t[1] = timediff(time, nav->pclk[index + 1].time);
    c[0] = nav->pclk[index].clk[sat - 1][0];
    c[1] = nav->pclk[index + 1].clk[sat - 1][0];

    if (t[0] <= 0.0) {
        if ((dts[0] = c[0]) == 0.0) return 0;
        std = nav->pclk[index].std[sat - 1][0] * CLIGHT - EXTERR_CLK * t[0];
    }
    else if (t[1] >= 0.0) {
        if ((dts[0] = c[1]) == 0.0) return 0;
        std = nav->pclk[index + 1].std[sat - 1][0] * CLIGHT + EXTERR_CLK * t[1];
    }
    else if (c[0] != 0.0 && c[1] != 0.0) {
        dts[0] = (c[1] * t[0] - c[0] * t[1]) / (t[0] - t[1]);
        i = t[0] < -t[1] ? 0 : 1;
        std = nav->pclk[index + i].std[sat - 1][0] * CLIGHT + EXTERR_CLK * fabs(t[i]);
    }
    else {
        trace(3, "prec clock outage %s sat=%2d\n", time2str(time, tstr, 0), sat);
        return 0;
    }
    if (varc) *varc = SQR(std);
    return 1;
}

/* =============================================================================
 * 功能：卫星天线相位中心偏移 (PCO) 改正
 * 说明：将卫星质心坐标改正为天线发射信号的相位中心坐标。
 * 这是精密定位必须的步骤，因为它处理了不同频率之间的频率相关偏差 (Phase Center Offset)。
 * ============================================================================= */
extern void satantoff(gtime_t time, const double* rs, int sat, const nav_t* nav,
    double* dant)
{
    const pcv_t* pcv = nav->pcvs + sat - 1;
    double ex[3], ey[3], ez[3], es[3], r[3], rsun[3], gmst, erpv[5] = { 0 }, freq[2];
    double C1, C2, dant1, dant2;
    int i, sys;

    char tstr[40];
    trace(4, "satantoff: time=%s sat=%2d\n", time2str(time, tstr, 3), sat);

    dant[0] = dant[1] = dant[2] = 0.0;

    /* 计算太阳在 ECEF 系中的位置，用于卫星姿态判断 */
    sunmoonpos(gpst2utc(time), erpv, rsun, NULL, &gmst);

    /* 建立卫星本体坐标系 (X, Y, Z) */
    for (i = 0; i < 3; i++) r[i] = -rs[i];
    if (!normv3(r, ez)) return;
    for (i = 0; i < 3; i++) r[i] = rsun[i] - rs[i];
    if (!normv3(r, es)) return;
    cross3(ez, es, r);
    if (!normv3(r, ey)) return;
    cross3(ey, ez, ex);

    /* 选择无电离层组合的基准频率，用于进行 LC 组合偏移改正 */
    sys = satsys(sat, NULL);
    if (sys == SYS_GPS || sys == SYS_QZS) { freq[0] = FREQL1; freq[1] = FREQL2; }
    else if (sys == SYS_GLO) { freq[0] = sat2freq(sat, CODE_L1C, nav); freq[1] = sat2freq(sat, CODE_L2C, nav); }
    else if (sys == SYS_GAL) { freq[0] = FREQL1; freq[1] = FREQE5b; }
    else if (sys == SYS_CMP) { freq[0] = FREQ1_CMP; freq[1] = FREQ2_CMP; }
    else if (sys == SYS_IRN) { freq[0] = FREQL5; freq[1] = FREQs; }
    else return;

    C1 = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
    C2 = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

    /* 对双频无电离层组合进行天线相位偏移改正 */
    for (i = 0; i < 3; i++) {
        dant1 = pcv->off[0][0] * ex[i] + pcv->off[0][1] * ey[i] + pcv->off[0][2] * ez[i];
        dant2 = pcv->off[1][0] * ex[i] + pcv->off[1][1] * ey[i] + pcv->off[1][2] * ez[i];
        dant[i] = C1 * dant1 + C2 * dant2;
    }
}

/* =============================================================================
 * 功能：peph2pos (卫星位置与钟差计算主接口)
 * 说明：这是上层调用计算卫星状态的主入口，整合了精密星历轨道、精密时钟校正、
 * 天线相位中心偏移以及相对论效应改正。
 * ============================================================================= */
extern int peph2pos(gtime_t time, int sat, const nav_t* nav, int opt,
    double* rs, double* dts, double* var)
{
    gtime_t time_tt;
    double rss[3], rst[3], dtss[1], dtst[1], dant[3] = { 0 }, vare = 0.0, varc = 0.0, tt = 1E-3;
    int i;

    char tstr[40];
    trace(4, "peph2pos: time=%s sat=%2d opt=%d\n", time2str(time, tstr, 3), sat, opt);

    if (sat <= 0 || MAXSAT < sat) return 0;

    /* 1. 计算当前时刻的精密坐标 rss 和精密钟差 dtss */
    if (!pephpos(time, sat, nav, rss, dtss, &vare, &varc) ||
        !pephclk(time, sat, nav, dtss, &varc)) return 0;

    /* 2. 微分估算卫星速度 */
    time_tt = timeadd(time, tt);
    if (!pephpos(time_tt, sat, nav, rst, dtst, NULL, NULL) ||
        !pephclk(time_tt, sat, nav, dtst, NULL)) return 0;

    /* 3. 应用天线相位偏移 (APC) 改正 */
    if (opt) {
        satantoff(time, rss, sat, nav, dant);
    }
    for (i = 0; i < 3; i++) {
        rs[i] = rss[i] + dant[i];
        rs[i + 3] = (rst[i] - rss[i]) / tt;
    }
    /* 4. 应用狭义相对论效应修正钟差 (s) */
    if (dtss[0] != 0.0) {
        dts[0] = dtss[0] - 2.0 * dot3(rs, rs + 3) / CLIGHT / CLIGHT;
        dts[1] = (dtst[0] - dtss[0]) / tt;
    }
    else {
        dts[0] = dts[1] = 0.0;
    }
    if (var) *var = vare + varc;

    return 1;
}
