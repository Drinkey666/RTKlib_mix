/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.0  moved from rtkcmn.c
*                           changed api:
*                               pntpos()
*                           deleted api:
*                               pntvel()
*           2011/01/12 1.1  add option to include unhealthy satellite
*                           reject duplicated observation data
*                           changed api: ionocorr()
*           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
*           2012/12/25 1.3  add variable snr mask
*           2014/05/26 1.4  support galileo and beidou
*           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
*           2018/10/10 1.6  support api change of satexclude()
*           2020/11/30 1.7  support NavIC/IRNSS in pntpos()
*                           no support IONOOPT_LEX option in ioncorr()
*                           improve handling of TGD correction for each system
*                           use E1-E5b for Galileo dual-freq iono-correction
*                           use API sat2freq() to get carrier frequency
*                           add output of velocity estimation error in estvel()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* Smartphone PPP V3.2: UTF-8, MSVC /TC compatible, no variable-length arrays. */

/* constants/macros ----------------------------------------------------------*/

#define SQR(x)      ((x)*(x))
#define MAX(x,y)    ((x)>=(y)?(x):(y))

#define QZSDT /* enable GPS-QZS time offset estimation */
#ifdef QZSDT
#define NX          (4+5)       /* # of estimated parameters */
#else
#define NX          (4+4)       /* # of estimated parameters */
#endif
#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
# define MAX_GDOP   30          /* max gdop for valid solution  */

/* pseudorange measurement error variance 计算单颗卫星的伪距误差方差------------------------------------*/
static double varerr(const prcopt_t* opt, const obsd_t* obs, double el, int sys)
{
    double fact = 1.0, varr, snr_rover;

    // 1. 根据不同卫星系统赋予不同的基础误差放大因子
    switch (sys) {
    case SYS_GPS: fact *= EFACT_GPS; break;
    case SYS_GLO: fact *= EFACT_GLO; break;
    case SYS_GAL: fact *= EFACT_GAL; break;
    case SYS_SBS: fact *= EFACT_SBS; break;
    case SYS_CMP: fact *= EFACT_CMP; break;
    case SYS_QZS: fact *= EFACT_QZS; break;
    case SYS_IRN: fact *= EFACT_IRN; break;
    default:      fact *= EFACT_GPS; break;
    }

    // 2. 物理极限保护：限制最小高度角，防止除以零或无穷大
    if (el < MIN_EL) el = MIN_EL;

    // 3. 核心数学模型：计算基础方差
    // var = a^2 + b^2 / sin(el)  (a为常数误差，b为随高度角变化的误差)
    varr = SQR(opt->err[1]) + SQR(opt->err[2]) / sin(el);

    // 4. 信噪比(SNR)惩罚项 (极其适合手机数据！)
    if (opt->err[6] > 0.0) {  /* 如果配置了信噪比相关项 */
        snr_rover = obs->SNR[0] > 0 ? SNR_UNIT * obs->SNR[0] : opt->err[5];
        // 惩罚指数衰减：当实际信噪比低于理想最大值时，误差按指数级急剧放大
        varr += SQR(opt->err[6]) * pow(10, 0.1 * MAX(opt->err[5] - snr_rover, 0));
    }

    // 5. 伪距/相位权重比例尺乘数 (比如手机配置 opt->eratio[0] = 400，就会在这里把方差放大几万倍)
    varr *= SQR(opt->eratio[0]);

    if (opt->err[7] > 0.0) {
        varr += SQR(opt->err[7] * 0.01 * (1 << (obs->Pstd[0] + 5)));  /* 引入接收机自带的伪距标准差指标 */
    }

    // 6. 无电离层组合(IF)的毒性：IF组合会将原本的噪声放大 3 倍左右！
    if (opt->ionoopt == IONOOPT_IFLC) varr *= SQR(3.0); /* iono-free */

    return SQR(fact) * varr;
}
/* get group delay parameter (m) 提取卫星的硬件延迟 ---------------------------------------------*/
static double gettgd(int sat, const nav_t* nav, int type)
{
    int i, sys = satsys(sat, NULL);

    // 1. GLONASS 系统的频分多址(FDMA)处理逻辑较为特殊
    if (sys == SYS_GLO) {
        for (i = 0; i < nav->ng; i++) {
            if (nav->geph[i].sat == sat) break; // 在广播星历中找到这颗卫星
        }
        // 返回 dtaun 参数并乘以光速 CLIGHT，将时间延迟转换为距离误差(米)
        return (i >= nav->ng) ? 0.0 : -nav->geph[i].dtaun * CLIGHT;
    }
    else {
        // 2. 其他系统 (GPS, Galileo, BDS等)
        for (i = 0; i < nav->n; i++) {
            if (nav->eph[i].sat == sat) break; // 查找对应星历
        }
        // 返回星历中记录的 TGD 参数，乘以光速转化为米
        return (i >= nav->n) ? 0.0 : nav->eph[i].tgd[type] * CLIGHT;
    }
}
/* test SNR mask 信噪比门限测试-------------------------------------------------------------*/
static int snrmask(const obsd_t* obs, const double* azel, const prcopt_t* opt)
{
    int f2;

    // 1. 检查第一频点 (通常是 L1) 是否通过了高度角和信噪比的双重掩码测试
    if (testsnr(0, 0, azel[1], obs->SNR[0] * SNR_UNIT, &opt->snrmask)) {
        return 0; // 测试未通过，返回 0 (不可用)
    }

    // 2. 如果开启了双频无电离层组合 (IFLC)
    if (opt->ionoopt == IONOOPT_IFLC) {
        f2 = seliflc(opt->nf, satsys(obs->sat, NULL)); // 找到第二频点的索引 (比如 L2 或 L5)
        // 既然要做双频组合，第二频点的信号也必须及格，否则整个卫星弃用
        if (testsnr(0, f2, azel[1], obs->SNR[f2] * SNR_UNIT, &opt->snrmask)) return 0;
    }
    return 1; // 全部通过，返回 1 (可用)
}
/* iono-free or "pseudo iono-free" pseudorange with code bias correction 生成修正后的伪距 (-----*/
static double prange(const obsd_t* obs, const nav_t* nav, const prcopt_t* opt, double* var)
{
    double P1, P2, gamma, b1, b2;
    int sat, sys, f2, bias_ix, osb1 = 0, osb2 = 0;

    sat = obs->sat;
    sys = satsys(sat, NULL);
    P1 = obs->P[0]; // 获取第一频点原始伪距
    f2 = seliflc(opt->nf, satsys(obs->sat, NULL)); // 决定第二频点用什么(L2还是L5)
    P2 = obs->P[f2]; // 获取第二频点原始伪距
    *var = 0.0;

    // 如果观测值缺失，直接返回 0
    if (P1 == 0.0 || (opt->ionoopt == IONOOPT_IFLC && P2 == 0.0)) return 0.0;

    // ==============================================
    // 步骤一：卫星端硬件延迟/码偏差 (Code Bias) 修正
    // ==============================================
    /* Bias-SINEX absolute OSB has priority over legacy relative DCB.
       Sign convention: corrected code = observed code - OSB. */
    if (obs->code[0] > CODE_NONE && obs->code[0] <= MAXCODE &&
        nav->osb_valid[sat - 1][obs->code[0]]) {
        P1 -= nav->osb[sat - 1][obs->code[0]];
        osb1 = 1;
    }
    else {
        bias_ix = code2bias_ix(sys, obs->code[0]);
        if (bias_ix > 0 && bias_ix <= MAX_CODE_BIASES)
            P1 += nav->cbias[sat - 1][0][bias_ix - 1];
    }

    if (f2 >= 0 && f2 < NFREQ && obs->code[f2] > CODE_NONE &&
        obs->code[f2] <= MAXCODE && nav->osb_valid[sat - 1][obs->code[f2]]) {
        P2 -= nav->osb[sat - 1][obs->code[f2]];
        osb2 = 1;
    }
    else if (f2 >= 0 && f2 < MAX_CODE_BIAS_FREQS) {
        bias_ix = code2bias_ix(sys, obs->code[f2]);
        if (bias_ix > 0 && bias_ix <= MAX_CODE_BIASES)
            P2 += nav->cbias[sat - 1][f2][bias_ix - 1];
    }

    // ==============================================
    // 步骤二：构建无电离层组合 (Dual-frequency Iono-Free)
    // ==============================================
    if (opt->ionoopt == IONOOPT_IFLC) {

        /* With two absolute OSBs, form IF directly from the corrected native
           codes. Broadcast TGD/BGD must not be applied a second time. */
        if (osb1 && osb2) {
            double f1hz = sat2freq(sat, obs->code[0], nav);
            double f2hz = sat2freq(sat, obs->code[f2], nav);
            if (f1hz > 0.0 && f2hz > 0.0 && fabs(f1hz - f2hz) > 1.0) {
                gamma = SQR(f1hz / f2hz);
                return (P2 - gamma * P1) / (1.0 - gamma);
            }
        }

        if (sys == SYS_GPS || sys == SYS_QZS) { /* GPS 和 QZSS 系统 */
            // gamma 是频率平方比值：(f1 / f2)^2
            gamma = f2 == 1 ? SQR(FREQL1 / FREQL2) : SQR(FREQL1 / FREQL5);
            // 核心公式：IF组合伪距 = (P2 - gamma*P1) / (1 - gamma)
            // 这样组合出来的伪距，其电离层延迟项会被完全抵消
            return (P2 - gamma * P1) / (1.0 - gamma);
        }
        else if (sys == SYS_GLO) { /* GLONASS 系统 (频率稍有不同) */
            gamma = f2 == 1 ? SQR(FREQ1_GLO / FREQ2_GLO) : SQR(FREQ1_GLO / FREQ3_GLO);
            return (P2 - gamma * P1) / (1.0 - gamma);
        }
        else if (sys == SYS_GAL) { /* Galileo 系统 */
            gamma = f2 == 1 ? SQR(FREQL1 / FREQE5b) : SQR(FREQL1 / FREQL5);
            if (f2 == 1 && getseleph(SYS_GAL)) {
                P2 -= gettgd(sat, nav, 0) - gettgd(sat, nav, 1); // 伽利略特有的BGD组延迟修正
            }
            return (P2 - gamma * P1) / (1.0 - gamma);
        }
        else if (sys == SYS_CMP) { /* 北斗系统 (BDS) */
            // 北斗的频率和 TGD/ISC 修正非常复杂，需要严格区分信号频段
            gamma = SQR(((obs->code[0] == CODE_L2I) ? FREQ1_CMP : FREQL1) / FREQ2_CMP);
            if (obs->code[0] == CODE_L2I) b1 = gettgd(sat, nav, 0);
            else if (obs->code[0] == CODE_L1P) b1 = gettgd(sat, nav, 2);
            else b1 = gettgd(sat, nav, 2) + gettgd(sat, nav, 4);
            b2 = gettgd(sat, nav, 1);
            // 组合公式同样消除了电离层，并补偿了北斗特有的 TGD/ISC 延迟
            return ((P2 - gamma * P1) - (b2 - gamma * b1)) / (1.0 - gamma);
        }
        else if (sys == SYS_IRN) { /* NavIC 系统 */
            gamma = SQR(FREQL5 / FREQs);
            return (P2 - gamma * P1) / (1.0 - gamma);
        }
    }
    // ==============================================
    // 步骤三：单频解算 (Single-frequency) 处理
    // ==============================================
    else {
        *var = SQR(ERR_CBIAS);

        /* Absolute OSB is clock-datum consistent; do not apply broadcast
           TGD/BGD again to that already corrected code observable. */
        if (osb1) return P1; // 赋予一个默认的码偏差方差

        // 单频无法抵消电离层，只能减去星历中播发的 TGD (时间群延迟)
        if (sys == SYS_GPS || sys == SYS_QZS) {
            b1 = gettgd(sat, nav, 0);
            return P1 - b1; // 单频伪距 - TGD修正
        }
        else if (sys == SYS_GLO) {
            gamma = SQR(FREQ1_GLO / FREQ2_GLO);
            b1 = gettgd(sat, nav, 0);
            return P1 - b1 / (gamma - 1.0);
        }
        else if (sys == SYS_GAL) {
            if (getseleph(SYS_GAL)) b1 = gettgd(sat, nav, 0);
            else                    b1 = gettgd(sat, nav, 1);
            return P1 - b1;
        }
        else if (sys == SYS_CMP) {
            if (obs->code[0] == CODE_L2I) b1 = gettgd(sat, nav, 0);
            else if (obs->code[0] == CODE_L1P) b1 = gettgd(sat, nav, 2);
            else b1 = gettgd(sat, nav, 2) + gettgd(sat, nav, 4);
            return P1 - b1;
        }
        else if (sys == SYS_IRN) {
            gamma = SQR(FREQs / FREQL5);
            b1 = gettgd(sat, nav, 0);
            return P1 - gamma * b1;
        }
    }
    return P1; // 最终兜底防错，返回原始伪距
}
/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t* nav, int sat, const double* pos,
    const double* azel, int ionoopt, double* ion, double* var)
{
    int err = 0;
    char tstr[40];

    // 打印调试信息
    trace(4, "ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
        time2str(time, tstr, 3), ionoopt, sat, pos[0] * R2D, pos[1] * R2D, azel[0] * R2D,
        azel[1] * R2D);

    /* 1. SBAS (星基增强系统) 电离层模型 */
    if (ionoopt == IONOOPT_SBAS) {
        if (sbsioncorr(time, nav, pos, azel, ion, var)) return 1;
        err = 1; // 如果 SBAS 数据不可用，标记错误并降级处理
    }

    /* 2. IONEX TEC 格网模型 (PPP 最常用) */
    // 如果你加载了 .26i 这类文件，就会走这里，精度很高
    if (ionoopt == IONOOPT_TEC) {
        if (iontec(time, nav, pos, azel, 1, ion, var)) return 1;
        err = 1;
    }

    /* 3. QZSS 广播电离层模型 (日本准天顶系统特有) */
    if (ionoopt == IONOOPT_QZS && norm(nav->ion_qzs, 8) > 0.0) {
        *ion = ionmodel(time, nav->ion_qzs, pos, azel);
        *var = SQR(*ion * ERR_BRDCI);
        return 1;
    }

    /* 4. GPS 广播电离层模型 (Klobuchar 模型) */
    // 这是单点定位最常见的保底方案。如果前面高级模型失败 (err==1)，也会退化到这里
    if (ionoopt == IONOOPT_BRDC || err == 1) {
        *ion = ionmodel(time, nav->ion_gps, pos, azel);
        *var = SQR(*ion * ERR_BRDCI); // 广播星历电离层方差通常较大
        return 1;
    }

    /* 5. 关闭电离层改正 (或者使用双频无电离层组合时，直接设为0) */
    *ion = 0.0;
    *var = ionoopt == IONOOPT_OFF ? SQR(ERR_ION) : 0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t* nav, const double* pos,
    const double* azel, int tropopt, double* trp, double* var)
{
    char tstr[40];
    trace(4, "tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
        time2str(time, tstr, 3), tropopt, pos[0] * R2D, pos[1] * R2D, azel[0] * R2D,
        azel[1] * R2D);

    /* 1. Saastamoinen 经验模型 (GNSS 界的绝对主流) */
    // 无论是单点定位(SAAS)还是PPP估计(EST/ESTG)，初始保底延迟都用此模型计算
    if (tropopt == TROPOPT_SAAS || tropopt == TROPOPT_EST || tropopt == TROPOPT_ESTG) {
        *trp = tropmodel(time, pos, azel, REL_HUMI); // 使用标准大气参数计算
        *var = SQR(ERR_SAAS / (sin(azel[1]) + 0.1));  // 仰角越低，通过的大气越厚，方差越大
        return 1;
    }

    /* 2. SBAS 对流层模型 */
    if (tropopt == TROPOPT_SBAS) {
        *trp = sbstropcorr(time, pos, azel, var);
        return 1;
    }

    /* 3. 关闭对流层改正 */
    *trp = 0.0;
    *var = tropopt == TROPOPT_OFF ? SQR(ERR_TROP) : 0.0;
    return 1;
}
/* pseudorange residuals -----------------------------------------------------*/
/* pseudorange residuals -----------------------------------------------------*/
static int rescode(int iter, const obsd_t* obs, int n, const double* rs,
    const double* dts, const double* vare, const int* svh,
    const nav_t* nav, const double* x, const prcopt_t* opt,
    const ssat_t* ssat, double* v, double* H, double* var,
    double* azel, int* vsat, double* resp, int* ns)
{
    gtime_t time;
    double r, freq, dion = 0.0, dtrp = 0.0, vmeas, vion = 0.0, vtrp = 0.0, rr[3], pos[3], dtr, e[3], P;
    int i, j, nv = 0, sat, sys, mask[NX - 3] = { 0 };

    // 1. 提取当前估计的接收机位置 (x, y, z) 和接收机钟差 dtr
    for (i = 0; i < 3; i++) rr[i] = x[i];
    dtr = x[3]; // x[3] 通常存放 GPS 系统的接收机钟差

    ecef2pos(rr, pos); // 将空间直角坐标 (XYZ) 转换为大地坐标 (经纬高)，用于算大气
    trace(3, "rescode: rr=%.3f %.3f %.3f\n", rr[0], rr[1], rr[2]);

    // 2. 开始遍历当前历元接收到的所有卫星观测值
    for (i = *ns = 0; i < n && i < MAXOBS; i++) {
        vsat[i] = 0; azel[i * 2] = azel[1 + i * 2] = resp[i] = 0.0;
        time = obs[i].time;
        sat = obs[i].sat;
        if (!(sys = satsys(sat, NULL))) continue; // 识别卫星所属系统 (GPS, BDS 等)

        /* 剔除重复的观测数据 */
        if (i < n - 1 && i < MAXOBS - 1 && sat == obs[i + 1].sat) {
            i++; continue;
        }

        /* 剔除被用户手动排除或健康状态(svh)异常的卫星 */
        if (satexclude(sat, vare[i], svh[i], opt)) continue;

        /* 3. 计算接收机到卫星的真实几何距离 r，并获取视线向量 e */
        if ((r = geodist(rs + i * 6, rr, e)) <= 0.0) continue;
        /* 计算并测试高度角掩码 (低于设置的 elmin 则剔除) */
        if (satazel(pos, e, azel + i * 2) < opt->elmin) continue;

        // 如果是第二次及以后的迭代 (iter>0)，才开始计算复杂大气模型
        if (iter > 0) {
            /* 测试信噪比掩码 (上一轮讲过的 SNR 过滤) */
            if (!snrmask(obs + i, azel + i * 2, opt)) continue;

            /* 计算电离层延迟 (调用刚刚的 ionocorr 函数) */
            if (!ionocorr(time, nav, sat, pos, azel + i * 2, opt->ionoopt, &dion, &vion)) continue;

            /* 频率转换：底层模型通常算的是 L1 频点延迟，需根据当前频点换算 */
            if ((freq = sat2freq(sat, obs[i].code[0], nav)) == 0.0) continue;
            dion *= SQR(FREQL1 / freq);
            vion *= SQR(SQR(FREQL1 / freq));

            /* 计算对流层延迟 (调用刚刚的 tropcorr 函数) */
            if (!tropcorr(time, nav, pos, azel + i * 2, opt->tropopt, &dtrp, &vtrp)) continue;
        }

        /* 4. 提取伪距 P，并应用硬件延迟偏差(Code Bias)修正 (调用你之前问过的 prange) */
        if ((P = prange(obs + i, nav, opt, &vmeas)) == 0.0) continue;

        /* 5. 核心：计算伪距残差 v (OMC: 观测值 - 计算值) */
        // P: 测量伪距
        // r: 几何距离, dtr: 接收机钟差, dts: 卫星钟差(已乘光速), dion: 电离层, dtrp: 对流层
        v[nv] = P - (r + dtr - CLIGHT * dts[i * 2] + dion + dtrp);

        /* 6. 构建设计矩阵 H (偏导数矩阵，用于最小二乘或滤波的增益更新) */
        // 前三个元素是位置的偏导数 (视线向量的反方向)
        for (j = 0; j < NX; j++) {
            H[j + nv * NX] = j < 3 ? -e[j] : (j == 3 ? 1.0 : 0.0); // 第4个元素对应GPS钟差，偏导为1
        }

        /* 7. 处理多系统间的系统时间偏差 (ISB, Inter-System Bias) */
        // 由于不同卫星系统的时间基准不同，接收机需要为 GLONASS, Galileo, BeiDou 额外估计一个相对于 GPS 的钟差
        if (sys == SYS_GLO) { v[nv] -= x[4]; H[4 + nv * NX] = 1.0; mask[1] = 1; }
        else if (sys == SYS_GAL) { v[nv] -= x[5]; H[5 + nv * NX] = 1.0; mask[2] = 1; }
        else if (sys == SYS_CMP) { v[nv] -= x[6]; H[6 + nv * NX] = 1.0; mask[3] = 1; }
        else if (sys == SYS_IRN) { v[nv] -= x[7]; H[7 + nv * NX] = 1.0; mask[4] = 1; }
#ifdef QZSDT
        else if (sys == SYS_QZS) { v[nv] -= x[8]; H[8 + nv * NX] = 1.0; mask[5] = 1; }
#endif
        else mask[0] = 1; // GPS 系统

        vsat[i] = 1; resp[i] = v[nv]; (*ns)++; // 标记有效卫星，记录残差

        /* 8. 计算该观测值的综合方差 (观测方差 + 轨道/钟方差 + 电离层方差 + 对流层方差) */
        var[nv] = vare[i] + vmeas + vion + vtrp;
        // 加上最底层的伪距测量噪声方差 (调用你之前问过的 varerr 函数)
        var[nv++] += varerr(opt, &obs[i], azel[1 + i * 2], sys);
    }

    /* 9. 约束处理以防止秩亏 (Rank-Deficient) */
    // 比如当前历元只有 GPS 和 BDS 卫星，没有 GLONASS 卫星。
    // 那么 GLONASS 的系统时差参数就无法估计。为了防止矩阵无法求逆导致崩溃，赋予它一个虚拟的极小方差和零残差。
    for (i = 0; i < NX - 3; i++) {
        if (mask[i]) continue; // 如果该系统有卫星，正常处理
        v[nv] = 0.0;             // 虚拟残差
        for (j = 0; j < NX; j++) H[j + nv * NX] = j == i + 3 ? 1.0 : 0.0;
        var[nv++] = 0.01;        // 给定一个小方差约束
    }
    return nv; // 返回有效观测方程的总个数
}
/* validate solution 解算结果质量检验---------------------------------------------------------*/
static int valsol(const double* azel, const int* vsat, int n,
    const prcopt_t* opt, const double* v, int nv, int nx, char* msg)
{
    double azels[MAXOBS * 2], dop[4], vv;
    int i, ns;

    trace(3, "valsol  : n=%d nv=%d\n", n, nv);

    /* 1. 残差的卡方 (Chi-square) 检验 */
    vv = dot(v, v, nv); // 计算残差向量 v 的平方和 (v^T * v)
    // nv是方程数(观测数)，nx是未知数。如果平方和大于统计学阈值 chisqr
    if (nv > nx && vv > chisqr[nv - nx - 1]) {
        sprintf(msg, "Warning: large chi-square error nv=%d vv=%.1f cs=%.1f", nv, vv, chisqr[nv - nx - 1]);
        /* 手机数据残差极大，RTKLIB 源码在这里把 return 0 注释掉了，只是报警但继续使用 */
    }

    /* 2. GDOP (几何精度因子) 检验 */
    for (i = ns = 0; i < n; i++) {
        if (!vsat[i]) continue; // 只提取参与解算的有效卫星
        azels[ns * 2] = azel[i * 2];   // 方位角
        azels[1 + ns * 2] = azel[1 + i * 2]; // 高度角
        ns++;
    }
    dops(ns, azels, opt->elmin, dop); // 计算 DOP 值矩阵

    // 如果几何分布极差（比如卫星排成一条线），拒接接受该结果
    if (dop[0] <= 0.0 || dop[0] > MAX_GDOP) {
        sprintf(msg, "gdop error nv=%d gdop=%.1f", nv, dop[0]);
        return 0; // 检验不通过
    }
    return 1; // 质检合格
}
/* estimate receiver position 接收机位置估计------------------------------------------------*/
static int estpos(const obsd_t* obs, int n, const double* rs, const double* dts,
    const double* vare, const int* svh, const nav_t* nav,
    const prcopt_t* opt, const ssat_t* ssat, sol_t* sol, double* azel,
    int* vsat, double* resp, char* msg)
{
    double x[NX] = { 0 }, dx[NX], Q[NX * NX], * v, * H, * var, sig;
    int i, j, k, info, stat, nv, ns;

    v = mat(n + NX - 3, 1); H = mat(NX, n + NX - 3); var = mat(n + NX - 3, 1);

    for (i = 0; i < 3; i++) x[i] = sol->rr[i]; // 取上一历元的位置作为迭代初值

    // 开始最大 MAXITR (通常为 10) 次的最小二乘迭代
    for (i = 0; i < MAXITR; i++) {

        /* 1. 构建观测方程，获取伪距残差 v, 设计矩阵 H, 方差 var (调用之前的 rescode) */
        nv = rescode(i, obs, n, rs, dts, vare, svh, nav, x, opt, ssat, v, H, var, azel, vsat, resp, &ns);

        // 自由度检查：有效观测数必须大于未知数（通常是3个坐标+系统时差）
        if (nv < NX) {
            sprintf(msg, "lack of valid sats ns=%d", nv);
            break;
        }

        /* 2. 赋予权重：用方差的平方根去除 残差 v 和 矩阵 H */
        for (j = 0; j < nv; j++) {
            sig = sqrt(var[j]); // 标准差
            v[j] /= sig;        // 误差大的卫星，这里除以一个大数，残差贡献变小 (降权)
            for (k = 0; k < NX; k++) H[k + j * NX] /= sig;
        }

        /* 3. 核心数学运算：最小二乘求解 (H^T * H * dx = H^T * v) */
        // dx 是本次迭代计算出的改正数，Q 是协方差矩阵 (用于评估精度)
        if ((info = lsq(H, v, NX, nv, dx, Q))) {
            sprintf(msg, "lsq error info=%d", info);
            break; // 矩阵奇异，无法求逆
        }

        // 更新未知数：当前状态 + 改正数
        for (j = 0; j < NX; j++) x[j] += dx[j];

        /* 4. 收敛判定 */
        // 如果位置改正数 dx 的模长小于 1E-4 米 (0.1毫米)，认为已经找到最优解！
        if (norm(dx, NX) < 1E-4) {
            sol->type = 0; // 0: XYZ坐标系
            sol->time = timeadd(obs[0].time, -x[3] / CLIGHT); // 修正接收机钟差引起的时间延迟

            // 保存接收机钟差 (GPS) 及其他系统的系统间时差 (ISB)
            sol->dtr[0] = x[3] / CLIGHT; /* GPS clock bias */
            sol->dtr[1] = x[4] / CLIGHT; /* GLO-GPS time offset */
            sol->dtr[2] = x[5] / CLIGHT; /* GAL-GPS time offset */
            sol->dtr[3] = x[6] / CLIGHT; /* BDS-GPS time offset */
            sol->dtr[4] = x[7] / CLIGHT; /* IRN-GPS time offset */
#ifdef QZSDT
            sol->dtr[5] = x[8] / CLIGHT; /* QZS-GPS time offset */
#else
            sol->dtr[5] = 0.0;
#endif

            for (j = 0; j < 6; j++) sol->rr[j] = j < 3 ? x[j] : 0.0; // 保存 XYZ 位置

            /* 保存协方差 Q 阵到解结构体中，用于后续 PPP 滤波的初始方差！ */
            sol->qr[0] = (float)Q[0];
            sol->qr[1] = (float)Q[1 + NX];
            sol->qr[2] = (float)Q[2 + 2 * NX];
            sol->qr[3] = (float)Q[1];
            sol->qr[4] = (float)Q[2 + NX];
            sol->qr[5] = (float)Q[2];
            sol->ns = (uint8_t)ns;

            /* 5. 调用 valsol 进行质量检验 */
            if ((stat = valsol(azel, vsat, n, opt, v, nv, NX, msg))) {
                sol->stat = opt->sateph == EPHOPT_SBAS ? SOLQ_SBAS : SOLQ_SINGLE; // 标记为单点定位解
            }
            free(v); free(H); free(var);
            return stat; // 成功返回
        }
    }
    if (i >= MAXITR) sprintf(msg, "iteration divergent i=%d", i); // 迭代发散，未能收敛

    free(v); free(H); free(var);
    return 0;
}
/* Startup-only post-fit check. The full clock/ISB design is retained; common
 * clock jumps are therefore not mistaken for individual code outliers.
 * Require two residual degrees of freedom and a 30 m AND 6-sigma fault.
 * Leverage correction prevents a fitted high-leverage error hiding itself.
 * Ordinary pntpos()/single/relative positioning are not changed. */
static int startup_metrics(const obsd_t *obs,int n,const double *rs,
    const double *dts,const double *vare,const int *svh,const nav_t *nav,
    const prcopt_t *opt,const ssat_t *ssat,const sol_t *sol,
    double *score,double *worst_z,double *worst_m,int *redundancy)
{
    double x[NX]={0},v[MAXOBS+NX],H[NX*(MAXOBS+NX)],var[MAXOBS+NX];
    double azel[MAXOBS*2],resp[MAXOBS],dx[NX],Q[NX*NX],sum=0.0;
    int used[MAXOBS]={0},i,j,k,nv,ns,row=0,candidate=-1;
    for(i=0;i<3;i++)x[i]=sol->rr[i];
    for(i=3;i<NX;i++)x[i]=sol->dtr[i-3]*CLIGHT;
    nv=rescode(1,obs,n,rs,dts,vare,svh,nav,x,opt,ssat,v,H,var,azel,used,resp,&ns);
    *score=1E99;*worst_z=*worst_m=0.0;*redundancy=nv-NX;
    if(*redundancy<2)return -1;
    for(i=0;i<nv;i++) {
        double sig=sqrt(var[i]);
        if(!isfinite(sig)||sig<=0.0||!isfinite(v[i]))return -1;
        v[i]/=sig;sum+=v[i]*v[i];
        for(j=0;j<NX;j++)H[j+i*NX]/=sig;
    }
    if(lsq(H,v,NX,nv,dx,Q))return -1;
    *score=sum/(*redundancy);
    for(i=0;i<n;i++)if(used[i]) {
        double leverage=0.0,z;
        for(j=0;j<NX;j++)for(k=0;k<NX;k++)
            leverage+=H[j+row*NX]*Q[j+k*NX]*H[k+row*NX];
        z=fabs(v[row])/sqrt(MAX(0.05,1.0-leverage));
        if(fabs(resp[i])>30.0&&z>6.0&&z>*worst_z) {
            candidate=i;*worst_z=z;*worst_m=resp[i];
        }
        row++;
    }
    return candidate;
}
static int startup_fde(const obsd_t *obs,int n,const double *rs,
    const double *dts,const double *vare,const int *svh,const nav_t *nav,
    const prcopt_t *opt,const ssat_t *ssat,sol_t *sol,double *azel,
    int *vsat,double *resp,int initial_ok,unsigned char *rejected,char *msg)
{
    obsd_t trial_obs[MAXOBS],work[MAXOBS];
    sol_t trial,best;
    double trial_azel[MAXOBS*2],trial_resp[MAXOBS],best_azel[MAXOBS*2],best_resp[MAXOBS];
    double score,z,metres,best_score;int i,round,dof,candidate,best_i,ok=initial_ok;
    int trial_used[MAXOBS],best_used[MAXOBS];char trial_msg[128],t[40],sid[8];
    memcpy(work,obs,n*sizeof(*obs));time2str(obs[0].time,t,3);
    for(round=0;round<3;round++) {
        candidate=ok?startup_metrics(work,n,rs,dts,vare,svh,nav,opt,ssat,sol,&score,&z,&metres,&dof):-1;
        if(ok&&candidate<0)return 1; /* sparse geometry: do not invent evidence */
        best_i=-1;best_score=ok?score*0.5:9.0;
        /* A trial must actually improve a redundant fit. Evaluate all used
         * hypotheses instead of blindly deleting the largest raw residual. */
        for(i=0;i<n;i++) {
            double trial_score,tz,tm;int tdof;
            if(work[i].P[0]==0.0||(ok&&!vsat[i]))continue;
            memcpy(trial_obs,work,n*sizeof(*obs));trial_obs[i].P[0]=0.0;
            trial=*sol;trial.stat=SOLQ_NONE;
            if(!estpos(trial_obs,n,rs,dts,vare,svh,nav,opt,ssat,&trial,trial_azel,trial_used,trial_resp,trial_msg))continue;
            startup_metrics(trial_obs,n,rs,dts,vare,svh,nav,opt,ssat,&trial,&trial_score,&tz,&tm,&tdof);
            if(tdof<2||!isfinite(trial_score)||trial_score>=best_score)continue;
            best_score=trial_score;best_i=i;best=trial;
            memcpy(best_azel,trial_azel,n*2*sizeof(double));
            memcpy(best_resp,trial_resp,n*sizeof(double));memcpy(best_used,trial_used,n*sizeof(int));
        }
        if(best_i<0) {
            trace(2,"$START_SPP_FAIL,%s,reason=NO_SAFE_FDE,initial_ok=%d\n",t,ok);
            strcpy(msg,"startup gross code: no safe exclusion");sol->stat=SOLQ_NONE;return 0;
        }
        satno2id(work[best_i].sat,sid);
        trace(2,"$START_SPP_REJECT,%s,sat=%s,sig=%s,res=%.3f,worst_z=%.2f,score_before=%.3f,score_after=%.3f\n",
            t,sid,code2obs(work[best_i].code[0]),ok?resp[best_i]:0.0,ok?z:0.0,ok?score:0.0,best_score);
        work[best_i].P[0]=0.0;rejected[best_i]=1;*sol=best;ok=1;
        memcpy(azel,best_azel,n*2*sizeof(double));memcpy(resp,best_resp,n*sizeof(double));memcpy(vsat,best_used,n*sizeof(int));
    }
    candidate=startup_metrics(work,n,rs,dts,vare,svh,nav,opt,ssat,sol,&score,&z,&metres,&dof);
    if(candidate>=0){strcpy(msg,"startup gross code: exclusion budget exhausted");sol->stat=SOLQ_NONE;return 0;}
    return ok;
}
/* RAIM FDE (failure detection and exclusion) 接收机自主完好性监测-------------------------------*/
static int raim_fde(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, const ssat_t *ssat, 
                    sol_t *sol, double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[40],name[8],msg_e[128];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;
    
    trace(3,"raim_fde: %s n=%2d\n",time2str(obs[0].time,tstr,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n); 
    
    for (i=0;i<n;i++) {
        
        /* satellite exclusion */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e +6*k,rs +6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,ssat,&sol_e,azel_e,
                    vsat_e,resp_e,msg_e)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                  obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);
        
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        sol_e.eventime = sol->eventime;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
        strcpy(msg,msg_e);
    }
#ifdef TRACE
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }
#endif
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}
/* range rate residuals ------------------------------------------------------*/
static int resdop(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double err, double *v,
                  double *H)
{
    double freq,rate,pos[3],E[9],a[3],e[3],vs[3],cosel,sig;
    int i,j,nv=0;
    
    trace(3,"resdop  : n=%d\n",n);
    
    ecef2pos(rr,pos); xyz2enu(pos,E);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        freq=sat2freq(obs[i].sat,obs[i].code[0],nav);
        
        if (obs[i].D[0]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
            continue;
        }
        /* LOS (line-of-sight) vector in ECEF */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,1,3,E,a,e);
        
        /* satellite velocity relative to receiver in ECEF */
        for (j=0;j<3;j++) {
            vs[j]=rs[j+3+i*6]-x[j];
        }
        /* range rate with earth rotation correction */
        rate=dot3(vs,e)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                     rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        
        /* Std of range rate error (m/s) */
        sig=(err<=0.0)?1.0:err*CLIGHT/freq;
        
        /* range rate residual (m/s) */
        v[nv]=(-obs[i].D[0]*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]))/sig;
        
        /* design matrix */
        for (j=0;j<4;j++) {
            H[j+nv*4]=((j<3)?-e[j]:1.0)/sig;
        }
        nv++;
    }
    return nv;
}
/* estimate receiver velocity ------------------------------------------------*/
static void estvel(const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat)
{
    double x[4]={0},dx[4],Q[16],*v,*H;
    double err=opt->err[4]; /* Doppler error (Hz) */
    int i,j,nv;
    
    v=mat(n,1); H=mat(4,n);
    
    for (i=0;i<MAXITR;i++) {
        
        /* range rate residuals (m/s) */
        if ((nv=resdop(obs,n,rs,dts,nav,sol->rr,x,azel,vsat,err,v,H))<4) {
            break;
        }
        /* least square estimation */
        if (lsq(H,v,4,nv,dx,Q)) break;
        
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        if (norm(dx,4)<1E-6) {
            trace(3,"estvel : vx=%.3f vy=%.3f vz=%.3f, n=%d\n",x[0],x[1],x[2],n);
            matcpy(sol->rr+3,x,3,1);
            sol->qv[0]=(float)Q[0];  /* xx */
            sol->qv[1]=(float)Q[5];  /* yy */
            sol->qv[2]=(float)Q[10]; /* zz */
            sol->qv[3]=(float)Q[1];  /* xy */
            sol->qv[4]=(float)Q[6];  /* yz */
            sol->qv[5]=(float)Q[2];  /* zx */
            break;
        }
    }
    free(v); free(H);
}
/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error) 总入口
*-----------------------------------------------------------------------------*/
static int pntpos_impl(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg,int startup)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
    int i,stat,vsat[MAXOBS]={0},svh[MAXOBS];
    unsigned char rejected[MAXOBS]={0};
    
    char tstr[40];
    if(!obs||n<=0||n>MAXOBS){sol->stat=SOLQ_NONE;strcpy(msg,"invalid observation count");return 0;}
    trace(3,"pntpos  : tobs=%s n=%d\n",time2str(obs[0].time,tstr,3),n);
    
    sol->stat=SOLQ_NONE;
    
    if (n<=0) {
        strcpy(msg,"no observation data");
        return 0;
    }
    sol->time=obs[0].time;
    msg[0]='\0';
    sol->eventime = obs[0].eventime;
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
    
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].snr_rover[0]=0;
            ssat[i].snr_base[0]=0;
        }
        for (i=0;i<n;i++)
            ssat[obs[i].sat-1].snr_rover[0]=obs[i].SNR[0];
    }
    // 如果不是纯 SPP 模式 (比如正在做 PPP 初始化)，强制使用广播星历电离层和 Saastamoinen 对流层作为底座
    if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }

    /* 步骤 2：核心动作，利用伪距估计接收机位置和钟差 */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* estimate receiver position and time with pseudorange */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    if(startup)stat=startup_fde(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,stat,rejected,msg);
    
    /* 步骤 3：如果开启了 RAIM 且卫星数足够，执行完好性检测剔除坏星 */
    if (!startup&&!stat&&n>=6&&opt->posopt[4]) {
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    }
    /* 步骤 4：位置算成功了，利用多普勒观测值估计接收机速度 */
    if (stat) {
        estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
        }
        for (i=0;i<n;i++) {
            if(startup)ssat[obs[i].sat-1].ppp_code_bad[0]=rejected[i];
            ssat[obs[i].sat-1].azel[0]=azel_[  i*2];
            ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1;
            ssat[obs[i].sat-1].resp[0]=resp[i];
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}
extern int pntpos(const obsd_t *obs,int n,const nav_t *nav,const prcopt_t *opt,
    sol_t *sol,double *azel,ssat_t *ssat,char *msg)
{ return pntpos_impl(obs,n,nav,opt,sol,azel,ssat,msg,0); }
extern int pntpos_startup(const obsd_t *obs,int n,const nav_t *nav,const prcopt_t *opt,
    sol_t *sol,double *azel,ssat_t *ssat,char *msg)
{ return pntpos_impl(obs,n,nav,opt,sol,azel,ssat,msg,1); }
