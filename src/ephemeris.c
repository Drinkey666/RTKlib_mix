/*------------------------------------------------------------------------------
* ephemeris.c : satellite ephemeris and clock functions
*
*          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
*
* references :
*     [1] IS-GPS-200K, Navstar GPS Space Segment/Navigation User Interfaces,
*         May 6, 2019
*     [2] Global Navigation Satellite System GLONASS, Interface Control Document
*         Navigational radiosignal In bands L1, L2, (Version 5.1), 2008
*     [3] RTCA/DO-229C, Minimum operational performance standards for global
*         positioning system/wide area augmentation system airborne equipment,
*         RTCA inc, November 28, 2001
*     [4] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*         Code Biases, URA
*     [5] RTCM Paper 012-2009-SC104-528, January 28, 2009 (previous ver of [4])
*     [6] RTCM Paper 012-2009-SC104-582, February 2, 2010 (previous ver of [4])
*     [7] European GNSS (Galileo) Open Service Signal In Space Interface Control
*         Document, Issue 1.3, December, 2016
*     [8] Quasi-Zenith Satellite System Interface Specification Satellite
*         Positioning, Navigation and Timing Service (IS-QZSS-PNT-003), Cabinet
*         Office, November 5, 2018
*     [9] BeiDou navigation satellite system signal in space interface control
*         document open service signal B1I (version 3.0), China Satellite
*         Navigation office, February, 2019
*     [10] RTCM Standard 10403.3, Differential GNSS (Global Navigation
*         Satellite Systems) Services - version 3, October 7, 2016
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.1  moved from rtkcmn.c
*                           added api:
*                               eph2clk(),geph2clk(),seph2clk(),satantoff()
*                               satposs()
*                           changed api:
*                               eph2pos(),geph2pos(),satpos()
*                           deleted api:
*                               satposv(),satposiode()
*           2010/08/26 1.2  add ephemeris option EPHOPT_LEX
*           2010/09/09 1.3  fix problem when precise clock outage
*           2011/01/12 1.4  add api alm2pos()
*                           change api satpos(),satposs()
*                           enable valid unhealthy satellites and output status
*                           fix bug on exception by glonass ephem computation
*           2013/01/10 1.5  support beidou (compass)
*                           use newton's method to solve kepler eq.
*                           update ssr correction algorithm
*           2013/03/20 1.6  fix problem on ssr clock relativistic correction
*           2013/09/01 1.7  support negative pseudorange
*                           fix bug on variance in case of ura ssr = 63
*           2013/11/11 1.8  change constant MAXAGESSR 70.0 -> 90.0
*           2014/10/24 1.9  fix bug on return of var_uraeph() if ura<0||15<ura
*           2014/12/07 1.10 modify MAXDTOE for qzss,gal and bds
*                           test max number of iteration for Kepler
*           2015/08/26 1.11 update RTOL_ELPLER 1E-14 -> 1E-13
*                           set MAX_ITER_KEPLER for alm2pos()
*           2017/04/11 1.12 fix bug on max number of obs data in satposs()
*           2018/10/10 1.13 update reference [7]
*                           support ura value in var_uraeph() for galileo
*                           test eph->flag to recognize beidou geo
*                           add api satseleph() for ephemeris selection
*           2020/11/30 1.14 update references [1],[2],[8],[9] and [10]
*                           add API getseleph()
*                           rename API satseleph() as setseleph()
*                           support NavIC/IRNSS by API satpos() and satposs()
*                           support BDS C59-63 as GEO satellites in eph2pos()
*                           default selection of I/NAV for Galileo ephemeris
*                           no support EPHOPT_LEX by API satpos() and satposs()
*                           unselect Galileo ephemeris with AOD<=0 in seleph()
*                           fix bug on clock iteration in eph2clk(), geph2clk()
*                           fix bug on clock reference time in satpos_ssr()
*                           fix bug on wrong value with ura=15 in var_ura()
*                           use integer types in stdint.h
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants and macros ------------------------------------------------------*/

#define SQR(x)   ((x)*(x))

#define RE_GLO   6378136.0        /* radius of earth (m)            ref [2] */
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */

#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */

#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */

#define ERREPH_GLO 5.0            /* error of glonass ephemeris (m) */
#define TSTEP    60.0             /* integration step glonass ephemeris (s) */
#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */

#define DEFURASSR 0.15            /* default accuracy of ssr corr (m) */
#define MAXECORSSR 10.0           /* max orbit correction of ssr (m) */
#define MAXCCORSSR (1E-6*CLIGHT)  /* max clock correction of ssr (m) */
#define MAXAGESSR 90.0            /* max age of ssr orbit and clock (s) */
#define MAXAGESSR_HRCLK 10.0      /* max age of ssr high-rate clock (s) */
#define STD_BRDCCLK 30.0          /* error of broadcast clock (m) */
#define STD_GAL_NAPA 500.0        /* error of galileo ephemeris for NAPA (m) */

#define MAX_ITER_KEPLER 30        /* max number of iteration of Kepler */

/* ephemeris selections ------------------------------------------------------*/
/* =============================================================================
 * 变量：eph_sel
 * 功能：星历选择标志数组，分别对应 GPS, GLO, GAL, QZS, BDS, IRN, SBS。
 * ============================================================================= */
static int eph_sel[] = { /* GPS,GLO,GAL,QZS,BDS,IRN,SBS */
    0,0,0,0,0,0,0
};

/* =============================================================================
 * 函数：var_uraeph
 * 功能：将广播星历中播发的 URA（User Range Accuracy，用户测距精度）指数转换为方差。
 * 解释：卫星会在导航电文中告诉接收机“我现在的轨道大概有多准”。GPS 和伽利略的定义不同，
 * 这里根据不同的系统将整数指数转换为具体的误差方差 (m^2)，用于卡尔曼滤波的 R 阵加权。
 * ============================================================================= */
static double var_uraeph(int sys, int ura)
{
    const double ura_value[] = {
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
        3072.0,6144.0
    };
    if (sys == SYS_GAL) { /* Galileo 系统的 SISA 精度指标计算法 */
        if (ura <= 49) return SQR(ura * 0.01);
        if (ura <= 74) return SQR(0.5 + (ura - 50) * 0.02);
        if (ura <= 99) return SQR(1.0 + (ura - 75) * 0.04);
        if (ura <= 125) return SQR(2.0 + (ura - 100) * 0.16);
        return SQR(STD_GAL_NAPA);
    }
    else { /* GPS 系统的 URA 计算法 */
        return ura < 0 || 14 < ura ? SQR(6144.0) : SQR(ura_value[ura]);
    }
}

/* =============================================================================
 * 函数：var_urassr
 * 功能：将 SSR（状态空间表示，通常用于星星地通讯或精密单点定位）的 URA 转换为方差。
 * ============================================================================= */
static double var_urassr(int ura)
{
    double std;
    if (ura <= 0) return SQR(DEFURASSR);
    if (ura >= 63) return SQR(5.4665);
    std = (pow(3.0, (ura >> 3) & 7) * (1.0 + (ura & 7) / 4.0) - 1.0) * 1E-3;
    return SQR(std);
}

/* =============================================================================
 * 函数：alm2pos
 * 功能：根据历书（Almanac）计算卫星的粗略位置和钟差。
 * 解释：历书是精简版的星历，精度很低（几千米误差），但有效期长达几个月。
 * 主要用于接收机刚开机时“盲搜”天空中有哪些卫星，不用于最终的高精度定位。
 * ============================================================================= */
extern void alm2pos(gtime_t time, const alm_t* alm, double* rs, double* dts)
{
    double tk, M, E, Ek, sinE, cosE, u, r, i, O, x, y, sinO, cosO, cosi, mu;
    int n;

    char tstr[40];
    trace(4, "alm2pos : time=%s sat=%2d\n", time2str(time, tstr, 3), alm->sat);

    tk = timediff(time, alm->toa); /* 计算当前时间与历书参考时间(toa)的差值 */

    if (alm->A <= 0.0) { /* 半长轴无效，防错保护 */
        rs[0] = rs[1] = rs[2] = *dts = 0.0;
        return;
    }
    mu = satsys(alm->sat, NULL) == SYS_GAL ? MU_GAL : MU_GPS; /* 地球引力常数 */

    /* 迭代求解开普勒方程 (Kepler's Equation: E = M + e*sin(E)) 得到偏近点角 E */
    M = alm->M0 + sqrt(mu / (alm->A * alm->A * alm->A)) * tk;
    for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER && n < MAX_ITER_KEPLER; n++) {
        Ek = E; E -= (E - alm->e * sin(E) - M) / (1.0 - alm->e * cos(E));
    }
    if (n >= MAX_ITER_KEPLER) {
        trace(2, "alm2pos: kepler iteration overflow sat=%2d\n", alm->sat);
    }
    /* 计算轨道面内的极坐标并转换到地固坐标系 (ECEF) */
    sinE = sin(E); cosE = cos(E);
    u = atan2(sqrt(1.0 - alm->e * alm->e) * sinE, cosE - alm->e) + alm->omg;
    r = alm->A * (1.0 - alm->e * cosE);
    i = alm->i0;
    O = alm->OMG0 + (alm->OMGd - OMGE) * tk - OMGE * alm->toas;
    x = r * cos(u); y = r * sin(u); sinO = sin(O); cosO = cos(O); cosi = cos(i);
    rs[0] = x * cosO - y * cosi * sinO;
    rs[1] = x * sinO + y * cosi * cosO;
    rs[2] = y * sin(i);
    *dts = alm->f0 + alm->f1 * tk; /* 历书提供的简单一阶钟差修正 */
}

/* =============================================================================
 * 函数：eph2clk
 * 功能：利用广播星历计算卫星钟差（不包含相对论修正）。
 * 解释：通过星历中播发的卫星钟偏移(f0)、钟速(f1)、钟漂(f2)三个多项式系数，
 * 计算当前时刻卫星钟与系统时的偏差。
 * ============================================================================= */
extern double eph2clk(gtime_t time, const eph_t* eph)
{
    double t, ts;
    int i;

    char tstr[40];
    trace(4, "eph2clk : time=%s sat=%2d\n", time2str(time, tstr, 3), eph->sat);

    t = ts = timediff(time, eph->toc); /* 当前时间与星历时钟参考时间(toc)差值 */

    /* 钟差本身会影响信号传播时间的计算，这里用迭代法消除时间闭合差 */
    for (i = 0; i < 2; i++) {
        t = ts - (eph->f0 + eph->f1 * t + eph->f2 * t * t);
    }
    trace(4, "ephclk: t=%.12f ts=%.12f dts=%.12f f0=%.12f f1=%.9f f2=%.9f\n", t, ts,
        eph->f0 + eph->f1 * t + eph->f2 * t * t, eph->f0, eph->f1, eph->f2);

    return eph->f0 + eph->f1 * t + eph->f2 * t * t;
}

/* =============================================================================
 * 函数：eph2pos
 * 功能：利用广播星历计算卫星精确的三维坐标 (XYZ) 和 钟差。
 * 解释：这是标准单点定位最底层的核心函数！通过 16 个开普勒与摄动参数，完全还原
 * GPS/Galileo/BDS 卫星在太空中的绝对位置，并加入了狭义相对论效应改正。
 * ============================================================================= */
extern void eph2pos(gtime_t time, const eph_t* eph, double* rs, double* dts,
    double* var)
{
    double tk, M, E, Ek, sinE, cosE, u, r, i, O, sin2u, cos2u, x, y, sinO, cosO, cosi, mu, omge;
    double xg, yg, zg, sino, coso;
    int n, sys, prn;

    char tstr[40];
    trace(4, "eph2pos : time=%s sat=%2d\n", time2str(time, tstr, 3), eph->sat);

    if (eph->A <= 0.0) { /* 容错：无效星历 */
        rs[0] = rs[1] = rs[2] = *dts = *var = 0.0;
        return;
    }
    tk = timediff(time, eph->toe); /* 与星历参考时间(toe)的差值 */

    /* 区分不同导航系统的地球引力常数(mu)和地球自转角速度(omge) */
    switch ((sys = satsys(eph->sat, &prn))) {
    case SYS_GAL: mu = MU_GAL; omge = OMGE_GAL; break;
    case SYS_CMP: mu = MU_CMP; omge = OMGE_CMP; break;
    default:      mu = MU_GPS; omge = OMGE;     break;
    }
    /* 计算平近点角 M，并加上星历给定的摄动修正项 deln */
    M = eph->M0 + (sqrt(mu / (eph->A * eph->A * eph->A)) + eph->deln) * tk;

    /* 迭代求解开普勒方程得到偏近点角 E */
    for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER && n < MAX_ITER_KEPLER; n++) {
        Ek = E; E -= (E - eph->e * sin(E) - M) / (1.0 - eph->e * cos(E));
    }
    if (n >= MAX_ITER_KEPLER) {
        trace(2, "eph2pos: kepler iteration overflow sat=%2d\n", eph->sat);
    }
    sinE = sin(E); cosE = cos(E);

    trace(4, "kepler: sat=%2d e=%8.5f n=%2d del=%10.3e\n", eph->sat, eph->e, n, E - Ek);

    /* 计算真近点角，结合近地点幅角(omg)得到升交角距 u */
    u = atan2(sqrt(1.0 - eph->e * eph->e) * sinE, cosE - eph->e) + eph->omg;
    r = eph->A * (1.0 - eph->e * cosE);
    i = eph->i0 + eph->idot * tk;
    sin2u = sin(2.0 * u); cos2u = cos(2.0 * u);

    /* 加入二阶调和摄动修正项（极度重要，这是精密定位区分于历书的关键） */
    u += eph->cus * sin2u + eph->cuc * cos2u; /* 纬度幅角摄动 */
    r += eph->crs * sin2u + eph->crc * cos2u; /* 轨道半径摄动 */
    i += eph->cis * sin2u + eph->cic * cos2u; /* 轨道倾角摄动 */
    x = r * cos(u); y = r * sin(u); cosi = cos(i); /* 卫星在轨道面内的二维坐标 */

    /* 处理北斗 GEO (地球静止轨道) 卫星的特殊坐标系转换 */
    /* BDS的GEO卫星是在CGCS2000坐标系下播发的，且轨道面有一个 5 度的倾角，必须单独处理 */
    if (sys == SYS_CMP && (prn <= 5 || prn >= 59)) { /* ref [9] table 4-1 */
        O = eph->OMG0 + eph->OMGd * tk - omge * eph->toes;
        sinO = sin(O); cosO = cos(O);
        xg = x * cosO - y * cosi * sinO;
        yg = x * sinO + y * cosi * cosO;
        zg = y * sin(i);
        sino = sin(omge * tk); coso = cos(omge * tk);
        /* 乘上针对 GEO 卫星的特定欧拉角旋转矩阵 (-5度) */
        rs[0] = xg * coso + yg * sino * COS_5 + zg * sino * SIN_5;
        rs[1] = -xg * sino + yg * coso * COS_5 + zg * coso * SIN_5;
        rs[2] = -yg * SIN_5 + zg * COS_5;
    }
    else { /* 其他常规 MEO/IGSO 卫星的地固系(ECEF)转换 */
        O = eph->OMG0 + (eph->OMGd - omge) * tk - omge * eph->toes;
        sinO = sin(O); cosO = cos(O);
        rs[0] = x * cosO - y * cosi * sinO;
        rs[1] = x * sinO + y * cosi * cosO;
        rs[2] = y * sin(i);
    }
    tk = timediff(time, eph->toc);
    *dts = eph->f0 + eph->f1 * tk + eph->f2 * tk * tk;

    /* 相对论修正 (Relativity Correction)
     * 因为卫星速度极快且处于高空(低重力)，时间流逝跟地球表面不一样 */
    *dts -= 2.0 * sqrt(mu * eph->A) * eph->e * sinE / SQR(CLIGHT);

    *var = var_uraeph(sys, eph->sva); /* 获取卫星坐标的估计方差 */
}

/* =============================================================================
 * 函数：deq (Differential Equation)
 * 功能：GLONASS 轨道推算的微分方程组。
 * 解释：构建了包含地球中心引力，以及地球非球形引力项 (J2) 摄动的加速度公式。
 * ============================================================================= */
static void deq(const double* x, double* xdot, const double* acc)
{
    double a, b, c, r2 = dot3(x, x), r3 = r2 * sqrt(r2), omg2 = SQR(OMGE_GLO);

    if (r2 <= 0.0) {
        xdot[0] = xdot[1] = xdot[2] = xdot[3] = xdot[4] = xdot[5] = 0.0;
        return;
    }
    /* GLONASS PZ-90 坐标系下的 J2 摄动加速度计算 */
    a = 1.5 * J2_GLO * MU_GLO * SQR(RE_GLO) / r2 / r3; /* 3/2*J2*mu*Ae^2/r^5 */
    b = 5.0 * x[2] * x[2] / r2;                    /* 5*z^2/r^2 */
    c = -MU_GLO / r3 - a * (1.0 - b);                /* -mu/r^3-a(1-b) */
    xdot[0] = x[3]; xdot[1] = x[4]; xdot[2] = x[5]; /* 速度即为位置的导数 */
    /* 速度的导数 (加速度) 包含引力、离心力、科里奥利力以及播发在星历中的日月摄动力(acc) */
    xdot[3] = (c + omg2) * x[0] + 2.0 * OMGE_GLO * x[4] + acc[0];
    xdot[4] = (c + omg2) * x[1] - 2.0 * OMGE_GLO * x[3] + acc[1];
    xdot[5] = (c - 2.0 * a) * x[2] + acc[2];
}

/* =============================================================================
 * 函数：glorbit
 * 功能：使用 四阶龙格-库塔法 (Runge-Kutta 4th order, RK4) 对 GLONASS 轨道进行数值积分。
 * 解释：已知当前时间的速度和加速度，积分推导极短时间 (t) 后的速度和位置。
 * ============================================================================= */
static void glorbit(double t, double* x, const double* acc)
{
    double k1[6], k2[6], k3[6], k4[6], w[6];
    int i;

    /* RK4 的四次采样预估 */
    deq(x, k1, acc); for (i = 0; i < 6; i++) w[i] = x[i] + k1[i] * t / 2.0;
    deq(w, k2, acc); for (i = 0; i < 6; i++) w[i] = x[i] + k2[i] * t / 2.0;
    deq(w, k3, acc); for (i = 0; i < 6; i++) w[i] = x[i] + k3[i] * t;
    deq(w, k4, acc);
    /* 加权平均推算最终状态 */
    for (i = 0; i < 6; i++) x[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * t / 6.0;
}

/* =============================================================================
 * 函数：geph2clk
 * 功能：利用 GLONASS 广播星历计算卫星钟差。
 * 解释：GLONASS 星历参数命名不同于 GPS (不叫 f0/f1，叫 taun/gamn)。
 * ============================================================================= */
extern double geph2clk(gtime_t time, const geph_t* geph)
{
    double t, ts;
    int i;

    char tstr[40];
    trace(4, "geph2clk: time=%s sat=%2d\n", time2str(time, tstr, 3), geph->sat);

    t = ts = timediff(time, geph->toe);

    /* 迭代法消除时间闭合差 */
    for (i = 0; i < 2; i++) {
        t = ts - (-geph->taun + geph->gamn * t);
    }
    trace(4, "geph2clk: t=%.12f ts=%.12f taun=%.12f gamn=%.12f\n", t, ts, geph->taun,
        geph->gamn);
    return -geph->taun + geph->gamn * t;
}

/* =============================================================================
 * 函数：geph2pos
 * 功能：利用 GLONASS 广播星历计算卫星精确的三维坐标 (XYZ) 和 钟差。
 * 解释：读取参考历元 toe 时的初始状态 (pos, vel)，通过逐步累加步长 (TSTEP)
 * 调用 RK4 数值积分函数，一直递推计算到当前时刻 time 的卫星坐标。
 * ============================================================================= */
extern void geph2pos(gtime_t time, const geph_t* geph, double* rs, double* dts,
    double* var)
{
    double t, tt, x[6];
    int i;

    char tstr[40];
    trace(4, "geph2pos: time=%s sat=%2d\n", time2str(time, tstr, 3), geph->sat);

    t = timediff(time, geph->toe); /* 目标时间与参考历元的总时长 */

    *dts = -geph->taun + geph->gamn * t; /* 获取钟差 */
    trace(4, "geph2pos: sat=%d\n", geph->sat);

    /* 装载参考历元 toe 的卫星初始位置(x,y,z)与速度(vx,vy,vz) */
    for (i = 0; i < 3; i++) {
        x[i] = geph->pos[i];
        x[i + 3] = geph->vel[i];
    }
    /* 核心：数值积分逼近。如果时间相差较大，需要切分为多个 TSTEP(如60秒)的小步长逐步积分 */
    for (tt = t < 0.0 ? -TSTEP : TSTEP; fabs(t) > 1E-9; t -= tt) {
        if (fabs(t) < TSTEP) tt = t; /* 积分最后一步的余数处理 */
        glorbit(tt, x, geph->acc); /* 调用 RK4 进行一次状态递推 */
    }
    for (i = 0; i < 3; i++) rs[i] = x[i]; /* 取出最终积分后的 XYZ 坐标 */

    *var = SQR(ERREPH_GLO); /* 赋予 GLONASS 广播星历默认的经验方差 */
}
/* =============================================================================
 * 函数：seph2clk
 * 功能：利用 SBAS（星基增强系统）的星历计算卫星钟差。
 * 解释：SBAS 卫星（如 WAAS, EGNOS, BDSBAS 等）播发的星历格式与 GPS 不同。
 * 它的钟差模型只有两个参数：af0 (时间偏移) 和 af1 (时间漂移)。
 * ============================================================================= */
extern double seph2clk(gtime_t time, const seph_t* seph)
{
    char tstr[40];
    trace(4, "seph2clk: time=%s sat=%2d\n", time2str(time, tstr, 3), seph->sat);

    // 计算当前时间与星历参考时间 (t0) 的时间差 ts
    double ts = timediff(time, seph->t0), t = ts;

    // 迭代 2 次以消除钟差引起的时间闭合差 (光速传播时间)
    for (int i = 0; i < 2; i++) {
        t = ts - (seph->af0 + seph->af1 * t);
    }
    return seph->af0 + seph->af1 * t; // 返回秒为单位的钟差
}

/* =============================================================================
 * 函数：seph2pos
 * 功能：利用 SBAS 星历计算卫星三维坐标 (XYZ) 和钟差。
 * 解释：SBAS 卫星本身就是地球同步轨道卫星(GEO)，它的轨道不像 GPS 播发复杂的开普勒参数，
 * 而是直接播发某一时刻的位置、速度和加速度，通过简单的牛顿运动学公式(匀加速直线运动模型)推算位置。
 * ============================================================================= */
extern void seph2pos(gtime_t time, const seph_t* seph, double* rs, double* dts,
    double* var)
{
    double t;
    int i;

    char tstr[40];
    trace(4, "seph2pos: time=%s sat=%2d\n", time2str(time, tstr, 3), seph->sat);

    t = timediff(time, seph->t0); // 距离参考历元的时间差

    // S = S0 + V*t + 1/2*A*t^2 (基于位置、速度、加速度的二次抛物线预报)
    for (i = 0; i < 3; i++) {
        rs[i] = seph->pos[i] + seph->vel[i] * t + seph->acc[i] * t * t / 2.0;
    }
    *dts = seph->af0 + seph->af1 * t; // 附带计算钟差

    *var = var_uraeph(SYS_SBS, seph->sva); // 转换 URA 为方差
}

/* =============================================================================
 * 函数：seleph
 * 功能：在所有广播星历中，为指定的 GPS/Galileo/BDS/QZSS 卫星挑选“最合适”的星历。
 * 解释：由于导航电文每隔 2 小时才更新一次，程序需要遍历内存中所有的星历块，
 * 寻找星历参考时间 (toe) 距离当前观测时间最近的那一块，并且判断是否过期 (tmax)。
 * ============================================================================= */
static eph_t* seleph(gtime_t time, int sat, int iode, const nav_t* nav)
{
    double t, tmax, tmin;
    int i, j = -1, sys, sel = 0;

    char tstr[40];
    trace(4, "seleph  : time=%s sat=%2d iode=%d\n", time2str(time, tstr, 3), sat, iode);

    sys = satsys(sat, NULL);
    // 1. 设置不同系统星历的“保质期”(tmax)。例如 GPS 星历通常有效期为 2 小时 (MAXDTOE)。
    switch (sys) {
    case SYS_GPS: tmax = MAXDTOE + 1.0; sel = eph_sel[0]; break;
    case SYS_GAL: tmax = MAXDTOE_GAL; sel = eph_sel[2]; break;
    case SYS_QZS: tmax = MAXDTOE_QZS + 1.0; sel = eph_sel[3]; break;
    case SYS_CMP: tmax = MAXDTOE_CMP + 1.0; sel = eph_sel[4]; break;
    case SYS_IRN: tmax = MAXDTOE_IRN + 1.0; sel = eph_sel[5]; break;
    default: tmax = MAXDTOE + 1.0; break;
    }
    tmin = tmax + 1.0;

    // 2. 遍历内存中已解析的所有星历数据 (nav->eph 数组)
    for (i = 0; i < nav->n; i++) {
        if (nav->eph[i].sat != sat) continue; // 卫星号不匹配，跳过
        if (iode >= 0 && nav->eph[i].iode != iode) continue; // 如果指定了星历期号(IODE)，必须严格匹配

        // Galileo 系统的特殊电文选择逻辑 (I/NAV 或 F/NAV)
        if (sys == SYS_GAL) {
            sel = getseleph(SYS_GAL);
            if (sel == 1 && !(nav->eph[i].code & (1 << 9))) continue; /* I/NAV */
            if (sel == 2 && !(nav->eph[i].code & (1 << 8))) continue; /* F/NAV */
            if (timediff(nav->eph[i].toe, time) >= 0.0) continue; /* AOD<=0 */
        }

        // 3. 检查星历是否过期
        if ((t = fabs(timediff(nav->eph[i].toe, time))) > tmax) continue;
        if (iode >= 0) return nav->eph + i;

        // 4. 核心：不断记录与当前观测时间差值 t 最小的星历索引 j
        if (t <= tmin) { j = i; tmin = t; } /* toe closest to time */
    }
    if (iode >= 0 || j < 0) {
        trace(2, "no broadcast ephemeris: %s sat=%2d iode=%3d\n", time2str(time, tstr, 0),
            sat, iode);
        return NULL; // 没找到可用星历！
    }
    trace(4, "seleph: sat=%d dt=%.0f\n", sat, tmin);
    return nav->eph + j; // 返回找到的最佳星历指针
}

/* =============================================================================
 * 函数：selgeph
 * 功能：为 GLONASS 卫星挑选最近的广播星历。
 * 解释：逻辑与 seleph 类似，但针对的是 nav->geph 数组 (GLONASS 星历结构体独立)。
 * ============================================================================= */
static geph_t* selgeph(gtime_t time, int sat, int iode, const nav_t* nav)
{
    double t, tmax = MAXDTOE_GLO, tmin = tmax + 1.0;
    int i, j = -1;

    char tstr[40];
    trace(4, "selgeph : time=%s sat=%2d iode=%2d\n", time2str(time, tstr, 3), sat, iode);

    for (i = 0; i < nav->ng; i++) {
        if (nav->geph[i].sat != sat) continue;
        if (iode >= 0 && nav->geph[i].iode != iode) continue;
        if ((t = fabs(timediff(nav->geph[i].toe, time))) > tmax) continue; // 判断过期
        if (iode >= 0) return nav->geph + i;
        if (t <= tmin) { j = i; tmin = t; } /* toe closest to time */
    }
    if (iode >= 0 || j < 0) {
        trace(3, "no glonass ephemeris  : %s sat=%2d iode=%2d\n", time2str(time, tstr, 0),
            sat, iode);
        return NULL;
    }
    trace(4, "selgeph: sat=%d dt=%.0f\n", sat, tmin);
    return nav->geph + j;
}

/* =============================================================================
 * 函数：selseph
 * 功能：为 SBAS 卫星挑选最近的广播星历。
 * 解释：逻辑同上，针对 nav->seph 数组。
 * ============================================================================= */
static seph_t* selseph(gtime_t time, int sat, const nav_t* nav)
{
    double t, tmax = MAXDTOE_SBS, tmin = tmax + 1.0;
    int i, j = -1;

    char tstr[40];
    trace(4, "selseph : time=%s sat=%2d\n", time2str(time, tstr, 3), sat);

    for (i = 0; i < nav->ns; i++) {
        if (nav->seph[i].sat != sat) continue;
        if ((t = fabs(timediff(nav->seph[i].t0, time))) > tmax) continue;
        if (t <= tmin) { j = i; tmin = t; } /* toe closest to time */
    }
    if (j < 0) {
        trace(3, "no sbas ephemeris     : %s sat=%2d\n", time2str(time, tstr, 0), sat);
        return NULL;
    }
    return nav->seph + j;
}

/* =============================================================================
 * 函数：ephclk
 * 功能：利用广播星历计算单颗卫星钟差的“总调度器”。
 * 解释：根据卫星系统(sys)，调用对应的星历选择函数(seleph/selgeph等)，
 * 并转交到底层特定的钟差计算公式(eph2clk/geph2clk等)。
 * ============================================================================= */
static int ephclk(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
    double* dts)
{
    eph_t* eph;
    geph_t* geph;
    seph_t* seph;
    int sys;

    char tstr[40];
    trace(4, "ephclk  : time=%s sat=%2d\n", time2str(time, tstr, 3), sat);

    sys = satsys(sat, NULL); // 获取卫星所属系统

    if (sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS || sys == SYS_CMP || sys == SYS_IRN) {
        if (!(eph = seleph(teph, sat, -1, nav))) return 0;
        *dts = eph2clk(time, eph);
    }
    else if (sys == SYS_GLO) {
        if (!(geph = selgeph(teph, sat, -1, nav))) return 0;
        if (fabs(geph->taun) > 1) return 0; /* 剔除非法数据防止浮点数异常 */
        *dts = geph2clk(time, geph);
    }
    else if (sys == SYS_SBS) {
        if (!(seph = selseph(teph, sat, nav))) return 0;
        *dts = seph2clk(time, seph);
    }
    else return 0; // 不支持的系统

    return 1;
}

/* =============================================================================
 * 函数：ephpos
 * 功能：利用广播星历计算卫星坐标与钟差的“总调度器”，并粗略估算卫星速度。
 * 解释：它除了计算位置，还会将时间人为向后拨 1 毫秒 (tt=1E-3) 再次计算坐标，
 * 进而利用微分的思想 (差分求导) 粗略估算出卫星的三维运动速度。
 * ============================================================================= */
static int ephpos(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
    int iode, double* rs, double* dts, double* var, int* svh)
{
    eph_t* eph;
    geph_t* geph;
    seph_t* seph;
    double rst[3], dtst[1], tt = 1E-3;
    int i, sys;

    char tstr[40];
    trace(4, "ephpos  : time=%s sat=%2d iode=%d\n", time2str(time, tstr, 3), sat, iode);

    sys = satsys(sat, NULL);
    *svh = -1; // 卫星健康状态初始化为非法

    if (sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS || sys == SYS_CMP || sys == SYS_IRN) {
        if (!(eph = seleph(teph, sat, iode, nav))) return 0;
        eph2pos(time, eph, rs, dts, var); // 计算当前时刻的坐标 rs 和 钟差 dts

        /* 差分求导：计算 1 毫秒后的位置，用于估算速度 */
        time = timeadd(time, tt);
        eph2pos(time, eph, rst, dtst, var);
        *svh = eph->svh; // 记录卫星健康标志
    }
    else if (sys == SYS_GLO) {
        if (!(geph = selgeph(teph, sat, iode, nav))) return 0;
        geph2pos(time, geph, rs, dts, var);
        time = timeadd(time, tt);
        geph2pos(time, geph, rst, dtst, var);
        *svh = geph->svh;
    }
    else if (sys == SYS_SBS) {
        if (!(seph = selseph(teph, sat, nav))) return 0;
        seph2pos(time, seph, rs, dts, var);
        time = timeadd(time, tt);
        seph2pos(time, seph, rst, dtst, var);
        *svh = seph->svh;
    }
    else return 0;

    /* satellite velocity and clock drift by differential approx */
    /* 通过 (下一时刻位置 - 当前位置) / 时间差 得到速度。
       之所以这么做，是因为广播星历不直接播发速度，且开普勒方程直接求导十分复杂 */
    for (i = 0; i < 3; i++) rs[i + 3] = (rst[i] - rs[i]) / tt;
    dts[1] = (dtst[0] - dts[0]) / tt; // 钟漂(钟速)

    return 1; // 成功
}

/* =============================================================================
 * 函数：satpos_sbas
 * 功能：利用 SBAS 改正数计算更高精度的卫星位置和钟差。
 * 解释：首先用常规广播星历算出基础坐标，如果当前卫星存在 SBAS（星基增强）下发的
 * 长/快周期改正数据，则将其叠加到基础坐标和钟差上。
 * ============================================================================= */
static int satpos_sbas(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
    double* rs, double* dts, double* var, int* svh)
{
    const sbssatp_t* sbs = NULL;
    int i;

    char tstr[40];
    trace(4, "satpos_sbas: time=%s sat=%2d\n", time2str(time, tstr, 3), sat);

    /* 1. 在 nav 结构体中搜索是否包含这颗卫星的 SBAS 改正信息 */
    for (i = 0; i < nav->sbssat.nsat; i++) {
        sbs = nav->sbssat.sat + i;
        if (sbs->sat == sat) break;
    }
    if (i >= nav->sbssat.nsat) {
        trace(2, "no sbas, use brdcast: %s sat=%2d\n", time2str(time, tstr, 0), sat);
        // 如果没有 SBAS 数据，退化回使用纯广播星历 (ephpos)
        if (!ephpos(time, teph, sat, nav, -1, rs, dts, var, svh)) return 0;
        return 1;
    }

    /* 2. 首先用 SBAS 改正数据指向的星历期号(IODE)，算出一个基础广播星历轨道 */
    if (!ephpos(time, teph, sat, nav, sbs->lcorr.iode, rs, dts, var, svh)) return 0;

    /* 3. 应用 SBAS 卫星改正模型 (叠加长期轨道改正和快变钟差改正) */
    if (sbssatcorr(time, sat, nav, rs, dts, var)) return 1;
    *svh = -1;
    return 0;
}

/* =============================================================================
 * 函数：satpos_ssr
 * 功能：利用 SSR（状态空间表示）改正数计算精密卫星位置和钟差（实时 PPP 核心！）。
 * 解释：SSR 改正数通常通过 NTRIP 网络流实时下发。这种模式允许我们在没有静态精密星历
 * (.SP3) 的情况下，也能利用广播星历 + SSR 实时算出厘米级的卫星轨道，从而实现实时 PPP。
 * SSR 播发的是径向(Radial)、切向(Along)、法向(Cross)的三维轨道残差，以及高频钟差残差。
 * ============================================================================= */
static int satpos_ssr(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
    int opt, double* rs, double* dts, double* var, int* svh)
{
    const ssr_t* ssr;
    eph_t* eph;
    double t1, t2, t3, er[3], ea[3], ec[3], rc[3], deph[3], dclk, dant[3] = { 0 }, tk;
    int i, sys;

    char tstr[40];
    trace(4, "satpos_ssr: time=%s sat=%2d\n", time2str(time, tstr, 3), sat);

    ssr = nav->ssr + sat - 1;

    // 检查是否有有效的 SSR 轨道和钟差更新
    if (!ssr->t0[0].time) {
        trace(2, "no ssr orbit correction: %s sat=%2d\n", time2str(time, tstr, 0), sat);
        return 0;
    }
    if (!ssr->t0[1].time) {
        trace(2, "no ssr clock correction: %s sat=%2d\n", time2str(time, tstr, 0), sat);
        return 0;
    }
    /* 保证轨道改正和钟差改正的星历期号 (IOD) 是配套匹配的 */
    if (ssr->iod[0] != ssr->iod[1]) {
        trace(2, "inconsist ssr correction: %s sat=%2d iod=%d %d\n",
            time2str(time, tstr, 0), sat, ssr->iod[0], ssr->iod[1]);
        *svh = -1;
        return 0;
    }

    // 计算当前时间与 SSR 历元的时差
    t1 = timediff(time, ssr->t0[0]);
    t2 = timediff(time, ssr->t0[1]);
    t3 = timediff(time, ssr->t0[2]);

    /* 如果数据超期 (Max Age)，则丢弃 */
    if (fabs(t1) > MAXAGESSR || fabs(t2) > MAXAGESSR) {
        trace(2, "age of ssr error: %s sat=%2d t=%.0f %.0f\n", time2str(time, tstr, 0),
            sat, t1, t2);
        *svh = -1;
        return 0;
    }
    if (ssr->udi[0] >= 1.0) t1 -= ssr->udi[0] / 2.0;
    if (ssr->udi[1] >= 1.0) t2 -= ssr->udi[1] / 2.0;

    // 1. 利用时间漂移率，推算当前精确的轨道改正量 (deph) 和钟差改正量 (dclk)
    for (i = 0; i < 3; i++) deph[i] = ssr->deph[i] + ssr->ddeph[i] * t1;
    dclk = ssr->dclk[0] + ssr->dclk[1] * t2 + ssr->dclk[2] * t2 * t2;

    /* 附加高频钟差改正 (HRCLK) */
    if (ssr->iod[0] == ssr->iod[2] && ssr->t0[2].time && fabs(t3) < MAXAGESSR_HRCLK) {
        dclk += ssr->hrclk;
    }
    // 粗差防爆检查
    if (norm(deph, 3) > MAXECORSSR || fabs(dclk) > MAXCCORSSR) {
        trace(3, "invalid ssr correction: %s deph=%.1f dclk=%.1f\n",
            time2str(time, tstr, 0), norm(deph, 3), dclk);
        *svh = -1;
        return 0;
    }

    /* 2. 调用 ephpos 计算基础的广播星历卫星位置 rs */
    if (!ephpos(time, teph, sat, nav, ssr->iode, rs, dts, var, svh)) return 0;

    /* 重新计算相对论效应 (抵消广播星历里自带的狭义相对论钟差) */
    sys = satsys(sat, NULL);
    if (sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS || sys == SYS_CMP) {
        if (!(eph = seleph(teph, sat, ssr->iode, nav))) return 0;

        tk = timediff(time, eph->toc);
        dts[0] = eph->f0 + eph->f1 * tk + eph->f2 * tk * tk;
        dts[1] = eph->f1 + 2.0 * eph->f2 * tk;

        dts[0] -= 2.0 * dot3(rs, rs + 3) / CLIGHT / CLIGHT;
    }

    /* 3. 坐标转换：将卫星随体系 (径向ea, 切向er, 法向ec) 转换为地固坐标系 (ECEF) */
    if (!normv3(rs + 3, ea)) return 0;
    cross3(rs, rs + 3, rc);
    if (!normv3(rc, ec)) {
        *svh = -1;
        return 0;
    }
    cross3(ea, ec, er);

    /* SSR 天线相位中心偏移 (PCO) 改正 (通常用精密星历或 SSR 时要加 APC 改正) */
    if (opt) {
        satantoff(time, rs, sat, nav, dant);
    }

    /* 4. 将计算好的三维轨道残差 (deph) 投影到地固系，并叠加到基础轨道 rs 上！ */
    for (i = 0; i < 3; i++) {
        rs[i] += -(er[i] * deph[0] + ea[i] * deph[1] + ec[i] * deph[2]) + dant[i];
    }
    /* t_corr = t_sv - (dts(brdc) + dclk(ssr) / CLIGHT) */
    dts[0] += dclk / CLIGHT; // 加上精密钟差改正

    /* 方差降级：用 SSR 指定的精度(URA)替代广播星历的粗劣精度 */
    *var = var_urassr(ssr->ura);

    trace(5, "satpos_ssr: %s sat=%2d deph=%6.3f %6.3f %6.3f er=%6.3f %6.3f %6.3f dclk=%6.3f var=%6.3f\n",
        time2str(time, tstr, 2), sat, deph[0], deph[1], deph[2], er[0], er[1], er[2], dclk, *var);

    return 1;
}
/* =============================================================================
 * 函数：satpos (Satellite Position)
 * 功能：单颗卫星坐标与钟差计算的“终极路由器”。
 * 解释：它不负责具体的数学计算，而是根据你在配置中设定的 ephopt（星历选项），
 * 决定把计算任务分发给哪个底层函数（广播星历、SBAS、SSR 还是精密星历 SP3）。
 * ============================================================================= */
extern int satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
    const nav_t* nav, double* rs, double* dts, double* var,
    int* svh)
{
    char tstr[40];
    trace(4, "satpos  : time=%s sat=%2d ephopt=%d\n", time2str(time, tstr, 3), sat, ephopt);

    *svh = 0;

    /* 核心路由分发器：根据设定的星历类型调用不同的底层算法 */
    switch (ephopt) {
    case EPHOPT_BRDC: return ephpos(time, teph, sat, nav, -1, rs, dts, var, svh); /* 广播星历 */
    case EPHOPT_SBAS: return satpos_sbas(time, teph, sat, nav, rs, dts, var, svh); /* SBAS增强星历 */
    case EPHOPT_SSRAPC: return satpos_ssr(time, teph, sat, nav, 0, rs, dts, var, svh); /* SSR状态空间改正(天线相位中心) */
    case EPHOPT_SSRCOM: return satpos_ssr(time, teph, sat, nav, 1, rs, dts, var, svh); /* SSR状态空间改正(卫星质心) */
    case EPHOPT_PREC:
        /* 精密星历 (PPP 模式的核心) */
        if (!peph2pos(time, sat, nav, 1, rs, dts, var)) break; else return 1;
    }
    *svh = -1;
    return 0; /* 如果都不匹配或计算失败，返回 0 */
}

/* =============================================================================
 * 函数：satposs (Satellite Positions - 注意末尾有s，代表复数)
 * 功能：批量计算当前历元所有可见卫星的位置和钟差（极其核心的物理时间回溯）。
 * 解释：接收机记录的观测时间是“信号到达接收机的时间”。但卫星一直以几公里每秒的速度在飞，
 * 我们必须知道信号“发射”时卫星到底在哪里。这个函数利用伪距除以光速，时光倒流推算出
 * 发射时刻 (Transmission Time)，再调用 satpos 算出该时刻的精确卫星位置。
 * ============================================================================= */
extern void satposs(gtime_t teph, const obsd_t* obs, int n, const nav_t* nav,
    int ephopt, double* rs, double* dts, double* var, int* svh)
{
    gtime_t time[2 * MAXOBS] = { {0} };
    double dt, pr;
    int i, j;

    char tstr[40];
    trace(3, "satposs : teph=%s n=%d ephopt=%d\n", time2str(teph, tstr, 3), n, ephopt);

    /* 遍历当前历元接收机观测到的所有卫星 (n 颗) */
    for (i = 0; i < n && i < 2 * MAXOBS; i++) {
        for (j = 0; j < 6; j++) rs[j + i * 6] = 0.0;
        for (j = 0; j < 2; j++) dts[j + i * 2] = 0.0;
        var[i] = 0.0; svh[i] = 0;

        /* 1. 搜索任意一个有效的伪距观测值 (用于计算信号飞行时间) */
        for (j = 0, pr = 0.0; j < NFREQ; j++) if ((pr = obs[i].P[j]) != 0.0) break;

        /* 如果连伪距都没有，无法计算飞行时间，直接跳过该卫星 */
        if (j >= NFREQ) {
            trace(2, "no pseudorange %s sat=%2d\n", time2str(obs[i].time, tstr, 3), obs[i].sat);
            continue;
        }

        /* 2. 时光倒流：发射时间 = 接收时间 - 伪距/光速 */
        time[i] = timeadd(obs[i].time, -pr / CLIGHT);

        /* 3. 修正钟差导致的发射时间偏差 (用广播星历粗略估算当前时刻的钟差) */
        if (!ephclk(time[i], teph, obs[i].sat, nav, &dt)) {
            trace(3, "no broadcast clock %s sat=%2d\n", time2str(time[i], tstr, 3), obs[i].sat);
            continue;
        }
        time[i] = timeadd(time[i], -dt); /* 得到真正精确的信号发射时间 (Transmission Time) */

        /* 4. 真正的高潮：调用上面的 satpos，算出该精确发射时刻下，卫星在宇宙中的三维绝对坐标 */
        if (!satpos(time[i], teph, obs[i].sat, ephopt, nav, rs + i * 6, dts + i * 2, var + i,
            svh + i)) {
            trace(3, "no ephemeris %s sat=%2d\n", time2str(time[i], tstr, 3), obs[i].sat);
            continue;
        }

        /* 5. 补丁机制：如果使用的是精密星历，但当前这颗卫星精密钟差缺失，退化使用广播星历的钟差保底 */
        if (dts[i * 2] == 0.0) {
            if (!ephclk(time[i], teph, obs[i].sat, nav, dts + i * 2)) continue;
            dts[1 + i * 2] = 0.0;
            *var = SQR(STD_BRDCCLK); /* 并且将方差放大，标记其精度不佳 */
        }
        trace(4, "satposs: %d,time=%.9f dt=%.9f pr=%.3f rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f\n",
            obs[i].sat, time[i].sec, dt, pr, rs[i * 6], rs[1 + i * 6], rs[2 + i * 6], dts[i * 2] * 1E9,
            var[i]);
    }

    /* 打印结果到日志 */
    for (i = 0; i < n && i < 2 * MAXOBS; i++) {
        trace(4, "%s sat=%2d rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f svh=%02X\n",
            time2str(time[i], tstr, 9), obs[i].sat, rs[i * 6], rs[1 + i * 6], rs[2 + i * 6],
            dts[i * 2] * 1E9, var[i], svh[i]);
    }
}

/* =============================================================================
 * 函数：setseleph
 * 功能：设置系统针对多种格式电文时的偏好策略。
 * 解释：现代导航系统往往播发多种电文（例如 GPS 的旧版 LNAV 和新版 CNAV，
 * 伽利略的 I/NAV 和 F/NAV）。这个函数用于全局设定遇到多类型电文时优先用哪个。
 * ============================================================================= */

// *-----------------------------------------------------------------------------*/
extern void setseleph(int sys, int sel)
{
    switch (sys) {
    case SYS_GPS: eph_sel[0] = sel; break;
    case SYS_GLO: eph_sel[1] = sel; break;
    case SYS_GAL: eph_sel[2] = sel; break;
    case SYS_QZS: eph_sel[3] = sel; break;
    case SYS_CMP: eph_sel[4] = sel; break;
    case SYS_IRN: eph_sel[5] = sel; break;
    case SYS_SBS: eph_sel[6] = sel; break;
    }
}

/* =============================================================================
 * 函数：getseleph
 * 功能：获取当前设定的电文类型偏好。
 * ============================================================================= */
 //*-----------------------------------------------------------------------------*/
extern int getseleph(int sys)
{
    switch (sys) {
    case SYS_GPS: return eph_sel[0];
    case SYS_GLO: return eph_sel[1];
    case SYS_GAL: return eph_sel[2];
    case SYS_QZS: return eph_sel[3];
    case SYS_CMP: return eph_sel[4];
    case SYS_IRN: return eph_sel[5];
    case SYS_SBS: return eph_sel[6];
    }
    return 0;
}