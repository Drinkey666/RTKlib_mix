/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL  use IERS tide model
*           -DOUTSTAT_AMB output ambiguity parameters to solution status
*
* references :
*    [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*    [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*        2003, November 2003
*    [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*        Space Technology Library, 2004
*    [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*        May 2009
*    [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*        Code Biases, URA
*    [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*        celestial reference frames, Geophys. Res. Let., 1997
*    [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*    [8] J.Kouba, A simplified yaw-attitude model for eclipsing GPS satellites,
*        GPS Solutions, 13:1-12, 2009
*    [9] F.Dilssner, GPS IIF-1 satellite antenna phase center and attitude
*        modeling, InsideGNSS, September, 2010
*    [10] F.Dilssner, The GLONASS-M satellite yaw-attitude model, Advances in
*        Space Research, 2010
*    [11] IGS MGEX (http://igs.org/mgex)
*
* version : $Revision:$ $Date:$
* history : 2010/07/20 1.0  new
*                           added api:
*                               tidedisp()
*           2010/12/11 1.1  enable exclusion of eclipsing satellite
*           2012/02/01 1.2  add gps-glonass h/w bias correction
*                           move windupcorr() to rtkcmn.c
*           2013/03/11 1.3  add otl and pole tides corrections
*                           involve iers model with -DIERS_MODEL
*                           change initial variances
*                           suppress acos domain error
*           2013/09/01 1.4  pole tide model by iers 2010
*                           add mode of ionosphere model off
*           2014/05/23 1.5  add output of trop gradient in solution status
*           2014/10/13 1.6  fix bug on P0(a[3]) computation in tide_oload()
*                           fix bug on m2 computation in tide_pole()
*           2015/03/19 1.7  fix bug on ionosphere correction for GLO and BDS
*           2015/05/10 1.8  add function to detect slip by MW-LC jump
*                           fix ppp solution problem with large clock variance
*           2015/06/08 1.9  add precise satellite yaw-models
*                           cope with day-boundary problem of satellite clock
*           2015/07/31 1.10 fix bug on nan-solution without glonass nav-data
*                           pppoutsolsat() -> pppoutstat()
*           2015/11/13 1.11 add L5-receiver-dcb estimation
*                           merge post-residual validation by rnx2rtkp_test
*                           support support option opt->pppopt=-GAP_RESION=nnnn
*           2016/01/22 1.12 delete support for yaw-model bug
*                           add support for ura of ephemeris
*           2018/10/10 1.13 support api change of satexclude()
*           2020/11/30 1.14 use sat2freq() to get carrier frequency
*                           use E1-E5b for Galileo iono-free LC
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* =============================================================================
 * 中文阅读指南：当前项目的手机多系统、多频非组合 PPP
 *
 * 一、输入和职责
 *   rtkpos() 先组织观测并执行单点定位，再按配置调用本文件的 pppos()。
 *   本文件不读取 Android TXT/RINEX 文件，也不下载产品；输入是同一历元的
 *   obsd_t[] 和已加载的 nav_t。TXT 与 RINEX 入口最终使用同一个 PPP 模型。
 *   rtk_t 必须跨历元保存：x/P 是滤波状态/协方差，ssat 保存周跳、平滑和
 *   信号质量历史。不要把“每历元初始化某个状态”误读成重新 rtkinit()。
 *
 * 二、常用单位和下标（阅读观测方程前先确认）
 *   obs.P：米；obs.L：周；obs.D：Hz；距离变化率为 -lambda*D。
 *   corr_meas() 后的 P、L 都是米，模糊度状态也是米而不是整周数。
 *   obs.time：GPST；pos[0/1]：纬度/经度弧度；pos[2]：椭球高米。
 *   obs.SNR*SNR_UNIT：C/N0，单位 dB-Hz；azel：方位角/高度角弧度。
 *   sat 从 1 开始，ssat[sat-1] 从 0 开始；f/frq 是从 0 开始的存储槽。
 *   槽号不是信号身份：必须同时看 constellation 和 obs.code[f]。
 *   P 矩阵按列存储，P[row+col*nx]；对角元素是方差，不是标准差。
 *
 * 三、建议阅读顺序
 *   1. 状态数量/索引宏及 ppp_ifb_index()：弄清每个参数的含义。
 *   2. pppos()：看一个历元的完整调用顺序。
 *   3. precheck_*、uddoppsm_ppp()：看粗差、周跳和因果伪距平滑。
 *   4. udstate_ppp()：看位置/钟差/大气/RCB/模糊度预测与初始化。
 *   5. corr_meas()、model_trop()/model_iono()：看单位和模型改正。
 *   6. ppp_res()：看 GNSS 残差、H/R、抗差及外部大气伪观测。
 *   7. update_stat()/ppp_diag_emit()：看 Q 值、输出及诊断含义。
 *
 * 四、解释当前实现时必须注意
 *   对流层状态存的是总 ZTD；实际估计湿延迟为 ZTD-ZHD。VMF3 用于
 *   映射、ZHD 和初始化；-VMF3ZWDSIG>0 时还追加 ZWD 软约束。
 *   电离层状态存的是 L1 参考频率的垂直延迟（米），不是 TECU；
 *   IONEX 是有方差的软先验，不能当作无误差真值重复使用。
 *   GPS 钟差每历元重置；其他系统存 GPS 参考下的 ISB，跨历元连续。
 *   RCB 是接收机“实际信号”的相对码偏差；OSB 是卫星产品，两者不同。
 *   Q6 表示载波参与的 PPP 状态，不保证厘米级；Q5 是单点/降级状态。
 *   坐标标准差反映模型内不确定度，不能代替与独立真值的误差评估。
 *
 * 本次中文说明只解释现有代码；不改变任何模型、参数或数值运算。
 * ============================================================================= */

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

/* 残差编辑最多试 16 次；每次从同一预测状态重新滤波，不是 16 个新历元。
 * 4 sigma 启动方差膨胀，8 sigma 进入验后硬拒绝候选；用膨胀前 sigma 判断。
 * VAR_* 是初始方差，prn[] 是过程噪声参数；二者不能混为观测定权。 */
#define MAX_ITER    16              /* max post-fit editing iterations (phone multi-frequency PPP) */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4              /* min satellite number for solution */
#define THRES_ROBUST 4.0            /* start robust down-weighting (sigma) */
#define THRES_REJECT 8.0            /* hard reject threshold after robust weighting (raw sigma) */

#define THRES_MW_JUMP 10.0

#define VAR_POS     SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_VEL     SQR(10.0)       /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)       /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) */
#define VAR_ZTD     SQR(0.15)       /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_IFB     SQR(10.0)       /* init variance common F3 code bias (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO SQR(10.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */

/* IONEX-constrained ionosphere estimation -----------------------------------
 * The PPP ionosphere state is vertical L1 delay (m). IONEX is used as a soft
 * prior, not as an error-free correction. -IONCONS= controls the 1-sigma
 * floor of the vertical prior in metres. -IONCONSINT= is the effective
 * correlation time (s) of successive constraints from the same GIM. */
#define ION_CONSTR_SIGMA_DEF 1.5
#define ION_CONSTR_SIGMA_MIN 0.20
#define ION_CONSTR_SIGMA_MAX 10.0
#define ION_CONSTR_INTERVAL_DEF 1.0
#define ION_CONSTR_ROBUST    3.0
#define ION_CONSTR_MAXPEN    25.0

/* VMF3 对流层辅助：提供 ZHD、映射系数及 ZTD 初值。
 * GNSS 方程中的湿延迟始终来自 PPP 状态；是否持续加入产品 ZWD 软约束，
 * 由 -VMF3ZWDSIG 决定（<=0 关闭）。不能只看已加载文件就认定约束生效。 */
#define VMF3_NLAT            36
#define VMF3_NLON            72
#define VMF3_NGRID           (VMF3_NLAT*VMF3_NLON)
#define VMF3_INIT_SIGMA_DEF  0.20
#define VMF3_INIT_SIGMA_MIN  0.05
#define VMF3_INIT_SIGMA_MAX  2.00
#define VMF3_ZWD_ROBUST      3.0
#define VMF3_ZWD_MAXPEN      25.0

#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */

/* 状态向量布局（所有索引由 opt 决定，禁止凭固定数字访问）：
 * [位置/速度/加速度][GPS 钟差及 ISB][ZTD/梯度][每星电离层][RCB][每星每槽模糊度]
 * NP：静态或无动态模型为 3，开启 dynamics 为 9；NC 固定 5。
 * NT：不估计为 0，EST 为 1，ESTG 为 3；NI 仅 IONOOPT_EST 时为 MAXSAT。
 * ND：nf>=4 为 5，nf>=3 为 4，否则为 0；NB 按有效频槽和 MAXSAT 分配。
 * 分配了 MAXSAT 个状态不代表所有卫星都参与更新；零值/无方差状态可不活跃。
 * IC/IT/II/ID/IB 分别定位钟差、对流层、电离层、RCB 和模糊度。
 * IFLC 模式只估一个组合模糊度，因此 NF(opt)=1，不能直接等同 opt->nf。 */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)

/*
 * PPP receiver-clock states use fixed time-system slots. Do not size this
 * block with NSYS: NSYS changes with compile-time ENAxxx flags, while the
 * measurement model historically used fixed GPS/GLO/GAL/BDS/IRN indices.
 * A compact NSYS block therefore aliases clock states with trop/iono states
 * whenever (for example) GLONASS is compiled out.
 */
#define NCLK_PPP    5               /* GPS/QZS, GLO, GAL, BDS, IRN */
#define NC(opt)     (NCLK_PPP)

#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)

/*
 * Receiver code biases are estimated only relative to each constellation's
 * reference code (GPS/GAL C1C, BDS C2I). Estimating a state for a reference
 * code would be indistinguishable from its receiver-clock state. The original
 * four states match the three-frequency phone profile. B1C
 * has its own state relative to the B1I reference code, not a replacement
 * for B1I and not another BDS receiver clock.
 */
enum {
    RCB_GPS_C5Q = 0,
    RCB_GAL_C5Q,
    RCB_BDS_C7I,
    RCB_BDS_C5P,
    RCB_BDS_C1P,
    NIFB_PPP
};
#define ND(opt)     ((opt)->nf>=4?NIFB_PPP:((opt)->nf>=3?RCB_BDS_C1P:0))

#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define ID(i,opt)   (NP(opt)+NC(opt)+NT(opt)+NI(opt)+(i))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

  /* PPP receiver-clock slot ----------------------------------------------------*/
static int ppp_clk_index(int sys)
{
    switch (sys) {
    case SYS_GPS:
    case SYS_QZS:
    case SYS_SBS: return 0;
    case SYS_GLO: return 1;
    case SYS_GAL: return 2;
    case SYS_CMP: return 3;
    case SYS_IRN: return 4;
    }
    return -1;
}

/* 按系统+实际信号+存储槽寻找接收机码偏差状态。
 * 当前参考码：GPS/GAL 1C、BDS 2I；参考码偏差与接收机钟差不可分离，
 * 不再单独估计。GPS/GAL 5Q、BDS 7I/5P/1P 分别使用自己的 RCB 状态。
 * 返回 -1 不一定是不支持：也可能是参考码；调用方还需判断 frq==0。
 * BDS B1I 与 B1C 是不同真实频率；本 PPP 的 B1C 放在第 4 存储槽。 */
static int ppp_ifb_index(int sys, int frq, uint8_t code, const prcopt_t* opt)
{
    const char *obs = code2obs(code);

    if (!opt || opt->nf < 3 || !obs || !*obs) return -1;

    if ((sys == SYS_GPS || sys == SYS_QZS) &&
        frq == 2 && !strcmp(obs, "5Q")) return RCB_GPS_C5Q;
    if (sys == SYS_GAL &&
        frq == 2 && !strcmp(obs, "5Q")) return RCB_GAL_C5Q;
    if (sys == SYS_CMP &&
        frq == 1 && !strcmp(obs, "7I")) return RCB_BDS_C7I;
    if (sys == SYS_CMP &&
        frq == 2 && !strcmp(obs, "5P")) return RCB_BDS_C5P;
    if (sys == SYS_CMP && opt->nf >= 4 &&
        frq == 3 && !strcmp(obs, "1P")) return RCB_BDS_C1P;

    return -1; /* reference code or a signal not configured for this data */
}

/* An unmodelled secondary receiver code bias is not measurement noise.
 * Never use its code in the PPP update. Its carrier can still be used with
 * a float ambiguity and is monitored by the phase-quality checks. */
static int ppp_code_modelled(int sys, int frq, uint8_t code,
                             const prcopt_t *opt)
{
    return frq == 0 || ppp_ifb_index(sys, frq, code, opt) >= 0;
}

/* PPP 层支持范围，不等同 Android adapter 的 ppp-safe 输入策略。
 * 前三个槽仍需结合其他模型检查；第 4 槽只允许 BDS 1P。
 * GAL 7Q 默认隔离，可用 -GALE5BPHASE 显式试验；识别信号并不代表使用它。
 * 本次仅添加解释，不调整任何准入规则。 */
static int ppp_signal_supported(int sys, int frq, uint8_t code,
                                const prcopt_t *opt)
{
    const char *obs = code2obs(code);
    const char *p;
    int gal_e5b_phase = 0;

    /* E5b/7Q currently has no receiver-code-bias state, so only its carrier
     * enters PPP. Three independent phone data sets repeatedly reject that
     * phase, while excluding it preserves Galileo E1/E5a and Q6 continuity.
     * Keep an explicit opt-in for testing a future calibrated E5b model. */
    if (opt && (p = strstr(opt->pppopt, "-GALE5BPHASE="))) {
        if (sscanf(p, "-GALE5BPHASE=%d", &gal_e5b_phase) != 1)
            gal_e5b_phase = 0;
    }
    if (!gal_e5b_phase && sys == SYS_GAL && frq == 1 && obs &&
        !strcmp(obs, "7Q")) return 0;
    if (frq < 3) return 1;
    return frq == 3 && sys == SYS_CMP && obs && !strcmp(obs, "1P");
}

/* forward declaration: ambiguity initialization uses the same ionosphere model
   as the measurement equation */
static int model_iono(gtime_t time, const double* pos, const double* azel,
    const prcopt_t* opt, int sat, const double* x,
    const nav_t* nav, double* dion, double* var);

/* VMF3 5x5 grid storage. The row order in the official files is north to
 * south, then west to east. Times in the VMF3 header are UTC; RTKLIB PPP
 * observation times are GPST, hence utc2gpst() at load time. */
typedef struct {
    int valid;
    gtime_t time;
    double ah [VMF3_NGRID];
    double aw [VMF3_NGRID];
    double zhd[VMF3_NGRID];
    double zwd[VMF3_NGRID];
} vmf3_grid_t;

static vmf3_grid_t vmf3_grid[2];
static double vmf3_orog[VMF3_NGRID];
static int vmf3_orog_valid = 0;

/* 读取一份官方 5x5 VMF3 文件：共 36*72 个格点，存 ah/aw 和 ZHD/ZWD。
 * 文件头 ! Epoch 是 UTC，加载时转 GPST；使用文件内容而不是文件名猜时间。
 * 格点纬经度虽然读入，但此实现按官方固定行序存储，不支持任意乱序格点。
 * 文件缺历元、数量错误或无法打开时返回 0，由调用层决定失败或模型回退。 */
static int readvmf3(const char *file, vmf3_grid_t *grid)
{
    FILE *fp;
    char buff[256], tstr[64], *p;
    double ep[6] = {0}, lat, lon, ah, aw, zhd, zwd;
    int n = 0, have_epoch = 0;

    if (!file || !*file || !(fp = fopen(file, "r"))) {
        trace(1, "vmf3: cannot open %s\n", file ? file : "(null)");
        return 0;
    }
    memset(grid, 0, sizeof(*grid));

    while (fgets(buff, sizeof(buff), fp)) {
        if ((p = strstr(buff, "! Epoch:"))) {
            if (sscanf(p + 8, "%lf%lf%lf%lf%lf%lf", ep, ep + 1, ep + 2,
                ep + 3, ep + 4, ep + 5) == 6) have_epoch = 1;
            continue;
        }
        if (buff[0] == '!') continue;
        if (sscanf(buff, "%lf%lf%lf%lf%lf%lf", &lat, &lon, &ah, &aw,
            &zhd, &zwd) != 6) continue;
        if (n >= VMF3_NGRID) {
            fclose(fp);
            trace(1, "vmf3: too many grid records in %s\n", file);
            return 0;
        }
        grid->ah [n] = ah;
        grid->aw [n] = aw;
        grid->zhd[n] = zhd;
        grid->zwd[n] = zwd;
        n++;
    }
    fclose(fp);
    if (!have_epoch || n != VMF3_NGRID) {
        trace(1, "vmf3: invalid file %s (epoch=%d records=%d)\n", file,
            have_epoch, n);
        return 0;
    }
    grid->time = utc2gpst(epoch2time(ep));
    grid->valid = 1;
    time2str(grid->time, tstr, 0);
    trace(2, "vmf3: loaded %s, epoch=%s, records=%d\n", file,
        tstr, n);
    return 1;
}

/* 配套格网高程 orography_ell_5x5：每个格点一个椭球高，不是新一期气象产品。
 * 行序必须与 VMF3 格点相同，用于把格点天顶延迟归算到接收机椭球高。 */
static int readvmf3orog(const char *file)
{
    FILE *fp;
    char buff[128];
    int n = 0;

    if (!file || !*file || !(fp = fopen(file, "r"))) {
        trace(1, "vmf3: cannot open orography file %s\n",
            file ? file : "(null)");
        return 0;
    }
    while (fgets(buff, sizeof(buff), fp)) {
        if (n < VMF3_NGRID && sscanf(buff, "%lf", vmf3_orog + n) == 1) n++;
    }
    fclose(fp);
    vmf3_orog_valid = n == VMF3_NGRID;
    if (!vmf3_orog_valid) {
        trace(1, "vmf3: invalid orography file %s (records=%d)\n", file, n);
    }
    return vmf3_orog_valid;
}

/* 启动时加载两期 VMF3 及格网高程；要求后一期时间严格大于前一期。
 * 存储是文件级 static 全局缓存，不属于每个 rtk_t；同一进程多个独立
 * 解算器若并发使用不同 VMF3 产品，需要另行隔离缓存，不能随意重载。 */
extern int pppvmf3load(const char *file0, const char *file1, const char *orog)
{
    memset(vmf3_grid, 0, sizeof(vmf3_grid));
    vmf3_orog_valid = 0;
    return readvmf3(file0, vmf3_grid) && readvmf3(file1, vmf3_grid + 1) &&
           readvmf3orog(orog) && timediff(vmf3_grid[1].time, vmf3_grid[0].time) > 0.0;
}

/* 两期之间线性时间插值，再取接收机周围四格点做双线性空间插值。
 * 先对四个格点的 ZHD/ZWD 分别做高程归算，再做空间加权。
 * pos 纬经度为弧度，高程为米；经度按 0..360 周期处理。
 * 仅在两期夹住的时间内可用（a 在 0..1），不对预报文件进行时间外推。
 * 返回 0 时 model_trop 等调用方会按其逻辑回退，不代表自动找下一期文件。 */
static int vmf3_grid_interp(gtime_t time, const double *pos, double *ah_out,
                             double *aw_out, double *zhd_out, double *zwd_out)
{
    double lat, lon, glat, glon, u, w, a, dt, den, zhd_t, zwd_t;
    int i0, i1, j0, j1, q, idx[4];
    double ah[4], aw[4], zhd[4], zwd[4], h[4], zhd4[4], zwd4[4];

    if (!ah_out || !aw_out || !zhd_out || !zwd_out ||
        !vmf3_grid[0].valid || !vmf3_grid[1].valid ||
        !vmf3_orog_valid) return 0;
    dt = timediff(vmf3_grid[1].time, vmf3_grid[0].time);
    a = timediff(time, vmf3_grid[0].time) / dt;
    if (a < 0.0 || a > 1.0) return 0; /* do not extrapolate forecast files */

    lat = pos[0] * R2D;
    lon = pos[1] * R2D;
    if (lon < 0.0) lon += 360.0;
    if (lat < -87.5 || lat > 87.5) return 0;

    glat = (87.5 - lat) / 5.0;
    i0 = (int)floor(glat);
    if (i0 < 0) i0 = 0;
    if (i0 >= VMF3_NLAT - 1) i0 = VMF3_NLAT - 2;
    i1 = i0 + 1;
    u = glat - i0;

    glon = (lon - 2.5) / 5.0;
    while (glon < 0.0) glon += VMF3_NLON;
    while (glon >= VMF3_NLON) glon -= VMF3_NLON;
    j0 = (int)floor(glon);
    j1 = (j0 + 1) % VMF3_NLON;
    w = glon - j0;

    idx[0] = i0 * VMF3_NLON + j0;
    idx[1] = i0 * VMF3_NLON + j1;
    idx[2] = i1 * VMF3_NLON + j0;
    idx[3] = i1 * VMF3_NLON + j1;
    for (q = 0; q < 4; q++) {
        ah[q] = vmf3_grid[0].ah[idx[q]] +
                a * (vmf3_grid[1].ah[idx[q]] - vmf3_grid[0].ah[idx[q]]);
        aw[q] = vmf3_grid[0].aw[idx[q]] +
                a * (vmf3_grid[1].aw[idx[q]] - vmf3_grid[0].aw[idx[q]]);
        zhd[q] = vmf3_grid[0].zhd[idx[q]] +
                 a * (vmf3_grid[1].zhd[idx[q]] - vmf3_grid[0].zhd[idx[q]]);
        zwd[q] = vmf3_grid[0].zwd[idx[q]] +
                 a * (vmf3_grid[1].zwd[idx[q]] - vmf3_grid[0].zwd[idx[q]]);
        h[q] = vmf3_orog[idx[q]];

        /* Kouba (2008): grid height -> receiver ellipsoidal height. */
        den = 1.0 - 0.00266 * cos(2.0 * pos[0]) - 0.00000028 * h[q];
        zhd_t = (zhd[q] / 0.0022768) * den;
        zhd_t *= pow(1.0 - 0.0000226 * (pos[2] - h[q]), 5.225);
        zhd_t = 0.0022768 * zhd_t /
                (1.0 - 0.00266 * cos(2.0 * pos[0]) - 0.00000028 * pos[2]);
        zwd_t = zwd[q] * exp(-(pos[2] - h[q]) / 2000.0);
        zhd4[q] = zhd_t;
        zwd4[q] = zwd_t;
    }
    *ah_out = (1.0 - u) * ((1.0 - w) * ah[0] + w * ah[1]) +
              u * ((1.0 - w) * ah[2] + w * ah[3]);
    *aw_out = (1.0 - u) * ((1.0 - w) * aw[0] + w * aw[1]) +
              u * ((1.0 - w) * aw[2] + w * aw[3]);
    *zhd_out = (1.0 - u) * ((1.0 - w) * zhd4[0] + w * zhd4[1]) +
               u * ((1.0 - w) * zhd4[2] + w * zhd4[3]);
    *zwd_out = (1.0 - u) * ((1.0 - w) * zwd4[0] + w * zwd4[1]) +
               u * ((1.0 - w) * zwd4[2] + w * zwd4[3]);
    return *zhd_out > 1.5 && *zhd_out < 3.0 &&
           *zwd_out >= 0.0 && *zwd_out < 1.5;
}

/* 连分式映射：把天顶延迟映射到卫星斜路径；el 为高度角，输出无量纲。
 * 必须由上层保证正高度角，低高度角放大会同时影响延迟和观测噪声。 */
static double vmf3_mapf(double el, double a, double b, double c)
{
    double sinel = sin(el);
    return (1.0 + a / (1.0 + b / (1.0 + c))) /
           (sinel + a / (sinel + b / (sinel + c)));
}

/* VMF3 a-coefficients, ZHD and mapping factors at the receiver ------------
 * VMF3 supplies the ray-traced a-coefficients. The empirical b/c terms are
 * the standard VMF continued-fraction terms; the small hydrostatic height
 * correction follows Niell (1996). */
static int vmf3_trop(gtime_t time, const double *pos, const double *azel,
                     double *mfh, double *mfw, double *zhd, double *zwd)
{
    const double bh = 0.0029, ch = 0.0620;
    const double bw = 0.00146, cw = 0.04391;
    const double aht = 2.53E-5, bht = 5.49E-3, cht = 1.14E-3;
    double ah, aw, dm;

    if (!mfh || !mfw || !zhd || !zwd || azel[1] <= 0.0 ||
        !vmf3_grid_interp(time, pos, &ah, &aw, zhd, zwd)) return 0;

    *mfh = vmf3_mapf(azel[1], ah, bh, ch);
    *mfw = vmf3_mapf(azel[1], aw, bw, cw);
    dm = (1.0 / sin(azel[1]) - vmf3_mapf(azel[1], aht, bht, cht)) *
         MAX(pos[2], 0.0) / 1E3;
    *mfh += dm;
    return 1;
}

/* 返回 ZHD+ZWD 作为对流层状态初值；本函数本身不追加滤波方程。
 * 持续的湿延迟约束由 ppp_res() 中独立的 VMF3 ZWD 分支实现。 */
static int vmf3_ztd_prior(gtime_t time, const double *pos, double *ztd)
{
    double ah, aw, zhd, zwd;
    if (!ztd || !vmf3_grid_interp(time, pos, &ah, &aw, &zhd, &zwd)) return 0;
    *ztd = zhd + zwd;
    return 1;
}

/* -VMF3SIG 是初始化 ZTD 的标准差（米），不是持续 ZWD 约束标准差。
 * 未设置时采用本文件默认值，并按上下限保护；平方后才作为初始方差。 */
static double vmf3_init_sigma(const prcopt_t *opt)
{
    const char *p;
    double sig = VMF3_INIT_SIGMA_DEF;

    if (opt && (p = strstr(opt->pppopt, "-VMF3SIG="))) {
        if (sscanf(p, "-VMF3SIG=%lf", &sig) != 1) sig = VMF3_INIT_SIGMA_DEF;
    }
    return MIN(VMF3_INIT_SIGMA_MAX, MAX(VMF3_INIT_SIGMA_MIN, sig));
}

/* VMF3 wet-delay soft constraint. A zero sigma disables it. The interval is
 * the decorrelation time of the forecast, not an on/off sampling interval.
 * Per-epoch variance is enlarged by interval/dt so that a 1 Hz solution does
 * not count the same slowly varying forecast as an independent measurement
 * hundreds of times. */
static double vmf3_zwd_constraint_sigma(const prcopt_t *opt)
{
    const char *p;
    double sig = 0.0;

    if (opt && (p = strstr(opt->pppopt, "-VMF3ZWDSIG="))) {
        if (sscanf(p, "-VMF3ZWDSIG=%lf", &sig) != 1 || !isfinite(sig)) {
            sig = 0.0;
        }
    }
    if (sig <= 0.0) return 0.0;
    return MIN(2.0, MAX(0.05, sig));
}

/* -VMF3ZWDINT 是重复产品先验的信息分配间隔（秒），不是每隔多少秒
 * 才读文件。默认 300 s；实际每历元方差由 ppp_res() 的时间记账计算。 */
static double vmf3_zwd_constraint_interval(const prcopt_t *opt)
{
    const char *p;
    double sec = 300.0;

    if (opt && (p = strstr(opt->pppopt, "-VMF3ZWDINT="))) {
        if (sscanf(p, "-VMF3ZWDINT=%lf", &sec) != 1 || !isfinite(sec)) {
            sec = 300.0;
        }
    }
    return MIN(3600.0, MAX(1.0, sec));
}

/* 标准差 sqrt(Pii)：固定解读取 Pa，其他解读取 P；内部矩阵存的是方差。
 * SQRT 宏把非正数/NaN 映射为 0 仅用于输出保护，不能据此断言状态精确。 */
static double STD(rtk_t* rtk, int i)
{
    if (rtk->sol.stat == SOLQ_FIX) return SQRT(rtk->Pa[i + i * rtk->na]);
    return SQRT(rtk->P[i + i * rtk->nx]);
}
/* =============================================================================
 * 功能：输出 PPP 解算状态到 .stat 文件中
 * 说明：本函数输出位置、钟差、RCB、大气及可选模糊度状态；逐星残差由其他输出路径负责。
 * $CLK 单位为 ns，顺序为 GPS 钟差及 GLO/GAL/BDS/IRN 的相对 ISB；
 * $RCB/$TROP/$ION/$AMB 单位为米，后面的 STD 是标准差。
 * $POS 是 ECEF 坐标，不是经纬高；$TROP 存总 ZTD，不是纯 ZWD。
 * OUTSTAT_AMB 是编译期开关；无有效解时本函数不输出状态。
 * ============================================================================= */
extern int pppoutstat(rtk_t* rtk, char* buff)
{
    ssat_t* ssat;
    double tow, pos[3], vel[3], acc[3], * x;
    int i, j, week;
    char id[8], * p = buff;

    if (!rtk->sol.stat) return 0;

    trace(3, "pppoutstat:\n");

    tow = time2gpst(rtk->sol.time, &week);

    x = rtk->sol.stat == SOLQ_FIX ? rtk->xa : rtk->x;

    /* receiver position */
    p += sprintf(p, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow,
        rtk->sol.stat, x[0], x[1], x[2], STD(rtk, 0), STD(rtk, 1), STD(rtk, 2));

    /* receiver velocity and acceleration */
    if (rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr, pos);
        ecef2enu(pos, rtk->x + 3, vel);
        ecef2enu(pos, rtk->x + 6, acc);
        p += sprintf(p, "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,"
            "%.4f,%.5f,%.5f,%.5f\n", week, tow, rtk->sol.stat, vel[0], vel[1],
            vel[2], acc[0], acc[1], acc[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    /* receiver clocks */
    i = IC(0, &rtk->opt);
    /* $CLK: GPS clock, then GLO/GAL/BDS/IRN ISBs, followed by their std-devs */
    p += sprintf(p, "$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
        week, tow, rtk->sol.stat, 1,
        x[i] * 1E9 / CLIGHT, x[i + 1] * 1E9 / CLIGHT,
        x[i + 2] * 1E9 / CLIGHT, x[i + 3] * 1E9 / CLIGHT,
        x[i + 4] * 1E9 / CLIGHT,
        STD(rtk, i) * 1E9 / CLIGHT, STD(rtk, i + 1) * 1E9 / CLIGHT,
        STD(rtk, i + 2) * 1E9 / CLIGHT, STD(rtk, i + 3) * 1E9 / CLIGHT,
        STD(rtk, i + 4) * 1E9 / CLIGHT);

    /* receiver code biases (m), relative to each constellation reference */
    if (ND(&rtk->opt) > 0) {
        static const char *name[NIFB_PPP] = {
            "GPS_C5Q", "GAL_C5Q", "BDS_C7I", "BDS_C5P", "BDS_C1P"
        };

        for (i = 0; i < ND(&rtk->opt) && i < NIFB_PPP; i++) {
            int idx = ID(i, &rtk->opt);

            if (idx >= 0 && idx < rtk->nx) {
                p += sprintf(p, "$RCB,%d,%.3f,%d,%s,%.4f,%.4f\n",
                    week, tow, rtk->sol.stat, name[i], x[idx], STD(rtk, idx));
            }
        }
    }

    /* tropospheric parameters */
    if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
        i = IT(&rtk->opt);
        p += sprintf(p, "$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow, rtk->sol.stat,
            1, x[i], STD(rtk, i));
    }
    if (rtk->opt.tropopt == TROPOPT_ESTG) {
        i = IT(&rtk->opt);
        p += sprintf(p, "$TRPG,%d,%.3f,%d,%d,%.5f,%.5f,%.5f,%.5f\n", week, tow,
            rtk->sol.stat, 1, x[i + 1], x[i + 2], STD(rtk, i + 1), STD(rtk, i + 2));
    }
    /* ionosphere parameters */
    if (rtk->opt.ionoopt == IONOOPT_EST) {
        for (i = 0; i < MAXSAT; i++) {
            ssat = rtk->ssat + i;
            if (!ssat->vs) continue;
            j = II(i + 1, &rtk->opt);
            if (rtk->x[j] == 0.0) continue;
            satno2id(i + 1, id);
            p += sprintf(p, "$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n", week, tow,
                rtk->sol.stat, id, rtk->ssat[i].azel[0] * R2D,
                rtk->ssat[i].azel[1] * R2D, x[j], STD(rtk, j));
        }
    }
#ifdef OUTSTAT_AMB
    /* ambiguity parameters */
    int k;
    for (i = 0; i < MAXSAT; i++) for (j = 0; j < NF(&rtk->opt); j++) {
        k = IB(i + 1, j, &rtk->opt);
        if (rtk->x[k] == 0.0) continue;
        satno2id(i + 1, id);
        p += sprintf(p, "$AMB,%d,%.3f,%d,%s,%d,%.4f,%.4f\n", week, tow,
            rtk->sol.stat, id, j + 1, x[k], STD(rtk, k));
    }
#endif
    return (int)(p - buff);
}
/* =============================================================================
 * 功能：剔除进入地球阴影区 (地影区) 的卫星
 * 说明：地影可能使卫星姿态及相位模型不可靠；此实现针对 BLOCK IIA
 *       （天线类型为空时也进入检查），不是把所有星座的地影卫星一律删除。
 *       通过把卫星位置置零使后续几何检查失败；由 posopt[3] 控制是否调用。
 * ============================================================================= */
static void testeclipse(const obsd_t* obs, int n, const nav_t* nav, double* rs)
{
    double rsun[3], esun[3], r, ang, erpv[5] = { 0 }, cosa;
    int i, j;
    const char* type;

    trace(3, "testeclipse:\n");

    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time), erpv, rsun, NULL, NULL);
    normv3(rsun, esun);

    for (i = 0; i < n; i++) {
        type = nav->pcvs[obs[i].sat - 1].type;

        if ((r = norm(rs + i * 6, 3)) <= 0.0) continue;

        /* only block IIA */
        if (*type && !strstr(type, "BLOCK IIA")) continue;

        /* sun-earth-satellite angle */
        cosa = dot3(rs + i * 6, esun) / r;
        cosa = cosa < -1.0 ? -1.0 : (cosa > 1.0 ? 1.0 : cosa);
        ang = acos(cosa);

        /* test eclipse */
        if (ang<PI / 2.0 || r * sin(ang)>RE_WGS84) continue;

        char tstr[40];
        trace(3, "eclipsing sat excluded %s sat=%2d\n", time2str(obs[0].time, tstr, 0),
            obs[i].sat);

        for (j = 0; j < 3; j++) rs[j + i * 6] = 0.0;
    }
}
/* 名义偏航角计算 ------------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta) < 1E-12 && fabs(mu) < 1E-12) return PI;
    return atan2(-tan(beta), sin(mu)) + PI;
}
/* 卫星偏航角计算：当前函数仅使用 yaw_nominal()，没有在这里实现
 * 各型号的完整地影/正午转弯模型；接口保留 sat/type/opt 供后续扩展。 */
extern int yaw_angle(int sat, const char* type, int opt, double beta, double mu,
    double* yaw)
{
    *yaw = yaw_nominal(beta, mu);
    return 1;
}
/* =============================================================================
 * 功能：卫星姿态偏航角模型
 * 说明：用于计算精确的卫星天线相位中心变化以及相位缠绕
 * ============================================================================= */
static int sat_yaw(gtime_t time, int sat, const char* type, int opt,
    const double* rs, double* exs, double* eys)
{
    double rsun[3], ri[6], es[3], esun[3], n[3], p[3], en[3], ep[3], ex[3], E, beta, mu;
    double yaw, cosy, siny, erpv[5] = { 0 };
    int i;

    sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);

    /* beta and orbit angle */
    matcpy(ri, rs, 6, 1);
    ri[3] -= OMGE * ri[1];
    ri[4] += OMGE * ri[0];
    cross3(ri, ri + 3, n);
    cross3(rsun, n, p);
    if (!normv3(rs, es) || !normv3(rsun, esun) || !normv3(n, en) ||
        !normv3(p, ep)) return 0;
    beta = PI / 2.0 - acos(dot3(esun, en));
    E = acos(dot3(es, ep));
    mu = PI / 2.0 + (dot3(es, esun) <= 0 ? -E : E);
    if (mu < -PI / 2.0) mu += 2.0 * PI;
    else if (mu >= PI / 2.0) mu -= 2.0 * PI;

    /* yaw-angle of satellite */
    if (!yaw_angle(sat, type, opt, beta, mu, &yaw)) return 0;

    /* satellite fixed x,y-vector */
    cross3(en, es, ex);
    cosy = cos(yaw);
    siny = sin(yaw);
    for (i = 0; i < 3; i++) {
        exs[i] = -siny * en[i] + cosy * ex[i];
        eys[i] = -cosy * en[i] - siny * ex[i];
    }
    return 1;
}
/* =============================================================================
 * 功能：相位缠绕 (Phase Windup) 误差模型
 * 说明：根据卫星/接收机天线方向计算载波极化造成的相位缠绕。
 *       phw 单位为周，通过与上次 phw 就近接续保持整周连续；不能每秒清零。
 *       此项只改正载波：在 corr_meas 中乘波长后从 L 扣除，不改正伪距。
 * ============================================================================= */
static int model_phw(gtime_t time, int sat, const char* type, int opt,
    const double* rs, const double* rr, double* phw)
{
    double exs[3], eys[3], ek[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
    double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph;
    int i;

    if (opt <= 0) return 1; /* no phase windup */

    /* satellite yaw attitude model */
    if (!sat_yaw(time, sat, type, opt, rs, exs, eys)) return 0;

    /* unit vector satellite to receiver */
    for (i = 0; i < 3; i++) r[i] = rr[i] - rs[i];
    if (!normv3(r, ek)) return 0;

    /* unit vectors of receiver antenna */
    ecef2pos(rr, pos);
    xyz2enu(pos, E);
    exr[0] = E[1]; exr[1] = E[4]; exr[2] = E[7]; /* x = north */
    eyr[0] = -E[0]; eyr[1] = -E[3]; eyr[2] = -E[6]; /* y = west  */

    /* phase windup effect */
    cross3(ek, eys, eks);
    cross3(ek, eyr, ekr);
    for (i = 0; i < 3; i++) {
        ds[i] = exs[i] - ek[i] * dot3(ek, exs) - eks[i];
        dr[i] = exr[i] - ek[i] * dot3(ek, exr) + ekr[i];
    }
    cosp = dot3(ds, dr) / norm(ds, 3) / norm(dr, 3);
    if (cosp < -1.0) cosp = -1.0;
    else if (cosp > 1.0) cosp = 1.0;
    ph = acos(cosp) / 2.0 / PI;
    cross3(ds, dr, drs);
    if (dot3(ek, drs) < 0.0) ph = -ph;

    *phw = ph + floor(*phw - ph + 0.5); /* in cycle */
    return 1;
}
/* =============================================================================
 * 功能：动态计算测量误差方差 (观测噪声 R 阵)
 * 说明：根据高度角和 C/N0 计算伪距、载波相位的观测方差。
 * -WGTELCN=0 保留现有分段 C/N0 模型；=1 使用 Li et al. (2022)
 * 式 (7)-(8)；=2 仅首频使用式 (8)，其余频点保留原模型。
 * ============================================================================= */
/* Storage slots are not calibration identities. Resolve the nominal noise
 * band from constellation+RINEX signal; retain the existing band defaults. */
/* 依据实际信号查标称噪声频带，而非把数组槽直接当作校准频带。
 * 例如 BDS 1P 虽放第 4 存储槽，仍需用其真实频带对应的噪声参数。 */
static int signal_noise_index_ppp(int sys, uint8_t signal, int fallback)
{
    int index = code2idx(sys, signal);
    return index >= 0 && index < NFREQ ? index :
        (fallback >= 0 && fallback < NFREQ ? fallback : 0);
}

/* Optional variance multipliers, phase and code independently. Example:
 * -SIGW=G5Q:2:1;E1C:3:1 means code variance x2/x3, phase unchanged.
 * Only down-weighting (1..100) is permitted. No uncalibrated signal-specific
 * constants or truth-dependent factors are applied by default. */
/* -SIGW 的两项分别是码/相位“方差倍率”，不是标准差倍率；type=0 相位，
 * type=1 伪距。只允许 1..100 的降权；未配置时返回 1，不自动偏爱某星座。 */
static double signal_weight_ppp(int sys, uint8_t signal, int type,
                               const prcopt_t *opt)
{
    const char *p = strstr(opt->pppopt, "-SIGW="), *obs = code2obs(signal);
    char system, sig[3], target = sys == SYS_GPS ? 'G' : sys == SYS_GAL ? 'E' :
        sys == SYS_CMP ? 'C' : sys == SYS_GLO ? 'R' : sys == SYS_QZS ? 'J' :
        sys == SYS_SBS ? 'S' : sys == SYS_IRN ? 'I' : '?';
    double code_factor, phase_factor;
    int consumed;
    if (!p || !obs || !*obs) return 1.0;
    p += 6;
    while (*p && *p != ' ') {
        consumed = 0;
        if (sscanf(p, "%c%2[0-9A-Z]:%lf:%lf%n", &system, sig,
                   &code_factor, &phase_factor, &consumed) != 4 || !consumed)
            break;
        if (p[consumed] && p[consumed] != ';' && p[consumed] != ' ') break;
        if (system == target && !strcmp(sig, obs) &&
            isfinite(code_factor) && isfinite(phase_factor) &&
            code_factor >= 1.0 && code_factor <= 100.0 &&
            phase_factor >= 1.0 && phase_factor <= 100.0)
            return type ? code_factor : phase_factor;
        p += consumed;
        if (*p != ';') break;
        p++;
    }
    return 1.0;
}

/* 返回观测方差（m^2），不是权本身：方差越大，对滤波的影响越弱。
 * 基线：a^2+b^2/sin(el)^2；a/b 由 err[1/2] 和星座系数得到，伪距再
 * 乘 eratio[f] 的标准差倍率。err[0] 不是此处伪距方差的直接来源。
 * 基线再乘分段 C/N0 方差系数；有效方差下限为码 1 m^2、相位 0.01^2 m^2。
 * -WGTELCN=1/2 才走可选联合模型，且受频带/CN0 条件限制；该分支
 * 成功后直接返回，不再叠加后面的基线 C/N0 分段倍率。
 * 轨道、大气、平滑、持续异常和抗差方差在 ppp_res 中另行叠加。
 * 拟合 C/N0 系数具有设备依赖性，本次不修改参数或启用试验模型。 */
static double varerr(unsigned char sat, int sys, double el, double snr_dbhz,
    int freq, int type, const prcopt_t* opt)
{
    double fact = 1.0, eratio = 1.0, a, b, var, sinel, snr_weight = 1.0;
    double cna = 35.0833, cnb = 0.1365, cnc = -0.0005;
    double eldeg, cnfit, cnmax, cnmin, cnvar, cnfitvar, scale, var0;
    const char* p;
    int combined = 0;
    int f = freq >= 0 && freq < NFREQ ? freq : 0;

    (void)sat;

    switch (sys) {
    case SYS_GLO: fact = EFACT_GLO; break;
    case SYS_SBS: fact = EFACT_SBS; break;
    case SYS_GPS: fact = EFACT_GPS; break;
    case SYS_QZS: fact = EFACT_QZS; break;
    case SYS_CMP: fact = EFACT_CMP; break;
    case SYS_GAL: fact = EFACT_GAL; break;
    case SYS_IRN: fact = EFACT_IRN; break;
    default:      fact = 1.0;       break;
    }

    /* RTKLIB semantics: err[1]/err[2] are phase base/elevation terms;
       pseudorange is scaled by eratio, not by err[0]. */
    if (type) {
        eratio = opt->eratio[f] > 0.0 ? opt->eratio[f] : opt->eratio[0];
        if (eratio <= 0.0) eratio = 100.0;
        fact *= eratio;
    }
    if (opt->ionoopt == IONOOPT_IFLC) fact *= 3.0;

    sinel = sin(MAX(el, 5.0 * D2R));
    a = fact * opt->err[1];
    b = fact * opt->err[2];
    var = a * a + b * b / (sinel * sinel);

    /* Li et al., Sensors 2022, 22, 2804, Eq. (7)-(8). The polynomial is
     * device-specific (Xiaomi MI8); allow calibration without recompiling.
     * Use the existing phase/code zenith variance and constellation factor
     * as sigma0^2 so code and phase stay on their original noise scales.
     * f=0 covers L1/G1/B1/E1; f=2 is L5/E5a only for GPS/Galileo.
     * Other bands and missing C/N0 retain the tested baseline model. */
    if ((p = strstr(opt->pppopt, "-WGTELCN="))) {
        if (sscanf(p, "-WGTELCN=%d", &combined) != 1) combined = 0;
    }
    if ((combined == 1 || combined == 2) && snr_dbhz > 0.0 &&
        snr_dbhz < 100.0 &&
        (f == 0 || (combined == 1 && f == 2 &&
            (sys == SYS_GPS || sys == SYS_GAL)))) {
        if ((p = strstr(opt->pppopt, "-WGTCNA="))) {
            if (sscanf(p, "-WGTCNA=%lf", &cna) != 1) cna = 35.0833;
        }
        if ((p = strstr(opt->pppopt, "-WGTCNB="))) {
            if (sscanf(p, "-WGTCNB=%lf", &cnb) != 1) cnb = 0.1365;
        }
        if ((p = strstr(opt->pppopt, "-WGTCNC="))) {
            if (sscanf(p, "-WGTCNC=%lf", &cnc) != 1) cnc = -0.0005;
        }
        if (fabs(cna) < 100.0 && fabs(cnb) < 10.0 && fabs(cnc) < 1.0) {
            eldeg = MIN(90.0, MAX(10.0, el * R2D));
            sinel = sin(eldeg * D2R);
            var0 = a * a + b * b;
            cnmax = f == 0 ? 45.0 : 40.0;
            cnmin = f == 0 ? 25.0 : 20.0;
            cnvar = var0 * pow(10.0, MAX(cnmax - snr_dbhz, 0.0) / 10.0);
            scale = (1.0 / SQR(sin(10.0 * D2R)) - 1.0) /
                (pow(10.0, (cnmax - cnmin) / 10.0) - 1.0);
            if (f == 0) {
                cnfit = cna + cnb * eldeg + cnc * eldeg * eldeg;
                cnfitvar = var0 * pow(10.0, MAX(cnmax - cnfit, 0.0) / 10.0);
                var = var0 / (sinel * sinel) + fabs(cnvar - cnfitvar) * scale;
            }
            else {
                var = cnvar * scale;
            }
            if (isfinite(var)) {
                return MAX(var, type ? SQR(1.0) : SQR(0.01));
            }
        }
    }

    /* Smartphone C/N0 weighting. This is deliberately moderate: bad data are
       down-weighted, but a hard variance cap never makes a bad observation
       look better than it really is. */
    if (snr_dbhz > 0.0) {
        if (snr_dbhz >= 40.0) snr_weight = 1.0;
        else if (snr_dbhz >= 35.0) snr_weight = 1.0 + (40.0 - snr_dbhz) * 0.2; /* 1..2 */
        else if (snr_dbhz >= 30.0) snr_weight = 2.0 + (35.0 - snr_dbhz) * 0.4; /* 2..4 */
        else if (snr_dbhz >= 25.0) snr_weight = 4.0 + (30.0 - snr_dbhz) * 0.8; /* 4..8 */
        else                     snr_weight = 12.0;
    }
    else snr_weight = 4.0;

    return MAX(var * snr_weight, type ? SQR(1.0) : SQR(0.01));
}
/* =============================================================================
 * 功能：初始化一个状态参数及其协方差，不是初始化整个滤波器。
 * 赋 x[i]=xi，清除 P 第 i 行/列与其他状态的旧相关，再设 Pii=var。
 * 周跳、换信号等场景必须同时清相关，不能只覆盖模糊度数值。
 * 0 在当前 RTKLIB 滤波实现中还用于不活跃状态；1E-6 常用来激活近零状态。
 * ============================================================================= */
static inline void initx(rtk_t* rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i] = xi;
    for (j = 0; j < rtk->nx; j++) rtk->P[i + j * rtk->nx] = 0.0;
    for (j = 0; j < rtk->nx; j++) rtk->P[j + i * rtk->nx] = 0.0;
    rtk->P[i + i * rtk->nx] = var;
}
/* =============================================================================
 * 功能：计算几何无关组合 (Geometry-Free, GF)
 * 说明：计算 lambda1*L1-lambdaf*Lf，输出米；几何/钟差被相消，
 *       但电离层和模糊度仍在，因此组合跳变不是纯伪距误差。
 *       f 必须是次频槽，任一相位/频率缺失返回 0，表示本次不可检测。
 * ============================================================================= */
static double gfmeas(const obsd_t* obs, const nav_t* nav, int f)
{
    double freq1, freq2;

    if (f <= 0 || f >= NFREQ) return 0.0;
    freq1 = sat2freq(obs->sat, obs->code[0], nav);
    freq2 = sat2freq(obs->sat, obs->code[f], nav);
    if (freq1 == 0.0 || freq2 == 0.0 || obs->L[0] == 0.0 || obs->L[f] == 0.0) return 0.0;
    return (obs->L[0] / freq1 - obs->L[f] / freq2) * CLIGHT;
}
/* =============================================================================
 * 功能：计算 MW 组合 (Melbourne-Wubbena)
 * 说明：宽巷相位减窄巷伪距，输出米，需两频均有 P/L。
 *       手机伪距噪声也进入 MW；MW 跳变不能无条件归因于两个载波同时周跳。
 * ============================================================================= */
static double mwmeas(const obsd_t* obs, const nav_t* nav, int f)
{
    double freq1, freq2;

    if (f <= 0 || f >= NFREQ) return 0.0;
    freq1 = sat2freq(obs->sat, obs->code[0], nav);
    freq2 = sat2freq(obs->sat, obs->code[f], nav);

    if (freq1 == 0.0 || freq2 == 0.0 || obs->L[0] == 0.0 || obs->L[f] == 0.0 ||
        obs->P[0] == 0.0 || obs->P[f] == 0.0) return 0.0;
    return (obs->L[0] - obs->L[f]) * CLIGHT / (freq1 - freq2) -
        (freq1 * obs->P[0] + freq2 * obs->P[f]) / (freq1 + freq2);
}
/* =============================================================================
 * 功能：对伪距和载波相位进行各类修正
 * 说明：将输入相位“周”转“米”，减天线项与相位缠绕；码减天线项。
 *       P/L 独立判断缺测，不因某个码缺失而一起丢掉同频载波。
 *       精密产品路径按当前 sat+code+time 查码 OSB：P_corrected=P_raw-OSB。
 *       OSB 已知但过期时禁用该码；从未提供时按现有规则尝试 legacy DCB。
 *       BDS 第 4 槽 C1P 无 legacy DCB 回退，没有时效内 OSB 时仅保留载波。
 *       本函数没有对 L 扣除相位 OSB；不要把“读取 BIA”当作已经实现 PPP-AR。
 *       此处不扣除电离层/对流层/钟差：它们统一进入后面的观测方程。
 *       Lc/Pc 为可选消电离层组合；非组合 PPP 仍使用逐频的 L[f]/P[f]。
 * ============================================================================= */
static void corr_meas(const obsd_t* obs, const nav_t* nav, const double* azel,
    const prcopt_t* opt, const double* dantr,
    const double* dants, double phw, double* L, double* P,
    double* Lc, double* Pc)
{
    double freq[NFREQ] = { 0 }, C1, C2, osb;
    int i, ix = 0, frq, frq2 = -1, bias_ix, osb_status;
    int sys = satsys(obs->sat, NULL);

    for (i = 0; i < opt->nf && i < NFREQ; i++) {
        L[i] = P[i] = 0.0;
        freq[i] = sat2freq(obs->sat, obs->code[i], nav);
        if (freq[i] == 0.0) continue;
        if (testsnr(0, i, azel[1], obs->SNR[i] * SNR_UNIT, &opt->snrmask)) continue;

        /* Phase and code can be intermittently missing on phones. Keep each
           observable independently instead of discarding both. */
        if (obs->L[i] != 0.0) {
            L[i] = obs->L[i] * CLIGHT / freq[i] - dants[i] - dantr[i]
                - phw * CLIGHT / freq[i];
        }
        if (obs->P[i] != 0.0) {
            P[i] = obs->P[i] - dants[i] - dantr[i];
        }
        if (P[i] == 0.0) continue;

        if (opt->sateph == EPHOPT_SSRAPC || opt->sateph == EPHOPT_SSRCOM) {
            if (sys == SYS_GPS)      ix = (i == 0 ? CODE_L1W - 1 : CODE_L2W - 1);
            else if (sys == SYS_GLO) ix = (i == 0 ? CODE_L1P - 1 : CODE_L2P - 1);
            else if (sys == SYS_GAL) ix = (i == 0 ? CODE_L1X - 1 : CODE_L7X - 1);
            else continue;
            if (obs->code[i] > 0) {
                P[i] += nav->ssr[obs->sat - 1].cbias[obs->code[i] - 1]
                    - nav->ssr[obs->sat - 1].cbias[ix];
            }
        }
        else {
            /* Prefer absolute Bias-SINEX code OSB when available.
             * Convention: corrected code = observed code - OSB.
             * Fall back to the legacy relative DCB table otherwise. */
            osb_status = codeosb_at(nav, obs->time, obs->sat,
                                    obs->code[i], &osb);
            if (osb_status > 0) {
                P[i] -= osb;
            }
            else if (osb_status < 0) {
                /* A known but expired OSB must not be silently reused or
                 * replaced with an uncalibrated code observation. */
                P[i] = 0.0;
            }
            else {
                frq = i;
                /* There is no legacy B1C DCB slot. Without a current C1P
                 * OSB, retain its carrier but never use uncalibrated code. */
                if (sys == SYS_CMP && frq == 3) {
                    P[i] = 0.0;
                    continue;
                }
                if (frq >= MAX_CODE_BIAS_FREQS) continue;
                bias_ix = code2bias_ix(sys, obs->code[i]);
                if (bias_ix > 0 && bias_ix <= MAX_CODE_BIASES) {
                    P[i] += nav->cbias[obs->sat - 1][frq][bias_ix - 1];
                }
            }
        }
    }
    for (; i < opt->nf; i++) { /* defensive for builds with opt->nf > NFREQ */
        /* caller arrays are NFREQ wide; nothing can be represented here */
    }

    *Lc = *Pc = 0.0;
    if (opt->nf < 2 || freq[0] == 0.0) return;

    /* Prefer L2 if present; otherwise use the native L5/E5 slot (index 2). */
    if (freq[1] != 0.0 && (L[1] != 0.0 || P[1] != 0.0)) frq2 = 1;
    else if (opt->nf >= 3 && NFREQ >= 3 && freq[2] != 0.0 && (L[2] != 0.0 || P[2] != 0.0)) frq2 = 2;
    if (frq2 < 0 || fabs(SQR(freq[0]) - SQR(freq[frq2])) < 1.0) return;

    C1 = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[frq2]));
    C2 = -SQR(freq[frq2]) / (SQR(freq[0]) - SQR(freq[frq2]));
    if (L[0] != 0.0 && L[frq2] != 0.0) *Lc = C1 * L[0] + C2 * L[frq2];
    if (P[0] != 0.0 && P[frq2] != 0.0) *Pc = C1 * P[0] + C2 * P[frq2];
}

/* causal pseudorange smoothing by Doppler ----------------------------------
 * Mirrors the phone profile in raPPPid: P(k) is blended with the previous
 * smoothed range propagated by integrated Doppler. It is deliberately
 * causal (current and previous epoch only), so it is usable in real time.
 * Enable with -DOPPSM=0.90; omit the option to retain the original RTKLIB
 * measurements exactly. */
/* -DOPPSM 是上一历元预测伪距的保留权重，0 关闭，允许范围 0..0.995。
 * 例如 0.90 表示标称 1 s 时旧预测占 90%，当前原始伪距占 10%。 */
static double doppsm_factor(const prcopt_t *opt)
{
    const char *p;
    double factor = 0.0;

    if (opt && (p = strstr(opt->pppopt, "-DOPPSM="))) {
        if (sscanf(p, "-DOPPSM=%lf", &factor) != 1 || !isfinite(factor)) factor = 0.0;
    }
    return MIN(0.995, MAX(0.0, factor));
}

/* duration of the code+Doppler bootstrap, in seconds ----------------------*/
static double doppsm_warmup_time(const prcopt_t *opt)
{
    const char *p;
    double sec = 0.0;

    if (opt && (p = strstr(opt->pppopt, "-DOPPWARM="))) {
        if (sscanf(p, "-DOPPWARM=%lf", &sec) != 1) sec = 0.0;
    }
    return MIN(300.0, MAX(0.0, sec));
}

/* read a bounded numeric switch from pppopt -------------------------------*/
static double pppopt_number(const prcopt_t *opt, const char *key,
                            double def, double minval, double maxval)
{
    const char *p;
    double value = def;

    if (opt && key && (p = strstr(opt->pppopt, key))) {
        if (sscanf(p + strlen(key), "%lf", &value) != 1 || !isfinite(value)) {
            value = def;
        }
    }
    return MIN(maxval, MAX(minval, value));
}

/* Statistics revision bits: 1=smoother covariance/time, 2=prior information
 * bookkeeping, 4=experimental causal correlation inflation (requires bit 1).
 * Zero is the exact legacy numerical path. No truth coordinates are used. */
/* -STATMOD 位掩码，默认 3=位1+位2：
 * 位1（数值1）：按真实 dt 调整平滑权重并传播热噪声协方差；
 * 位2（数值2）：记录产品先验的成功使用时间，避免初始化/重复历元双计信息；
 * 位3（数值4）：试验性平滑时间相关方差膨胀。0 保留旧数值路径。
 * 这里是近似统计处理，不是完整的有色噪声滤波或产品误差相关模型。 */
static int statistical_model_ppp(const prcopt_t *opt)
{
    return (int)pppopt_number(opt, "-STATMOD=", 3.0, 0.0, 7.0);
}

/* Fraction of one independent prior's information. Initialization already
 * consumes a prior. Duplicate/reversed epochs get no additional information;
 * outages never grant more than one full prior. This is an information-budget
 * approximation, NOT a fitted product-error correlation time. */
/* information=min(1,dt/interval)：本次先验最多占一份独立先验的信息量。
 * 已有记录且重复/倒序历元返回 0；无记录返回 1；长中断最多仍为 1。
 * 调用方用 R_base/information 扩大方差，而不是把残差强制改成更小。
 * 该间隔是当前信息预算的控制参数，不声称由此测得了产品真实相关时间。 */
static double atmosphere_information_ppp(gtime_t time, gtime_t last, double interval)
{
    double dt;
    if (!last.time) return 1.0;
    dt = timediff(time, last);
    if (!isfinite(dt) || dt <= 0.0 || !isfinite(interval) || interval <= 0.0) return 0.0;
    return MIN(1.0, dt / interval);
}

/* Shared trapezoidal Doppler sample induces a cross term. Rate noise is m/s,
 * not Hz. White input-noise assumptions describe only receiver thermal noise;
 * atmospheric/product/multipath errors are NOT reduced by this propagation. */
/* 梯形积分的相邻历元共享一个 Doppler 样本，不能当作完全独立噪声。
 * variance 保存平滑伪距热噪声方差，ratecov 为与当前距离变化率的协方差，
 * pastcov 记录对过去平滑值的累计相关项；rawvar 为米平方，ratevar 为
 * (m/s)^2。这里只传播白输入噪声，不代表多路径/大气误差也被平均消除。 */
static void smoother_covariance_ppp(double alpha, double dt, double rawvar,
    double ratevar, double *variance, double *old_ratevar, double *ratecov,
    double *pastcov)
{
    double b = dt * 0.5, v = *variance, c = *ratecov;
    *variance = SQR(alpha) * (v + SQR(b) * (*old_ratevar + ratevar) - 2.0 * b * c)
              + SQR(1.0 - alpha) * rawvar;
    *pastcov = alpha * (*pastcov + v - b * c);
    *ratecov = -alpha * b * ratevar;
    *old_ratevar = ratevar;
}

/* zero ambiguity states at the code-to-phase handover ---------------------
 * udbias_ppp() will initialize them from the already-smoothed code at the
 * following epoch. Clearing covariance rows as well is essential: otherwise
 * dormant ambiguities from the bootstrap remain correlated with position. */
static void reset_ambiguities_ppp(rtk_t *rtk)
{
    int i, f, sat, j;

    for (f = 0; f < NF(&rtk->opt); f++) for (sat = 1; sat <= MAXSAT; sat++) {
        j = IB(sat, f, &rtk->opt);
        rtk->x[j] = 0.0;
        for (i = 0; i < rtk->nx; i++) {
            rtk->P[i + j * rtk->nx] = 0.0;
            rtk->P[j + i * rtk->nx] = 0.0;
        }
    }
}

/* select the observation stage for this epoch -----------------------------*/
/* -DOPPWARM>0 且开启平滑时，启动阶段暂不把相位行送入测量更新。
 * 用观测 GPST 计算持续时间，不用电脑墙钟；切入相位时重置旧模糊度。
 * warmup=0 不启用此阶段；Q5 启动不能仅凭该函数判断，还需看实际验后状态。 */
static int update_doppsm_warmup(rtk_t *rtk, gtime_t time)
{
    double warmup = doppsm_warmup_time(&rtk->opt), elapsed;

    rtk->ppp_code_warmup = 0;
    if (warmup <= 0.0 || doppsm_factor(&rtk->opt) <= 0.0) return 0;

    if (!rtk->ppp_start_valid) {
        rtk->ppp_start_time = time;
        rtk->ppp_start_valid = 1;
    }
    elapsed = timediff(time, rtk->ppp_start_time);
    if (elapsed < warmup) {
        rtk->ppp_code_warmup = 1;
        return 1;
    }
    if (!rtk->ppp_phase_started) {
        reset_ambiguities_ppp(rtk);
        rtk->ppp_phase_started = 1;
        trace(3, "$DOPP_HANDOVER,elapsed=%.1f,phase=enabled\n", elapsed);
    }
    return 0;
}

/* update one causal code smoother per satellite/frequency once per epoch ---*/
/* 每星每频槽的因果平滑，一历元执行一次：
 *   P_pred=P_smooth_old-lambda*(D_old+D_now)/2*dt
 *   P_smooth=(1-alpha)*P_raw+alpha*P_pred
 * STATMOD 位1启用时 alpha=factor^dt，适应非 1 Hz 的采样间隔。
 * 缺码/缺 Doppler、周跳、预检坏码使缓存无效；首次有效观测或异常时间间隔
 * 用原始伪距重建缓存，不把历史预测跨长中断外推。
 * 平滑缓存是原始码域；产品/天线改正保留到 corr_meas 之后再应用。
 * obsd_t 没有 Android Doppler uncertainty，此处使用 err[4] 的 Hz 标准差
 * （无效时回退 1 Hz），转米/秒后平方，不直接沿用 ADR uncertainty。 */
static int uddoppsm_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
    const double *rs)
{
    double factor = doppsm_factor(&rtk->opt), freq, dt, pred;
    double alpha, rawvar, ratevar, pos[3], e[3], azel[2], el, dop_sigma;
    double sum_raw=0.0,sum_sm=0.0,sum_rate=0.0,max_corr_ratio=1.0;
    int model = statistical_model_ppp(&rtk->opt);
    int i, f, sat, nvalid = 0;

    if (factor <= 0.0) return 0;

    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        if (sat < 1 || sat > MAXSAT) continue;

        for (f = 0; f < rtk->opt.nf && f < NFREQ; f++) {
            ssat_t *ssat = rtk->ssat + sat - 1;

            if (!ppp_signal_supported(satsys(sat, NULL), f,
                                      obs[i].code[f], &rtk->opt)) continue;

            freq = sat2freq(sat, obs[i].code[f], nav);
            if (obs[i].P[f] == 0.0 || !isfinite(obs[i].P[f]) || freq == 0.0 ||
                obs[i].D[f] == 0.0 || !isfinite(obs[i].D[f]) ||
                fabs(obs[i].D[f]) > 20000.0 ||
                (obs[i].LLI[f] & LLI_SLIP) ||
                (ssat->slip[f] & LLI_SLIP) || ssat->ppp_code_bad[f]) {
                ssat->psmvalid[f] = 0;
                continue;
            }
            dt = ssat->psmvalid[f] ? timediff(obs[i].time, ssat->psmt[f]) : 0.0;
            if ((model & 1) && ssat->psmvalid[f] && dt == 0.0) continue;
            el = ssat->azel[1];
            if (rs && norm(rtk->sol.rr, 3) > 1E6 && geodist(rs+i*6,rtk->sol.rr,e)>0.0) {
                ecef2pos(rtk->sol.rr,pos);satazel(pos,e,azel);el=azel[1];
            }
            rawvar = varerr(sat,satsys(sat,NULL),MAX(el,rtk->opt.elmin),
                obs[i].SNR[f]*SNR_UNIT,
                signal_noise_index_ppp(satsys(sat,NULL),obs[i].code[f],f),1,&rtk->opt) *
                signal_weight_ppp(satsys(sat,NULL),obs[i].code[f],1,&rtk->opt);
            /* obsd_t has no Android rate uncertainty: use configured Hz sigma,
               with 1 Hz fallback, and expose it in diagnostics. */
            dop_sigma=rtk->opt.err[4];
            if (!isfinite(dop_sigma) || dop_sigma<=0.0) dop_sigma=1.0;
            ratevar=SQR(CLIGHT/freq*dop_sigma);
            if (!ssat->psmvalid[f] || ((model&1)?dt<0.0:dt<0.2) || dt > 5.0) {
                ssat->psmP[f] = obs[i].P[f];
                ssat->psmD[f] = obs[i].D[f];
                ssat->psmt[f] = obs[i].time;
                ssat->psmvalid[f] = 1;
                ssat->psm_var[f]=rawvar;ssat->psm_rate_var[f]=ratevar;
                ssat->psm_rate_cov[f]=ssat->psm_past_cov[f]=0.0;
                continue;
            }

            /* RINEX Doppler is cycles/s; range-rate is -lambda*D. */
            pred = ssat->psmP[f] - CLIGHT / freq *
                   0.5 * (ssat->psmD[f] + obs[i].D[f]) * dt;
            alpha=(model & 1)?pow(factor,dt):factor; /* factor is retained weight at 1 s */
            if(model&1)ssat->psmP[f] = (1.0 - alpha) * obs[i].P[f] + alpha * pred;
            else ssat->psmP[f] = (1.0 - factor) * obs[i].P[f] + factor * pred;
            smoother_covariance_ppp(alpha,dt,rawvar,ratevar,&ssat->psm_var[f],
                &ssat->psm_rate_var[f],&ssat->psm_rate_cov[f],&ssat->psm_past_cov[f]);
            trace(4,"$SMOOTH_STAT,sat=%d,sig=%s,dt=%.6f,alpha=%.6f,rawvar=%.6f,var=%.6f,ratevar=%.6f,ratecov=%.6f,pastcov=%.6f,mode=%d\n",
                sat,code2obs(obs[i].code[f]),dt,alpha,rawvar,ssat->psm_var[f],
                ratevar,ssat->psm_rate_cov[f],ssat->psm_past_cov[f],model);
            ssat->psmD[f] = obs[i].D[f];
            ssat->psmt[f] = obs[i].time;
            sum_raw+=rawvar;sum_sm+=ssat->psm_var[f];sum_rate+=ratevar;
            max_corr_ratio=MAX(max_corr_ratio,
                (ssat->psm_var[f]+2.0*ssat->psm_past_cov[f])/MAX(rawvar,1E-12));
            nvalid++;
        }
    }
    if(nvalid>0)trace(3,"$SMOOTH_EPOCH_STAT,mode=%d,n=%d,raw_var_mean=%.6f,thermal_var_mean=%.6f,rate_var_mean=%.6f,correlation_ratio_max=%.6f,correlation_inflation=%s\n",
        model,nvalid,sum_raw/nvalid,sum_sm/nvalid,sum_rate/nvalid,max_corr_ratio,
        (model&4)?"EXPERIMENTAL":"OFF_RAW_FLOOR");
    return nvalid;
}

/* 改正后的码 P 加上“平滑原始码-原始码”的差值。
 * 这样不重复做 OSB/天线改正，也不在平滑缓存内混入跨时变化的产品改正。
 * 已被 corr_meas 禁用的码 P==0 不因平滑缓存存在而被复活。 */
static int applydoppsm_ppp(const rtk_t *rtk, const obsd_t *obs, double *P)
{
    int f, nvalid = 0, sat = obs->sat;
    const ssat_t *ssat;

    if (doppsm_factor(&rtk->opt) <= 0.0 || sat < 1 || sat > MAXSAT) return 0;
    ssat = rtk->ssat + sat - 1;
    for (f = 0; f < rtk->opt.nf && f < NFREQ; f++) {
        if (!ssat->psmvalid[f] || ssat->ppp_code_bad[f] ||
            obs->P[f] == 0.0 || P[f] == 0.0) continue;
        P[f] += ssat->psmP[f] - obs->P[f];
        nvalid++;
    }
    return nvalid;
}

/* Keep a bad carrier signal out of the filter without discarding its code or
 * other frequencies. The quarantine is disabled unless -PPPQUAR=N,SEC is set. */
/* 查询该星该槽载波是否仍在隔离时间内；隔离不自动删除同频伪距。
 * 隔离历史绑定实际 signal，换信号时由 prepare_signal_quality_ppp 清理。 */
static int ppp_phase_quarantined(const ssat_t *ssat, int f, gtime_t time)
{
    return f >= 0 && f < NFREQ && ssat->ppp_phase_block_until[f].time &&
        timediff(time, ssat->ppp_phase_block_until[f]) < 0.0;
}

/* A frequency slot is not a signal identity: RINEX selection can change its
 * code within a session. Never carry an ambiguity, Doppler smoother or
 * quarantine history from one actual signal into another. */
/* 先识别存储槽内的实际 code 是否变化，再开始本历元预处理。
 * 换信号必须清掉旧模糊度、GF/MW、相位/Doppler 历史、平滑、隔离及质量
 * 状态；否则两个硬件/频率不同的信号会被错误拼成一条连续相位弧。
 * 同一信号因较长中断/禁用 SIGQC 时还会清除持续异常判别的历史。
 * 这里只管理历史和输出噪声配置日志，不重新挑选 RINEX/Android 输入信号。 */
static void prepare_signal_quality_ppp(rtk_t *rtk, const obsd_t *obs, int n)
{
    int i, f, sat, j, k;
    char sid[8], str[40];

    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        if (sat < 1 || sat > MAXSAT) continue;
        for (f = 0; f < rtk->opt.nf && f < NFREQ; f++) {
            ssat_t *ss = rtk->ssat + sat - 1;
            uint8_t old_code = ss->ppp_track_code[f];
            uint8_t new_code = obs[i].code[f];

            /* Sparse/unmodelled observations must not retain a permanent
             * penalty. Never share history with another actual signal. */
            for (k = 0; k < 2; k++) {
                ppp_signal_quality_t *q = &ss->ppp_res_quality[f][k];
                double gap = q->time.time ? timediff(obs[i].time, q->time) : 0.0;
                if (pppopt_number(&rtk->opt,"-SIGQC=",1.0,0.0,1.0)<0.5 ||
                    (q->time.time && (gap < 0.0 || gap > 10.0))) {
                    if (q->weak) {
                        satno2id(sat,sid);time2str(obs[i].time,str,2);
                        trace(2,"$SIG_QUALITY_RESET,%s,sat=%s,sig=%s,type=%s,reason=GAP_OR_DISABLED\n",
                              str,sid,code2obs(old_code),k?"CODE":"PHASE");
                    }
                    memset(q,0,sizeof(*q));
                }
            }

            if (new_code == CODE_NONE || new_code == old_code) continue;
            satno2id(sat,sid);time2str(obs[i].time,str,2);
            trace(2,"$SIG_NOISE_PROFILE,%s,sat=%s,sig=%s,F%d,noise_band=%d,code_var=%.3f,phase_var=%.3f,quality=%d,phase_auto=%d\n",
                str,sid,code2obs(new_code),f+1,
                signal_noise_index_ppp(satsys(sat,NULL),new_code,f)+1,
                signal_weight_ppp(satsys(sat,NULL),new_code,1,&rtk->opt),
                signal_weight_ppp(satsys(sat,NULL),new_code,0,&rtk->opt),
                pppopt_number(&rtk->opt,"-SIGQC=",1.0,0.0,1.0)>=0.5,
                pppopt_number(&rtk->opt,"-SIGQCPH=",0.0,0.0,1.0)>=0.5);
            if (old_code != CODE_NONE) {
                j = IB(sat, f, &rtk->opt);
                if (j >= 0 && j < rtk->nx) initx(rtk, 0.0, 0.0, j);
                ss->outc[f] = 0;
                ss->slip[f] = 0;
                ss->psmvalid[f] = 0;
                ss->ppp_dop_valid[f] = 0;
                ss->ppp_half_seen[f] = 0;
                ss->ppp_phase_reject_streak[f] = 0;
                ss->ppp_phase_quar_count[f] = 0;
                ss->ppp_phase_block_until[f].time = 0;
                ss->ppp_phase_reject_time[f].time = 0;
                ss->ppp_phase_good_since[f].time = 0;
                ss->ppp_diag_cmc_time[f].time = 0;
                ss->ppp_diag_cmc_code[f] = CODE_NONE;
                ss->ppp_cmc_quality_time[f].time = 0;
                ss->ppp_cmc_quality_count[f] = 0;
                ss->ppp_cmc_code_weak[f] = 0;
                memset(ss->ppp_res_quality[f],0,sizeof(ss->ppp_res_quality[f]));
                ss->pt[0][f].time = ss->pt[1][f].time = 0;
                ss->ph[0][f] = ss->ph[1][f] = 0.0;
                if (f == 0) {
                    for (k = 0; k < NFREQ - 1; k++) {
                        ss->gf[k] = ss->mw[k] = 0.0;
                    }
                }
                else {
                    ss->gf[f - 1] = ss->mw[f - 1] = 0.0;
                }
                satno2id(sat, sid);
                time2str(obs[i].time, str, 2);
                trace(2, "$SIG_SWITCH,%s,sat=%s,F%d,old=%s,new=%s,reset=AMB+HISTORY\n",
                      str, sid, f + 1, code2obs(old_code), code2obs(new_code));
            }
            ss->ppp_track_code[f] = new_code;
        }
    }
}
/* =============================================================================
 * 功能：利用 LLI (Loss of Lock Indicator) 标志位探测周跳
 * 说明：读取 LLI_SLIP，同时记录 LLI_HALFC 的有效状态。
 *       半周不确定状态发生变化时也标记周跳；持续半周不确定不每秒重复重置。
 *       本函数只设置 ssat.slip，模糊度的真正重建在 udbias_ppp 中执行。
 *       有 RESET/SLIP 的当前相位可以保留，不能把标记周跳等同删除相位值。
 * ============================================================================= */
static void detslp_ll(rtk_t* rtk, const obsd_t* obs, int n)
{
    int i, j, nf = rtk->opt.nf;

    trace(3, "detslp_ll: n=%d\n", n);

    for (i = 0; i < n && i < MAXOBS; i++) for (j = 0; j < nf; j++) {
        ssat_t *ss = &rtk->ssat[obs[i].sat - 1];
        int half_invalid, slip;
        if (!ppp_signal_supported(satsys(obs[i].sat, NULL), j,
                                  obs[i].code[j], &rtk->opt)) continue;
        if (obs[i].L[j] == 0.0) continue;

        half_invalid = (obs[i].LLI[j] & LLI_HALFC) != 0;
        slip = (obs[i].LLI[j] & LLI_SLIP) != 0 ||
            (ss->ppp_half_seen[j] &&
             ss->ppp_half_invalid[j] != half_invalid);
        ss->ppp_half_seen[j] = 1;
        ss->ppp_half_invalid[j] = (uint8_t)half_invalid;
        if (!slip) continue;

        trace(3, "detslp_ll: slip detected sat=%2d f=%d lli=%d half=%d\n",
            obs[i].sat, j + 1, obs[i].LLI[j], half_invalid);

        ss->slip[j] |= LLI_SLIP;
    }
}
/* =============================================================================
 * 功能：利用 GF 组合 (Geometry-Free) 探测周跳，阈值为 thresslip（米）。
 * 对每个次频与首频分别检测；无合适的单信号检测依据时标记组合的两频。
 * 若 PREPROC/DOPPSLIP 已开启且两频有可用 Doppler 历史，GF 跳变交给
 * 单信号连续性检测归因，避免一条坏次频把正常首频模糊度一起重置。
 * ============================================================================= */
static void detslp_gf(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
    double g0, g1, thres = rtk->opt.thresslip;
    int i, f, sat;

    if (thres <= 0.0) return;

    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        for (f = 1; f < rtk->opt.nf && f < NFREQ; f++) {
            if (!ppp_signal_supported(satsys(sat, NULL), f,
                                      obs[i].code[f], &rtk->opt)) continue;
            if ((g1 = gfmeas(obs + i, nav, f)) == 0.0) continue;
            g0 = rtk->ssat[sat - 1].gf[f - 1];
            rtk->ssat[sat - 1].gf[f - 1] = g1;
            if (g0 != 0.0 && fabs(g1 - g0) > thres) {
                ssat_t *ss = rtk->ssat + sat - 1;
                double dt0 = timediff(obs[i].time, ss->ppp_dop_time[0]);
                double dtf = timediff(obs[i].time, ss->ppp_dop_time[f]);

                /* When both single-signal Doppler continuity tests can run,
                 * let them attribute the jump. A noisy secondary signal must
                 * not automatically reset the good primary ambiguity. */
                if (pppopt_number(&rtk->opt, "-PREPROC=", 0.0, 0.0, 1.0) >= 0.5 &&
                    pppopt_number(&rtk->opt, "-DOPPSLIP=", 0.0, 0.0, 10.0) > 0.0 &&
                    ss->ppp_dop_valid[0] && ss->ppp_dop_valid[f] &&
                    obs[i].D[0] != 0.0 && obs[i].D[f] != 0.0 &&
                    fabs(obs[i].D[0]) < 20000.0 &&
                    fabs(obs[i].D[f]) < 20000.0 &&
                    dt0 >= 0.2 && dt0 <= 5.0 && dtf >= 0.2 && dtf <= 5.0) {
                    trace(3, "$GF_DEFER,sat=%d,F%d,dGF=%.3f,reason=DOPPLER_READY\n",
                          sat, f + 1, g1 - g0);
                    continue;
                }
                rtk->ssat[sat - 1].slip[0] |= LLI_SLIP;
                rtk->ssat[sat - 1].slip[f] |= LLI_SLIP;
                trace(2, "detslp_gf: slip sat=%2d L1-L%d dGF=%.3f m\n",
                    sat, f + 1, g1 - g0);
            }
        }
    }
}
/* =============================================================================
 * 功能：利用 MW 组合 (Melbourne-Wubbena) 探测宽巷周跳。
 * -MWTHRES>0 才启用，阈值米；坏码、时间中断、已有周跳均不重复归因。
 * tested 两频均为真时表示本历元两条相位/Doppler 检测实际执行过，
 * 此时 MW 仅报告而不再同时重置两频；不把手机伪距跳变当成载波周跳。
 * ============================================================================= */
static void detslp_mw(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav,
                     const unsigned char tested[MAXSAT][NFREQ])
{
    double g0, g1, dt;
    double thres = pppopt_number(&rtk->opt, "-MWTHRES=", 0.0, 0.0, 100.0);
    int i, f, sat;

    if (thres <= 0.0) return;
    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        if (sat < 1 || sat > MAXSAT) continue;
        for (f = 1; f < rtk->opt.nf && f < NFREQ; f++) {
            ssat_t *ss = rtk->ssat + sat - 1;
            if (!ppp_signal_supported(satsys(sat, NULL), f,
                                      obs[i].code[f], &rtk->opt)) continue;
            if (ss->ppp_code_bad[0] || ss->ppp_code_bad[f]) continue;
            if ((g1 = mwmeas(obs + i, nav, f)) == 0.0) continue;
            g0 = ss->mw[f - 1];
            dt = timediff(obs[i].time, ss->pt[0][0]);
            ss->mw[f - 1] = g1;
            if (g0 == 0.0 || dt < 0.2 || dt > 5.0 ||
                (ss->slip[0] & LLI_SLIP) || (ss->slip[f] & LLI_SLIP)) {
                continue;
            }
            if (fabs(g1 - g0) > thres) {
                /* MW also contains noisy phone code. If BOTH independent
                 * carrier/Doppler tests ran, use their per-signal attribution
                 * instead of resetting two good ambiguities for a code jump. */
                if(tested && tested[sat-1][0] && tested[sat-1][f]) {
                    trace(3,"$MW_DEFER,sat=%d,L1-L%d,dMW=%.3f,reason=DOPPLER_TESTED\n",sat,f+1,g1-g0);
                    continue;
                }
                ss->slip[0] |= LLI_SLIP;
                ss->slip[f] |= LLI_SLIP;
                trace(2, "$PRE_MW_SLIP,sat=%d,L1-L%d,dMW=%.3f,thres=%.3f\n",
                    sat, f + 1, g1 - g0, thres);
            }
        }
    }
}

/* median for a small per-epoch work array ---------------------------------*/
/* 小数组中位数：在原数组内排序，所以需传工作副本，不可传原始观测数组。 */
static double ppp_median(double *a, int n)
{
    double v;
    int i, j;

    if (n <= 0) return 0.0;
    for (i = 1; i < n; i++) {
        v = a[i];
        for (j = i; j > 0 && a[j - 1] > v; j--) a[j] = a[j - 1];
        a[j] = v;
    }
    return n & 1 ? a[n / 2] : 0.5 * (a[n / 2 - 1] + a[n / 2]);
}

/* Causal raw-observation precheck -----------------------------------------
 * Carrier change and Doppler-integrated range change should agree apart from
 * receiver-clock motion. Remove the per-system median so a common phone clock
 * jump is not labelled as a slip on every satellite. The same prediction is
 * used to reject only very large raw-code jumps before code smoothing. */
/* 有历史的历元预检，先检测后更新历史：
 * 1. 重置本历元 slip/code_bad，再合并 LLI 和 GF 检测结果。
 * 2. 相位创新=(L_now-L_old)*lambda+lambda*D_avg*dt；
 *    码创新=P_now-(P_smooth_old-lambda*D_avg*dt)。
 * 3. 按系统去中位共同项，防止手机钟跳把整组卫星误报为周跳/粗差；
 *    用 distinct satellites 计同伴数量，不把同一星多频当作多颗星。
 * 4. 单信号相位超过 DOPPSLIP 标记周跳；码超过 CODEJUMP 标坏且清平滑。
 *    同组码相共同大跳变则清平滑的旧钟差基准，不硬改 GNSS 观测值。
 * 5. MW 作为补充；最后才推进 Doppler 历史，避免“当前减当前”的伪检测。
 * 仅用当前/过去观测，不访问未来历元或参考真值；缺历史不能完成此类检测。
 * STARTQC 启用时 30 s 之后清理旧 pt/ph 的规则按当前代码保留，
 * 此处只描述其行为，不在加注释任务中调整历史采样策略。 */
static void precheck_obs_ppp(rtk_t *rtk, const obsd_t *obs, int n,
                             const nav_t *nav)
{
    enum { NWORK = MAXOBS * NFREQ };
    double piv[NWORK], civ[NWORK], work[NWORK], medp[NCLK_PPP] = { 0 };
    double medc[NCLK_PPP] = { 0 }, phase_thres, code_thres;
    int poi[NWORK], pof[NWORK], pgrp[NWORK], coi[NWORK], cof[NWORK], cgrp[NWORK];
    int pcnt[NCLK_PPP] = { 0 }, ccnt[NCLK_PPP] = { 0 };
    unsigned char pseen[NCLK_PPP][MAXSAT]={{0}},cseen[NCLK_PPP][MAXSAT]={{0}};
    unsigned char dop_tested[MAXSAT][NFREQ]={{0}};
    int i, f, q, sat, sys, grp, np = 0, nc = 0;
    int history_done=pppopt_number(&rtk->opt,"-STARTQC=",1.0,0.0,1.0)>=0.5 &&
        rtk->ppp_start_time.time && fabs(timediff(obs[0].time,rtk->ppp_start_time))>30.0;

    phase_thres = pppopt_number(&rtk->opt, "-DOPPSLIP=", 0.0, 0.0, 10.0);
    code_thres = pppopt_number(&rtk->opt, "-CODEJUMP=", 0.0, 0.0, 1000.0);

    for (sat = 0; sat < MAXSAT; sat++) for (f = 0; f < NFREQ; f++) {
        rtk->ssat[sat].slip[f] = 0;
        rtk->ssat[sat].ppp_code_bad[f] = 0;
        if(history_done) {
            rtk->ssat[sat].pt[0][f].time=0;rtk->ssat[sat].pt[0][f].sec=0.0;
            rtk->ssat[sat].ph[0][f]=0.0;
        }
    }
    detslp_ll(rtk, obs, n);
    detslp_gf(rtk, obs, n, nav);

    if (phase_thres > 0.0 || code_thres > 0.0) {
        for (i = 0; i < n && i < MAXOBS; i++) {
            sat = obs[i].sat;
            if (sat < 1 || sat > MAXSAT ||
                (grp = ppp_clk_index(sys = satsys(sat, NULL))) < 0) continue;
            for (f = 0; f < rtk->opt.nf && f < NFREQ; f++) {
                ssat_t *ss = rtk->ssat + sat - 1;
                double freq, lambda, dt, avgd, innov;

                if (!ppp_signal_supported(sys, f, obs[i].code[f],
                                          &rtk->opt)) continue;

                freq = sat2freq(sat, obs[i].code[f], nav);
                if (freq == 0.0 || obs[i].D[f] == 0.0 ||
                    fabs(obs[i].D[f]) > 20000.0 || !ss->ppp_dop_valid[f] ||
                    fabs(ss->ppp_doppler[f]) > 20000.0) continue;
                dt = timediff(obs[i].time, ss->ppp_dop_time[f]);
                if (dt < 0.2 || dt > 5.0) continue;
                lambda = CLIGHT / freq;
                avgd = 0.5 * (ss->ppp_doppler[f] + obs[i].D[f]);

                if (phase_thres > 0.0 && obs[i].L[f] != 0.0 &&
                    ss->ph[0][f] != 0.0 &&
                    fabs(timediff(ss->ppp_dop_time[f], ss->pt[0][f])) < 0.05 &&
                    !(ss->slip[f] & LLI_SLIP) && np < NWORK) {
                    innov = (obs[i].L[f] - ss->ph[0][f]) * lambda +
                        lambda * avgd * dt;
                    piv[np] = innov; poi[np] = i; pof[np] = f; pgrp[np] = grp;
                    if(!pseen[grp][sat-1]){pcnt[grp]++;pseen[grp][sat-1]=1;} np++;
                }
                if (code_thres > 0.0 && obs[i].P[f] != 0.0 &&
                    ss->psmvalid[f] && nc < NWORK) {
                    innov = obs[i].P[f] -
                        (ss->psmP[f] - lambda * avgd * dt);
                    civ[nc] = innov; coi[nc] = i; cof[nc] = f; cgrp[nc] = grp;
                    if(!cseen[grp][sat-1]){ccnt[grp]++;cseen[grp][sat-1]=1;} nc++;
                }
            }
        }
        for (grp = 0; grp < NCLK_PPP; grp++) {
            for (i = q = 0; i < np; i++) if (pgrp[i] == grp) work[q++] = piv[i];
            if (q >= 3) medp[grp] = ppp_median(work, q);
            for (i = q = 0; i < nc; i++) if (cgrp[i] == grp) work[q++] = civ[i];
            if (q >= 3) medc[grp] = ppp_median(work, q);
            /* A code+phase common clock step is not a constellation-wide
             * cycle slip. Reset the code smoother's clock datum instead of
             * propagating last epoch's offset into current code for seconds. */
            if(pcnt[grp]>=3&&ccnt[grp]>=3&&code_thres>0.0&&
                fabs(medc[grp])>code_thres&&
                fabs(medp[grp]-medc[grp])<MAX(5.0,code_thres)) {
                for(i=0;i<n&&i<MAXOBS;i++)if(ppp_clk_index(satsys(obs[i].sat,NULL))==grp)
                    for(f=0;f<rtk->opt.nf&&f<NFREQ;f++)rtk->ssat[obs[i].sat-1].psmvalid[f]=0;
                trace(2,"$PRE_CLOCK_RESET,group=%d,code_common=%.3f,phase_common=%.3f,action=RESET_SMOOTHER\n",
                    grp,medc[grp],medp[grp]);
            }
        }
        for (i = 0; i < np; i++) if (pcnt[pgrp[i]] >= 3) {
            sat = obs[poi[i]].sat; f = pof[i];
            dop_tested[sat-1][f]=1;
            if(fabs(piv[i]-medp[pgrp[i]])<=phase_thres)continue;
            rtk->ssat[sat - 1].slip[f] |= LLI_SLIP;
            rtk->ssat[sat - 1].psmvalid[f] = 0;
            trace(2, "$PRE_DOP_SLIP,sat=%d,F%d,res=%.3f,common=%.3f,thres=%.3f\n",
                sat, f + 1, piv[i] - medp[pgrp[i]], medp[pgrp[i]], phase_thres);
        }
        for (i = 0; i < nc; i++) if (ccnt[cgrp[i]] >= 3 &&
            fabs(civ[i] - medc[cgrp[i]]) > code_thres) {
            sat = obs[coi[i]].sat; f = cof[i];
            rtk->ssat[sat - 1].ppp_code_bad[f] = 1;
            rtk->ssat[sat - 1].psmvalid[f] = 0;
            trace(2, "$PRE_CODE_REJECT,sat=%d,F%d,res=%.3f,common=%.3f,thres=%.3f\n",
                sat, f + 1, civ[i] - medc[cgrp[i]], medc[cgrp[i]], code_thres);
        }
    }
    detslp_mw(rtk, obs, n, nav,dop_tested);

    /* Advance Doppler history only after all innovations have been formed. */
    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        if (sat < 1 || sat > MAXSAT) continue;
        for (f = 0; f < rtk->opt.nf && f < NFREQ; f++) {
            ssat_t *ss = rtk->ssat + sat - 1;
            if (obs[i].D[f] == 0.0 || !isfinite(obs[i].D[f]) || fabs(obs[i].D[f]) > 20000.0) {
                ss->ppp_dop_valid[f] = 0;
                continue;
            }
            ss->ppp_dop_time[f] = obs[i].time;
            ss->ppp_doppler[f] = obs[i].D[f];
            ss->ppp_dop_valid[f] = 1;
        }
    }
}
/* No-history code check before smoother/ambiguity seeding. Native code OSBs
 * and an atmospheric prior are applied, then clock/ISB/receiver-code common
 * offsets are removed by grouping on constellation AND actual signal code.
 * Conservative MAD and coarse-position uncertainty gates are code-only.
 * An unknown first-epoch carrier offset is an ambiguity, not evidence of a
 * cycle slip. No absolute phase/code consistency gate is used here. */
static double ppp_iono_mapf(const nav_t *nav,const double *pos,const double *azel);
static int ionex_vertical_prior(gtime_t time,const nav_t *nav,const double *pos,
    const double *azel,const prcopt_t *opt,double *ion_v,double *var_v);
/* 无历史的启动伪距检查：STARTQC 开启且使用 EST 电离层时才执行。
 * 利用 SPP 粗位置、卫星钟差、OSB 和宽松大气模型形成检查量；首频检查
 * P-(几何-卫星钟+大气)，有参考码的次频则检查“次频码-参考码”的差。
 * 每组按系统+实际码+参考码分开，去中位共同偏移以容纳钟差/ISB/RCB。
 * 普通检查至少 5 个同伴，码差检查至少 4 个；MAD 与位置不确定度给保守
 * 门限，下限 30 m。稀疏组不强行判粗差；只标 code_bad，不删卫星/相位。
 * 未知首历元载波偏移主要是模糊度，不能拿绝对 L-P 大小做首历元周跳判据。
 * 此函数与后续 PPP 验后拒绝不同：要在平滑和模糊度初值建立前阻止坏码播种。 */
static void precheck_startup_code_ppp(rtk_t *rtk,const obsd_t *obs,int n,
    const nav_t *nav,const double *rs,const double *dts,const double *var_rs,
    const int *svh)
{
    enum { NWORK=MAXOBS*NFREQ };
    double residual[NWORK],work[NWORK],pos[3],rr[3],e[3],azel[2],zero[NFREQ]={0};
    double L[NFREQ],P[NFREQ],Lc,Pc,position_gate;
    int oi[NWORK],of[NWORK],group[NWORK],code[NWORK],refcode[NWORK],done[NWORK]={0};
    int count=0,i,f,j,k,q;const prcopt_t *opt=&rtk->opt;
    if(pppopt_number(opt,"-STARTQC=",1.0,0.0,1.0)<0.5 ||
        opt->ionoopt!=IONOOPT_EST || (norm(rtk->x,3)>0.0 &&rtk->ppp_start_time.time &&
        fabs(timediff(obs[0].time,rtk->ppp_start_time))>30.0))return;
    for(i=0;i<3;i++)rr[i]=rtk->sol.rr[i];
    if(!isfinite(norm(rr,3))||norm(rr,3)<1E6)return;
    ecef2pos(rr,pos);
    position_gate=4.0*sqrt(MAX(0.0,rtk->sol.qr[0]+rtk->sol.qr[1]+rtk->sol.qr[2]));
    for(i=0;i<n&&i<MAXOBS;i++) {
        int sat=obs[i].sat,sys=satsys(sat,NULL),grp=ppp_clk_index(sys);
        double r,ion,vi,dtrop,vtr,zhd,zwd,mh,mw;
        if(grp<0||satexclude(sat,var_rs[i],svh[i],opt)||
            (r=geodist(rs+i*6,rr,e))<=0.0||!isfinite(r)||
            satazel(pos,e,azel)<opt->elmin)continue;
        if(!ionex_vertical_prior(obs[i].time,nav,pos,azel,opt,&ion,&vi))
            ion=ionmodel(obs[i].time,nav->ion_gps,pos,azel)/MAX(ppp_iono_mapf(nav,pos,azel),1.0);
        ion*=ppp_iono_mapf(nav,pos,azel);
        if(vmf3_trop(obs[i].time,pos,azel,&mh,&mw,&zhd,&zwd))dtrop=mh*zhd+mw*zwd;
        else if(!tropcorr(obs[i].time,nav,pos,azel,TROPOPT_SAAS,&dtrop,&vtr))continue;
        corr_meas(obs+i,nav,azel,opt,zero,zero,0.0,L,P,&Lc,&Pc);
        for(f=0;f<opt->nf&&f<NFREQ;f++) {
            double freq,y,freq0;int reference=0;
            if(P[f]==0.0||rtk->ssat[sat-1].ppp_code_bad[f]||
                !ppp_signal_supported(sys,f,obs[i].code[f],opt)||
                !ppp_code_modelled(sys,f,obs[i].code[f],opt)||
                (freq=sat2freq(sat,obs[i].code[f],nav))<=0.0)continue;
            y=P[f]-(r-CLIGHT*dts[i*2]+dtrop+SQR(FREQL1/freq)*ion);
            /* A secondary-minus-reference CODE pair cancels geometry, clock
             * and trop exactly. Its signal-specific receiver bias is removed
             * by the across-satellite median, not assumed to be zero. */
            if(f>0&&P[0]!=0.0&&!rtk->ssat[sat-1].ppp_code_bad[0]&&
                ppp_code_modelled(sys,0,obs[i].code[0],opt)&&
                (freq0=sat2freq(sat,obs[i].code[0],nav))>0.0) {
                y=P[f]-P[0]-(SQR(FREQL1/freq)-SQR(FREQL1/freq0))*ion;
                reference=obs[i].code[0];
            }
            if(!isfinite(y)||count>=NWORK)continue;
            residual[count]=y;oi[count]=i;of[count]=f;group[count]=grp;
            code[count]=obs[i].code[f];refcode[count]=reference;count++;
        }
    }
    for(i=0;i<count;i++) {
        double centre,mad,gate;
        if(done[i])continue;
        for(j=q=0;j<count;j++)if(group[j]==group[i]&&code[j]==code[i]&&refcode[j]==refcode[i]) {
            done[j]=1;work[q++]=residual[j];
        }
        /* A sparse constellation/signal cannot provide a robust consensus. */
        if(q<(refcode[i]?4:5))continue;
        centre=ppp_median(work,q);
        for(j=k=0;j<count;j++)if(group[j]==group[i]&&code[j]==code[i]&&refcode[j]==refcode[i])work[k++]=fabs(residual[j]-centre);
        mad=1.4826*ppp_median(work,k);
        gate=MAX(30.0,MAX(6.0*mad,refcode[i]?0.0:position_gate));
        for(j=0;j<count;j++)if(group[j]==group[i]&&code[j]==code[i]&&refcode[j]==refcode[i]&&fabs(residual[j]-centre)>gate) {
            int sat=obs[oi[j]].sat;char t[40],sid[8];f=of[j];
            if(refcode[j]&&rtk->ssat[sat-1].ppp_code_bad[0])continue;
            rtk->ssat[sat-1].ppp_code_bad[f]=1;rtk->ssat[sat-1].psmvalid[f]=0;
            time2str(obs[oi[j]].time,t,3);satno2id(sat,sid);
            trace(2,"$START_CODE_REJECT,%s,sat=%s,F%d,sig=%s,res=%.3f,common=%.3f,mad=%.3f,gate=%.3f,peers=%d,ref=%s\n",
                t,sid,f+1,code2obs(code[j]),residual[j]-centre,centre,mad,gate,q,refcode[j]?code2obs(refcode[j]):"GEOMETRY");
        }
    }
}
/* =============================================================================
 * 功能：卡尔曼滤波时间更新 (预测) -> 接收机位置、速度、加速度。
 * FIXED 模式使用 opt.ru；STATIC 保持上历元位置并增加 prn[5] 过程方差。
 * 无 dynamics 的运动模式每历元用 SPP 粗位置重建位置块；有 dynamics
 * 使用位置/速度/加速度转移 F，以实际 tt 传播 x 与 P，再加加速度噪声。
 * 初始位置来自 rtk->sol.rr，不来自分析用参考点；不同 mode 行为不能混读。
 * ============================================================================= */
static void udpos_ppp(rtk_t* rtk)
{
    double* F, * P, * FP, * x, * xp, pos[3], Q[9] = { 0 }, Qv[9], var = 0.0;
    int i, j, * ix, nx;

    trace(3, "udpos_ppp:\n");

    /* fixed mode */
    if (rtk->opt.mode == PMODE_PPP_FIXED) {
        for (i = 0; i < 3; i++) initx(rtk, rtk->opt.ru[i], 1E-8, i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x, 3) <= 0.0) {
        for (i = 0; i < 3; i++) initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        if (rtk->opt.dynamics) {
            for (i = 3; i < 6; i++) initx(rtk, rtk->sol.rr[i], VAR_VEL, i);
            for (i = 6; i < 9; i++) initx(rtk, 1E-6, VAR_ACC, i);
        }
    }
    /* static ppp mode */
    if (rtk->opt.mode == PMODE_PPP_STATIC) {
        for (i = 0; i < 3; i++) {
            rtk->P[i * (1 + rtk->nx)] += SQR(rtk->opt.prn[5]) * fabs(rtk->tt);
        }
        return;
    }
    /* kinematic mode without dynamics */
    if (!rtk->opt.dynamics) {
        for (i = 0; i < 3; i++) {
            initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        }
        return;
    }
    /* check variance of estimated position */
    for (i = 0; i < 3; i++) var += rtk->P[i + i * rtk->nx];
    var /= 3.0;

    if (var > VAR_POS) {
        /* reset position with large variance */
        for (i = 0; i < 3; i++) initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        for (i = 3; i < 6; i++) initx(rtk, rtk->sol.rr[i], VAR_VEL, i);
        for (i = 6; i < 9; i++) initx(rtk, 1E-6, VAR_ACC, i);
        trace(2, "reset rtk position due to large variance: var=%.3f\n", var);
        return;
    }
    /* generate valid state index */
    ix = imat(rtk->nx, 1);
    for (i = nx = 0; i < rtk->nx; i++) {
        if (i < 9 || (rtk->x[i] != 0.0 && rtk->P[i + i * rtk->nx] > 0.0)) ix[nx++] = i;
    }
    /* state transition of position/velocity/acceleration */
    F = eye(nx); P = mat(nx, nx); FP = mat(nx, nx); x = mat(nx, 1); xp = mat(nx, 1);

    for (i = 0; i < 6; i++) {
        F[i + (i + 3) * nx] = rtk->tt;
    }
    /* include accel terms if filter is converged */
    if (var < rtk->opt.thresar[1]) {
        for (i = 0; i < 3; i++) {
            F[i + (i + 6) * nx] = SQR(rtk->tt) / 2.0;
        }
    }
    else trace(3, "pos var too high for accel term: %.4f,%.4f\n", var, rtk->opt.thresar[1]);
    for (i = 0; i < nx; i++) {
        x[i] = rtk->x[ix[i]];
        for (j = 0; j < nx; j++) {
            P[i + j * nx] = rtk->P[ix[i] + ix[j] * rtk->nx];
        }
    }
    /* x=F*x, P=F*P*F+Q */
    matmul("NN", nx, 1, nx, F, x, xp);
    matmul("NN", nx, nx, nx, F, P, FP);
    matmul("NT", nx, nx, nx, FP, F, P);

    for (i = 0; i < nx; i++) {
        rtk->x[ix[i]] = xp[i];
        for (j = 0; j < nx; j++) {
            rtk->P[ix[i] + ix[j] * rtk->nx] = P[i + j * nx];
        }
    }
    /* process noise added to only acceleration */
    Q[0] = Q[4] = SQR(rtk->opt.prn[3]) * fabs(rtk->tt);
    Q[8] = SQR(rtk->opt.prn[4]) * fabs(rtk->tt);
    ecef2pos(rtk->x, pos);
    covecef(pos, Q, Qv);
    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
        rtk->P[i + 6 + (j + 6) * rtk->nx] += Qv[i + j * 3];
    }
    free(ix); free(F); free(P); free(FP); free(x); free(xp);
}
/* =============================================================================
 * 功能：卡尔曼滤波时间更新 (预测) -> 接收机钟差 (以及多系统系统时差)
 * 说明：仅 GPS/QZS 公共钟差每历元按白噪声重建，初值来自 SPP dtr[0]。
 *       GLO/GAL/BDS/IRN 存的是相对 GPS 的 ISB，保留跨历元连续性并加
 *       (0.10 m/sqrt(s))^2*abs(tt) 过程方差，不是每秒全部重置。
 *       内部钟差单位米，输出 dtr 是秒；非 GPS 方程依赖 GPS 钟差+本系统 ISB。
 * ============================================================================= */
static void udclk_ppp(rtk_t* rtk)
{
    int i;

    trace(3, "udclk_ppp:\n");

    /* GPS/QZS receiver clock: white-noise state */
    initx(rtk, CLIGHT * rtk->sol.dtr[0], VAR_CLK, IC(0, &rtk->opt));

    /* Inter-system biases: keep temporal continuity */
    for (i = 1; i < NCLK_PPP; i++) {
        int j = IC(i, &rtk->opt);

        if (rtk->x[j] == 0.0) {
            /*
             * RTKLIB uses x==0.0 as the inactive-state sentinel in the
             * reduced Kalman update. A literal zero would therefore leave an
             * ISB permanently unestimated and force ZTD/ambiguities to absorb
             * the constellation-dependent code offset.
             */
            initx(rtk, 1E-6, VAR_CLK, j);
        }
        else {
            /* 0.10 m/sqrt(s): deliberately loose for phone receiver */
            rtk->P[j + j * rtk->nx] += SQR(0.10) * fabs(rtk->tt);
        }
    }
}
/* =============================================================================
 * 功能：卡尔曼滤波时间更新 (预测) -> 对流层参数 (天顶延迟ZTD + 水平梯度)
 * 说明：第一个状态是总 ZTD（米），不是单独存储的 ZWD。
 *       初始化优先 VMF3 的 ZHD+ZWD，否则使用 SBAS 天顶延迟；ESTG 还建
 *       两个梯度。随后保持状态并以 prn[2]^2*abs(tt) 增加 ZTD 方差。
 *       STATMOD 位2启用时先验记账使用 obs GPST，不使用 SPP 改钟后的 sol.time；
 *       用过 VMF3 初始化后记录时间，避免同历元再次当作独立产品观测强化。
 * ============================================================================= */
static void udtrop_ppp(rtk_t* rtk, gtime_t epoch)
{
    double pos[3], azel[] = { 0.0,PI / 2.0 }, ztd, var;
    char tstr[64];
    int i = IT(&rtk->opt), j;
    /* sol.time can be clock-corrected by SPP: prior bookkeeping must use
       the same observation timestamp as subsequent pseudo-observations. */
    gtime_t prior_time=(statistical_model_ppp(&rtk->opt)&2)?epoch:rtk->sol.time;

    trace(3, "udtrop_ppp:\n");

    if (rtk->x[i] == 0.0) {
        ecef2pos(rtk->sol.rr, pos);
        time2str(prior_time, tstr, 0);
        if (vmf3_ztd_prior(prior_time, pos, &ztd)) {
            /* Start at the external forecast, but retain a conservative
             * variance because this is a 5-degree forecast grid. */
            var = SQR(vmf3_init_sigma(&rtk->opt));
            if (statistical_model_ppp(&rtk->opt)&2)
                rtk->ppp_zwd_prior_time=epoch;
            trace(2, "$VMF3_INIT,time=%s,ZTD=%.4f,sig=%.3f\n",
                tstr, ztd, sqrt(var));
        }
        else {
            rtk->ppp_zwd_prior_time.time=0;rtk->ppp_zwd_prior_time.sec=0.0;
            ztd = sbstropcorr(prior_time, pos, azel, &var);
            trace(2, "$VMF3_INIT,time=%s,fallback=SBAS,ZTD=%.4f\n",
                tstr, ztd);
        }
        initx(rtk, ztd, var, i);

        if (rtk->opt.tropopt >= TROPOPT_ESTG) {
            for (j = i + 1; j < i + 3; j++) initx(rtk, 1E-6, VAR_GRA, j);
        }
    }
    else {
        rtk->P[i + i * rtk->nx] += SQR(rtk->opt.prn[2]) * fabs(rtk->tt);

        if (rtk->opt.tropopt >= TROPOPT_ESTG) {
            for (j = i + 1; j < i + 3; j++) {
                rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[2] * 0.1) * fabs(rtk->tt);
            }
        }
    }
}
/* -IONCONS 是 IONEX 垂直 L1 延迟先验的标准差下限（米）。
 * 产品 RMS 缺失时不能给零方差；此参数不是 TECU，也不是斜路径标准差。 */
static double ion_constraint_sigma(const prcopt_t* opt)
{
    const char* p;
    double sig = ION_CONSTR_SIGMA_DEF;

    if (opt && (p = strstr(opt->pppopt, "-IONCONS="))) {
        if (sscanf(p, "-IONCONS=%lf", &sig) != 1 || !isfinite(sig)) {
            sig = ION_CONSTR_SIGMA_DEF;
        }
    }
    if (sig < ION_CONSTR_SIGMA_MIN) sig = ION_CONSTR_SIGMA_MIN;
    if (sig > ION_CONSTR_SIGMA_MAX) sig = ION_CONSTR_SIGMA_MAX;
    return sig;
}

/* Do not count a slowly varying GIM as a new independent observation at
 * every receiver epoch. Scaling R by interval/dt approximately preserves
 * one prior's information over each correlation interval. A value of 1 s
 * reproduces the former behavior for 1 Hz observations. */
/* -IONCONSINT 控制同一缓变 GIM 被重复应用的信息间隔；1 Hz 且设置
 * 1 s 时后续历元仍接近旧的每秒约束强度，并不自动推算 GIM 的实际相关时间。 */
static double ion_constraint_interval(const prcopt_t* opt)
{
    const char* p;
    double sec = ION_CONSTR_INTERVAL_DEF;

    if (opt && (p = strstr(opt->pppopt, "-IONCONSINT="))) {
        if (sscanf(p, "-IONCONSINT=%lf", &sec) != 1 || !isfinite(sec)) {
            sec = ION_CONSTR_INTERVAL_DEF;
        }
    }
    return MIN(86400.0, MAX(1.0, sec));
}

/* Keep the estimated vertical TEC state, its IONEX prior, and its design
 * matrix on the same single-layer geometry. iontec() uses the IONEX grid's
 * earth radius and shell height; ionmapf() instead hard-codes HION=350 km.
 * IONEX grids normally share this geometry across their map epochs. */
/* 统一垂直/斜路径几何：优先用 IONEX 的地球半径和单层高度求映射因子。
 * 初始化、GNSS 电离层项与 H 偏导必须采用同一映射，不能各自用不同壳层。
 * 产品不满足单层条件或不可用时回退 ionmapf；输出为无量纲。 */
static double ppp_iono_mapf(const nav_t* nav, const double* pos,
    const double* azel)
{
    double posp[3], mapf;
    const tec_t* tec;

    if (nav && nav->nt > 0 && nav->tec) {
        tec = nav->tec;
        if (tec->ndata[2] == 1 && tec->rb > 0.0 && tec->hgts[0] > 0.0) {
            mapf = ionppp(pos, azel, tec->rb, tec->hgts[0], posp);
            if (isfinite(mapf) && mapf > 0.0) return mapf;
        }
    }
    return ionmapf(pos, azel);
}

/* iontec() returns slant L1 delay. Convert it to the vertical L1 state used by
 * this PPP and impose a conservative variance floor because many GIM IONEX
 * products (including the current CODE file) contain TEC maps without RMS maps. */
/* iontec 给 L1 斜路径延迟和方差：除 mapf/mapf^2 得到垂直先验及方差。
 * 若产品没有 RMS，仍至少使用 IONCONS^2 的方差；异常数值/高度角返回失败。
 * 本函数只算外部先验，不更新 PPP 状态、不做重复信息时间记账。 */
static int ionex_vertical_prior(gtime_t time, const nav_t* nav,
    const double* pos, const double* azel, const prcopt_t* opt,
    double* ion_v, double* var_v)
{
    double ion_s = 0.0, var_s = 0.0, mapf, sig;

    if (!nav || !pos || !azel || azel[1] <= 0.0) return 0;
    if ((mapf = ppp_iono_mapf(nav, pos, azel)) <= 0.0) return 0;
    if (!iontec(time, nav, pos, azel, 1, &ion_s, &var_s)) return 0;

    *ion_v = ion_s / mapf;
    *var_v = var_s > 0.0 ? var_s / SQR(mapf) : 0.0;

    sig = ion_constraint_sigma(opt);
    if (*var_v < SQR(sig)) *var_v = SQR(sig);

    if (*ion_v != *ion_v || *var_v != *var_v || fabs(*ion_v) > 1E4 || *var_v <= 0.0) {
        return 0;
    }
    return 1;
}

/* =============================================================================
 * 功能：卡尔曼滤波时间更新 (预测) -> 电离层参数 (每颗卫星单独估计)
 * 说明：每星估一个参考 L1 的垂直延迟，不要求该星此刻两频码都齐全。
 *       初始化优先采用 IONEX，缺失时用广播电离层作宽松初值；不人为制造
 *       GNSS 测量。随后增加 (prn[1]/sin(el))^2*abs(tt) 的过程方差。
 *       仅所有可用频槽均超 GAP_RESION 中断计数时重置，避免首频暂缺误清状态。
 * ============================================================================= */
static void udiono_ppp(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
    double ion, vari, sinel, pos[3], * azel;
    char* p;
    int i, j, f, gap_resion = GAP_RESION, sat, outage;

    trace(3, "udiono_ppp (IONEX constrained):\n");

    if ((p = strstr(rtk->opt.pppopt, "-GAP_RESION="))) {
        sscanf(p, "-GAP_RESION=%d", &gap_resion);
    }

    /* Reset an ionosphere state only if every usable frequency of the satellite
     * has been absent for the configured gap. Using outc[0] alone incorrectly
     * resets L5/E5-only phone measurements when L1 is temporarily missing. */
    for (i = 0; i < MAXSAT; i++) {
        j = II(i + 1, &rtk->opt);
        outage = 1000000000;
        for (f = 0; f < NF(&rtk->opt) && f < NFREQ; f++) {
            if ((int)rtk->ssat[i].outc[f] < outage) outage = (int)rtk->ssat[i].outc[f];
        }
        if (rtk->x[j] != 0.0 && outage > gap_resion) {
            initx(rtk, 0.0, 0.0, j);
            rtk->ssat[i].ppp_ion_prior_time.time=0;
            rtk->ssat[i].ppp_ion_prior_time.sec=0.0;
            trace(3, "ion reset: sat=%d outage=%d\n", i + 1, outage);
        }
    }

    ecef2pos(rtk->sol.rr, pos);

    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        j = II(sat, &rtk->opt);
        azel = rtk->ssat[sat - 1].azel;

        if (rtk->x[j] == 0.0) {
            rtk->ssat[sat-1].ppp_ion_prior_time.time=0;
            rtk->ssat[sat-1].ppp_ion_prior_time.sec=0.0;
            /* First choice: initialize from the same IONEX product that will be
             * used as a soft constraint. This avoids noisy phone code-difference
             * initialization being contaminated by receiver IFB/code noise. */
            if (ionex_vertical_prior(obs[i].time, nav, pos, azel, &rtk->opt, &ion, &vari)) {
                if (fabs(ion) < 1E-8) ion = 1E-6; /* RTKLIB filter active-state sentinel */
                initx(rtk, ion, vari, j);
                if (statistical_model_ppp(&rtk->opt)&2)
                    rtk->ssat[sat-1].ppp_ion_prior_time=obs[i].time;
                trace(3, "ion init IONEX: sat=%d Iv=%.4f sig=%.3f\n",
                    sat, ion, sqrt(vari));
            }
            else {
                /* IONEX unavailable: use broadcast ionosphere only as a very
                 * loose initialization, never as a hard constraint. */
                double mapf = ppp_iono_mapf(nav, pos, azel);
                if (mapf <= 0.0) continue;
                ion = ionmodel(obs[i].time, nav->ion_gps, pos, azel) / mapf;
                if (ion != ion || fabs(ion) > 1E4) ion = 1E-6;
                if (fabs(ion) < 1E-8) ion = 1E-6;
                initx(rtk, ion, VAR_IONO, j);
                trace(3, "ion init BRDC: sat=%d Iv=%.4f sig=%.3f\n",
                    sat, ion, sqrt(VAR_IONO));
            }
        }
        else {
            sinel = sin(MAX(azel[1], 5.0 * D2R));
            rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[1] / sinel) * fabs(rtk->tt);
        }
    }
}
/* =============================================================================
 * 功能：初始化/预测按实际信号划分的接收机相对码偏差 RCB（代码沿用 uddcb 名称）。
 * 初值 1E-6 用于激活状态，方差 VAR_IFB；随后以 0.01^2*abs(tt) 随机游走。
 * 这里没有读“接收机偏差产品”；RCB 是利用观测估计的未知数，不是卫星 OSB。
 * ============================================================================= */
 /* 状态只对应前面枚举的信号，不是为所有次频自动建立通用公共偏差。 */
static void uddcb_ppp(rtk_t* rtk)
{
    int i, j;

    if (ND(&rtk->opt) <= 0) return;

    for (i = 0; i < ND(&rtk->opt); i++) {
        j = ID(i, &rtk->opt);
        if (j < 0 || j >= rtk->nx) continue;

        if (rtk->x[j] == 0.0) {
            /* A nonzero value activates the state in RTKLIB's reduced filter. */
            initx(rtk, 1E-6, VAR_IFB, j);
        }
        else {
            /* Receiver code hardware delay is nearly constant but not perfectly
               fixed over a phone session. */
            rtk->P[j + j * rtk->nx] += SQR(0.01) * fabs(rtk->tt);
        }
    }
}
/* =============================================================================
 * 功能：卡尔曼滤波时间更新 (预测) -> 载波相位模糊度 (Phase Biases)
 * 说明：非组合 PPP 的模糊度以米估计，为浮点相位偏置，不在这里固定整周。
 *       PREPROC 开启时沿用之前 LLI/GF/Doppler/MW 的检测结果；关闭时在这里
 *       使用旧的 LLI/GF 路径。MW 是否启用取决于 PREPROC 和 MWTHRES。
 *       超 maxout 中断、瞬时固定模式或配置的钟跳处理可清掉旧模糊度。
 *       初始化同频组合：B=L_f-(P_f-RCB_f)+2*(fL1/f)^2*I_slant。
 *       没有建模次频码时，可用有效参考码与两频电离层系数之和建立载波初值。
 *       码已经过同一 OSB 和平滑路径，不能用未改正码与已改正相位混配播种。
 *       有确认周跳时用当前历元重建该信号模糊度，并清相关；没有有效初值
 *       就暂不初始化，不通过删当前相位或写真值来“加速收敛”。
 * ============================================================================= */
static void udbias_ppp(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
    double L[NFREQ], P[NFREQ], Lc, Pc, bias[MAXOBS], offset = 0.0, pos[3] = { 0 };
    double dion = 0.0, vari = 0.0, gamma, freq, P_clean;
    double dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
    int i, j, k, f, f2, sat, sys, slip[MAXOBS] = { 0 }, clk_jump = 0, ifb_slot, id_ifb;
    char tstr[40];

    trace(3, "udbias  : n=%d\n", n);

    if (rtk->opt.posopt[5]) {
        clk_jump = ROUND(time2gpst(obs[0].time, NULL) * 10) % 864000 == 0;
    }
    if (pppopt_number(&rtk->opt, "-PREPROC=", 0.0, 0.0, 1.0) < 0.5) {
        /* Legacy path retained for exact A/B comparison. */
        for (i = 0; i < MAXSAT; i++) for (j = 0; j < rtk->opt.nf; j++) {
            rtk->ssat[i].slip[j] = 0;
            rtk->ssat[i].ppp_code_bad[j] = 0;
        }
        detslp_ll(rtk, obs, n);
        detslp_gf(rtk, obs, n, nav);
    }
    ecef2pos(rtk->sol.rr, pos);
    id_ifb = -1;

    for (f = 0; f < NF(&rtk->opt); f++) {
        offset = 0.0;

        for (i = 0; i < MAXSAT; i++) {
            uint32_t outc = ++rtk->ssat[i].outc[f];
            int reset_outage = outc > (uint32_t)rtk->opt.maxout;
            int reset_inst = rtk->opt.modear == ARMODE_INST;

            if (reset_outage || reset_inst || clk_jump) {
                char sid[8] = "";
                j = IB(i + 1, f, &rtk->opt);

                if (j >= 0 && j < rtk->nx && rtk->x[j] != 0.0) {
                    satno2id(i + 1, sid);
                    trace(2,
                        "$AMBRESET,%s,%s,F%d,%s,outc=%u,old=%.4f,std=%.4f\n",
                        time2str(rtk->sol.time, tstr, 2), sid, f + 1,
                        clk_jump ? "CLOCKJUMP" :
                        (reset_inst ? "INSTANT_AR" : "OUTAGE"),
                        (unsigned int)outc, rtk->x[j],
                        SQRT(rtk->P[j + j * rtk->nx]));
                }
                if (j >= 0 && j < rtk->nx) initx(rtk, 0.0, 0.0, j);
            }
        }

        for (i = k = 0; i < n && i < MAXOBS; i++) {
            sat = obs[i].sat;
            bias[i] = 0.0;
            slip[i] = 0;
            if (!ppp_signal_supported(satsys(sat, NULL), f,
                                      obs[i].code[f], &rtk->opt)) continue;
            if (ppp_phase_quarantined(&rtk->ssat[sat - 1], f, obs[i].time)) continue;
            j = IB(sat, f, &rtk->opt);
            slip[i] = rtk->ssat[sat - 1].slip[f] & LLI_SLIP;

            if (j < 0 || j >= rtk->nx) {
                trace(2, "udbias_ppp: invalid ambiguity index j=%d nx=%d sat=%d f=%d\n",
                    j, rtk->nx, sat, f);
                continue;
            }

            corr_meas(obs + i, nav, rtk->ssat[sat - 1].azel, &rtk->opt, dantr, dants,
                0.0, L, P, &Lc, &Pc);
            applydoppsm_ppp(rtk, obs + i, P);
            if (rtk->ssat[sat - 1].ppp_code_bad[f]) P[f] = 0.0;

            if (rtk->opt.ionoopt == IONOOPT_IFLC) {
                if (Lc == 0.0 || Pc == 0.0) continue;
                bias[i] = Lc - Pc;
                f2 = seliflc(rtk->opt.nf, satsys(sat, NULL));
                if (f2 > 0 && f2 < NFREQ) {
                    slip[i] = (rtk->ssat[sat - 1].slip[0] |
                        rtk->ssat[sat - 1].slip[f2]) & LLI_SLIP;
                }
            }
            else {
                if (f >= NFREQ || L[f] == 0.0) continue;
                if ((freq = sat2freq(sat, obs[i].code[f], nav)) == 0.0) continue;
                sys = satsys(sat, NULL);
                if (f > 0 && !ppp_code_modelled(sys, f, obs[i].code[f],
                                                &rtk->opt) &&
                    P[0] != 0.0 && !rtk->ssat[sat - 1].ppp_code_bad[0]) {
                    double freq0 = sat2freq(sat, obs[i].code[0], nav);
                    if (freq0 == 0.0) continue;
                    P_clean = P[0];
                    gamma = SQR(FREQL1 / freq) + SQR(FREQL1 / freq0);
                }
                else {
                    if (P[f] == 0.0) continue;
                    P_clean = P[f];
                    gamma = 2.0 * SQR(FREQL1 / freq);
                }

                /*
                 * Receiver code IFB for every non-reference frequency.
                 * Subtract it before forming the phase-code ambiguity seed.
                 */
                ifb_slot = ppp_ifb_index(sys, f, obs[i].code[f], &rtk->opt);
                id_ifb = (ifb_slot >= 0 && ifb_slot < ND(&rtk->opt)) ?
                    ID(ifb_slot, &rtk->opt) : -1;
                if (id_ifb >= 0 && id_ifb < rtk->nx) {
                    P_clean -= rtk->x[id_ifb];
                }

                dion = vari = 0.0;
                if (!model_iono(obs[i].time, pos, rtk->ssat[sat - 1].azel,
                    &rtk->opt, sat, rtk->x, nav, &dion, &vari)) {
                    continue;
                }
                /* L_f - P_ref needs (gamma_f + gamma_ref)*I; same-band
                 * code uses 2*gamma_f*I. */
                bias[i] = L[f] - P_clean + gamma * dion;
                if (bias[i] == 0.0) bias[i] = 1E-4;
            }

            if (rtk->x[j] == 0.0 || slip[i] || bias[i] == 0.0) continue;
            offset += bias[i] - rtk->x[j];
            k++;
        }

        /* 传统码相共同大跳变保护：需至少两个已有模糊度，门限为 0.5 ms*c。
         * 这是异常量级的共同偏移处理，不是把普通米级残差平均后改坐标。 */
        if (k >= 2 && fabs(offset / k) > 0.0005 * CLIGHT) {
            for (i = 0; i < MAXSAT; i++) {
                j = IB(i + 1, f, &rtk->opt);
                if (j >= 0 && j < rtk->nx && rtk->x[j] != 0.0) rtk->x[j] += offset / k;
            }
            trace(2, "phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
                time2str(rtk->sol.time, tstr, 0), k, offset / k / CLIGHT);
        }

        for (i = 0; i < n && i < MAXOBS; i++) {
            sat = obs[i].sat;
            if (ppp_phase_quarantined(&rtk->ssat[sat - 1], f, obs[i].time)) continue;
            j = IB(sat, f, &rtk->opt);
            if (j < 0 || j >= rtk->nx) continue;

            if (rtk->P[j + j * rtk->nx] > 0.0) {
                rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[0]) * fabs(rtk->tt);
            }
            if (bias[i] == 0.0) continue;

            /* A confirmed cycle slip changes the ambiguity itself. Reinitialize
               from this epoch instead of retaining the old value. */
            if (rtk->x[j] == 0.0 || slip[i]) {
                char sid[8] = "";
                double old_bias = rtk->x[j];
                double old_std = SQRT(rtk->P[j + j * rtk->nx]);

                satno2id(sat, sid);
                if (slip[i]) {
                    trace(2,
                        "$AMBRESET,%s,%s,F%d,SLIP,outc=%u,old=%.4f,std=%.4f,new=%.4f\n",
                        time2str(rtk->sol.time, tstr, 2), sid, f + 1,
                        (unsigned int)rtk->ssat[sat - 1].outc[f],
                        old_bias, old_std, bias[i]);
                    /* udbias_ppp clears the slip flag below, so count it here. */
                    rtk->ssat[sat - 1].slipc[f]++;
                }
                else {
                    trace(4, "$AMBINIT,%s,%s,F%d,bias=%.4f\n",
                        time2str(rtk->sol.time, tstr, 2), sid, f + 1, bias[i]);
                }

                initx(rtk, bias[i], VAR_BIAS, j);
                rtk->ssat[sat - 1].slip[f] &= (uint8_t)~LLI_SLIP;
                for (k = 0; k < MAXSAT; k++) rtk->ambc[sat - 1].flags[k] = 0;
            }
        }
    }
}
/* =============================================================================
 * 功能：所有状态量的时间更新/初始化总入口，一历元只调用一次。
 * 顺序是位置 -> GPS钟差/ISB -> ZTD/梯度 -> 电离层 -> RCB -> 模糊度。
 * 模糊度初值依赖当前电离层和 RCB，因此不能随意把 udbias 移到它们之前。
 * “时间更新”以旧状态预测新历元并加过程噪声，与后面的观测测量更新不同。
 * ============================================================================= */
static void udstate_ppp(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
    trace(3, "udstate_ppp: n=%d\n", n);

    /* temporal update of position */
    udpos_ppp(rtk);

    /* temporal update of clock */
    udclk_ppp(rtk);

    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
        udtrop_ppp(rtk,obs[0].time);
    }
    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt == IONOOPT_EST) {
        udiono_ppp(rtk, obs, n, nav);
    }
    /* 预测当前已配置的各系统/实际信号 RCB，不仅仅是 GPS L5。 */
    if (rtk->opt.nf >= 3) {
        uddcb_ppp(rtk);
    }
    /* temporal update of phase-bias */
    udbias_ppp(rtk, obs, n, nav);
}
/* =============================================================================
 * 功能：根据卫星到接收机方向与卫星天底方向的夹角，查询卫星天线 PCV。
 * 结果 dant 是各频槽的距离改正；精密轨道参考点/PCO 处理还涉及其他模块，
 * 不应把此处 PCV 查询当作全部天线改正。天线型号来自已加载的 ANTEX。
 * ============================================================================= */
static void satantpcv(const double* rs, const double* rr, const pcv_t* pcv,
    double* dant)
{
    double ru[3], rz[3], eu[3], ez[3], nadir, cosa;
    int i;

    for (i = 0; i < 3; i++) {
        ru[i] = rr[i] - rs[i];
        rz[i] = -rs[i];
    }
    if (!normv3(ru, eu) || !normv3(rz, ez)) return;

    cosa = dot3(eu, ez);
    cosa = cosa < -1.0 ? -1.0 : (cosa > 1.0 ? 1.0 : cosa);
    nadir = acos(cosa);

    antmodel_s(pcv, nadir, dant);
}
/* =============================================================================
 * 功能：精密对流层投影函数与偏导数模型计算
 * 说明：x[0] 为总 ZTD，zhd 为模型干延迟，所以 ZWD_est=x[0]-zhd。
 *       斜路径延迟 T=m_h*ZHD+m_w*ZWD_est，两者对伪距/载波同号。
 *       VMF3 可用时取其 ZHD/映射，不直接用产品 ZWD 替换 GNSS 湿延迟；
 *       不可用时用 Saastamoinen 干延迟及 tropmapf 映射。
 *       x[1/2] 为当前实现的梯度参数，通过修正湿映射进入模型。
 *       返回 dtdx 供 H 使用：ZTD 偏导是修正后的 m_w，梯度偏导还乘 ZWD。
 *       model_trop 在 EST 模式把梯度设 0，ESTG 才读取三维状态。
 * ============================================================================= */
static double trop_model_prec(gtime_t time, const double* pos,
    const double* azel, const double* x, double* dtdx,
    double* var)
{
    const double zazel[] = { 0.0,PI / 2.0 };
    double zhd, zwd_init, m_h, m_w, cotz, grad_n, grad_e;

    /*
     * VMF3 provides the hydrostatic delay and mapping factors. The forecast
     * ZWD is intentionally not used here: (x[0]-zhd) is the locally estimated
     * PPP wet zenith delay and remains free to evolve as a random walk.
     */
    if (!vmf3_trop(time, pos, azel, &m_h, &m_w, &zhd, &zwd_init)) {
        zhd = tropmodel(time, pos, zazel, 0.0);
        m_h = tropmapf(time, pos, azel, &m_w);
    }

    if (azel[1] > 0.0) {

        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz = 1.0 / tan(azel[1]);
        grad_n = m_w * cotz * cos(azel[0]);
        grad_e = m_w * cotz * sin(azel[0]);
        m_w += grad_n * x[1] + grad_e * x[2];
        dtdx[1] = grad_n * (x[0] - zhd);
        dtdx[2] = grad_e * (x[0] - zhd);
    }
    dtdx[0] = m_w;
    *var = SQR(0.01);
    return m_h * zhd + m_w * (x[0] - zhd);
}
/* 按 tropopt 分派：SAAS/SBAS 是纯模型改正；EST/ESTG 从滤波 x 取状态。
 * EST 下本地 trp[3] 的梯度初值为 0，避免误读不存在的两个梯度状态。
 * 返回值是模型可用标志，dtrp 是米，var 是模型残差方差（m^2）。 */
static int model_trop(gtime_t time, const double* pos, const double* azel,
    const prcopt_t* opt, const double* x, double* dtdx,
    const nav_t* nav, double* dtrp, double* var)
{
    double trp[3] = { 0 };

    if (opt->tropopt == TROPOPT_SAAS) {
        *dtrp = tropmodel(time, pos, azel, REL_HUMI);
        *var = SQR(ERR_SAAS);
        return 1;
    }
    if (opt->tropopt == TROPOPT_SBAS) {
        *dtrp = sbstropcorr(time, pos, azel, var);
        return 1;
    }
    if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
        matcpy(trp, x + IT(opt), opt->tropopt == TROPOPT_EST ? 1 : 3, 1);
        *dtrp = trop_model_prec(time, pos, azel, trp, dtdx, var);
        return 1;
    }
    return 0;
}
/* =============================================================================
 * 功能：返回参考 L1 的斜路径电离层延迟（米）及模型方差。
 * EST 时 dion=垂直状态*映射因子，频率缩放/码相符号在 ppp_res 中处理。
 * EST 的 var=0 不是说电离层没有误差，而是其不确定度由状态协方差 P 传播。
 * TEC/BRDC/SBAS 是相应产品或模型路径；IFLC 返回 0，因为组合已消一阶电离层。
 * ============================================================================= */
static int model_iono(gtime_t time, const double* pos, const double* azel,
    const prcopt_t* opt, int sat, const double* x,
    const nav_t* nav, double* dion, double* var)
{
    if (opt->ionoopt == IONOOPT_SBAS) {
        return sbsioncorr(time, nav, pos, azel, dion, var);
    }
    if (opt->ionoopt == IONOOPT_TEC) {
        if (iontec(time, nav, pos, azel, 1, dion, var)) return 1;

        /*
         * IONEX gaps must not remove the whole satellite from a phone PPP
         * epoch. Fall back to broadcast ionosphere with deliberately loose
         * variance; the stochastic model will down-weight it.
         */
        *dion = ionmodel(time, nav->ion_gps, pos, azel);
        *var = SQR(MAX(fabs(*dion) * ERR_BRDCI, 2.0));
        trace(3, "IONEX unavailable: fallback to broadcast ionosphere\n");
        return 1;
    }
    if (opt->ionoopt == IONOOPT_BRDC) {
        *dion = ionmodel(time, nav->ion_gps, pos, azel);
        *var = SQR(*dion * ERR_BRDCI);
        return 1;
    }
    if (opt->ionoopt == IONOOPT_EST) {
        /* Estimated delay is a vertical delay, apply the mapping function. */
        *dion = x[II(sat, opt)] * ppp_iono_mapf(nav, pos, azel);
        *var = 0.0;
        return 1;
    }
    if (opt->ionoopt == IONOOPT_IFLC) {
        *dion = *var = 0.0;
        return 1;
    }
    return 0;
}
/* 验后诊断快照：每次 ppp_res 调用重新填充，最终通过编辑的一份用于输出。
 * pppos 另存第一次验后的 evidence（也含之后可能拒绝的候选）用于持续异常
 * 判别，避免一次坏观测被剔除后反而完全失去其质量证据。
 * phase_ok/code_ok 表示本次残差方程已构建，并不意味着所有试算快照都已接受。
 * physical_sigma 保存适应/抗差膨胀之前的 sigma；CMC 用未平滑的改正观测。 */
typedef struct {
    double phase_res, code_res;
    double physical_sigma[2]; /* before temporal/adaptive/robust penalties */
    double cmc, iono, snr, el;
    uint8_t phase_ok, code_ok, cmc_ok, signal;
} ppp_diag_signal_t;
typedef struct {
    ppp_diag_signal_t obs[MAXOBS][NFREQ];
} ppp_diag_epoch_t;

/* Five bad tested epochs spanning >=4 s, ten good epochs spanning >=9 s.
 * Neutral/untested epochs are not evidence of recovery. A single outlier
 * cannot trigger a persistent penalty. Called at most once per input epoch. */
/* 持续异常的滞回状态机：5 个坏检测且跨度>=4 s 进入 weak；10 个好检测
 * 且跨度>=9 s 恢复。重复时间不计第二次，倒序/大中断清历史。
 * neutral 会清当前计数但不直接解除已有 weak，防止未检测被当作恢复。
 * 返回是否切换 weak，供日志使用；本函数不直接改 PPP 观测值。 */
static int signal_quality_step_ppp(ppp_signal_quality_t *q, gtime_t time,
                                  int bad, int good)
{
    int old = q->weak;
    double dt = q->time.time ? timediff(time,q->time) : 0.0;
    if (q->time.time && dt == 0.0) return 0;
    if (q->time.time && (dt < 0.0 || dt > 5.0)) memset(q,0,sizeof(*q));
    q->time = time;
    if (bad) {
        q->good_count=0;q->good_since.time=0;
        if (!q->bad_count) q->bad_since=time;
        if (q->bad_count<255) q->bad_count++;
        if (q->bad_count>=5 && timediff(time,q->bad_since)>=4.0) q->weak=1;
    }
    else if (good) {
        q->bad_count=0;q->bad_since.time=0;
        if (!q->good_count) q->good_since=time;
        if (q->good_count<255) q->good_count++;
        if (q->good_count>=10 && timediff(time,q->good_since)>=9.0) q->weak=0;
    }
    else {
        q->bad_count=q->good_count=0;
        q->bad_since.time=q->good_since.time=0;
    }
    return old != q->weak;
}

/* weak 时码方差乘 4；相位默认只监测，需 SIGQCPH 显式启用才自动乘 4。
 * SIGQC 关闭则不施加惩罚；严重载波的隔离逻辑是另一个独立机制。 */
static int persistent_variance_factor_ppp(const prcopt_t *opt,
                                          const ppp_signal_quality_t *q,int type)
{
    if(pppopt_number(opt,"-SIGQC=",1.0,0.0,1.0)<0.5 || !q->weak) return 1;
    /* Carrier residuals are entangled with ambiguity/ionosphere states.
     * Monitor moderate persistent phase errors by default; automatic phase
     * reweighting is opt-in. Existing severe-phase quarantine remains on. */
    if(!type && pppopt_number(opt,"-SIGQCPH=",0.0,0.0,1.0)<0.5) return 1;
    return 4;
}

/* Spatial common-mode removal is diagnostic only: never edit residuals used
 * by the EKF. Compare like constellation+actual code+observation type, with
 * >=5 distinct satellites. Robust group scatter prevents common atmosphere,
 * clock/RCB errors from being labelled as individual signal failures. */
/* 用同系统+同实际 signal+同观测类型的验后残差去空间共同中位项。
 * 至少 5 颗不同卫星才检测；这样不把共同钟差/RCB/大气误差全归因于单星。
 * 依据物理 sigma 和组内 MAD 判持续异常，把本历元证据记入下历元的质量状态。
 * 去中位数仅用于检测，不修改送进 EKF 的 v，不能用此操作消掉真实模型偏差。
 * pppos 传入第一次验后快照，可能包含后来被硬拒绝的候选，并只计一次历元。 */
static void update_signal_quality_ppp(rtk_t *rtk, const obsd_t *obs, int n,
                                      const ppp_diag_epoch_t *evidence)
{
    int i,f,t,k,j,sys,peers,seen[MAXSAT];
    double residuals[MAXOBS],deviations[MAXOBS],center,mad,res,sigma,z;
    char sid[8],str[40];
    if(pppopt_number(&rtk->opt,"-SIGQC=",1.0,0.0,1.0)<0.5) return;
    for(i=0;i<n && i<MAXOBS;i++) for(f=0;f<NF(&rtk->opt)&&f<NFREQ;f++) {
        int sat=obs[i].sat;
        const ppp_diag_signal_t *d=&evidence->obs[i][f];
        ssat_t *ss=&rtk->ssat[sat-1];
        if(!d->signal || d->signal!=ss->ppp_track_code[f])continue;
        sys=satsys(sat,NULL);
        for(t=0;t<2;t++) {
            if(!t && ((obs[i].LLI[f]|ss->slip[f]) & LLI_SLIP)) {
                memset(&ss->ppp_res_quality[f][t],0,sizeof(ss->ppp_res_quality[f][t]));
                continue;
            }
            if(!(t?d->code_ok:d->phase_ok))continue;
            peers=0;memset(seen,0,sizeof(seen));
            for(k=0;k<n&&k<MAXOBS;k++) {
                int s=obs[k].sat;
                if(satsys(s,NULL)!=sys||seen[s-1])continue;
                for(j=0;j<NF(&rtk->opt)&&j<NFREQ;j++) {
                    const ppp_diag_signal_t *v=&evidence->obs[k][j];
                    if(v->signal!=d->signal||!(t?v->code_ok:v->phase_ok))continue;
                    residuals[peers++]=t?v->code_res:v->phase_res;
                    seen[s-1]=1;break;
                }
            }
            if(peers<5)continue;
            center=ppp_median(residuals,peers);
            for(k=0;k<peers;k++)deviations[k]=fabs(residuals[k]-center);
            mad=1.4826*ppp_median(deviations,peers);
            sigma=d->physical_sigma[t];
            if(!isfinite(sigma)||sigma<=0.0)continue;
            res=(t?d->code_res:d->phase_res)-center;z=fabs(res)/sigma;
            if(signal_quality_step_ppp(&ss->ppp_res_quality[f][t],obs[i].time,
                fabs(res)>MAX(4.0*sigma,6.0*mad),z<2.0)) {
                satno2id(sat,sid);time2str(obs[i].time,str,2);
                trace(2,"$SIG_QUALITY,%s,sat=%s,sig=%s,F%d,type=%s,centered_res=%.4f,sigma=%.4f,z=%.2f,MAD=%.4f,peers=%d,weak=%d,var_factor=%d\n",
                    str,sid,code2obs(d->signal),f+1,t?"CODE":"PHASE",
                    res,sigma,z,mad,peers,ss->ppp_res_quality[f][t].weak,
                    persistent_variance_factor_ppp(&rtk->opt,&ss->ppp_res_quality[f][t],t));
            }
        }
    }
}

/* Use only accepted, corrected, unsmoothed code/phase pairs from the final
 * filter update. This is a one-epoch-late, causal quality signal; it cannot
 * force a position toward a known reference or use future observations. */
/* 可选码相漂移降权，SIGCMCDRIFT 默认 0（关闭）。
 * CMC=P_corrected-L_corrected-2*I_est，仍含模糊度/接收机码偏差及
 * 剩余电离层/多路径；只能同一连续弧内比较，不能宣称它是纯码误差。
 * 换信号/周跳/中断重建基准；成熟弧相对基准的漂移用于下一历元码方差。
 * 本次不打开该研究开关，不改滤波模型。 */
static void update_cmc_quality_ppp(rtk_t *rtk, const obsd_t *obs, int n,
                                   const ppp_diag_epoch_t *diag)
{
    /* Code-minus-carrier still depends on the estimated ionosphere. Keep
     * this experimental reweighting opt-in until it is validated against
     * independent antenna coordinates on multiple receivers. */
    double threshold = pppopt_number(&rtk->opt, "-SIGCMCDRIFT=",
                                      0.0, 0.0, 30.0);
    int i, f, sat;

    if (!diag) return;
    if (threshold <= 0.0) {
        /* A live option change must stop applying any old CMC penalties. */
        for (i = 0; i < n && i < MAXOBS; i++) {
            sat = obs[i].sat;
            if (sat >= 1 && sat <= MAXSAT) {
                memset(rtk->ssat[sat - 1].ppp_cmc_code_weak, 0,
                       sizeof(rtk->ssat[sat - 1].ppp_cmc_code_weak));
            }
        }
        return;
    }
    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        if (sat < 1 || sat > MAXSAT) continue;
        for (f = 0; f < NF(&rtk->opt) && f < NFREQ; f++) {
            const ppp_diag_signal_t *d = &diag->obs[i][f];
            ssat_t *ss = rtk->ssat + sat - 1;
            double gap, drift;
            int weak;
            char sid[8], str[40];

            if (!d->cmc_ok || !d->phase_ok || !d->code_ok ||
                d->signal != ss->ppp_track_code[f]) continue;
            gap = ss->ppp_cmc_quality_time[f].time ?
                timediff(obs[i].time, ss->ppp_cmc_quality_time[f]) : 0.0;
            if (ss->ppp_cmc_quality_count[f] == 0 || gap <= 0.0 ||
                gap > MAX(3.0, 2.5 * fabs(rtk->tt)) ||
                (obs[i].LLI[f] & LLI_SLIP)) {
                ss->ppp_cmc_quality_mean[f] = d->cmc;
                ss->ppp_cmc_quality_ref[f] = d->cmc;
                ss->ppp_cmc_quality_count[f] = 1;
                ss->ppp_cmc_code_weak[f] = 0;
            }
            else {
                ss->ppp_cmc_quality_mean[f] +=
                    0.10 * (d->cmc - ss->ppp_cmc_quality_mean[f]);
                if (ss->ppp_cmc_quality_count[f] < 255) {
                    ss->ppp_cmc_quality_count[f]++;
                }
                if (ss->ppp_cmc_quality_count[f] == 10) {
                    ss->ppp_cmc_quality_ref[f] = ss->ppp_cmc_quality_mean[f];
                }
                drift = ss->ppp_cmc_quality_mean[f] -
                        ss->ppp_cmc_quality_ref[f];
                weak = ss->ppp_cmc_code_weak[f];
                if (ss->ppp_cmc_quality_count[f] >= 30 &&
                    fabs(drift) > threshold) weak = 1;
                else if (fabs(drift) < threshold * 0.5) weak = 0;
                if (weak != ss->ppp_cmc_code_weak[f]) {
                    satno2id(sat, sid);
                    time2str(obs[i].time, str, 2);
                    trace(2, "$SIG_CMC_QUALITY,%s,sat=%s,sig=%s,F%d,drift=%.3f,thres=%.3f,code_weight=%s\n",
                          str, sid, code2obs(d->signal), f + 1, drift,
                          threshold, weak ? "DOWN_4X_VAR" : "NORMAL");
                    ss->ppp_cmc_code_weak[f] = (uint8_t)weak;
                }
            }
            ss->ppp_cmc_quality_time[f] = obs[i].time;
        }
    }
}

/* =============================================================================
 * 功能：构建 GNSS 观测和大气伪观测的 v/H/R（核心测量模型）。
 *
 * 非组合模型，所有距离项均为米：
 *   P_f=r+c*dt_GPS+ISB_sys-c*dt_sat+T+gamma_f*I_slant+RCB_f+noise
 *   L_f=r+c*dt_GPS+ISB_sys-c*dt_sat+T-gamma_f*I_slant+B_f+noise
 * gamma_f=(FREQL1/f)^2；L_f 已在 corr_meas 中由周转米并去天线/缠绕项。
 * 码/相位的几何、钟差、对流层同号；电离层反号；RCB 只进码，B 只进相位。
 * v=观测-计算值，H 存“计算值对状态”的偏导，不是 v 对状态的偏导。
 *
 * post=0：验前构造滤波方程，返回有效方程数 nv，并追加 IONEX/VMF3 软约束。
 * post>0：用试更新 xp 重算验后残差，返回 1 通过/0 需继续编辑；
 *         不追加外部伪观测，不把产品偏离作为整颗卫星拒绝依据。
 * H 是 nx*nv 按列存储，第 nv 列在 H[k+nx*nv]；R 为 nv*nv 对角方差阵。
 * exc[i] 是整星屏蔽；mexc[i*2*NFREQ+j] 是单条码/相位屏蔽，不能混用。
 * 观测序列为 {L1,P1,L2,P2,...}，j/2 是频槽，j%2 是观测类型。
 * 物理 sigma、持续异常倍率、时间相关处理、抗差膨胀分开保存，避免“降权后
 * sigma 变大”掩盖粗差证据；验后按原物理标准化残差选最坏的一条编辑。
 * ============================================================================= */
 /* phase and code residuals --------------------------------------------------*/
static int ppp_res(int post, const obsd_t* obs, int n, const double* rs,
    const double* dts, const double* var_rs, const int* svh,
    const double* dr, int* exc, unsigned char* mexc, const nav_t* nav,
    const double* x, rtk_t* rtk, double* v, double* H, double* R,
    double* azel, ppp_diag_epoch_t* diag, int *product_use)
{
    prcopt_t* opt = &rtk->opt;
    double y, r, cdtr, bias, rr[3], pos[3], e[3], dtdx[3], L[NFREQ], P[NFREQ], Lc, Pc;
    double var[MAXOBS * 2 * NFREQ + MAXSAT + 3], dtrp = 0.0, dion = 0.0, vart = 0.0, vari = 0.0, dcb, freq;
    double dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
    double ve[MAXOBS * 2 * NFREQ] = { 0 };
    double ze[MAXOBS * 2 * NFREQ] = { 0 };
    double se[MAXOBS * 2 * NFREQ] = { 0 };
    double snre[MAXOBS * 2 * NFREQ] = { 0 };
    double ele[MAXOBS * 2 * NFREQ] = { 0 };
    double ifbe[MAXOBS * 2 * NFREQ] = { 0 };
    char str[40], sid[8];
    const char *code_only_opt, *bds_code_opt, *bds_code_warm_opt;
    double code_only_factor = 1.0;
    double bds_code_factor = 1.0, bds_code_warm = 0.0, elapsed;
    int ne = 0, obsi[MAXOBS * 2 * NFREQ] = { 0 }, frqi[MAXOBS * 2 * NFREQ] = { 0 };
    int maxobs = -1, maxfrq = -1, rej = -1;
    int i, j, k, sat, sys, nv = 0, nx = rtk->nx, stat = 1, frq, code;
    int ionex_constraints = 0, vmf3_zwd_constraints = 0;
    int statmodel=statistical_model_ppp(opt);
    double res, C;

    time2str(obs[0].time, str, 2);
    if (!post && product_use) memset(product_use,0,sizeof(int)*(MAXSAT+3));
    if ((code_only_opt = strstr(opt->pppopt, "-CODENOPH=")) &&
        sscanf(code_only_opt, "-CODENOPH=%lf", &code_only_factor) != 1) {
        code_only_factor = 1.0;
    }
    code_only_factor = MIN(100.0, MAX(1.0, code_only_factor));
    if ((bds_code_opt = strstr(opt->pppopt, "-BDSCODEVAR=")) &&
        (sscanf(bds_code_opt, "-BDSCODEVAR=%lf", &bds_code_factor) != 1 ||
            !isfinite(bds_code_factor))) {
        bds_code_factor = 1.0;
    }
    bds_code_factor = MIN(100.0, MAX(1.0, bds_code_factor));
    if ((bds_code_warm_opt = strstr(opt->pppopt, "-BDSCODEWARM=")) &&
        (sscanf(bds_code_warm_opt, "-BDSCODEWARM=%lf", &bds_code_warm) != 1 ||
            !isfinite(bds_code_warm))) {
        bds_code_warm = 0.0;
    }
    bds_code_warm = MIN(600.0, MAX(0.0, bds_code_warm));
    /* Code is valuable while position and ambiguities are being initialized.
     * Increase the BeiDou code variance smoothly only after phase arcs have
     * begun to form, avoiding an abrupt information change in the EKF. */
    if (bds_code_warm > 0.0 && rtk->ppp_start_valid) {
        elapsed = MAX(0.0, timediff(obs[0].time, rtk->ppp_start_time));
        bds_code_factor = 1.0 + (bds_code_factor - 1.0) *
            MIN(1.0, elapsed / bds_code_warm);
    }
    if (diag) memset(diag, 0, sizeof(*diag));

    for (i = 0; i < MAXSAT; i++) {
        for (j = 0; j < opt->nf; j++) rtk->ssat[i].vsat[j] = 0;
    }

    /* dr 为本历元潮汐位移：参与几何模型，不直接把位移写进位置状态。
     * x 的位置按所选参考约定保留，当前观测位置 rr=x_position+dr。 */
    for (i = 0; i < 3; i++) rr[i] = x[i] + dr[i];
    ecef2pos(rr, pos);

    for (i = 0; i < n && i < MAXOBS; i++) {
        int sat_has_phase = 0;
        sat = obs[i].sat;

        for (j = 0; j < NF(opt) && j < NFREQ; j++) {
            if (obs[i].L[j] != 0.0) sat_has_phase = 1;
        }

        if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 ||
            satazel(pos, e, azel + i * 2) < opt->elmin) {
            exc[i] = 1;
            continue;
        }

        /*
         * Satellite-level exclusion is reserved for geometry/health/product
         * failures. Frequency-specific residual outliers are handled by mexc[]
         * so that one bad E5b/L5 measurement does not delete good L1/E5a data
         * from the same satellite.
         */
        if (!(sys = satsys(sat, NULL)) ||
            satexclude(sat, var_rs[i], svh[i], opt) || exc[i]) {
            exc[i] = 1;
            continue;
        }

        if (!model_trop(obs[i].time, pos, azel + i * 2, opt, x, dtdx, nav, &dtrp, &vart) ||
            !model_iono(obs[i].time, pos, azel + i * 2, opt, sat, x, nav, &dion, &vari)) {
            continue;
        }

        if (opt->posopt[0]) satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, dants);
        antmodel(opt->pcvr, opt->antdel[0], azel + i * 2, opt->posopt[1], dantr);

        if (!model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type,
            opt->posopt[2] ? 2 : 0, rs + i * 6, rr, &rtk->ssat[sat - 1].phw)) {
            continue;
        }

        /* 先统一码/相位单位与产品/天线改正，再保存未平滑 CMC，最后应用
         * 原始码域的平滑差值。诊断漂移与滤波用平滑码不应混成同一个量。 */
        corr_meas(obs + i, nav, azel + i * 2, &rtk->opt, dantr, dants,
            rtk->ssat[sat - 1].phw, L, P, &Lc, &Pc);
        if (diag && opt->ionoopt != IONOOPT_IFLC) {
            for (j = 0; j < NF(opt) && j < NFREQ; j++) {
                ppp_diag_signal_t* d = &diag->obs[i][j];
                double signal_freq, iono_scale;

                if (L[j] == 0.0 || P[j] == 0.0 ||
                    rtk->ssat[sat - 1].ppp_code_bad[j] ||
                    !obs[i].code[j] || obs[i].code[j] > MAXCODE) continue;
                if ((signal_freq = sat2freq(sat, obs[i].code[j], nav)) == 0.0) continue;
                iono_scale = SQR(FREQL1 / signal_freq) * dion;
                d->cmc = P[j] - L[j] - 2.0 * iono_scale;
                d->iono = iono_scale;
                d->snr = obs[i].SNR[j] * SNR_UNIT;
                d->el = azel[1 + i * 2] * R2D;
                d->signal = obs[i].code[j];
                d->cmc_ok = 1;
            }
        }
        applydoppsm_ppp(rtk, obs + i, P);

        /* stack phase and code residuals {L1,P1,L2,P2,...} */
        for (j = 0; j < 2 * NF(opt); j++) {
            int midx = i * 2 * NFREQ + j;
            double snr_real, sigma0, z0, penalty = 1.0;
            double ifb_used = 0.0;

            if (mexc && midx >= 0 && midx < MAXOBS * 2 * NFREQ && mexc[midx]) {
                continue;
            }

            C = 0.0;
            dcb = bias = 0.0;
            code = j % 2; /* 0=phase, 1=code */
            frq = j / 2;

            if (!ppp_signal_supported(sys, frq, obs[i].code[frq], opt)) continue;

            if (code && frq < NFREQ &&
                rtk->ssat[sat - 1].ppp_code_bad[frq]) continue;

            if (code && frq > 0 &&
                !ppp_code_modelled(sys, frq, obs[i].code[frq], opt)) {
                trace(4, "$SIG_MODEL_SKIP,%s,sat=%d,sig=%s,reason=NO_RCB_STATE\n",
                      str, sat, code2obs(obs[i].code[frq]));
                continue;
            }

            if (!code && ppp_phase_quarantined(&rtk->ssat[sat - 1], frq,
                obs[i].time)) continue;

            /* raPPPid-style startup: first solve a code+Doppler-smoothed
             * position. Introducing carrier ambiguities only after this
             * short bootstrap prevents their arbitrary initial values from
             * dominating the first few seconds of a phone PPP solution. */
            if (!code && rtk->ppp_code_warmup) continue;

            if (opt->ionoopt == IONOOPT_IFLC) {
                if ((y = code == 0 ? Lc : Pc) == 0.0) continue;
            }
            else {
                if (frq >= NFREQ || (y = code == 0 ? L[frq] : P[frq]) == 0.0) continue;
                if ((freq = sat2freq(sat, obs[i].code[frq], nav)) == 0.0) continue;
                /* 频率和正负号一起进入系数 C；dion 已是参考 L1 斜路径延迟。
                 * 后面 H 的电离层偏导还需乘 mapf，因状态本身是垂直延迟。 */
                C = SQR(FREQL1 / freq) * (code == 0 ? -1.0 : 1.0);
            }

            if (H) {
                for (k = 0; k < nx; k++) H[k + nx * nv] = 0.0;
                /* e 是接收机指向卫星的视线单位向量；接收机朝卫星移动，
                 * 几何距离变短，故计算距离对 ECEF 位置的偏导为 -e。 */
                for (k = 0; k < 3; k++) H[k + nx * nv] = -e[k];
            }

            k = ppp_clk_index(sys);
            if (k < 0 || k >= NCLK_PPP) continue;
            /*
             * State IC(0) is the GPS/QZS receiver clock. The remaining
             * clock states are inter-system biases relative to GPS.
             * Keep the computed observation and its H column consistent:
             * non-GPS observations depend on BOTH IC(0) and IC(k).
             */
            cdtr = x[IC(0, opt)];

            if (H) H[IC(0, opt) + nx * nv] = 1.0;

            if (k > 0) {
                cdtr += x[IC(k, opt)];
                if (H) H[IC(k, opt) + nx * nv] = 1.0;
            }

            if (H) {

                if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
                    for (k = 0; k < (opt->tropopt >= TROPOPT_ESTG ? 3 : 1); k++) {
                        H[IT(opt) + k + nx * nv] = dtdx[k];
                    }
                }
            }

            if (opt->ionoopt == IONOOPT_EST) {
                if (x[II(sat, opt)] == 0.0) continue;
                if (H) H[II(sat, opt) + nx * nv] = C * ppp_iono_mapf(nav, pos, azel + i * 2);
            }

            /* Receiver code bias for this exact non-reference signal. */
            if (code == 1 && ND(opt) > 0) {
                int ifb_slot = ppp_ifb_index(sys, frq, obs[i].code[frq], opt);
                int id_ifb = (ifb_slot >= 0 && ifb_slot < ND(opt)) ?
                    ID(ifb_slot, opt) : -1;

                if (id_ifb >= 0 && id_ifb < nx) {
                    ifb_used = x[id_ifb];
                    dcb += ifb_used;
                    if (H) H[id_ifb + nx * nv] = 1.0;
                }
            }

            if (code == 0) {
                if ((bias = x[IB(sat, frq, opt)]) == 0.0) continue;
                if (H) H[IB(sat, frq, opt) + nx * nv] = 1.0;
            }

            /* dts[i*2] 是卫星钟差秒，乘 c 后以负号进模型；cdtr 已为米。
             * dcb 只在码行非零，bias 只在相位行非零；残差 res 单位米。
             * 初始化模糊度时采用的 RCB/电离层符号必须与这一行一致。 */
            res = y - (r + cdtr - CLIGHT * dts[i * 2] + dtrp + C * dion + dcb + bias);
            if (v) v[nv] = res;
            if (diag && frq < NFREQ) {
                ppp_diag_signal_t* d = &diag->obs[i][frq];
                d->signal = obs[i].code[frq];
                if (code) {
                    d->code_res = res;
                    d->code_ok = 1;
                }
                else {
                    d->phase_res = res;
                    d->phase_ok = 1;
                }
            }

            if (code == 0) rtk->ssat[sat - 1].resc[frq] = res;
            else           rtk->ssat[sat - 1].resp[frq] = res;

            /*
             * Keep the physical (pre-robust) sigma for outlier diagnosis.
             * Robust down-weighting is used in the EKF R matrix, but must not
             * erase the evidence used by the post-fit outlier detector.
             */
            snr_real = (frq < NFREQ) ? obs[i].SNR[frq] * SNR_UNIT : 35.0;
            var[nv] = varerr(sat, sys, azel[1 + i * 2], snr_real,
                signal_noise_index_ppp(sys,obs[i].code[frq],frq),code,opt) *
                signal_weight_ppp(sys,obs[i].code[frq],code,opt);
            if (sys == SYS_CMP && code) var[nv] *= bds_code_factor;
            if (code && !sat_has_phase && !rtk->ppp_code_warmup) {
                var[nv] *= code_only_factor;
            }
            /* 模型误差合入 R：对流层+频率缩放后的电离层+卫星位置/钟差产品。
             * EST 电离层的不确定度由状态 P 处理，不在此重复加一份先验方差。
             * 当前 R 仍是对角近似，未显式建立同星码相/跨星产品误差相关。 */
            var[nv] += vart + SQR(C) * vari + var_rs[i];
            if (sys == SYS_GLO && code == 1) var[nv] += VAR_GLO_IFB;

            /* 此时冻结原物理 sigma0 和 z0。后面即使 R 被持续惩罚/抗差放大，
             * 硬拒绝仍使用 z0，不允许同一粗差靠不断膨胀 sigma 逃过检测。 */
            sigma0 = sqrt(MAX(var[nv], 1E-12));
            z0 = fabs(res) / sigma0;
            if (diag && frq < NFREQ)
                diag->obs[i][frq].physical_sigma[code]=sigma0;
            /* Never claim the white-noise variance reduction as independent
               PPP information: retain the raw-code floor. An experimental
               causal row-sum accounts for temporal covariance approximately,
               but is not equivalent to augmented coloured-noise filtering. */
            if (code && opt->ionoopt!=IONOOPT_IFLC && (statmodel&1) &&
                doppsm_factor(opt)>0.0 && rtk->ssat[sat-1].psmvalid[frq]) {
                const ssat_t *ss=rtk->ssat+sat-1;
                double raw=varerr(sat,sys,azel[1+i*2],snr_real,
                    signal_noise_index_ppp(sys,obs[i].code[frq],frq),1,opt)*
                    signal_weight_ppp(sys,obs[i].code[frq],1,opt);
                double effective=ss->psm_var[frq];
                if(statmodel&4)effective+=2.0*ss->psm_past_cov[frq];
                if(isfinite(effective)&&effective>raw) {
                    double scale=1.0;
                    if(sys==SYS_CMP)scale*=bds_code_factor;
                    if(!sat_has_phase&&!rtk->ppp_code_warmup)scale*=code_only_factor;
                    var[nv]+=(effective-raw)*scale;
                }
            }
            /* Persistent penalties affect R, never erase the physical
             * standardized residual used for detection/hard rejection. */
            if (code && rtk->ssat[sat - 1].ppp_cmc_code_weak[frq]) var[nv]*=4.0;
            var[nv]*=persistent_variance_factor_ppp(opt,
                &rtk->ssat[sat-1].ppp_res_quality[frq][code],code);
            trace(4,"$SIG_WEIGHT,%s,sat=%d,sig=%s,F%d,type=%s,noise_band=%d,cal_var=%.3f,persistent_var=%d,physical_sigma=%.4f\n",
                str,sat,code2obs(obs[i].code[frq]),frq+1,code?"CODE":"PHASE",
                signal_noise_index_ppp(sys,obs[i].code[frq],frq)+1,
                signal_weight_ppp(sys,obs[i].code[frq],code,opt),
                persistent_variance_factor_ppp(opt,&rtk->ssat[sat-1].ppp_res_quality[frq][code],code),sigma0);

            if (z0 > THRES_ROBUST) {
                penalty = SQR(z0 / THRES_ROBUST);
                if (penalty > 25.0) penalty = 25.0;
                var[nv] *= penalty;

                trace(3,
                    "robust down-weight: sat=%d %s%d z0=%.2f penalty=%.2f sigma0=%.3f sigma=%.3f\n",
                    sat, code ? "P" : "L", frq + 1, z0, penalty,
                    sigma0, sqrt(var[nv]));
            }

            trace(3, "%s sat=%2d %s%d res=%9.4f sig0=%8.4f sig=%8.4f z0=%6.2f el=%4.1f snr=%4.1f ifb=%8.3f\n",
                str, sat, code ? "P" : "L", frq + 1, res, sigma0,
                sqrt(var[nv]), z0, azel[1 + i * 2] * R2D, snr_real, ifb_used);

            /*
             * Post-fit hard outlier detection uses the uninflated sigma.
             * Candidate ranking also uses standardized residual z0 instead of
             * raw metres, otherwise code measurements always dominate phase.
             */
            if (post && z0 > THRES_REJECT && ne < MAXOBS * 2 * NFREQ) {
                obsi[ne] = i;
                frqi[ne] = j;
                ve[ne] = res;
                ze[ne] = z0;
                se[ne] = sigma0;
                snre[ne] = snr_real;
                ele[ne] = azel[1 + i * 2] * R2D;
                ifbe[ne] = ifb_used;
                ne++;
            }

            /* vsat/lock/outc describe accepted carrier tracking only. A
             * code-only filter update is labelled separately below. */
            if (code == 0) rtk->ssat[sat - 1].vsat[frq] = 1;
            nv++;
        }
    }

    /* IONEX soft constraints -------------------------------------------------
     * For IONOOPT_EST, GNSS code/phase observations estimate one vertical L1
     * ionosphere state per satellite. IONEX is added as a loose pseudo-
     * observation Iv_IONEX - Iv_state, rather than being directly subtracted
     * from every measurement. This preserves local information in the phone
     * observations while preventing the ionosphere state from drifting freely.
     *
     * Constraints are added only in the pre-fit filter equation. Post-fit
     * rejection remains reserved for real GNSS measurements. */
    /* IONEX 伪观测是一颗卫星一条垂直状态约束，不是每个频率一条。
     * 设计矩阵仅该星 II 列为 1；单位与状态一致：v=Iv_product-Iv_state。
     * STATMOD 位2使初始化已用过的本历元先验不再加入；后续按 dt 分配信息。
     * product_use 只记试算里用了哪些先验，成功提交 xp/Pp 后才更新使用时间。 */
    if (!post && opt->ionoopt == IONOOPT_EST && v && H) {
        unsigned char used[MAXSAT] = { 0 };
        int logged_interval = 0;

        for (i = 0; i < n && i < MAXOBS; i++) {
            double ion_prior, var_prior, sigma_prior, zc, penalty = 1.0;
            double interval, dt, temporal_factor, information=1.0;
            int ii;

            sat = obs[i].sat;
            if (sat < 1 || sat > MAXSAT || used[sat - 1] || exc[i]) continue;
            used[sat - 1] = 1;
            if (azel[1 + i * 2] < opt->elmin) continue;

            ii = II(sat, opt);
            if (ii < 0 || ii >= nx || x[ii] == 0.0) continue;

            if (!ionex_vertical_prior(obs[i].time, nav, pos, azel + i * 2,
                opt, &ion_prior, &var_prior)) {
                continue;
            }
            interval = ion_constraint_interval(opt);
            if (statmodel&2) {
                information=atmosphere_information_ppp(obs[i].time,
                    rtk->ssat[sat-1].ppp_ion_prior_time,interval);
                if(information<=0.0)continue;
            }

            if (nv >= MAXOBS * 2 * NFREQ + MAXSAT + 3) break;

            for (k = 0; k < nx; k++) H[k + nx * nv] = 0.0;
            H[ii + nx * nv] = 1.0;
            v[nv] = ion_prior - x[ii];

            /* If the local GNSS solution strongly disagrees with the global GIM,
             * weaken the external constraint instead of forcing the state back
             * to IONEX. This is essential for disturbed/local ionosphere. */
            sigma_prior = sqrt(MAX(var_prior, 1E-12));
            zc = fabs(v[nv]) / sigma_prior;
            if (zc > ION_CONSTR_ROBUST) {
                penalty = SQR(zc / ION_CONSTR_ROBUST);
                if (penalty > ION_CONSTR_MAXPEN) penalty = ION_CONSTR_MAXPEN;
                var_prior *= penalty;
            }
            dt = (statmodel&2)&&rtk->ssat[sat-1].ppp_ion_prior_time.time ?
                timediff(obs[i].time,rtk->ssat[sat-1].ppp_ion_prior_time):MAX(0.1,fabs(rtk->tt));
            temporal_factor = (statmodel&2)?1.0/information:MAX(1.0,interval/dt);
            var[nv] = var_prior * temporal_factor;
            if(product_use)product_use[3+sat-1]=1;
            trace(3,"$ATM_STAT,%s,type=IONEX,sat=%d,dt=%.6f,information=%.9f,base_sigma=%.6f,robust_penalty=%.6f,effective_variance=%.9f,mode=%d\n",
                str,sat,dt,1.0/temporal_factor,sigma_prior,penalty,var[nv],statmodel);

            if (!logged_interval) {
                trace(3, "$IONCON_TIME,%s,dt=%.3f,interval=%.0f,factor=%.1f,sig=%.3f,effsig=%.3f\n",
                    str, dt, interval, temporal_factor, sigma_prior,
                    sqrt(var[nv]));
                logged_interval = 1;
            }

            if (zc > ION_CONSTR_ROBUST) {
                satno2id(sat, sid);
                trace(2, "$IONCON_WARN,%s,%s,Iest=%.4f,Iionex=%.4f,res=%.4f,sig=%.3f,effsig=%.3f,interval=%.0f,z=%.2f,pen=%.2f\n",
                    str, sid, x[ii], ion_prior, v[nv], sigma_prior,
                    sqrt(var[nv]), interval, zc, penalty);
            }
            trace(4, "$IONCON,%s,sat=%d,Iest=%.4f,Iionex=%.4f,res=%.4f,sig=%.3f,effsig=%.3f,interval=%.0f\n",
                str, sat, x[ii], ion_prior, v[nv], sigma_prior,
                sqrt(var[nv]), interval);
            nv++;
            ionex_constraints++;
        }
    }

    /* VMF3 ZWD soft constraint ----------------------------------------------
     * ZHD remains a deterministic VMF3 model correction in trop_model_prec().
     * The estimated state is total ZTD, so constraining
     *   ZWD = ZTD_state - ZHD_VMF3
     * is equivalent to a pseudo-observation of ZHD_VMF3 + ZWD_VMF3. The
     * product is deliberately weak and time-decorrelated; large disagreement
     * inflates its variance instead of forcing the PPP state to the forecast. */
    /* VMF3 约束的是未知湿延迟，而内部仍沿用总 ZTD 状态：
     * v=ZWD_product-(ZTD_state-ZHD_product)，所以 H[IT]=1。
     * 这与约束 ZTD 到 ZHD+ZWD 在固定产品 ZHD 时等价；不另建一个独立 ZWD。
     * VMF3ZWDSIG<=0 或产品时间/空间不可用时不追加此行，不强行改写状态。
     * 有效 R 还乘重复信息与抗差因子；产品偏离大时软化，不当作无误差真值。 */
    if (!post && v && H &&
        (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG)) {
        double sig_zwd = vmf3_zwd_constraint_sigma(opt);
        double ah, aw, zhd_vmf, zwd_vmf, dt, interval, zc, penalty = 1.0;
        double information=1.0;
        int it = IT(opt);
        interval=vmf3_zwd_constraint_interval(opt);
        if(statmodel&2)information=atmosphere_information_ppp(obs[0].time,
            rtk->ppp_zwd_prior_time,interval);

        if (information>0.0 && sig_zwd > 0.0 && it >= 0 && it < nx && x[it] != 0.0 &&
            vmf3_grid_interp(obs[0].time, pos, &ah, &aw,
                &zhd_vmf, &zwd_vmf) &&
            nv < MAXOBS * 2 * NFREQ + MAXSAT + 3) {
            for (k = 0; k < nx; k++) H[k + nx * nv] = 0.0;
            H[it + nx * nv] = 1.0;
            v[nv] = zwd_vmf - (x[it] - zhd_vmf);

            zc = fabs(v[nv]) / sig_zwd;
            if (zc > VMF3_ZWD_ROBUST) {
                penalty = SQR(zc / VMF3_ZWD_ROBUST);
                penalty = MIN(VMF3_ZWD_MAXPEN, penalty);
            }
            dt = (statmodel&2)&&rtk->ppp_zwd_prior_time.time?
                timediff(obs[0].time,rtk->ppp_zwd_prior_time):MAX(0.1,fabs(rtk->tt));
            if(statmodel&2)var[nv] = SQR(sig_zwd) / information * penalty;
            else var[nv] = SQR(sig_zwd) * interval / dt * penalty;
            if(product_use)product_use[2]=1;
            trace(3,"$ATM_STAT,%s,type=VMF3_ZWD,sat=0,dt=%.6f,information=%.9f,base_sigma=%.6f,robust_penalty=%.6f,effective_variance=%.9f,mode=%d\n",
                str,dt,(statmodel&2)?information:dt/interval,sig_zwd,penalty,var[nv],statmodel);

            trace(zc > VMF3_ZWD_ROBUST ? 2 : 4,
                "$VMF3_ZWD_CONSTR,%s,ZWDest=%.4f,ZWDvmf=%.4f,res=%.4f,"
                "sig=%.3f,interval=%.0f,effsig=%.3f,z=%.2f,pen=%.2f\n",
                str, x[it] - zhd_vmf, zwd_vmf, v[nv], sig_zwd,
                interval, sqrt(var[nv]), zc, penalty);
            nv++;
            vmf3_zwd_constraints++;
        }
    }

    if (!post && product_use) {
        product_use[0] = ionex_constraints;
        product_use[1] = vmf3_zwd_constraints;
    }

    /* 一次编辑只屏蔽验后 z0 最大的单条观测，并返回 stat=0 请求重算。
     * 伪距噪声比相位大，不能直接按残差米数排序，否则总会优先删码。
     * 不把 mexc 换成 exc：同星其他正常频率/另一种观测应有机会保留。 */
    if (post && ne > 0) {
        rej = 0;
        for (j = 1; j < ne; j++) {
            if (ze[rej] >= ze[j]) continue;
            rej = j;
        }

        maxobs = obsi[rej];
        maxfrq = frqi[rej];
        sat = obs[maxobs].sat;
        frq = maxfrq / 2;
        code = maxfrq % 2;

        if (mexc) {
            int midx = maxobs * 2 * NFREQ + maxfrq;
            if (midx >= 0 && midx < MAXOBS * 2 * NFREQ) mexc[midx] = 1;
        }

        if (frq >= 0 && frq < NFREQ) rtk->ssat[sat - 1].rejc[frq]++;

        satno2id(sat, sid);
        trace(2,
            "$PPP_REJECT,%s,iter=%d,sat=%s,F%d,%s,res=%.4f,sigma=%.4f,z=%.2f,snr=%.1f,el=%.1f,ifb=%.4f,candidates=%d\n",
            str, post, sid, frq + 1, code ? "CODE" : "PHASE",
            ve[rej], se[rej], ze[rej], snre[rej], ele[rej], ifbe[rej], ne);

        /*
         * Request another EKF iteration after excluding only this measurement.
         * Do NOT set exc[maxobs]: one bad E5b/L5 observable must not delete the
         * whole satellite and its otherwise good frequencies.
         */
        stat = 0;
    }

    /* 构造实际送 filter() 的对角协方差矩阵；var[] 是方差不是权，
     * 外部伪观测行也在同一个 v/H/R 中，矩阵大小按最终 nv 而非原始 n。 */
    if (R) {
        for (j = 0; j < nv; j++) {
            for (i = 0; i < nv; i++) R[i + j * nv] = 0.0;
        }
        for (i = 0; i < nv; i++) R[i + i * nv] = var[i];
    }

    return post ? stat : nv;
}

/* Count rejections across epochs, never across editing iterations. Only a
 * repeatedly rejected carrier is blocked; the corresponding code survives. */
/* -PPPQUAR=N,SEC：跨不同输入历元连续拒绝 N 次才隔离该载波。
 * 同一历元 16 次残差编辑最多只计一次；未建模/没测试不当作好观测奖励。
 * 重复坏弧逐级延长隔离（最长 600 s），连续有效跟踪 120 s 才清退避等级。
 * 进入隔离清掉模糊度及相关；同频码/其他频率仍可用，重新进入时重建相位弧。
 * 此处是严重异常阻断，与 SIGQC 的持续温和降权不同。 */
static void update_phase_quarantine_ppp(rtk_t *rtk, const obsd_t *obs, int n,
    const unsigned char *mexc, const ppp_diag_epoch_t *diag)
{
    const char *p = strstr(rtk->opt.pppopt, "-PPPQUAR=");
    int i, f, sat, j, trigger = 0;
    double hold = 0.0;
    char sid[8] = "", str[40];

    if (!p || sscanf(p, "-PPPQUAR=%d,%lf", &trigger, &hold) != 2 ||
        trigger < 2 || hold <= 0.0) return;
    trigger = MIN(trigger, 20);
    hold = MIN(hold, 120.0);
    time2str(obs[0].time, str, 2);

    for (i = 0; i < n && i < MAXOBS; i++) {
        ssat_t *ssat;
        double actual_hold;
        sat = obs[i].sat;
        if (sat < 1 || sat > MAXSAT) continue;
        ssat = &rtk->ssat[sat - 1];
        for (f = 0; f < NF(&rtk->opt) && f < NFREQ; f++) {
            if (obs[i].L[f] == 0.0 ||
                ppp_phase_quarantined(ssat, f, obs[i].time)) continue;
            j = i * 2 * NFREQ + 2 * f;
            if (!mexc[j]) {
                /* Absent/unmodelled/low-elevation data are not accepted
                 * evidence. Do not clear rejection history or reward them. */
                if (!diag || !diag->obs[i][f].phase_ok) continue;
                ssat->ppp_phase_reject_streak[f] = 0;
                if (ssat->vsat[f]) {
                    if (!ssat->ppp_phase_good_since[f].time) {
                        ssat->ppp_phase_good_since[f] = obs[i].time;
                    }
                    else if (timediff(obs[i].time,
                             ssat->ppp_phase_good_since[f]) >= 120.0) {
                        ssat->ppp_phase_quar_count[f] = 0;
                    }
                }
                continue;
            }
            ssat->ppp_phase_good_since[f].time = 0;
            if (!ssat->ppp_phase_reject_time[f].time ||
                timediff(obs[i].time, ssat->ppp_phase_reject_time[f]) > 5.0) {
                ssat->ppp_phase_reject_streak[f] = 0;
            }
            ssat->ppp_phase_reject_time[f] = obs[i].time;
            if (ssat->ppp_phase_reject_streak[f]<255)
                ssat->ppp_phase_reject_streak[f]++;
            if (ssat->ppp_phase_reject_streak[f] < trigger) continue;

            ssat->ppp_phase_reject_streak[f] = 0;
            if (ssat->ppp_phase_quar_count[f] < 7) {
                ssat->ppp_phase_quar_count[f]++;
            }
            /* Repeatedly rejected arcs must not re-enter every two minutes.
             * Accepted tracking for 120 s still clears this backoff. */
            actual_hold = MIN(600.0, hold * (1 <<
                              (ssat->ppp_phase_quar_count[f] - 1)));
            ssat->ppp_phase_block_until[f] = timeadd(obs[i].time, actual_hold);
            j = IB(sat, f, &rtk->opt);
            if (j >= 0 && j < rtk->nx) initx(rtk, 0.0, 0.0, j);
            satno2id(sat, sid);
            trace(2, "$PPP_QUAR,%s,sat=%s,sig=%s,F%d,hold=%.1f,level=%d,reason=%d_CONSEC_REJECTS\n",
                str, sid, code2obs(obs[i].code[f]), f + 1, actual_hold,
                ssat->ppp_phase_quar_count[f], trigger);
        }
    }
}

typedef struct {
    char sys;
    uint8_t signal;
    int np, nl, nd, resets;
    double sum_p, sum_p2, sum_l, sum_l2, sum_d, sum_d2;
} ppp_diag_group_t;

/* Emit one ZTD record per epoch, then accepted residuals grouped by actual
 * constellation/signal. CMC is P-L-2*I on the corrected but UNSMOOTHED
 * measurements; its arc-relative change is diagnostic, not a pure hardware
 * bias estimate (residual ionosphere/multipath can contribute). */
/* 每历元诊断输出（不用于真值反馈或坐标改正）：
 * $PPP_DIAG_ZTD：ZTD/ZHD/ZWD_est/产品ZWD/标准差，vmf=1 才表示本次插值有效。
 * $PPP_DIAG_SIG：按真实星座+code 汇总最终残差的数量/均值/RMS，不按槽混组。
 * $PPP_DIAG_CMC：未平滑改正码相组合的连续弧变化，周跳/换信号/半周变化
 * 和中断重置基准。绝对 CMC 包含模糊度/RCB，不能把漂移等同纯码偏差。
 * mean 是有符号平均，RMS=sqrt(mean(res^2))，不是去均值后的标准差。
 * ppp_ok 与 Q 一起读；滤波失败仍可报告预测 ZTD，但不能叫它有效定位结果。 */
static void ppp_diag_emit(rtk_t* rtk, const obsd_t* obs, int n,
    const ppp_diag_epoch_t* diag, int ppp_ok)
{
    ppp_diag_group_t groups[MAXOBS * NFREQ] = { 0 };
    const prcopt_t* opt = &rtk->opt;
    const double zenith[2] = { 0.0, PI / 2.0 };
    double pos[3], zhd = 0.0, zwd_vmf = 0.0, ah, aw, ztd, sigma;
    double step, drift, gap, allowed_gap;
    char str[40], sid[8];
    int i, f, g, ng = 0, sat, vmf_ok = 0, reset;

    if (!obs || n <= 0) return;
    time2str(obs[0].time, str, 2);
    if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
        ztd = rtk->x[IT(opt)];
        sigma = sqrt(MAX(rtk->P[IT(opt) + IT(opt) * rtk->nx], 0.0));
        if (norm(rtk->x, 3) > 1E6) {
            ecef2pos(rtk->x, pos);
            vmf_ok = vmf3_grid_interp(obs[0].time, pos, &ah, &aw,
                &zhd, &zwd_vmf);
            if (!vmf_ok) zhd = tropmodel(obs[0].time, pos, zenith, 0.0);
        }
        trace(2,
            "$PPP_DIAG_ZTD,time=%s,ppp=%d,Q=%d,vmf=%d,ZTD=%.4f,ZHD=%.4f,ZWDest=%.4f,ZWDvmf=%.4f,dZTDvmf=%.4f,sig=%.4f\n",
            str, ppp_ok, rtk->sol.stat, vmf_ok, ztd, zhd, ztd - zhd, zwd_vmf,
            vmf_ok ? ztd - zhd - zwd_vmf : 0.0, sigma);
    }
    if (!diag) return;

    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        if (sat < 1 || sat > MAXSAT) continue;
        satno2id(sat, sid);
        for (f = 0; f < NF(opt) && f < NFREQ; f++) {
            const ppp_diag_signal_t* d = &diag->obs[i][f];
            ssat_t* ss = &rtk->ssat[sat - 1];
            ppp_diag_group_t* group;

            if (!d->signal || d->signal > MAXCODE ||
                (!d->phase_ok && !d->code_ok)) continue;
            for (g = 0; g < ng; g++) {
                if (groups[g].sys == sid[0] && groups[g].signal == d->signal) break;
            }
            if (g == ng) {
                if (ng >= MAXOBS * NFREQ) continue;
                groups[ng].sys = sid[0];
                groups[ng].signal = d->signal;
                ng++;
            }
            group = &groups[g];
            if (d->code_ok) {
                group->np++;
                group->sum_p += d->code_res;
                group->sum_p2 += SQR(d->code_res);
            }
            if (d->phase_ok) {
                group->nl++;
                group->sum_l += d->phase_res;
                group->sum_l2 += SQR(d->phase_res);
            }
            if (!d->cmc_ok || !d->phase_ok || !d->code_ok) continue;

            gap = ss->ppp_diag_cmc_time[f].time ?
                timediff(obs[i].time, ss->ppp_diag_cmc_time[f]) : 0.0;
            allowed_gap = MAX(3.0, 2.5 * fabs(rtk->tt));
            reset = !ss->ppp_diag_cmc_time[f].time ||
                gap <= 0.0 || gap > allowed_gap ||
                ss->ppp_diag_cmc_code[f] != d->signal ||
                (ss->slip[f] & LLI_SLIP) ||
                (obs[i].LLI[f] & LLI_SLIP) ||
                (ss->ppp_diag_cmc_half_invalid[f] !=
                    ((obs[i].LLI[f] & LLI_HALFC) != 0));
            if (reset) {
                ss->ppp_diag_cmc_start[f] = obs[i].time;
                ss->ppp_diag_cmc_first[f] = d->cmc;
                step = drift = 0.0;
                group->resets++;
            }
            else {
                step = d->cmc - ss->ppp_diag_cmc_last[f];
                drift = d->cmc - ss->ppp_diag_cmc_first[f];
                group->nd++;
                group->sum_d += drift;
                group->sum_d2 += SQR(drift);
            }
            ss->ppp_diag_cmc_time[f] = obs[i].time;
            ss->ppp_diag_cmc_last[f] = d->cmc;
            ss->ppp_diag_cmc_code[f] = d->signal;
            ss->ppp_diag_cmc_half_invalid[f] =
                (obs[i].LLI[f] & LLI_HALFC) != 0;
            trace(3,
                "$PPP_DIAG_CMC,time=%s,sat=%s,sig=%s,F%d,arc_s=%.0f,reset=%d,CMC=%.3f,iono=%.3f,step=%.3f,drift=%.3f,snr=%.1f,el=%.1f\n",
                str, sid, code2obs(d->signal), f + 1,
                timediff(obs[i].time, ss->ppp_diag_cmc_start[f]), reset,
                d->cmc, d->iono, step, drift, d->snr, d->el);
        }
    }
    for (g = 0; g < ng; g++) {
        const ppp_diag_group_t* group = &groups[g];
        trace(2,
            "$PPP_DIAG_SIG,time=%s,sys=%c,sig=%s,nP=%d,meanP=%.4f,rmsP=%.4f,nL=%d,meanL=%.4f,rmsL=%.4f,nDrift=%d,meanDrift=%.4f,rmsDrift=%.4f,resets=%d\n",
            str, group->sys, code2obs(group->signal),
            group->np,
            group->np ? group->sum_p / group->np : 0.0,
            group->np ? sqrt(group->sum_p2 / group->np) : 0.0,
            group->nl,
            group->nl ? group->sum_l / group->nl : 0.0,
            group->nl ? sqrt(group->sum_l2 / group->nl) : 0.0,
            group->nd,
            group->nd ? group->sum_d / group->nd : 0.0,
            group->nd ? sqrt(group->sum_d2 / group->nd) : 0.0,
            group->resets);
    }
}

/* rtkinit 等调用方通过此接口获知 PPP 状态总数；索引宏和分配大小要一致。
 * 这是分配维数，不是每历元的观测数/活跃状态数，也不是卫星数量。 */
extern int pppnx(const prcopt_t* opt)
{
    return NX(opt);
}

/* 向状态输出/诊断层公开模糊度索引，f 为零基槽；非法输入返回 -1。
 * 调用方不可自己复制固定偏移，nf/ionoopt/dynamics 会改变前面状态块大小。 */
extern int pppambidx(int sat, int f, const prcopt_t* opt)
{
    if (!opt || sat < 1 || sat > MAXSAT || f < 0 || f >= NF(opt)) return -1;
    return IB(sat, f, opt);
}
/* =============================================================================
 * 功能：把已接受滤波状态复制到外部解 rtk->sol，并更新有效跟踪计数。
 * 任一载波频槽通过编辑即计该星有效，一颗星只计一次；不是只数首频 L1。
 * 少于 4 颗载波有效卫星先设 Q0；pppos 随后可对足够码更新标记为 Q5 降级。
 * 位置 rr 为 ECEF，qr 存位置协方差（不是 std）；dtr 由米除光速转秒。
 * Q6 仅表示载波 PPP，Q1 才是固定解；Q 值不直接表示绝对误差或收敛阈值。
 * ============================================================================= */
static void update_stat(rtk_t* rtk, const obsd_t* obs, int n, int stat)
{
    const prcopt_t* opt = &rtk->opt;
    int i, j;

    /* test # of valid satellites
     * Count a satellite once if ANY carrier frequency survives PPP editing.
     * Phones may have usable E5/L5 while L1 is temporarily absent; counting only
     * j==0 can incorrectly turn an otherwise valid PPP epoch into SOLQ_NONE.
     */
    rtk->sol.ns = 0;
    for (i = 0; i < n && i < MAXOBS; i++) {
        int sat_valid = 0;

        for (j = 0; j < opt->nf; j++) {
            if (!rtk->ssat[obs[i].sat - 1].vsat[j]) continue;
            rtk->ssat[obs[i].sat - 1].lock[j]++;
            rtk->ssat[obs[i].sat - 1].outc[j] = 0;
            sat_valid = 1;
        }
        if (sat_valid) rtk->sol.ns++;
    }
    rtk->sol.stat = rtk->sol.ns < MIN_NSAT_SOL ? SOLQ_NONE : stat;

    if (rtk->sol.stat == SOLQ_FIX) {
        for (i = 0; i < 3; i++) {
            rtk->sol.rr[i] = rtk->xa[i];
            rtk->sol.qr[i] = (float)rtk->Pa[i + i * rtk->na];
        }
        rtk->sol.qr[3] = (float)rtk->Pa[1];
        rtk->sol.qr[4] = (float)rtk->Pa[1 + 2 * rtk->na];
        rtk->sol.qr[5] = (float)rtk->Pa[2];
    }
    else {
        for (i = 0; i < 3; i++) {
            rtk->sol.rr[i] = rtk->x[i];
            rtk->sol.qr[i] = (float)rtk->P[i + i * rtk->nx];
        }
        rtk->sol.qr[3] = (float)rtk->P[1];
        rtk->sol.qr[4] = (float)rtk->P[2 + rtk->nx];
        rtk->sol.qr[5] = (float)rtk->P[2];

        if (rtk->opt.dynamics) { /* velocity and covariance */
            for (i = 3; i < 6; i++) {
                rtk->sol.rr[i] = rtk->x[i];
                rtk->sol.qv[i - 3] = (float)rtk->P[i + i * rtk->nx];
            }
            rtk->sol.qv[3] = (float)rtk->P[4 + 3 * rtk->nx];
            rtk->sol.qv[4] = (float)rtk->P[5 + 4 * rtk->nx];
            rtk->sol.qv[5] = (float)rtk->P[5 + 3 * rtk->nx];
        }
    }
    rtk->sol.dtr[0] = rtk->x[IC(0, opt)] / CLIGHT;
    rtk->sol.dtr[1] = rtk->x[IC(1, opt)] / CLIGHT;
    rtk->sol.dtr[2] = rtk->x[IC(2, opt)] / CLIGHT;
    rtk->sol.dtr[3] = rtk->x[IC(3, opt)] / CLIGHT;
    rtk->sol.dtr[4] = rtk->x[IC(4, opt)] / CLIGHT;
    rtk->sol.dtr[5] = 0.0; /* QZS shares GPST */

    for (i = 0; i < n && i < MAXOBS; i++) for (j = 0; j < opt->nf; j++) {
        rtk->ssat[obs[i].sat - 1].snr_rover[j] = obs[i].SNR[j];
        rtk->ssat[obs[i].sat - 1].snr_base[j] = 0;
    }
    for (i = 0; i < MAXSAT; i++) for (j = 0; j < opt->nf; j++) {
        if (rtk->ssat[i].slip[j] & (LLI_SLIP | LLI_HALFC)) rtk->ssat[i].slipc[j]++;
        if (rtk->ssat[i].fix[j] == 2 && stat != SOLQ_FIX) rtk->ssat[i].fix[j] = 1;
    }
}
/* 固定解连续保持测试：仅 modear=FIXHOLD 时参与，新模糊度组合会清 nfix。
 * 本文件保留 AR 接口不表示默认启用；当前 modear=OFF 时不做整周固定。
 * 本次没有修改 ppp_ar 或固定/保持策略。 */
static int test_hold_amb(rtk_t* rtk)
{
    int i, j, stat = 0;

    /* no fix-and-hold mode */
    if (rtk->opt.modear != ARMODE_FIXHOLD) return 0;

    /* reset # of continuous fixed if new ambiguity introduced */
    for (i = 0; i < MAXSAT; i++) {
        if (rtk->ssat[i].fix[0] != 2 && rtk->ssat[i].fix[1] != 2) continue;
        for (j = 0; j < MAXSAT; j++) {
            if (rtk->ssat[j].fix[0] != 2 && rtk->ssat[j].fix[1] != 2) continue;
            if (!rtk->ambc[j].flags[i] || !rtk->ambc[i].flags[j]) stat = 1;
            rtk->ambc[j].flags[i] = rtk->ambc[i].flags[j] = 1;
        }
    }
    if (stat) {
        rtk->nfix = 0;
        return 0;
    }
    /* test # of continuous fixed */
    return ++rtk->nfix >= rtk->opt.minfix;
}
/* =============================================================================
 * 功能：对一个有效输入历元执行完整 PPP，rtk_t 状态跨历元持续存在。
 * 当前实现顺序：信号历史准备 -> 卫星位置/钟差 -> 预检 -> 平滑/启动切换
 * -> 状态预测 -> 产品诊断/地影/潮汐 -> 测量更新与残差编辑 -> 解状态/日志。
 * 调用前 obs/n/nav 必须有效；本函数使用 obs[0]，不是空事件的处理入口。
 * 卫星产品在外部启动时加载，本函数不每秒读产品、更不每秒 rtkinit。
 * 多次循环是同一历元的粗差编辑试算：每次从相同预测 x/P 出发，只提交
 * 最终验后通过的 xp/Pp 一次；否则会重复把同一秒数据当作新观测强化。
 * ============================================================================= */
extern void pppos(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
    const prcopt_t* opt = &rtk->opt;
    double* rs, * dts, * var, * v, * H, * R, * azel, * xp, * Pp, dr[3] = { 0 }, std[3];
    char str[40];
    const char* diag_opt;
    ppp_diag_epoch_t diag, evidence;
    int i, j, nv, info, nsmooth, warmup, diag_enabled = 0,
        product_use[MAXSAT+3] = {0},
        svh[MAXOBS], exc[MAXOBS] = { 0 }, stat = SOLQ_SINGLE;
    unsigned char mexc[MAXOBS * 2 * NFREQ] = { 0 };
    memset(&evidence,0,sizeof(evidence));

    time2str(obs[0].time, str, 2);
    if ((diag_opt = strstr(opt->pppopt, "-PPPDIAG="))) {
        sscanf(diag_opt, "-PPPDIAG=%d", &diag_enabled);
    }
    trace(3, "pppos   : time=%s nx=%d n=%d\n", str, rtk->nx, n);
    trace(3,"$STAT_MODEL,%s,mode=%d,doppler_sigma_hz=%.3f,ion_interval=%.1f,zwd_interval=%.1f\n",
        str,statistical_model_ppp(opt),opt->err[4],ion_constraint_interval(opt),
        vmf3_zwd_constraint_interval(opt));
    rs = mat(6, n); dts = mat(2, n); var = mat(1, n); azel = zeros(2, n);

    for (i = 0; i < MAXSAT; i++) for (j = 0; j < opt->nf; j++) rtk->ssat[i].fix[j] = 0;
    for (i = 0; i < n && i < MAXOBS; i++) for (j = 0; j < opt->nf; j++) {
        rtk->ssat[obs[i].sat - 1].snr_rover[j] = obs[i].SNR[j];
        rtk->ssat[obs[i].sat - 1].snr_base[j] = 0;
    }

    /* 第 1 步：同槽换 signal 时先断开历史；然后根据本历元 nav/观测取
     * 卫星位置速度 rs（6*n）、钟差/漂移 dts（2*n）、产品方差和健康标志。
     * satposs 的结果也是启动粗差预检所需的输入，故它在 udstate 之前。 */
    prepare_signal_quality_ppp(rtk, obs, n);
    /* Ephemerides depend only on the input epoch/nav, not predicted states. */
    satposs(obs[0].time,obs,n,nav,opt->sateph,rs,dts,var,svh);

    /* Update the causal code smoother before ambiguity initialization and
     * residual formation. This runs once per epoch, never inside the
     * post-fit editing iterations. */
    /* Causal raw-observation quality control runs before smoothing and state
     * prediction, so rejected code and detected slips never seed PPP states. */
    /* 第 2 步：只依赖当前/历史数据的粗差和周跳预检；先标坏再建立平滑。
     * LLI 等 slip 标记通知重建模糊度；code_bad 阻止坏码参与初值/滤波。 */
    if (pppopt_number(opt, "-PREPROC=", 0.0, 0.0, 1.0) >= 0.5) {
        precheck_obs_ppp(rtk, obs, n, nav);
        precheck_startup_code_ppp(rtk,obs,n,nav,rs,dts,var,svh);
    }
    /* 第 3 步：平滑及可选码阶段选择，各执行一次；严禁移进残差编辑循环，
     * 否则同一 epoch 会重复平滑和增加历史计数，fast/realtime 也可能不等价。 */
    nsmooth = uddoppsm_ppp(rtk, obs, n, nav, rs);
    warmup = update_doppsm_warmup(rtk, obs[0].time);
    if (doppsm_factor(opt) > 0.0) {
        trace(3, "$DOPPSM,time=%s,alpha=%.3f,propagated=%d,stage=%s\n", str,
            doppsm_factor(opt), nsmooth, warmup ? "CODE_DOPPLER" : "PPP_PHASE");
    }

    /* 第 4 步：形成本历元的预测状态 rtk->x/P。
     * 模糊度最后初始化，使用已平滑码及一致的电离层/RCB 模型。 */
    udstate_ppp(rtk, obs, n, nav);

    /* Satellite positions/clocks were computed before startup code checks. */

    /* Product coverage is reported at the observation epoch, not inferred
     * from the mere presence of input files. Satellite-specific interpolation
     * can still fail; the valid satellite count is reported separately. */
    /* 第 5 步：产品诊断只输出，不通过“文件存在”认定产品有效。
     * span=IN_RANGE 是总体时间覆盖；某卫星缺记录仍可能插值失败。
     * code_osb_in_time 必须匹配 epoch+sat+signal，而非只数 BIA 文件记录。
     * $PPP_PRODUCT_USE 稍后才报告实际加入的 IONEX/VMF3 伪观测数量。 */
    if (opt->sateph == EPHOPT_PREC) {
        int oi, fi, osb_ok = 0, osb_expired = 0, osb_missing = 0;
        int valid_satpos = 0, vmf3_model = 0;
        const char *orbit = "MISSING", *clock = "MISSING";
        const char *ionex = "MISSING", *vmf3 = "MISSING";
        gtime_t epoch = obs[0].time;
        double osb, model_pos[3], ah, aw, zhd, zwd;

        if (nav->ne > 0 && nav->peph) {
            orbit = timediff(epoch, nav->peph[0].time) >= 0.0 &&
                    timediff(epoch, nav->peph[nav->ne - 1].time) <= 0.0 ?
                    "IN_RANGE" : "OUT_OF_RANGE";
        }
        if (nav->nc > 0 && nav->pclk) {
            clock = timediff(epoch, nav->pclk[0].time) >= 0.0 &&
                    timediff(epoch, nav->pclk[nav->nc - 1].time) <= 0.0 ?
                    "IN_RANGE" : "OUT_OF_RANGE";
        }
        if (nav->nt > 0 && nav->tec) {
            ionex = timediff(epoch, nav->tec[0].time) >= 0.0 &&
                    timediff(epoch, nav->tec[nav->nt - 1].time) <= 0.0 ?
                    "IN_RANGE" : "OUT_OF_RANGE";
        }
        if (vmf3_grid[0].valid && vmf3_grid[1].valid && vmf3_orog_valid) {
            vmf3 = timediff(epoch, vmf3_grid[0].time) >= 0.0 &&
                   timediff(epoch, vmf3_grid[1].time) <= 0.0 ?
                   "IN_RANGE" : "OUT_OF_RANGE";
        }
        if (norm(rtk->x, 3) > 0.0) {
            ecef2pos(rtk->x, model_pos);
            vmf3_model = vmf3_grid_interp(epoch, model_pos, &ah, &aw,
                                          &zhd, &zwd);
        }
        for (oi = 0; oi < n; oi++) {
            if (norm(rs + oi * 6, 3) > 0.0) valid_satpos++;
            for (fi = 0; fi < opt->nf && fi < NFREQ; fi++) {
                int status;
                if (obs[oi].P[fi] == 0.0 || obs[oi].code[fi] == CODE_NONE) continue;
                status = codeosb_at(nav, epoch, obs[oi].sat,
                                    obs[oi].code[fi], &osb);
                if (status > 0) osb_ok++;
                else if (status < 0) osb_expired++;
                else osb_missing++;
            }
        }
        trace(2, "$PPP_PRODUCT,%s,orbit_span=%s,clock_span=%s,"
              "ionex_span=%s,vmf3_span=%s,vmf3_model=%s,"
              "satpos=%d/%d,code_osb_in_time=%d,code_osb_expired=%d,"
              "code_osb_absent=%d\n", str, orbit, clock, ionex, vmf3,
              vmf3_model ? "VMF3" : "FALLBACK",
              valid_satpos, n, osb_ok, osb_expired, osb_missing);
    }

    /* 第 6 步：按开关处理地影和潮汐。潮汐接口使用 UTC，观测仍保持 GPST。
     * dr 仅用于几何位置，不能把这一位移再重复写进最终解坐标。 */
    /* exclude measurements of eclipsing satellite (block IIA) */
    if (rtk->opt.posopt[3]) {
        testeclipse(obs, n, nav, rs);
    }
    /* earth tides correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time), rtk->x, opt->tidecorr, &nav->erp,
            opt->odisp[0], dr);
    }
    /* 为 GNSS 的码相行和外部大气约束预留矩阵容量。
     * xp/Pp 是本次试算副本；rtk->x/P 暂保留同一预测起点，待验后接受再提交。
     * 当前 nx 包含 MAXSAT 预留状态，因此内存/耗时不只取决于可见卫星数。 */
    nv = n * rtk->opt.nf * 2 + MAXSAT + 3;
    xp = mat(rtk->nx, 1); Pp = zeros(rtk->nx, rtk->nx);
    v = mat(nv, 1); H = mat(rtk->nx, nv); R = mat(nv, nv);

    /* 第 7 步：测量更新 + 单观测粗差编辑，最多 MAX_ITER 次。
     * mexc 保留本历元累计屏蔽项；x/P 每次却回到同一预测状态重新求解。
     * 这是重新选观测子集，不是把同一组观测连续滤波 16 次。 */
    for (i = 0; i < MAX_ITER; i++) {

        matcpy(xp, rtk->x, rtk->nx, 1);
        matcpy(Pp, rtk->P, rtk->nx, rtk->nx);

        /* prefit residuals */
        if (!(nv = ppp_res(0, obs, n, rs, dts, var, svh, dr, exc, mexc, nav, xp, rtk, v, H, R, azel, NULL, product_use))) {
            trace(2, "%s ppp (%d) no valid obs data\n", str, i + 1);
            break;
        }
        /* filter 在 rtkcmn 中完成 EKF 测量更新：根据 v/H/R 算增益，更新
         * xp/Pp。info!=0 为数值更新失败，不能把该试算状态提交成有效 PPP。 */
        if ((info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {
            trace(2, "%s ppp (%d) filter error info=%d\n", str, i + 1, info);
            break;
        }
        /* 用更新后的 xp 重算残差，而不是复用验前残差判断滤波是否通过。
         * post_ok=0 时仅保留新增 mexc 屏蔽，然后下一轮从原预测状态重算。 */
        {
        int post_ok=ppp_res(i + 1, obs, n, rs, dts, var, svh, dr, exc, mexc, nav, xp, rtk, NULL, NULL, NULL, azel,
            &diag, NULL);
        /* First post-fit residuals include measurements later edited out.
         * Count their evidence once, never once per editing iteration. */
        /* 第一份验后证据只用于持续质量检测，包含后续可能被剔除的观测。
         * 最终输出/CMC 使用验后通过的 diag；两份快照的含义不能混淆。 */
        if(i==0)evidence=diag;
        if (post_ok) {
            /* 仅验后通过才提交测量更新；成功使用的产品先验时间也在此提交。
             * 初值已消耗的先验和编辑失败的试算不会被当作额外独立信息。 */
            matcpy(rtk->x, xp, rtk->nx, 1);
            matcpy(rtk->P, Pp, rtk->nx, rtk->nx);
            stat = SOLQ_PPP;
            /* Trial editing iterations start from the same prior P: they are
               not extra measurements. Commit consumption only on success. */
            if(statistical_model_ppp(opt)&2) {
                for(j=0;j<MAXSAT;j++)if(product_use[3+j])
                    rtk->ssat[j].ppp_ion_prior_time=obs[0].time;
                if(product_use[2])rtk->ppp_zwd_prior_time=obs[0].time;
            }
            break;
        }
        }
    }
    if (i >= MAX_ITER) {
        trace(2, "%s ppp (%d) iteration overflows\n", str, i);
    }
    trace(2, "$PPP_PRODUCT_USE,%s,ionex_constraints=%d,vmf3_zwd_constraints=%d\n",
          str, product_use[0], product_use[1]);

    /* Compact epoch summary for targeted trace diagnosis. */
    {
        int n_sat_exc = 0, n_meas_exc = 0, ii;

        for (ii = 0; ii < n && ii < MAXOBS; ii++) {
            if (exc[ii]) n_sat_exc++;
        }
        for (ii = 0; ii < MAXOBS * 2 * NFREQ; ii++) {
            if (mexc[ii]) n_meas_exc++;
        }
        trace(2, "$PPP_EPOCH,%s,status=%s,iter=%d,sat_exc=%d,meas_exc=%d\n",
            str, stat == SOLQ_PPP ? "OK" : "FAIL", i + 1,
            n_sat_exc, n_meas_exc);
    }

    /* 第 8 步：仅成功载波 PPP 进入可选 AR/状态输出。
     * stat 初始是 SOLQ_SINGLE；未通过滤波时不把预测坐标伪装成 Q6。
     * 外层 rtkpos 对失败/单点回退的处理需同时查看，不能只看本地 stat。 */
    if (stat == SOLQ_PPP) {

        if (opt->modear != ARMODE_OFF &&
            ppp_ar(rtk, obs, n, exc, nav, azel, xp, Pp) &&
            ppp_res(9, obs, n, rs, dts, var, svh, dr, exc, mexc, nav, xp, rtk, NULL, NULL, NULL, azel, NULL, NULL)) {

            matcpy(rtk->xa, xp, rtk->nx, 1);
            matcpy(rtk->Pa, Pp, rtk->nx, rtk->nx);

            for (i = 0; i < 3; i++) std[i] = sqrt(Pp[i + i * rtk->nx]);
            if (norm(std, 3) < MAX_STD_FIX) stat = SOLQ_FIX;
        }
        else {
            rtk->nfix = 0;
        }
        /* 把状态复制为 rtk->sol，按实际通过的载波卫星数决定最终 Q/ns；
         * 隔离只跨历元更新一次，不随前面的试算次数增加。 */
        update_stat(rtk, obs, n, stat);
        update_phase_quarantine_ppp(rtk, obs, n, mexc, &diag);

        /* A successful code-only/mixed update is useful for continuity but
         * is not a carrier PPP solution. Label it degraded (Q=5), never Q=6.
         * This also covers the intentional code+Doppler warm-up. */
        /* 这一 Q5 是成功码更新但载波不足时的真实降级标记，不是为消除
         * 开头 Q5 人为改标签。足够码卫星与足够载波卫星是不同条件。 */
        if (rtk->sol.stat == SOLQ_NONE) {
            int code_sat = 0, phase_sat = 0, oi, fi;
            for (oi = 0; oi < n && oi < MAXOBS; oi++) {
                int has_code = 0, has_phase = 0;
                for (fi = 0; fi < NF(opt) && fi < NFREQ; fi++) {
                    has_code |= diag.obs[oi][fi].code_ok != 0;
                    has_phase |= diag.obs[oi][fi].phase_ok != 0;
                }
                code_sat += has_code;
                phase_sat += has_phase;
            }
            if (code_sat >= MIN_NSAT_SOL) {
                rtk->sol.stat = SOLQ_SINGLE;
                rtk->sol.ns = (uint8_t)code_sat;
                trace(2, "$PPP_DEGRADED,%s,Q=5,code_sat=%d,phase_sat=%d,warmup=%d\n",
                    str, code_sat, phase_sat, rtk->ppp_code_warmup);
            }
        }

        if (stat == SOLQ_FIX && test_hold_amb(rtk)) {
            matcpy(rtk->x, xp, rtk->nx, 1);
            matcpy(rtk->P, Pp, rtk->nx, rtk->nx);
            trace(2, "%s hold ambiguity\n", str);
            rtk->nfix = 0;
        }
    }
    /* 第 9 步：在已成功的滤波历元末更新质量历史，影响后续历元；
     * CMC 基于最终通过的码相对，持续残差检测基于首次验后的 evidence。 */
    if (stat == SOLQ_PPP || stat == SOLQ_FIX) {
        update_cmc_quality_ppp(rtk, obs, n, &diag);
        update_signal_quality_ppp(rtk, obs, n, &evidence);
    }
    if (diag_enabled) {
        ppp_diag_emit(rtk, obs, n, stat == SOLQ_PPP || stat == SOLQ_FIX ? &diag : NULL,
            rtk->sol.stat == SOLQ_PPP || rtk->sol.stat == SOLQ_FIX);
    }
    /* 第 10 步：仅释放本历元临时数组；rtk_t 的 x/P/ssat 和 nav 产品
     * 都由调用层持续保留，整个会话结束才 rtkfree/释放产品。 */
    free(rs); free(dts); free(var); free(azel);
    free(xp); free(Pp); free(v); free(H); free(R);
}
