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
/*==============================================================================
 * 本项目的 PC 程序入口：先选择文件和参数，再把任务交给对应的解算模块。
 *
 * 【先看整体流程】
 *   1. 读取命令行，决定使用文件选择窗口，还是环境变量指定的自动运行方式。
 *   2. 观测输入二选一：RINEX OBS 或 Android 原始 TXT，不能同时作为输入。
 *   3. 检查文件路径，整理产品路径，生成一份统一的 PPP 处理参数。
 *   4. 按输入类型分流：
 *      RINEX -> postpos() -> RINEX 观测读取 -> 按历元解算 -> POS/stat/trace。
 *      TXT   -> pc_txt_replay / pc_live_view -> gnss_replay -> gnss_adapter
 *            -> ppp-safe 信号选择 -> 持续同一个 rtk_t -> rtkpos()。
 *   5. 检查结果、关闭日志，结束程序。
 *
 * 【这个文件负责什么，不负责什么】
 *   负责启动、输入选择、参数传递、模块调用和完成/错误提示。
 *   不在这里实现卡尔曼滤波、周跳检测、定权、IONEX/VMF3观测方程。
 *   上述模型主要在 ppp.c / rtkpos.c 等核心文件中；调整窗口不会重置模型。
 *   TXT 的解析、时间转换、相位/多普勒转换在 replay/adapter 中实现。
 *
 * 【相关文件索引】
 *   smartphone_ppp_config.c : 共用的默认处理参数、配置检查和配置签名日志。
 *   pc_file_picker.c        : 中文文件选择窗口，记住上次选择，检查文件路径。
 *   pc_txt_replay.c         : 把选定 TXT/产品/参数交给既有回放器的桥接层。
 *   pc_live_view.c          : 后台回放 + 只读坐标窗口，可安全停止回放。
 *   gnss_replay.c           : 按原始 GNSS 时间回放，加载产品，逐历元调用 rtkpos。
 *   postpos.c               : RINEX 后处理任务及内部历元处理流程。
 *
 * 【运行方式】
 *   正常双击：打开窗口，选择输入类型、观测、产品及结果保存位置。
 *   自动测试：--no-gui，并用 RTK_PPP_* 环境变量指定文件/可选参数。
 *   TXT 模拟实时按原始历元时间 1×回放，不等于 PC 已连接手机在线采集。
 *   本次中文注释仅用于说明现有行为，不改默认参数、不改解算逻辑。
 *==============================================================================*/
#include <stdarg.h>
#include "./rtklib.h"                  /* RTKLIB 公共类型、常量和解算接口。 */
#include "smartphone_ppp_config.h"      /* PC/replay/JNI 共享默认配置接口。 */
#include "pc_file_picker.h"             /* 文件列表、输入类型和中文文件窗口。 */
#include "pc_txt_replay.h"              /* TXT 无窗口回放入口。 */
#include "pc_live_view.h"               /* TXT 坐标窗口及后台回放入口。 */

#define PROGNAME    "rnx2rtkp"          /* 保留原程序名宏，本入口未据此生成帮助标题。 */
#define MAXFILE     16                  /* infile数组容量；RINEX实际传入5个文件。 */

/* Optional overrides keep controlled A/B runs from overwriting the normal
 * observation, solution and trace files. Empty variables retain defaults. */
/* 从一个环境变量安全地复制字符串到已有缓冲区。
 * dst      : 接收路径/参数的缓冲区；原值就是未覆盖时的默认值。
 * capacity : 缓冲区总容量，包含最后的 '\0'。
 * name     : 环境变量名，例如 RTK_PPP_OUTPUT。
 * 返回1表示可继续，0表示变量过长，应停止而不是截断路径后继续解算。
 * 注意：变量不存在或为空时“不改dst”；它不是清空字段的指令。
 * 此处只读取进程环境，不会写入系统环境或修改源码中的默认值。 */
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

/* 帮助只描述当前入口接受的选项。
 * 原版rnx2rtkp的所有历史选项并未在本入口完整保留，不能仅依据上面的
 * 历史记录就认为-k/-h等选项可用。路径环境变量主要用于可复现的自动测试。 */
/* help text -----------------------------------------------------------------*/
static const char* help[] = {
"RTKLIB PC PPP",
"  No arguments: open the file-selection window (Windows).",
"  --select-files: explicitly open the window, even with environment overrides.",
"  --no-gui: unattended processing with existing defaults/environment variables.",
"  --check-files: validate paths only; do not run PPP.",
"  --source rinex|txt-realtime|txt-fast: unattended observation source.",
"  RTK_PPP_TXT selects Android raw TXT; never set it together with RTK_PPP_OBS.",
"  TXT uses ppp-safe, ADR uncertainty limit 1.0m, no RINEX OBS dependency.",
"  --max-epochs N: TXT diagnostic limit only (normally replay the whole file).",
"  --live-view: display TXT epoch coordinates during an unattended-path run.",
"  --help: display this help.",
"Products can be overridden for unattended runs using RTK_PPP_OBS,",
"RTK_PPP_NAV, RTK_PPP_SP3, RTK_PPP_CLK, RTK_PPP_BIA, RTK_PPP_IONEX,",
"RTK_PPP_VMF_FIRST, RTK_PPP_VMF_SECOND, RTK_PPP_OROGRAPHY, RTK_PPP_ATX.",
"RTK_PPP_OUTPUT and RTK_PPP_TRACE select diagnostic output files.",
"Existing RTK_PPP_NF/SYSTEMS/GLO/EXSATS/OPTS options are unchanged.",
"An explicit RTK_PPP_OBS retains legacy unattended behavior.",
};
/* show message --------------------------------------------------------------*/
/* 后处理需要的进度回调：format/...类似printf，由内部模块传入。
 * 信息写到标准错误输出；'\r'回到当前行开头，适合刷新进度，不是换行。
 * 返回0表示本回调不请求中断。TXT窗口的停止走独立的回放取消接口，
 * 并不是通过此showmsg回调实现。 */
extern int showmsg(const char* format, ...)
{
    va_list arg;
    va_start(arg, format); vfprintf(stderr, format, arg); va_end(arg);
    fprintf(stderr, "\r");
    return 0;
}
/* 保留后处理要求的两个回调符号。目前是空实现：
 * settspan可通知界面总起止时间，settime可通知当前处理时间。
 * 本项目的TXT坐标刷新使用独立结果快照，不在这里更新。
 * settime只是RTKLIB的进度回调，不会修改Windows系统时间。 */
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}

/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i = 0; i < (int)(sizeof(help) / sizeof(*help)); i++) fprintf(stderr, "%s\n", help[i]);
    exit(0); /* 请求帮助不是一次解算，打印完立即结束，不读取产品。 */
}
/* 统计RINEX分支生成的POS中“非空且非头部”的数据行。
 * 跳过行首空白、%/#说明和空行，用于发现没有新结果的情况。
 * 不是精度/收敛验证，也没有检查每一行是否Q6；有行不代表厘米级定位。
 * TXT有效历元由replay内部统计，不用本函数统计CSV（避免把CSV表头当结果）。 */
static long solution_rows(const char *path)
{
    FILE *fp=fopen(path,"r");
    char line[4096],*s;
    long rows=0;
    if(!fp) return 0;
    while(fgets(line,sizeof(line),fp)) {
        for(s=line;*s==' '||*s=='\t';s++); /* 空循环，仅移动s到第一个非空白字符。 */
        if(*s && *s!='%' && *s!='#' && *s!='\r' && *s!='\n') rows++;
    }
    fclose(fp); return rows;
}
/*============================ 主程序入口 ===================================*/
int main(int argc, char **argv) {
    /* 第1部分：启动方式与文件容器。
     * argc/argv来自命令行，argv[0]是程序名，后面才是用户选项。
     * getenv返回进程环境中的只读字符串，不需要也不能free。
     * 观测环境变量非空时默认自动运行；否则默认打开文件窗口。
     * 后面的启动方式选项还可覆盖这个gui标志。 */
    const char *obs_env=getenv("RTK_PPP_OBS"), *txt_env=getenv("RTK_PPP_TXT");
    int i, n=5, ret, gui=!((obs_env && *obs_env)||(txt_env && *txt_env));
    int check_only=0, requested_source=-1, realtime=1, max_epochs=0, live_view=0;
    /* gui：是否显示文件选择窗口，不等于是否显示坐标窗口。
     * check_only：只检查路径，不加载产品、不初始化滤波器。
     * requested_source：-1为未指定；0/1由PC_SOURCE_RINEX/TXT定义。
     * realtime：TXT默认1×回放；txt-fast置0，只跳过等待，不改PPP参数。
     * max_epochs：0不人为截断TXT；正数用于短时回放测试。
     * live_view：不用文件选择窗口时，也可显式要求显示TXT坐标窗口。 */
    double tint=0.0;
    gtime_t ts={0},te={0};
    /* ts/te为零表示本入口不指定后处理起止时间，使用输入可用时段。
     * tint=0表示不指定后处理抽样间隔，不是把观测时间改成0。
     * 这三个变量只传给RINEX postpos，不控制TXT的单调时钟调度。 */
    const char *infile[MAXFILE];
    /* 下列路径仅为历史“无窗口自动运行”的兼容默认值。
     * 窗口启动时会将files清零，不用旧日期文件自动填入表单。
     * 换日期不必改这里：窗口选文件或用环境变量覆盖即可。
     * path顺序必须与pc_file_picker.h的PC_*枚举对应：观测、NAV、SP3、
     * CLK、BIA、IONEX、VMF前/后、格网高程、ATX、输出。
     * 未显式初始化的source为0，即RINEX；稍后按选定输入重新确定。 */
    pc_ppp_files_t files={{
        "E:\\RTKLIB_Data\\OBS\\K80-GNSS00GEO_R_20260760533_20M_01S_MO.rnx",
        "E:\\RTKLIB_Data\\NAV\\BRDC00IGS_R_20260760000_01D_MN.rnx",
        "E:\\RTKLIB_Data\\SP3\\WUM0MGXFIN_20260760000_01D_05M_ORB.SP3",
        "E:\\RTKLIB_Data\\CLK\\WUM0MGXFIN_20260760000_01D_30S_CLK.CLK",
        "E:\\RTKLIB_Data\\BIA\\WUM0MGXFIN_20260760000_01D_01D_OSB.BIA",
        "E:\\RTKLIB_Data\\IONEX\\COD0OPSFIN_20260760000_01D_01H_GIM.INX",
        "E:\\RTKLIB_Data\\tro\\VMF3_20260317.H00",
        "E:\\RTKLIB_Data\\tro\\VMF3_20260317.H06",
        "E:\\RTKLIB_Data\\tro\\orography_ell_5x5",
        "E:\\RTKLIB_Data\\Tables\\igs20.atx",
        "E:\\RTKLIB_Data\\PPP_Result\\K80-GNSS00GEO_R_20260760533_20M_01S_MO.pos"
    }};
    /* 环境变量名与path字段一一对应。
     * PC_OBS其实是“唯一观测入口槽”，也可放TXT路径；不是同时存两种观测。
     * TXT单独读取RTK_PPP_TXT，所以后面产品覆盖从索引1开始，
     * 不会把RTK_PPP_OBS的内容误复制给TXT入口。 */
    const char *file_env[PC_FILE_COUNT]={
        "RTK_PPP_OBS","RTK_PPP_NAV","RTK_PPP_SP3","RTK_PPP_CLK","RTK_PPP_BIA",
        "RTK_PPP_IONEX","RTK_PPP_VMF_FIRST","RTK_PPP_VMF_SECOND",
        "RTK_PPP_OROGRAPHY","RTK_PPP_ATX","RTK_PPP_OUTPUT"
    };
    char outfile[MAXSTRPATH],trace_file[MAXSTRPATH],file_error[2048];
    filopt_t filopt={0}; /* RINEX辅助文件选项，区别于PPP处理参数prcopt。 */
    /* before/after比较文件修改时间，防止把旧POS当成本次新结果。 */
    WIN32_FILE_ATTRIBUTE_DATA before={0},after={0};
    int output_existed;
    long rows;
    /* 第2部分：解析当前入口支持的命令行选项。
     * 从左到右解析，后面的启动方式选项可覆盖前面的gui标志。
     * --source会关闭文件窗口；其后加--select-files可重新打开窗口。
     * --live-view只要求坐标窗口，不自行选文件或改变输入类型。
     * 未知参数直接报错，不静默忽略，避免误以为某设置已经生效。 */
    for(i=1;i<argc;i++) {
        if(!strcmp(argv[i],"--help") || !strcmp(argv[i],"-?")) printhelp();
        else if(!strcmp(argv[i],"--select-files")) gui=1;
        else if(!strcmp(argv[i],"--no-gui")) gui=0;
        else if(!strcmp(argv[i],"--check-files")) {gui=0;check_only=1;}
        else if(!strcmp(argv[i],"--live-view")) live_view=1;
        else if(!strcmp(argv[i],"--source") && i+1<argc) {
            /* ++i消耗紧跟选项的值，避免下一轮把这个值当成新选项。
             * txt-realtime/fast共用同一个回放入口，区别仅为是否按时间等待。 */
            const char *source=argv[++i];gui=0;
            if(!strcmp(source,"rinex")) requested_source=PC_SOURCE_RINEX;
            else if(!strcmp(source,"txt-realtime")) requested_source=PC_SOURCE_TXT,realtime=1;
            else if(!strcmp(source,"txt-fast")) requested_source=PC_SOURCE_TXT,realtime=0;
            else {fprintf(stderr,"Unknown source: %s\n",source);return 2;}
        }
        else if(!strcmp(argv[i],"--max-epochs") && i+1<argc) {
            /* strtol检查整串是否有效正整数；end非空表示有非法尾随字符。 */
            char *end;long value=strtol(argv[++i],&end,10);
            if(!argv[i][0] || *end || value<=0 || value>2147483647L) {
                fprintf(stderr,"--max-epochs must be a positive integer\n");return 2;
            }
            max_epochs=(int)value;
        }
        else {fprintf(stderr,"Unknown argument: %s (use --help)\n",argv[i]);return 2;}
    }
    /* 第3部分：从窗口或环境变量取得唯一的观测输入与共用产品。 */
    if(gui) {
        /* First GUI run starts blank; only saved choices populate the form.
           Historical hardcoded products are for legacy automated runs only. */
        memset(&files,0,sizeof(files));
        /* pc_file_picker只选文件，不运行PPP。
         * 返回1：点击开始且检查通过；0：用户取消；负数：窗口等错误。
         * 取消不会运行历史默认观测，也不会覆盖结果。 */
        ret=pc_file_picker(&files,file_error,sizeof(file_error));
        if(ret<=0) {
            if(ret<0) fprintf(stderr,"%s\n",file_error);
            return ret<0 ? 2 : 0;
        }
    }
    else {
        /* 自动模式RINEX/TXT只能二选一，两种环境输入同时非空就拒绝。
         * 窗口模式不执行此else，因此--select-files可摆脱VS中旧的路径环境。
         * 但下面的PPP参数环境覆盖仍会执行，不是清除所有环境参数。 */
        if(obs_env && *obs_env && txt_env && *txt_env) {
            fprintf(stderr,"INPUT_CONFLICT: choose RTK_PPP_OBS or RTK_PPP_TXT, never both\n");return 2;
        }
        files.source=requested_source>=0 ? requested_source :
                     (txt_env && *txt_env ? PC_SOURCE_TXT : PC_SOURCE_RINEX);
        /* 命令行显式source优先；未指定时非空TXT选TXT，其余选RINEX。
         * 命令行与观测环境输入类型不一致时停止，不擅自替用户选一种。 */
        if((files.source==PC_SOURCE_TXT && obs_env && *obs_env) ||
           (files.source==PC_SOURCE_RINEX && txt_env && *txt_env)) {
            fprintf(stderr,"INPUT_CONFLICT: source disagrees with observation environment\n");return 2;
        }
        for(i=1;i<PC_FILE_COUNT;i++)
            if(!env_copy(files.path[i],sizeof(files.path[i]),file_env[i])) return 2;
        if(files.source==PC_SOURCE_TXT) {
            /* env_copy遇到空变量会保留旧值，因此先清空唯一观测槽。
             * TXT未指定时后面检查失败，绝不回退到旧RINEX。
             * 自动TXT还必须指定输出路径，防止覆盖历史POS。 */
            files.path[PC_OBS][0]=0; /* Never fall back to the historical RINEX. */
            if(!env_copy(files.path[PC_OBS],sizeof(files.path[PC_OBS]),"RTK_PPP_TXT")) return 2;
            if(!(getenv("RTK_PPP_OUTPUT") && *getenv("RTK_PPP_OUTPUT"))) {
                /* No legacy result can be overwritten by a new TXT run. */
                fprintf(stderr,"TXT unattended mode requires RTK_PPP_OUTPUT (CSV)\n");return 2;
            }
        }
        else if(!env_copy(files.path[PC_OBS],sizeof(files.path[PC_OBS]),"RTK_PPP_OBS")) return 2;
    }
    /* 以下两项仅适用于TXT，不假装给RINEX后处理也设置了这些功能。 */
    if(max_epochs && files.source!=PC_SOURCE_TXT) {
        fprintf(stderr,"--max-epochs is only supported for TXT replay\n");return 2;
    }
    if(live_view && files.source!=PC_SOURCE_TXT) {
        fprintf(stderr,"--live-view is only supported for TXT replay\n");return 2;
    }
    /* 只检查文件存在/可读、输出目录和不覆盖所选输入等路径条件。
     * 不证明产品时间覆盖、内容正确、OSB支持或最终精度；
     * 同名文件也可能属于旧日期，格式/时间检查仍由产品读取器完成。 */
    if(!pc_files_validate(&files,file_error,sizeof(file_error))) {
        fprintf(stderr,"FILE_SELECTION_ERROR: %s\n",file_error);return 2;
    }
    /* 打印实际选定路径，便于发现误用旧日期产品，不是打印默认路径列表。 */
    printf("PPP_INPUT,SOURCE=%s\n",files.source==PC_SOURCE_TXT ? "ANDROID_TXT" : "RINEX_OBS");
    for(i=0;i<PC_FILE_COUNT;i++)
        printf("PPP_INPUT,%s=%s\n",i==PC_OBS && files.source==PC_SOURCE_TXT ? "TXT" : pc_file_keys[i],files.path[i]);
    if(check_only) {puts("FILE_CHECK_OK (existence/readability only; not product coverage)");return 0;}
    /* 第4部分：整理输出与RINEX辅助路径。
     * 这里只复制字符串，不读取观测或产品。
     * TXT也执行这些公共赋值，但在分流处提前return，不调用postpos。 */
    strcpy(outfile,files.path[PC_OUTPUT]);
    /* Preserve the original input ordering and bias loading route. */
    infile[0]=files.path[PC_OBS]; infile[1]=files.path[PC_NAV];
    infile[2]=files.path[PC_SP3]; infile[3]=files.path[PC_CLK];
    infile[4]=files.path[PC_IONEX];
    /* infile是postpos输入列表；BIA/ATX等辅助路径由filopt传入。
     * filopt.dcb沿用历史字段名，不代表只能读旧式DCB；本项目偏差读取器
     * 支持已接入的BIA/OSB。这里不直接解析产品内容。
     * VMF3和格网高程走后面的pppvmf3load，不计入这5个infile条目。 */
    strcpy(filopt.satantp,files.path[PC_ATX]);
    strcpy(filopt.iono,files.path[PC_IONEX]);
    strcpy(filopt.dcb,files.path[PC_BIA]);
    /* 窗口模式/TXT默认将trace放在结果旁，如result.csv.trace。
     * 纯自动RINEX保留历史ppp_debug.trace；无文件窗口时可用RTK_PPP_TRACE
     * 覆盖，即使有live_view坐标窗口也可使用这一覆盖。 */
    if(gui || files.source==PC_SOURCE_TXT) {
        snprintf(trace_file,sizeof(trace_file),"%s.trace",outfile);
        if(!gui && !env_copy(trace_file,sizeof(trace_file),"RTK_PPP_TRACE")) return 2;
    }
    else {
        strcpy(trace_file,"E:\\RTKLIB_Data\\ppp_debug.trace");
        if(!env_copy(trace_file,sizeof(trace_file),"RTK_PPP_TRACE")) return 2;
    }
    /* 留16字节给.stat/.trace/.obs.csv等伴随文件后缀。 */
    if(strlen(outfile)+16>=sizeof(outfile)) {fprintf(stderr,"Output path too long\n");return 2;}
    prcopt_t prcopt = prcopt_default;
    solopt_t solopt = solopt_default;
    /* 第5部分：生成一份处理参数，再交给选中的解算链。
     * prcopt_t：系统/频槽/模型/随机模型/PPP扩展选项等处理参数。
     * solopt_t：输出时间、坐标格式等，不是滤波状态。
     * filopt_t：辅助路径，不是已加载产品nav_t。
     * 先取RTKLIB默认值，再调用共享配置统一为当前手机PPP基线。
     * 目前基线为静态PPP、GPS+Galileo+BDS、nf=4、精密星历，估计电离层
     * 和对流层状态，不固定模糊度；具体数值以共享配置为准。
     * “实时节奏回放”不会自动切换成动态定位模型，这里未另改mode。 */
    /* All entry points share the current PC baseline. */
    smartphone_ppp_configure(&prcopt, &solopt);
    {
        /* Reproducible B1C ablation without changing any other PPP option. */
        const char *nf_env = getenv("RTK_PPP_NF");
        /* nf是允许使用的RTKLIB频率槽数，不代表每颗卫星都有4个频率。
         * 这里只接受3/4，用于已有B1C第四槽的对照；3覆盖为3，4保留基线。
         * 具体有何信号仍取决于观测读取器/adapter/信号策略。 */
        if (nf_env && *nf_env) {
            if (!strcmp(nf_env, "3")) prcopt.nf = 3;
            else if (strcmp(nf_env, "4")) {
                fprintf(stderr, "RTK_PPP_NF must be 3 or 4\n");
                return -1;
            }
        }
    }
    /* 非空RTK_PPP_OPTS会整体替换pppopt，而不是追加一个选项。
     * pppopt包含本项目IONEX、VMF3、预处理等自定义参数，核心解释其含义，
     * 这里仅传递字符串。只填一个选项可能丢失其余基线选项，
     * 做对照时应保留需要的整套选项并核对PPP_CONFIG日志。 */
    if (!env_copy(prcopt.pppopt, sizeof(prcopt.pppopt), "RTK_PPP_OPTS")) return -1;
    {
        const char *systems = getenv("RTK_PPP_SYSTEMS");
        /* 系统按位组合，不是卫星颗数：G=GPS，R=GLONASS，E=Galileo，
         * C=BeiDou，J=QZSS；例如GEC表示GPS+Galileo+BDS。
         * 非空字符串整体替换navsys，未知字符直接报错。
         * 选择某系统不保证adapter支持它或精密产品覆盖其全部观测。 */
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
        /* 在SYSTEMS之后执行，所以GLO=1/0优先于SYSTEMS中的R。
         * |=打开GLONASS位，&=~关闭该位，不改其他星座。
         * 此开关仅用于受控实验，不解决GLONASS偏差/产品支持问题。 */
        if (glo && *glo) {
            if (!strcmp(glo, "1")) prcopt.navsys |= SYS_GLO;
            else if (!strcmp(glo, "0")) prcopt.navsys &= ~SYS_GLO;
            else {
                fprintf(stderr, "RTK_PPP_GLO must be 0 or 1\n");
                return -1;
            }
        }
    }
    /* 运行时选GLONASS还需编译时ENAGLO支持。NSATGLO为0时没有该星座
     * 的卫星容量，不能靠运行配置字符串补出编译时未启用的能力。 */
    if ((prcopt.navsys & SYS_GLO) && NSATGLO == 0) {
        fprintf(stderr, "GLONASS requested but ENAGLO is not enabled in this build\n");
        return -1;
    }
    {
        /* Controlled satellite leave-one-out tests; no default exclusion. */
        const char *excluded = getenv("RTK_PPP_EXSATS");
        /* 可选手动排星诊断，例如G12,C19，逗号/空格均可分隔。
         * satid2no将卫星ID转内部编号；sat从1开始，数组从0开始，故用sat-1。
         * exsats置1表示排除，不是强制使用；默认不在此入口额外排星。
         * 这是人为实验设置，不等同于核心残差粗差剔除或持续异常检测。 */
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
        /* 对所有覆盖后的最终参数做结构/取值检查，不只检查默认配置。
         * 检查通过不代表观测质量/精度合格，仍需查看真实结果及trace。 */
        if (!smartphone_ppp_validate(&prcopt, config_error, sizeof(config_error))) {
            fprintf(stderr, "PPP configuration error: %s\n", config_error);
            return -1;
        }
    }

    /* 第6部分：TXT分流，执行后return，不继续运行RINEX后处理。
     * 同一份prcopt直接传入回放器，不另设一套“结果更好看”的参数。
     * gui/live_view任一为真：坐标窗口的后台线程启动同一回放器，
     * 主线程响应窗口；界面只读快照，不能直接修改滤波器状态。
     * 两者均假：无窗口回放，适合自动验证。
     * 产品加载、一次rtkinit、逐历元rtkpos、结束rtkfree在回放器内部，
     * 不是此main每秒重新加载/初始化。
     * realtime仅控制时间等待；max_epochs仅控制实验截断，不进入观测方程。 */
    if(files.source==PC_SOURCE_TXT) {
        /* Exactly one runner; postpos() and RINEX OBS loading are not invoked. */
        ret=(gui || live_view) ? pc_live_view_run(&files,&prcopt,trace_file,realtime,max_epochs) :
                  pc_txt_replay(&files,&prcopt,trace_file,realtime,max_epochs);
        /* TXT返回0正常完成，3用户安全停止，其余为错误。
         * CSV逐历元保存位置/Q/ns，stat保存状态，obs.csv保存观测转换字段。
         * 完成不表示所有历元均Q6，也不表示达到某个精度阈值。 */
        if(!ret) printf("\nTXT replay completed. Result: %s\nState: %s.stat\nObservations: %s.obs.csv\nTrace: %s\n",
            outfile,outfile,outfile,trace_file);
        else if(ret==3) printf("TXT replay stopped by user. Completed epochs saved to: %s\n",outfile);
        else fprintf(stderr,"TXT replay failed (status=%d).\n",ret);
        return ret;
    }
    /* 第7部分：只有RINEX输入走到这里。trace是诊断日志，不是POS/stat。
     * 默认等级3，等级越高通常日志越详细、文件越大；不是定位质量等级，
     * 不会通过提高tracelevel把单点Q5变成PPP Q6。 */
    printf("Starting Smartphone Un-combined PPP Processing...\n");
    {
        const char *level=getenv("RTK_PPP_TRACE_LEVEL");
        traceopen(trace_file);
        tracelevel(level ? atoi(level) : 3);
    }
    /* 记录最终配置和signature，核对不同入口/实验是否真的用了相同参数。 */
    smartphone_ppp_log("pc-postprocess",&prcopt);
    /* 启动后处理前加载VMF3前后两期及独立格网高程。
     * 文件加载失败明确停止，不在本入口静默回退其他对流层模型。
     * 是否插值/约束ZWD/具体权重由既有核心模型和pppopt决定，不在这里设置。
     * 两期应覆盖观测；能读入文件不代表每个历元都在产品有效时段内。 */
    if(!pppvmf3load(files.path[PC_VMF_FIRST],files.path[PC_VMF_SECOND],
                    files.path[PC_OROGRAPHY])) {
        fprintf(stderr,"ERROR: selected VMF3 files could not be loaded; no PPP run started.\n");
        traceclose();return 2;
    }
    /* 第8部分：记录原输出文件状态，再执行RINEX后处理。
     * 同名文件已有时，记录其修改时间以供后面检查结果新鲜度。
     * clock计时包含postpos读取/解算，不是逐历元processing_ms，也不是
     * TXT调度用的单调时钟，不能拿此总耗时验证1×回放是否发生时间漂移。 */
    output_existed=GetFileAttributesExA(outfile,GetFileExInfoStandard,&before);
    long t1=clock();
    /* postpos各参数：
     * ts/te：起止时间；tint：抽样间隔；0.0：本入口不指定分段时长。
     * prcopt：处理参数；solopt：输出格式；filopt：偏差/天线等辅助文件路径。
     * infile/n：输入列表/实际数量；outfile：结果位置。
     * 最后两空字符串：不在此指定流动站/基准站ID筛选列表。
     * 真正的观测读取、产品读取及按历元滤波由postpos内部完成。
     * 此main不会逐颗卫星读取后单独调用一次rtkpos。 */
    ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,"","");
    long t2=clock();
    /* 第9部分：防止“内部读入失败但返回0”被误报成功。
     * 部分历史后处理失败分支也返回0，所以还要求至少一行结果、文件存在，
     * 已有结果文件的修改时间发生变化。不满足则ret改为2，表示没有可确认
     * 的新结果。这是入口级检查，不是完整POS解析或精度验证。
     * 同时间戳更新也可能保守判作未刷新，应结合文件内容和trace进一步查验。 */
    rows=solution_rows(outfile);
    if(!ret && (!rows || !GetFileAttributesExA(outfile,GetFileExInfoStandard,&after) ||
        (output_existed && !CompareFileTime(&before.ftLastWriteTime,&after.ftLastWriteTime)))) ret=2;
    if(!ret) {
        printf("\nProcessing completed: %ld solution rows.\n",rows);
        printf("* Result saved to: %s\n* Trace saved to: %s\n",outfile,trace_file);
    }
    else fprintf(stderr,"Processing failed or produced no fresh solution rows (status=%d).\n",ret);
    printf("* The total processing time: %6.3f seconds\n",(double)(t2-t1)/CLOCKS_PER_SEC);
    /* 第10部分：关闭RINEX trace并退出。
     * 仅“窗口选文件的RINEX”等待回车，让用户看清控制台提示。
     * TXT上面已经返回，坐标窗口负责其关闭，不执行此处getchar。
     * 常见返回值：0正常完成/取消文件选择，2输入或结果检查错误，
     * -1本入口部分配置错误，3仅表示TXT用户安全停止。 */
    traceclose();
    if(gui) {printf("Press 'Enter' key to exit...\n");getchar();}
    return ret;
}
