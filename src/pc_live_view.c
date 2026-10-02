/* Read-only PPP result display. No navigation products/filter state in the UI. */
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <math.h>
#include <stdio.h>
#include <wchar.h>
#include "pc_live_view.h"

#define WM_EPOCH (WM_APP+10)
#define WM_FINISHED (WM_APP+11)
#define STOP_ID 400
enum { VALUE_TIME, VALUE_Q, VALUE_NS, VALUE_LAT, VALUE_LON, VALUE_HEIGHT,
       VALUE_MS, VALUE_COUNT };
typedef struct {
    const pc_ppp_files_t *files;
    const prcopt_t *options;
    const char *trace_path;
    int realtime, max_epochs;
    HWND window, value[VALUE_COUNT], banner, history, button;
    HFONT font, large_font;
    HANDLE worker;
    CRITICAL_SECTION lock;
    gnss_replay_epoch_t latest;
    int epochs, drawn_epochs, rc, close_requested;
    volatile LONG cancelled, finished, update_pending;
    double scale;
} live_view_t;

static int cancelled(void *context)
{
    live_view_t *p=(live_view_t *)context;
    return InterlockedCompareExchange(&p->cancelled,0,0)!=0;
}
static void on_epoch(void *context,const gnss_replay_epoch_t *epoch)
{
    live_view_t *p=(live_view_t *)context;
    EnterCriticalSection(&p->lock);
    p->latest=*epoch;p->epochs++;
    LeaveCriticalSection(&p->lock);
    /* Coalesce fast replay updates; every epoch remains in the output CSV. */
    if(!InterlockedExchange(&p->update_pending,1) &&
       !PostMessageW(p->window,WM_EPOCH,0,0))
        InterlockedExchange(&p->update_pending,0);
}
static DWORD WINAPI replay_worker(void *context)
{
    live_view_t *p=(live_view_t *)context;
    gnss_replay_observer_t observer={on_epoch,cancelled,p};
    p->rc=pc_txt_replay_observed(p->files,p->options,p->trace_path,
                                p->realtime,p->max_epochs,&observer);
    InterlockedExchange(&p->finished,1);
    PostMessageW(p->window,WM_FINISHED,0,0);
    return 0;
}
static HWND control(live_view_t *p,const wchar_t *type,const wchar_t *text,
                    DWORD style,int x,int y,int w,int h,int id)
{
    HWND child=CreateWindowExW(!wcscmp(type,L"EDIT")?WS_EX_CLIENTEDGE:0,
        type,text,WS_CHILD|WS_VISIBLE|style,(int)(x*p->scale),(int)(y*p->scale),
        (int)(w*p->scale),(int)(h*p->scale),p->window,(HMENU)(INT_PTR)id,
        GetModuleHandleW(NULL),NULL);
    SendMessageW(child,WM_SETFONT,(WPARAM)p->font,TRUE);
    return child;
}
static void render(live_view_t *p)
{
    gnss_replay_epoch_t r;
    double pos[3],tow;
    int week,count,q,valid;
    char time_text[64];
    wchar_t text[512],clock_text[64];
    EnterCriticalSection(&p->lock);r=p->latest;count=p->epochs;LeaveCriticalSection(&p->lock);
    if(!count || count==p->drawn_epochs) return;
    p->drawn_epochs=count;q=r.solution.stat;
    valid=q!=SOLQ_NONE && isfinite(r.solution.rr[0]) &&
          isfinite(r.solution.rr[1]) && isfinite(r.solution.rr[2]) &&
          norm(r.solution.rr,3)>1.0;
    tow=time2gpst(r.time,&week);time2str(r.time,time_text,3);
    MultiByteToWideChar(CP_ACP,0,time_text,-1,clock_text,64);
    swprintf_s(text,512,L"%ls  |  week %d / TOW %.3f",clock_text,week,tow);
    SetWindowTextW(p->value[VALUE_TIME],text);
    swprintf_s(text,512,L"Q=%d  %ls",q,q==SOLQ_PPP ? L"PPP 浮点解" :
        q==SOLQ_SINGLE ? L"单点解（尚未进入 PPP）" : q==SOLQ_NONE ? L"无有效解" : L"其他解算状态");
    SetWindowTextW(p->value[VALUE_Q],text);
    swprintf_s(text,512,L"解算使用 %d  /  输入 %d",r.solution.ns,r.nobs);
    SetWindowTextW(p->value[VALUE_NS],text);
    if(valid) {
        ecef2pos(r.solution.rr,pos);
        swprintf_s(text,512,L"%.9f°",pos[0]*R2D);SetWindowTextW(p->value[VALUE_LAT],text);
        swprintf_s(text,512,L"%.9f°",pos[1]*R2D);SetWindowTextW(p->value[VALUE_LON],text);
        swprintf_s(text,512,L"%.4f m（椭球高）",pos[2]);SetWindowTextW(p->value[VALUE_HEIGHT],text);
    }
    else {
        SetWindowTextW(p->value[VALUE_LAT],L"暂无有效坐标");
        SetWindowTextW(p->value[VALUE_LON],L"暂无有效坐标");
        SetWindowTextW(p->value[VALUE_HEIGHT],L"暂无有效坐标");
    }
    swprintf_s(text,512,L"%.3f ms  |  已处理 %d 历元",r.processing_ms,count);
    SetWindowTextW(p->value[VALUE_MS],text);
    if(!InterlockedCompareExchange(&p->finished,0,0) && !cancelled(p))
        SetWindowTextW(p->banner,p->realtime ? L"运行中：按原始时间 1×回放；产品已加载，滤波状态持续保留。" : L"运行中：快速诊断回放；界面显示最新历元。");
    if(GetWindowTextLengthW(p->history)>24000) SetWindowTextW(p->history,L"");
    if(valid) swprintf_s(text,512,L"%ls  Q=%d ns=%d  %.9f  %.9f  %.4f m  %.3f ms\r\n",
        clock_text,q,r.solution.ns,pos[0]*R2D,pos[1]*R2D,pos[2],r.processing_ms);
    else swprintf_s(text,512,L"%ls  Q=%d ns=%d  无有效坐标  %.3f ms\r\n",clock_text,q,r.solution.ns,r.processing_ms);
    SendMessageW(p->history,EM_SETSEL,(WPARAM)-1,(LPARAM)-1);
    SendMessageW(p->history,EM_REPLACESEL,FALSE,(LPARAM)text);
    InvalidateRect(p->value[VALUE_Q],NULL,TRUE);
}
static void stop(live_view_t *p)
{
    InterlockedExchange(&p->cancelled,1);
    SetWindowTextW(p->banner,L"正在停止，等待当前加载/历元结束并关闭结果文件，请稍候…");
    EnableWindow(p->button,FALSE);
}
static LRESULT CALLBACK window_proc(HWND window,UINT msg,WPARAM wp,LPARAM lp)
{
    live_view_t *p=(live_view_t *)GetWindowLongPtrW(window,GWLP_USERDATA);
    if(msg==WM_NCCREATE) {
        p=(live_view_t *)((CREATESTRUCTW *)lp)->lpCreateParams;p->window=window;
        SetWindowLongPtrW(window,GWLP_USERDATA,(LONG_PTR)p);
    }
    if(msg==WM_CREATE && p) {
        const wchar_t *labels[]={L"GNSS 时间（GPST）",L"解算状态",L"卫星数",L"纬度 Latitude",L"经度 Longitude",L"高程 Height",L"处理耗时"};
        wchar_t output[PC_PPP_PATH_CAP];int i;
        p->font=CreateFontW((int)(-16*p->scale),0,0,0,FW_NORMAL,0,0,0,DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Microsoft YaHei UI");
        p->large_font=CreateFontW((int)(-21*p->scale),0,0,0,FW_MEDIUM,0,0,0,DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Microsoft YaHei UI");
        p->banner=control(p,L"STATIC",L"正在加载本地产品，尚未开始回放…",0,20,15,785,36,0);
        for(i=0;i<VALUE_COUNT;i++) {
            control(p,L"STATIC",labels[i],0,20,65+i*40,175,32,0);
            p->value[i]=control(p,L"STATIC",L"—",SS_LEFT,200,62+i*40,600,34,500+i);
            if(i>=VALUE_LAT && i<=VALUE_HEIGHT) SendMessageW(p->value[i],WM_SETFONT,(WPARAM)p->large_font,TRUE);
        }
        control(p,L"STATIC",L"Q=5 为单点解，Q=6 为 PPP 浮点解；Q值不代表已达到厘米级。\n这是 TXT 模拟实时输入，不是从手机在线接收数据。",0,20,354,785,50,0);
        p->history=control(p,L"EDIT",L"",ES_MULTILINE|ES_READONLY|ES_AUTOVSCROLL|WS_VSCROLL,20,413,785,155,0);
        SendMessageW(p->history,EM_SETLIMITTEXT,30000,0);
        control(p,L"STATIC",L"完整逐历元结果保存位置：",0,20,579,600,26,0);
        MultiByteToWideChar(CP_ACP,0,p->files->path[PC_OUTPUT],-1,output,PC_PPP_PATH_CAP);
        control(p,L"EDIT",output,ES_READONLY|ES_AUTOHSCROLL,20,607,785,27,0);
        p->button=control(p,L"BUTTON",L"停止回放",WS_TABSTOP|BS_DEFPUSHBUTTON,665,649,140,34,STOP_ID);
        return 0;
    }
    if(msg==WM_EPOCH && p) {InterlockedExchange(&p->update_pending,0);render(p);return 0;}
    if(msg==WM_FINISHED && p) {
        render(p);
        SetWindowTextW(p->banner,p->rc==0 ? L"回放完成，结果已保存。窗口保留最后历元坐标，可关闭。" :
            p->rc==3 ? L"回放已停止，已处理历元的结果已保存。" : L"回放失败，请查看控制台和 trace；当前坐标仅为最后已处理历元。");
        SetWindowTextW(p->button,L"关闭窗口");EnableWindow(p->button,TRUE);
        if(p->close_requested) DestroyWindow(window);
        return 0;
    }
    if((msg==WM_COMMAND && LOWORD(wp)==STOP_ID) || msg==WM_CLOSE) {
        if(p && !InterlockedCompareExchange(&p->finished,0,0)) {
            if(msg==WM_CLOSE) p->close_requested=1;
            stop(p);
        }
        else DestroyWindow(window);
        return 0;
    }
    if(msg==WM_CTLCOLORSTATIC && p && (HWND)lp==p->value[VALUE_Q]) {
        int q;EnterCriticalSection(&p->lock);q=p->latest.solution.stat;LeaveCriticalSection(&p->lock);
        SetTextColor((HDC)wp,q==SOLQ_PPP ? RGB(0,110,50) : q==SOLQ_SINGLE ? RGB(150,85,0) : RGB(170,25,25));
        SetBkMode((HDC)wp,TRANSPARENT);return (LRESULT)GetSysColorBrush(COLOR_BTNFACE);
    }
    if(msg==WM_DESTROY) {PostQuitMessage(0);return 0;}
    return DefWindowProcW(window,msg,wp,lp);
}
int pc_live_view_run(const pc_ppp_files_t *files,const prcopt_t *options,
                     const char *trace_path,int realtime,int max_epochs)
{
    WNDCLASSW cls={0};live_view_t p={0};MSG msg;RECT work,r;HDC dc;int gm;
    p.files=files;p.options=options;p.trace_path=trace_path;p.realtime=realtime;p.max_epochs=max_epochs;p.rc=2;p.scale=1.0;
    InitializeCriticalSection(&p.lock);
    dc=GetDC(NULL);if(dc){p.scale=GetDeviceCaps(dc,LOGPIXELSX)/96.0;ReleaseDC(NULL,dc);}
    SystemParametersInfoW(SPI_GETWORKAREA,0,&work,0);
    if(p.scale*860>work.right-work.left) p.scale=(work.right-work.left-40)/860.0;
    if(p.scale*750>work.bottom-work.top) p.scale=(work.bottom-work.top-40)/750.0;
    cls.lpfnWndProc=window_proc;cls.hInstance=GetModuleHandleW(NULL);cls.lpszClassName=L"RTKLIB_PPP_LiveView";
    cls.hCursor=LoadCursorW(NULL,MAKEINTRESOURCEW(32512));cls.hbrBackground=(HBRUSH)(COLOR_BTNFACE+1);
    if(!RegisterClassW(&cls) && GetLastError()!=ERROR_CLASS_ALREADY_EXISTS) goto done;
    r.left=r.top=0;r.right=(int)(830*p.scale);r.bottom=(int)(704*p.scale);
    AdjustWindowRectEx(&r,WS_OVERLAPPED|WS_CAPTION|WS_SYSMENU|WS_MINIMIZEBOX,FALSE,0);
    p.window=CreateWindowExW(0,cls.lpszClassName,L"RTKLIB — 模拟实时 PPP 坐标",
        WS_OVERLAPPED|WS_CAPTION|WS_SYSMENU|WS_MINIMIZEBOX,
        work.left+(work.right-work.left-(r.right-r.left))/2,work.top+(work.bottom-work.top-(r.bottom-r.top))/2,
        r.right-r.left,r.bottom-r.top,NULL,NULL,cls.hInstance,&p);
    if(!p.window) goto done;
    ShowWindow(p.window,SW_SHOW);UpdateWindow(p.window);
    p.worker=CreateThread(NULL,0,replay_worker,&p,0,NULL);
    if(!p.worker) {DestroyWindow(p.window);goto done;}
    while((gm=GetMessageW(&msg,NULL,0,0))>0) {
        if(!IsDialogMessageW(p.window,&msg)) {TranslateMessage(&msg);DispatchMessageW(&msg);}
    }
    /* Never free the observer/UI context while the replay worker is using it. */
    if(gm<0) InterlockedExchange(&p.cancelled,1);
    WaitForSingleObject(p.worker,INFINITE);CloseHandle(p.worker);
    if(IsWindow(p.window)) DestroyWindow(p.window);
    if(gm<0) p.rc=2;
done:
    if(p.font) DeleteObject(p.font);if(p.large_font) DeleteObject(p.large_font);
    DeleteCriticalSection(&p.lock);
    return p.rc;
}
