/* Native Windows file-selection front end. No PPP options or state here. */
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <commdlg.h>
#include <stdio.h>
#include <string.h>
#include <wchar.h>
#include "pc_file_picker.h"

#define WPATH_CAP 4096
#define EDIT_ID 200
#define BROWSE_ID 100
#define SOURCE_RINEX_ID 300
#define SOURCE_TXT_ID 301
const char *pc_file_keys[PC_FILE_COUNT] = {
    "OBS", "NAV", "SP3", "CLK", "BIA", "IONEX",
    "VMF_FIRST", "VMF_SECOND", "OROGRAPHY", "ATX", "OUTPUT"
};
static const wchar_t *labels[PC_FILE_COUNT] = {
    L"观测文件 OBS", L"广播星历 NAV", L"精密轨道 SP3", L"精密钟差 CLK",
    L"码偏差 BIA / OSB", L"电离层 IONEX", L"VMF3 前一期", L"VMF3 后一期",
    L"VMF3 格网高程", L"天线文件 ATX", L"结果保存位置 POS"
};
static const wchar_t *filters[PC_FILE_COUNT] = {
    L"RINEX观测 (*.rnx;*.obs;*.??o)\0*.rnx;*.obs;*.??o\0所有文件\0*.*\0",
    L"RINEX星历 (*.rnx;*.nav;*.??n;*.??p)\0*.rnx;*.nav;*.??n;*.??p\0所有文件\0*.*\0",
    L"精密轨道 (*.sp3;*.eph)\0*.sp3;*.eph\0所有文件\0*.*\0",
    L"精密钟差 (*.clk;*.rnx)\0*.clk;*.rnx\0所有文件\0*.*\0",
    L"偏差产品 (*.bia;*.bsx;*.dcb)\0*.bia;*.bsx;*.dcb\0所有文件\0*.*\0",
    L"IONEX (*.inx;*.ion;*.??i)\0*.inx;*.ion;*.??i\0所有文件\0*.*\0",
    L"VMF3 (*.H??;*.vmf3)\0*.H??;*.vmf3\0所有文件\0*.*\0",
    L"VMF3 (*.H??;*.vmf3)\0*.H??;*.vmf3\0所有文件\0*.*\0",
    L"格网高程 (orography*)\0orography*\0所有文件\0*.*\0",
    L"ANTEX (*.atx)\0*.atx\0所有文件\0*.*\0",
    L"定位结果 (*.pos)\0*.pos\0所有文件\0*.*\0"
};
typedef struct {
    pc_ppp_files_t *files;
    HWND edit[PC_FILE_COUNT];
    HWND label[PC_FILE_COUNT];
    HWND radio[2];
    int source;
    HFONT font;
    wchar_t settings[WPATH_CAP];
    int accepted;
    double scale;
} picker_t;

static int to_wide(const char *s, wchar_t *out, int count)
{
    return MultiByteToWideChar(CP_ACP, 0, s, -1, out, count) != 0;
}
static int to_native(const wchar_t *s, char *out, int count)
{
    BOOL substituted = FALSE;
    UINT cp = GetACP();
    int n = WideCharToMultiByte(cp, cp == CP_UTF8 ? 0 : WC_NO_BEST_FIT_CHARS,
                              s, -1, out, count, NULL,
                              cp == CP_UTF8 ? NULL : &substituted);
    return n && !substituted;
}
static int regular_file(const char *s)
{
    DWORD a = GetFileAttributesA(s);
    return a != INVALID_FILE_ATTRIBUTES && !(a & FILE_ATTRIBUTE_DIRECTORY);
}
int pc_files_validate(const pc_ppp_files_t *files, char *error, size_t capacity)
{
    char full[PC_PPP_PATH_CAP], other[PC_PPP_PATH_CAP], parent[PC_PPP_PATH_CAP];
    char *slash;
    DWORD len, attr;
    int i;
    if (error && capacity) error[0] = 0;
    if (!files) { if(error && capacity) snprintf(error,capacity,"No file configuration"); return 0; }
    if(files->source!=PC_SOURCE_RINEX && files->source!=PC_SOURCE_TXT) {
        if(error && capacity) snprintf(error,capacity,"Select exactly one observation source"); return 0;
    }
    for (i=0;i<PC_FILE_COUNT;i++) {
        if (!memchr(files->path[i],0,PC_PPP_PATH_CAP) || !files->path[i][0]) {
            if(error && capacity) snprintf(error,capacity,"%s: select a file first",pc_file_keys[i]);
            return 0;
        }
        if(i<PC_OUTPUT && !regular_file(files->path[i])) {
            if(error && capacity) snprintf(error,capacity,"%s: file does not exist: %s",pc_file_keys[i],files->path[i]);
            return 0;
        }
        if(i<PC_OUTPUT) {
            FILE *fp=fopen(files->path[i],"rb");
            if(!fp) { if(error && capacity) snprintf(error,capacity,"%s: cannot read file",pc_file_keys[i]); return 0; }
            fclose(fp);
        }
    }
    len=GetFullPathNameA(files->path[PC_OUTPUT],sizeof(full),full,NULL);
    if(!len || len>=sizeof(full) || strlen(full)+16>=PC_PPP_PATH_CAP) {
        if(error && capacity) snprintf(error,capacity,"OUTPUT: path is too long or invalid"); return 0;
    }
    for(i=0;i<PC_OUTPUT;i++) {
        len=GetFullPathNameA(files->path[i],sizeof(other),other,NULL);
        if(!len || len>=sizeof(other)) { if(error && capacity) snprintf(error,capacity,"%s: invalid path",pc_file_keys[i]); return 0; }
        if(!_stricmp(full,other)) {
            if(error && capacity) snprintf(error,capacity,"OUTPUT must not overwrite an input file (%s)",pc_file_keys[i]); return 0;
        }
        {
            char companion[PC_PPP_PATH_CAP];
            const char *suffix[]={".stat",".trace",".obs.csv"};
            int j;
            for(j=0;j<3;j++) {
                snprintf(companion,sizeof(companion),"%s%s",full,suffix[j]);
                if(!_stricmp(companion,other)) {
                    if(error && capacity) snprintf(error,capacity,"Output diagnostics must not overwrite %s",pc_file_keys[i]);return 0;
                }
            }
            strcpy(companion,full);
            { char *dot=strrchr(companion,'.'), *sep=strrchr(companion,'\\');
              if(!dot || (sep && dot<sep)) dot=companion+strlen(companion);
              memmove(dot+7,dot,strlen(dot)+1);memcpy(dot,"_events",7); }
            if(!_stricmp(companion,other)) {
                if(error && capacity) snprintf(error,capacity,"Output events must not overwrite %s",pc_file_keys[i]);return 0;
            }
        }
    }
    strcpy(parent,full); slash=strrchr(parent,'\\');
    if(!slash || !slash[1]) { if(error && capacity) snprintf(error,capacity,"OUTPUT: choose a file, not a folder"); return 0; }
    if(slash==parent+2 && parent[1]==':') slash[1]=0; else *slash=0;
    attr=GetFileAttributesA(parent);
    if(attr==INVALID_FILE_ATTRIBUTES || !(attr & FILE_ATTRIBUTE_DIRECTORY)) {
        if(error && capacity) snprintf(error,capacity,"OUTPUT: parent folder does not exist"); return 0;
    }
    attr=GetFileAttributesA(full);
    if(attr!=INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY)) {
        if(error && capacity) snprintf(error,capacity,"OUTPUT: select a file, not a folder"); return 0;
    }
    return 1;
}

static void settings_path(picker_t *p)
{
    wchar_t folder[WPATH_CAP];
    DWORD n=GetEnvironmentVariableW(L"LOCALAPPDATA",folder,WPATH_CAP);
    if(!n || n+64>=WPATH_CAP) return;
    wcscat_s(folder,WPATH_CAP,L"\\RTKLIB_PPP");
    if(!CreateDirectoryW(folder,NULL) && GetLastError()!=ERROR_ALREADY_EXISTS) return;
    swprintf_s(p->settings,WPATH_CAP,L"%ls\\file_selection.ini",folder);
}
static void load_settings(picker_t *p)
{
    int i;
    if(!p->settings[0]) return;
    p->source=GetPrivateProfileIntW(L"Files",L"Source",0,p->settings)==1 ? PC_SOURCE_TXT : PC_SOURCE_RINEX;
    for(i=0;i<PC_FILE_COUNT;i++) {
        wchar_t key[32],value[WPATH_CAP];
        to_wide(pc_file_keys[i],key,32);
        GetPrivateProfileStringW(L"Files",key,L"",value,WPATH_CAP,p->settings);
        SetWindowTextW(p->edit[i],value);
    }
}
static int save_settings(picker_t *p)
{
    HANDLE h;
    int i;
    if(!p->settings[0]) return 0;
    /* Profile APIs preserve Unicode when a new INI starts with UTF-16 BOM. */
    h=CreateFileW(p->settings,GENERIC_WRITE,0,NULL,CREATE_NEW,FILE_ATTRIBUTE_NORMAL,NULL);
    if(h!=INVALID_HANDLE_VALUE) {
        WORD bom=0xFEFF; DWORD written;
        int ok=WriteFile(h,&bom,sizeof(bom),&written,NULL) && written==sizeof(bom);
        CloseHandle(h); if(!ok) return 0;
    }
    else if(GetLastError()!=ERROR_FILE_EXISTS) return 0;
    if(!WritePrivateProfileStringW(L"Files",L"Source",
        p->source==PC_SOURCE_TXT ? L"1" : L"0",p->settings)) return 0;
    for(i=0;i<PC_FILE_COUNT;i++) {
        wchar_t key[32],value[WPATH_CAP];
        to_wide(pc_file_keys[i],key,32);
        GetWindowTextW(p->edit[i],value,WPATH_CAP);
        if(!WritePrivateProfileStringW(L"Files",key,value,p->settings)) return 0;
    }
    return 1;
}
static void suggest_output(picker_t *p)
{
    wchar_t path[WPATH_CAP], output[WPATH_CAP], *ext, *slash;
    SYSTEMTIME st;
    if(GetWindowTextLengthW(p->edit[PC_OUTPUT])) return;
    GetWindowTextW(p->edit[PC_OBS],path,WPATH_CAP);
    if(!path[0] || wcslen(path)+64>=WPATH_CAP) return;
    slash=wcsrchr(path,L'\\'); ext=wcsrchr(path,L'.');
    if(ext && (!slash || ext>slash)) *ext=0;
    GetLocalTime(&st);
    swprintf_s(output,WPATH_CAP,L"%ls_ppp_%04d%02d%02d_%02d%02d%02d.%ls",path,
               st.wYear,st.wMonth,st.wDay,st.wHour,st.wMinute,st.wSecond,
               p->source==PC_SOURCE_TXT ? L"csv" : L"pos");
    SetWindowTextW(p->edit[PC_OUTPUT],output);
}
static void browse(HWND owner,picker_t *p,int row)
{
    OPENFILENAMEW dialog={0};
    wchar_t path[WPATH_CAP];
    DWORD attr;
    if(row==PC_OUTPUT) suggest_output(p);
    GetWindowTextW(p->edit[row],path,WPATH_CAP);
    attr=GetFileAttributesW(path);
    if(row!=PC_OUTPUT && (attr==INVALID_FILE_ATTRIBUTES || (attr & FILE_ATTRIBUTE_DIRECTORY))) path[0]=0;
    dialog.lStructSize=sizeof(dialog); dialog.hwndOwner=owner;
    dialog.lpstrFilter=filters[row]; dialog.lpstrFile=path;
    dialog.nMaxFile=WPATH_CAP; dialog.lpstrTitle=labels[row];
    if(p->source==PC_SOURCE_TXT && row==PC_OBS) {
        dialog.lpstrFilter=L"Android 原始观测 (*.txt)\0*.txt\0所有文件\0*.*\0";
        dialog.lpstrTitle=L"选择 Android 原始 GNSS TXT（不是 RINEX）";
    }
    if(p->source==PC_SOURCE_TXT && row==PC_OUTPUT)
        dialog.lpstrFilter=L"逐历元定位结果 (*.csv)\0*.csv\0所有文件\0*.*\0";
    dialog.Flags=OFN_EXPLORER|OFN_NOCHANGEDIR|OFN_PATHMUSTEXIST;
    if(row==PC_OUTPUT) {dialog.Flags|=OFN_OVERWRITEPROMPT;dialog.lpstrDefExt=p->source==PC_SOURCE_TXT ? L"csv" : L"pos";}
    else dialog.Flags|=OFN_FILEMUSTEXIST;
    if(row==PC_OUTPUT ? GetSaveFileNameW(&dialog) : GetOpenFileNameW(&dialog)) {
        SetWindowTextW(p->edit[row],path);
        if(row==PC_OBS) suggest_output(p);
    }
    else if(CommDlgExtendedError())
        MessageBoxW(owner,L"文件对话框未能打开，请检查路径长度或系统环境。",L"文件选择",MB_OK|MB_ICONERROR);
}
static int output_exists(const pc_ppp_files_t *files)
{
    char name[PC_PPP_PATH_CAP];
    const char *extensions[]={"", ".stat", ".trace", ".obs.csv"};
    int i;
    for(i=0;i<4;i++) {
        snprintf(name,sizeof(name),"%s%s",files->path[PC_OUTPUT],extensions[i]);
        if(GetFileAttributesA(name)!=INVALID_FILE_ATTRIBUTES) return 1;
    }
    strcpy(name,files->path[PC_OUTPUT]);
    { char *dot=strrchr(name,'.'), *slash=strrchr(name,'\\');
      if(!dot || (slash && dot<slash)) dot=name+strlen(name);
      memmove(dot+7,dot,strlen(dot)+1);memcpy(dot,"_events",7); }
    return GetFileAttributesA(name)!=INVALID_FILE_ATTRIBUTES;
}
static int start(HWND owner,picker_t *p)
{
    pc_ppp_files_t selected={0};
    wchar_t value[WPATH_CAP],message[WPATH_CAP];
    char error[PC_PPP_PATH_CAP];
    int i;
    selected.source=p->source;
    for(i=0;i<PC_FILE_COUNT;i++) {
        GetWindowTextW(p->edit[i],value,WPATH_CAP);
        if(!to_native(value,selected.path[i],PC_PPP_PATH_CAP)) {
            MessageBoxW(owner,L"路径过长，或含当前 RTKLIB 文件接口不能表示的字符。请将文件放到较短、可表示的目录（例如英文目录）再选择。",L"路径无法使用",MB_OK|MB_ICONERROR);
            SetFocus(p->edit[i]); return 0;
        }
    }
    if(!pc_files_validate(&selected,error,sizeof(error))) {
        to_wide(error,message,WPATH_CAP);
        MessageBoxW(owner,message,L"请检查文件选择",MB_OK|MB_ICONWARNING); return 0;
    }
    if(output_exists(&selected) && MessageBoxW(owner,L"此保存位置已有结果或诊断文件。继续可能覆盖它们。\n\n建议选择新的结果文件名。仍要继续吗？",L"确认覆盖结果",MB_YESNO|MB_DEFBUTTON2|MB_ICONWARNING)!=IDYES) return 0;
    *p->files=selected;
    if(!save_settings(p)) MessageBoxW(owner,L"文件已选择，但无法记住本次选择。解算仍可继续。",L"设置保存提示",MB_OK|MB_ICONWARNING);
    p->accepted=1; DestroyWindow(owner); return 1;
}
static HWND control(HWND parent,picker_t *p,const wchar_t *type,const wchar_t *text,
                    DWORD style,int x,int y,int w,int h,int id)
{
    HWND child=CreateWindowExW(!wcscmp(type,L"EDIT")?WS_EX_CLIENTEDGE:0,type,text,
        WS_CHILD|WS_VISIBLE|style,(int)(x*p->scale),(int)(y*p->scale),
        (int)(w*p->scale),(int)(h*p->scale),parent,(HMENU)(INT_PTR)id,
        GetModuleHandleW(NULL),NULL);
    SendMessageW(child,WM_SETFONT,(WPARAM)p->font,TRUE);
    return child;
}
static void source_labels(picker_t *p)
{
    SendMessageW(p->radio[0],BM_SETCHECK,p->source==PC_SOURCE_RINEX ? BST_CHECKED : BST_UNCHECKED,0);
    SendMessageW(p->radio[1],BM_SETCHECK,p->source==PC_SOURCE_TXT ? BST_CHECKED : BST_UNCHECKED,0);
    SetWindowTextW(p->label[PC_OBS],p->source==PC_SOURCE_TXT ? L"Android 原始 TXT" : labels[PC_OBS]);
    SetWindowTextW(p->label[PC_OUTPUT],p->source==PC_SOURCE_TXT ? L"结果保存位置 CSV" : labels[PC_OUTPUT]);
}
static LRESULT CALLBACK picker_proc(HWND window,UINT msg,WPARAM wp,LPARAM lp)
{
    picker_t *p=(picker_t *)GetWindowLongPtrW(window,GWLP_USERDATA);
    if(msg==WM_NCCREATE) { p=(picker_t *)((CREATESTRUCTW *)lp)->lpCreateParams;
        SetWindowLongPtrW(window,GWLP_USERDATA,(LONG_PTR)p); }
    if(msg==WM_CREATE && p) {
        int i;
        p->font=CreateFontW((int)(-16*p->scale),0,0,0,FW_NORMAL,0,0,0,
                            DEFAULT_CHARSET,0,0,CLEARTYPE_QUALITY,0,L"Microsoft YaHei UI");
        control(window,p,L"STATIC",L"选择同一采集时段的观测和产品文件；不会自动匹配日期。PPP 算法和参数不变。",0,20,14,855,26,0);
        p->radio[0]=control(window,p,L"BUTTON",L"RINEX OBS 后处理",WS_TABSTOP|WS_GROUP|BS_AUTORADIOBUTTON,20,43,230,28,SOURCE_RINEX_ID);
        p->radio[1]=control(window,p,L"BUTTON",L"Android 原始 TXT（1×模拟实时）",WS_TABSTOP|BS_AUTORADIOBUTTON,270,43,420,28,SOURCE_TXT_ID);
        for(i=0;i<PC_FILE_COUNT;i++) {
            p->label[i]=control(window,p,L"STATIC",labels[i],0,20,91+i*38,166,26,0);
            p->edit[i]=control(window,p,L"EDIT",L"",WS_TABSTOP|ES_AUTOHSCROLL,190,87+i*38,590,28,EDIT_ID+i);
            SendMessageW(p->edit[i],EM_SETLIMITTEXT,WPATH_CAP-1,0);
            control(window,p,L"BUTTON",L"浏览…",WS_TABSTOP|BS_PUSHBUTTON,790,87+i*38,85,28,BROWSE_ID+i);
        }
        control(window,p,L"STATIC",L"观测输入二选一；切换类型会清空观测和输出路径，产品路径保留。\nTXT 按原始时间 1×回放（10分钟数据约需10分钟），不生成或读取 RINEX OBS。\nVMF3 两期应夹住观测时间。结果旁保存 stat、trace；TXT另存观测字段CSV。",0,20,519,855,65,0);
        control(window,p,L"BUTTON",L"开始解算",WS_TABSTOP|BS_DEFPUSHBUTTON,637,588,130,34,IDOK);
        control(window,p,L"BUTTON",L"取消",WS_TABSTOP|BS_PUSHBUTTON,785,588,90,34,IDCANCEL);
        settings_path(p); load_settings(p); source_labels(p);
        return 0;
    }
    if(msg==WM_COMMAND && p) {
        int id=LOWORD(wp);
        if(id==SOURCE_RINEX_ID || id==SOURCE_TXT_ID) {
            int source=id==SOURCE_TXT_ID ? PC_SOURCE_TXT : PC_SOURCE_RINEX;
            if(p->source!=source) {
                p->source=source;
                SetWindowTextW(p->edit[PC_OBS],L"");
                SetWindowTextW(p->edit[PC_OUTPUT],L"");
            }
            source_labels(p);return 0;
        }
        if(id>=BROWSE_ID && id<BROWSE_ID+PC_FILE_COUNT) {browse(window,p,id-BROWSE_ID);return 0;}
        if(id==IDOK) {start(window,p);return 0;}
        if(id==IDCANCEL) {DestroyWindow(window);return 0;}
    }
    if(msg==WM_CLOSE) {DestroyWindow(window);return 0;}
    if(msg==WM_DESTROY) {PostQuitMessage(0);return 0;}
    return DefWindowProcW(window,msg,wp,lp);
}
int pc_file_picker(pc_ppp_files_t *files,char *error,size_t capacity)
{
    WNDCLASSW cls={0}; picker_t p={0}; HWND window; MSG msg; RECT r,work;
    HDC dc; int gm;
    if(error && capacity) error[0]=0;
    p.files=files; p.scale=1.0;
    dc=GetDC(NULL);if(dc){p.scale=GetDeviceCaps(dc,LOGPIXELSX)/96.0;ReleaseDC(NULL,dc);}
    SystemParametersInfoW(SPI_GETWORKAREA,0,&work,0);
    if(p.scale*920>work.right-work.left) p.scale=(work.right-work.left-40)/920.0;
    if(p.scale*687>work.bottom-work.top) p.scale=(work.bottom-work.top-40)/687.0;
    cls.lpfnWndProc=picker_proc;cls.hInstance=GetModuleHandleW(NULL);
    cls.hCursor=LoadCursorW(NULL,MAKEINTRESOURCEW(32512));cls.hbrBackground=(HBRUSH)(COLOR_BTNFACE+1);
    cls.lpszClassName=L"RTKLIB_PPP_FilePicker";
    if(!RegisterClassW(&cls) && GetLastError()!=ERROR_CLASS_ALREADY_EXISTS) {
        if(error && capacity) snprintf(error,capacity,"Cannot register file-picker window"); return -1;
    }
    r.left=r.top=0;r.right=(int)(900*p.scale);r.bottom=(int)(642*p.scale);
    AdjustWindowRectEx(&r,WS_OVERLAPPED|WS_CAPTION|WS_SYSMENU|WS_MINIMIZEBOX,FALSE,0);
    window=CreateWindowExW(0,cls.lpszClassName,L"RTKLIB — PPP 文件选择",
        WS_OVERLAPPED|WS_CAPTION|WS_SYSMENU|WS_MINIMIZEBOX,
        work.left+(work.right-work.left-(r.right-r.left))/2,
        work.top+(work.bottom-work.top-(r.bottom-r.top))/2,
        r.right-r.left,r.bottom-r.top,NULL,NULL,cls.hInstance,&p);
    if(!window) {if(p.font)DeleteObject(p.font);if(error && capacity)snprintf(error,capacity,"Cannot create file-picker window");return -1;}
    ShowWindow(window,SW_SHOW);UpdateWindow(window);SetForegroundWindow(window);
    while((gm=GetMessageW(&msg,NULL,0,0))>0) {
        if(!IsDialogMessageW(window,&msg)) {TranslateMessage(&msg);DispatchMessageW(&msg);}
    }
    if(p.font) DeleteObject(p.font);
    if(gm<0) {if(IsWindow(window))DestroyWindow(window);if(error && capacity)snprintf(error,capacity,"File-picker message loop failed");return -1;}
    return p.accepted;
}
