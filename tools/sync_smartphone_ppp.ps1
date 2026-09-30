# Verify by default. -UpdateAndroid explicitly refreshes the vendored snapshot.
param(
    [Parameter(Mandatory=$true)][string]$AndroidRoot,
    [string]$PcRoot = (Split-Path -Parent $PSScriptRoot),
    [switch]$UpdateAndroid
)
$ErrorActionPreference = 'Stop'
$PcRoot = (Resolve-Path -LiteralPath $PcRoot).Path
$AndroidRoot = (Resolve-Path -LiteralPath $AndroidRoot).Path
$sourceDir = Join-Path $PcRoot 'src'
$targetDir = Join-Path $AndroidRoot 'app/src/main/cpp/ppp_live/rtklib'
$manifestPath = Join-Path $targetDir 'pc_profile.sha256'
$records = [IO.File]::ReadAllLines($manifestPath) |
    Where-Object { $_ -match '^[0-9a-f]{64} [A-Za-z0-9_.]+$' }
if ($records.Count -lt 20) { throw 'Incomplete PPP source manifest' }
$sha = [Security.Cryptography.SHA256]::Create()
function TextHash([string]$Path) {
    $text = [IO.File]::ReadAllText($Path).Replace("`r`n", "`n").TrimEnd("`r", "`n")
    return ([BitConverter]::ToString($sha.ComputeHash(
        [Text.Encoding]::UTF8.GetBytes($text)))).Replace('-', '').ToLowerInvariant()
}
$changed = @()
$newRecords = @()
foreach ($record in $records) {
    $expected, $name = $record -split ' ', 2
    $source = Join-Path $sourceDir $name
    $target = Join-Path $targetDir $name
    $sourceHash = TextHash $source
    $targetHash = TextHash $target
    $newRecords += "$sourceHash $name"
    if ($sourceHash -cne $targetHash -or $targetHash -cne $expected) {
        $changed += $name
    }
}
if ($changed.Count -eq 0) {
    Write-Output "PPP_PROFILE_OK: $($records.Count) files; PC and Android match."
    exit 0
}
if (-not $UpdateAndroid) {
    throw "PC/Android snapshot differs: $($changed -join ', '). Use -UpdateAndroid only after validating the PC baseline."
}
# Literal filenames come from a restricted manifest, never shell expressions.
$backupDir = Join-Path $AndroidRoot (
    'tools/ppp_sync_backups/' + (Get-Date -Format 'yyyyMMdd_HHmmss_fff'))
New-Item -ItemType Directory -Path $backupDir -ErrorAction Stop | Out-Null
Copy-Item -LiteralPath $manifestPath -Destination $backupDir
foreach ($name in $changed) {
    Copy-Item -LiteralPath (Join-Path $targetDir $name) -Destination $backupDir
    Copy-Item -LiteralPath (Join-Path $sourceDir $name) -Destination (Join-Path $targetDir $name)
}
$lines = @('# PC smartphone PPP snapshot; generated from validated PC src.',
    '# SHA256: UTF-8 text, normalized CRLF, trailing newlines removed.') + $newRecords
[IO.File]::WriteAllText($manifestPath, ($lines -join "`n") + "`n",
    (New-Object Text.UTF8Encoding($false)))
Write-Output "Synced $($changed.Count) files. Previous files preserved in $backupDir"
