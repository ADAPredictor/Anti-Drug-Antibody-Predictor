# Provision an isolated Python 3.10 runtime + mhcflurry (real MHC-I neural
# models) without touching the project's Python 3.7 application environment.
# Usage (from repo root or anywhere):
#   powershell -ExecutionPolicy Bypass -File tools\setup_mhcflurry.ps1
# Optional: $env:ADA_PIP_INDEX = "https://pypi.tuna.tsinghua.edu.cn/simple"
$ErrorActionPreference = 'Stop'

$root = Join-Path $PSScriptRoot '.'
$root = (Resolve-Path $root).Path
$zip = Join-Path $root 'py310.zip'
$pydir = Join-Path $root 'py310'
$py = Join-Path $pydir 'python.exe'

function Download-File([string]$Url, [string]$OutPath) {
    Write-Host "Downloading $Url ..."
    # Prefer curl.exe when available (more reliable on Windows); else WebClient.
    if (Get-Command curl.exe -ErrorAction SilentlyContinue) {
        & curl.exe -L --fail --retry 3 -o $OutPath $Url
        if ($LASTEXITCODE -ne 0) { throw "download failed: $Url" }
    } else {
        $wc = New-Object System.Net.WebClient
        $wc.DownloadFile($Url, $OutPath)
    }
}

New-Item -ItemType Directory -Force -Path $root | Out-Null

if (-not (Test-Path $py)) {
    Write-Host 'Downloading embeddable Python 3.10.11 ...'
    Download-File 'https://www.python.org/ftp/python/3.10.11/python-3.10.11-embed-amd64.zip' $zip
    Expand-Archive $zip $pydir -Force
    $pth = Join-Path $pydir 'python310._pth'
    if (Test-Path $pth) {
        (Get-Content $pth) -replace '#\s*import site', 'import site' | Set-Content $pth
    }
}

$getpip = Join-Path $root 'get-pip.py'
if (-not (Test-Path (Join-Path $pydir 'Scripts\pip.exe'))) {
    Write-Host 'Bootstrapping pip ...'
    Download-File 'https://bootstrap.pypa.io/get-pip.py' $getpip
    & $py $getpip --no-warn-script-location
    if ($LASTEXITCODE -ne 0) { throw 'get-pip failed' }
}

Write-Host 'Installing mhcflurry (PyTorch CPU backend, large wheel) ...'
$pipIndex = if ($env:ADA_PIP_INDEX) { @('-i', $env:ADA_PIP_INDEX) } else { @() }
& $py -m pip install --timeout 180 --retries 5 @pipIndex mhcflurry==2.2.1
if ($LASTEXITCODE -ne 0) { throw 'pip install mhcflurry failed' }

Write-Host 'Fetching mhcflurry class-I presentation models ...'
$env:MHCFLURRY_DATA_DIR = Join-Path $root 'mhcflurry_data'
New-Item -ItemType Directory -Force -Path $env:MHCFLURRY_DATA_DIR | Out-Null
& (Join-Path $pydir 'Scripts\mhcflurry-downloads.exe') fetch models_class1_presentation
if ($LASTEXITCODE -ne 0) { throw 'mhcflurry model fetch failed' }

Write-Host 'DONE'
& $py -c "import mhcflurry; print('mhcflurry', mhcflurry.__version__)"
