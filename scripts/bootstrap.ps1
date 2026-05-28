param(
    [switch]$SkipDownloads
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
$Python = Join-Path $Root ".venv\Scripts\python.exe"
$Pip = Join-Path $Root ".venv\Scripts\pip.exe"

if (-not (Test-Path $Python)) {
    py -3 -m venv (Join-Path $Root ".venv")
}

& $Pip install --upgrade pip
& $Pip install -e $Root

$Tools = Join-Path $Root "tools"
$Downloads = Join-Path $Tools "downloads"
New-Item -ItemType Directory -Force -Path $Downloads | Out-Null

if (-not $SkipDownloads) {
    $ViennaPath = Join-Path $Tools "python-packages\RNA"
    if (-not (Test-Path $ViennaPath)) {
        & $Pip install --target (Join-Path $Tools "python-packages") ViennaRNA==2.7.2
    }

    $Bowtie2Exe = Join-Path $Tools "bowtie2\bowtie2-2.5.0-mingw-x86_64\bowtie2-align-s.exe"
    if (-not (Test-Path $Bowtie2Exe)) {
        $Zip = Join-Path $Downloads "bowtie2-2.5.0-mingw-x86_64.zip"
        curl.exe -L --fail --retry 3 -o $Zip "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.0/bowtie2-2.5.0-mingw-x86_64.zip/download"
        Expand-Archive -LiteralPath $Zip -DestinationPath (Join-Path $Tools "bowtie2") -Force
    }
}

& $Python -m rnai_designer.cli --check-deps
Write-Host ""
Write-Host "Ready. Use scripts\run.ps1 to launch RNAi Designer."
