param(
    [string]$Name = "rnai-designer-windows",
    [switch]$IncludeTools = $true
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
$Dist = Join-Path $Root "dist"
$Stage = Join-Path $Dist $Name
$Zip = Join-Path $Dist "$Name.zip"

if (Test-Path $Stage) { Remove-Item -LiteralPath $Stage -Recurse -Force }
if (Test-Path $Zip) { Remove-Item -LiteralPath $Zip -Force }
New-Item -ItemType Directory -Force -Path $Stage | Out-Null

$Items = @("src", "tests", "examples", "scripts", "README.md", "pyproject.toml", ".gitignore")
foreach ($Item in $Items) {
    $Source = Join-Path $Root $Item
    if (Test-Path $Source) {
        Copy-Item -LiteralPath $Source -Destination $Stage -Recurse -Force
    }
}

if ($IncludeTools -and (Test-Path (Join-Path $Root "tools"))) {
    Copy-Item -LiteralPath (Join-Path $Root "tools") -Destination $Stage -Recurse -Force
}

Compress-Archive -LiteralPath $Stage -DestinationPath $Zip -Force
Write-Host $Zip
