# SWINGERsim installer for Windows: SWINGERsim runs inside WSL (Windows Subsystem for Linux).
#
#   powershell -ExecutionPolicy Bypass -File install.ps1 [options of install.sh]
#
# Checks that WSL with a Linux distribution is present, then runs install.sh in this folder
# inside WSL. Without WSL it prints how to install it (one command, then a restart).
param([Parameter(ValueFromRemainingArguments = $true)][string[]]$Rest)
$ErrorActionPreference = "Stop"

function Fail($msg) { Write-Host "ERROR: $msg" -ForegroundColor Red; exit 1 }

if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    Fail "WSL is not installed. Open PowerShell as administrator, run 'wsl --install', restart Windows, open 'Ubuntu' once to create a user, then run this script again."
}
$distros = (wsl.exe --list --quiet 2>$null) -replace "`0", "" | Where-Object { $_.Trim() -ne "" }
if (-not $distros) {
    Fail "WSL has no Linux distribution. Run 'wsl --install -d Ubuntu', open 'Ubuntu' once to create a user, then run this script again."
}

$here = (Get-Location).Path
$linuxPath = (wsl.exe wslpath -a "$here").Trim()
Write-Host "==> Running install.sh inside WSL in $linuxPath"
Write-Host "    (if Python 3.11+ or venv is missing, install it in Ubuntu: sudo apt install python3 python3-venv python3-dev build-essential)"
$argsLine = ($Rest | ForEach-Object { "'$_'" }) -join " "
wsl.exe bash -lc "cd '$linuxPath' && chmod +x install.sh && ./install.sh $argsLine"
if ($LASTEXITCODE -ne 0) { Fail "install.sh failed (see the messages above)" }

Write-Host ""
Write-Host "SWINGERsim is installed inside WSL. To use it, open Ubuntu (or run 'wsl' here) and type:"
Write-Host "    cd '$linuxPath'"
Write-Host "    source .venv/bin/activate"
Write-Host "    swingersim gui        # then open the printed address in your Windows browser"
Write-Host ""
Write-Host "Tip: the folder is faster inside the Linux file system, e.g. clone into ~/SWINGERsim in Ubuntu."
