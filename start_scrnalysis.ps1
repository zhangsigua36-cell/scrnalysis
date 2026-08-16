$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$pythonPath = 'C:\Users\86135\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$nodeBin = 'C:\Users\86135\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin'

$backendListening = Get-NetTCPConnection -LocalPort 8000 -State Listen -ErrorAction SilentlyContinue
if (-not $backendListening) {
  Start-Process -FilePath 'powershell.exe' -WindowStyle Hidden -WorkingDirectory $projectRoot -ArgumentList @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-Command', "& '$pythonPath' 'backend/server.py'")
}

$frontendListening = Get-NetTCPConnection -LocalPort 3000 -State Listen -ErrorAction SilentlyContinue
if (-not $frontendListening) {
  $frontendCommand = "`$env:Path = '$nodeBin;' + `$env:Path; `$env:WRANGLER_LOG_PATH = '.wrangler/wrangler.log'; & './node_modules/.bin/vinext.CMD' dev --hostname 0.0.0.0 --port 3000"
  Start-Process -FilePath 'powershell.exe' -WindowStyle Hidden -WorkingDirectory $projectRoot -ArgumentList @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-Command', $frontendCommand)
}

Start-Sleep -Seconds 3
Write-Host 'Scrnalysis 已启动。请打开 http://localhost:3000 或本机 IP:3000。'
Write-Host '另一台电脑请使用运行 Scrnalysis 这台电脑的 IP 地址。'
