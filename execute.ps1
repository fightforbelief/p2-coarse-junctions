$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot
$py = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
& $py -VV
if ($LASTEXITCODE -ne 0) { throw "Environment audit failed" }
& $py -c "import sys; from pathlib import Path; assert Path(sys.prefix).resolve()==(Path.cwd()/'.venv').resolve(); assert sys.version_info >= (3,12)"
if ($LASTEXITCODE -ne 0) { throw "Wrong virtual environment" }
& $py -m unittest discover -v -s tests
if ($LASTEXITCODE -ne 0) { throw "Test failure" }
& $py benchmark.py
if ($LASTEXITCODE -ne 0) { throw "Benchmark failure" }
