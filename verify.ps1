$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot
wsl --cd $PWD -- bash -lc '~/venvs/fragsim/bin/python -m pytest -q'
