@echo off
setlocal
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0reconstruct_all_genomes.ps1" %*
exit /b %ERRORLEVEL%
