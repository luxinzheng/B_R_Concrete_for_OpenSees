@echo off
chcp 65001 >nul
echo Building forumat.f90 and test_adina_concrete.f90 ...
gfortran -O2 -c forumat.f90 -o forumat.o
if errorlevel 1 goto err
gfortran -O2 -c test_adina_concrete.f90 -o test_adina_concrete.o
if errorlevel 1 goto err
gfortran -O2 -o test_adina_concrete.exe forumat.o test_adina_concrete.o
if errorlevel 1 goto err
echo Running tests...
test_adina_concrete.exe
echo Done. Check output *.txt files.
goto end
:err
echo Build or run failed.
exit /b 1
:end
exit /b 0
