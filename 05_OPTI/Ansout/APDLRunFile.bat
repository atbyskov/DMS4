@echo off
rem This batch file is placed in your working directory
SET ANSWAIT=1
set ANSYS_LOCK=OFF
rem set ANS_CONSEC=YES
"C:\Program files\ANSYS Inc\v251\ANSYS\bin\winx64\ansys251" -b -p ansys -smp -np 6 -i Ansout\APDL_Input.txt -dir "Ansout" -o Ansout\AnsysOutputWindow.txt 
