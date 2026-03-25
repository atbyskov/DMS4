@echo off
rem This batch file is placed in your working directory
SET ANSWAIT=1
set ANSYS_LOCK=OFF
rem set ANS_CONSEC=YES
"C:\Program files\ANSYS Inc\v251\ANSYS\bin\winx64\ansys251" -b -p ansys -smp -np 8 -i AnsoutNonlin\APDL_Nonlin_Input.txt -dir "AnsoutNonlin" -o AnsoutNonlin\AnsysOutputWindow.txt 
