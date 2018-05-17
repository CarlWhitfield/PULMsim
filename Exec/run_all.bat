@echo off 
SET exec_name=PULMsim.exe

FOR %%G in (*.lung) DO (
	FOR /F "tokens=1 delims=." %%a in ('ECHO %%G') DO (
		MKDIR %%a
		MOVE %%G %%a
		COPY %exec_name% %%a
		CD %%a
		.\%exec_name% %%G
		CD ../
	)
)