@ECHO OFF
"%WIN_SDK_ROOT%\%WINDOWS_SDK_VERSION%\Setup\WindowsSdkVer.exe" -q -version:%WINDOWS_SDK_VERSION%
"%WIN_SDK_ROOT%\%WINDOWS_SDK_VERSION%\Bin\SetEnv.cmd"

IF %PYTHON_ARCH% == 64 (
    appveyor DownloadFile "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4-dll64.zip"
    SET FFTW_DLL_FILENAME=fftw-3.3.4-dll64.zip
) ELSE (
    appveyor DownloadFile "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4-dll32.zip"
    SET FFTW_DLL_FILENAME=fftw-3.3.4-dll32.zip
)

call 7z.exe e %FFTW_DLL_FILENAME% -opyfftw *.dll
call 7z.exe e %FFTW_DLL_FILENAME% -opyfftw *.def
call lib /def:pyfftw\libfftw3-3.def
call lib /def:pyfftw\libfftw3f-3.def
call lib /def:pyfftw\libfftw3l-3.def