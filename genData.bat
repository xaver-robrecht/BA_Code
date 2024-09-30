setlocal enabledelayedexpansion

set "rs_values=4 6 8"
set "lambda_values=1e-6 0.001 0.01 0.1 0.2"

for %%r in (%rs_values%) do (
    for %%l in (%lambda_values%) do (
        echo Running: ./bin/main.exe %%r 0.01 0 %%l
        C:\Users\Xaver\Desktop\Bachelorarbeit\Programme\bin\main.exe %%r 0.01 0 %%l
    )
)

echo All runs completed.