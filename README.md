# whale algorithm optimization (cmake project)

## Usage

- woa.cpp

## Input

- set the parameters of the woa agent

## Output

- the extremum of the function


## Include

- woa_config.hpp

## how to use

-set the solver parameters in the macro definition
    -CONSTRAINS stands for solving range constraints
    -OUT_PRECISION stands for output position accuracy
    -FUNCTION stands for the function to be solved，The default parameter is schaffer，eggholder，booth，matyas，cross_in_tray and levi
    -MAX_MIN_STATUS solving the maximum or minimum value，The default parameter is MAX_FIT and MIN_FIT.

-then cmake 
    ```bash
    mkdir build
    cd build
    cmake ..
    make
    ```

-run ./woa to get the results


## supplement
  For the given function，schaffer function is unstable and Levi function cannot be solved temporarily. The bug still needs to be modified

## contact
  If you have any questions, please send an email to qxy1020@163.com


