# observation_planner
observation planner for MuSCAT2

![observation_planner_demo](/img/observation_planner_demo.gif)

## usage:

    python observation_planner.py --obsdate=yymmdd --minp=int --all

this python script allows you to browse through upcoming transit events with MuSCAT2 to plan observations. <br/>
it can also suggest observation plans so we can save some time!<br/>

specify the observation date with --obsdate and minimum priority with --minp

### 1. retrieve airmass plots, gantt plots and sky view (eg. Dec 25 2022)

    python observation_planner.py --obsdate=221225 --minp=2 --all

![airmass_plots](/img/observation_planner_demo.gif)

### 2. suggest observation plans

    python observation_planner.py --obsdate=230120 --minp=2

![sample_plan](/img/sample_plan.png)
