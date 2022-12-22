# observation_planner
observation planner for MuSCAT2

## usage:

    python observation_planner.py --obsdate=yymmdd --minp=int --all

this python script allows you to browse through upcoming transit events with MuSCAT2 to plan observations. <br/>
it can also suggest observation plans so we can save some time!

### 1. retrieve airmass and gantt plots (eg. Jan 20 2023)

    python observation_planner.py --obsdate=230120 --minp=2 --all

![airmass_plots](/img/airmass_plots.png)

### 2. suggest observation plans

    python observation_planner.py --obsdate=230120 --minp=2

![sample_plan](/img/sample_plan.png)
