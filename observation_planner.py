import argparse
import datetime
import csv
import io
import warnings
import json
import os
import getpass
import time
import sys
import math
try:
    from astropy import units as u
    from bs4 import BeautifulSoup
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib import gridspec
    from matplotlib.animation import FuncAnimation
    import mplcursors
    import numpy as np
    import pandas as pd
    import requests
    import tqdm
except ModuleNotFoundError:
    print('_________________________________________________\n')
    print("Some required modules are not installed.\n\n    pip install -r requirements.txt\n\nPlease run the above command and try again.")
    print('_________________________________________________\n')
    sys.exit(1)

warnings.filterwarnings("ignore")
#warnings.filterwarnings('ignore', message="posx and posy should be finite values")

parser = argparse.ArgumentParser(description=\
'## obslog formatter ver. 2022 Nov. 19 ##')

parser.add_argument('--obsdate', type=int, help='observation date in yymmdd format')
parser.add_argument('--minp', type=int, help="minimum priority of objects")
parser.add_argument('--all', help='if provided, plots all visible objects', action='store_true')

args = parser.parse_args()

def obsdate_to_date(date):
    year = str(date)[0:2]
    month = str(date)[2:4]
    day = str(date)[4:6]
    return f'20{year}-{month}-{day}'

data = {
        "date": obsdate_to_date(args.obsdate),
        "observatory":"OT",
        "minimum_priority": args.minp,
        "maximum_priority": 1,
        "filler": "",
}

def dms_to_deg(deg):
    d = int(deg[0:3])
    m = int(deg[4:6])
    s = float(deg[7:])
    return d + m/60 + s/3600 

def hms_to_deg(deg):
    h = int(deg[0:2])
    m = int(deg[3:5])
    s = float(deg[6:])
    return h*15 + m*15/60 + s*15/3600

def deg_to_dms(deg):
    d = int(deg)
    m = int((deg - d) * 60)
    s = np.round((((deg - d) * 60) - m) * 60, decimals=2)
    return f'{d}:{np.abs(m)}:{np.abs(s)}'

def deg_to_hms(deg):
    d = int(deg/15)
    m = int((deg/15 - d) * 60)
    s = np.round((((deg/15 - d) * 60) - m) * 60, decimals=2)
    return f'{d}:{np.abs(m)}:{np.abs(s)}'

def replace_header(df):
    new_header = df.iloc[0] #grab the first row for the header
    df = df[1:] #get dataframe without the first row
    df.columns = new_header 
    return df

def time_to_datetime(str):
    try:
        dt = datetime.datetime.strptime(str, "%Y-%m-%d %H:%M:%S")
    except ValueError:
        dt = ""
    return dt

def timedelta_in_days(td):
    days = td.days
    seconds_in_hours = td.seconds/3600
    hours_in_days =  seconds_in_hours/24
    days = days + hours_in_days
    return days

def datetime_in_hours(dt):
    hours = dt.hour
    minutes_in_hours = dt.minute/60
    seconds_in_hours = dt.minute/3600
    return hours + minutes_in_hours + seconds_in_hours

def local_sidereal_time(d,long,ut):
    time = datetime_in_hours(ut)
    lst = 100.46 + 0.985647*d + long + 15*time
    while lst > 360:
        lst -= 360
    return lst

def hour_angle(lst,ra):
    ha = lst - ra
    while ha < 0:
        ha += 360
    return ha

def alt_az(lat,dec,ha):
    lat = lat*2*np.pi/360
    dec = dec*2*np.pi/360
    ha  = ha*2*np.pi/360
    sin_alt = np.sin(dec) * np.sin(lat) + np.cos(dec)*np.cos(lat)*np.cos(ha)
    alt = np.arcsin(sin_alt)
    cos_a = (np.sin(dec) - np.sin(alt)*np.sin(lat)) / (np.cos(alt)*np.cos(lat))
    a = np.arccos(cos_a)

    if np.sin(ha) < 0:
        az = a
    else:
        az = 2*np.pi - a
    return {"altitude": alt*360/(2*np.pi), "azimuth": az*360/(2*np.pi)}

def adjust_lst(float, date):
    if float > 12:
        lst = date + datetime.timedelta(hours=float)
    else:
        lst = date + datetime.timedelta(days=1,hours=float)
    return lst

def alt_at_time(df, longitude, latitude, mode):
    if mode == "start":
        time = df['Obs begin DT']
    elif mode == "end":
        time = df['Obs end DT']
    j2000 = datetime.datetime(2000,1,1,12)

    df['day_since_j2000'] = [timedelta_in_days(item - j2000) for item in time]

    df['Local sidereal time'] = [local_sidereal_time(day,longitude,ut) for day, ut in zip(df['day_since_j2000'], time)]
    df['Local sidereal time in hours'] = df['Local sidereal time']/15

    df['Hour angle'] = [hour_angle(lst,ra) for lst, ra in zip(df['Local sidereal time'], df['RA in deg'])]

    df['Alt'] = [alt_az(latitude, dec, ha)['altitude'] for dec, ha in zip(df['Dec in deg'],df['Hour angle'])]
    df['Az'] = [alt_az(latitude, dec, ha)['azimuth'] for dec, ha in zip(df['Dec in deg'],df['Hour angle'])]
    return df['Alt']

def moon_position(longitude,latitude,lst):
    # Get the Moon's position in the ICRS (International Celestial Reference System) frame
    julian_date = (lst - longitude + (2451545.0 * 360.98564736629) - 280.46061837)/360.98564736629 
    ecliptic_longitude = (218.316 + 13.176396 * julian_date)*2*np.pi/360
    ecliptic_latitude = 5.13 * np.sin(ecliptic_longitude*2*np.pi/360)
    right_ascension = math.degrees(math.atan2(np.cos(ecliptic_longitude*2*np.pi/360) * np.sin(ecliptic_latitude*2*np.pi/360), np.cos(ecliptic_latitude*2*np.pi/360)))
    declination = math.degrees(math.asin(math.sin(ecliptic_longitude*2*np.pi/360) * math.sin(ecliptic_latitude*2*np.pi/360)))
    return alt_az(latitude,declination,lst-right_ascension)

def moon_separation(alt,az,moon_alt, moon_az):
    separation = np.sin(moon_alt*2*np.pi/360)*np.sin(alt*2*np.pi/360) + np.cos(moon_alt*2*np.pi/360)*np.cos(alt*2*np.pi/360)*np.cos(np.abs((az*2*np.pi/360) - moon_az*2*np.pi/360))
    return np.arccos(separation)*360/(2*np.pi)

def darken_color(color,value):
    color = color - value
    color = [item if item > 0 else 0 for item in color]
    color = [item if item < 1 else 1 for item in color]
    return np.array(color)

with requests.Session() as s:

    day = datetime.datetime.strptime(data["date"], "%Y-%m-%d")
    j2000 = datetime.datetime(2000,1,1,12)
    resolution = 10

    latitude, longitude = 28.3, -16.5097

    try:
        with open("twilights.json", 'r') as openfile:
            twilights = json.load(openfile)
    except FileNotFoundError:
        twilights ={
            "date": "",
            "minimum_priority": "",
            "night": "",
            "morning": ""
        }
        with open('twilights.json', 'w') as f:
            json.dump(twilights, f, ensure_ascii=False)   

    if twilights["date"] == data["date"] and twilights["minimum_priority"] >= data['minimum_priority']:
        targets_df = pd.read_csv("targets.csv")
        df = pd.read_csv("obs_plan.csv")
        df = df[df["Priority"] <= data['minimum_priority']]

        with open("twilights.json", 'r') as openfile:
            twilights = json.load(openfile)
    else:
        dirname = os.path.dirname(__file__)
        credentials = os.path.join(dirname, 'cred.json')

        try: 
            with open(credentials, 'r') as openfile:
                payload = json.load(openfile)
        except FileNotFoundError:
            print('\n_____MuSCAT2_wiki_login__________________________\n')
            print("***you can also provide a cred.json file to bypass login***\n")
            username = input('username: ')
            password = getpass.getpass(prompt='password: ')
            print('_________________________________________________\n')
            payload = {"username": username, "password": password}

        print('_________________________________________________\n')
        print('Authenticating... (takes about 10-15 seconds)')
        print('_________________________________________________')
        
        time_start = time.time()
        p = s.post('https://research.iac.es/proyecto/muscat/users/login', data=payload)
        time_end = time.time()
        elapsed = time_end - time_start
        if elapsed < 2:
            print("\nlogin failed: wrong username/password\n")
            sys.exit(1)
        obs_path = 'https://research.iac.es/proyecto/muscat/observations/export'
        targets_path = 'https://research.iac.es/proyecto/muscat/stars/export'
        path_list = [obs_path,targets_path]
        r = [s.get(path) for path in tqdm.tqdm(path_list, desc=f'Downloading past observation and registered targets data... (about 15 seconds)')]
        obs_text = r[0].text#.encode('utf-8')
        targets_text = r[1].text
        #column = ["id","telescope","observer","weather","seeing","temperature","start_time","end_time","flats","bias","exposure_time","lightcurve","quicklook","comments","files","star_id","code"]
        obs_reader = csv.reader(io.StringIO(obs_text), delimiter=';',quotechar='"',
                            quoting=csv.QUOTE_ALL, skipinitialspace=True)#,index_col=column)
        targets_reader = csv.reader(io.StringIO(targets_text), delimiter=';',quotechar='"',
                            quoting=csv.QUOTE_ALL, skipinitialspace=True)
        observations_df = pd.DataFrame([row for row in obs_reader])
        targets_df = pd.DataFrame([row for row in targets_reader])
        targets_df = replace_header(targets_df)

        #print(registration.content)
        print('_________________________________________________\n')
        print(f'Generating airmass plots for {day.strftime("%B %d, %Y")}')
        print('_________________________________________________')
        
        registration = s.post('https://research.iac.es/proyecto/muscat/stars/scheduler', data=data)
        
        mat = []  # 保存先の行列
        soup = BeautifulSoup(registration.content, 'html.parser')

        # tableの取得
        special_targets = []
        table = soup.find('table')
        lists = soup.find_all("ul")
        special_observations = lists[-1].find_all("li")
        for item in special_observations:
            extracted_date = item.text.split(" ")[:5]
            special_date = [extracted_date[0]] + extracted_date[3:]
            special_date = datetime.datetime.strptime(" ".join(special_date),"%Y %d %B")
            if day == special_date:
                special_targets.append(item.text)
        
        if len(special_targets) != 0:
            for item in special_targets:
                print(f'\n{item}')
            input(f'\nWARNING: The above special observations are detected on {day.strftime("%B %d, %Y")}. Press y to continue: ')

        # theadの解析
        r = []  # 保存先の行
        thead = table.find('thead')  # theadタグを探す
        ths = thead.tr.find_all('th')
        for th in ths:  # thead -> trからthタグを探す
            r.append(th.text)  # thタグのテキストを保存

        mat.append(r)  # 行をテーブルに保存
        # tbodyの解析
        tbody = table.find('tbody')  # tbodyタグを探す
        trs = tbody.find_all('tr')  # tbodyからtrタグを探す
        for tr in trs:
            r = []  # 保存先の行
            for td in tr.find_all('td'):  # trタグからtdタグを探す
                r.append(td.text)  # tdタグのテキストを保存
            mat.append(r)
        # 出力
        #for r in mat:
            #print(','.join(r))  # カンマ（,）で列を結合して表示
        df = pd.DataFrame(mat)
        df = replace_header(df)

        night_info = soup.find('div', class_='night_info_box').text
        night_twilight = night_info.find("Night twilight nautical: ")
        morning_twilight = night_info.find("Morning twilight nautical: ")

        night_twilight = night_info[night_twilight + len("Night twilight nautical: "): night_twilight + len("Night twilight nautical: ")+5]
        morning_twilight = night_info[morning_twilight + len("Morning twilight nautical: "):morning_twilight + len("Morning twilight nautical: ")+5]

        twilights ={
            "date": data["date"],
            "minimum_priority": data["minimum_priority"],
            "night": night_twilight,
            "morning": morning_twilight
        }

        df['Acc period error'] = ["+ 00:00:00" if item == "" else item for item in df['Acc period error']]

        targets_df.to_csv("targets.csv")
        ra_list = []
        dec_list = []
        for index, row in df.iterrows():
            print(df)
            print(row)
            meta = targets_df[targets_df["name"] == row["Name"]]
            ra_list.append(float(meta['RA']))
            dec_list.append(float(meta['Decl']))
        print(ra_list)
        print(dec_list)
        df['RA'] = ra_list
        df['Dec'] = dec_list

        df.to_csv("obs_plan.csv")

        with open('twilights.json', 'w') as f:
            json.dump(twilights, f, ensure_ascii=False)
 
    night_twilight = day + datetime.timedelta(hours=int(twilights["night"][0:2]),minutes=int(twilights["night"][3:5]))
    morning_twilight = day + datetime.timedelta(days=1,hours=int(twilights["morning"][0:2]),minutes=int(twilights["morning"][3:5]))

    df['Moon'] = [float(item[:-2]) for item in df['Moon']]
    df = df[df['Moon'] > 30]

    ##RA Decからminimum&maximum
    ##whichever the direction in which the star is facing in the ra direction, 
    # we are only interested in the change in dec over time, because a telescope can point to any direction horizontally 
    # but not vertically

    #First, there's an easy check: find the declination of your object and the latitude of the observatory, 
    # both in degrees. Any object with a declination larger than If 90 minus the latitude of the observatory 
    # is definitely going to be visible from that observatory all the time

    # 1: Preapre RA, Dec, Time, Latitude, Longitude
    '''
    df['RA in deg'] = [hms_to_deg(item) for item in df['RA']]
    df['Dec in deg'] = [dms_to_deg(item) for item in df['Dec']]
    '''
    df['RA in deg'] = df['RA']
    df['Dec in deg'] = df['Dec']

    df['Ephem error TD'] = [datetime.timedelta(hours=int(item[2:4]),minutes=int(item[5:7])) for item in df['Acc period error']]

    df['Transit begin DT'] = [time_to_datetime(str(item)) for item in df["Transit begin"]]
    df['Obs begin DT'] = [begin - error - datetime.timedelta(minutes=75) for begin, error in zip(df['Transit begin DT'], df['Ephem error TD'])]

    df['Transit middle DT'] = [time_to_datetime(str(item)) for item in df["Transit middle"]]
    
    df['Transit end DT'] = [time_to_datetime(str(item)) for item in df["Transit end"]]
    df['Obs end DT'] = [end + error + datetime.timedelta(minutes=60) for end, error in zip(df['Transit end DT'], df['Ephem error TD'])]

    df['day_since_j2000'] = [timedelta_in_days(item - j2000) for item in df['Obs begin DT']]

    df['Local sidereal time'] = [local_sidereal_time(day,longitude,ut) for day, ut in zip(df['day_since_j2000'], df['Obs begin DT'])]
    df['Local sidereal time in hours'] = df['Local sidereal time']/15

    df['Hour angle'] = [hour_angle(lst,ra) for lst, ra in zip(df['Local sidereal time'], df['RA in deg'])]

    df['Alt'] = [alt_az(latitude, dec, ha)['altitude'] for dec, ha in zip(df['Dec in deg'],df['Hour angle'])]
    df['Az'] = [alt_az(latitude, dec, ha)['azimuth'] for dec, ha in zip(df['Dec in deg'],df['Hour angle'])]

    df_start = alt_at_time(df,longitude,latitude,"start")
    df_end = alt_at_time(df,longitude,latitude,"end")
    #df = df[df_start > 30]
    #df = df[df_end> 30]

    filler = df[~df['Filler'].str.contains("No")]

    #df = df[df['Obs begin DT'] > night_twilight]
    #df = df[df['Obs end DT'] < morning_twilight]

    ##ここに最適化アルゴリズムを噛ませて、それ以外の天体をフィルタリングする df[df['Name'] in {アルゴリズムの出力}]]

    #print(df[['Alt','Az']])

    # 2: Convert them all to angles except for time
    # 3: Convert day with respect to J2000
    # 4: Calculate the local sidereal time in degrees: (100.46 * 0.985647*d + longitude(deg) + 15*UT(hours))
    # 5: Calculate the hour angle: LST - RA

    #上のコードは全ての天体のトランジット開始時刻でのaltitudeの計算
    #ここからは全ての天体について任意の時刻でできるように書き直し

    if not args.all:
        df_sorted = df.sort_values('Obs begin DT')

        plans = []

        for i in range(len(df_sorted)):
            plan = []
            df_next_gen = df_sorted.iloc[i:,:]
            #print(df_next_gen)
            while len(df_next_gen) > 0:
                #一つ目
                plan.append(df_next_gen.iloc[0])
                #時間が被っていない
                df_next_gen = df_next_gen[df_next_gen['Obs begin DT'] > df_next_gen['Obs end DT'].iloc[0]]
                #一番はやく終わるやつ
                df_next_gen = df_next_gen.sort_values('Obs end DT')
            plan = pd.DataFrame(plan)
            plan = pd.concat([plan, filler], sort=False)
            plans.append(plan)

        #print([len(plan) for plan in plans])
        plans = [plan for plan in plans if (len(plan) != 0)]
        plans = [plan for plan in plans if (len(plan) != 1)]
        if len(plans) > 9:
            plans = plans[:9]
    else:
        #plans = [df.sort_values('Obs begin DT',ascending=False)]
        plans = [df]

    constants = pd.DataFrame()
    constants['UT'] = np.arange(night_twilight-datetime.timedelta(minutes=30),morning_twilight+datetime.timedelta(minutes=30), datetime.timedelta(minutes=resolution)).astype(datetime.datetime)
    constants['JST'] = [item + datetime.timedelta(hours=9) for item in constants['UT']]
    constants['day_since_j2000'] = [timedelta_in_days(item - j2000) for item in constants['UT']]

    constants['Local sidereal time'] = [local_sidereal_time(day,longitude,ut) for day, ut in zip(constants['day_since_j2000'], constants['UT'])]
    constants['Local sidereal time in hours'] = constants['Local sidereal time']/15
    constants['Local sidereal time DT'] = [adjust_lst(lst,day) for lst in constants['Local sidereal time in hours']]
    constants['Moon altitude'] = [moon_position(longitude,latitude,lst)['altitude'] for lst in constants['Local sidereal time']]
    constants['Moon azimuth'] = [moon_position(longitude,latitude,lst)['azimuth'] for lst in constants['Local sidereal time']]

    for i, plan in enumerate(plans):
        plan = plan.sort_values(['Filler','Priority','Name'],ascending=[False,False,False])

        fig = plt.figure(figsize=(15,8))
        gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1]) 
        ax_airmass_plot = plt.subplot(gs[0,0])
        ax_gantt_plot   = plt.subplot(gs[1,0])
        ax_polar_plot   = plt.subplot(gs[0,1],polar=True)

        ax_airmass_plot.axvline(mdates.date2num(morning_twilight + datetime.timedelta(hours=9)),alpha=0.5)
        ax_airmass_plot.axvline(mdates.date2num(night_twilight + datetime.timedelta(hours=9)),alpha=0.5)
        ax_gantt_plot.axvline(mdates.date2num(morning_twilight),alpha=0.5)
        ax_gantt_plot.axvline(mdates.date2num(night_twilight),alpha=0.5)

        print(f'______Plan {i+1}/{len(plans)}___________________________________')

        #print(constants['Moon position'])
        #ax_airmass_plot.scatter(mdates.date2num(constants['JST']), constants['Moon altitude'], marker='D',color='black',s=2)
        ax_polar_plot.bar(np.linspace(0,360,50),1,bottom=np.cos(30*2*np.pi/360), color='pink', alpha=0.4)
        #ax_polar_plot.scatter(constants['Moon azimuth']*2*np.pi/360, np.cos(constants['Moon altitude']*2*np.pi/360), marker='D',color='black',s=2)

        object_df_list = []
        object_info_list = []
        altitude_plot_list = []
        altitude_marker_list = []
        observation_plot_list = []
        moon_separation_list = []
        polar_plot_list = []
        polar_plot_trajectory_list = []

        #print(plans[0].columns)

        for index, object in plan.iterrows():
            meta = targets_df[targets_df["name"] == object["Name"]]
            color = np.random.uniform(low=0.42, high=0.95, size=(3,))
            text_color = darken_color(color, 0.3)

            df_altitude_plot = pd.DataFrame()

            df_altitude_plot['UT'] = constants['UT']
            df_altitude_plot['JST'] = constants['JST']

            df_altitude_plot['Hour angle'] = [hour_angle(lst,object['RA in deg']) for lst in constants['Local sidereal time']]

            df_altitude_plot['Alt'] = [alt_az(latitude, object['Dec in deg'], ha)['altitude'] for ha in df_altitude_plot['Hour angle']]
            df_altitude_plot['Az'] = [alt_az(latitude, object['Dec in deg'], ha)['azimuth'] for ha in df_altitude_plot['Hour angle']]

            df_altitude_plot['Moon separation'] = [moon_separation(alt,az,moon_alt,moon_az) for alt,az,moon_alt, moon_az in zip(df_altitude_plot['Alt'],df_altitude_plot['Az'],constants['Moon altitude'],constants['Moon azimuth'])]

            try:
                object_info = f'{object["Name"]} (Priority {object["Priority"]})\nRA, Dec: {deg_to_hms(float(meta["RA"]))} {deg_to_dms(float(meta["Decl"]))}\nTransit time: {object["Transit begin DT"].strftime("%H:%M")} - {object["Transit end DT"].strftime("%H:%M")} ({object["Acc period error"][0:7]})\nObs time: {object["Obs begin DT"].strftime("%H:%M")} - {object["Obs end DT"].strftime("%H:%M")}\nMoon: {object["Moon"]} (min)\nVmag: {np.round(float(meta["V_mag"]),1) if meta["V_mag"].iloc[0] != "" else "N/A"}\nComments: {meta["comments"].iloc[0][:30] if type(meta["comments"].iloc[0]) != float else "None"}\n                  {meta["comments"].iloc[0][30:60] + " ..." if type(meta["comments"].iloc[0]) != float and len(meta["comments"].iloc[0]) > 29 else ""}\n\n\n\n\n\n'
            except ValueError:
                object_info = f'{object["Name"]} (Priority {object["Priority"]}|{object["Filler"]})\nRA, Dec: {deg_to_hms(float(meta["RA"]))} {deg_to_dms(float(meta["Decl"]))}\nMoon: {object["Moon"]} (min)\nVmag: {np.round(float(meta["V_mag"]),1) if meta["V_mag"].iloc[0] != "" else "N/A"}\nComments: {meta["comments"].iloc[0][:30] if type(meta["comments"].iloc[0]) != float else "None"}\n                  {meta["comments"].iloc[0][30:60] + " ..." if type(meta["comments"].iloc[0]) != float and len(meta["comments"].iloc[0]) > 29 else ""}\n\n\n\n\n\n'

            transit_filter = (df_altitude_plot['UT'] > object['Transit begin DT']) & (df_altitude_plot['UT'] < object['Transit end DT'])
            altitude_filter = (df_altitude_plot['Alt'] > 0) & (df_altitude_plot['Alt'] < 90)
            obs_lim_filter = (df_altitude_plot['Alt'] > 30) & (df_altitude_plot['Alt'] < 90)

            intransit = df_altitude_plot[transit_filter]
            ootransit = df_altitude_plot[~transit_filter][altitude_filter]

            altitude_plot, = ax_airmass_plot.plot(mdates.date2num(df_altitude_plot['JST']), df_altitude_plot['Alt'], color=color, alpha=0.)
            ax_airmass_plot.plot(mdates.date2num(intransit['JST']), intransit['Alt'], color=color, label=object['Name'],linestyle="solid")
            ax_airmass_plot.scatter(mdates.date2num(ootransit['JST']), ootransit['Alt'], color=color,s=2)
            altitude_marker, = ax_airmass_plot.plot(mdates.date2num(df_altitude_plot['JST']), np.full(len(df_altitude_plot),0) , color="pink",alpha=0.0,)#, left=df_altitude_plot['JST'])

            obs_max_duration = mdates.date2num(df_altitude_plot['UT'][obs_lim_filter].iloc[0]) - mdates.date2num(df_altitude_plot['UT'][obs_lim_filter].iloc[-1])

            if type(object['Obs begin DT']) != pd._libs.tslibs.nattype.NaTType:
                transit_duration = mdates.date2num(object['Transit end DT']) - mdates.date2num(object['Transit begin DT'])#mdates.date2num(object['Transit end DT']) - mdates.date2num(object['Transit begin DT'])
                obs_duration = mdates.date2num(object['Obs end DT']) - mdates.date2num(object['Obs begin DT'])
                transit_duration_werror = mdates.date2num(object['Transit end DT'] + object['Ephem error TD']) - mdates.date2num(object['Transit begin DT'] - object['Ephem error TD'])#mdates.date2num(object['Transit end DT']) - mdates.date2num(object['Transit begin DT'])

                ax_gantt_plot.barh(object['Name'], left=mdates.date2num(object['Obs begin DT']), width=obs_duration, color=color,alpha=0.4,height=1)#, left=df_altitude_plot['JST'])
                ax_gantt_plot.barh(object['Name'], left=mdates.date2num(object['Transit begin DT'] - object['Ephem error TD']), width=transit_duration_werror, color=color, alpha=0.5, height=1,)#, left=df_altitude_plot['JST']) 
                observation_plot, = ax_gantt_plot.barh(object['Name'], left=mdates.date2num(object['Transit begin DT']), width=transit_duration, color=color,height=1)#, left=df_altitude_plot['JST'])
                ax_gantt_plot.text(mdates.date2num(object['Transit begin DT']) + transit_duration/2, object['Name'], f'{object["Name"]} [{str(object["Priority"])}]', va='center' ,ha='center', fontsize=10, color=text_color,weight='bold')
            else:
                observation_plot, = ax_gantt_plot.barh(object['Name'], left=mdates.date2num(df_altitude_plot['UT'][obs_lim_filter].iloc[-1]), width=obs_max_duration, color="gray",height=0.5,alpha=0.5)#, left=df_altitude_plot['JST'])
                ax_gantt_plot.text(mdates.date2num(df_altitude_plot['UT'][obs_lim_filter].iloc[-1]) + obs_max_duration/2, object['Name'], f'{object["Name"]} [{str(object["Priority"])}|{object["Filler"]}]', va='center' ,ha='center', fontsize=10, color=darken_color(text_color,0.05),weight='bold')
            ax_gantt_plot.hlines(object['Name'], df_altitude_plot['UT'][obs_lim_filter].iloc[0],df_altitude_plot['UT'][obs_lim_filter].iloc[-1] , color=color,alpha=0.2)#, left=df_altitude_plot['JST'])
            print(object["Name"])
            print(df_altitude_plot)
            #sys.exit(1)

            polar_plot_trajectory, = ax_polar_plot.plot(df_altitude_plot['Az'][obs_lim_filter]*2*np.pi/360 + (np.pi/2),np.cos(df_altitude_plot['Alt'][obs_lim_filter]*2*np.pi/360),color=color,alpha=0,linestyle="dotted")
            polar_plot, = ax_polar_plot.plot(df_altitude_plot['Az'][obs_lim_filter]*2*np.pi/360 + (np.pi/2),np.cos(df_altitude_plot['Alt'][obs_lim_filter]*2*np.pi/360), color=color,alpha=0)

            object_df_list.append([df_altitude_plot[['Alt','Az']],object["Name"],text_color])
            altitude_plot_list.append(altitude_plot)
            altitude_marker_list.append(altitude_marker)
            observation_plot_list.append(observation_plot)
            object_info_list.append(object_info)
            moon_separation_list.append(df_altitude_plot['Moon separation'])
            polar_plot_list.append(polar_plot)
            polar_plot_trajectory_list.append(polar_plot_trajectory)

            print(f'\n{object["Name"]} (Priority {object["Priority"]})')
            print(f'RA, Dec: {deg_to_hms(float(meta["RA"]))} {deg_to_dms(float(meta["Decl"]))}')
            try:
                print(f'Transit time: {object["Transit begin DT"].strftime("%H:%M")} - {object["Transit end DT"].strftime("%H:%M")} ({object["Acc period error"][0:7]})')
                print(f'Obs time: {object["Obs begin DT"].strftime("%H:%M")} - {object["Obs end DT"].strftime("%H:%M")}')
            except ValueError:
                pass
            if meta["V_mag"].iloc[0] != '':
                print(f'Vmag: {float(meta["V_mag"])}')
            else:
                print(f'Vmag: N/A')
            if type(meta["comments"].iloc[0]) != float:
                print(f'Comments: {meta["comments"].iloc[0]}')
            if type(meta["comments_sg1"].iloc[0]) != float:
                print(f'SG1 comments: {meta["comments_sg1"].iloc[0]}')
            
        print('_________________________________________________')

        cursor = mplcursors.cursor(
                observation_plot_list,
                hover=True,  # Transient
                annotation_kwargs=dict(
                    bbox=dict(
                        boxstyle="square,pad=0.5",
                        facecolor="white",
                        edgecolor="#ddd",
                        linewidth=0.,
                    ),
                    linespacing=1.5,
                    arrowprops=None,
                ),
                highlight=True,
                highlight_kwargs=dict(alpha=0.5),
            )

        pairs = dict(zip(observation_plot_list, altitude_plot_list))
        pairs.update(zip(observation_plot_list,altitude_plot_list))
        pairs_2 = dict(zip(observation_plot_list, object_info_list))
        pairs_2.update(zip(observation_plot_list,object_info_list))
        pairs_3 = dict(zip(observation_plot_list, polar_plot_list))
        pairs_3.update(zip(observation_plot_list,polar_plot_list))
        pairs_4 = dict(zip(observation_plot_list, polar_plot_trajectory_list))
        pairs_4.update(zip(observation_plot_list,polar_plot_trajectory_list))
        pairs_5 = dict(zip(observation_plot_list, altitude_marker_list))
        pairs_5.update(zip(observation_plot_list,altitude_marker_list))

        @cursor.connect("add")
        def on_add(sel):
            sel.annotation.set_text(pairs_2[sel.artist])
            sel.annotation.set(position=(mdates.date2num(constants['UT'].iloc[-1]+datetime.timedelta(hours=0.5)), 0))
            sel.extras.append(cursor.add_highlight(pairs[sel.artist]))
            sel.extras.append(cursor.add_highlight(pairs_3[sel.artist]))
            sel.extras.append(cursor.add_highlight(pairs_4[sel.artist]))
            sel.extras.append(cursor.add_highlight(pairs_5[sel.artist]))

        live_plot = ax_polar_plot.scatter([],[])
        live_time_jst = ax_airmass_plot.axvline(0)
        live_time_ut = ax_gantt_plot.axvline(0)
        live_annotations = [ax_polar_plot.annotate("",[0,0]) for item in plan]
        
        def animate_planets(i):
            now = datetime.datetime.utcnow()
            past_midnight = day + datetime.timedelta(days=1)
            today = (day.year == now.year and day.month == now.month and day.day == now.day and now.hour > 12) or (past_midnight.year == now.year and past_midnight.month == now.month and past_midnight.day == now.day and now.hour < 12)
            if today:
                pass
            else:
                now = constants['UT'].iloc[0] + datetime.timedelta(minutes=(i*5))
            if constants['UT'].iloc[0] < now and constants['UT'].iloc[-1] > now:
                live_time_jst.set_xdata([mdates.date2num(now+datetime.timedelta(hours=9))])# = ax_airmass_plot.axvline(mdates.date2num(now+datetime.timedelta(hours=9)))
                live_time_ut.set_xdata([mdates.date2num(now)])# = ax_gantt_plot.axvline(mdates.date2num(now))
                live_time_ut.set_color("pink")
                live_time_jst.set_color("pink")
                if (i==0 or i%(resolution*60) == 0) or (not today):
                    current_time = constants['UT'][constants['UT'] > now].index[0]
                    live_alt = [df[0]['Alt'].iloc[current_time] *2*np.pi/360 if df[0]['Alt'].iloc[current_time] > 0 else -90 for df in object_df_list] 
                    live_az = [df[0]['Az'].iloc[current_time] *2*np.pi/360 + (np.pi/2) for df in object_df_list]
                    live_color = [df[2] for df in object_df_list]
                    live_name = [df[1] for df in object_df_list]
                    live_plot.set_offsets(np.transpose(np.asarray((live_az, np.cos(live_alt)))))#, color=live_color,s=2)
                    live_plot.set_color(live_color)
                    for alt, altitude_marker in zip(live_alt,altitude_marker_list):
                        altitude_marker.set_data([mdates.date2num(constants['JST'].iloc[0]),mdates.date2num(constants['JST'].iloc[-1])],[alt*360/(2*np.pi),alt*360/(2*np.pi)])
                    for live_annotation, text, alt, az, color in zip(live_annotations,live_name,live_alt,live_az,live_color):
                            live_annotation.set_position((az,np.cos(alt)))
                            live_annotation.set_text(text)
                            live_annotation.set_color(color)
                    for polar_plot, polar_plot_trajectory, df in zip(polar_plot_list,polar_plot_trajectory_list, object_df_list):
                        trajectory = df[0][current_time+1:][df[0]['Alt'] > 0]*2*np.pi/360
                        path = df[0][:current_time][df[0]['Alt'] > 0]*2*np.pi/360
                        polar_plot_trajectory.set_data(trajectory['Az'] + (np.pi/2),np.cos(trajectory['Alt']))#, color=live_color,s=2)
                        polar_plot.set_data(path['Az'] + (np.pi/2),np.cos(path['Alt']))                 
        animation = FuncAnimation(fig, animate_planets, interval = 1000,)
        
        #ax_airmass_plot.plot(np.linspace(mdates.date2num(constants['JST'].iloc[0]),mdates.date2num(constants['JST'].iloc[-1]),100),np.full(100,30),color="red", linestyle="dashed",alpha=0.2)
        ax_airmass_plot.fill_between(np.linspace(mdates.date2num(constants['JST'].iloc[0]),mdates.date2num(constants['JST'].iloc[-1]),100), 0, 30,color='pink', alpha=0.7)
        ax_airmass_plot.set_title(f'Observation Plan {i+1}/{len(plans)} on {twilights["date"]}')
        ax_airmass_plot.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax_airmass_plot.xaxis.set_major_locator(mdates.HourLocator(interval=1))
        ax_airmass_plot.xaxis.tick_top()
        ax_airmass_plot.set_xlim(mdates.date2num(constants['JST'].iloc[0]),mdates.date2num(constants['JST'].iloc[-1]))
        ax_airmass_plot.set_ylim(0,90)
        ax_airmass_plot.tick_params(labelbottom=False,labeltop=True)
        ax_airmass_plot.set_xlabel("Time (JST)")
        ax_airmass_plot.xaxis.set_label_position('top')
        ax_airmass_plot.set_ylabel("Elevation")
        ax_airmass_plot.set_yticks(np.linspace(0,90,10))
        ax_airmass_plot.set_yticklabels([f'${int(angle)}^\circ$' for angle in np.linspace(0,90,10)])

        ax_gantt_plot.xaxis.set_major_locator(mdates.HourLocator(interval=1))
        ax_gantt_plot.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax_gantt_plot.set_xlim(mdates.date2num(constants['UT'].iloc[0]),mdates.date2num(constants['UT'].iloc[-1]))

        ax_gantt_plot.set_xlabel("Time (UT)")
        ax_gantt_plot.set_yticks([])
        
        ax_polar_plot.set_xticklabels(["E", "NE",f'N (0$^\circ$)', "NW", "W", "SW", "S", "SE", ])
        ax_polar_plot.set_ylim(0,1)

        ax_polar_plot.set_yticks([np.cos(0*2*np.pi/360),np.cos(30*2*np.pi/360),np.cos(60*2*np.pi/360),np.cos(90*2*np.pi/360)])
        ax_polar_plot.set_yticklabels([f'$0^\circ$',f'$30^\circ$',f'$60^\circ$',f'$90^\circ$'],color="gray")
        ax_polar_plot.set_title(f'Sky view')
        
        fig.tight_layout()
        plt.subplots_adjust(hspace=0.)
        plt.show()
        plt.savefig("image.png")
        # note: altitude = 0 になる時間を解ける？
