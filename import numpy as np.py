import numpy as np
import requests
import matplotlib.dates as mdates
from bs4 import BeautifulSoup
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import tqdm
import csv
import io
import warnings
import json
import math

warnings.filterwarnings("ignore")

data = {
        "date": "2022-12-17",
        "observatory":"OT",
        "minimum_priority": 2,
        "maximum_priority": 1,
        "filler": "",
}

payload = {
    "username": "observer",
    "password": "observer-2018"
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
    cos_a = (np.sin(dec) - np.sin(alt)*np.sin(lat)) / np.cos(alt)*np.cos(lat)
    a = np.arccos(cos_a)

    if np.sin(ha) < 0:
        az = a
    else:
        az = 360 - a
    return {"altitude": alt*360/(2*np.pi), "azimuth":az/(2*np.pi)}

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

with requests.Session() as s:

    day = datetime.datetime.strptime(data["date"], "%Y-%m-%d")
    j2000 = datetime.datetime(2000,1,1,12)

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
        

    if twilights["date"] == data["date"] and twilights["minimum_priority"] == data['minimum_priority']:
        targets_df = pd.read_csv("targets.csv")
        df = pd.read_csv("obs_plan.csv")

        with open("twilights.json", 'r') as openfile:
            twilights = json.load(openfile)
    else:
        print('_________________________________________________\n')
        print('Authenticating... (takes about 10-15 seconds)')
        print('_________________________________________________')
        
        p = s.post('http://research.iac.es/proyecto/muscat/users/login', data=payload)
        obs_path = 'http://research.iac.es/proyecto/muscat/observations/export'
        targets_path = 'http://research.iac.es/proyecto/muscat/stars/export'
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
        
        registration = s.post('http://research.iac.es/proyecto/muscat/stars/scheduler', data=data)
        
        mat = []  # 保存先の行列
        soup = BeautifulSoup(registration.content, 'html.parser')

        # tableの取得
        table = soup.find('table')

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
        df.to_csv("obs_plan.csv")

        with open('twilights.json', 'w') as f:
            json.dump(twilights, f, ensure_ascii=False)
 
    night_twilight = day + datetime.timedelta(hours=int(twilights["night"][0:2]),minutes=int(twilights["night"][3:5]))
    morning_twilight = day + datetime.timedelta(days=1,hours=int(twilights["morning"][0:2]),minutes=int(twilights["morning"][3:5]))

    df['Moon'] = [float(item[:-2]) for item in df['Moon']]
    df = df[df['Moon'] > 30]

    #print(df)
    ##RA Decからminimum&maximum
    ##whichever the direction in which the star is facing in the ra direction, 
    # we are only interested in the change in dec over time, because a telescope can point to any direction horizontally 
    # but not vertically

    #First, there's an easy check: find the declination of your object and the latitude of the observatory, 
    # both in degrees. Any object with a declination larger than If 90 minus the latitude of the observatory 
    # is definitely going to be visible from that observatory all the time

    # 1: Preapre RA, Dec, Time, Latitude, Longitude
    
    df['RA in deg'] = [hms_to_deg(item) for item in df['RA']]
    df['Dec in deg'] = [dms_to_deg(item) for item in df['Dec']]

    df['Ephem error TD'] = [datetime.timedelta(hours=int(item[2:4]),minutes=int(item[5:7])) for item in df['Acc period error']]

    df['Transit begin DT'] = [time_to_datetime(str(item)) for item in df["Transit begin"]]
    df['Obs begin DT'] = [begin - error - datetime.timedelta(minutes=45) for begin, error in zip(df['Transit begin DT'], df['Ephem error TD'])]

    df['Transit middle DT'] = [time_to_datetime(str(item)) for item in df["Transit middle"]]
    
    df['Transit end DT'] = [time_to_datetime(str(item)) for item in df["Transit end"]]
    df['Obs end DT'] = [end + error + datetime.timedelta(minutes=30) for end, error in zip(df['Transit end DT'], df['Ephem error TD'])]

    df['day_since_j2000'] = [timedelta_in_days(item - j2000) for item in df['Obs begin DT']]

    df['Local sidereal time'] = [local_sidereal_time(day,longitude,ut) for day, ut in zip(df['day_since_j2000'], df['Obs begin DT'])]
    df['Local sidereal time in hours'] = df['Local sidereal time']/15

    df['Hour angle'] = [hour_angle(lst,ra) for lst, ra in zip(df['Local sidereal time'], df['RA in deg'])]

    df['Alt'] = [alt_az(latitude, dec, ha)['altitude'] for dec, ha in zip(df['Dec in deg'],df['Hour angle'])]
    df['Az'] = [alt_az(latitude, dec, ha)['azimuth'] for dec, ha in zip(df['Dec in deg'],df['Hour angle'])]

    df_start = alt_at_time(df,longitude,latitude,"start")
    df_end = alt_at_time(df,longitude,latitude,"end")
    df = df[df_start > 30]
    df = df[df_end> 30]

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
    
    time = np.arange(night_twilight-datetime.timedelta(minutes=30),morning_twilight+datetime.timedelta(minutes=30), datetime.timedelta(minutes=10)).astype(datetime.datetime)
    jst  = [item + datetime.timedelta(hours=9) for item in time]

    #targets_filter = np.array([item in df['Name'].tolist() for item in targets_df['name']])

    df_sorted = df.sort_values('Obs begin DT')
    
    plans = []

    for i in range(len(df_sorted)):
        plan = []
        df_sorted = df_sorted.iloc[i:,:]
        df_next_gen = df_sorted
        while len(df_next_gen) > 0:
            #一つ目
            plan.append(df_next_gen.iloc[0])
            #時間が被っていない、一番はやく始まるやつ
            df_next_gen = df_next_gen[df_next_gen['Obs begin DT'] > df_next_gen['Obs end DT'].iloc[0]]
        plan = pd.DataFrame(plan)
        plans.append(plan)

    plans = [plan for plan in plans if len(plan) != 0]

    for plan in plans:

        fig, ax = plt.subplots(2,1,gridspec_kw={'height_ratios': [2,2]},figsize=(10,8))

        for index, object in plan.iterrows():
            meta = targets_df[targets_df["name"] == object["Name"]]
            df_altitude_plot = pd.DataFrame()

            df_altitude_plot['day_since_j2000'] = [timedelta_in_days(item - j2000) for item in time]
            df_altitude_plot['UT'] = time
            df_altitude_plot['JST'] = jst

            df_altitude_plot['Local sidereal time'] = [local_sidereal_time(day,longitude,ut) for day, ut in zip(df_altitude_plot['day_since_j2000'], time)]
            df_altitude_plot['Local sidereal time in hours'] = df_altitude_plot['Local sidereal time']/15
            df_altitude_plot['Local sidereal time DT'] = [adjust_lst(lst,day) for lst in df_altitude_plot['Local sidereal time in hours']]

            df_altitude_plot['Hour angle'] = [hour_angle(lst,object['RA in deg']) for lst in df_altitude_plot['Local sidereal time']]

            df_altitude_plot['Alt'] = [alt_az(latitude, object['Dec in deg'], ha)['altitude'] for ha in df_altitude_plot['Hour angle']]
            df_altitude_plot['Az'] = [alt_az(latitude, object['Dec in deg'], ha)['azimuth'] for ha in df_altitude_plot['Hour angle']]

            transit_duration = mdates.date2num(object['Transit end DT']) - mdates.date2num(object['Transit begin DT'])#mdates.date2num(object['Transit end DT']) - mdates.date2num(object['Transit begin DT'])
            obs_duration = mdates.date2num(object['Obs end DT']) - mdates.date2num(object['Obs begin DT'])
            transit_duration_werror = mdates.date2num(object['Transit end DT'] + object['Ephem error TD']) - mdates.date2num(object['Transit begin DT'] - object['Ephem error TD'])#mdates.date2num(object['Transit end DT']) - mdates.date2num(object['Transit begin DT'])

            color = np.random.rand(2,)
            color = np.append(color,0.3)
            transit_filter = (df_altitude_plot['UT'] > object['Transit begin DT']) & (df_altitude_plot['UT'] < object['Transit end DT'])
            altitude_filter = (df_altitude_plot['Alt'] > 0) & (df_altitude_plot['Alt'] < 90)
            intransit = df_altitude_plot[transit_filter]
            ootransit = df_altitude_plot[~transit_filter][altitude_filter]

            jst_plt = object['Transit begin DT'] + datetime.timedelta(hours=9)

            ax[0].plot(mdates.date2num(intransit['JST']), intransit['Alt'], color=color, label=object['Name'],linestyle="solid")
            ax[0].scatter(mdates.date2num(ootransit['JST']), ootransit['Alt'], color=color,s=2)

            ax[1].barh(object['Name'], left=mdates.date2num(object['Obs begin DT']), width=obs_duration, color=color,alpha=0.4,height=1,)#, left=df_altitude_plot['JST'])
            ax[1].barh(object['Name'], left=mdates.date2num(object['Transit begin DT'] - object['Ephem error TD']), width=transit_duration_werror, color=color,alpha=0.5,height=1,)#, left=df_altitude_plot['JST']) 
            ax[1].barh(object['Name'], left=mdates.date2num(object['Transit begin DT']), width=transit_duration, color=color,height=1,label=object['Name'])#, left=df_altitude_plot['JST'])
            ax[1].text(mdates.date2num(object['Transit begin DT']) + transit_duration/2, object['Name'], object['Name'], va='center' ,ha='center', fontsize=10, color='white',weight='bold')

            print('_________________________________________________')
            print(f'{object["Name"]} (Priority {object["Priority"]})')
            print(f'RA, Dec: {deg_to_hms(float(meta["RA"]))} {deg_to_dms(float(meta["Decl"]))}')
            print(f'Transit time: {object["Transit begin DT"].strftime("%H:%M")} - {object["Transit end DT"].strftime("%H:%M")} ({object["Acc period error"]})')
            print(f'Obs time: {object["Obs begin DT"].strftime("%H:%M")} - {object["Obs end DT"].strftime("%H:%M")}')
            print(f'Vmag: {float(meta["V_mag"])}')
            if type(meta["comments"].iloc[0]) != float:
                print(f'Comments: {meta["comments"].iloc[0]}')
            if type(meta["comments_sg1"].iloc[0]) != float:
                print(f'SG1 comments: {meta["comments_sg1"].iloc[0]}')

                #time.sleep(2)
        print('_________________________________________________')
            #plt.plot(time, df_altitude_plot['Alt'], color='red')    
        #plt.xlim(0,90)

        ax[0].axvline(mdates.date2num(morning_twilight + datetime.timedelta(hours=9)))
        ax[0].axvline(mdates.date2num(night_twilight + datetime.timedelta(hours=9)))
        ax[1].axvline(mdates.date2num(morning_twilight))
        ax[1].axvline(mdates.date2num(night_twilight))

        ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax[0].xaxis.set_major_locator(mdates.HourLocator(interval=1))
        ax[0].xaxis.tick_top()
        ax[0].set_xlim(mdates.date2num(jst[0]),mdates.date2num(jst[-1]))
        ax[0].set_ylim(0,90)
        ax[0].tick_params(labelbottom=False,labeltop=True)
        ax[0].set_xlabel("Time (JST)")
        ax[0].xaxis.set_label_position('top')

        ax[0].set_ylabel("Elevation")

        ax[1].xaxis.set_major_locator(mdates.HourLocator(interval=1))
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax[1].set_xlim(mdates.date2num(time[0]),mdates.date2num(time[-1]))
        #ax[1].set_xlabel("Time (UT)")

        ax[1].set_xlabel("Time (UT)")
        ax[1].set_yticks([])

        fig.tight_layout()

    plt.show()
        # note: altitude = 0 になる時間を解ける？