## import modules

import numpy as np
import math
import pandas as pd
import matplotlib
from datetime import datetime, timezone
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import dates, rcParams, colormaps, colors
import colorcet as cc
from obspy.clients.fdsn import Client
from obspy import Stream, Trace, UTCDateTime, read
from compute_metrics import compute_metrics, process_waveform
import sys
from waveform_collection import gather_waveforms
from scipy.signal import spectrogram, medfilt
from matplotlib.patches import Rectangle
from scipy import stats
import datetime
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.dates import DateFormatter, AutoDateLocator
import matplotlib.dates as mdates
import pytz

# user input
# TOP_DIR = '/Users/dfee/Documents/AVO/shishaldin/2023/'
# MSEED_DIR = f'{TOP_DIR}mseed/'
# chron_path = f'{TOP_DIR}Shishaldin_2023_chronology_review.xlsx'

SVDIR = '/Users/kamerritt/Desktop/Shishaldin/seis_inf_plume/'
#
TOP_DIR = '/Users/kamerritt/Documents/'
MSEED_DIR = '/Users/kamerritt/PycharmProjects/shishaldin_research/'
chron_path_inf = f'{MSEED_DIR}Shishaldin_2023_chronology_submitted.xlsx'
chron_path_vrp = f'{MSEED_DIR}Shishaldin_2023_chronology_submitted.xlsx'
#chron_path_vrp = f'VRP_AllEvents_Cleaned.xlsx'

LOAD_MSEED = True

FONT_S = 18

SOURCE = 'IRIS'
NET = 'AV'
LOC = ''
CHAN = 'BHZ,BDF'

"""dict_keys = ['Event 1', 'Event 2', 'Event 3', 'Event 4', 'Event 5',
             'Event 6', 'Event 8', 'Event 10',
             'Event 11', 'Event 12', 'Event 13']""" #for vrp_noaa

"""dict_keys = ['Event 1', 'Event 2', 'Event 3', 'Event 4', 'Event 5',
             'Event 6', 'Event 8', 'Event 10',
             'Event 11', 'Event 12', 'Event 13'] #for vrp_goes -- no data for event 2"""
dict_keys = ['Event 3', 'Event 6', 'Event 12']
shish_dict = dict.fromkeys(dict_keys)
# shish_dict = {}

# just doing it for a single event now so simplifying a bit
#evt_num = {1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 13}
evt_num = [3, 6, 12]
# evt_num = {7,9}
SEIS_STA = 'SSLS'
INF_STA = 'SSLS'

# select timescale: utc or dr
timescale = 'utc'

# Define time to read data in before RSAM peak
START_OFF = 10 * 3600  # 30
TIME_END = 10 * 3600  # 15

FILT = [1, 10]
FILT_INF = [0.5, 10]
DR_WIN_LEN = 300

PAD = DR_WIN_LEN * 4  # padding length in seconds for edges of data
SAVE = True
SAVE_ALL = False  # scatter plot of all events

VLATLON = (54.7554, -163.9711)
colorm = cc.cm.rainbow  # colormap for plotting
LW = 0.9
# MS = 40 # scatter plot marker size
MS = 55
olap = 0.5
SPECWIN = 50
spec_yl = [.1, 10]
inf_co = 'gray'

timescale = 'utc'

seis_co = 'black'

# plot timeseries of plumes, DR and RMS pressure, and spetrograms
PLOT_TIMESERIES = True

# plot scatter plot of ash height vs. RMS pressure, either for single event
# or all events in a single plot
PLOT_SCATTER = 'SINGLE'  # 'ALL' or 'SINGLE'
DR_AXIS_LIM = [0, 22]
RMS_AXIS_LIM = [0, 1.5]

# rcParams.update({'font.size': FONT_S})
rcParams.update({'font.size': 14})

# %% read in chronology

fig, axes = plt.subplots(nrows=4, ncols=len(evt_num), figsize=(18, 12), sharex='col')


# if only one event, axes won’t be 2D -> ensure consistency
if len(evt_num) == 1:
    axes = np.array(axes).reshape(4, 1)

for col_idx, evt_tmp in enumerate(sorted(evt_num)):
    evtsum = pd.read_excel(chron_path_inf, sheet_name="EventSummary")

    #choose which vrp data to read in
    vrp_h = pd.read_excel(chron_path_inf, sheet_name="VRP_HotLink")
    vrp_g = pd.read_excel(chron_path_inf, sheet_name="VRP_GOES_MIROVA")
    #vrp = pd.read_excel(chron_path_vrp, sheet_name="Sheet1")

    # loop over the events in evt_num
    ct = 1
     # start at line 4 in chron because of 2020
    evt_tmp_ind = evt_tmp + 2  # 14
     # start at line 178 in ash_GOES

    inf_start = evtsum.loc[evt_tmp_ind, 'infrasound_start']
    e_start_dnum = dates.date2num(inf_start)
    inf_end = evtsum.loc[evt_tmp_ind, 'infrasound_end']
    e_end_dnum = dates.date2num(inf_end)

    # evt_num = evtsum['Event']
    rsam_start = evtsum['RSAM_start']
    rsam_peak_time = evtsum['RSAM_peak_time']
    rsam_start_dnum = dates.date2num(evtsum['RSAM_start'])
    rsam_peak_time_dnum = dates.date2num(evtsum['RSAM_peak_time'])

    vrp_g['DateTime'] = pd.to_datetime(vrp_g['DateTime'], errors='coerce')
    vrp_g_start = vrp_g['DateTime']
    vrp_g_start_dnum = dates.date2num(vrp_g_start)
    vrp_g_end_shifted = vrp_g_start.shift(-1).dropna()
    vrp_g_end_dnum = dates.date2num(vrp_g_end_shifted)
    vrp_h['DateTime'] = pd.to_datetime(vrp_h['DateTime'], errors='coerce')
    vrp_h_start = vrp_h['DateTime']
    vrp_h_start_dnum = dates.date2num(vrp_h_start)
    vrp_h_end_shifted = vrp_h_start.shift(-1).dropna()
    vrp_h_end_dnum = dates.date2num(vrp_h_end_shifted)
    ash = pd.read_excel(chron_path_inf, sheet_name="Ash_GOES")
    ash_start = ash['DateTime']
    ash_start_dnum = dates.date2num(ash_start)
    ash_end_dnum = dates.date2num(ash_start.shift(-1))

    if timescale == 'dr':
        inf_start_hr = (e_start_dnum - rsam_peak_time_dnum[evt_tmp_ind]) * 24
        inf_end_hr = (e_end_dnum - rsam_peak_time_dnum[evt_tmp_ind]) * 24

        ash_start_off = (ash_start_dnum - rsam_peak_time_dnum[evt_tmp_ind]) * 24
        ash_end_off = (ash_end_dnum - rsam_peak_time_dnum[evt_tmp_ind]) * 24

    elif timescale == 'utc':
        if evt_tmp == 7:
                inf_start_hr = UTCDateTime('2023-8-14T08:00:00')
                inf_end_hr = UTCDateTime('2023-8-16T10:44:00')
        elif evt_tmp == 9:
            inf_start_hr = UTCDateTime('2023-9-04T23:00:00')
            inf_end_hr = UTCDateTime('2023-9-5T19:03:00')

        else:
            inf_start_hr = UTCDateTime(e_start_dnum)
            inf_end_hr = UTCDateTime(e_end_dnum)

            ash_start_off = ash_start_dnum
            ash_end_off = ash_end_dnum

        # # Timeseries plot individual - Ash Plumes

        # rsam_SSBA = pd.read_excel(chron_path, sheet_name="RSAM_SSBA")
        # DR = rsam_SSBA['DR']
        # DR_tvec = dates.date2num(rsam_SSBA['DateTime'])

    max_t_off = rsam_peak_time - rsam_start

        # %% read in data
    client = Client(SOURCE)
    st_raw_bak = Stream()

    print(f'Reading data from mseed file for Event {evt_tmp}')
    st_inf = read(f'{MSEED_DIR}event{evt_tmp}_{INF_STA}-BDF.mseed', format="MSEED")
    st_seis = read(f'{MSEED_DIR}event{evt_tmp}_{SEIS_STA}-BHZ.mseed', format="MSEED")
    inv = client.get_stations(network=NET, station='SSLS,SSLN,SSBA',
                                  channel='BDF,BHZ', level="response")
    # st_raw_bak.attach_response(inv)

    """for evt_num_tmp in evt_num:
            i=evt_num_tmp-1
            print(f'Running for event {evt_num_tmp}')
    
            st_seis = gather_waveforms('IRIS', NET, SEIS_STA, LOC, 'BHZ', UTCDateTime(rsam_peak_time[i])-START_OFF-PAD,
                                        UTCDateTime(rsam_peak_time[i])+TIME_END+PAD, remove_response=False)
            st_inf = gather_waveforms('IRIS', NET, INF_STA, LOC, 'BDF', UTCDateTime(rsam_peak_time[i])-START_OFF-PAD,
                                        UTCDateTime(rsam_peak_time[i])+TIME_END+PAD, remove_response=False)"""


        # convert seismic to displacement
    st_seis.integrate()

    # combine streams
    st_raw = st_seis + st_inf

    """st_raw.trim(starttime=UTCDateTime(rsam_peak_time[evt_tmp_ind]) - START_OFF - PAD,
                    endtime=UTCDateTime(rsam_peak_time[evt_tmp_ind]) + TIME_END + PAD)"""
        # st_raw_bak = st_seis+st_inf

        # plot raw waveforms as a check
        # fig1=plt.figure()
        # st_raw.plot(fig=fig1)

    for tr in st_raw:
        coords = inv.get_coordinates(tr.id)
        tr.stats.longitude = coords['longitude']
        tr.stats.latitude = coords['latitude']
        tr.stats.elevation = coords['elevation']
        # st_raw = Stream(st_raw_bak[i])
    print(f'Computing metrics for {st_raw[0].stats.starttime} - {st_raw[0].stats.endtime}')

    tmpl_s, Dr, pe, fc, fd, fsd = compute_metrics(Stream(st_raw[0]), process_taper=True, metric_taper=None,
                                                      filter_band=FILT,
                                                      window_length=DR_WIN_LEN, overlap=0.75, vlatlon=VLATLON)
        # print(Dr)
    print('Done with seis')

    tmpl_i, RMS_p, pe, fc, fd, fsd = compute_metrics(Stream(st_raw[1]), process_taper=True, metric_taper=None,
                                                         filter_band=FILT_INF,
                                                         window_length=DR_WIN_LEN, overlap=0.75,
                                                         vlatlon=VLATLON)
    print('Done with infra')

    # %%
    tr = Trace()
    tr.stats.network = st_raw[0].stats.network
    tr.stats.station = st_raw[0].stats.station
    tr.stats.location = st_raw[0].stats.location
    tr.stats.channel = st_raw[0].stats.channel
    tr.data = Dr[0, :]
    tr.stats.starttime = UTCDateTime(dates.num2date(tmpl_s[0, 1]))
    tr.stats.sampling_rate = 1 / np.round(86400 * (tmpl_s[0, 1] - tmpl_s[0, 0]))

    # correct DR and RMSP or different stations (updated scaling values 7/17)
    if tr.stats.station == 'SSLS' and tr.stats.channel == 'BHZ':
        tr.data = tr.data / 2.02
    elif tr.stats.station == 'SSLN' and tr.stats.channel == 'BHZ':
        tr.data = tr.data / 3.63

    if tr.stats.station == 'SSLS' and tr.stats.channel == 'BDF':
        tr.data = tr.data / 2.76
    elif tr.stats.station == 'SSLN' and tr.stats.channel == 'BDF':
        tr.data = tr.data / 4.58

    st_Dr = Stream()
    st_Dr.append(tr)

    # build infrasound
    tr = Trace()
    tr.stats.network = st_raw[1].stats.network
    tr.stats.station = st_raw[1].stats.station
    tr.stats.location = st_raw[1].stats.location
    tr.stats.channel = st_raw[1].stats.channel
    if evt_tmp == 7 or evt_tmp == 9:
        print(f'Event {evt_tmp}: filling with nans')
        tr.data = np.nan * np.zeros(len(RMS_p[0, :]))
    else:
        tr.data = RMS_p[0, :]
    tr.stats.starttime = UTCDateTime(dates.num2date(tmpl_i[0, 1]))
    tr.stats.sampling_rate = 1 / np.round(86400 * (tmpl_i[0, 1] - tmpl_i[0, 0]))

    if timescale == 'dr':
        # time vector in hours
        tvec = st_Dr[0].times()  # /3600
        tvec_off = (tvec - START_OFF - PAD) / 3600
    elif timescale == 'utc':
        tvec_off = st_Dr[0].times("matplotlib")

    st_Dr.append(tr)
    st_Dr.taper(.01)

    # st_Dr_10Hz = st_Dr.copy()
    if timescale == 'dr':
        # plume times in hours
        vrp_g_start_off = (vrp_g_start_dnum - rsam_peak_time_dnum[
            evt_tmp_ind]) * 24
        vrp_g_end_off = (vrp_g_end_dnum - rsam_peak_time_dnum[evt_tmp_ind]) * 24

        vrp_h_start_off = (vrp_h_start_dnum - rsam_peak_time_dnum[
            evt_tmp_ind]) * 24
        vrp_h_end_off = (vrp_h_end_dnum - rsam_peak_time_dnum[evt_tmp_ind]) * 24

    elif timescale == 'utc':
        vrp_g_start_off = vrp_g_start_dnum
        vrp_g_end_off = vrp_g_end_dnum

        vrp_h_start_off = vrp_h_start_dnum
        vrp_h_end_off = vrp_h_end_dnum

        rms_times = [dates.num2date(t) for t in tmpl_i[0, :]]
        rms_times = pd.to_datetime(rms_times, utc=True)
        rms_vals = RMS_p[0, :]

        dr_times = [dates.num2date(t) for t in tmpl_s[0, :]]
        dr_times = pd.to_datetime(dr_times, utc=True)
        dr_vals = Dr[0, :]

        vrp_g['DateTime'] = pd.to_datetime(vrp_g['DateTime'], errors='coerce', utc=True)
        vrp_g = vrp_g.dropna(subset=['DateTime'])

        vrp_h['DateTime'] = pd.to_datetime(vrp_h['DateTime'], errors='coerce', utc=True)
        vrp_h = vrp_h.dropna(subset=['DateTime'])

        t_start = pd.to_datetime(rms_times[0])
        t_end = pd.to_datetime(rms_times[-1])

            # Filter VRP within time window
            #print("RMS window:", t_start, "to", t_end)
            #print("VRP times min:", vrp['DateTime'].min(), "max:", vrp['DateTime'].max())

        vrp_g_window = vrp_g[(vrp_g['DateTime'] >= t_start) & (vrp_g['DateTime'] <= t_end)]

        vrp_g_times = vrp_g_window['DateTime'].tolist()
        vrp_g_vals = vrp_g_window['TotalRP'].values

        vrp_h_window = vrp_h[(vrp_h['DateTime'] >= t_start) & (vrp_h['DateTime'] <= t_end)]

        vrp_h_times = vrp_h_window['DateTime'].tolist()
        vrp_h_vals = vrp_h_window['VRP'].values

            #print(vrp_window[['DateTime', 'TotalRP']])
            #print(vrp_h_window['VRP'])

            # %% trigger counts

        sta_lta_data = np.load(f'{MSEED_DIR}dict_shish_INF_STA_LTA_THRESH2.5_0903.npy', allow_pickle=True).item()

        evt_key = f'Event {evt_tmp}'
        if evt_key in sta_lta_data and sta_lta_data[evt_key] is not None:
            counts = np.array(sta_lta_data[evt_key]['trigger_counts'])
            bins = np.array(sta_lta_data[evt_key]['bins'])

            if len(bins) == len(counts) + 1:
                bins = (bins[:-1] + (bins[1:] - bins[:-1]) / 2)
            sta_lta_times = [b.datetime.replace(tzinfo=pytz.UTC) for b in bins]

            if len(sta_lta_times) > 1:
                bin_width = np.median(np.diff([t.timestamp() for t in sta_lta_times])) / 86400.0

            else:
                bin_width = 60.0 / 86400.0

        else:
            counts = []
            sta_lta_times = []
            bin_width = 0.001

        # %% plotting

    ax1 = axes[0, col_idx]
    ax2 = axes[1, col_idx]
    ax3 = axes[2, col_idx]
    ax4 = axes[3, col_idx]

    # RMS subplot
    ax1.plot(rms_times, rms_vals, color='gray', lw=2)
    ax1.set_ylabel("RMSP [Pa]", color='gray')
    ax1.set_ylim(0, 10)

    ax1b = ax1.twinx()
    ax1b.plot(tvec_off, st_Dr[0].data, color='black', lw=2)
    if evt_tmp in [2, 6, 10, 13]:
        ax1b.set_ylim(0, 25)
    ax1b.set_ylabel("DR [cm$^2$]", color='black')
    if evt_tmp == 6:
        ax1.set_title(f"Seismoacoustics, VRP, Explosions, & Plumes for Events 3, 6, and 12")
    ax1.axvline(pd.to_datetime(inf_start, utc=True), color='gray', ls='--', lw=2, label='Infrasound Start/End')
    ax1.axvline(pd.to_datetime(inf_end, utc=True), color='gray', ls='--', lw=2)
    ax1.grid(True)
    ax1.legend(loc='upper right')

    # VRP subplot
    ax2.plot([], [])
    ax2.set_xlabel("Time (UTC)")
    ax2.set_ylabel("GOES VRP [W]", color='firebrick')
    ax2.scatter(vrp_g_times, vrp_g_vals, color='firebrick', marker='o', s=50)
    #ax2.set_ylim(1e5, 2e10)
    ax2.set_yscale('log')
    ax2.set_ylim(1e6, 2e10)

    """ax2b = ax2.twinx()
    ax2b.scatter(vrp_h_times, vrp_h_vals, color='orange', marker='*', s=50)
    ax2b.set_ylabel("HotLINK VRP [W]", color='orange')
    #ax2b.set_ylim(1e3, 1e9)
    ax2b.set_yscale('log')"""


    ax2.axvline(pd.to_datetime(inf_start, utc=True), color='gray', ls='--', lw=2)
    ax2.axvline(pd.to_datetime(inf_end, utc=True), color='gray', ls='--', lw=2)

    ax2.xaxis.set_major_formatter(DateFormatter("%d %H:%M"))
    ax2.tick_params(axis='x', rotation=45)
    ax2.grid(True)

    # --- STA/LTA trigger counts subplot ---
    ax3.bar(sta_lta_times, counts, width=bin_width, align='center', edgecolor='black')
    ax3.set_ylabel("Trigger Counts")
    ax3.set_xlabel("Time (UTC)")
    ax3.grid(True)
    ax3.set_ylim(0, 25)

    ax3.axvline(pd.to_datetime(inf_start, utc=True), color='gray', ls='--', lw=2)
    ax3.axvline(pd.to_datetime(inf_end, utc=True), color='gray', ls='--', lw=2)

    ax3.xaxis.set_major_formatter(DateFormatter("%d %H:%M"))
    ax3.tick_params(axis='x', rotation=45)

    ash_times = pd.to_datetime(ash['DateTime'], errors='coerce', utc=True)

    ash_starts = ash_times[:-1].reset_index(drop=True)
    ash_ends = ash_times.shift(-1)[:-1].reset_index(drop=True)
    heights_raw = ash['Height_km'][:-1].reset_index(drop=True)

    ash_starts = pd.to_datetime(ash_starts, errors='coerce', utc=True)
    ash_ends = pd.to_datetime(ash_ends, errors='coerce', utc=True)

    mask = (~ash_starts.isna()) & (~ash_ends.isna())
    ash_starts = ash_starts[mask].reset_index(drop=True)
    ash_ends = ash_ends[mask].reset_index(drop=True)
    heights_raw = heights_raw[mask].reset_index(drop=True)

    ash_start_py = [t.to_pydatetime() for t in ash_starts if pd.notna(t)]
    ash_end_py = [t.to_pydatetime() for t in ash_ends if pd.notna(t)]

    if not isinstance(ash_start_py, list):
        ash_start_py = [ash_start_py]
    if not isinstance(ash_end_py, list):
        ash_end_py = [ash_end_py]

    ash_start_md = mdates.date2num(ash_start_py)
    ash_end_md = mdates.date2num(ash_end_py)
    widths_md = ash_end_md-ash_start_md

    if len(ash_start_md) > 0:
        for j in range(len(ash_start_md)):
            if pd.isna(heights_raw[j]):
                h = 15
                color = 'gray'
            else:
                h = heights_raw[j]
                color = 'darkred'
            rect = Rectangle((ash_start_md[j], 0), widths_md[j], h, color=color, alpha=0.5)
            ax4.add_patch(rect)

    ax4.set_ylabel("Height [km]")
    ax4.set_xlabel("Time (UTC)")
    ax4.set_ylim(0, 15)
    ax4.grid(True)
    ax4.axvline(pd.to_datetime(inf_start, utc=True), color='gray', ls='--', lw=2)
    ax4.axvline(pd.to_datetime(inf_end, utc=True), color='gray', ls='--', lw=2)
    ax4.xaxis.set_major_formatter(DateFormatter("%d %H:%M"))
    ax4.tick_params(axis='x', rotation=45)

for r in range(4):
    for c in range(len(evt_num)):
        axes[r, c].xaxis.set_major_formatter(DateFormatter("%d %H:%M"))
        axes[r, c].tick_params(axis='x', rotation=45)

# Y labels only on leftmost column
axes[0,0].set_ylabel("RMSP [Pa]", color='gray')
axes[1,0].set_ylabel("GOES VRP [W]", color='firebrick')
axes[2,0].set_ylabel("Trigger Counts")
axes[3,0].set_ylabel("Height [km]")


# X labels only on bottom row
for c in range(len(evt_num)):
    axes[3, c].set_xlabel("Time (UTC)")

plt.tight_layout()
plt.show()