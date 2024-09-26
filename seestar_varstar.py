import socket
import json
import time
from datetime import datetime, timezone
import threading
import sys
import argparse
import pandas as pd
import logging
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import ephem
import pytz
import json
import seestar_varstar_params as params

def CreateLogger():
    # Create a custom logger 
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Create handlers
    console_handler = logging.StreamHandler(sys.stdout)
    file_handler = logging.FileHandler('seestar_varstar.log')

    # Set levels for handlers
    console_handler.setLevel(logging.INFO)
    file_handler.setLevel(logging.DEBUG)

    # Create formatters and add them to handlers
    console_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    console_handler.setFormatter(console_format)
    file_handler.setFormatter(file_format)

    # Add handlers to the logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    return(logger)

def heartbeat(): #I noticed a lot of pairs of test_connection followed by a get if nothing was going on
    json_message("test_connection")
#    json_message("scope_get_equ_coord")

def send_message(data):
    global s
    try:
        s.sendall(data.encode())  # TODO: would utf-8 or unicode_escaped help here
    except socket.error as e:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((HOST, PORT))
        send_message(data)

def get_socket_msg():
    global s
    try:
        data = s.recv(1024 * 60)  # comet data is >50kb
    except socket.error as e:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((HOST, PORT))
        data = s.recv(1024 * 60)
    data = data.decode("utf-8")
    if is_debug:
        logger.debug(f"Socket received : {data}")
    return data
    
def receieve_message_thread_fn():
    global is_watch_events
    global op_state
    global s
        
    msg_remainder = ""
    while is_watch_events:
        #print("checking for msg")
        data = get_socket_msg()
        if data:
            msg_remainder += data
            first_index = msg_remainder.find("\r\n")
            
            while first_index >= 0:
                first_msg = msg_remainder[0:first_index]
                msg_remainder = msg_remainder[first_index+2:]            
                parsed_data = json.loads(first_msg)
                
                if 'Event' in parsed_data and parsed_data['Event'] == "AutoGoto":
                    state = parsed_data['state']
                    logger.debug("AutoGoto state: %s" % state)
                    if state == "complete" or state == "fail":
                        op_state = state
                
                if is_debug:
                    logger.debug(parsed_data)
                    
                first_index = msg_remainder.find("\r\n")
        else:
            logger.debug(f'no data from socket')
            #sys.exit(1)
        time.sleep(1)

def json_message(instruction):
    global cmdid
    data = {"id": cmdid, "method": instruction}
    cmdid += 1
    json_data = json.dumps(data)
    if is_debug:
        logger.debug("Sending %s" % json_data)
    send_message(json_data+"\r\n")

def json_message2(data):
    if data:
        json_data = json.dumps(data)
        if is_debug:
            logger.debug("Sending2 %s" % json_data)
        resp = send_message(json_data + "\r\n")


def goto_target(ra, dec, target_name, exp_time, exp_cont, is_lp_filter):
    global cmdid
    """
    First we set the exposure parameters for this target
    """
    logger.info(f'Setting parameters for {target_name}')
    data = {}
    data['id'] = cmdid
    cmdid += 1
    data['method'] = 'set_setting'
    params = {}
    params['exp_ms'] = {}
    params['exp_ms']['stack_l']=int(exp_time)*1000
    params['exp_ms']['continous']=int(exp_cont)*1000
    data['params'] = params
    logger.debug(f'Exposure Settings: {data}')
    json_message2(data)
    """
    Then we slew to this target
    """
    logger.info(f"going to target {target_name}...")
    data = {}
    data['id'] = cmdid
    cmdid += 1
    data['method'] = 'iscope_start_view'
    params = {}
    params['mode'] = 'star'
    ra_dec = [ra, dec]
    params['target_ra_dec'] = ra_dec
    params['target_name'] = target_name
    params['lp_filter'] = int(is_lp_filter)
    data['params'] = params
    logger.debug(f'Target Settings: {data}')
    json_message2(data)
    
def start_stack():
    global cmdid
    logger.info("starting to stack...")
    data = {}
    data['id'] = cmdid
    cmdid += 1
    data['method'] = 'iscope_start_stack'
    params = {}
    params['restart'] = True
    data['params'] = params
    json_message2(data)

def stop_stack():
    global cmdid
    logger.info("stop stacking...")
    data = {}
    data['id'] = cmdid
    cmdid += 1
    data['method'] = 'iscope_stop_view'
    params = {}
    params['stage'] = 'Stack'
    data['params'] = params
    json_message2(data)

def wait_end_op():
    global op_state
    op_state = "working"
    heartbeat_timer = 0
    while op_state == "working":
        heartbeat_timer += 1
        if heartbeat_timer > 5:
            heartbeat_timer = 0
            json_message("test_connection")
        time.sleep(1)

    
def sleep_with_heartbeat(session_time):
    stacking_timer = 0
    while stacking_timer < session_time:         # stacking time for each object
        stacking_timer += 1
        if stacking_timer % 5 == 0:
            json_message("test_connection")
        time.sleep(1)

def parse_ra_to_float(ra_string):
    # Split the RA string into hours, minutes, and seconds
    ra_string = ra_string[1:-1]
    hours, minutes, seconds = map(float, ra_string.split(':'))

    # Convert to decimal degrees
    ra_decimal = hours + minutes / 60 + seconds / 3600

    return ra_decimal
    
def parse_dec_to_float(dec_string):
    # Split the Dec string into degrees, minutes, and seconds
    dec_string = dec_string[1:-1]
    if dec_string[0] == '-':
        sign = -1
        dec_string = dec_string[1:]
    else:
        sign = 1
    degrees, minutes, seconds = map(float, dec_string.split(':'))

    # Convert to decimal degrees
    dec_decimal = sign * (degrees + minutes / 60 + seconds / 3600)

    return dec_decimal

def set_stack_settings():
    global cmdid
    logger.debug("set stack setting to record individual frames")
    data = {}
    data['id'] = cmdid
    cmdid += 1
    data['method'] = 'set_stack_setting'
    params = {}
    params['save_discrete_frame'] = True
    data['params'] = params
    return(json_message2(data))

def shutdown_seestar():
    global cmdid
    data = {}
    data['id'] = cmdid
    cmdid+=1
    data['method'] = 'pi_shutdown'
    json_message2(data)

def get_coord_object(target_names):
    result_table = Simbad.query_objects(target_names)
    object_ra = result_table['RA'].data  # Right Ascension
    object_dec = result_table['DEC'].data  # Declination
    coord=SkyCoord(object_ra, object_dec,unit=(u.deg))
    return coord.ra.deg, coord.dec.deg

def calc_twilight():
    S50Observer = ephem.Observer()
    # Set the date and time 
    S50Observer.date = datetime.today().strftime('%Y-%m-%d')#"2024-07-14"

    # Location 
    S50Observer.lat = str(params.Longitude)
    S50Observer.lon = str(params.Latitude)

    # Elevation 
    S50Observer.elev = params.Elevation

    # To get U.S. Naval Astronomical Almanac values, use these settings
    S50Observer.pressure = 0
    S50Observer.horizon = '-0:34'

    # Calculate sunrise, solar noon, and sunset
    S50Observer.horizon = '-12'  # -6=civil twilight, -12=nautical, -18=astronomical
    beg_twilight = S50Observer.next_rising(ephem.Sun(), use_center=True)
    end_twilight = S50Observer.next_setting(ephem.Sun(), use_center=True)
    return(beg_twilight.datetime().replace(tzinfo=pytz.utc), end_twilight.datetime().replace(tzinfo=pytz.utc))

    
is_watch_events = True
    
def main():
    global HOST
    global PORT
    global session_time
    global s
    global cmdid
    global is_watch_events
    global is_debug, logger
    
    logger = CreateLogger()
    version_string = params.S50_run_version
    logger.info(f"Seestar Control for Variable Star Observing Version: {version_string}")
    
    parser = setup_argparse()
    args = parser.parse_args()
    HOST = params.ip
    is_debug = args.is_debug
    is_use_LP_filter = 0 # do not use this filter for photometric observations
    PORT = params.port
    cmdid = 999

    # Get the schedule of targets
    try:
        target_df = pd.read_csv(args.schedule_file)
    except Exception as e:
         logger.error(f'Unable to load schedule - {e}')
         sys.exit(1)
    target_names=target_df['Name'].values
    target_stack_times = target_df['TotalExp'].values
    target_exptimes = target_df['ExpTime'].values
    
    # Get object coordinate from simbad query and convert to Jnow
    ras,decs = get_coord_object(target_names)
    logger.info(f'list of targets {target_names, ras, decs}')

    # TARGET OBSERVATIONS
    # determine safe limits for observing
    (twilight_begin, twilght_end) = calc_twilight()

    logger.info(f'Nautical Twilight is from {twilight_begin} to {twilght_end}')

    S50safe = False
    while not S50safe:
        if (twilight_begin <= datetime.now(timezone.utc) <= twilght_end):
            S50safe = True
        else:
            logger.info('Outside observing times - waiting...')
            logger.info(f'Current UTC: {datetime.now(timezone.utc)}')
            time.sleep(1)
            break #testing

    # determine the repetition pattern requested
    if (args.target_seq_mode not in ['repeat', 'single']):
        logger.error(f'target sequence mode not known: {args.target_seq_mode}')
        raise RuntimeError()
    elif (args.target_seq_mode == 'repeat'):
            logger.info(f'Targets ({len(ras)}) will be cycled repeatedly until dawn - mode {args.target_seq_mode}')
            repeat = True
    elif (args.target_seq_mode == 'single'):
            logger.info(f'Targets ({len(ras)}) will be observed in order - mode {args.target_seq_mode}')
            repeat = False        

    # preliminary overall settings
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        s.connect((HOST, PORT))
        logger.debug('Seestar connection is good')
    except Exception as e:
            logger.error(f'Unable to connect to Seestar - check powered on')
            sys.exit(1)
    set_stack_settings() # set to save frames in the stack
    ras = [10.5, 12.5]
    decs = [-40, -60]
    loop = True
    with s:
        # flush the socket input stream for garbage
        get_socket_msg()
        iloop = 0
        nloops = 0
        while loop:
            # print input requests
            logger.info("Received parameters:")
            logger.debug(f"  ip address    :  {HOST}")
            logger.info(f"  target        :  {target_names[iloop]}")
            logger.debug(f"  RA            :  {ras[iloop]}")
            logger.debug(f"  Dec           :  {decs[iloop]}")
            logger.info(f"  exp time      :  {target_exptimes[iloop]}")
            logger.info(f"  session time  :  {target_stack_times[iloop]}")
            logger.debug(f"  debug         :  {is_debug}")

        
            get_msg_thread = threading.Thread(target=receieve_message_thread_fn)
            get_msg_thread.start()
        
            goto_target(ras[iloop], decs[iloop], target_names[iloop], target_exptimes[iloop], target_stack_times[iloop], is_lp_filter=0)
            wait_end_op()
            logger.info("Goto operation finished")
            time.sleep(3)
            if op_state == "complete":
                start_stack()    
                sleep_with_heartbeat(target_stack_times[iloop])
                stop_stack()
                logger.info("Stacking operation finished" + target_names[iloop])
            else:
                logger.warning("Goto failed.")
            # test if we are safe to continue with observing
            if not (twilight_begin <= datetime.now(timezone.utc) <= twilght_end):
                logger.warning('Now outside observing time... shutting down')
                logger.info(f'Current UTC: {datetime.now(timezone.utc)}')
                if not is_debug:
                    shutdown_seestar()
                    break
            iloop +=1
            if (iloop>=len(target_names) and repeat):
                    iloop = 0
                    nloops +=1
                    logger.info(f'Loop {nloops} executed')
            elif (iloop>=len(target_names) and not repeat):
                    loop = False    
    logger.info("Finished seestar_run")
    is_watch_events = False
    get_msg_thread.join(timeout=10)
    time.sleep(10)
    s.close()
#    if not is_debug:
    shutdown_seestar()

    
def setup_argparse():
    parser = argparse.ArgumentParser(description='Seestar Run_VarStar')
    parser.add_argument('schedule_file', type=str, help="Observation Target Schedule")
    parser.add_argument('target_seq_mode', type=str, choices=['repeat', 'single'], default='single', help='Schedule Mode: repeat or single loop through targets')
    parser.add_argument('is_debug', type=str, default=False, nargs='?', help="optional - print debug logs while running.")

    return parser
    
if __name__ == "__main__":
    main()
    

 
 #conda activate wombat
 #cd C:\Users\Admin\Documents\GitHub\seestar_run-main\seestar_run-main
 #python seestar_varstar.py schedule_20240925.dat single False
 #"""