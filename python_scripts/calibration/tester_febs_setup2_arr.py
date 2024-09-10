import time
import sys
import os
import logging
import pprint
import uhal
import subprocess
from datetime import datetime
sys.path.append('../autogen/agwb/python/')
import agwb
#import test_setup as setup
from smx_tester import *
import msts_defs as smc
#from opm import Opm
import smx_oper as smxoper

# def rawfifo_to_smxframes(rfifo_dat):
#     smxframes = {}
#     for u in rfifo_dat.keys():
#         smxframes[u] = []
#         for d in rfifo_dat[u]:
#             smxframes[u].append(format.smx.determinator(d.smxframe))
#             # systime is not a part of smxframe, but is is convenient to add it to python class
#             smxframes[u][-1].systime = d.systime
#     return smxframes
            

# --------------------------- SETTING LOGGING DETAILS ------------------------------------------


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s,%(msecs)03d:%(module)s:%(levelname)s:%(message)s",
    datefmt="%H:%M:%S",
)

log = logging.getLogger()
def set_logging_details(module_dir = "module_files"):    
    date  = time.strftime("%y%m%d_%H%M")
    logFileName = sys.argv[0].replace('.py', '')+"_"+ date +".log" 
    logFilePath = module_dir + "/" + logFileName 
    fh = logging.FileHandler(logFilePath, 'w')
    fh.setLevel(logging.DEBUG)
    fmt = logging.Formatter('[%(levelname)s] %(message)s')
    fh.setFormatter(fmt)
    log.addHandler(fh)
    uhal.setLogLevelTo(uhal.LogLevel.WARNING)

# -----------------------------------------------------------
# --- User tests start here 
# -----------------------------------------------------------

# define list of ASICs used for tests
#  this is a sublist of the detected Asics 'smxes' 
#smx_l = smxes[0:3]   # example: first three(!) Asics 
#smx_l = smxes   #[0:1]      # all Asics

#opm = Opm(smx_l)
 
# configure Asics with default configuration and readback config registers 
#  - configuration from 'default' dictionary in file smx_tester/mcbm_msts_active_params.py



# all below are just examples andare optional
# --- Do some tests --------
'''
log.info("Reading SMX chip address")
for smx in smx_l:
    addr = smx.read(192, 22)
    log.info("Asic address (emu %d   downlink %d   hw_addr %d): \t %d", smx.rob, smx.downlink, smx.address, addr)

for smx in smx_l:
	smx.reset_asic()

for smx in smx_l:
    smx.set_channel_mask_enable_all()

#  set selected parameters  
#  - example: set the CSA currents
#      - list of available Asic config setting names: see the above file
for smx in smx_l:
    smx.conf_func("csa_front",31)
    smx.conf_func("csa_back", 31)
    log.info(" CSA settings: front %d   back %d", smx.read(130,0)&0xff, smx.read(130,13)&0xff)

for smx in smx_l:
    v_diag, dac_diag = smx.read_diag("Temp")
    log.info("Asic %d: Temp %f mV (ADC %d)", smx.address, v_diag, dac_diag)
    log.info("Asic %d: Temp calib. %f ", smx.address, smx.read_temp())
    v_diag, dac_diag = smx.read_diag("Vddm")
    log.info("Asic %d: Vddm %f mV (ADC %d)", smx.address, v_diag, dac_diag)
log.info("")
'''

#  Get some specific smx by (downlink, hw_address) or (downlink, aseq)
#   "aseq" is "asic_number"
#smx = opm.smx_dl_aseq(1,0) #original
# try:
#     #smx = opm.smx_dl_hw(1,4) #my
#     #smx = smx_l[0]
#     addr = smx.read(192, 22)
#     log.info("Asic address (emu %d   downlink %d   aseq %d): \t %d", smx.rob, smx.downlink, smx.aseq, addr)
# except:
#     pass



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ------------------------ CONNECTING EMUs, FINDING FEBs AND SYNC ASICs -----------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def general_sync(emu_list, active_downlinks = [1,2]):
    # Function to establish connection with EMU boards, finding FEBs and synchronizing ASICs 
    manager = uhal.ConnectionManager("file://devices.xml")
    # Open xml file with list of EMU devices
    
    # Arrays to organize setup elements (FEBs)
    setup_elements = []
    emu_elements = []

    link_scan = True
    device_scan = True
    link_sync = True

    for emu in emu_list:
        ipbus_interface = IPbusInterface(manager, emu)
        agwb_top = agwb.top(ipbus_interface, 0)
        smx_tester = SmxTester(agwb_top, CLK_160)
        setup_elements_tmp = smx_tester.scan_setup()
        log.info("setup_elements_tmp: %d", len(setup_elements_tmp))
        setup_elements.extend(setup_elements_tmp)
        emu_elements.extend( [emu for i in range(len(setup_elements_tmp))])

    #setup_elements = [setup_elements[0]]
    #active_downlinks = active_downlinks
    #active_downlinks = [0,3]
    #active_downlinks = [1,2]
    #active_downlinks = [3]

    remove_list_01 = []
    remove_list_23 = []
    
    for se in list(setup_elements):
        log.info(f"SE:\t{se.downlink}\t\t{se.uplinks}  ({len(se.uplinks)})")
        if not se.downlink in active_downlinks:
            setup_elements.remove(se)
            log.info(f"Remove setup element for downlink no. {se.downlink} with uplinks {se.uplinks}")

    for se in setup_elements:
        log.info(f"Cleaned - SE:\t{se.downlink}\t\t{se.uplinks}  ({len(se.uplinks)})")

        for uplink in se.uplinks:
            if se.downlink == 0 or se.downlink == 1:
                if uplink > 15:
                    remove_list_01.append(uplink)
                    log.info("Uplink appended to remove list: {}".format(uplink))
            else:
                if uplink < 16:
                    remove_list_23.append(uplink)
                    log.info("Uplink appended to remove list: {}".format(uplink))
        if se.downlink == 0 or se.downlink == 1:
            se.uplinks = [i for i in se.uplinks if i not in remove_list_01]
        else:
            se.uplinks = [i for i in se.uplinks if i not in remove_list_23]
        log.info(f"Final List UPLINKS SE:\t{se.downlink}\t\t{se.uplinks}  ({len(se.uplinks)})")
        
    # Temporary workaround becase there are no independent elink clocks.
    # setup_elements = [setup_elements[0]]
    # smx_tester.elinks.disable_downlink_clock(4)

    if link_scan:
        for se in setup_elements:
            se.characterize_clock_phase()
            log.info(f"post char clk phase\n {se}")
            se.initialize_clock_phase()
            log.info(f"post set\n {se}")

        for se in setup_elements:
            se.characterize_data_phases()
            log.info(f"post char data phase\n {se}")
            se.initialize_data_phases()
            log.info(f"post set\n {se}")

    if device_scan:        
        for se in setup_elements:
            se.scan_smx_asics_map()
            log.info(f"post scan map|n {se}")

    if link_sync:        
        for se in setup_elements:
            se.synchronize_elink()

        for se in setup_elements:
            se.write_smx_elink_masks()

    smxes = []
    for se, emu in zip(setup_elements, emu_elements):
        smxes_tmp = smxes_from_setup_element(se)
        log.info("smxes_tmp: %d", len(smxes_tmp))
        for smx in smxes_tmp:
            smx.rob = int(emu[4:])
        smxes.extend(smxes_tmp)
    
    print(emu_elements) 
    log.info("Number of setup elements: %d", len(setup_elements))
    log.info("Number of smx: %d", len(smxes))
    log.info("")
    return smxes

def scanning_asics(smx_l):
    # Function to find the number of ASICs per dwnlink or devices    
    n_asics = 0
    p_asics = 0
    for smx in smx_l:
        log.info("Asic (emu %d   downlink %d   hw_addr %d): ", smx.rob, smx.downlink, smx.address)
        smx.func_to_reg(smc.R_ACT)
        smx.write_reg_all()
        smx.read_reg_all(compFlag = False)
        if (smx.downlink == 1):
            n_asics+=1
        elif (smx.downlink == 0):
            n_asics+=1
        elif (smx.downlink == 2):
            p_asics+=1
        else:
            p_asics+=1       
    return (n_asics, p_asics)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ----------------------------- MODULE TESTING: LV & HV CONTROLS ------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def reading_lv(pol_usr = 'N'):
    # Function to read the voltage and current of the FEB
    # LV mapping                                                                                                                     
    # LV channels vs FEB                                                                                                             
    # u904 ................ 1.2 V LDOs N-side                                                                                        
    # u905 ................ 1.8 V LDOs N-side                                                                                        
    # u906 ................ 1.2 V LDOs P-side                                                                                        
    # u907 ................ 1.8 V LDOs P-side                                                                                        

    pol = pol_usr

    if (pol == 'N' or pol == '0'):
        #reading LV in channels u904 and u905                                                                                         
        print("Reading LV for N-side (electrons polarity):")
        #Reading Voltage and current values and writting them in a file                                                               
        #os.system("ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_u700 voltage | grep -Eo '[0-9]+([.][0-9]+)' > lv_iv.txt")
        os.system("./../../../lv/mpod_hhlab_2024/lab24_mpod02_u500 voltage | grep -Eo '[0-9]+([.][0-9]+)' > lv_iv.txt")
        os.system("./../../../lv/mpod_hhlab_2024/lab24_mpod02_u501 voltage | grep -Eo '[0-9]+([.][0-9]+)' >> lv_iv.txt")
    elif(pol == 'P' or pol == '1'):
        print("Reading LV for P-side (holes polarity):")
        #Reading Voltage and current values and writting them in a file                                                               
        os.system("./../../../lv/mpod_hhlab_2024/lab24_mpod02_u502 voltage | grep -Eo '[0-9]+([.][0-9]+)' > lv_iv.txt")
        os.system("./../../../lv/mpod_hhlab_2024/lab24_mpod02_u503 voltage | grep -Eo '[0-9]+([.][0-9]+)' >> lv_iv.txt")
    else:
        log.error("Please, indicate a polarity for reading the corresponding LV potentials, in the following way:")
        log.error("N or 0 for n-side, 0 for electrons polarity")
        log.error("P or 1 for p-side, 1 for holes polarity")
        sys.exit()

    #Open filename to read values                                                                                                    
    filename_tmp_lv =  open("lv_iv.txt",'r')
    #perform file operations                                                                                                         
    lv_12 = filename_tmp_lv.readline().strip()
    i_12   = filename_tmp_lv.readline().strip()
    lv_18 = filename_tmp_lv.readline().strip()
    i_18   = filename_tmp_lv.readline().strip()
    # closing temporary LV file                                                                                                      
    filename_tmp_lv.close()

    print("LV values: ")
    if (pol == 'N' or pol == '0'):
        print("LV_12_N [V]: " + str(lv_12) + " I_12_N [A]: " +str(i_12))
        print("LV_18_N [V]: " + str(lv_18) + " I_18_N [A]: " +str(i_18))
    else:
        print("LV_12_P [V]: " + str(lv_12) + " I_12_P [A]: " +str(i_12))
        print("LV_18_P [V]: " + str(lv_18) + " I_18_P [A]: " +str(i_18))

    return (lv_12, i_12, lv_18, i_18)

def powerOn_lv(pol_usr = 'N'):
    # Function to Turn On LV. According to the argument, the corresponding side is powered ON
    # Function to read the voltage and current of the FEB
    # LV mapping                                                                                                                     
    # LV channels vs FEB                                                                                                             
    # u904 ................ 1.2 V LDOs N-side                                                                                        
    # u905 ................ 1.8 V LDOs N-side                                                                                        
    # u906 ................ 1.2 V LDOs P-side                                                                                        
    # u907 ................ 1.8 V LDOs P-side                                                                                        

    pol = pol_usr

    if (pol == 'N' or pol == '0'):
        #reading LV in channels u904 and u905                                                                                         
        print("Turning ON for N-side (electrons polarity):")
        #Reading Voltage and current values and writting them in a file                                                               
        os.system("ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_u700 on")
        os.system("ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_u701 on")
        time.sleep(10)
        reading_lv(pol)
    elif(pol == 'P' or pol == '1'):
        print("Turning ON LV for P-side (holes polarity):")
        #Reading Voltage and current values and writting them in a file                                                               
        os.system("ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_u702 on")
        os.system("ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_u703 on")
        time.sleep(10)
        reading_lv(pol)
    else:
        log.error("Please, indicate a polarity for reading the corresponding LV potentials, in the following way:")
        log.error("N or 0 for n-side, 0 for electrons polarity")
        log.error("P or 1 for p-side, 1 for holes polarity")
        sys.exit()

def powerOn_EMU(lv_channel):
    # Function to Turn On EMU. According to the argument, the corresponding side is powered ON
    # Function to read the voltage and current of the FEB
    # LV mapping                                                                                                                     
    # LV channels vs FEB                                                                                                             
    # u800 ................ EMU_213                                                                                        
    # u801 ................ EMU_234                                                                                         
    # u802 ................ EMU_238
        
    log.info("Turning ON LV channel %s for EMU", lv_channel)
    #Reading Voltage and current values and writting them in a file
    command = "ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_{} on".format(lv_channel)
    os.system(command)
    time.sleep(15)
    command = "ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_{} status".format(lv_channel)
    os.system(command)
    
    
def powerOff_lv(pol_usr = 'N'):
    # Function to Turn Off LV. According to the argument, the corresponding side is powered ON
    # Function to read the voltage and current of the FEB
    # LV mapping                                                                                                                     
    # LV channels vs FEB                                                                                                             
    # u904 ................ 1.2 V LDOs N-side                                                                                        
    # u905 ................ 1.8 V LDOs N-side                                                                                        
    # u906 ................ 1.2 V LDOs P-side                                                                                        
    # u907 ................ 1.8 V LDOs P-side                                                                                        

    pol = pol_usr

    if (pol == 'N' or pol == '0'):
        #reading LV in channels u904 and u905                                                                                         
        print("Turning OFF for N-side (electrons polarity):")
        #Reading Voltage and current values and writting them in a file                                                               
        os.system("ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_u700 on")
        os.system("ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_u701 on")
        time.sleep(10)
        reading_lv(pol)
    elif(pol == 'P' or pol == '1'):
        print("Turning OFF LV for P-side (holes polarity):")
        #Reading Voltage and current values and writting them in a file                                                               
        os.system("ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_u702 on")
        os.system("ssh cbm@cbmflib01 ./mpod_lab2021/lab21_mpod02_u703 on")
        time.sleep(10)
        reading_lv(pol)
    else:
        log.error("Please, indicate a polarity for reading the corresponding LV potentials, in the following way:")
        log.error("N or 0 for n-side, 0 for electrons polarity")
        log.error("P or 1 for p-side, 1 for holes polarity")
        sys.exit()  

def powerON_hv(hv_n_channel, hv_p_channel, bias_voltage):
    # Function to Turn ON HV. According to the argument, the corresponding voltage is provided
    # HV mapping                                                                                                                     
    # HV channels vs SETUP                                                                                                             
    # u118 ................ N-side (Positive)    SETUP_0                                                                                        
    # u119 ................ P-side (Negativo)                                                                                        

    # u120 ................ N-side (Positive)    SETUP_1                                                                             
    # u121 ................ P-side (Negativo)
    
    # u118 ................ N-side (Positive)    SETUP_2                                                                                
    # u119 ................ P-side (Negativo)                                                                                        

    # Initizalizing the HV channels. Setting Output voltage = 0
    command = "ssh cbm@cbmflib01 snmpset -v 2c -m +WIENER-CRATE-MIB -c guru 10.203.0.64 outputVoltage.{} F {}".format(hv_n_channel, 0)
    os.system(command)
    command = "ssh cbm@cbmflib01 snmpset -v 2c -m +WIENER-CRATE-MIB -c guru 10.203.0.64 outputVoltage.{} F {}".format(hv_p_channel, 0)
    os.system(command)
    # Turning on the channels:
    command = "ssh cbm@cbmflib01 snmpset -v 2c -m +WIENER-CRATE-MIB -c guru 10.203.0.64 outputSwitch.{} i 1".format(hv_n_channel)
    os.system(command)
    command = "ssh cbm@cbmflib01 snmpset -v 2c -m +WIENER-CRATE-MIB -c guru 10.203.0.64 outputSwitch.{} i 1".format(hv_p_channel)
    os.system(command)
    #Ramping up the HV, up to the bias_voltage.
    #Symmetric mode:
    for hv in range(0, bias_voltage, 5):
        # n-side
        command = "ssh cbm@cbmflib01 snmpset -v 2c -m +WIENER-CRATE-MIB -c guru 10.203.0.64 outputVoltage.{} F {}".format(hv_n_channel, hv)
        os.system(command)
        # p-side
        hv *=(-1)
        command = "ssh cbm@cbmflib01 snmpset -v 2c -m +WIENER-CRATE-MIB -c guru 10.203.0.64 outputVoltage.{} F {}".format(hv_p_channel, hv)
        os.system(command)
        time.sleep(10)
    #reading HV currents @ the biasing voltage
    command = "ssh cbm@cbmflib01 snmpget -v 2c -m +WIENER-CRATE-MIB -c public 10.203.0.64 outputMeasurementCurrent.{} | grep -Eo '[0-9]+([.][0-9]+)'".format(hv_n_channel)
    hv_i_nside = subprocess.check_output(command, shell = True)
    hv_i_nside = hv_i_nside.decode("utf-8").strip()
    hv_i_nside = 1000000*float(hv_i_nside)
    command = "ssh cbm@cbmflib01 snmpget -v 2c -m +WIENER-CRATE-MIB -c public 10.203.0.64 outputMeasurementCurrent.{} | grep -Eo '[0-9]+([.][0-9]+)'".format(hv_p_channel)
    hv_i_pside = subprocess.check_output(command, shell = True)
    hv_i_pside = hv_i_pside.decode("utf-8").strip()
    hv_i_pside = 1000000*float(hv_i_pside)
    
    log.info("N-side: {} uA".format(hv_i_nside))
    log.info("P-side: {} uA".format(hv_i_pside))

    return (hv_i_nside, hv_i_pside) 
        
def powerOff_hv(hv_n_channel, hv_p_channel):
    # Function to Turn ON HV. According to the argument, the corresponding voltage is provided
    # HV mapping                                                                                                                     
    # HV channels vs SETUP                                                                                                             
    # u118 ................ N-side (Positive)    SETUP_0                                                                                        
    # u119 ................ P-side (Negativo)                                                                                        

    # u120 ................ N-side (Positive)    SETUP_1                                                                             
    # u121 ................ P-side (Negativo)
    
    # u118 ................ N-side (Positive)    SETUP_2                                                                                
    # u119 ................ P-side (Negativo)                                                                                        
    thr_v_off = 5                                # Voltage threshold to consider a channel off
    
    # Turning off the channels:
    command = "ssh cbm@cbmflib01 snmpset -v 2c -m +WIENER-CRATE-MIB -c guru 10.203.0.64 outputSwitch.{} i 0".format(hv_n_channel)
    os.system(command)
    command = "ssh cbm@cbmflib01 snmpset -v 2c -m +WIENER-CRATE-MIB -c guru 10.203.0.64 outputSwitch.{} i 0".format(hv_p_channel)
    os.system(command)
    time.sleep(180)
    
    # reading voltage after OFF 
    command = "ssh cbm@cbmflib01 snmpget -v 2c -m +WIENER-CRATE-MIB -c public 10.203.0.64 outputMeasurementTerminalVoltage.{} | grep -Eo '[0-9]+([.][0-9]+)'".format(hv_n_channel)
    hv_v_nside = subprocess.check_output(command, shell = True)
    hv_v_nside = hv_v_nside.decode("utf-8").strip()
    hv_v_nside = float(hv_v_nside)
    command = "ssh cbm@cbmflib01 snmpget -v 2c -m +WIENER-CRATE-MIB -c public 10.203.0.64 outputMeasurementTerminalVoltage.{} | grep -Eo '[0-9]+([.][0-9]+)'".format(hv_p_channel)
    hv_v_pside = subprocess.check_output(command, shell = True)
    hv_v_pside = hv_v_pside.decode("utf-8").strip()
    hv_v_pside = float(hv_v_pside)
    
    if (hv_v_nside < thr_v_off and hv_v_pside < thr_v_off):
        return True
    else:
        return False
        

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ----------------------------- MODULE TESTING: DIRECTORY AND FILES ---------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Module testing directory & files                                                                                                                                                                                                           
def read_moduleId():
    module_sn = input("INTRODUCE THE MODULE NAME: ")
    #module_sn = "L3DL3 01Tr-6-PA2 42 C"
    return module_sn
    
def check_moduleId(module_str):
    module_id = ""
    sensor_size = ""
    sensor_qgrade = ""
    module_sn = module_str.split()
    if (len(module_sn) == 0):
        print("MODULE NAME should not be left empty")
        sys.exit()
    elif (len(module_sn) == 1):
        module_id = module_sn[0]
        if (module_id[-2:-1] == 'A' or module_id[-2:-1] == 'B'):
            sensor_size = input("INTRODUCE SENSOR SIZE [mm]:")
            sensor_qgrade = input("INTRODUCE SENSOR QUALITY GRADE (ex: A, B, C, D): ")
            return (module_id, sensor_size, sensor_qgrade)
        else:
            print("MODULE ID must contain the FEB type (A/B) in the second-to-last position")
            return ['na']
    else:
        module_id = module_sn[1]
        sensor_size = module_sn[2]
        sensor_qgrade = module_sn[3]
        if (module_id[-2:-1] == 'A' or module_id[-2:-1] == 'B'):
            sensor_size = module_sn[2]
            sensor_qgrade = module_sn[3]
            return (module_id, sensor_size, sensor_qgrade)
        else:
            print("MODULE ID must contain the FEB type (A/B) in the second-to-last position")
            return ['na']
    
def read_test_center():
    #test_center = input("INTRODUCE THE TEST CENTER ([1 or GSI, 2 or KIT): ")
    test_center = "GSI"
    test_center_list = ['1', '2', 'GSI', 'KIT','gsi', 'kit']
    if (test_center in test_center_list):
        if (test_center =='1' or test_center == 'gsi'):
            test_center = 'GSI'
        if (test_center =='2' or test_center == 'kit'):
            test_center = 'KIT'
        return test_center
    else:
        log.error("TEST CENTER NAME is not included in the list")
        return 'na'

def read_operator_id():
    #operator_id = input("Insert the Operator ID (example: Jane Doe): ")
    operator_id = "Dairon Rodriguez Garces"
    operator_id_list = []
    with open("operators_file.txt") as file:
        for line in file:
            operator_id_list.append(line.rstrip())
    file.close()
    
    if (operator_id in operator_id_list):
        return operator_id
    else:
        #add_to_list = input("OPERATOR'S ID is not found. Should be added to the list (Y/N):")
        add_to_list = "Y"
        if (add_to_list == 'Y' or 'y' or 'yes'):
           operator_id_file = open("operators_file.txt", "a")
           operator_id_file.write(operator_id)
        else:
            return 'na'

def initWorkingDirectory(module_dir = "module_files/M00", module_sn ='M00'):
    # Creating the directory of module testing 
    # importing the module name  
    # module_sn = "M" + module_sn
    # module_path = module_path
    # module_dir  = module_path + '/' + module_sn
    module_dir  = module_dir
    # other variables:                                   
    date  = time.strftime("%y%m%d-%H:%M:%S")
   # checking if directory exists                       
    if (os.path.isdir(module_dir) == True):
        print("This module directory already exists")
    else:
        # creating the module files directory             
        try:
            os.makedirs(module_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
            else:
                print("Succesfully created Module directory %s:" + (module_dir))

    # creating log file                                  
    logfile = open(module_dir + "/" + module_sn + "_log.log","a")
    
    # checking if data file exists                       
    if (os.path.isfile(module_dir + "/" + module_sn + "_data.dat") == True):
        print("Data file already exists. Would you like to re-write it (Y/N):")
        #writing_flag = input()
        writing_flag = "N"
        if ( writing_flag == 'N' or writing_flag == 'No' or writing_flag == 'n' or writing_flag == 'no'):
            module_sn = module_sn + "_" + time.strftime("%y%m%d_%H%M")       # changing file name with different date                                        
            print("Data file will be created with the following name: %s" %(module_sn + "_data.dat"))
        else:
            module_sn = module_sn
    else:
        module_sn = module_sn
        print("Data file does not exists & it will be created with the following name: %s"  %( module_sn + "_data.dat"))

    # creating data file and initializing it
    datafile = open(module_dir + "/" + module_sn + "_data.dat","w+")
    datafile.write("FILENAME_DATA: ")
    datafile.write(module_sn)
    datafile.write("\n")
    datafile.write("TEST_DATE: ")
    datafile.write(date)
    datafile.write("\n")
    datafile.close()

    # initializing log file                        
    logfile.write(date)
    logfile.write("\t")
    logfile.write("MODULE_ID:")
    logfile.write("\t")
    logfile.write(module_sn)
    logfile.write("\n")    
    logfile.close()

    return module_sn

def close_log_file(module_sn):
    # ending testing & closing data and log files
    logfile.write(date)
    logfile.write("\t")
    logfile.write("Ending test sequence")
    logfile.write("\n")
    logfile.close()
    datafile.close()                                                                    

def write_data_file(module_dir, module_sn = "M00", info = ""):
    datafile = open(module_dir + "/" + str(module_sn) + "_data.dat","a")
    datafile.write(str(info))
    datafile.write("\n")
    datafile.close()

def write_log_file(module_dir, module_sn = "M00", info = ""):
    date = time.strftime("%y%m%d-%H:%M:%S")
    logfile = open(module_dir + "/" + str(module_sn) + "_log.log","a")
    logfile.write(date)
    logfile.write("\t")
    logfile.write("[INFO]")
    logfile.write("\t")
    logfile.write(str(info))
    logfile.write("\n")
    logfile.close()
 
def making_pscan_dir(module_dir):
    # Creating the pscan_files directory
    # checking if directory exists
    print("Creating Pscan files directory")
    pscan_dir = module_dir + "/pscan_files/"
    if (os.path.isdir(pscan_dir) == True):
        print("This module directory already exists")
        return pscan_dir
    else:
        # creating the module files directory             
        try:
            os.makedirs(pscan_dir)
            return pscan_dir
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
            else:
                print("Succesfully created directory: %s" + (pscan_dir))


def making_trim_dir(module_dir):
    # Creating the trim_files directory
    # checking if directory exists
    print("Creating Trim files directory")
    trim_dir = module_dir + "/trim_files/"
    if (os.path.isdir(trim_dir) == True):
        print("This module directory already exists")
        return trim_dir
    else:
        # creating the module files directory             
        try:
            os.makedirs(trim_dir)
            return trim_dir
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
            else:
                print("Succesfully created directory: %s" + (trim_dir))

def making_conn_check_dir(module_dir):
    # Creating the connection_files directory                                                                                                                        
    # checking if directory exists                                                                                                                                                                  
    print("Creating Connection files directory")
    conn_check_dir = module_dir + "/conn_check_files/"
    if (os.path.isdir(conn_check_dir) == True):
        print("This module directory already exists")
        return conn_check_dir
    else:
        # creating the module files directory                                                                                                                                                       
        try:
            os.makedirs(conn_check_dir)
            return conn_check_dir
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
            else:
                print("Succesfully created directory: %s" + (conn_check_dir))


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ----------------------------- MODULE TESTING: OPERATING FUNCTIONS ---------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def read_FebSN_Nside(module_sn = 'M00'):
    # Fucntion to read the FEB serial number
    feb_sn = input("--> INTRODUCE THE SERIAL NUMBER for the N-side FEB. (XXXType<A/B>NumberOfUplinks ex: 138A2): ")
    #feb_sn = "2099B2"
    if (feb_sn[-2:-1] == 'A' or feb_sn[-2:-1] == 'B'):
        if (feb_sn[-2:-1] != module_sn[-2:-1]):
            return feb_sn
        else:
            log.error("Please, check that the SN of the N-side FEB is the right one. It contradicts the module ID")
            return 'na'
    else:
        log.error("FEB ID must contain the FEB type (A/B) in the second-to-last position")
        return 'na'

def read_FebSN_Pside(module_sn = 'M00'):
    # Fucntion to read the FEB serial number
    feb_sn = input("--> INTRODUCE THE SERIAL NUMBER for the P-side FEB. (XXXType<A/B>NumberOfUplinks ex: 138B2): ")
    #feb_sn = "1100A2"
    if (feb_sn[-2:-1] == 'A' or feb_sn[-2:-1] == 'B'):
        if (feb_sn[-2:-1] == module_sn[-2:-1]):
            return feb_sn
        else:
            log.error("Please, check the FEB_SN is the right one. It contradicts the module ID")
            return 'na'
    else:
        log.error("FEB ID must contain the FEB type (A/B) in the second-to-last position")
        return 'na'
        
def read_asicIDs_FEB(smx_l_side, pol, feb_type):
    # Function to read the overall ASIC ID (FEB where ASIC belongs, polarity, HW address, ASIC e-fuse ID(string) and (int)) 
    info = ""
    feb_type_sw = []
    write_data_file(module_dir, module_sn_tmp, info)
    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
        pol_calib = 0
        info = 'EFUSE_ID_N'
    else:
        pol_str = 'P-side'
        pol_calib = 1
        info = 'EFUSE_ID_P'
    write_data_file(module_dir, module_sn_tmp, info)
    header_line  = "FEB-ID_\t\t_POLARITY_\t\t_HW-ADDR_\t\t_EFUSE-ID-(STR)_\t\t_EFUSE-ID-(INT)_"
    log.info(header_line)
    write_data_file(module_dir, module_sn_tmp, header_line)

    if (feb_type[-2:-1] == "A"):
        feb_type_sw.extend(feb_type_sw_A)
    else:
        feb_type_sw.extend(feb_type_sw_B)

    smx_counter =0
    for asic_sw in feb_type_sw:
        for smx in smx_l_side:
            if (smx.address == asic_sw):
                addr = smx.address
                asic_id_int = smx.read_efuse()
                asic_id_str = smx.read_efuse_str()
                info  = "{} \t\t {} \t\t {} \t\t {} \t\t {}".format(feb_type, pol_str, addr, asic_id_str, asic_id_int)
                write_data_file(module_dir, module_sn_tmp, info)
                log.info(info)
            else:
                pass
    return 0

                
def read_VDDM_TEMP_FEB(smx_l_side, pol, feb_type):
    # Function to read the VDDM and TEMP of the ASICs
    info = ""
    feb_type_sw = []
    write_data_file(module_dir, module_sn_tmp, info)
    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
        info = 'VDDM_TEMP_N'
    else:
        pol_str = 'P-side'
        info = 'VDDM_TEMP_P'
    write_data_file(module_dir, module_sn_tmp, info)
    header_line  = "FEB-ID_\t\t_POLARITY_\t\t_HW-ADDR_\t_VDDM POTENTIAL [LSB] | [mV]_\t\t_TEMP [mV] | [C]_"    
    log.info(header_line)
    write_data_file(module_dir, module_sn_tmp,header_line)

    if (feb_type[-2:-1] == "A"):
        feb_type_sw.extend(feb_type_sw_A)
    else:
        feb_type_sw.extend(feb_type_sw_B)

    for asic_sw in feb_type_sw:
        for smx in smx_l_side:
            if (smx.address == asic_sw):
                asic_vddm = smx.read_vddm()
                #asic_vddm_src = smx.read_diag("Vddm")
                asic_temp = smx.read_temp()
                #asic_temp_src = smx.read_diag("Temp")
                info = "{} \t\t {} \t\t {} \t\t\t {} \t {:.1f} \t\t\t {:.1f} \t {:.1f}".format(feb_type, pol_str, smx.address, asic_vddm[0], asic_vddm[1], asic_temp[0], asic_temp[1])
                write_data_file(module_dir, module_sn_tmp,info)
                log.info(info)
            else:
                pass
    return 0

def load_STD_Config(smx_l_side, pol, feb_type):
    # Function to write on each ASIC the default settings
    feb_type_sw = []
    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
        pol_calib = 0 
    else:
        pol_str = 'P-side'
        pol_calib = 1

    if (feb_type[-2:-1] == "A"):
        feb_type_sw.extend(feb_type_sw_A)
    else:
        feb_type_sw.extend(feb_type_sw_B)

    for asic_sw in feb_type_sw:
        for smx in smx_l_side:
            if (smx.address == asic_sw):
                header_line  = "--> SETTING STANDARD CONFIGURATION for ASIC with HW address {} and polarity {}".format(smx.address,pol_str)
                log.info(header_line)
                smx.write_def_ana_reg(smx.address, pol_calib )
                smx.read_reg_all(compFlag = False)
            else:
                pass
    return 0
            
def set_Trim_default(smx_l_side, pol, feb_type, cal_asic_list):
    feb_type_sw = []
    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
    else:
        pol_str = 'P-side'

    if (feb_type[-2:-1] == "A"):
        feb_type_sw.extend(feb_type_sw_A)
    else:
        feb_type_sw.extend(feb_type_sw_B)

    for asic_sw in feb_type_sw:
        for smx in smx_l_side:
            if ((smx.address == asic_sw) and (smx.address in cal_asic_list)):
                header_line  = "--> SETTING DEFAULT TRIM for ASIC with HW address {} and polarity {}".format(smx.address, pol_str)
                log.info(header_line)
                smx.set_trim_default(128,36)
            else:
                #info = "--> NO SETTING DEFAULT TRIM for ASIC with HW ADDRESS {} in {}".format(asic_hw_addr, pol_str)
                #log.info(info)
                pass
    return 0
        
def scan_VrefP_N_Thr2glb(smx_l_side, pol, feb_type, cal_asic_list, npulses = 100, test_ch = 64, amp_cal_min = 30, amp_cal_max = 247, amp_cal_fast = 30, vref_t = 118):
    smx_cnt = 0
    feb_type_sw = []
    #cal_set_asic = [0 for n_asics in len(smx_l_side)[0 for vref in range(3)]]
    cal_set_asic = []
    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
        pol_calib = 0 
    else:
        pol_str = 'P-side'
        pol_calib = 1

    if (feb_type[-2:-1] == "A"):
        feb_type_sw.extend(feb_type_sw_A)
    else:
        feb_type_sw.extend(feb_type_sw_B)

    for asic_sw in feb_type_sw:
        for smx in smx_l_side:
            if ((smx.address == asic_sw) and (smx.address in cal_asic_list)):
                header_line  = "--> SCANNING VREF_P,N & THR2_GLB for ASIC with HW address {} and polarity {}".format(smx.address, pol_str)
                log.info(header_line)
                cal_set_asic.append(smx.vrefpn_scan(pol_calib, test_ch, npulses, amp_cal_min, amp_cal_max, amp_cal_fast, vref_t))
            else:
                #info = "--> NO SCANNING VREF_P/N & THR2_GLB FOR ASIC with HW ADDRESS {} in {}".format(asic_hw_addr, pol_str)
                #log.info(info)
                pass
    return cal_set_asic

def print_cal_settings(pol, cal_set_side, cal_asic_list, vref_t_calib = 118):
    # Publishing the VREF_P_N & Thr2_glb results
    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
    else:
        pol_str = 'P-side'
    info = "--> CALIBRATION SETTINGS for {}".format(pol_str)
    log.info(info)
    header_line = "HW_address \t VRef_P \t VRef_N \t VRef_T \t Thr2_glb"
    log.info(header_line)
    
    smx_cnt = 0
    for asic in cal_set_side:
        info = "{}\t{}\t{}\t{}\t{}".format(cal_asic_list[smx_cnt], asic[0], asic[1], vref_t_calib,  asic[2])
        log.info(info)
        smx_cnt +=1
    return 0

def writing_cal_settings(smx_l_side, pol, feb_type,  cal_set_side, cal_asic_list, vref_t_calib = 118):
    vref_n_arr = []
    vref_p_arr = []
    thr2_glb_arr = []
    feb_type_sw = []
    #if (len(smx_l_side)!= len(cal_set_side)):
    #log.error("Length of the calibration settings array and the number of ASICs does not coincide")
    #else:
    info = ""
    write_data_file(module_dir, module_sn_tmp, info)
    for asic in cal_set_side:
        vref_p_arr.append(asic[0])
        vref_n_arr.append(asic[1])
        thr2_glb_arr.append(asic[2])
    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
        info = "CAL_SETTINGS_N"
    else:
        pol_str = 'P-side'
        info = "CAL_SETTINGS_P"
    log.info(info)
    write_data_file(module_dir, module_sn_tmp, info)
    info = "INDEX: \t\t HW_ADDR: \t VRef_P: \t VRef_N: \t VRef_T: \t VRef_T_range \t Thr2_glb:"
    write_data_file(module_dir, module_sn_tmp, info)

    if (feb_type[-2:-1] == "A"):
        feb_type_sw.extend(feb_type_sw_A)
    else:
        feb_type_sw.extend(feb_type_sw_B)

    smx_cnt = 0
    for asic_sw in feb_type_sw:
        for smx in smx_l_side:
            if ((smx.address == asic_sw) and (smx.address in cal_asic_list)):
                smx.write(130, 9, vref_p_arr[smx_cnt])
                smx.write(130, 8, vref_n_arr[smx_cnt])
                smx.write(130, 18, vref_t_calib)
                smx.write(130, 7, thr2_glb_arr[smx_cnt])
                vreft_range  = (smx.read(130,10)&64)>>4 or (smx.read(130,18)&192)>>6
                # -----------------------------------------------------------------------
                info = "{}\t\t {}\t\t {}\t\t {}\t\t {}\t\t {}\t\t {}".format(smx_cnt, smx.address, smx.read(130,9)&0xff, smx.read(130,8)&0xff, smx.read(130,18)&0xff, vreft_range, smx.read(130,7)&0xff)
                log.info(info)
                write_data_file(module_dir, module_sn_tmp, info)
                smx_cnt +=1
            else:
                pass
            
# To check calibration function
def calib_FEB(smx_l_side, trim_dir, pol, feb_type, cal_asic_list, npulses = 40, amp_cal_min = 30, amp_cal_max = 247, amp_cal_fast = 30, much_mode_on = 0):
    # Function to calibrate the ADC and FAST discriminator of each ASIC according to the given polarity
    info = ""
    feb_type_sw = []
    write_data_file(module_dir, module_sn_tmp, info)
    filename_trim = trim_dir
    pol_calib = 0
    if (pol == 'N' or pol == '0'):
        pol_str = 'elect'
        pol_calib = 0
        info = 'TRIM_FILE_N'
    elif (pol == 'P' or pol == '1'):
        pol_str = 'holes'
        pol_calib = 1
        info = 'TRIM_FILE_P'
    else:
        log.error("Please check the polarity is correct: (ex: string 'N' or '0', 'P' or '1')")
    write_data_file(module_dir, module_sn_tmp, info)
    # definition of the final trim array
    trim_final = [[0 for d in range(32)] for ch in range(128)]    
    if (feb_type[-2:-1] == "A"):
        feb_type_sw.extend(feb_type_sw_A)
    else:
        feb_type_sw.extend(feb_type_sw_B)

    for asic_sw in feb_type_sw:
        for smx in smx_l_side:
            if ((smx.address == asic_sw) and (smx.address in cal_asic_list)):
                asic_hw_addr = smx.address
                # Launching the calibration with the corresponding
                info = "--> RUNNING CALIBRATION FOR ASIC with HW ADDRESS {} in {}".format(asic_hw_addr, pol_str)
                log.info(info)
                # Elements for the filename
                asic_id_str = smx.read_efuse_str()
                vref_n = smx.read(130,8)&0xff
                vref_p = smx.read(130,9)&0xff
                vref_t = smx.read(130,18)&0xff
                thr2_glb = smx.read(130,7)&0xff
                date_fw = datetime.now().strftime("%y%m%d_%H%M")
                filename_str = 'ftrim_' + asic_id_str + '_HW_' + str(asic_hw_addr) + '_SET_'+ str(vref_p) + '_' + str(vref_n) + '_' + str(vref_t) + '_' + str(thr2_glb) + '_R_' + str(amp_cal_min) + '_' + str(amp_cal_max) + '_' + pol_str
                filename_trim = trim_dir + filename_str
                # Executing calibration
                log.info("Calib file name: {}".format(filename_trim))
                smx.get_trim_adc_SA(pol_calib, trim_final, 40, amp_cal_min, amp_cal_max, much_mode_on)
                #smx.get_trim_adc(pol_calib, trim_final, 40, amp_cal_min, amp_cal_max,vref_t,  much_mode_on)
                smx.get_trim_fast(pol_calib, trim_final, npulses, amp_cal_fast, much_mode_on)
                # Writing calibration file
                smx.write_trim_file(filename_trim, pol_calib, trim_final, amp_cal_min, amp_cal_max, amp_cal_fast, much_mode_on)
                info = "CAL_ASIC_HW_ADDR_{}: {}.txt".format(asic_hw_addr,filename_str)
                write_data_file(module_dir, module_sn_tmp, info)
                write_log_file(module_dir, module_sn, info)            
            else:
                #info = "SKIP_ASIC_HW_ADDR_{}".format(asic_hw_addr)
                #log.info(info)
                #write_log_file(module_dir, module_sn, info)
                pass
    return 0

def set_trim_calib(smx_l_side, trim_dir, pol, feb_type, cal_asic_list, much_mode_on = 0):
    # Function to set the calibration values the ADC and FAST discriminator of each ASIC according to the given polarity and the ASIC ID                            
    filename_trim = trim_dir
    feb_type_sw = []
    pol_calib = 0
    if (pol == 'N' or pol == '0'):
        pol_str = 'elect'
        pol_calib = 0
    elif (pol == 'P' or pol == '1'):
        pol_str = 'holes'
        pol_calib = 1
    else:
        log.error("Please check the polarity is correct: (ex: string 'N' or '0', 'P' or '1')")

    if (feb_type[-2:-1] == "A"):
        feb_type_sw.extend(feb_type_sw_A)
    else:
        feb_type_sw.extend(feb_type_sw_B)

    for asic_sw in feb_type_sw:
        for smx in smx_l_side:
            if ((smx.address == asic_sw) and (smx.address in cal_asic_list)):
                asic_hw_addr = smx.address
                # Setting the TRIM calibration values                                                                                                                  
                info = "--> SETTING TRIM CALIBRATION VALUES FOR ASIC with HW ADDRESS {} in {}".format(asic_hw_addr, pol_str)
                log.info(info)
                # Elements for the trim file
                asic_id_str = smx.read_efuse_str()
                smx.set_trim(trim_dir, pol_calib, asic_id_str)
            else:
                pass
    return 0

def check_trim(smx_l_side, pscan_dir, pol, feb_type, cal_asic_list, disc_list = [5,10,16,24,30,31], vp_min = 0, vp_max = 255, vp_step = 1, npulses = 100):
    # Function to measure the ADC and FAST discriminator response function for a  given number of discriminators of an ASIC                                                          
    info = ""
    feb_type_sw = []
    write_data_file(module_dir, module_sn_tmp, info)
    pol_calib = 0
    disc_list = disc_list
    npulses = npulses
    vp_min = vp_min
    vp_max = vp_max
    vp_step = vp_step
    if (pol == 'N' or pol == '0'):
        pol_str = 'elect'
        pol_calib = 0
        info = 'PSCAN_FILE_N'
    elif (pol == 'P' or pol == '1'):
        pol_str = 'holes'
        pol_calib = 1
        info = 'PSCAN_FILE_P'
    else:
        log.error("Please check the polarity is correct: (ex: string 'N' or '0', 'P' or '1')")
    write_data_file(module_dir, module_sn_tmp, info)

    if (feb_type[-2:-1] == "A"):
        feb_type_sw.extend(feb_type_sw_A)
    else:
        feb_type_sw.extend(feb_type_sw_B)

    for asic_sw in feb_type_sw:
        for smx in smx_l_side:
            asic_hw_addr = smx.address
            if ((smx.address == asic_sw) and (smx.address in cal_asic_list)):
                # Checking ENC and calibration results for an ASIC                                                                                                     
                info = "PSCAN_ASIC_HW_ADDR_{}: {}".format(asic_hw_addr, pol_str)
                log.info(info)
                # Elements for the pscan file                                                                                                                          
                asic_id_str = smx.read_efuse_str()
                pscan_filename = smx.check_trim_red(pscan_dir, pol_calib, asic_id_str, disc_list, vp_min, vp_max, vp_step, npulses)
                info = "PSCAN_ASIC_HW_ADDR_{}: {}".format(asic_hw_addr, pscan_filename)
                write_data_file(module_dir, module_sn_tmp, info)
                write_log_file(module_dir, module_sn, info)
            else:
                pass
                #info = "SKIP_PSCAN_ASIC_HW_ADDR_{}".format(asic_hw_addr)
                #log.info(info)
                #write_data_file(module_dir, module_sn_tmp, info)
                #write_log_file(module_dir, module_sn, info)
    return 0


def connection_check(smx_l_side, conn_check_dir, pol, feb_type, cal_asic_list, nloops = 5, vref_t = 108 ):
    # Function to check the connectivy of each channel by counting noise hits at a lower threshold
    info = ""
    feb_type_sw = []
    write_data_file(module_dir, module_sn_tmp, info)
    if (pol == 'N' or pol == '0'):
        pol_str = 'elect'
        pol_calib = 0
        info = 'CH_CONN_CHECK_FILE_N'
    elif (pol == 'P' or pol == '1'):
        pol_str = 'holes'
        pol_calib = 1
        info = 'CH_CONN_CHECK_FILE_P'
    else:
        log.error("Please check the polarity is correct: (ex: string 'N' or '0', 'P' or '1')")
    write_data_file(module_dir, module_sn_tmp, info)
    
    if (feb_type[-2:-1] == "A"):
        feb_type_sw.extend(feb_type_sw_A)
    else:
        feb_type_sw.extend(feb_type_sw_B)

    log.info("FEB type: {}".format(feb_type[-2:-1]))
    for asic_sw in feb_type_sw:
        for smx in smx_l_side:
            #asic_hw_addr = smx.address
            if ((smx.address == asic_sw) and (smx.address in cal_asic_list)):
                # Checking ENC and calibration results for an ASIC     
                info = "CONN-CHECK_ASIC_HW_ADDR_{}: {}".format(smx.address, pol_str)
                log.info(info)
                write_data_file(module_dir, module_sn_tmp, info)
                write_log_file(module_dir, module_sn, info)
                # Elements for the connection check file                                                                                                                     
                #asic_id_str = smx.read_efuse_str()
                smx.connection_check(conn_check_dir, pol_calib, nloops, vref_t)
            else:
                pass
            #    info = "SKIP_CONN_ASIC_HW_ADDR_{}".format(asic_hw_addr)
            #    log.info(info)
            #    write_data_file(module_dir, module_sn_tmp, info)
            #    write_log_file(module_dir, module_sn, info)
    return 0
            

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ----------------------------- MODULE TESTING: VARIABLES DEFINITION --------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
emu_list = ["EMU_234"]                           #["EMU_213"]     # List of EMU boards used during the test
module_path  = "module_files"                    # Path to store all the information concerning module testing
module_sn = " "                                  # Initiialization of the module SN

# List of ASICs HW address to be tested in case the 
#cal_asic_list = [0, 1, 2, 3, 4, 5, 6, 7]
cal_asic_list_nside = [0, 1, 2, 3, 4, 5, 6, 7]
#cal_asic_list_nside = []
cal_asic_list_pside = [0, 1, 2, 3, 4, 5, 6, 7]
#cal_asic_list_pside = []
#smx_l_short = smx_l[0:1]
#smx_cnt = 0
lv_channel_emu = 'u800'
hv_n_channel = 'u118'
hv_p_channel = 'u119'
bias_voltage = 80                                # Symmetric bias voltage. +/- 75V 


# Operating variables GetVrefs paramaters
npulses = 100
test_ch = 64
# Global variables used everywhere
feb_type_sw_B = [1, 0, 3, 2, 5, 4, 7, 6]
feb_type_sw_A = [7, 6, 5, 4, 3, 2, 1, 0]
# Global variables, used only during the get_trim function
# Calibration settings:  ADC Disc 30 -> amp_cal_min, ADC Disc 0 -> amp_cal_max, FAST Disc -> amp_cal_fast
amp_cal_min = 30
amp_cal_max = 247
amp_cal_fast = 30
vref_t = 118                    # Vrf_T value  = 54 in the largest VRef_T range. To detemrine calib par and to calibrate ADC
# Check Trim Parameters
disc_list = [5,10,16,24,30,31]
#disc_list = [15,18,20,25,30,31]
vp_min = 0
vp_max = 255
vp_step = 1
# Connection check Parameters
vref_t_low = 102                # Low ADC threshold to count noise hits in the discriminators 
nloops = 5                      # Number of loops to count noise hits in discriminators
        
#Lists of subsequences for module tests
# ------------------------------------------------------------------------------------------------------------------
# Possible test sequences.

test_list_init = ["#power_on_emu", "full_sync","#turn_hv_on", "read_lv_bc"]
test_list_comm = ["std_config","read_asic_id","set_trim_default", "read_lv_ac", "check_vddm_temp"]
test_list_calib = ["set_trim_default", "get_vrefs", "set_calib_par","get_trim"]
test_list_check = ["set_trim_calib", "check_trim", "#turn_hv_off", "#conn_check"]
test_list_stress = ["#reg_config_stress", "#iv_meas", "#set_mbias", "long_run"]
#test_list_rdata = []
#test_list_sysoff = ["#set_lv_off", "#set_hv_off"]


test_list = []
# Add the list of test to perfrom
test_list.extend(test_list_init)
test_list.extend(test_list_comm)
test_list.extend(test_list_calib)
test_list.extend(test_list_check)
#test_list.extend(test_list_stress)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ---------------------------------- MODULE TESTING: RUNNING  ---------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Step 0: -------------- Creating the working directory -------------------
module_sn = 'na'
module_str = []
nfails = 0
while(module_sn == 'na' and nfails < 3):
    module_str.extend(check_moduleId(str(read_moduleId())))
    module_sn = module_str[0]
    print("Module_str: {}".format(module_str))
    print("Module_SN: {}".format(module_sn))
    if (module_sn == 'na'):
        module_str.clear()
        nfails+=1       
    if (nfails ==3):
        print("Multiple fails on Writing Module ID. It should contain A or B in the second-to-last position. Please check the Module's information")
        sys.exit()
    else:
        pass

module_dir = module_path + "/" + str(module_sn)
module_sn_tmp = initWorkingDirectory(module_dir, module_sn)
pscan_dir = making_pscan_dir(module_dir)
trim_dir = making_trim_dir(module_dir)
conn_check_dir = making_conn_check_dir(module_dir)
# Setting logging directory
set_logging_details(module_dir)


# Step 0: --------------- Initiazlizing data & log files -----------------

info = "TEST_CENTER: {}".format(read_test_center())
write_data_file(module_dir, module_sn_tmp, info)
write_log_file(module_dir, module_sn, info)
info = "OPERATOR_ID: {}".format(read_operator_id())
write_data_file(module_dir, module_sn_tmp, info)
write_log_file(module_dir, module_sn, info)

info = "MODULE ID: \t{}".format(module_sn)
write_data_file(module_dir, module_sn_tmp, info)
write_log_file(module_dir, module_sn, info)
info = "SENSOR_SIZE [mm]: \t{}".format(module_str[1])
write_data_file(module_dir, module_sn_tmp, info)
write_log_file(module_dir, module_sn, info)
info = "SENSOR_QGRADE: \t{}".format(module_str[2])
write_data_file(module_dir, module_sn_tmp, info)
write_log_file(module_dir, module_sn, info)

# Step 1.1: -------------- Turning ON LV for N & P-side -------------------
#powerOn_lv('N')
#powerOn_lv('P')
#time.sleep(5)
#lv_nside_bc = reading_lv('N')
#lv_pside_bc = reading_lv('P')

# Step 1.2: ---------------- Setting active downlinks ---------------------
# Setting active downlinks according to the module type.
# Since N-side should be connected on the slot 0 of the EMU_FMC, the active downlinks
# can be determined for this configuration by knowing the Module ID. The Module ID
# contains the P-side FEB type

active_downlinks = []

if (module_sn[-2:-1] == 'A'):
    active_downlinks.append(0)
    active_downlinks.append(3)
elif (module_sn[-2:-1] == 'B'):
    active_downlinks.append(1)
    active_downlinks.append(2)
else:
    for i in range(0,4):
        active_downlinks.append(i)

# Setp 2: ------------------ Checking the FEBs ID -------------------------
# System exit if nfails = 3
feb_nside ='na'
nfails = 0
if (nfails < 3):
    while(feb_nside == 'na'):
        feb_nside = read_FebSN_Nside(module_sn)
        if(feb_nside == 'na'):
            nfails+=1
else:
    sys.exit()
info = "FEB_SN_N:\t {}".format(feb_nside) 
write_data_file(module_dir, module_sn_tmp, info)

feb_pside ='na'
nfails = 0
if (nfails < 3):
    while(feb_pside =='na'):
        feb_pside = read_FebSN_Pside(module_sn)
        if (feb_pside =='na'):
            nfails+=1
else:
    sys.exit()    
info = "FEB_SN_P:\t {}".format(feb_pside)
write_data_file(module_dir, module_sn_tmp, info)

# TEMP_ARR. Fix HV readout
hv_current_n = input("Enter sensor leakage current N-side [uA]: ")
info = "I_SENSOR_150V_N: {} [uA]".format(hv_current_n)
write_data_file(module_dir, module_sn_tmp, info)

hv_current_p = input("Enter sensor leakage current P-side [uA]: ")
info = "I_SENSOR_150V_P: {} [uA]".format(hv_current_p)
write_data_file(module_dir, module_sn_tmp, info)


# RUNNING TEST SEQUENCE
# -------------------------------------------
for test_step in test_list:
    if (test_step == "power_on_emu"):
        log.info(" ------------------- POWERING UP EMU BOARD ---------------------- ")
        info = "-->> POWERING UP EMU BOARD"
        write_log_file(module_dir, module_sn, info)
        # Step 1: -------------- Turning ON LV for EMU board -------------------
        powerOn_EMU("u801")
        time.sleep(10)

    elif (test_step =="full_sync"):
        log.info(" ------------------- RUNNING SYNCHRONIZATION ---------------------- ")
        info = "-->> RUNNING SYNCHRONIZATION"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Step 2: -------- Starting with connection and synchronization ----------- 
        info = "EMU_BOARD_SN: {}".format(emu_list[0])
        write_data_file(module_dir, module_sn_tmp, info)
        smx_l = general_sync(emu_list, active_downlinks)
        # 2.1 Determining the number of ASICs per side
        n_asic_all = scanning_asics(smx_l)
        # 2.2 Assigning the ASICs according to polarities
        n_asics = n_asic_all[0]
        p_asics = n_asic_all[1]
        write_data_file(module_dir, module_sn_tmp, "")
        info = "No_SYNC_N:\t {}".format(n_asics)
        write_data_file(module_dir, module_sn_tmp, info)
        info = "No_SYNC_P:\t {}".format(p_asics)
        write_data_file(module_dir, module_sn_tmp, info)
        if (n_asics!=0):
            info = "DWN_LINK_N: {}".format(active_downlinks[0])
        else:
            info = "DWN_LINK_N: -1"
        write_data_file(module_dir, module_sn_tmp, info)
        if (p_asics!=0):
            info = "DWN_LINK_P: {}".format(active_downlinks[1])
        else:
            info = "DWN_LINK_P: -1"
        write_data_file(module_dir, module_sn_tmp, info)
        # 2.3 Distributing the ASICs in arrays according to operational polarity 
        smx_l_nside = smx_l[0:n_asics]
        smx_l_pside = smx_l[n_asics:n_asics + p_asics]
        info = "<<-- FINISHED SYNCHRONIZATION"
        write_log_file(module_dir, module_sn, info)

    elif(test_step == "turn_hv_on"):
        log.info(" ---------------------------- TURNING ON HIGH VOLTAGE -------------------------------- ")
        info = "-->> TURNING ON HIGH VOLTAGE FOR THE Si SENSOR"
        write_log_file(module_dir, module_sn, info)
        # Function
        hv_current_n, hv_current_p = powerON_hv(hv_n_channel, hv_p_channel, bias_voltage)
        info = "I_SENSOR_150V_N: {} [uA]".format(hv_current_n)
        write_data_file(module_dir, module_sn_tmp, info)
        info = "I_SENSOR_150V_P: {} [uA]".format(hv_current_p)
        write_data_file(module_dir, module_sn_tmp, info)
        info = "<<-- FINISHED BIASING THE Si SENSOR"
        write_log_file(module_dir, module_sn, info)
    
    elif (test_step =="read_lv_bc"):
        log.info(" ---------------------- READING LV VALUES BEFORE CONFIGURATION------------------------- ")
        info = "-->> READING LV VALUES BEFORE CONFIGURATION"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Step 4: --------------- Measuring LV before configuration ---------------
        info = ""
        write_data_file(module_dir, module_sn_tmp, info)
        info = "LV_BEF_CONFIG_N"
        write_data_file(module_dir, module_sn_tmp, info)
        lv_nside_bc = reading_lv('N')
        info = "FEB N-side:\t{}\tLV_1.2_BC_N [V]:\t{}\tI_1.2_BC_N [A]:\t{}".format(feb_nside, lv_nside_bc[0], lv_nside_bc[1])
        write_data_file(module_dir, module_sn_tmp, info)
        info = "FEB N-side:\t{}\tLV_1.8_BC_N [V]:\t{}\tI_1.8_BC_N [A]:\t{}".format(feb_nside, lv_nside_bc[2], lv_nside_bc[3])
        write_data_file(module_dir, module_sn_tmp, info)
        info = "LV_BEF_CONFIG_P"
        write_data_file(module_dir, module_sn_tmp, info)
        lv_pside_bc = reading_lv('P')
        info = "FEB P-side:\t{}\tLV_1.2_BC_P [V]:\t{}\tI_1.2_BC_P [A]:\t{}".format(feb_pside, lv_pside_bc[0], lv_pside_bc[1])
        write_data_file(module_dir, module_sn_tmp, info)
        info = "FEB P-side:\t{}\tLV_1.8_BC_P [V]:\t{}\tI_1.8_BC_P [A]:\t{}".format(feb_pside, lv_pside_bc[2], lv_pside_bc[3])
        write_data_file(module_dir, module_sn_tmp, info)
        info = "<<-- FINISHED READING LV VALUES BEFORE CONFIGURATION"
        write_log_file(module_dir, module_sn, info)

    elif (test_step =="std_config"):
        log.info(" ------------- LOADING STANDARD CONFIGURATION --------------------- ")
        info = "-->> LOADING STANDARD CONFIGURATION"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Setp 5: --------------- Setting the standard ASIC configuration ---------
        load_STD_Config(smx_l_nside, 'N', feb_nside)
        load_STD_Config(smx_l_pside, 'P', feb_pside)
        info = "<<-- FINISHED LOADING STANDARD CONFIGURATION"
        write_log_file(module_dir, module_sn, info)

    elif (test_step =="read_asic_id"):
        log.info("----------------------- READING ASICs ID -------------------------- ")
        info = "-->> READING ASICs ID"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Step 8: ----------------------- Reading ASIC ID -------------------------
        read_asicIDs_FEB(smx_l_nside, 'N', feb_nside)
        read_asicIDs_FEB(smx_l_pside, 'P', feb_pside)
        info = "<<-- FINISHED READING ASICs ID"
        write_log_file(module_dir, module_sn, info)

    elif (test_step =="set_trim_default"):
        log.info("---------------- LOADING DEFAULT TRIM VALUES ---------------------- ")
        info = "-->> LOADING DEFAULT TRIM VALUES"
        write_log_file(module_dir, module_sn, info)
        #Function
        # Step 6: --------------- Loading the default trim on the ASICs -----------
        set_Trim_default(smx_l_nside, 'N', feb_nside, cal_asic_list_nside)
        set_Trim_default(smx_l_pside, 'P', feb_pside, cal_asic_list_pside)
        info = "<<-- FINISHED LOADING DEFAULT TRIM VALUES"
        write_log_file(module_dir, module_sn, info)

    elif (test_step =="read_lv_ac"):
        log.info(" ---------------------- READING LV VALUES AFTER CONFIGURATION ------------------------- ")
        info = "-->> READING LV VALUES AFTER CONFIGURATION"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Step 7: --------------- Measuring LV after  configuration ---------------
        info = ""
        write_data_file(module_dir, module_sn_tmp, info)
        info = "LV_AFT_CONFIG_N"
        write_data_file(module_dir, module_sn_tmp, info)
        time.sleep(10)
        lv_nside_ac = reading_lv('N')
        info = "FEB N-side:\t{}\tLV_1.2_AC_N [V]:\t{}\tI_1.2_AC_N [A]:\t{}".format(feb_nside, lv_nside_ac[0], lv_nside_ac[1])
        write_data_file(module_dir, module_sn_tmp, info)
        info = "FEB N-side:\t{}\tLV_1.8_AC_N [V]:\t{}\tI_1.8_AC_N [A]:\t{}".format(feb_nside, lv_nside_ac[2], lv_nside_ac[3])
        write_data_file(module_dir, module_sn_tmp, info)
        info = "LV_AFT_CONFIG_P"
        write_data_file(module_dir, module_sn_tmp, info)
        lv_pside_ac = reading_lv('P')
        info = "FEB P-side:\t{}\tLV_1.2_AC_P [V]:\t{}\tI_1.2_AC_P [A]:\t{}".format(feb_pside, lv_pside_ac[0], lv_pside_ac[1])
        write_data_file(module_dir, module_sn_tmp, info)
        info = "FEB P-side:\t{}\tLV_1.8_AC_P [V]:\t{}\tI_1.8_AC_P [A]:\t{}".format(feb_pside, lv_pside_ac[2], lv_pside_ac[3])
        write_data_file(module_dir, module_sn_tmp, info)
        info = "<<-- FINISHED READING LV VALUES AFTER CONFIGURATION"
        write_log_file(module_dir, module_sn, info)

    elif (test_step =="check_vddm_temp"):
        log.info("-------------- READING VDDM & TEMPERATURE ------------------------- ")
        info = "-->> READING VDDM & TEMPERATURE"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Step 9: -------------------- Reading VVDM & TEMP ------------------------
        read_VDDM_TEMP_FEB(smx_l_nside, 'N', feb_nside)
        read_VDDM_TEMP_FEB(smx_l_pside, 'P', feb_pside)
        info = "<<-- FINISHED READING VDDM & TEMPERATURE"
        write_log_file(module_dir, module_sn, info)
    
    elif (test_step =="get_vrefs"):
        log.info("-------------- FINDING VREF_P, VREF_N & THR@_GLB FOR CALIBRATION ------------------------- ")
        info = "-->> FINDING VREF_P, VREF_N & THR@_GLB FOR CALIBRATION"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Step 10: ----------------- SCAN VREF_P, N & THR2_GLB --------------------
        # 10.1 Scanning ADC and FAST discriminator potentials beforte calibration
        # define default cal settings just in case
        cal_set_nside = scan_VrefP_N_Thr2glb(smx_l_nside, 'N', feb_nside, cal_asic_list_nside, npulses, test_ch, amp_cal_min, amp_cal_max, amp_cal_fast, vref_t)
        cal_set_pside = scan_VrefP_N_Thr2glb(smx_l_pside, 'P', feb_pside, cal_asic_list_pside, npulses, test_ch, amp_cal_min, amp_cal_max, amp_cal_fast, vref_t)
        info = "<<-- FINISHED FINDING VREF_P, VREF_N & THR@_GLB FOR CALIBRATION"
        write_log_file(module_dir, module_sn, info)
    
    elif (test_step =="set_calib_par"):
        log.info("-------------- SETTING THE CALIBRATION PARAMETERS ------------------------- ")
        info = "--> SETTING THE CALIBRATION PARAMETERS"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Step 10.2: --------- Printing calibration settings ----------------------
        print_cal_settings('N', cal_set_nside, cal_asic_list_nside)
        print_cal_settings('P', cal_set_pside, cal_asic_list_pside)
        # 10.3 Writing calibration settings on the ASIC
        writing_cal_settings(smx_l_nside, 'N', feb_nside, cal_set_nside, cal_asic_list_nside)
        writing_cal_settings(smx_l_pside, 'P', feb_pside, cal_set_pside, cal_asic_list_pside)
        info = "<<-- FINISHED SETTING THE CALIBRATION PARAMETERS"
        write_log_file(module_dir, module_sn, info)

    elif (test_step =="get_trim"):
        log.info("-------------- CALIBRATING THE MODULE ------------------------- ")
        info = "-->> CALIBRATING THE MODULE"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Step 11: ------------------ Calibrating the Module ---------------------
        calib_FEB(smx_l_nside, trim_dir, 'N', feb_nside, cal_asic_list_nside, npulses, amp_cal_min, amp_cal_max, amp_cal_fast, much_mode_on = 0)
        calib_FEB(smx_l_pside, trim_dir, 'P', feb_pside, cal_asic_list_pside, npulses, amp_cal_min, amp_cal_max, amp_cal_fast, much_mode_on = 0)
        info = "<<-- FINISHED CALIBRATING THE MODULE"
        write_log_file(module_dir, module_sn, info)
            
    elif (test_step =="set_trim_calib"):
        log.info("---------------- SETTING THE CALIBRATION TRIM ----------------------------- ")
        info = "--> SETTING THE CALIBRATION TRIM"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Step 12: ----------------- Setting Calibration Trim --------------------   
        set_trim_calib(smx_l_nside, trim_dir, 'N', feb_nside, cal_asic_list_nside, much_mode_on = 0)
        set_trim_calib(smx_l_pside, trim_dir, 'P', feb_pside, cal_asic_list_pside, much_mode_on = 0)
        info = "<<-- FINISHED SETTING THE CALIBRATION TRIM"
        write_log_file(module_dir, module_sn, info)
    
    elif (test_step =="check_trim"):
        log.info("-------------- MEASURING ENC AND CHECKING CALIBRATION ------------------------- ")
        info = "--> MEASURING ENC AND CHECKING CALIBRATION"
        write_log_file(module_dir, module_sn, info)
        # Function
        # Step 12: ------------------- ENC & Checking Trim -----------------------   
        check_trim(smx_l_nside, pscan_dir, 'N', feb_nside, cal_asic_list_nside, disc_list, vp_min, vp_max, vp_step, npulses)
        check_trim(smx_l_pside, pscan_dir, 'P', feb_pside, cal_asic_list_pside, disc_list, vp_min, vp_max, vp_step, npulses)
        info = "<<-- FINISHED MEASURING ENC AND CHECKING CALIBRATION"
        write_log_file(module_dir, module_sn, info)

    elif (test_step == "turn_hv_off"):
        log.info(" --------------------------- TURNING OFF HIGH VOLTAGE -------------------------------- ")
        info = "-->> TURNING OFF HIGH VOLTAGE FOR THE Si SENSOR"
        write_log_file(module_dir, module_sn, info)
        # Function
        flag_hv_off = 0
        while( flag_hv_off == 0):
            if (powerOff_hv(hv_n_channel, hv_p_channel) == True):
                info = "<<-- HV Si SENSOR is OFF"
                write_log_file(module_dir, module_sn, info)
                flag_hv_off = 1
            else:
                log.info("MODULE HV STILL ON")
                time.sleep(60)
        info = "<<-- FINISHED TURNING OFF HV FOR Si SENSOR"
        write_log_file(module_dir, module_sn, info)

    elif (test_step =="conn_check"):
        log.info("-------------- CHECKING CHANNELS' CONNECTIVITY ------------------------- ")
        info = "--> CHECKING CHANNELS' CONNECTIVITY"
        write_log_file(module_dir, module_sn, info)
        # Function
        connection_check(smx_l_nside, conn_check_dir, 'N', feb_nside, cal_asic_list_nside, nloops, vref_t_low)
        connection_check(smx_l_pside, conn_check_dir, 'P', feb_pside, cal_asic_list_pside, nloops, vref_t_low)
        info = "<<-- FINISHED CHECKING CHANNELS' CONNECTIVITY"
        write_log_file(module_dir, module_sn, info)

    elif(test_step == "long_run"):
        log.info("-------------- CHECKING ENC LONG RUN STABILITYY ------------------------- ")
        info = "--> CHECKING ENC LONG RUN STABILITY "
        write_log_file(module_dir, module_sn, info)
        # Function                                                                                                                                                                                 
        vref_t_low = 118
        nloops = 30
        cal_asic_list_nside = [2]
        cal_asic_list_pside = [2]
        for n in range(0,nloops):
            read_VDDM_TEMP_FEB(smx_l_nside, 'N', feb_nside )
            check_trim(smx_l_nside, pscan_dir, 'N', feb_nside, cal_asic_list_nside, disc_list, vp_min, vp_max, vp_step, npulses)
            read_VDDM_TEMP_FEB(smx_l_pside, 'P', feb_pside)
            check_trim(smx_l_pside, pscan_dir, 'P', feb_pside, cal_asic_list_pside, disc_list, vp_min, vp_max, vp_step, npulses)
            info = "--> RUN {} of CHECKING ENC LONG RUN STABILITY ".format(n)
            log.info(info)
            write_log_file(module_dir, module_sn, info)
        info = "<<-- FINISHED ENC LONG RUN STABILITYY"
        write_log_file(module_dir, module_sn, info)
    elif(test_step == "set_lv_off"):
        log.info("-------------- SETTING LV OFF. FINALIZING TESTS ------------------------- ")
        info = "--> SETTING OFF LV N"
        write_log_file(module_dir, module_sn, info)
        powerOff_lv('P')
        time.sleep(15)
        flag_off = 0
        while(flag_off!= 1):
            powerOff_lv('N')
            time.sleep(15)
            lv_nside_end = reading_lv('N')
            if (lv_nside_end[0] == 0 and lv_nside_end[2] == 0):
                info = "<<-- N-side FEB is OFF. LV_12_END_N: \t{}\t LV_18_END_N: \t{}".format(lv_nside_end[0], lv_nside_end[2])
                log.info(info)
                write_log_file(module_dir, module_sn, info)
                flag_off = 1

            
                
            
    else:
        if (test_step[0] == '#'):
            log.info("----------------------- SKIPPING TEST: %s -------------------------- ", test_step)
