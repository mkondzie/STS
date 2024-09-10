import os
import sys
import json
import argparse
import time
import logging
import epics
from epics import PV
from datetime import datetime
import subprocess
import pprint

for root, dirs, files in os.walk("/home/cbm/cbmsoft/emu_test_module/"):
    for dir in dirs:
        path = os.path.join(root, dir)
        sys.path.append(path)

import uhal
import agwb
from smx_tester import *
import msts_defs as smc
import smx_oper as smxoper
#import test_setup as setup
#from opm import Opm


parser = argparse.ArgumentParser(description='A simple Python script with argument parsing.')

parser.add_argument('--config_file',  type=str, default="module_test_config.json", help='Toggle switch for Connection check task')
parser.add_argument('--conn_check',  type=bool, default=False, help='Toggle switch for Connection check task')
parser.add_argument('--feb_n_side',  type=str, default="", help='Toggle switch for Connection check task')
parser.add_argument('--feb_p_side',  type=str, default="", help='Toggle switch for Connection check task')
parser.add_argument('--module_id',  type=str, default="", help='Toggle switch for Connection check task')
parser.add_argument('--output_path', type=str,  default="",    help='Output folder path')

# Parsing the arguments
args = parser.parse_args()

# Accessing the values
config_file = args.config_file
conn_check = args.conn_check
feb_n_side = args.feb_n_side
feb_p_side = args.feb_p_side
module_id = args.module_id
output_path = args.output_path


#Lists of subsequences for module tests
# ------------------------------------------------------------------------------------------------------------------
# Possible test sequences.

test_list_init = ["full_sync"]
test_list_comm = ["read_lv_ac","std_config","read_asic_id", "set_trim_default", "read_lv_ac", "check_vddm_temp"]
test_list_check = ["conn_check"]

test_list = []

if (config_file != "None" and os.path.isfile(config_file)):
    with open(config_file, "r") as json_file:
        loaded_data = json.load(json_file)
    for test in loaded_data:
        print (test, loaded_data[test])
        if (loaded_data[test]):
            test_list.append(test)

# Your code logic here
print(f"module_id path: {module_id}")
print(f"feb_n_side path: {feb_n_side}")
print(f"feb_p_side path: {feb_p_side}")
print(f"Output path: {output_path}")
print(f"run_conn_check: {conn_check}")

summary_file = open(("{}/{}_last_test.info").format(output_path, module_id), "w+")
print("saving file to ", "{}/{}_last_test.info".format(output_path, module_id))
header_line  = "# FEB-ID_\t_POLARITY_\t_ASIC-HW-ADDRESS_\t_VDDM POTENTIAL (LSB) (mV)_\t_TEMPERATURE (LSB) (C)_\n"
summary_file.write(str(header_line))

# --------------------------- SETTING LOGGING DETAILS ------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s,%(msecs)03d:%(module)s:%(levelname)s:%(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger()
fh = logging.FileHandler(sys.argv[0].replace('py', 'log'), 'a+')
fh.setLevel(logging.DEBUG)
fmt = logging.Formatter('[%(levelname)s] %(message)s')
fh.setFormatter(fmt)
log.addHandler(fh)
uhal.setLogLevelTo(uhal.LogLevel.WARNING)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ------------------------ CONNECTING EMUs, FINDING FEBs AND SYNC ASICs -----------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def general_sync(emu_list):
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
    #active_downlinks = [1,2] #For modules PB
    active_downlinks = [0,3]  #For modules PA

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
    pol = pol_usr

    if (pol == 'N' or pol == '0'):
        print("Reading LV for N-side (electrons polarity):")
        #Reading Voltage and current values and writting them in a file
        lv_12 = PV("CBM:STS:WIENLV_U100:vmon").value
        lv_18 = PV("CBM:STS:WIENLV_U101:vmon").value
        i_12 = PV("CBM:STS:WIENLV_U100:imon").value
        i_18 = PV("CBM:STS:WIENLV_U101:imon").value

    elif(pol == 'P' or pol == '1'):
        print("Reading LV for P-side (holes polarity):")
        #Reading Voltage and current values and writting them in a file
        lv_12 = PV("CBM:STS:WIENLV_U102:vmon").value
        lv_18 = PV("CBM:STS:WIENLV_U103:vmon").value
        i_12 = PV("CBM:STS:WIENLV_U102:imon").value
        i_18 = PV("CBM:STS:WIENLV_U103:imon").value

    else:
        log.error("Please, indicate a polarity for reading the corresponding LV potentials, in the following way:")
        log.error("N or 0 for n-side, 0 for electrons polarity")
        log.error("P or 1 for p-side, 1 for holes polarity")
        sys.exit()

    print("LV values: ")
    if (pol == 'N' or pol == '0'):
        print("LV_12_N [V]: {}  I_12_N [A]: {}".format(lv_12, i_12))
        print("LV_18_N [V]: {}  I_18_N [A]: {}".format(lv_18, i_18))
    else:
        print("LV_12_P [V]: {}  I_12_P [A]: {}".format(lv_12, i_12))
        print("LV_18_P [V]: {}  I_18_P [A]: {}".format(lv_18, i_18))

    return (lv_12, i_12, lv_18, i_18)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ----------------------------- MODULE TESTING: DIRECTORY AND FILES ---------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
feb_list = [feb_p_side, feb_n_side]
connection_dir = output_path

# Module testing directory & files
print (feb_list)

def initWorkingDirectory():
    # creating log file

    logfile = open(output_path + "/" + str(module_id) + "/" + str(module_id) +".log", "a+")
    datafile = open(output_path + "/" + str(module_id) + "/" + str(module_id) +".dat","a+")
    date = time.strftime("%y%m%d-%H:%M:%S")

    # creating data file and initializing it
    datafile.write("MODULE_ID: ")
    datafile.write(str(module_id))
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
    logfile.write(str(module_id))
    logfile.write("\n")
    logfile.close()

    return module_id

def close_log_file(module_id):
    logfile = open(output_path + str(module_id) + "/" + str(module_id) +".log", "a+")
    datafile = open(output_path + str(module_id) + "/" + str(module_id) +".dat","a+")
    # ending testing & closing data and log files
    logfile.write(datetime)
    logfile.write("\t")
    logfile.write("Ending test sequence")
    logfile.write("\n")
    logfile.close()
    datafile.close()

def write_data_file(output_path, module_id, info = ""):
    datafile = open(output_path + str(module_id) + "/" + str(module_id) +".dat","a+")
    datafile.write(str(info))
    datafile.write("\n")
    datafile.close()

def lv_file(output_path, module_id, info = ""):
    date = time.strftime("%y%m%d-%H:%M:%S")
    lv_info = open(output_path + str(module_id) + "/" + "lv_" + str(module_id) +".dat","a+")
    lv_info.write(date)
    lv_info.write(";")
    lv_info.write(str(info))
    lv_info.write("\n")
    lv_info.close()

def write_log_file(output_path, module_id, info = ""):
    date = time.strftime("%y%m%d-%H:%M:%S")
    logfile = open(output_path + str(module_id) + "/" + str(module_id) +".log", "a+")
    logfile.write(date)
    logfile.write("\t")
    logfile.write("[INFO]")
    logfile.write("\t")
    logfile.write(str(info))
    logfile.write("\n")
    logfile.close()

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ----------------------------- MODULE TESTING: OPERATING FUNCTIONS ---------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def read_asicIDs_FEB(feb_sn, smx_l_side, pol):
    # Function to read the overall ASIC ID (FEB where ASIC belongs, polarity, HW address, ASIC e-fuse ID(string) and (int))
    write_data_file(output_path, module_sn_tmp, " ")
    header_line  = "--> FEB-ID_\t_POLARITY_\t_ASIC-HW-ADDRESS_\t\t_ASIC-EFUSE-ID-(STR)_\t\t_ASIC-EFUSE-ID-(INT)_"
    log.info(header_line)
    write_data_file(output_path, module_sn_tmp, header_line)

    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
        pol_calib = 0
    else:
        pol_str = 'P-side'
        pol_calib = 1

    for smx in smx_l_side:
        addr = smx.address
        asic_id_int = smx.read_efuse()
        asic_id_str = smx.read_efuse_str()
        info  = "{} \t\t {} \t\t {} \t\t {} \t\t {}".format(feb_sn, pol_str, addr, asic_id_str, asic_id_int)
        write_data_file(output_path, module_sn_tmp, info)
        log.info(info)

def read_VDDM_TEMP_FEB(feb_sn, smx_l_side, pol):
    # Function to read the VDDM and TEMP of the ASICs
    header_line  = "--> FEB-ID_\t_POLARITY_\t_ASIC-HW-ADDRESS_\t_VDDM POTENTIAL (LSB) (mV)_\t\t_TEMPERATURE (LSB) (C)_\n"
    write_data_file(output_path, module_sn_tmp, " ")

    log.info(header_line)

    write_data_file(output_path, module_sn_tmp,header_line)


    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
    else:
        pol_str = 'P-side'
    for smx in smx_l_side:
        asic_vddm = smx.read_vddm()
        #asic_vddm_src = smx.read_diag("Vddm")
        asic_temp = smx.read_temp()
        #asic_temp_src = smx.read_diag("Temp")
        info = "{}\t{}\t{}\t{}\t{:.1f}\t{:.1f}\t{}\n".format(feb_sn, pol_str, smx.address, asic_vddm[0], asic_vddm[1], asic_temp[0], asic_temp[1])
        write_data_file(output_path, module_sn_tmp,info)
        summary_file.write(str(info))
        log.info(info)

def load_STD_Config(smx_l_side, pol):
    # Function to write on each ASIC the default settings
    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
        pol_calib = 0
    else:
        pol_str = 'P-side'
        pol_calib = 1
    for smx in smx_l_side:
        header_line  = "--> SETTING STANDARD CONFIGURATION for ASIC with HW address {} and polarity {}".format(smx.address,pol_str)
        log.info(header_line)
        smx.write_def_ana_reg(smx.address, pol_calib )
        smx.read_reg_all(compFlag = False)

def set_Trim_default(smx_l_side, pol, cal_asic_list):
    if (pol == 'N' or pol == '0'):
        pol_str = 'N-side'
    else:
        pol_str = 'P-side'
    for smx in smx_l_side:
        asic_hw_addr = smx.address
        if (asic_hw_addr in cal_asic_list):
            header_line  = "--> SETTING DEFAULT TRIM for ASIC with HW address {} and polarity {}".format(smx.address, pol_str)
            log.info(header_line)
            smx.set_trim_default(128,36)
        else:
            info = "--> NO SETTING DEFAULT TRIM for ASIC with HW ADDRESS {} in {}".format(asic_hw_addr, pol_str)
            log.info(info)
            pass

def connection_check(smx_l_side, conn_dir, pol, cal_asic_list, nloops = 10, vref_t = 118 ):
    if (pol == 'N' or pol == '0'):
        pol_str = 'elect'
        pol_calib = 0
    elif (pol == 'P' or pol == '1'):
        pol_str = 'holes'
        pol_calib = 1
    else:
        log.error("Please check the polarity is correct: (ex: string 'N' or '0', 'P' or '1')")
    conn_dir = connection_dir + str(module_id)
    for smx in smx_l_side:
        asic_hw_addr = smx.address
        if (asic_hw_addr in cal_asic_list):
            # Checking ENC and calibration results for an ASIC
            info = "--> CONNECTION-CHECK SCAN for ASIC with HW ADDRESS {} in {}".format(asic_hw_addr, pol_str)
            log.info(info)
            # Elements for the pscan file
            #asic_id_str = smx.read_efuse_str()
            smx.connection_check(conn_dir, pol_calib, nloops, vref_t)
    return 0

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ----------------------------- MODULE TESTING: VARIABLES DEFINITION --------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
emu_list = ["EMU_239"]                           #["EMU_213"]     # List of EMU boards used during the test

# List of ASICs HW address to be tested in case the
#cal_asic_list = [0, 1, 2, 3, 4, 5, 6, 7]
#cal_asic_list_nside = []
cal_asic_list_nside = [0, 1, 2, 3, 4, 5, 6, 7]
cal_asic_list_pside = [0, 1, 2, 3, 4, 5, 6, 7]


# Operating variables GetVrefs parameters
vref_t = 118                    # Vrf_T value  = 54 in the largest VRef_T range. To determine calib par and to calibrate ADC


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ---------------------------------- MODULE TESTING: RUNNING  ---------------------------------
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Step 0: -------------- Creating the working directory -------------------
module_sn_tmp = initWorkingDirectory()
info = "MODULE ID: \t{}".format(module_id)
write_data_file(output_path, module_sn_tmp, info)
write_log_file(output_path, module_id, info)

# Setp 1: ------------------ FEBs ID -------------------------
feb_nside =feb_n_side
feb_pside = feb_p_side
# RUNNING TEST SEQUENCE
# -------------------------------------------
for test_step in test_list:
    if (test_step == "full_sync"):
        log.info(" ------------------- RUNNING SYNCHRONIZATION ---------------------- ")
        info = "-->> RUNNING SYNCHRONIZATION"
        write_log_file(output_path, module_id, info)
        # Function
        # Step 1: -------- Starting with connection and synchronization -----------
        info = "--> EMU board: {}".format(emu_list)
        write_data_file(output_path, module_sn_tmp, info)
        smx_l = general_sync(emu_list)
        # 2.1 Determining the number of ASICs per side
        n_asic_all = scanning_asics(smx_l)
        # 2.2 Assigning the ASICs according to polarities
        n_asics = n_asic_all[0]
        p_asics = n_asic_all[1]
        write_data_file(output_path, module_sn_tmp, "")
        info = "--> Number of sync ASICs for N-side:\t {}".format(n_asics)
        write_data_file(output_path, module_sn_tmp, info)
        info = "--> Number of sync ASICs for P-side:\t {}".format(p_asics)
        write_data_file(output_path, module_sn_tmp, info)
        # 2.3 Distributing the ASICs in arrays according to operational polarity
        smx_l_nside = smx_l[0:n_asics]
        smx_l_pside = smx_l[n_asics:n_asics + p_asics]
        info = "<<-- FINISHED SYNCHRONIZATION"
        write_log_file(output_path, module_id, info)

    elif (test_step =="std_config"):
        log.info(" ------------- LOADING STANDARD CONFIGURATION --------------------- ")
        info = "-->> LOADING STANDARD CONFIGURATION"
        write_log_file(output_path, module_id, info)
        # Function
        # Setp 2: --------------- Setting the standard ASIC configuration ---------
        load_STD_Config(smx_l_nside, 'N')
        load_STD_Config(smx_l_pside, 'P')
        info = "<<-- FINISHED LOADING STANDARD CONFIGURATION"
        write_log_file(output_path, module_id, info)
        #------------------------ QA setting --------------------------------------
        #qa_write_read(smx_l_nside, 'N')
        #qa_write_read(smx_l_pside, 'P')

    elif (test_step =="read_asic_id"):
        log.info("----------------------- READING ASICs ID -------------------------- ")
        info = "-->> READING ASICs ID"
        write_log_file(output_path, module_id, info)
        # Function
        # Step 3: ----------------------- Reading ASIC ID -------------------------
        read_asicIDs_FEB(feb_nside, smx_l_nside, 'N')
        read_asicIDs_FEB(feb_pside, smx_l_pside, 'P')
        info = "<<-- FINISHED READING ASICs ID"
        write_log_file(output_path, module_id, info)

    elif (test_step =="read_lv_ac"):
        log.info(" ---------------------- READING LV VALUES ------------------------- ")
        info = "-->> READING LV VALUES"
        write_log_file(output_path, module_id, info)
        # Function
        # Step 4: --------------- Measuring LV after  configuration ---------------
        info = "--> VOLTAGE AND CURRENT"
        write_data_file(output_path, module_sn_tmp, info)
        time.sleep(10)
        lv_nside_ac = reading_lv('N')
        info = "N-side;{};LV_1.2_N[V];{};I_1.2_N[A];{}".format(feb_nside, lv_nside_ac[0], lv_nside_ac[1])
        write_data_file(output_path, module_sn_tmp, info)
        lv_file(output_path, module_sn_tmp, info)
        info = "N-side;{};LV_1.8_N[V];{};I_1.8_N[A];{}".format(feb_nside, lv_nside_ac[2], lv_nside_ac[3])
        write_data_file(output_path, module_sn_tmp, info)
        lv_file(output_path, module_sn_tmp, info)
        lv_pside_ac = reading_lv('P')
        info = "P-side;{};LV_1.2_P[V];{};I_1.2_P[A];{}".format(feb_pside, lv_pside_ac[0], lv_pside_ac[1])
        write_data_file(output_path, module_sn_tmp, info)
        lv_file(output_path, module_sn_tmp, info)
        info = "P-side;{};LV_1.8_P[V];{};I_1.8_P[A];{}".format(feb_pside, lv_pside_ac[2], lv_pside_ac[3])
        write_data_file(output_path, module_sn_tmp, info)
        lv_file(output_path, module_sn_tmp, info)
        info = "<<-- FINISHED READING LV VALUES"
        write_log_file(output_path, module_id, info)

    elif (test_step =="check_vddm_temp"):
        log.info("-------------- READING VDDM & TEMPERATURE ------------------------- ")
        info = "-->> READING VDDM & TEMPERATURE"
        write_log_file(output_path, module_id, info)
        # Function
        # Step 5: -------------------- Reading VVDM & TEMP ------------------------
        read_VDDM_TEMP_FEB(feb_nside, smx_l_nside, 'N')
        read_VDDM_TEMP_FEB(feb_pside, smx_l_pside, 'P')
        info = "<<-- FINISHED READING VDDM & TEMPERATURE"
        write_log_file(output_path, module_id, info)

    elif (test_step =="set_trim_default"):
        log.info("---------------- LOADING DEFAULT TRIM VALUES ---------------------- ")
        info = "-->> LOADING DEFAULT TRIM VALUES"
        write_log_file(output_path, module_id, info)
        #Function
        # Step 6: --------------- Loading the default trim on the ASICs -----------
        set_Trim_default(smx_l_nside, 'N', cal_asic_list_nside)
        set_Trim_default(smx_l_pside, 'P', cal_asic_list_pside)
        info = "<<-- FINISHED LOADING DEFAULT TRIM VALUES"
        write_log_file(output_path, module_id, info)

    elif (test_step =="conn_check"):
        log.info("-------------- CHECKING CHANNELS' CONNECTIVITY ------------------------- ")
        info = "--> CHECKING CHANNELS' CONNECTIVITY"
        write_log_file(connection_dir, module_id, info)

        vref_t_low = 102
        nloops = 10
        connection_check(smx_l_nside, output_path, 'N', cal_asic_list_nside, nloops, vref_t_low)
        connection_check(smx_l_pside, output_path, 'P', cal_asic_list_pside, nloops, vref_t_low)
        info = "<<-- FINISHED CHECKING CHANNELS' CONNECTIVITY"
        write_log_file(connection_dir, module_id, info)

    else:
        if (test_step[0] == '#'):
            log.info("----------------------- SKIPPING TEST: %s -------------------------- ", test_step)
