#!/usr/bin/python
import re
import sys
import os
import tkinter as tk
from graphics import *

HEIGHT = 900
WIDTH  = 1200


####If numbers cut out, need to be doubled
######################################
def Get_input(x):
    input = ""
    with open(x,"r") as f:
        for line in f:
            input += line

    return input
######################################

def Get_dictionaries(x):
    """
    takes in IgG in SMILES format and identifies variables for dynamically rendering image
    """
    y       = re.sub("\s","",x)
    if "|" in y:
        splitx  = y.split("|")
    else:
        splitx = [y]
    chains  = []
    ADCs    = []
    VHa     = {}
    VHb     = {}
    VLa     = {}
    VLb     = {}
    fragment= {}
    Salt_bridges   = ""
    VHa_VLa_bonds  = ""
    VHb_VLb_bonds  = ""
    CH1a_CL1a_bonds=""
    CH1b_CL1b_bonds=""
    #get chains and watch out for ADCs
    for i in splitx:
        if len(splitx) > 1:
            if i[0] != "X":
                chain = i.split("(")[0]
                chains.append(chain)
            elif i[0] == "X":
                chain = i.split("-")[1]
                chain = chain.split("(")[0]
                chains.append(chain)
        elif len(splitx) == 1:
            chains.append("fragment")

    for i in range(len(chains)):
        dict = str(re.sub("\.","",str(chains[i])))
        chain = splitx[i].split("-")
        for j in range(len(chain)):
            domain   =  re.findall("^CH[1-9]|^VL.[ab]|^VH.[ab]|^VL[1-9]|^VH[1-9]|^CL|^VH|^VL|^H|^X", str(chain[j]))
            domain   =  str(re.sub("\[|\'|\]","", str(domain)))

            #if chain[i] != fragment:
            #    if domain == ""
            ##Get Bond numbers
            if domain == "H" and Salt_bridges == "":
                Salt_bridges = re.findall("\{.*?\}", str(chain[j]))
                Salt_bridges = str(Salt_bridges)
                Salt_bridges = int(re.sub("\{|\'|\}|\[|\]","", Salt_bridges))
            elif domain == "VH":
                if dict == "VHa" and VHa_VLa_bonds == "":
                    VHa_VLa_bonds_match = re.findall("\{.*?\}", str(chain[j]))
                    if VHa_VLa_bonds_match != []:
                        VHa_VLa_bonds = str(VHa_VLa_bonds_match)
                        VHa_VLa_bonds = int(re.sub("\{|\'|\}|\[|\]","", VHa_VLa_bonds))
                elif dict == "VHb" and VHb_VLb_bonds == "":
                    VHb_VLb_bonds_match = re.findall("\{.*?\}", str(chain[j]))
                    if VHb_VLb_bonds_match != []:
                        VHb_VLb_bonds = str(VHb_VLb_bonds_match)
                        VHb_VLb_bonds = int(re.sub("\{|\'|\}|\[|\]","", VHb_VLb_bonds))
            elif domain == "CH1":
                if dict == "VHa" and CH1a_CL1a_bonds == "":
                    CH1a_CL1a_bonds_match = re.findall("\{.*?\}", str(chain[j]))
                    if CH1a_CL1a_bonds_match != []:
                        CH1a_CL1a_bonds = str(CH1a_CL1a_bonds_match)
                        CH1a_CL1a_bonds = int(re.sub("\{|\'|\}|\[|\]","", CH1a_CL1a_bonds))
                elif dict == "VHb" and CH1b_CL1b_bonds == "":
                    CH1b_CL1b_bonds_match = re.findall("\{.*?\}", str(chain[j]))
                    if CH1b_CL1b_bonds_match != []:
                        CH1b_CL1b_bonds = str(CH1b_CL1b_bonds_match)
                        CH1b_CL1b_bonds = int(re.sub("\{|\'|\}|\[|\]","", CH1b_CL1b_bonds))

            location    = []
            locationstr =  re.findall("\((.*?)\)", str(chain[j]))
            locationstr =  str(re.sub("\[|\'|\]","", str(locationstr)))
            locationstr =  str(re.sub(":",",", str(locationstr)))
            locationstr = locationstr.split(",")

            if domain != "X" and domain != "":
                for x in range(len(locationstr)):
                    location.append(int(locationstr[x]))
            elif domain == "X":
                location = re.findall("\[(.*?)\]", str(chain[j]))

            if dict == "VHa" and domain !="":
                VHa[domain] = location
            elif dict == "VHb" and domain !="":
                VHb[domain] = location
            elif dict == "VLa" and domain !="":
                VLa[domain] = location
            elif dict == "VLb" and domain !="":
                VLb[domain] = location
            elif dict == "fragment" and domain !="":
                fragment[domain] = location
            else:
                continue

    return(VHa,VLa,VHb,VLb,Salt_bridges,VHa_VLa_bonds,VHb_VLb_bonds,CH1a_CL1a_bonds,CH1b_CL1b_bonds,fragment)


######################################
def Check_interactions(chains_list):
#########Set Variables################

    VHa       = chains_list[0]
    VLa       = chains_list[1]
    VHb       = chains_list[2]
    VLb       = chains_list[3]
    print(chains_list[0])
    print(chains_list[1])
    print(chains_list[2])
    print(chains_list[3])
    Salt_bridge_count     = chains_list[4]
    VHa_VLa_bond_count    = chains_list[5]
    VHb_VLb_bond_count    = chains_list[6]
    CH1a_CL1a_bond_count  = chains_list[7]
    CH1b_CL1b_bond_count  = chains_list[8]
    fragment              = chains_list[9]
    Heavy_chain_a_domains = []
    Light_chain_a_domains = []
    Heavy_chain_b_domains = []
    Light_chain_b_domains = []
    Salt_bridges          = []
    Label_Locations       = []
    Label_Text            = []

    VHa_VLa   = False
    CH1a_CLa  = False
    VHb_VLb   = False
    CH1b_CLb  = False
    CH2a_CH2b = False
    CH3a_CH3b = False

#########VLa Chain coordinates################


    VL1a                       = [100,100,120,100,132,120,150,120,180,174,140,174]
    VL1a_no_V                  = [100,100,140,110,180,174,140,174]
    VL1a_label_location        = [80,140]
    VLa1_top_bond_location     = [120,100]
    VL1a_bottom_bond_location  = [160,174]


    CL1a                       = [140,180,180,180,220,254,180,254]
    CL1a_label_location        = [120,220]
    CL1a_top_bond_location     = [160,180]
    CL1a_bottom_bond_location  = [200,254]

    #scFV_VeryL_IgGa = [0,10,20,10,40,50,30,50]
    #tandem_scFca    = [10,25,20,25,35,45,20,45]
    #Fab_scFV_Fca    = [20,45,30,45,40,65,30,65]
    #scFV_VeryLscFVa = [30,65,40,65,40,85,40,85]
    #scFV_L_IgGa     = [15,5,25,5,35,25,25,25]
    #scFV_LscFVa     = [45,65,55,65,55,85,45,85]

#########VHa Chain coordinates################

    VH1a                      = [170,120,190,120,180,100,200,100,240,174,200,174]
    VH1a_no_V                 = [160,100,200,100,240,174,200,174]
    VH1a_label_location       = [240,140]
    VH1a_top_bond_location    = [180,100]
    VH1a_bottom_bond_location = [220,174]

    CH1a                      = [200,180,240,180,280,254,240,254]
    CH1a_label_location       = [280,220]
    CH1a_top_bond_location    = [220,180]
    CH1a_bottom_bond_location = [260,254]

    CH2a                      = [280,280,320,280,320,354,280,354]
    CH2a_label_location       = [260,320]
    CH2a_top_bond_location    = [300,280]
    CH2a_bottom_bond_location = [300,354]

    CH3a                      = [280,360,320,360,320,434,280,434]
    CH3a_notch                = [280,360,320,360,340,480,320,434,280,434]
    CH3a_label_location       = [260,400]
    CH3a_top_bond_location    = [300,360]
    CH3a_bottom_bond_location = [300,434]


    #scFV_H_IgGa_outer = [60,10,80,10,100,50,80,50]
    #scFV_H_IgGa_inner = [90,10,110,10,130,50,110,50]
    #scFv4_IgG         = [110,50,130,50,150,90,130,90]

    IgG_2scFV1a                      = [280,440,320,440,320,514,280,514]
    IgG_2scFV1a_label_location       = [300,480]
    IgG_2scFV1a_top_bond_location    = [300,440]
    IgG_2scFV1a_bottom_bond_location = [300,514]

    IgG_2scFV2a                      = [280,320,360,520,360,594,280,594]
    IgG_2scFV2a_label_location       = [300,460]
    IgG_2scFV2a_top_bond_location    = [320,540]
    IgG_2scFV2a_bottom_bond_location = [320,694]

    IgG_2scFV3a                      = [220,440,260,440,260,514,220,514]
    IgG_2scFV3a_label_location       = [240,480]
    IgG_2scFV3a_top_bond_location    = [240,440]
    IgG_2scFV3a_bottom_bond_location = [240,514]

    IgG_2scFV4a                      = [220,520,260,520,260,594,220,594]
    IgG_2scFV4a_label_location       = [240,560]
    IgG_2scFV4a_top_bond_location    = [240,520]
    IgG_2scFV4a_bottom_bond_location = [240,594]



#########VLb Chain coordinates################

    VL1b                      = [560,100,540,100,530,120,510,120,480,174,520,174]
    VL1b_no_V                 = [560,100,520,50,480,87,520,87]
    VL1b_label_location       = [560,140]
    VL1b_top_bond_location    = [540,100]
    VL1b_bottom_bond_location = [500,174]


    CL1b                      = [520,180,480,180,440,254,480,254]
    CL1b_label_location       = [520,220]
    CL1b_top_bond_location    = [500,180]
    CL1b_bottom_bond_location = [460,254]

    #scFV_VeryL_IgGa = [165,5,155,5,145,25,165,5]
    #tandem_scFca= [155,25, 145,25,135,45,155,45]
    #Fab_scFV_Fca= [145,45,135,45,125,65,135,65]
    #scFV_VeryLscFVa = [135,65,125,65,125,85,135,85]
    #scFV_L_IgGb = [150,5,140,5,130,25,140,25]
    #scFV_LscFVb = [120,60,110,60,110,80,120,80]

#########VHb Chain coordinates################

    VH1b                      = [460,100,480,100,470,120,490,120,460,174,420,174]
    VH1b_no_v                 = [500,100,460,100,420,174,460,174]
    VH1b_label_location       = [420,140]
    VH1b_top_bond_location    = [480,100]
    VH1b_bottom_bond_location = [440,174]

    CH1b                      = [460,180,420,180,380,254,420,254]
    CH1b_label_location       = [380,220]
    CH1b_top_bond_location    = [440,180]
    CH1b_bottom_bond_location = [400,254]


    CH2b                      = [380,280,340,280,340,354,380,354]
    CH2b_label_location       = [400,320]
    CH2b_top_bond_location    = [360,280]
    CH2b_bottom_bond_location = [360,354]


    CH3b                      = [380,360,340,360,340,434,380,434]
    CH3b_label_location       = [400,400]
    CH3b_top_bond_location    = [360,360]
    CH3b_bottom_bond_location = [360,434]


    #scFv_H_IgGb_outer                = [270,10,250,10,130,50,250,50]
    #scFv_H_IgGb_inner                = [240,10,220,10,200,50,220,50]
    #scFv4_IgGb                       = [220,50,200,50,180,90,200,90]

    IgG_2scFV1b                      = [380,440,340,440,340,514,380,514]
    IgG_2scFV1b_label_location       = [360,480]
    IgG_2scFV1b_top_bond_location    = [360,440]
    IgG_2scFV1b_bottom_bond_location = [360,514]

    IgG_2scFV2b                      = [380,320,340,520,340,594,280,594]
    IgG_2scFV2b_label_location       = [360,600]
    IgG_2scFV2b_top_bond_location    = [360,420]
    IgG_2scFV2b_bottom_bond_location = [360,594]

    IgG_2scFV3b                      = [440,440,400,440,400,514,440,514]
    IgG_2scFV3b_label_location       = [420,480]
    IgG_2scFV3b_top_bond_location    = [420,440]
    IgG_2scFV3b_bottom_bond_location = [420,514]

    IgG_2scFV4b                      = [440,260,400,520,400,594,440,594]
    IgG_2scFV4b_label_location       = [420,600]
    IgG_2scFV4b_top_bond_location    = [420,420]
    IgG_2scFV4b_bottom_bond_location = [420,594]



#########Check Salt Bridge ################

    CH1_CH2bonda = CH1a_bottom_bond_location+CH2a_top_bond_location
    CH1_CH2bondb = CH1b_bottom_bond_location+CH2b_top_bond_location
    Saltbridge_a  =[270,262,388,262]
    Saltbridge_labela = [260,274]
    Saltbridge_b  = [288,272,374,272]
    Saltbridge_labelb = [400,274]

    if Salt_bridge_count   == 0:
        Salt_bridges = [CH1_CH2bonda, CH1_CH2bondb]
    elif Salt_bridge_count == 1:
        Salt_bridges = [CH1_CH2bonda, CH1_CH2bondb, Saltbridge_a]
    elif Salt_bridge_count == 2:
        Salt_bridges = [CH1_CH2bonda, CH1_CH2bondb, Saltbridge_a, Saltbridge_b]

#########Check Bonds between VH and VL Chains################

    VHa_VLa_bond_location   = VH1a_bottom_bond_location+VL1a_bottom_bond_location
    VHb_VLb_bond_location   = VH1b_bottom_bond_location+VL1b_bottom_bond_location
    CH1a_CL1a_bond_location = CH1a_bottom_bond_location+CL1a_bottom_bond_location
    CH1b_CL1b_bond_location = CH1b_bottom_bond_location+CL1b_bottom_bond_location
    if VHa_VLa_bond_count != "":
        Salt_bridges.append(VHa_VLa_bond_location)
    if VHb_VLb_bond_count!= "":
        Salt_bridges.append(VHb_VLb_bond_location)
    if CH1a_CL1a_bond_count != "":
        Salt_bridges.append(CH1a_CL1a_bond_location)
    if CH1b_CL1b_bond_count !="":
        Salt_bridges.append(CH1b_CL1b_bond_location)

#########Check VH/VL ADC Coordinates################
    VLa_ADC                    = [100,60,80,40,100,200,120,40,100,60]
    VLa_ADC_label_location     = [100,10]
    VLa_ADC_bottom_bond_location      = [100,60]
    VHa_ADC                    = [160,60,140,40,160,20,180,40,160,60]
    VHa_ADC_label_location     = [140,10]
    VHa_ADC_bottom_bond_location      = [160,60]
    VLb_ADC                    = [480,60,600,40,580,20,560,40,480,60]
    VLb_ADC_label_location     = [480,10]
    VLb_ADC_bottom_bond_location      = [480,60]
    VHb_ADC                    = [500,60,480,40,500,20,520,40,500,60]
    VHb_ADC_label_location     = [500,10]
    VHb_ADC_bottom_bond_location      = [500,60]
    VHa_VHb_ADC                = [300,180,280,140,300,100,360,100,380,140,360,180,300,180]
    VHa_VHb_ADC_label_location = [320,80]
    VHa_ADC_bond_location      = [300,100]
    VHb_ADC_bond_location      = [360,180]


#########Make bonds################
    def bondmaker(x,y):
        x_bond      = str(x)+"_bottom_bond_location"
        y_bond      = str(y)+   "_top_bond_location"
        return(x_bond,y_bond)


#######Heavy chain A check what domains are there################
    for i in range(len(VHa)):
        keyslist = list(VHa.keys())
        if i==0:
            if keyslist[i] == "VH.a":
                Heavy_chain_a_domains.append(VH1a)
                location = VHa.get(keyslist[i])[0]
                Label_Locations.append((VH1a_label_location))
                Label_Text.append(str(location))
            elif keyslist[i] == "X":
                Heavy_chain_a_domains.append(VHa_ADC)
                TYPE = str(re.sub("\[\'TYPE:|\]","",VHa.get(keyslist[i])[0]))
                Label_Locations.append((VHa_ADC_label_location))
                Label_Text.append(str(TYPE))




        elif i == 1:
                if keyslist[i] == "VH.a" and keyslist[i-1] == "X":
                    Heavy_chain_a_domains.append(VH1a)
                    location = VHa.get(keyslist[i])[0]
                    Label_Locations.append((VH1a_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker("VHa_ADC","VH1a")
                    Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))
                elif keyslist[i] == "VH.b" and keyslist[i-1] == "VL.a":
                    Heavy_chain_a_domains.append(VH1b)
                    location = VHa.get(keyslist[i])[0]
                    Label_Locations.append((VH1b_label_location))
                    Label_Text.append(str(location))
                elif keyslist[i] == "VL.b" and keyslist[i-1] == "VL.a":
                    Heavy_chain_a_domains.append(VL1b)
                    location = VHa.get(keyslist[i])[0]
                    Label_Locations.append((VL1b_label_location))
                    Label_Text.append(str(location))
                elif keyslist[i] == "CH1":
                    Heavy_chain_a_domains.append(CH1a)
                    location = VHa.get(keyslist[i])[0]
                    Label_Locations.append((CH1a_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker("VH1a",str(keyslist[i])+"a")
                    Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

        elif i >1:
            if VHa.get(keyslist[i])[0] == (int(VHa.get(keyslist[i-1])[0])+1) and VHa.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
                if keyslist[i] == "CH1":
                    Heavy_chain_a_domains.append(CH1a)
                    location = VHa.get(keyslist[i])[0]
                    Label_Locations.append((CH1a_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker("VH1a",str(keyslist[i])+"a")
                    Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

                elif keyslist[i] == "H":
                    location = VHa.get(keyslist[i])[0]
                    Label_Locations.append((Saltbridge_labela))
                    Label_Text.append(str(location))


                elif "CH2" in keyslist[i]:
                    Heavy_chain_a_domains.append(CH2a)
                    if VHa.get('CH2')[0] == VHb.get('CH2')[1] and VHa.get('CH2')[1] == VHb.get('CH2')[0]:
                        CH2a_CH2b = True
                        location = VHa.get(keyslist[i])[0]
                        Label_Locations.append((CH2a_label_location))
                        Label_Text.append(str(location))


                elif "CH3" in keyslist[i]:
                    Heavy_chain_a_domains.append(CH3a)
                    if VHa.get('CH3')[0] == VHb.get('CH3')[1] and VHa.get('CH3')[1] == VHb.get('CH3')[0]:
                        CH3a_CH3b = True
                        location = VHa.get(keyslist[i])[0]
                        Label_Locations.append((CH3a_label_location))
                        Label_Text.append(str(location))
                        bonds = bondmaker(str(keyslist[i-1])+"a",str(keyslist[i])+"a")
                        Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))


                #checking for extra chains at the end
                elif keyslist[i] == "VH2" or "CH4" and keyslist[i-1] == "CH3":
                    Heavy_chain_a_domains.append(IgG_2scFV1a)
                    location = VHa.get(keyslist[i])[0]
                    Label_Locations.append((IgG_2scFV1a_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker("CH3a","IgG_2scFV1a")
                    Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

                elif keyslist[i] == "VH2" and keyslist[i-1] == "CH4":
                    Heavy_chain_a_domains.append(IgG_2scFV2a)
                    location = VHa.get(keyslist[i])[0]
                    Label_Locations.append((IgG_2scFV2a_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker("IgG_2scFV1a","IgG_2scFV2a")
                    Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))



            elif VHa.get(keyslist[i])[0] == VHa.get(keyslist[i-1])[1] and  VHa.get(keyslist[i])[1] == VHa.get(keyslist[i-1])[0]:
                if keyslist[i] == "VL2":
                    if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH3" :
                        Heavy_chain_a_domains.append(IgG_2scFV3a)
                        location = VHa.get(keyslist[i])[0]
                        Label_Locations.append((IgG_2scFV3a_label_location))
                        Label_Text.append(str(location))
                        bonds = bondmaker("IgG_2scFV1a","IgG_2scFV3a")
                        Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

                    elif keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                        Heavy_chain_a_domains.append(IgG_2scFV4a)
                        location = VHa.get(keyslist[i])[0]
                        Label_Locations.append((IgG_4scFV1a_label_location))
                        Label_Text.append(str(location))
                        bonds = bondmaker("IgG_2scFV3a","IgG_2scFV4a")
                        Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

            elif VHa.get(keyslist[i])[0] != VHa.get(keyslist[i-1])[1] and  VHa.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
                if keyslist[i] == "VL2":
                    if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                        Heavy_chain_a_domains.append(IgG_2scFV3a)
                        location = VHa.get(keyslist[i])[0]
                        Label_Locations.append((IgG_2scFV3a_label_location))
                        Label_Text.append(str(location))
                        bonds = bondmaker("IgG_2scFV2a","IgG_2scFV4a")
                        Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))



#########Light chain A check what domains are there################
    for i in range(len(VLa)):
        keyslist = list(VLa.keys())
        if i==0:
            if keyslist[i] == "VL.a":
                Light_chain_a_domains.append(VL1a)
                location = VLa.get(keyslist[i])[0]
                Label_Locations.append((VL1a_label_location))
                Label_Text.append(str(location))

            elif keyslist[i] == "X":
                Salt_bridges.append(VLa_ADC)
                TYPE = str(re.sub("\[\'TYPE:|\]","",VLa.get(keyslist[i])[0]))
                Label_Locations.append(VLa_ADC_label_location)
                Label_Text.appen(str(TYPE))



        elif i ==1:
            if keyslist[i] == "VL.a":
                Light_chain_a_domains.append(VL1a)
                location = VLa.get(keyslist[i])[0]
                Label_Locations.append((VL1a_label_location))
                Label_Text.append(str(location))
            if keyslist[i] == "VH.a":
                Light_chain_a_domains.append(VH1a)
                location = VLa.get(keyslist[i])[0]
                Label_Locations.append((VH1a_label_location))
                Label_Text.append(str(location))
            if keyslist[i] == "VH.b":
                Light_chain_a_domains.append(VH1b)
                location = VLa.get(keyslist[i])[0]
                Label_Locations.append((VH1b_label_location))
                Label_Text.append(str(location))
            if keyslist[i] == "VL.b":
                Light_chain_a_domains.append(VL1b)
                location = VLa.get(keyslist[i])[0]
                Label_Locations.append((VL1b_label_location))
                Label_Text.append(str(location))
                bonds = bondmaker("VL1a","VL1b")
                Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

            if keyslist[i] == "CL" :
                if VLa.get(keyslist[i])[0] == (int(VLa.get(keyslist[i-1])[0])+1):
                    Light_chain_a_domains.append(CL1a)
                    location = VLa.get(keyslist[i])[0]
                    Label_Locations.append((CL1a_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker("VL1a",str(keyslist[i])+"1a")
                    Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

        elif i >1:
            if VLa.get(keyslist[i])[0] == (int(VLa.get(keyslist[i-1])[0])+1):
                if keyslist[i] == "CL":
                    if keyslist[i-1] == "VL.a":
                        Light_chain_a_domains.append(CL1a)
                        location = VLa.get(keyslist[i])[0]
                        Label_Locations.append((CL1a_label_location))
                        Label_Text.append(str(location))
                        bonds = bondmaker("VL1a",str(keyslist[i])+"1a")
                        Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))
                    elif keyslist[i-1] == "VL.b":
                        Light_chain_a_domains.append(CL1a)
                        location = VLa.get(keyslist[i])[0]
                        Label_Locations.append((CL1a_label_location))
                        Label_Text.append(str(location))
                        bonds = bondmaker("VL1b",str(keyslist[i])+"1a")
                        Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))



#########Heavy chain B check what domains are there################

    for i in range(len(VHb)):
            keyslist = list(VHb.keys())
            if i==0:
                if keyslist[i] == "VH.b":
                    Heavy_chain_b_domains.append(VH1b)
                    location = VHb.get(keyslist[i])[0]
                    Label_Locations.append((VH1b_label_location))
                    Label_Text.append(str(location))
                elif keyslist[i] == "X":
                    Heavy_chain_b_domains.append(VHb_ADC)
                    TYPE = str(re.sub("\[\'TYPE:|\]","",VHb.get(keyslist[i])[0]))
                    Label_Locations.append(VHb_ADC_label_location)
                    Label_Text.append(str(TYPE))




            elif i == 1:
                if keyslist[i] == "VH.b" and keyslist[i-1] == "X" :
                    Heavy_chain_b_domains.append(VH1b)
                    location = VHb.get(keyslist[i])[0]
                    Label_Locations.append((VH1b_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker("VHb_ADC","VH1b")
                    Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

                elif keyslist[i] == "VH.a" and keyslist[i-1] == "VH.b":
                    Heavy_chain_b_domains.append(VH1a)
                    location = VHb.get(keyslist[i])[0]
                    Label_Locations.append((VH1a_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker("VH1b","VH1a")
                    Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

                elif keyslist[i] == "CH1":
                    Heavy_chain_b_domains.append(CH1b)
                    location = VHb.get(keyslist[i])[0]
                    Label_Locations.append((CH1b_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker("VH1b",str(keyslist[i])+"b")
                    Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))


            elif i > 1:
                if VHb.get(keyslist[i])[0] == (int(VHb.get(keyslist[i-1])[0])+1) and VHb.get(keyslist[i])[1] != VHb.get(keyslist[i-1])[0]:
                    if keyslist[i] == "CH1":
                        if  keyslist[i-1] == "VH.b":
                            Heavy_chain_b_domains.append(CH1b)
                            location = VHb.get(keyslist[i])[0]
                            Label_Locations.append((CH1b_label_location))
                            Label_Text.append(str(location))
                            bonds = bondmaker("VH1b","CH1b")
                            Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))
                        elif keyslist[i-1] == "VH.a":
                            Heavy_chain_b_domains.append(CH1a)
                            location = VHb.get(keyslist[i])[0]
                            Label_Locations.append((CH1a_label_location))
                            Label_Text.append(str(location))
                            bonds = bondmaker("VH1a","CH1a")
                            Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

                    elif keyslist[i] == "H":
                        location = VHb.get(keyslist[i])[0]
                        Label_Locations.append((Saltbridge_labelb))
                        Label_Text.append(str(location))

                    elif keyslist[i] == "CH2":
                        Heavy_chain_b_domains.append(CH2b)
                        if VHb.get('CH2')[0] == VHa.get('CH2')[1] and VHb.get('CH2')[1] == VHa.get('CH2')[0]:
                            location = VHb.get(keyslist[i])[0]
                            Label_Locations.append((CH2b_label_location))
                            Label_Text.append(str(location))

                    elif keyslist[i] == "CH3":
                        Heavy_chain_b_domains.append(CH3b)
                        if VHb.get('CH3')[0] == VHa.get('CH3')[1] and VHb.get('CH3')[1] == VHa.get('CH3')[0]:
                            location = VHb.get(keyslist[i])[0]
                            Label_Locations.append((CH3b_label_location))
                            Label_Text.append(str(location))
                            bonds = bondmaker(str(keyslist[i-1])+"b",str(keyslist[i])+"b")
                            Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

                    #checking for extra chains at the end
                    elif keyslist[i] == "VH2" or "CH4" and keyslist[i-1] == "CH3":
                        Heavy_chain_b_domains.append(IgG_2scFV1b)
                        location = VHb.get(keyslist[i])[0]
                        Label_Locations.append((IgG_2scFV1b_label_location))
                        Label_Text.append(str(location))
                        bonds = bondmaker("CH3b","IgG_2scFV1b")
                        Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

                    elif keyslist[i] == "VH2" and keyslist[i-1] == "CH4":
                        Heavy_chain_b_domains.append(IgG_2scFV2b)
                        location = VHb.get(keyslist[i])[0]
                        Label_Locations.append((IgG_2scFV2b_label_location))
                        Label_Text.append(str(location))
                        bonds = bondmaker("IgG_2scFV1b","IgG_2scFV2b")
                        Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))



                elif VHb.get(keyslist[i])[0] == VHb.get(keyslist[i-1])[1] and  VHb.get(keyslist[i])[1] == VHb.get(keyslist[i-1])[0]:
                    if keyslist[i] == "VL2":
                        if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH3" :
                            Heavy_chain_b_domains.append(IgG_2scFV3b)
                            location = VHb.get(keyslist[i])[0]
                            Label_Locations.append((IgG_2scFV3b_label_location))
                            Label_Text.append(str(location))
                            bonds = bondmaker("IgG_2scFV1b","IgG_2scFV3b")
                            Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))


                        elif keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                            Heavy_chain_b_domains.append(IgG_2scFV4b)
                            location = VHb.get(keyslist[i])[0]
                            Label_Locations.append((IgG_4scFV4b_label_location))
                            Label_Text.append(str(location))
                            bonds = bondmaker("IgG_2scFV3b","IgG_2scFV4b")
                            Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))


                elif VHb.get(keyslist[i])[0] != VHb.get(keyslist[i-1])[1] and  VHb.get(keyslist[i])[1] != VHb.get(keyslist[i-1])[0]:
                    if keyslist[i] == "VL2":
                        if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                            Heavy_chain_b_domains.append(IgG_2scFV3b)
                            location = VHb.get(keyslist[i])[0]
                            Label_Locations.append((IgG_2scFV4b_label_location))
                            Label_Text.append(str(location))
                            bonds = bondmaker("IgG_2scFV2b","IgG_2scFV4b")
                            Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))



#########Light chain B check what domains are there################
    for i in range(len(VLb)):
        keyslist = list(VLb.keys())
        if i==0:
            if keyslist[i] == "VL.b":
                Light_chain_b_domains.append(VL1b)
                location = VLb.get(keyslist[i])[0]
                Label_Locations.append((VL1b_label_location))
                Label_Text.append(str(location))

            elif keyslist[i] == "X":
                Salt_bridges.append(VLb_ADC)
                TYPE = str(re.sub("\[\'TYPE:|\]","",VLb.get(keyslist[i])[0]))
                Label_Locations.append(VLb_ADC_label_location)
                Label_Text.append(str(TYPE))
                bonds = bondmaker("VLb_ADC","VL1b")
                Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

        elif i ==1:
            if keyslist[i] == "VL.b" and keyslist[i-1] == "X":
                Light_chain_b_domains.append(VL1b)
                location = VLb.get(keyslist[i])[0]
                Label_Locations.append((VL1b_label_location))
                Label_Text.append(str(location))

            if keyslist[i] == "VH.b" and keyslist[i-1] == "VL.b":
                Light_chain_b_domains.append(VH1b)
                location = VLb.get(keyslist[i])[0]
                Label_Locations.append((VH1b_label_location))
                Label_Text.append(str(location))
            if keyslist[i] == "VH.a" and keyslist[i-1] == "VL.b":
                Light_chain_b_domains.append(VH1a)
                location = VLb.get(keyslist[i])[0]
                Label_Locations.append((VH1a_label_location))
                Label_Text.append(str(location))

            if keyslist[i] == "CL" :
                if VLb.get(keyslist[i])[0] == (int(VLb.get(keyslist[i-1])[0])+1):
                    Light_chain_b_domains.append(CL1b)
                    location = VLb.get(keyslist[i])[0]
                    Label_Locations.append((CL1b_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker("VL1b",str(keyslist[i])+"1b")
                    Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))

        elif i >1:
            if VLb.get(keyslist[i])[0] == (int(VLb.get(keyslist[i-1])[0])+1):
                if keyslist[i] == "CL":
                    Light_chain_b_domains.append(CL1b)
                    location = VLb.get(keyslist[i])[0]
                    Label_Locations.append((CL1a_label_location))
                    Label_Text.append(str(location))
                if keyslist[i] == "CH1":
                    Light_chain_b_domains.append(CH1a)
                    location = VLb.get(keyslist[i])[0]
                    Label_Locations.append((CH1a_label_location))
                    Label_Text.append(str(location))

#########Check Fragment################

    if fragment != {}:
        for i in range(len(fragment)):
            keyslist = list(fragment.keys())
            if keyslist[i] == "VL.a":
                Light_chain_a_domains.append(VL1a)
                location = fragment.get(keyslist[i])[0]
                Label_Locations.append((VL1a_label_location))
                Label_Text.append(str(location))


            if keyslist[i] == "VH.a":
                Heavy_chain_a_domains.append(VH1a)
                location = fragment.get(keyslist[i])[0]
                Label_Locations.append((VH1a_label_location))
                Label_Text.append(str(location))
                if keyslist[i-1] == "VL.a":
                    if i > 0:
                        bonds = bondmaker("VL1a","VH1a")
                        Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))


            if keyslist[i] == "VH.b":
                Heavy_chain_b_domains.append(VH1b)
                location = fragment.get(keyslist[i])[0]
                Label_Locations.append((VH1b_label_location))
                Label_Text.append(str(location))
                if keyslist[i+1] == "VL.b":
                    Salt_bridges.append(VHb_VLb_bond_location)


            if keyslist[i] == "VL.b":
                Light_chain_b_domains.append(VL1b)
                location = fragment.get(keyslist[i])[0]
                Label_Locations.append((VL1b_label_location))
                Label_Text.append(str(location))
                bonds = bondmaker("VH1b","VL1b")
                Salt_bridges.append(eval(bonds[0])+eval(bonds[1]))


            if keyslist[i] == "X":
                Heavy_chain_a_domains.append(VHa_VHb_ADC)
                TYPE = str(re.sub("\[\'TYPE:|\]","",fragment.get(keyslist[i])[0]))
                Label_Locations.append((VHa_VHb_ADC_label_location))
                Label_Text.append(str(TYPE))
                if keyslist[i-1] == "VH.a":
                    Salt_bridges.append((VH1a_bottom_bond_location,VHa_ADC_bond_location))

                if keyslist[i+1] == "VH.b":
                    Salt_bridges.append((VHb_ADC_bond_location,VH1b_top_bond_location))


    return (Heavy_chain_a_domains,Light_chain_a_domains,Heavy_chain_b_domains,Light_chain_b_domains, Salt_bridges, Label_Locations,Label_Text)



def render(chains_list,canvas):
    canvas.delete("all")
    Polygons = [chains_list[0],chains_list[1],chains_list[2],chains_list[3]]
    Bonds    = chains_list[4]
    Label_positions = chains_list[5]
    Label_Text      = chains_list[6]


    for i in range(len(Polygons)):
        for j in range(len(Polygons[i])):
            canvas.create_polygon(Polygons[i][j], outline='#f11',fill='#1f1', width=2)
    for i in range(len(Bonds)):
        canvas.create_line(Bonds[i], fill='#f11', width = 2)
    for i in range(len(Label_positions)):
        canvas.create_text(Label_positions[i], text=Label_Text[i])
    #canvas.pack(fill=BOTH, expand=1)


def render_pipeline(entry,canvas):
    split_chains = Get_dictionaries(entry)
    coordinates  = Check_interactions(split_chains)
    render(coordinates, canvas)

root = tk.Tk()

canvas = tk.Canvas(root, height=HEIGHT, width=WIDTH, bg='#E7E0E6')
canvas.pack()

frame = tk.Frame(root, bg = '#80c1ff',bd=5)
frame.place(relx=0.05,rely=0.1,relwidth=0.4,relheight=0.8,)

entry = tk.Entry(frame, font=40)
entry.place(relwidth=1,relheight=0.5)

button = tk.Button(frame, text = "Get Structure", font=40, command=lambda: render_pipeline(entry.get(),lower_canvas))
button.place(relx=0.35,rely=0.55,relheight=0.1, relwidth=0.3)

lower_frame = tk.Frame(root, bg = '#80c1ff', bd=10)
lower_frame.place(relx=0.45, rely=0.1, relwidth=0.55,relheight=0.8)
lower_canvas = tk.Canvas(lower_frame)
lower_canvas.place(relheight=1,relwidth=1)




root.mainloop()
