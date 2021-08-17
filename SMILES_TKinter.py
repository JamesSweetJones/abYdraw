#!/usr/bin/python
import re
import sys
import os
import tkinter as tk
from graphics import *

HEIGHT = 700
WIDTH  = 800

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


    VL1a                       = [50,50,60,50,66,60,75,60,90,87,70,87]
    VL1a_no_V                  = [50,50,70,55,90,87,70,87]
    VL1a_label_location        = [40,70]
    VLa1_top_bond_location     = [60,50]
    VL1a_bottom_bond_location  = [80,87]


    CL1a                       = [70,90,90,90,110,127,90,127]
    CL1a_label_location        = [60,110]
    CL1a_top_bond_location     = [80,90]
    CL1a_bottom_bond_location  = [100,127]

    scFV_VeryL_IgGa = [0,10,20,10,40,50,30,50]
    tandem_scFca    = [10,25,20,25,35,45,20,45]
    Fab_scFV_Fca    = [20,45,30,45,40,65,30,65]
    scFV_VeryLscFVa = [30,65,40,65,40,85,40,85]
    scFV_L_IgGa     = [15,5,25,5,35,25,25,25]
    scFV_LscFVa     = [45,65,55,65,55,85,45,85]

#########VHa Chain coordinates################

    VH1a                      = [85,60,95,60,90,50,100,50,120,87,100,87]
    VH1a_no_V                 = [80,50,100,50,120,87,100,87]
    VH1a_label_location       = [120,70]
    VH1a_top_bond_location    = [90,50]
    VH1a_bottom_bond_location = [110,87]

    CH1a                      = [100,90,120,90,140,127,120,127]
    CH1a_label_location       = [140,110]
    CH1a_top_bond_location    = [110,90]
    CH1a_bottom_bond_location = [130,127]

    CH2a                      = [ 140,140,160,140,160,177,140,177]
    CH2a_label_location       = [130,160]
    CH2a_top_bond_location    = [150,140]
    CH2a_bottom_bond_location = [ 150,177]

    CH3a                      = [ 140,180,160,180,160,217,140,217]
    CH3a_notch                = [ 140,180,160,180,170,240,160,217,140,217]
    CH3a_label_location       = [130,200]
    CH3a_top_bond_location    = [150,180]
    CH3a_bottom_bond_location = [ 150,217]


    scFV_H_IgGa_outer = [ 60,10,80,10,100,50,80,50]
    scFV_H_IgGa_inner = [ 90,10,110,10,130,50,110,50]
    scFv4_IgG    = [ 110,50,130,50,150,90,130,90]

    IgG_2scFV1a                      = [ 140,220,160,220,160,257,140,257]
    IgG_2scFV1a_label_location       = [150,240]
    IgG_2scFV1a_top_bond_location    = [150,220]
    IgG_2scFV1a_bottom_bond_location = [ 150,257]

    IgG_2scFV2a                      = [ 140,260,180,260,180,297,140,297]
    IgG_2scFV2a_label_location       = [150,280]
    IgG_2scFV2a_top_bond_location    = [160,260]
    IgG_2scFV2a_bottom_bond_location = [ 160,297]

    IgG_2scFV3a                      = [ 110,220,130,220,130,257,110,257]
    IgG_2scFV3a_label_location       = [120,240]
    IgG_2scFV3a_top_bond_location    = [120,220]
    IgG_2scFV3a_bottom_bond_location = [ 120,257]

    IgG_2scFV4a                      = [ 110,260,130,260,130,297,110,297]
    IgG_2scFV4a_label_location       = [120,280]
    IgG_2scFV4a_top_bond_location    = [120,260]
    IgG_2scFV4a_bottom_bond_location = [ 120,297]



#########VLb Chain coordinates################

    VL1b                      = [ 280,50,270,50,265,60,255,60,240,87,260,87]
    VL1b_no_V                 = [ 280,50,260,50,240,87,260,87]
    VL1b_label_location       = [280,70]
    VL1b_top_bond_location    = [270,50]
    VL1b_bottom_bond_location = [ 250,87]


    CL1b                      = [ 260,90,240,90,220,127,240,127]
    CL1b_label_location       = [260,110]
    CL1b_top_bond_location    = [250,90]
    CL1b_bottom_bond_location = [ 230,127]

    scFV_VeryL_IgGa = [ 165,5,155,5,145,25,165,5]
    tandem_scFca= [ 155,25, 145,25,135,45,155,45]
    Fab_scFV_Fca= [ 145,45,135,45,125,65,135,65]
    scFV_VeryLscFVa = [ 135,65,125,65,125,85,135,85]
    scFV_L_IgGb = [ 150,5,140,5,130,25,140,25]
    scFV_LscFVb = [ 120,60,110,60,110,80,120,80]

#########VHb Chain coordinates################

    VH1b                      = [ 230,50, 240,50, 235,60, 245,60,230,87,210,87]
    VH1b_no_v                 = [ 250,50,230,50,210,87,230,87]
    VH1b_label_location       = [210,70]
    VH1b_top_bond_location    = [240,50]
    VH1b_bottom_bond_location = [ 220,87]

    CH1b                      = [ 230,90,210,90,190,127,210,127]
    CH1b_label_location       = [190,110]
    CH1b_top_bond_location    = [220,90]
    CH1b_bottom_bond_location = [200,127]


    CH2b                      = [190,140, 170,140,170,177,190,177]
    CH2b_label_location       = [200,160]
    CH2b_top_bond_location    = [180,140]
    CH2b_bottom_bond_location = [180,177]


    CH3b                      = [ 190,180, 170,180,170,217,190,217]
    CH3b_label_location       = [200,200]
    CH3b_top_bond_location    = [180,180]
    CH3b_bottom_bond_location = [ 180,217]


    scFv_H_IgGb_outer                = [270,10,250,10,130,50,250,50]
    scFv_H_IgGb_inner                = [ 240,10,220,10,200,50,220,50]
    scFv4_IgGb                       = [ 220,50,200,50,180,90,200,90]

    IgG_2scFV1b                      = [ 190,220,170,220,170,257,190,257]
    IgG_2scFV1b_label_location       = [180,240]
    IgG_2scFV1b_top_bond_location    = [180,220]
    IgG_2scFV1b_bottom_bond_location = [ 180,257]

    IgG_2scFV2b                      = [ 190,260,170,260,170,297,190,297]
    IgG_2scFV2b_label_location       = [180,300]
    IgG_2scFV2b_top_bond_location    = [180,260]
    IgG_2scFV2b_bottom_bond_location = [180,297]

    IgG_2scFV3b                      = [ 220,220,200,220,200,257,220,257]
    IgG_2scFV3b_label_location       = [210,240]
    IgG_2scFV3b_top_bond_location    = [210,220]
    IgG_2scFV3b_bottom_bond_location = [ 210,257]

    IgG_2scFV4b                      = [ 220,260,200,260,200,297,220,297]
    IgG_2scFV4b_label_location       = [210,300]
    IgG_2scFV4b_top_bond_location    = [210,260]
    IgG_2scFV4b_bottom_bond_location = [ 210,297]



#########Check Salt Bridge ################

    CH1_CH2bonda = CH1a_bottom_bond_location+CH2a_top_bond_location
    CH1_CH2bondb = CH1b_bottom_bond_location+CH2b_top_bond_location
    Saltbridge_a  =[134,131,194,131]
    Saltbridge_labela = [130,137]
    Saltbridge_b  = [144,136,187,136]
    Saltbridge_labelb = [200,137]

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
    VLa_ADC                    = [60,50,50,30,40,20,50,10,60,20,50,30]
    VLa_ADC_label_location     = [ 50,10]
    VHa_ADC                    = [90,50,80,30,70,20,80,10,90,20,80,30]
    VHa_ADC_label_location     = [ 80,10]
    VLb_ADC                    = [280,50,290,30,300,20,290,10,280,20,290,30]
    VLb_ADC_label_location     = [ 230,10]
    VHb_ADC                    = [240,50,250,30,240,20,250,10,260,20,250,30]
    VHb_ADC_label_location     = [ 200,10]
    VHa_VHb_ADC                = [150,90,140,70,150,50,180,50,190,70,180,90,150,90]
    VHa_VHb_ADC_label_location = [ 160,40]
    VHa_ADC_bond_location      = [ 150,50]
    VHb_ADC_bond_location      = [ 180,90]


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
                Salt_bridges.append(VHa_ADC)
                TYPE = str(re.sub("\[\'TYPE:|\]","",VHa.get(keyslist[i])[0]))
                Label_Locations.append((VHa_ADC_label_location,str(TYPE)))




        elif i == 1:
                if keyslist[i] == "VH.a" and keyslist[i-1] == "X":
                    Heavy_chain_a_domains.append(VH1a)
                    location = VHa.get(keyslist[i])[0]
                    Label_Locations.append((VH1a_label_location))
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
                Label_Locations.append((VLa_ADC_label_location,str(TYPE)))


        elif i ==1:
            if keyslist[i] == "VL.a":
                Light_chain_a_domains.append(VL1a)
                location = VLa.get(keyslist[i])[0]
                Label_Locations.append((VL1a_label_location))
                Label_Text.append(str(location))

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
                    Light_chain_a_domains.append(CL1a)
                    location = VLa.get(keyslist[i])[0]
                    Label_Locations.append((CL1a_label_location))
                    Label_Text.append(str(location))
                    bonds = bondmaker(str(VL1a),str(keyslist[i])+"1a")
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
                    Salt_bridges.append(VHb_ADC)
                    TYPE = str(re.sub("\[\'TYPE:|\]","",VHb.get(keyslist[i])[0]))
                    Label_Locations.append((VHb_ADC_label_location,str(TYPE)))




            elif i == 1:
                if keyslist[i] == "VH.b" and keyslist[i-1] == "X" :
                    Heavy_chain_b_domains.append(VH1b)
                    location = VHb.get(keyslist[i])[0]
                    Label_Locations.append((VH1b_label_location))
                    Label_Text.append(str(location))

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
                        Heavy_chain_b_domains.append(CH1b)
                        location = VHb.get(keyslist[i])[0]
                        Label_Locations.append((CH1b_label_location))
                        Label_Text.append(str(location))

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
                Label_Locations.append((VLb_ADC_label_location,str(TYPE)))


        elif i ==1:
            if keyslist[i] == "VL.b":
                Light_chain_b_domains.append(VL1b)
                location = VLb.get(keyslist[i])[0]
                Label_Locations.append((VL1b_label_location))
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
                Salt_bridges.append(VHa_VHb_ADC)
                TYPE = str(re.sub("\[\'TYPE:|\]","",fragment.get(keyslist[i])[0]))
                Label_Locations.append((VHa_VHb_ADC_label_location,str(TYPE)))
                if keyslist[i-1] == "VH.a":
                    Salt_bridges.append(Line(VH1a_bottom_bond_location,VHa_ADC_bond_location))

                if keyslist[i+1] == "VH.b":
                    Salt_bridges.append(Line(VHb_ADC_bond_location,VH1b_top_bond_location))


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
        canvas.create_line(Bonds[i])
    for i in range(len(Label_positions)):
        canvas.create_text(Label_positions[i], text=Label_Text[i])
    #canvas.pack(fill=BOTH, expand=1)


def render_pipeline(entry,canvas):
    split_chains = Get_dictionaries(entry)
    coordinates  = Check_interactions(split_chains)
    render(coordinates, canvas)

root = tk.Tk()

canvas = tk.Canvas(root, height=HEIGHT, width=WIDTH, bg='#FFC433')
canvas.pack()

frame = tk.Frame(root, bg = '#80c1ff',bd=5)
frame.place(relx=0.5,rely=0.1,relwidth=0.75,relheight=0.1,anchor='n')

entry = tk.Entry(frame, font=40)
entry.place(relwidth=0.65,relheight=1)

button = tk.Button(frame, text = "Get Structure", font=40, command=lambda: render_pipeline(entry.get(),lower_canvas))
button.place(relx=0.7, relheight=1, relwidth=0.3)

lower_frame = tk.Frame(root, bg = '#80c1ff', bd=10)
lower_frame.place(relx=0.5, rely=0.25, relwidth=0.75,relheight=0.6,anchor='n')
lower_canvas = tk.Canvas(lower_frame)
lower_canvas.place(relheight=1,relwidth=1)




root.mainloop()
