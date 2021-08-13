#!/usr/bin/python



"""
Graphics module found at https://mcsp.wartburg.edu/zelle/python/graphics.py

NEED TO DO
>Heavy chains pre VH domain
>Light chain adaptations
"""
import re
import sys
from graphics import *

#normal mAb
input1="VH.b(1:6)-CH1(2:7){1}-H(3:10){2}-CH2(4:11) - CH3(5:12) | VL.b(6:1)-CL(7:2){1} | VH.a(8:13)- CH1(9:14){1}-H(10:3){2}-CH2(11:4)-CH3(12:5) | VL.a(13:8)-CL(14:9){1}"
#IgG(H)-scFv
input2="VH.b(1:8)-CH1(2:9){1}-H(3:12){2}-CH2(4:13)-CH3(5:14)-L-VH2(6:7)-L-VL2(7:6)|VL.b(8:1)-CL(9:2){1}|VH.a(10:17)-CH1(11:18){1}-H(12:3){2}-CH2(13:4)-CH3(14:5)-L-VH2(15:16)-L-VL2(16:15)|VL.a(17:10)-CL(18:11){1}"

input3="VL.a(1:2)-L-VH.a(2:1)-L-X[TYPE:FUSION][NOTE:human serum albumin]-L-VH.b(3:4)-L-VL.b(4:3)"
#BsAb Fragment
input4="VH.b(1:2)-L-VL.b(2:1)-H*(3:6){1}[MOD: engineered disulphide bond]-X[TYPE:LEUCINE] | VH.a(4:5)-L-VL.a(5:4)-H*(6:3){1}[MOD: engineered disulphide bond]-X[TYPE:LEUCINE]"
#ADC
input5="X[TYPE: pharmacophore peptide heterodimer]-VH.b(1:6)- CH1(2:7){1}-H(3:10){2}-CH2(4:11)-CH3(5:12) | VL.b(6:1)-CL(7:2){1} | X[TYPE: pharmacophore peptide heterodimer]-VH.a(8:14)-CH1(9:14){1}-H(10:3){2}- CH2(11:4)-CH3 (12:5) | VL.a(13:8)-CL(14:9){1}"
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
        #print(dict,chain)
        for j in range(len(chain)):
            domain   =  re.findall("^CH[1-9]|^VL.[ab]|^VH.[ab]|^VL[1-9]|^VH[1-9]|^CL|^VH|^VL|^H|^X", str(chain[j]))
            domain   =  str(re.sub("\[|\'|\]","", str(domain)))
            print(domain)
            print(chain[j])
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


    print(fragment)
    return(VHa,VLa,VHb,VLb,Salt_bridges,VHa_VLa_bonds,VHb_VLb_bonds,CH1a_CL1a_bonds,CH1b_CL1b_bonds,fragment)



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
    Labels                = []


    VHa_VLa   = False
    CH1a_CLa  = False
    VHb_VLb   = False
    CH1b_CLb  = False
    CH2a_CH2b = False
    CH3a_CH3b = False

#########VLa Chain coordinates################

    VL1a = Polygon(Point(50,50),Point(70,50),Point(90,90),Point(70,90))
    VL1a_label_location = Point(40,70)
    CL1a = Polygon(Point(70,90),Point(90,90),Point(110,130),Point(90,130))
    CL1a_label_location = Point(60,110)
    scFV_VeryL_IgGa = Polygon( Point(0,10),Point(20,10),Point(40,50),Point(30,50))
    tandem_scFca= Polygon( Point(10,25),Point(20,25),Point(35,45),Point(20,45))
    Fab_scFV_Fca= Polygon( Point(20,45),Point(30,45),Point(40,65),Point(30,65))
    scFV_VeryLscFVa = Polygon( Point(30,65),Point(40,65),Point(40,85),Point(40,85))
    scFV_L_IgGa = Polygon( Point(15,5),Point(25,5),Point(35,25),Point(25,25))
    scFV_LscFVa = Polygon( Point(45,65),Point(55,65),Point(55,85),Point(45,85))

#########VHa Chain coordinates################

    VH1a = Polygon(Point(80,50),Point(100,50),Point(120,90),Point(100,90))
    VH1a_label_location = Point(120,70)
    CH1a = Polygon( Point(100,90),Point(120,90),Point(140,130),Point(120,130))
    CH1a_label_location = Point(140,110)
    CH2a = Polygon( Point(140,140),Point(160,140),Point(160,180),Point(140,180))
    CH2a_label_location = Point(130,160)
    CH3a = Polygon( Point(140,180),Point(160,180),Point(160,220),Point(140,220))
    CH3a_label_location = Point(130,200)
    scFV_H_IgGa_outer = Polygon( Point(60,10),Point(80,10),Point(100,50),Point(80,50))
    scFV_H_IgGa_inner = Polygon( Point(90,10),Point(110,10),Point(130,50),Point(110,50))
    scFv4_IgG    = Polygon( Point(110,50),Point(130,50),Point(150,90),Point(130,90))
    IgG_2scFV1a  = Polygon( Point(110,220),Point(130,220),Point(130,260),Point(110,260))
    IgG_2scFV1a_label_location = Point(120,240)
    IgG_2scFV2a  = Polygon( Point(110,260),Point(130,260),Point(130,300),Point(110,300))
    IgG_2scFV2a_label_location = Point(120,280)
    IgG_2scFV3a  = Polygon( Point(140,220),Point(160,220),Point(160,260),Point(140,260))
    IgG_2scFV3a_label_location = Point(150,240)
    IgG_2scFV4a  = Polygon( Point(140,260),Point(180,260),Point(180,300),Point(140,300))
    IgG_2scFV4a_label_location = Point(150,280)

#########VLb Chain coordinates################

    VL1b = Polygon( Point(280,50),Point(260,50),Point(240,90),Point(260,90))
    VL1b_label_location = Point(280,70)
    CL1b = Polygon( Point(260,90),Point(240,90),Point(220,130),Point(240,130))
    CL1b_label_location = Point(260,110)
    scFV_VeryL_IgGa = Polygon( Point(165,5),Point(155,5),Point(145,25),Point(165,5))
    tandem_scFca= Polygon( Point(155,25), Point(145,25),Point(135,45),Point(155,45))
    Fab_scFV_Fca= Polygon( Point(145,45),Point(135,45),Point(125,65),Point(135,65))
    scFV_VeryLscFVa = Polygon( Point(135,65),Point(125,65),Point(125,85),Point(135,85))
    scFV_L_IgGb = Polygon( Point(150,5),Point(140,5),Point(130,25),Point(140,25))
    scFV_LscFVb = Polygon( Point(120,60),Point(110,60),Point(110,80),Point(120,80))

#########VHb Chain coordinates################


    VH1b = Polygon( Point(250,50),Point(230,50),Point(210,90),Point(230,90))
    VH1b_label_location = Point(210,70)
    CH1b= Polygon( Point(230,90),Point(210,90),Point(190,130),Point(210,130))
    CH1b_label_location = Point(190,110)
    CH2b= Polygon( Point(190,140), Point(170,140),Point(170,180),Point(190,180))
    CH2b_label_location = Point(200,160)
    CH3b= Polygon( Point(190,180), Point(170,180),Point(170,220),Point(190,220))
    CH3b_label_location = Point(200,200)
    scFv_H_IgGb_outer = Polygon(Point(270,10),Point(250,10),Point(130,50),Point(250,50))
    scFv_H_IgGb_inner= Polygon( Point(240,10),Point(220,10),Point(200,50),Point(220,50))
    scFv4_IgGb   = Polygon( Point(220,50),Point(200,50),Point(180,90),Point(200,90))
    IgG_2scFV1b  = Polygon( Point(220,220),Point(200,220),Point(200,260),Point(220,260))
    IgG_2scFV1b_label_location = Point(180,240)
    IgG_2scFV2b  = Polygon( Point(220,260),Point(200,260),Point(200,300),Point(220,300))
    IgG_2scFV2b_label_location = Point(180,300)
    IgG_2scFV3b  = Polygon( Point(190,220),Point(170,220),Point(170,260),Point(190,260))
    IgG_2scFV3b_label_location = Point(210,240)
    IgG_2scFV4b  = Polygon( Point(190,260),Point(170,260),Point(170,300),Point(190,300))
    IgG_2scFV4b_label_location = Point(210,300)

#########Check Salt Bridge ################

    CH1_CH2bonda = Line(Point(130,130),Point(150,140))
    CH1_CH2bondb = Line(Point(200,130),Point(175,140))
    Saltbridge_a  = Line(Point(133,134),Point(194,134))
    Saltbridge_labela = Point(130,140)
    Saltbridge_b  = Line(Point(144,136),Point(190,136))
    Saltbridge_labelb = Point(200,140)

    if Salt_bridge_count   == 0:
        Salt_bridges = [CH1_CH2bonda, CH1_CH2bondb]
    elif Salt_bridge_count == 1:
        Salt_bridges = [CH1_CH2bonda, CH1_CH2bondb, Saltbridge_a]
    elif Salt_bridge_count == 2:
        Salt_bridges = [CH1_CH2bonda, CH1_CH2bondb, Saltbridge_a, Saltbridge_b]

#########Check Bonds between VH and VL Chains################

    VHa_VLa_bond_location   = Line(Point(90,90),Point(100,90))
    VHb_VLb_bond_location   = Line(Point(230,90),Point(240,90))
    CH1a_CL1a_bond_location = Line(Point(110,130),Point(120,130))
    CH1b_CL1b_bond_location = Line(Point(220,130),Point(210,130))
    if VHa_VLa_bond_count != "":
        Salt_bridges.append(VHa_VLa_bond_location)
    if VHb_VLb_bond_count!= "":
        Salt_bridges.append(VHb_VLb_bond_location)
    if CH1a_CL1a_bond_count != "":
        Salt_bridges.append(CH1a_CL1a_bond_location)
    if CH1b_CL1b_bond_count !="":
        Salt_bridges.append(CH1b_CL1b_bond_location)

#########Check VH/VL ADC Coordinates################
    VLa_ADC = Polygon(Point(60,50),Point(50,30),Point(40,20),Point(50,10),Point(60,20),Point(50,30))
    VLa_ADC_label_location = Point(50,10)
    VHa_ADC = Polygon(Point(90,50),Point(80,30),Point(70,20),Point(80,10),Point(90,20),Point(80,30))
    VHa_ADC_label_location = Point(80,10)
    VLb_ADC = Polygon(Point(280,50),Point(290,30),Point(300,20),Point(290,10),Point(280,20),Point(290,30))
    VLb_ADC_label_location = Point(230,10)
    VHb_ADC = Polygon(Point(240,50),Point(250,30),Point(240,20),Point(250,10),Point(260,20),Point(250,30))
    VHb_ADC_label_location = Point(200,10)
    VHa_VHb_ADC = Polygon(Point(150,90),Point(140,70),Point(150,50),Point(180,50),Point(190,70),Point(180,90),Point(150,90))
    VHa_VHb_ADC_label_location = Point(160,40)
    VHa_ADC_bond_location = Line(Point(110,90),Point(140,70))
    VHb_ADC_bond_location = Line(Point(190,70),Point(210,90))
#######Heavy chain A check what domains are there################
    for i in range(len(VHa)):
        keyslist = list(VHa.keys())
        if i==0:
            if keyslist[i] == "VH.a":
                Heavy_chain_a_domains.append(VH1a)
                location = VHa.get(keyslist[i])[0]
                Labels.append(Text(VH1a_label_location,str(location)))
            elif keyslist[i] == "X":
                Salt_bridges.append(VHa_ADC)
                TYPE = str(re.sub("\[\'TYPE:|\]","",VHa.get(keyslist[i])[0]))
                Labels.append(Text(VHa_ADC_label_location,str(TYPE)))




        elif i == 1:
                if keyslist[i] == "VH.a" and keyslist[i-1] == "X":
                    Heavy_chain_a_domains.append(VH1a)
                    location = VHa.get(keyslist[i])[0]
                    Labels.append(Text(VH1a_label_location,str(location)))
                elif keyslist[i] == "CH1":
                    Heavy_chain_a_domains.append(CH1a)
                    location = VHa.get(keyslist[i])[0]
                    Labels.append(Text(CH1a_label_location,str(location)))

        elif i >1:
            if VHa.get(keyslist[i])[0] == (int(VHa.get(keyslist[i-1])[0])+1) and VHa.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
                if keyslist[i] == "CH1":
                    Heavy_chain_a_domains.append(CH1a)
                    location = VHa.get(keyslist[i])[0]
                    Labels.append(Text(CH1a_label_location,str(location)))

                elif keyslist[i] == "H":
                    location = VHa.get(keyslist[i])[0]
                    Labels.append(Text(Saltbridge_labela,str(location)))

                elif keyslist[i] == "CH2":
                    Heavy_chain_a_domains.append(CH2a)
                    if VHa.get('CH2')[0] == VHb.get('CH2')[1] and VHa.get('CH2')[1] == VHb.get('CH2')[0]:
                        CH2a_CH2b = True
                        location = VHa.get(keyslist[i])[0]
                        Labels.append(Text(CH2a_label_location,str(location)))

                elif keyslist[i] == "CH3":
                    Heavy_chain_a_domains.append(CH3a)
                    if VHa.get('CH3')[0] == VHb.get('CH3')[1] and VHa.get('CH3')[1] == VHb.get('CH3')[0]:
                        CH3a_CH3b = True
                        location = VHa.get(keyslist[i])[0]
                        Labels.append(Text(CH3a_label_location,str(location)))

                #checking for extra chains at the end
                elif keyslist[i] == "VH2" or "CH4" and keyslist[i-1] == "CH3":
                    Heavy_chain_a_domains.append(IgG_2scFV1a)
                    location = VHa.get(keyslist[i])[0]
                    Labels.append(Text(IgG_2scFV1a_label_location,str(location)))
                elif keyslist[i] == "VH2" and keyslist[i-1] == "CH4":
                    Heavy_chain_a_domains.append(IgG_2scFV2a)
                    location = VHa.get(keyslist[i])[0]
                    Labels.append(Text(IgG_2scFV2a_label_location,str(location)))



            elif VHa.get(keyslist[i])[0] == VHa.get(keyslist[i-1])[1] and  VHa.get(keyslist[i])[1] == VHa.get(keyslist[i-1])[0]:
                if keyslist[i] == "VL2":
                    if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH3" :
                        Heavy_chain_a_domains.append(IgG_2scFV3a)
                        location = VHa.get(keyslist[i])[0]
                        Labels.append(Text(IgG_2scFV3a_label_location,str(location)))
                    elif keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                        Heavy_chain_a_domains.append(IgG_2scFV4a)
                        location = VHa.get(keyslist[i])[0]
                        Labels.append(Text(IgG_4scFV1a_label_location,str(location)))
            elif VHa.get(keyslist[i])[0] != VHa.get(keyslist[i-1])[1] and  VHa.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
                if keyslist[i] == "VL2":
                    if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                        Heavy_chain_a_domains.append(IgG_2scFV3a)
                        location = VHa.get(keyslist[i])[0]
                        Labels.append(Text(IgG_2scFV3a_label_location,str(location)))


#########Light chain A check what domains are there################
    for i in range(len(VLa)):
        keyslist = list(VLa.keys())
        if i==0:
            if keyslist[i] == "VL.a":
                Light_chain_a_domains.append(VL1a)
                location = VLa.get(keyslist[i])[0]
                Labels.append(Text(VL1a_label_location,str(location)))

            elif keyslist[i] == "X":
                Salt_bridges.append(VLa_ADC)
                TYPE = str(re.sub("\[\'TYPE:|\]","",VLa.get(keyslist[i])[0]))
                Labels.append(Text(VLa_ADC_label_location,str(TYPE)))


        elif i ==1:
            if keyslist[i] == "VL.a":
                Light_chain_a_domains.append(VL1a)
                location = VLa.get(keyslist[i])[0]
                Labels.append(Text(VL1a_label_location,str(location)))

            if keyslist[i] == "CL" :
                if VLa.get(keyslist[i])[0] == (int(VLa.get(keyslist[i-1])[0])+1):
                    Light_chain_a_domains.append(CL1a)
                    location = VLa.get(keyslist[i])[0]
                    Labels.append(Text(CL1a_label_location,str(location)))

        elif i >1:
            if VLa.get(keyslist[i])[0] == (int(VLa.get(keyslist[i-1])[0])+1):
                if keyslist[i] == "CL":
                    Light_chain_a_domains.append(CL1a)
                    location = VLa.get(keyslist[i])[0]
                    Labels.append(Text(CL1a_label_location,str(location)))



#########Heavy chain B check what domains are there################

    for i in range(len(VHb)):
            keyslist = list(VHb.keys())
            if i==0:
                if keyslist[i] == "VH.b":
                    Heavy_chain_b_domains.append(VH1b)
                    location = VHb.get(keyslist[i])[0]
                    Labels.append(Text(VH1b_label_location,str(location)))
                elif keyslist[i] == "X":
                    Salt_bridges.append(VHb_ADC)
                    TYPE = str(re.sub("\[\'TYPE:|\]","",VHb.get(keyslist[i])[0]))
                    Labels.append(Text(VHb_ADC_label_location,str(TYPE)))




            elif i == 1:
                if keyslist[i] == "VH.b" and keyslist[i-1] == "X" :
                    Heavy_chain_b_domains.append(VH1b)
                    location = VHb.get(keyslist[i])[0]
                    Labels.append(Text(VH1b_label_location,str(location)))

                elif keyslist[i] == "CH1":
                    Heavy_chain_b_domains.append(CH1b)
                    location = VHb.get(keyslist[i])[0]
                    Labels.append(Text(CH1b_label_location,str(location)))

            elif i > 1:
                if VHb.get(keyslist[i])[0] == (int(VHb.get(keyslist[i-1])[0])+1) and VHb.get(keyslist[i])[1] != VHb.get(keyslist[i-1])[0]:
                    if keyslist[i] == "CH1":
                        Heavy_chain_b_domains.append(CH1b)
                        location = VHb.get(keyslist[i])[0]
                        Labels.append(Text(CH1b_label_location,str(location)))

                    elif keyslist[i] == "H":
                        location = VHb.get(keyslist[i])[0]
                        Labels.append(Text(Saltbridge_labelb,str(location)))

                    elif keyslist[i] == "CH2":
                        Heavy_chain_b_domains.append(CH2b)
                        if VHb.get('CH2')[0] == VHa.get('CH2')[1] and VHb.get('CH2')[1] == VHa.get('CH2')[0]:
                            location = VHb.get(keyslist[i])[0]
                            Labels.append(Text(CH2b_label_location,str(location)))

                    elif keyslist[i] == "CH3":
                        Heavy_chain_b_domains.append(CH3b)
                        if VHb.get('CH3')[0] == VHa.get('CH3')[1] and VHb.get('CH3')[1] == VHa.get('CH3')[0]:
                            location = VHb.get(keyslist[i])[0]
                            Labels.append(Text(CH3b_label_location,str(location)))

                    #checking for extra chains at the end
                    elif keyslist[i] == "VH2" or "CH4" and keyslist[i-1] == "CH3":
                        Heavy_chain_b_domains.append(IgG_2scFV1b)
                        location = VHb.get(keyslist[i])[0]
                        Labels.append(Text(IgG_2scFV1b_label_location,str(location)))
                    elif keyslist[i] == "VH2" and keyslist[i-1] == "CH4":
                        Heavy_chain_b_domains.append(IgG_2scFV2b)
                        location = VHb.get(keyslist[i])[0]
                        Labels.append(Text(IgG_2scFV2b_label_location,str(location)))



                elif VHb.get(keyslist[i])[0] == VHb.get(keyslist[i-1])[1] and  VHb.get(keyslist[i])[1] == VHb.get(keyslist[i-1])[0]:
                    if keyslist[i] == "VL2":
                        if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH3" :
                            Heavy_chain_b_domains.append(IgG_2scFV3b)
                            location = VHb.get(keyslist[i])[0]
                            Labels.append(Text(IgG_2scFV3b_label_location,str(location)))
                        elif keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                            Heavy_chain_b_domains.append(IgG_2scFV4b)
                            location = VHb.get(keyslist[i])[0]
                            Labels.append(Text(IgG_4scFV1b_label_location,str(location)))
                            print("yay")
                elif VHb.get(keyslist[i])[0] != VHb.get(keyslist[i-1])[1] and  VHb.get(keyslist[i])[1] != VHb.get(keyslist[i-1])[0]:
                    if keyslist[i] == "VL2":
                        if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                            Heavy_chain_b_domains.append(IgG_2scFV3b)
                            location = VHb.get(keyslist[i])[0]
                            Labels.append(Text(IgG_2scFV3b_label_location,str(location)))


#########Light chain B check what domains are there################
    for i in range(len(VLb)):
        keyslist = list(VLb.keys())
        if i==0:
            if keyslist[i] == "VL.b":
                Light_chain_b_domains.append(VL1b)
                location = VLb.get(keyslist[i])[0]
                Labels.append(Text(VL1b_label_location,str(location)))

            elif keyslist[i] == "X":
                Salt_bridges.append(VLb_ADC)
                TYPE = str(re.sub("\[\'TYPE:|\]","",VLb.get(keyslist[i])[0]))
                Labels.append(Text(VLb_ADC_label_location,str(TYPE)))


        elif i ==1:
            if keyslist[i] == "VL.b":
                Light_chain_b_domains.append(VL1b)
                location = VLb.get(keyslist[i])[0]
                Labels.append(Text(VL1b_label_location,str(location)))

            if keyslist[i] == "CL" :
                if VLb.get(keyslist[i])[0] == (int(VLb.get(keyslist[i-1])[0])+1):
                    Light_chain_b_domains.append(CL1b)
                    location = VLb.get(keyslist[i])[0]
                    Labels.append(Text(CL1b_label_location,str(location)))

        elif i >1:
            if VLb.get(keyslist[i])[0] == (int(VLb.get(keyslist[i-1])[0])+1):
                if keyslist[i] == "CL":
                    Light_chain_b_domains.append(CL1b)
                    location = VLb.get(keyslist[i])[0]
                    Labels.append(Text(CL1a_label_location,str(location)))

#########Check Fragment################

    if fragment != {}:
        for i in range(len(fragment)):
            keyslist = list(fragment.keys())
            print(keyslist[i])
            if keyslist[i] == "VL.a":
                Light_chain_a_domains.append(VL1a)
                location = fragment.get(keyslist[i])[0]
                Labels.append(Text(VL1a_label_location,str(location)))

            if keyslist[i] == "VH.a":
                Heavy_chain_a_domains.append(VH1a)
                location = fragment.get(keyslist[i])[0]
                Labels.append(Text(VH1a_label_location,str(location)))
                if keyslist[i-1] == "VL.a":
                    Salt_bridges.append(VHa_VLa_bond_location)

            if keyslist[i] == "VH.b":
                Heavy_chain_b_domains.append(VH1b)
                location = fragment.get(keyslist[i])[0]
                Labels.append(Text(VH1b_label_location,str(location)))
                if keyslist[i+1] == "VL.b":
                    Salt_bridges.append(VHb_VLb_bond_location)

            if keyslist[i] == "VL.b":
                Light_chain_b_domains.append(VL1b)
                location = fragment.get(keyslist[i])[0]
                Labels.append(Text(VL1b_label_location,str(location)))

            if keyslist[i] == "X":
                Salt_bridges.append(VHa_VHb_ADC)
                TYPE = str(re.sub("\[\'TYPE:|\]","",fragment.get(keyslist[i])[0]))
                Labels.append(Text(VHa_VHb_ADC_label_location,str(TYPE)))
                if keyslist[i-1] == "VH.a":
                    Salt_bridges.append(VHa_ADC_bond_location)
                if keyslist[i+1] == "VH.b":
                    Salt_bridges.append(VHb_ADC_bond_location)

    return (Heavy_chain_a_domains,Light_chain_a_domains,Heavy_chain_b_domains,Light_chain_b_domains, Salt_bridges, Labels)



def render(chains_list):
  win  = GraphWin("My Window",350,350)
  print(chains_list)
  for i in range(len(chains_list)):
      for j in range(len(chains_list[i])):
        chains_list[i][j].draw(win)

  win.getMouse()
  win.close()






split_chains = Get_dictionaries(input2)
print(split_chains)
coordinates  = Check_interactions(split_chains)
render(coordinates)
