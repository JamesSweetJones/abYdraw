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

input1="VH.b(1:6)-CH1(2:7){1}-H(3:10){2}-CH2(4:11) - CH3(5:12) | VL.b(6:1)-CL(7:2){1} | VH.a(8:13)- CH1(9:14){1}-H(10:3){2}-CH2(11:4)-CH3(12:5) | VL.a(13:8)-CL(14:9){1}"
input2="VH.b(1:8)-CH1(2:9){1}-H(3:12){2}-CH2(4:13)-CH3(5:14)-L-VH2(6:7)-L-VL2(7:6)|VL.b(8:1)-CL(9:2){1}|VH.a(10:17)-CH1(11:18){1}-H(12:3){2}-CH2(13:4)-CH3(14:5)-L-VH2(15:16)-L-VL2(16:15)|VL.a(17:10)-CL(18:11){1}"
input3="VL.a(1:2)|VH.a(2:1)|VL.b(3:4)|VH.b(4:3)"
input4="VH.b(1:2)-L-VL.b(2:1)-H*(3:6){1}[MOD: engineered disulphide bond]-X[TYPE:LEUCINE] | VH.a(4:5)-L-VL.a(5:4)-H*(6:3){1}[MOD: engineered disulphide bond]-X[TYPE:LEUCINE]"
def Get_dictionaries(x):
    """
    takes in IgG in SMILES format and identifies variables for dynamically rendering image
    """
    y       = re.sub("\s","",x)
    splitx  = y.split("|")
    chains  = []
    ADCs    = []
    VHa     = {}
    VHb     = {}
    VLa     = {}
    VLb     = {}
    Salt_bridges = ""
    VHa_VLa_bonds= ""
    VHb_VLb_bonds= ""
    CH1a_CL1a_bonds=""
    CH1b_CL1b_bonds=""
    #get chains
    for i in splitx:
        chain = i.split("(")[0]
        chains.append(chain)
    for i in range(len(chains)):
        dict = str(re.sub("\.","",str(chains[i])))
        chain = splitx[i].split("-")
        #print(dict,chain)
        for j in range(len(chain)):
            #print(chain[j])
            domain   =  re.findall("^CH[1-9]|^VL[1-9]|^VH[1-9]|^CL|^VH|^VL|^H", str(chain[j]))
            domain   =  str(re.sub("\[|\'|\]","", str(domain)))
            ##Get Bond numbers
            if domain == "H" and Salt_bridges == "":
                Salt_bridges = re.findall("\{.*?\}", str(chain[j]))
                Salt_bridges = str(Salt_bridges)
                Salt_bridges = int(re.sub("\{|\'|\}|\[|\]","", Salt_bridges))
            elif domain != "VH":
                if dict == "VHa" and VHa_VLa_bonds == "":
                    VHa_VLa_bonds_match = re.findall("\{.*?\}", str(chain[j]))
                    if VHa_VLa_bonds_match != []:
                        VHa_VLa_bonds = str(VHa_VLa_bonds_match)
                        VHa_VLa_bonds = int(re.sub("\{|\'|\}|\[|\]","", VHa_VLa_bonds))
                elif dict == "VHb" and VHb_VLb_bonds == "":
                    VHb_VLb_bonds_match = re.findall("\{.*?\}", str(chain[j]))
                    if VHb_VLb_bonds_match != []:
                        print(VHb_VLb_bonds_match)
                        VHb_VLb_bonds = str(VHb_VLb_bonds_match)
                        VHb_VLb_bonds = int(re.sub("\{|\'|\}|\[|\]","", VHb_VLb_bonds))
            elif domain != "CH1":
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

            location = []
            locationstr =  re.findall("\((.*?)\)", str(chain[j]))
            locationstr =  str(re.sub("\[|\'|\]","", str(locationstr)))
            locationstr =  str(re.sub(":",",", str(locationstr)))
            locationstr = locationstr.split(",")
            if domain != "":
                for x in range(len(locationstr)):
                    location.append(int(locationstr[x]))
            if dict == "VHa" and domain !="":
                VHa[domain] = location
            elif dict == "VHb" and domain !="":
                VHb[domain] = location
            elif dict == "VLa" and domain !="":
                VLa[domain] = location
            elif dict == "VLb" and domain !="":
                VLb[domain] = location
            else:
                continue
    print(CH1a_CL1a_bonds)
    return(VHa,VLa,VHb,VLb,Salt_bridges,VHa_VLa_bonds,VHb_VLb_bonds,CH1a_CL1a_bonds,CH1b_CL1b_bonds)



def Check_interactions(chains_list):

    VHa       = chains_list[0]
    VLa       = chains_list[1]
    VHb       = chains_list[2]
    VLb       = chains_list[3]
    Salt_bridge_count     = chains_list[4]
    VHa_VLa_bond_count    = chains_list[5]
    VHb_VLb_bond_count    = chains_list[6]
    CH1a_CL1a_bond_count  = chains_list[7]
    CH1b_CL1b_bond_count  = chains_list[8]
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


    VL1a = Polygon(Point(50,50),Point(70,50),Point(90,90),Point(70,90))
    VL1a_label_location = Point(40,70)
    CL1a = Polygon(Point(70,90),Point(90,90),Point(110,130),Point(90,130))
    CL1a_label_location = Point(60,110)
    scFV_VeryL_IgGa = Polygon( Point(0,10),Point(20,10),Point(40,50),Point(30,50))
    #tandem_scFca= Polygon( 10,25,20,25,35,45,20,45))
    #Fab_scFV_Fca= Polygon( 20,45,30,45,40,65,30,65))
    #scFV_VeryLscFVa = Polygon( 30,65,40,65,40,85,40,85))
    #scFV_L_IgGa = Polygon( 15,5,25,5,35,25,25,25))
    #scFV_LscFVa = Polygon( 45,65,55,65,55,85,45,85))


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


    VL1b = Polygon( Point(280,50),Point(260,50),Point(240,90),Point(260,90))
    VL1b_label_location = Point(280,70)
    CL1b = Polygon( Point(260,90),Point(240,90),Point(220,130),Point(240,130))
    CL1b_label_location = Point(260,110)
    #scFV_VeryL_IgGa = Polygon( 165,5,155,5,145,25,165,5))
    #tandem_scFca= Polygon( 155,25, 145,25,135,45,155,45))
    #Fab_scFV_Fca= Polygon( 145,45,135,45,125,65,135,65))
    #scFV_VeryLscFVa = Polygon( 135,65,125,65,125,85,135,85))
    #scFV_L_IgGb = Polygon( 150,5,140,5,130,25,140,25))
    #scFV_LscFVb = Polygon( 120,60,110,60,110,80,120,80))


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

    VHa_VLa_bond_location   = Line(Point(90,90),Point(100,90))
    VHb_VLb_bond_location   = Line(Point(230,90),Point(210,90))
    CH1a_CL1a_bond_location = Line(Point(110,130),Point(120,130))
    CH1b_CL1b_bond_location = Line(Point(230,130),Point(210,130))
    if VHa_VLa_bond_count != "":
        Salt_bridges.append(VHa_VLa_bond_location)
    if VHb_VLb_bond_count!= "":
        Salt_bridges.append(VHb_VLb_bond_location)
    if CH1a_CL1a_bond_count != "":
        Salt_bridges.append(CH1a_CL1a_bond_location)
    if CH1b_CL1b_bond_count !="":
        Salt_bridges.append(CH1b_CL1b_bond_location)
#########Heavy chain A check what domains are there################
    for i in range(len(VHa)):
        keyslist = list(VHa.keys())
        if i==0:
            if keyslist[i] == "VH":
                Heavy_chain_a_domains.append(VH1a)
                location = VHa.get(keyslist[i])[0]
                Labels.append(Text(VH1a_label_location,str(location)))



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
            if keyslist[i] == "VL":
                Light_chain_a_domains.append(VL1a)
                location = VLa.get(keyslist[i])[0]
                Labels.append(Text(VL1a_label_location,str(location)))

        if VLa.get(keyslist[i])[0] == (int(VLa.get(keyslist[i-1])[0])+1):
            if keyslist[i] == "CL":
                Light_chain_a_domains.append(CL1a)
                location = VLa.get(keyslist[i])[0]
                Labels.append(Text(CL1a_label_location,str(location)))



#########Heavy chain B check what domains are there################

    for i in range(len(VHb)):
        keyslist = list(VHb.keys())
        if i==0:
            if keyslist[i] == "VH":
                Heavy_chain_b_domains.append(VH1b)
                location = VHb.get(keyslist[i])[0]
                Labels.append(Text(VH1b_label_location,str(location)))

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
                    CH2a_CH2b = True
                    location = VHb.get(keyslist[i])[0]
                    Labels.append(Text(CH2b_label_location,str(location)))

            elif keyslist[i] == "CH3":
                Heavy_chain_b_domains.append(CH3b)
                if VHb.get('CH3')[0] == VHa.get('CH3')[1] and VHb.get('CH3')[1] == VHa.get('CH3')[0]:
                    CH3a_CH3b = True
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
                    Labels.append(Text(IgG_2scFV4b_label_location,str(location)))
        elif VHb.get(keyslist[i])[0] != VHb.get(keyslist[i-1])[1] and  VHb.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
            if keyslist[i] == "VL2":
                if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                    Heavy_chain_b_domains.append(IgG_2scFV3b)
                    location = VHb.get(keyslist[i])[0]
                    Labels.append(Text(IgG_2scFV3b_label_location,str(location)))


#########Light chain B check what domains are there################
    for i in range(len(VLb)):
        keyslist = list(VLb.keys())
        if i==0:
            if keyslist[i] == "VL":
                Light_chain_b_domains.append(VL1b)
                location = VLb.get(keyslist[i])[0]
                Labels.append(Text(VL1b_label_location,str(location)))

        if VLb.get(keyslist[i])[0] == (int(VLb.get(keyslist[i-1])[0])+1):
            if keyslist[i] == "CL":
                Light_chain_b_domains.append(CL1b)
                location = VLb.get(keyslist[i])[0]
                Labels.append(Text(CL1b_label_location,str(location)))
    #Light chain a
    #if VHa.get('VH')[0] == VLa.get('VL')[1] and VHa.get('VH')[1] == VLa.get('VL')[0]:
    #    VHa_VLa = True
    #    VL1a = int(VLa.get('VL')[0])
    #if VHa.get('CH1')[0] == VLa.get('CL')[1] and VHa.get('CH1')[1] == VLa.get('CL')[0]:
    #    CH1a_CLa = True
    #    CL1a = int(VLa.get('CL')[0])

    #check light-heavy chain interactions


    #if VHb.get('VH')[0] == VLb.get('VL')[1] and VHb.get('VH')[1] == VLb.get('VL')[0]:
    #    VHb_VLb = True
    #if VHb.get('CH1')[0] == VLb.get('CL')[1] and VHb.get('CH1')[1] == VLb.get('CL')[0]:
    #    CH1b_CLb = True




    standard_interactions = [VHa_VLa,CH1a_CLa,VHb_VLb,CH1b_CLb,CH2a_CH2b,CH3a_CH3b,Salt_bridges]
    print(Labels)
    return (Heavy_chain_a_domains,Light_chain_a_domains,Heavy_chain_b_domains,Light_chain_b_domains, Salt_bridges, Labels)



def render(chains_list):
  win  = GraphWin("My Window",300,300)

  for i in range(len(chains_list)):
      for j in range(len(chains_list[i])):
        chains_list[i][j].draw(win)

  win.getMouse()
  win.close()





#    ADCs                    = re.findall("X\[TYPE\: (.*?)]-(.*?)\)", str(splitx))
#    Notch_mod               = re.findall("CH[1-3]@\(.*?\)", str(splitx))
#    Notch_mod_domain        = re.match("\([1-9]|[1-9]:", Notch_mod)
#    Notch_mod_domain        = int(re.sub("\(|:", "", Notch_mod_domain))
#    salt_bridges_location   = re.findall("H\(.*?\){[1-9]}", str(splitx))
#    Number_of_salt_bridges  = re.match("{[1-9]}",salt_bridges_location)
#    Number_of_salt_bridges  = int(re.sub("\{|\}","",Number_of_salt_bridges))

split_chains = Get_dictionaries(input1)
coordinates  = Check_interactions(split_chains)
render(coordinates)
