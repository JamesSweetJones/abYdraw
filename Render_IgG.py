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

def Get_dictionaries(x):
    """
    takes in IgG in SMILES format and identifies variables for dynamically rendering image
    """
    y       = re.sub("\s","",x)
    splitx  = y.split("|")
    chains  = []
    domains = {}
    ADCs    = []
    VHa     = {}
    VHb     = {}
    VLa     = {}
    VLb     = {}
    VHLa    = {}
    VHLb    = {}
    VLLa    = {}
    VLLb    = {}
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


    return(VHa,VLa,VHb,VLb)



def Check_interactions(chains_list):

    VHa       = chains_list[0]
    VLa       = chains_list[1]
    VHb       = chains_list[2]
    VLb       = chains_list[3]
    Heavy_chain_a_domains = []
    Light_chain_a_domains = []
    Heavy_chain_b_domains = []
    Light_chain_b_domains = []
    VHa_VLa   = False
    CH1a_CLa  = False
    VHb_VLb   = False
    CH1b_CLb  = False
    CH2a_CH2b = False
    CH3a_CH3b = False


    VL1a = Polygon(Point(25,25),Point(35,25),Point(45,45),Point(35,45))
    CL1a = Polygon(Point(35,45),Point(45,45),Point(55,65),Point(45,65))
    scFV_VeryL_IgGa = Polygon( Point(0,5),Point(10,5),Point(20,25),Point(15,25))
    #tandem_scFca= Polygon( 10,25,20,25,35,45,20,45))
    #Fab_scFV_Fca= Polygon( 20,45,30,45,40,65,30,65))
    #scFV_VeryLscFVa = Polygon( 30,65,40,65,40,85,40,85))
    #scFV_L_IgGa = Polygon( 15,5,25,5,35,25,25,25))
    #scFV_LscFVa = Polygon( 45,65,55,65,55,85,45,85))


    VH1a = Polygon(Point(40,25),Point(50,25),Point(60,45),Point(50,45))
    CH1a = Polygon( Point(50,45),Point(60,45),Point(70,65),Point(60,65))
    CH2a = Polygon( Point(70,70),Point(80,70),Point(80,90),Point(70,90))
    CH3a = Polygon( Point(70,90),Point(80,90),Point(80,110),Point(70,110))
    scFV_H_IgGa_outer = Polygon( Point(30,5),Point(40,5),Point(50,25),Point(40,25))
    scFV_H_IgGa_inner = Polygon( Point(45,5),Point(55,5),Point(65,25),Point(55,25))
    scFv4_IgG   = Polygon( Point(55,25),Point(65,25),Point(75,45),Point(65,45))
    IgG_2scFV1a  = Polygon( Point(55,110),Point(65,110),Point(65,130),Point(55,130))
    IgG_2scFV2a  = Polygon( Point(55,130),Point(65,130),Point(65,150),Point(55,150))
    IgG_2scFV3a  = Polygon( Point(70,110),Point(80,110),Point(80,130),Point(70,130))
    IgG_2scFV4a  = Polygon( Point(70,130),Point(90,130),Point(90,150),Point(70,150))


    VL1b = Polygon( Point(140,25),Point(130,25),Point(120,45),Point(130,45))
    CL1b = Polygon( Point(130,45),Point(120,45),Point(110,65),Point(120,65))
    #scFV_VeryL_IgGa = Polygon( 165,5,155,5,145,25,165,5))
    #tandem_scFca= Polygon( 155,25, 145,25,135,45,155,45))
    #Fab_scFV_Fca= Polygon( 145,45,135,45,125,65,135,65))
    #scFV_VeryLscFVa = Polygon( 135,65,125,65,125,85,135,85))
    #scFV_L_IgGb = Polygon( 150,5,140,5,130,25,140,25))
    #scFV_LscFVb = Polygon( 120,60,110,60,110,80,120,80))


    VH1b = Polygon( Point(125,25),Point(115,25),Point(105,45),Point(115,45))
    CH1b= Polygon( Point(115,45),Point(105,45),Point(95,65),Point(105,65))
    CH2b= Polygon( Point(95,70), Point(85,70),Point(85,90),Point(95,90))
    CH3b= Polygon( Point(95,90), Point(85,90),Point(85,110),Point(95,110))
    scFv_H_IgGb_outer = Polygon(Point(135,5),Point(125,5),Point(115,25),Point(125,25))
    scFv_H_IgGb_inner= Polygon( Point(120,5),Point(110,5),Point(100,25),Point(110,25))
    scFv4_IgGb   = Polygon( Point(110,25),Point(100,25),Point(90,45),Point(100,45))
    IgG_2scFV1b  = Polygon( Point(110,110),Point(100,110),Point(100,130),Point(110,130))
    IgG_2scFV2b  = Polygon( Point(110,130),Point(100,130),Point(100,150),Point(110,150))
    IgG_2scFV3b  = Polygon( Point(95,110),Point(85,110),Point(85,130),Point(95,130))
    IgG_2scFV4b  = Polygon( Point(95,130),Point(85,130),Point(85,150),Point(95,150))




#########Heavy chain A check what domains are there################
    for i in range(len(VHa)):
        keyslist = list(VHa.keys())
        if i==0:
            if keyslist[i] == "VH":
                Heavy_chain_a_domains.append(VH1a)


        if VHa.get(keyslist[i])[0] == (int(VHa.get(keyslist[i-1])[0])+1) and VHa.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
            if keyslist[i] == "CH1":
                Heavy_chain_a_domains.append(CH1a)


            elif keyslist[i] == "CH2":
                Heavy_chain_a_domains.append(CH2a)
                if VHa.get('CH2')[0] == VHb.get('CH2')[1] and VHa.get('CH2')[1] == VHb.get('CH2')[0]:
                    CH2a_CH2b = True

            elif keyslist[i] == "CH3":
                Heavy_chain_a_domains.append(CH3a)
                if VHa.get('CH3')[0] == VHb.get('CH3')[1] and VHa.get('CH3')[1] == VHb.get('CH3')[0]:
                    CH3a_CH3b = True
            #checking for extra chains at the end
            elif keyslist[i] == "VH2" or "CH4" and keyslist[i-1] == "CH3":
                Heavy_chain_a_domains.append(IgG_2scFV1a)
            elif keyslist[i] == "VH2" and keyslist[i-1] == "CH4":
                Heavy_chain_a_domains.append(IgG_2scFV2a)



        elif VHa.get(keyslist[i])[0] == VHa.get(keyslist[i-1])[1] and  VHa.get(keyslist[i])[1] == VHa.get(keyslist[i-1])[0]:
            if keyslist[i] == "VL2":
                if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH3" :
                    Heavy_chain_a_domains.append(IgG_2scFV3a)
                elif keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                    Heavy_chain_a_domains.append(IgG_2scFV4a)
        elif VHa.get(keyslist[i])[0] != VHa.get(keyslist[i-1])[1] and  VHa.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
            if keyslist[i] == "VL2":
                if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                    Heavy_chain_a_domains.append(IgG_2scFV3a)


#########Light chain A check what domains are there################
    for i in range(len(VLa)):
        keyslist = list(VLa.keys())
        if i==0:
            if keyslist[i] == "VL":
                Light_chain_a_domains.append(VL1a)

        if VLa.get(keyslist[i])[0] == (int(VLa.get(keyslist[i-1])[0])+1):
            if keyslist[i] == "CL":
                Light_chain_a_domains.append(CL1a)



#########Heavy chain B check what domains are there################

    for i in range(len(VHb)):
        keyslist = list(VHb.keys())
        if i==0:
            if keyslist[i] == "VH":
                Heavy_chain_b_domains.append(VH1b)

        if VHb.get(keyslist[i])[0] == (int(VHb.get(keyslist[i-1])[0])+1) and VHb.get(keyslist[i])[1] != VHb.get(keyslist[i-1])[0]:
            if keyslist[i] == "CH1":
                Heavy_chain_b_domains.append(CH1b)


            elif keyslist[i] == "CH2":
                Heavy_chain_b_domains.append(CH2b)
                if VHb.get('CH2')[0] == VHa.get('CH2')[1] and VHb.get('CH2')[1] == VHa.get('CH2')[0]:
                    CH2a_CH2b = True

            elif keyslist[i] == "CH3":
                Heavy_chain_b_domains.append(CH3b)
                if VHb.get('CH3')[0] == VHa.get('CH3')[1] and VHb.get('CH3')[1] == VHa.get('CH3')[0]:
                    CH3a_CH3b = True
            #checking for extra chains at the end
            elif keyslist[i] == "VH2" or "CH4" and keyslist[i-1] == "CH3":
                Heavy_chain_b_domains.append(IgG_2scFV1b)
            elif keyslist[i] == "VH2" and keyslist[i-1] == "CH4":
                Heavy_chain_b_domains.append(IgG_2scFV2b)




        elif VHb.get(keyslist[i])[0] == VHb.get(keyslist[i-1])[1] and  VHb.get(keyslist[i])[1] == VHb.get(keyslist[i-1])[0]:
            if keyslist[i] == "VL2":
                if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH3" :
                    Heavy_chain_b_domains.append(IgG_2scFV3b)
                elif keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                    Heavy_chain_b_domains.append(IgG_2scFV4b)
        elif VHb.get(keyslist[i])[0] != VHb.get(keyslist[i-1])[1] and  VHb.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
            if keyslist[i] == "VL2":
                if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                    Heavy_chain_b_domains.append(IgG_2scFV3b)


#########Light chain B check what domains are there################
    for i in range(len(VLb)):
        keyslist = list(VLb.keys())
        if i==0:
            if keyslist[i] == "VL":
                Light_chain_b_domains.append(VL1b)

        if VLb.get(keyslist[i])[0] == (int(VLb.get(keyslist[i-1])[0])+1):
            if keyslist[i] == "CL":
                Light_chain_b_domains.append(CL1b)

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




    standard_interactions = [VHa_VLa,CH1a_CLa,VHb_VLb,CH1b_CLb,CH2a_CH2b,CH3a_CH3b]
    return (Heavy_chain_a_domains,Light_chain_a_domains,Heavy_chain_b_domains,Light_chain_b_domains)



def render(chains_list):
  win  = GraphWin("My Window",250,250)

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

split_chains = Get_dictionaries(input2)
coordinates  = Check_interactions(split_chains)
render(coordinates)
