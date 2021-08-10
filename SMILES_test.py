#!/usr/bin/python


"""
NEED TO DO
>Heavy chains pre VH domain
>Light chain adaptations
"""
import re
import sys

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
    VHa_VLa   = False
    CH1a_CLa  = False
    VHb_VLb   = False
    CH1b_CLb  = False
    CH2a_CH2b = False
    CH3a_CH3b = False

    #Heavy chaina
    VH1a =0
    CH1a=0
    CH2a=0
    CH3a=0
    CH4a=0
    scFV_H_IgGa_outer=0
    scFV_H_IgGa_inner=0
    scFv4_IgGa=0
    IgG_2scFV1a=0
    IgG_2scFV2a=0
    IgG_2scFV3a=0
    IgG_2scFV4a=0

    #Light chaina
    VL1a=0
    CL1a=0
    CL2a=0

    #Heavy
    VH1b=0
    CH1b=0
    CH2b=0
    CH3b=0
    CH4b=0
    scFV_H_IgGb_outer=0
    scFV_H_IgGb_inner=0
    scFv4_IgGb=0
    IgG_2scFV1b=0
    IgG_2scFV2b=0
    IgG_2scFV3b=0
    IgG_2scFV4b=0

    #Light chainsb
    VL1b=0
    CL1b=0
    CL2b=0

    

#########Heavy chain A check what domains are there################
    for i in range(len(VHa)):
        keyslist = list(VHa.keys())
        if i==0:
            if keyslist[i] == "VH":
                VH1a = int(VHa.get('VH')[0])


        if VHa.get(keyslist[i])[0] == (int(VHa.get(keyslist[i-1])[0])+1) and VHa.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
            if keyslist[i] == "CH1":
                CH1a = int(VHa.get('CH1')[0])


            elif keyslist[i] == "CH2":
                CH2a = int(VHa.get('CH2')[0])
                if VHa.get('CH2')[0] == VHb.get('CH2')[1] and VHa.get('CH2')[1] == VHb.get('CH2')[0]:
                    CH2a_CH2b = True

            elif keyslist[i] == "CH3":
                CH3a = int(VHa.get('CH3')[0])
                if VHa.get('CH3')[0] == VHb.get('CH3')[1] and VHa.get('CH3')[1] == VHb.get('CH3')[0]:
                    CH3a_CH3b = True
            #checking for extra chains at the end
            elif keyslist[i] == "VH2" or "CH4" and keyslist[i-1] == "CH3":
                IgG_2scFV1a = int(VHa.get(keyslist[i])[0])
            elif keyslist[i] == "VH2" and keyslist[i-1] == "CH4":
                IgG_2scFV2a = int(VHa.get(keyslist[i])[0])



        elif VHa.get(keyslist[i])[0] == VHa.get(keyslist[i-1])[1] and  VHa.get(keyslist[i])[1] == VHa.get(keyslist[i-1])[0]:
            if keyslist[i] == "VL2":
                if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH3" :
                    IgG_2scFV3a = int(VHa.get(keyslist[i])[0])
                elif keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                    IgG_2scFV4a = int(VHa.get(keyslist[i])[0])
        elif VHa.get(keyslist[i])[0] != VHa.get(keyslist[i-1])[1] and  VHa.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
            if keyslist[i] == "VL2":
                if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                    IgG_2scFV3a = int(VHa.get(keyslist[i])[0])


#########Light chain A check what domains are there################
    for i in range(len(VLa)):
        keyslist = list(VLa.keys())
        if i==0:
            if keyslist[i] == "VL":
                VL1a = int(VLa.get('VL')[0])

        if VLa.get(keyslist[i])[0] == (int(VLa.get(keyslist[i-1])[0])+1):
            if keyslist[i] == "CL":
                CL1a = int(VLa.get('CL')[0])



#########Heavy chain B check what domains are there################

    for i in range(len(VHb)):
        keyslist = list(VHb.keys())
        if i==0:
            if keyslist[i] == "VH":
                VH1b = int(VHb.get('VH')[0])

        if VHb.get(keyslist[i])[0] == (int(VHb.get(keyslist[i-1])[0])+1) and VHb.get(keyslist[i])[1] != VHb.get(keyslist[i-1])[0]:
            if keyslist[i] == "CH1":
                CH1b = int(VHb.get('CH1')[0])


            elif keyslist[i] == "CH2":
                CH2b = int(VHb.get('CH2')[0])
                if VHb.get('CH2')[0] == VHa.get('CH2')[1] and VHb.get('CH2')[1] == VHa.get('CH2')[0]:
                    CH2a_CH2b = True

            elif keyslist[i] == "CH3":
                CH3b = int(VHb.get('CH3')[0])
                if VHb.get('CH3')[0] == VHa.get('CH3')[1] and VHb.get('CH3')[1] == VHa.get('CH3')[0]:
                    CH3a_CH3b = True
            #checking for extra chains at the end
            elif keyslist[i] == "VH2" or "CH4" and keyslist[i-1] == "CH3":
                IgG_2scFV1b = int(VHb.get(keyslist[i])[0])
            elif keyslist[i] == "VH2" and keyslist[i-1] == "CH4":
                IgG_2scFV2b = int(VHb.get(keyslist[i])[0])




        elif VHb.get(keyslist[i])[0] == VHb.get(keyslist[i-1])[1] and  VHb.get(keyslist[i])[1] == VHb.get(keyslist[i-1])[0]:
            if keyslist[i] == "VL2":
                if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH3" :
                    IgG_2scFV3b = int(VHb.get(keyslist[i])[0])
                elif keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                    IgG_2scFV4b = int(VHb.get(keyslist[i])[0])
        elif VHb.get(keyslist[i])[0] != VHb.get(keyslist[i-1])[1] and  VHb.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
            if keyslist[i] == "VL2":
                if keyslist[i-1] == "VH2" and keyslist[i-2] == "CH4":
                    IgG_2scFV3b = int(VHb.get(keyslist[i])[0])


#########Light chain B check what domains are there################
    for i in range(len(VLb)):
        keyslist = list(VLb.keys())
        if i==0:
            if keyslist[i] == "VL":
                VL1b = int(VLb.get('VL')[0])

        if VLb.get(keyslist[i])[0] == (int(VLb.get(keyslist[i-1])[0])+1):
            if keyslist[i] == "CL":
                CL1b = int(VLb.get('CL')[0])

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
    Heavy_chain_a_domains = [VH1a,CH1a,CH2a,CH3a,scFV_H_IgGa_outer, scFV_H_IgGa_inner, scFv4_IgGa,IgG_2scFV1a,IgG_2scFV2a,IgG_2scFV3a,IgG_2scFV4a]
    Light_chain_a_domains = [VL1a,CL1a,CL2a]
    Heavy_chain_b_domains = [VH1b,CH1b,CH2b,CH3b,scFV_H_IgGb_outer, scFV_H_IgGb_inner, scFv4_IgGb,IgG_2scFV1b,IgG_2scFV2b,IgG_2scFV3b,IgG_2scFV4b]
    Light_chain_b_domains = [VL1b,CL1b,CL2b]

    for i in range(len(Heavy_chain_b_domains)):
        if Heavy_chain_b_domains[i] > 0:
            print(Heavy_chain_b_domains[i])
    for i in range(len(Light_chain_b_domains)):
        if Light_chain_b_domains[i] > 0:
            print(Light_chain_b_domains[i])
    for i in range(len(Heavy_chain_a_domains)):
        if Heavy_chain_a_domains[i] > 0:
            print(Heavy_chain_a_domains[i])
    for i in range(len(Light_chain_a_domains)):
        if Light_chain_a_domains[i] > 0:
            print(Light_chain_a_domains[i])





    #for i in range(len(VHa)):
    #    keyslist = list(VHa.keys())
    #    if i==0:
    #        print("yes")
    #        print(VHa.get(keyslist[i]))
    #    if VHa.get(keyslist[i])[0] == (int(VHa.get(keyslist[i-1])[0])+1) and VHa.get(keyslist[i])[1] != VHa.get(keyslist[i-1])[0]:
    #        print("yes")
    #        print(VHa.get(keyslist[i]))
    #    elif VHa.get(keyslist[i])[0] == VHa.get(keyslist[i-1])[1] and  VHa.get(keyslist[i])[1] == VHa.get(keyslist[i-1])[0]:
    #        print("no")
    #        print(VHa.get(keyslist[i]))






#    ADCs                    = re.findall("X\[TYPE\: (.*?)]-(.*?)\)", str(splitx))
#    Notch_mod               = re.findall("CH[1-3]@\(.*?\)", str(splitx))
#    Notch_mod_domain        = re.match("\([1-9]|[1-9]:", Notch_mod)
#    Notch_mod_domain        = int(re.sub("\(|:", "", Notch_mod_domain))
#    salt_bridges_location   = re.findall("H\(.*?\){[1-9]}", str(splitx))
#    Number_of_salt_bridges  = re.match("{[1-9]}",salt_bridges_location)
#    Number_of_salt_bridges  = int(re.sub("\{|\}","",Number_of_salt_bridges))

split_chains = Get_dictionaries(input1)
#print(split_chains)
Check_interactions(split_chains)
