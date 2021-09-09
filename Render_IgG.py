#!/usr/bin/python
import re
import sys
import os
import tkinter as tk

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
        dict = str(re.sub("\.|\+|\_","",str(chains[i])))
        chain = splitx[i].split("-")
        for j in range(len(chain)):
            domain   =  re.findall("^CH[0-9][@+>_!*]\(.*\)\[.*\]|^CH[0-9]\(.*\)\[.*\]|^CH[0-9][@+>_!*]|^CL[0-9][@+>_!*]\(.*\)\[.*\]|^CL[0-9]\(.*\)\[.*\]|^CL[0-9][@+>_!*]|^CL[@+>_!*]|^CH[0-9]|^VL.[ab]|^VH.[ab]|^VL[1-9]|^VH[1-9]|VH[+_*]\.[ab]|VL[+_*]\.[ab]|^CL[0-9]|^CL|^VH|^VL|^H|^X", str(chain[j]))
            for i in range(len(domain)):
                if "*" in domain[i]:
                    domain = str(re.sub("\(.*\)\[MOD\:","",str(domain)))
                elif ("VH+") in domain[i] or ("VL+") in domain[i]:
                    domain = str(re.sub("\+","",domain[i]))
                    domain = domain+"+"
                elif ("VH_") in domain[i] or ("VL_") in domain[i]:
                    domain = str(re.sub("\_","",domain[i]))
                    domain = domain+"_"
            else:
                domain = str(re.sub("\[|\'|\]|\.","", str(domain)))


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

    VHa_chain       = chains_list[0]
    VLa_chain       = chains_list[1]
    VHb_chain       = chains_list[2]
    VLb_chain       = chains_list[3]
    print(chains_list[0])
    print(chains_list[1])
    print(chains_list[2])
    print(chains_list[3])
    Salt_bridge_count     = chains_list[4]
    VHa_VLa_bond_count    = chains_list[5]
    VHb_VLb_bond_count    = chains_list[6]
    CH1a_CLa_bond_count  = chains_list[7]
    CH1b_CLb_bond_count  = chains_list[8]
    fragment              = chains_list[9]

    if fragment != {}:
        VLa_chain = fragment

    print(fragment)



#########Make bonds################
    def bondmaker(x,y,side):


        if re.search("[a-b]",x) == None and re.search("[a-b]",y) == None:
            x_bond      = str(x)+str(side)+"_bottom_bond_location"
            y_bond      = str(y)+str(side)+   "_top_bond_location"
        elif re.search("[a-b]",x) != "" and re.search("[a-b]",y) == None:
            x_bond      = str(x)+"_bottom_bond_location"
            y_bond      = str(y)+str(side)+   "_top_bond_location"
        elif re.search("[a-b]",x) == None and re.search("[a-b]",y) != "" :
            x_bond      = str(x)+str(side)+"_bottom_bond_location"
            y_bond      = str(y)          +   "_top_bond_location"
        elif re.search("[a-b]",x) != "" and re.search("[a-b]",y) != "":
            x_bond      = str(x)+"_bottom_bond_location"
            y_bond      = str(y)+   "_top_bond_location"

        return(x_bond,y_bond)




#######Heavy chain A check what domains are there################
    def GetChainCoordinates(dictionary,chain):
#########VLa Chain coordinates################


        VLa                       = [97,100,117,100,129,120,147,120,177,174,137,174]
        VLa_no_V                  = [100,100,140,110,180,174,140,174]
        VLa_label_location        = [80,140]
        VLa_top_bond_location     = [120,100]
        VLa_bottom_bond_location  = [157,174]


        CLa                       = [140,180,180,180,220,254,180,254]
        CLa_notch                 = [140,180,180,180,220,220,220,254,180,254]
        CLa_inward_notch          = [140,180,180,180,180,220,220,254,180,254]
        CLa_label_location        = [120,220]
        CLa_top_bond_location     = [160,180]
        CLa_bottom_bond_location  = [200,254]

        scFV_VeryL_IgGa = [0,20,40,20,80,100,60,100]
        tandem_scFca    = [20,50,40,50,70,90,40,90]
        Fab_scFV_Fca    = [40,90,60,90,80,130,60,130]
        scFV_VeryLscFVa = [60,130,80,130,80,170,80,90]
        scFV_L_IgGa     = [30,10,50,10,70,50,50,50]
        scFV_LscFVa     = [90,130,110,130,110,170,90,170]

    #########VHa Chain coordinates################

        VHa                      = [167,120,187,120,177,100,197,100,237,174,197,174]
        VHa_no_V                 = [160,100,200,100,240,174,200,174]
        VHa_label_location       = [240,140]
        VHa_top_bond_location    = [180,100]
        VHa_bottom_bond_location = [217,174]

        CH1a                      = [200,180,240,180,280,254,240,254]
        CH1a_notch                = [200,180,240,180,280,254,240,254,200,220]
        CH1a_inward_notch         = [200,180,240,180,280,254,240,254,240,220]
        CH1a_label_location       = [280,220]
        CH1a_top_bond_location    = [220,180]
        CH1a_bottom_bond_location = [260,254]

        CH2a                      = [280,280,320,280,320,354,280,354]
        CH2a_notch                = [280,280,320,280,320,340,320,354,280,354]
        CH2a_inward_notch         = [280,280,320,280,320,300,320,354,280,354]
        CH2a_label_location       = [260,320]
        CH2a_top_bond_location    = [300,280]
        CH2a_bottom_bond_location = [300,354]

        CH3a                      = [280,360,320,360,320,434,280,434]
        CH3a_notch                = [280,360,320,360,340,380,320,434,280,434]
        CH3a_inward_notch         = [280,360,320,360,300,400,320,434,280,434]
        CH3a_label_location       = [260,400]
        CH3a_top_bond_location    = [300,360]
        CH3a_bottom_bond_location = [300,434]


        scFV_H_IgGa_outer = [120,20,160,20,200,100,160,100]
        scFV_H_IgGa_inner = [180,20,220,20,260,100,220,100]
        scFv4_IgG         = [220,100,260,100,300,180,260,180]

        VH2a                      = [300,440,320,440,320,514,280,514,280,460,300,440]
        CH4a                      = [280,440,320,440,320,514,280,514]
        VH2a_label_location       = CH4a_label_location       = [300,480]
        VH2a_top_bond_location    = CH4a_top_bond_location    = [300,440]
        VH2a_bottom_bond_location = CH4a_bottom_bond_location = [300,514]

        VL2a                      = [220,440,240,440,260,460,260,514,220,514]
        CL2a                      = [220,440,260,440,260,514,220,514]
        VL2a_label_location       = CL2a_label_location       = [240,480]
        VL2a_top_bond_location    = CL2a_top_bond_location    = [240,440]
        VL2a_bottom_bond_location = CL2a_bottom_bond_location = [240,514]


        VH2_2a                      = [280,520,320,520,320,594,280,594]
        VH2_2a_label_location       = CH5a_label_location       = [300,560]
        VH2_2a_top_bond_location    = CH5a_top_bond_location    = [300,540]
        VH2_2a_bottom_bond_location = CH5a_bottom_bond_location = [300,694]


        VL2_2a                      = [220,520,260,520,260,594,220,594]
        VL2_2a_label_location       = [240,560]
        VL2_2a_top_bond_location    = [240,520]
        VL2_2a_bottom_bond_location = [240,594]



    #########VLb Chain coordinates################

        VLb                      = [563,100,543,100,533,120,513,120,483,174,523,174]
        VLb_no_V                 = [560,100,520,50,480,87,520,87]
        VLb_label_location       = [560,140]
        VLb_top_bond_location    = [540,100]
        VLb_bottom_bond_location = [503,174]


        CLb                      = [520,180,480,180,440,254,480,254]
        CLb_notch                = [520,180,480,180,440,220,440,254,480,254]
        CLb_inward_notch         = [520,180,480,180,480,220,440,254,480,254]
        CLb_label_location       = [520,220]
        CLb_top_bond_location    = [500,180]
        CLb_bottom_bond_location = [460,254]

        scFV_VeryL_IgGa = [330,10,310,10,290,50,330,10]
        tandem_scFca= [330,50,290,50,270,90,330,90]
        Fab_scFV_Fca= [290,90,270,90,250,130,270,130]
        scFV_VeryLscFVa = [270,130,250,130,250,170,270,170]
        scFV_L_IgGb = [300,10,280,10,260,50,280,50]
        scFV_LscFVb = [240,120,220,120,220,160,240,160]

    #########VHb Chain coordinates################

        VHb                      = [463,100,483,100,473,120,493,120,463,174,423,174]
        VHb_no_v                 = [500,100,460,100,420,174,460,174]
        VHb_label_location       = [420,140]
        VHb_top_bond_location    = [480,100]
        VHb_bottom_bond_location = [443,174]

        CH1b                      = [460,180,420,180,380,254,420,254]
        CH1b_notch                = [460,180,420,180,380,254,420,254,460,220]
        CH1b_inward_notch         = [460,180,420,180,380,254,420,254,420,220]
        CH1b_label_location       = [380,220]
        CH1b_top_bond_location    = [440,180]
        CH1b_bottom_bond_location = [400,254]


        CH2b                      = [380,280,340,280,340,354,380,354]
        CH2b_notch                = [380,280,340,280,320,320,340,354,380,354]
        CH2b_inward_notch         = [380,280,340,280,320,360,340,354,380,354]
        CH2b_label_location       = [400,320]
        CH2b_top_bond_location    = [360,280]
        CH2b_bottom_bond_location = [360,354]


        CH3b                      = [380,360,340,360,340,434,380,434]
        CH3b_notch                = [380,360,340,360,320,400,340,434,380,434]
        CH3b_inward_notch         = [380,360,340,360,360,400,340,434,380,434]
        CH3b_label_location       = [400,400]
        CH3b_top_bond_location    = [360,360]
        CH3b_bottom_bond_location = [360,434]


        scFv_H_IgGb_outer                = [540,20,500,20,260,100,500,100]
        scFv_H_IgGb_inner                = [280,20,440,20,400,100,440,100]
        scFv4_IgGb                       = [440,100,400,100,360,180,400,180]

        VH2b                      = [380,460,360,440,340,440,340,514,380,514]
        CH4b                      = [380,440,340,440,340,514,380,514]
        VH2b_label_location       = CH4b_label_location       = [360,480]
        VH2b_top_bond_location    = CH4b_top_bond_location    = [360,440]
        VH2b_bottom_bond_location = CH4b_bottom_bond_location = [360,514]

        VL2b                      = [440,440,420,440,400,460,400,514,440,514]
        CL2b                      = [440,440,400,440,400,514,440,514]
        VL2b_label_location       = CL2b_label_location       = [420,480]
        VL2b_top_bond_location    = CL2b_top_bond_location    = [420,440]
        VL2b_bottom_bond_location = CL2b_bottom_bond_location = [420,514]


        VH2_2b                    = [380,520,340,520,340,594,380,594]
        CH5b                      = [380,520,340,520,340,594,380,594]
        CH5b_label_location       = VH2_2b_label_location       = [360,560]
        CH5b_top_bond_location    = VH2_2b_top_bond_location    = [360,520]
        CH5b_bottom_bond_location = VH2_2b_bottom_bond_location = [360,594]


        VL2_2b                      = [440,520,400,520,400,594,440,594]
        VL2_2b_label_location       = [420,560]
        VL2_2b_top_bond_location    = [420,520]
        VL2_2b_bottom_bond_location = [420,594]




    #########Check VH/VL ADC Coordinates################
        X_VLa                    = [100,60,80,40,100,200,120,40,100,60]
        X_VLa_label_location     = [100,10]
        X_VLa_bottom_bond_location      = [100,60]

        X_VHa                    = [160,60,140,40,160,20,180,40,160,60]
        X_VHa_label_location     = [140,10]
        X_VHa_bottom_bond_location      = [160,60]

        X_VLb                    = [480,60,600,40,580,20,560,40,480,60]
        X_VLb_label_location     = [480,10]
        X_VLb_bottom_bond_location      = [480,60]

        X_VHb                    = [500,60,480,40,500,20,520,40,500,60]
        X_VHb_label_location     = [500,10]
        X_VHb_bottom_bond_location      = [500,60]

        VHa_X_VHb                = [300,180,280,140,300,100,360,100,380,140,360,180,300,180]
        VHa_X_VHb_label_location = [320,80]
        VHa_X_VHb_top_bond_location      = [300,100]
        VHa_X_VHb_bottom_bond_location      = [360,180]

        X_CH3a                   = [280,480,300,500,280,520,260,500]
        X_CH3a_label_location    = [290,520]
        X_CH3a_top_bond_location = [290,480]

        X_CH3b                   = [380,480,400,500,380,520,360,500]
        X_CH3b_label_location    = [380,520]
        X_CH3b_top_bond_location = [380,480]


        CH1_CH2bonda = CH1a_bottom_bond_location+CH2a_top_bond_location
        CH1_CH2bondb = CH1b_bottom_bond_location+CH2b_top_bond_location
        Saltbridge_a  =[271,262,388,262]
        Saltbridge_labela = [260,274]
        Saltbridge_b  = [288,272,374,272]
        Saltbridge_labelb = [400,274]

#########Check VH/VL ADC Coordinates################



        coordinates_list      = []
        Bonds                 = []
        Label_Locations       = []
        Label_Text            = []
        domain_names          = []
        labels = []
        keyslist = list(dictionary.keys())
        for i in range(len(dictionary)):
            mod = ""
            if i == 0:
                if "V" in keyslist[i]:
                    if "+" in keyslist[i]:
                        domain = keyslist[i]
                        domain = re.sub("\+","",str(domain))
                        coordinates= eval(str(domain))
                        mod = " +"
                    elif "_" in keyslist[i]:
                        domain = keyslist[i]
                        domain = re.sub("\_","",str(domain))
                        coordinates= eval(str(domain))
                        mod = " -"
                    else:
                        domain = keyslist[i]
                        coordinates = eval(str(keyslist[i]))

                    coordinates_list.append(coordinates)
                    Label_Locations.append((eval(str(domain+"_label_location"))))
                    location = dictionary.get(keyslist[i])[0]
                    Label_Text.append(str(location)+mod)
                    domain_names.append(keyslist[i])

                elif "X" in keyslist[i]:
                    coordinates = eval(str(keyslist[i]+"_"+keyslist[i+1]))
                    coordinates_list.append(coordinates)
                    location = dictionary.get(keyslist[i])[0]
                    Label_Text.append(str(location))
                    Label_Locations.append((eval(str(keyslist[i]+"_"+keyslist[i+1]+"_label_location"))))
                    domain_names.append(keyslist[i])


            elif i > 0 and "V" in keyslist[i]:
                if re.search('\d', keyslist[i]) is None:
                    if keyslist[i-1] == "X":
                        domain = keyslist[i]
                        coordinates = eval(str(keyslist[i]))
                        coordinates_list.append(coordinates)
                        Label_Locations.append((eval(str(keyslist[i]+"_label_location"))))
                        if i == 1:
                            bonds = bondmaker(str("X_"+domain),domain,chain)
                        elif i > 1:
                            bonds = bondmaker(str("VHa_X_VHb"),domain,chain)
                        Bonds.append(eval(bonds[0])+eval(bonds[1]))
                        location = dictionary.get(keyslist[i])[0]
                        Label_Text.append(str(location))
                        domain_names.append(keyslist[i])

                    else:
                        domain = keyslist[i]
                        coordinates = eval(str(keyslist[i]))
                        coordinates_list.append(coordinates)
                        Label_Locations.append((eval(str(keyslist[i]+"_label_location"))))
                        bonds = bondmaker(str(keyslist[i-1]),domain,chain)
                        Bonds.append(eval(bonds[0])+eval(bonds[1]))
                        location = dictionary.get(keyslist[i])[0]
                        Label_Text.append(str(location))
                        domain_names.append(keyslist[i])

                elif "VH2" in keyslist[i]:
                    if "CH3" in keyslist[i-1]:
                        domain = "CH4"
                        coordinates = eval(str(keyslist[i]+chain))
                        coordinates_list.append(coordinates)
                        Label_Locations.append((eval(str(domain+chain+"_label_location"))))
                        bonds = bondmaker(str(keyslist[i-1]),domain,chain)
                        Bonds.append(eval(bonds[0])+eval(bonds[1]))
                        location = dictionary.get(keyslist[i])[0]
                        Label_Text.append(str(location))
                        domain_names.append(keyslist[i])

                    elif "CH4" in keyslist[i-1]:
                        domain = "CH5"
                        coordinates = eval(str("VH2_2"+chain))
                        coordinates_list.append(coordinates)
                        Label_Locations.append((eval(str(domain+chain+"_label_location"))))
                        bonds = bondmaker(str(keyslist[i-1]),domain,chain)
                        Bonds.append(eval(bonds[0])+eval(bonds[1]))
                        location = dictionary.get(keyslist[i])[0]
                        Label_Text.append(str(location))
                        domain_names.append(keyslist[i])

                elif "VL2" in keyslist[i]:
                    if  keyslist[i-1] == "VH2" and keyslist[i-2] == "CH3"  :
                        domain = "VL2"
                        coordinates = eval(str("VL2"+chain))
                        coordinates_list.append(coordinates)
                        Label_Locations.append((eval(str(domain+chain+"_label_location"))))
                        bonds = bondmaker("CH4",domain,chain)
                        Bonds.append(eval(bonds[0])+eval(bonds[1]))
                        location = dictionary.get(keyslist[i])[0]
                        Label_Text.append(str(location))
                        domain_names.append(keyslist[i])

                    elif "VH2" in keyslist[i-1] and "CH4" in keyslist[i-2]:
                        domain = "VL2_2"
                        coordinates = eval(str("VL2_2"+chain))
                        coordinates_list.append(coordinates)
                        Label_Locations.append((eval(str(domain+chain+"_label_location"))))
                        bonds = bondmaker("CH5",domain,chain)
                        Bonds.append(eval(bonds[0])+eval(bonds[1]))
                        location = dictionary.get(keyslist[i])[0]
                        Label_Text.append(str(location))
                        domain_names.append(keyslist[i])

                    elif "CL2" in keyslist[i-1]:
                        domain = "VL2_2"
                        coordinates = eval(str("VL2_2"+chain))
                        coordinates_list.append(coordinates)
                        Label_Locations.append((eval(str(domain+chain+"_label_location"))))
                        bonds = bondmaker(keyslist[i-1],domain,chain)
                        Bonds.append(eval(bonds[0])+eval(bonds[1]))
                        location = dictionary.get(keyslist[i])[0]
                        Label_Text.append(str(location))
                        domain_names.append(keyslist[i])

                elif "CL2" in keyslist[i]:
                    if dictionary == VHa_chain or  dictionary == VHb_chain:
                        domain = "CL2"
                        coordinates = eval(str("CL2"+chain))
                        coordinates_list.append(coordinates)
                        Label_Locations.append((eval(str(domain+chain+"_label_location"))))
                        bonds = bondmaker(keyslist[i-1],domain,chain)
                        Bonds.append(eval(bonds[0])+eval(bonds[1]))
                        location = dictionary.get(keyslist[i])[0]
                        domain_names.append(keyslist[i])



            elif keyslist[i] == "H":
                location = dictionary.get(keyslist[i])[0]
                Label_Locations.append((eval(str("Saltbridge_label"+chain))))
                Label_Text.append(str(location))

            elif i > 0 and "V" not in keyslist[i]:
                domain = keyslist[i]
                previous_domain = str(keyslist[i-1])
                previous_domain = re.sub("\*|\+|\_|\>|\@","",previous_domain)
                previous_domain_2 = str(keyslist[i-2])
                previous_domain_2 = re.sub("\*|\+|\_|\>|\@","",previous_domain_2)
                mod    = ""
                print(domain, previous_domain,previous_domain_2)
                if "@" in domain:
                    domain = re.sub("@","",str(domain))
                    coordinates= eval(str(domain+chain+"_notch"))
                    mod = " notch"

                elif ">" in keyslist[i]:
                    domain = re.sub(">","",str(domain))
                    coordinates= eval(str(domain+chain+"_inward_notch"))
                    mod = " inward notch"

                elif "+" in keyslist[i]:
                    domain = re.sub("\+","",str(domain))
                    coordinates= eval(str(domain+chain))
                    mod = " +"

                elif "_" in keyslist[i]:
                    domain = re.sub("\_","",str(domain))
                    coordinates= eval(str(domain+chain))
                    mod = " -"

                elif "*" in keyslist[i]:
                    domain = re.sub("\*.*","",str(domain))
                    coordinates= eval(str(domain+chain))
                    try:
                        mod = str(keyslist[i].split("*")[1])
                    except IndexError:
                        mod = ""


                elif keyslist[i] == "X":
                    if keyslist[i-1] == "VHa":
                        if keyslist[i+1] == "VHb":
                            coordinates = VHa_X_VHb
                            domain = "VHa_X_VHb"


                            chain = ""
                    elif "CH3" in keyslist[i-1]:
                        domain = "X_CH3"
                        coordinates = eval(domain+chain)


                else:
                    coordinates = eval(str(domain+chain))

                if previous_domain == "H":
                    print("yes")
                    bonds = bondmaker(previous_domain_2,domain,chain)
                else:
                    bonds = bondmaker(previous_domain,domain,chain)
                Bonds.append(eval(bonds[0])+eval(bonds[1]))
                coordinates_list.append(coordinates)
                Label_Locations.append((eval(str(domain+chain+"_label_location"))))
                location = dictionary.get(keyslist[i])[0]
                Label_Text.append(str(location)+mod)
                domain_names.append(keyslist[i])

        return(coordinates_list, Bonds,Label_Text, Label_Locations, domain_names)

    def Count_Salt_bridges(Salt_bridge_count):
        Disulphide_bonds = []
        Saltbridge_a  =[271,262,388,262]
        Saltbridge_b  = [288,272,374,272]
        if Salt_bridge_count   == 0:
            pass
        if Salt_bridge_count == 1:
            Disulphide_bonds = [Saltbridge_a]
        elif Salt_bridge_count == 2:
            Disulphide_bonds = [Saltbridge_a, Saltbridge_b]
        return Disulphide_bonds

    def Horizontal_bonds(VHa_VLa_bond_count,VHb_VLb_bond_count,CH1a_CL1a_bond_count,CH1b_CLb_bond_count):
        Horizontal_bonds = []
        VLa_bottom_bond_location  = [157,174]
        VHa_bottom_bond_location = [217,174]
        CLa_bottom_bond_location  = [200,254]
        CH1a_bottom_bond_location = [260,254]
        VLb_bottom_bond_location = [503,174]
        VHb_bottom_bond_location = [443,174]
        CLb_bottom_bond_location = [460,254]
        CH1b_bottom_bond_location = [400,254]
        VHa_VLa_bond_location   = VHa_bottom_bond_location+VLa_bottom_bond_location
        VHb_VLb_bond_location   = VHb_bottom_bond_location+VLb_bottom_bond_location
        CH1a_CL1a_bond_location = CH1a_bottom_bond_location+CLa_bottom_bond_location
        CH1b_CL1b_bond_location = CH1b_bottom_bond_location+CLb_bottom_bond_location
        if VHa_VLa_bond_count != "":
            Horizontal_bonds.append(VHa_VLa_bond_location)
        if VHb_VLb_bond_count!= "":
            Horizontal_bonds.append(VHb_VLb_bond_location)
        if CH1a_CLa_bond_count != "":
            Horizontal_bonds.append(CH1a_CL1a_bond_location)
        if CH1b_CLb_bond_count !="":
            Horizontal_bonds.append(CH1b_CL1b_bond_location)
        return Horizontal_bonds



    H_Coordinates   = GetChainCoordinates(VHa_chain,"a")[0] , GetChainCoordinates(VHb_chain,"b")[0]
    L_Coordinates   = GetChainCoordinates(VLa_chain,"a")[0] , GetChainCoordinates(VLb_chain,"b")[0]
    Bonds           = GetChainCoordinates(VLa_chain,"a")[1] + GetChainCoordinates(VHa_chain,"a")[1] + GetChainCoordinates(VLb_chain,"b")[1] + GetChainCoordinates(VHb_chain,"b")[1]
    Horizontal_bonds= Count_Salt_bridges(Salt_bridge_count)+ Horizontal_bonds(VHa_VLa_bond_count,VHb_VLb_bond_count,CH1a_CLa_bond_count,CH1b_CLb_bond_count)
    Label_Text      = GetChainCoordinates(VLa_chain,"a")[2] + GetChainCoordinates(VHa_chain,"a")[2] + GetChainCoordinates(VLb_chain,"b")[2] + GetChainCoordinates(VHb_chain,"b")[2]
    Label_spot      = GetChainCoordinates(VLa_chain,"a")[3] + GetChainCoordinates(VHa_chain,"a")[3] + GetChainCoordinates(VLb_chain,"b")[3] + GetChainCoordinates(VHb_chain,"b")[3]
    H_domain_names  = GetChainCoordinates(VHa_chain,"a")[4] , GetChainCoordinates(VHb_chain,"b")[4]
    L_domain_names  = GetChainCoordinates(VLb_chain,"b")[4] + GetChainCoordinates(VLb_chain,"b")[4]
    return(H_Coordinates,L_Coordinates,Bonds,Horizontal_bonds,Label_Text,Label_spot, H_domain_names,L_domain_names)



def render(chains_list,canvas):
    canvas.delete("all")
    H_Polygons      = chains_list[0]
    L_Polygons      = chains_list[1]
    Bonds           = chains_list[2]
    Horizontal_bonds= chains_list[3]
    Label_Text      = chains_list[4]
    Label_positions = chains_list[5]
    H_domain_names  = chains_list[6]
    L_domain_names  = chains_list[7]
    print(H_Polygons)
    print(L_Polygons)


    for i in range(len(Bonds)):
        canvas.create_line(Bonds[i], fill='#000000', width = 2)
    for i in range(len(Horizontal_bonds)):
        canvas.create_line(Horizontal_bonds[i], fill='#000000', width = 2)

    for i in range(len(H_Polygons)):
        for j in range(len(H_Polygons[i])):
            #tag = str(H_domain_names[i][j])
            canvas.create_polygon(H_Polygons[i][j], outline='#000000',fill='#006400', width=2, tags = str(H_domain_names[i][j]))
            #canvas.tag_bind(str(H_domain_names[i][j]),"<Enter>", button_hover_polygon(H_domain_names[i][j]))
            #canvas.tag_bind(str(H_domain_names[i][j]),"<Leave>", button_hover_polygon_leave)
    for i in range(len(L_Polygons)):
        for j in range(len(L_Polygons[i])):
            #tag = str(L_domain_names[i][j])
            canvas.create_polygon(L_Polygons[i][j], outline='#000000',fill='#1f1', width=2)
            #canvas.tag_bind(tag,"<Enter>", button_hover_polygon(tag))
            #canvas.tag_bind(tag,"<Leave>", button_hover_polygon_leave)
    for i in range(len(Label_positions)):
        canvas.create_text(Label_positions[i], text=Label_Text[i])
    #canvas.pack(fill=BOTH, expand=1)


def render_pipeline(canvas):
    entry=textBox.get("1.0","end-1c")
    split_chains = Get_dictionaries(entry)
    coordinates  = Check_interactions(split_chains)
    render(coordinates, canvas)

def button_hover(e):
    button["bg"] = "white"
    status_label.config(text="Get Input")
def button_hover_leave(e):
    button["bg"] = "grey"
    status_label.config(text="")
def button_hover_polygon(label):
    status_label.config(text=label)
def button_hover_polygon_leave(e):
    status_label.config(text="")




###############Main programme#######################

root = tk.Tk()

canvas = tk.Canvas(root, height=HEIGHT, width=WIDTH, bg='#E7E0E6')
canvas.pack()

frame = tk.Frame(root, bg = '#80c1ff',bd=5)
frame.place(relx=0.05,rely=0.1,relwidth=0.4,relheight=0.8,)

###Input box
textBox = tk.Text(frame, font=40)
textBox.place(rely = 0.05, relwidth=1,relheight=0.3)

###Option box
frame2 = tk.Frame(frame, bg = '#D3D3D3')
frame2.place(rely = 0.4, relheight = 0.35, relwidth = 1)

status_label = tk.Label(root, text='test', bd=1)
status_label.place(rely = 0.98, relheight = 0.02, relwidth = 1)

##Big button
button = tk.Button(frame, text = "Get Structure", bg = "grey", font=40, command=lambda: render_pipeline(lower_canvas))
button.place(relx=0.35,rely=0.85,relheight=0.1, relwidth=0.3)
button.bind("<Enter>", button_hover)
button.bind("<Leave>", button_hover_leave)

###Results canvas
lower_frame = tk.Frame(root, bg = '#80c1ff', bd=10)
lower_frame.place(relx=0.45, rely=0.1, relwidth=0.55,relheight=0.8)

lower_canvas = tk.Canvas(lower_frame)
lower_canvas.place(relheight=1,relwidth=1)




root.mainloop()
