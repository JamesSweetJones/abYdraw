#!/usr/bin/python
import re
import sys
import os
import tkinter as tk

HEIGHT = 1200
WIDTH  = 1600


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

    ###Split chains into dictionaries
    y       = re.sub("\s","",x)
    if "|" in y:
        splitx  = y.split("|")
    else:
        splitx = [y]
    if len(splitx) == 4:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a']
    elif len(splitx) == 2:
        chains  = ['VH.b','VH.a']
    elif len(splitx) == 1:
        chains  = ['VH.a']
    elif len(splitx) ==5:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a', 'fragment1']
    elif len(splitx) ==6:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a', 'fragment1', 'fragment2']


    chain_count = len(splitx)
    VHa     = {}
    VHb     = {}
    VLa     = {}
    VLb     = {}
    disulphide_bridges   = ""
    VHa_VLa_bonds  = ""
    VHb_VLb_bonds  = ""
    CH1a_CL1a_bonds=""
    CH1b_CL1b_bonds=""
    fragment1={}
    fragment2={}
    fragment3={}
    fragment4={}



    for i in range(len(chains)):
        chain = splitx[i].split("-")
        if i < 4:
            if  "VH" in chain[0]:
                if chain[1] != "L":
                    dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                elif chain[1] == "L":
                    if "VH" in chain[2]:
                        dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                    elif "VL" in chain[2]:
                        dict = str(re.sub("\.|\+|\_","",str(chains[i])))

            elif "VL" in chain[0]:
                if  chain[1] != "L":
                    dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                elif chain[1] == "L":
                    if len(chain) > 2:
                        if "VH" in chain[2]:
                            dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                        elif "VL" in chain[2]:
                            dict = str(re.sub("\.|\+|\_","",str(chains[i])))

            elif "X" in chain[0] and "VH" in chain[1]:
                 dict = str(re.sub("\.|\+|\_","",str(chains[i])))
        elif i >= 4:
            dict_number = i-3
            dict = str("fragment"+str(dict_number))


        for j in range(len(chain)):
            domain   =  re.findall("^CH[0-9][@+>_!*]\(.*\)\[.*\]|^CH[0-9]\(.*\)\[.*\]|^CH[0-9][@+>_!*]|^CL[0-9][@+>_!*]\(.*\)\[.*\]|^CL[0-9]\(.*\)\[.*\]|^CL[0-9][@+>_!*]|^CL[@+>_!*]|^CH[0-9]|^VL[1-9].[abcd]|^VL.[abcd]|^VH\.[abcd]|^VL[1-9]|^VH[1-9]\.[abcd]|^VH[1-9]|VH[+_*]\.[abcd]|VL[+_*]\.[abcd]|^CL[0-9]|^CL|^VH|^VL|^H|^X|^L", str(chain[j]))
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


            location    = []
            locationstr =  re.findall("\((.*?)\)", str(chain[j]))
            locationstr =  str(re.sub("\[|\'|\]","", str(locationstr)))
            locationstr =  str(re.sub(":",",", str(locationstr)))
            locationstr = locationstr.split(",")

            if domain != "X" and domain !="L" and domain != "":
                for x in range(len(locationstr)):
                    location.append(int(locationstr[x]))
            elif domain == "X":
                location = re.findall("\[(.*?)\]", str(chain[j]))
            elif domain == "L":
                domain = "L"+str([j])
                location.append("")

            if re.findall("\{.*?\}", str(chain[j])) != []:
                disulphide_bridges = re.findall("\{.*?\}", str(chain[j]))
                disulphide_bridges = re.sub("\{|\'|\}|\[|\]","", str(disulphide_bridges))
                disulphide_bridges = int(disulphide_bridges)
                location = [location,disulphide_bridges]


            if dict == "VHa" and domain !="":
                VHa[domain] = location
            elif dict == "VHb" and domain !="":
                VHb[domain] = location
            elif dict == "VLa" and domain !="":
                VLa[domain] = location
            elif dict == "VLb" and domain !="":
                VLb[domain] = location
            elif dict == "fragment1" and domain !="":
                fragment1[domain] = location
            elif dict == "fragment2" and domain !="":
                fragment2[domain] = location
            elif dict == "fragment3" and domain !="":
                fragment4[domain] = location
            elif dict == "fragment4" and domain !="":
                fragment1[domain] = location
            else:
                continue

    return(VHa,VLa,VHb,VLb,chain_count,fragment1,fragment2,fragment3,fragment4)

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
    chain_count     = chains_list[4]
    fragment1       = chains_list[5]
    fragment2       = chains_list[6]
    fragment3       = chains_list[7]
    fragment4       = chains_list[8]
    All_positions_and_chains    ={}
    extra_disulphide_bridges          ={}
    completed_disulphidebridges=[]



    def innie_or_outie(chain,Light_chain_check, chain_count):
        innie_or_outie_list = []
        keyslist = list(chain.keys())
        for n in range(len(chain)):

            if "V" in keyslist[n]:
                default = ""
                if Light_chain_check == False:
                    default = "outie"
                elif Light_chain_check == True:
                    default = "innie"
                if n == 0:
                    try:
                        if "L[" in keyslist[n+1]:
                            if len(chain.get(keyslist[n])) == 1:
                                innie_or_outie_list.append("Single_Fv_Chain")
                            try:
                                if  chain.get(keyslist[n])[0] == (chain.get(keyslist[n+2])[1]):
                                    if Light_chain_check == False:
                                        innie_or_outie_list.append("innie")
                                    elif Light_chain_check == True:
                                        innie_or_outie_list.append("innie")
                                elif chain.get(keyslist[n])[0] != (chain.get(keyslist[n+2])[1]):
                                    if chain_count == 2:
                                        if Light_chain_check == False:
                                            innie_or_outie_list.append("innie")
                                        elif Light_chain_check == True:
                                            innie_or_outie_list.append("innie")
                                    else:
                                        innie_or_outie_list.append(default)
                            except IndexError:
                                innie_or_outie_list.append("Single_Fv_Chain")
                        elif "L[" not in keyslist[n+1]:
                            if len(chain.get(keyslist[n])) == 1:
                                innie_or_outie_list.append("Single_Fv_Chain")
                            elif chain_count == 2 or chain_count == 1:
                                if Light_chain_check == False:
                                    innie_or_outie_list.append("innie")
                            else:
                                innie_or_outie_list.append(default)
                    except IndexError:
                        continue



                elif n > 0 and n < len(chain):
                    if "L[" in keyslist[n-1]:
                        if len(chain.get(keyslist[n])) == 1:
                            innie_or_outie_list.append("Single_Fv_Chain")
                        elif chain.get(keyslist[n-2])[0] == (chain.get(keyslist[n])[1]):

                            if innie_or_outie_list[-2] == "outie":
                                innie_or_outie_list.append("innie")
                            elif innie_or_outie_list[-2] == "innie":
                                innie_or_outie_list.append("outie")
                            else:
                                innie_or_outie_list.append(default)

                        elif  chain.get(keyslist[n])[0] != (chain.get(keyslist[n])[1]):
                            if "CL" in keyslist[n-2] and Light_chain_check == True:
                                innie_or_outie_list.append("outie")
                            elif innie_or_outie_list[-2] == "outie":
                                innie_or_outie_list.append("outie")
                            elif innie_or_outie_list[-2] == "innie":
                                innie_or_outie_list.append("innie")
                            else:
                                innie_or_outie_list.append(default)
                        else:
                            innie_or_outie_list.append(default)

                    else:
                        innie_or_outie_list.append(default)
                elif n == len(chain):
                    innie_or_outie_list.append("innie")
            elif "V" not in keyslist[n]:
                innie_or_outie_list.append("constant")


        return(innie_or_outie_list)

    def VH_VL_check(current_chain_pos,current_dictionary,Light_chain_check,Heavy_chain,Light_chain,have,lookingfor):
        interaction = ""
        if Light_chain_check == False:
            keyslist = list(Light_chain.keys())
            candidate_chain = Light_chain
        elif Light_chain_check == True:
            keyslist = list(Heavy_chain.keys())
            candidate_chain = Heavy_chain
        current_chain_loc =  current_dictionary.get(current_chain_pos)[0]
        for x in range(len(keyslist)):
            if lookingfor in keyslist[x]:
                candidate_chain_loc = candidate_chain.get(keyslist[x])[1]
                if current_chain_loc == candidate_chain_loc:
                    interaction = "H_L"
                elif current_chain_loc != candidate_chain_loc:
                    interaction = "H_H"

        return(interaction)

    def find_the_fragment(fragment_start,All_positions_and_chains):
        try:
            start_coordinates = All_positions_and_chains.get(int(fragment_start))
            return(start_coordinates)
        except IndexError:
            pass

    def domainmaker(startx,starty,righthanded,slant,V,direction,X,mod,interaction):
        if V == False and  direction == "constant" and mod == "" and X == False:

            firstx      = startx
            firsty      = starty
            if righthanded == False:
                secondx = firstx - 20
            elif righthanded == True:
                secondx = firstx + 20
            secondy     = firsty
            if slant == True and righthanded == False:
                thirdx      = secondx + 40
            elif slant == True and righthanded == True:
                thirdx      = secondx - 40
            else:
                thirdx      = secondx
            thirdy          = secondy+ 75
            if righthanded == False:
                fourthx     = thirdx + 20
            elif righthanded == True:
                fourthx     = thirdx - 20
            fourthy         = thirdy
            if righthanded == False:
                fifthx     = fourthx + 20
            elif righthanded == True:
                fifthx     = fourthx - 20
            fifthy         = fourthy
            if slant and righthanded == False:
                sixthx     = fifthx - 40
            elif slant and righthanded == True:
                sixthx     = fifthx + 40
            else:
                sixthx     = fifthx
            sixthy     = fifthy-75
            coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy]
        elif V == True and direction == "innie":
            firstx      = startx
            firsty      = starty
            if righthanded == False:
                secondx = firstx - 20
            elif righthanded == True:
                secondx = firstx + 20
            secondy     = firsty
            if slant == True and righthanded == False:
                thirdx      = secondx + 40
            elif slant == True and righthanded == True:
                thirdx      = secondx - 40
            else:
                thirdx      = secondx
            thirdy          = secondy+ 75
            if righthanded == False:
                fourthx     = thirdx + 20
            elif righthanded == True:
                fourthx     = thirdx - 20
            fourthy         = thirdy
            if righthanded == False:
                fifthx     = fourthx + 20
            elif righthanded == True:
                fifthx     = fourthx - 20
            fifthy         = fourthy
            if slant and righthanded == False:
                sixthx     = fifthx - 30
            elif slant and righthanded == True:
                sixthx     = fifthx + 30
            else:
                sixthx     = fifthx
            sixthy     = fifthy-55
            if righthanded == False:
                seventhx = sixthx - 20
            elif righthanded == True:
                seventhx = sixthx + 20
            seventhy    = sixthy
            coordinates = [firstx,firsty,secondx,secondy,thirdx,thirdy,fourthx,fourthy,fifthx,fifthy,sixthx,sixthy,seventhx,seventhy]
        elif V == True and direction == "outie":
            firstx      = startx
            firsty      = starty
            if righthanded == False:
                secondx = firstx + 20
            elif righthanded == True:
                secondx = firstx - 20
            secondy     = firsty
            if slant == True and righthanded == False:
                thirdx      = secondx + 40
            elif slant == True and righthanded == True:
                thirdx      = secondx - 40
            else:
                thirdx      = secondx
            thirdy          = secondy+ 75
            if righthanded == False:
                fourthx     = thirdx - 20
            elif righthanded == True:
                fourthx     = thirdx + 20
            fourthy         = thirdy
            if righthanded == False:
                fifthx     = fourthx - 20
            elif righthanded == True:
                fifthx     = fourthx + 20
            fifthy         = fourthy
            if slant and righthanded == False:
                sixthx     = fifthx - 30
            elif slant and righthanded == True:
                sixthx     = fifthx + 30
            else:
                sixthx     = fifthx
            sixthy     = fifthy-55
            if righthanded == False:
                seventhx = sixthx + 20
            elif righthanded == True:
                seventhx = sixthx - 20
            seventhy    = sixthy
            coordinates = [firstx,firsty,secondx,secondy,thirdx,thirdy,fourthx,fourthy,fifthx,fifthy,sixthx,sixthy,seventhx,seventhy]
        elif V == True and direction == "Single_Fv_Chain":
            firstx      = startx
            firsty      = starty
            if righthanded == False:
                secondx = firstx + 20
            elif righthanded == True:
                secondx = firstx - 20
            secondy     = firsty
            if slant == True and righthanded == False:
                thirdx      = secondx + 40
            elif slant == True and righthanded == True:
                thirdx      = secondx - 40
            else:
                thirdx      = secondx
            thirdy          = secondy+ 75
            if righthanded == False:
                fourthx     = thirdx - 20
            elif righthanded == True:
                fourthx     = thirdx + 20
            fourthy         = thirdy
            if righthanded == False:
                fifthx     = fourthx - 20
            elif righthanded == True:
                fifthx     = fourthx + 20
            fifthy         = fourthy
            if slant and righthanded == False:
                sixthx     = fifthx - 40
            elif slant and righthanded == True:
                sixthx     = fifthx + 40
            else:
                sixthx     = fifthx
            sixthy       = firsty
            if slant == False:
                seventhx    = firstx
            elif slant == True and righthanded == False:
                seventhx    = firstx+10
            elif slant == True and righthanded == True:
                seventhx    = firstx-10
            seventhy    = firsty+20

            eighthx     = secondx
            eighthy     = secondy
            coordinates = [secondx,secondy,secondx,secondy,thirdx,thirdy,fourthx,fourthy,fifthx,fifthy,sixthx,sixthy,seventhx,seventhy]

        elif V == False  and mod == "@" and direction == "outie" and X == False:
            firstx      = startx
            firsty      = starty
            if righthanded == False:
                secondx = firstx - 20
            elif righthanded == True:
                secondx = firstx + 20
            secondy     = firsty
            if slant == True and righthanded == False:
                thirdx      = secondx + 40
            elif slant == True and righthanded == True:
                thirdx      = secondx - 40
            else:
                thirdx      = secondx
            thirdy          = secondy+ 75
            if righthanded == False:
                fourthx     = thirdx + 20
            elif righthanded == True:
                fourthx     = thirdx - 20
            fourthy         = thirdy
            if righthanded == False:
                fifthx     = fourthx + 20
            elif righthanded == True:
                fifthx     = fourthx - 20
            fifthy         = fourthy
            if slant== True:
                sixthx     = fourthx-20
            else:
                sixthx     = fourthx
            sixthy     = fifthy-37
            if righthanded == False:
                seventhx   =  firstx + 20
            elif righthanded == True:
                seventhx   =  firstx - 20
            seventhy   =  firsty
            coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]
        elif V == False  and mod == "@" and direction == "innie" and X == False:
            firstx      = startx
            firsty      = starty
            if righthanded == False:
                secondx = firstx + 20
            elif righthanded == True:
                secondx = firstx - 20
            secondy     = firsty
            if slant == True and righthanded == False:
                thirdx      = secondx + 40
            elif slant == True and righthanded == True:
                thirdx      = secondx - 40
            else:
                thirdx      = secondx
            thirdy          = secondy+ 75
            if righthanded == False:
                fourthx     = thirdx - 20
            elif righthanded == True:
                fourthx     = thirdx + 20
            fourthy         = thirdy
            if righthanded == False:
                fifthx     = fourthx - 20
            elif righthanded == True:
                fifthx     = fourthx + 20
            fifthy         = fourthy
            if slant == True:
                sixthx     = fifthx
            elif slant == False and righthanded == False:
                sixthx     = fifthx + 20
            elif slant == False and righthanded == True:
                sixthx     = fifthx - 20
            sixthy     = fifthy-37
            if righthanded == False:
                seventhx   =  firstx - 20
            elif righthanded == True:
                seventhx   =  firstx + 20
            seventhy   =  firsty
            coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]
        elif V == False  and mod == ">" and direction == "outie" and X == False:
            firstx      = startx
            firsty      = starty
            if righthanded == False:
                secondx = firstx + 20
            elif righthanded == True:
                secondx = firstx - 20
            secondy     = firsty
            if slant == True and righthanded == False:
                thirdx      = secondx + 40
            elif slant == True and righthanded == True:
                thirdx      = secondx - 40
            else:
                thirdx      = secondx
            thirdy          = secondy+ 75
            if righthanded == False:
                fourthx     = thirdx - 20
            elif righthanded == True:
                fourthx     = thirdx + 20
            fourthy         = thirdy
            if righthanded == False:
                fifthx     = fourthx - 20
            elif righthanded == True:
                fifthx     = fourthx + 20
            fifthy         = fourthy
            if slant == True and righthanded == False:
                sixthx     =  firstx - 20
            elif slant == True and righthanded == True:
                seventhx   =  firstx + 20
            elif slant == False and righthanded == False:
                sixthx     = fifthx + 40
            elif slant == False and righthanded == True:
                sixthx     = fifthx - 40
            sixthy     = fifthy-37
            if righthanded == False:
                seventhx   =  firstx - 20
            elif righthanded == True:
                seventhx   =  firstx + 20
            seventhy   =  firsty
            coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]
        elif V == False  and mod == ">" and direction == "innie" and X == False:
            firstx      = startx
            firsty      = starty
            if righthanded == False:
                secondx = firstx - 20
            elif righthanded == True:
                secondx = firstx + 20
            secondy     = firsty
            if slant == True and righthanded == False:
                thirdx      = secondx + 40
            elif slant == True and righthanded == True:
                thirdx      = secondx - 40
            else:
                thirdx      = secondx
            thirdy          = secondy+ 75
            if righthanded == False:
                fourthx     = thirdx + 20
            elif righthanded == True:
                fourthx     = thirdx - 20
            fourthy         = thirdy
            if righthanded == False:
                fifthx     = fourthx + 20
            elif righthanded == True:
                fifthx     = fourthx - 20
            fifthy         = fourthy
            if slant == True:
                sixthx     = fifthx
            elif slant == False and righthanded == False:
                sixthx     = fifthx + 40
            elif slant == False and righthanded == True:
                sixthx     = fifthx - 40
            sixthy     = fifthy-37
            if righthanded == False:
                seventhx   =  firstx + 20
            elif righthanded == True:
                seventhx   =  firstx - 20
            seventhy   =  firsty
            coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]

        elif X == True:
            firstx      = startx
            firsty      = starty
            if righthanded == False:
                secondx     =   startx-30
            elif righthanded== True:
                secondx     =   startx+30
            secondy         =   firsty+37
            if righthanded == False:
                thirdx     =   secondx+30
            elif righthanded== True:
                thirdx     =   secondx-30
            thirdy         =   secondy+ 38
            fourthx         =   thirdx
            fourthy         =   thirdy
            if righthanded == False:
                fifthx      =   fourthx+30
            elif righthanded== True:
                fifthx     =    fourthx-30
            fifthy         =    fourthy-38

            coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy]

        top_bond    = [firstx,firsty+20]
        bottom_bond = [fourthx, fourthy]
        if slant == True and righthanded == False:
            Labelbond   = [firstx+20, (secondy+thirdy)/2]
        elif slant == True and righthanded == True:
            Labelbond   = [firstx-20, (secondy+thirdy)/2]
        else:
            Labelbond   = [firstx, (secondy+thirdy)/2]
        return(coordinates, bottom_bond,top_bond,Labelbond)



    def renderchains(dictionary,startx,starty):
        chain                   = []
        coordinates_list_heavy_a= []
        coordinates_list_light_a= []
        coordinates_list_heavy_b= []
        coordinates_list_light_b= []
        coordinates_list_heavy_c= []
        coordinates_list_light_c= []
        coordinates_list_heavy_d= []
        coordinates_list_light_d= []
        bonds = []
        disulphidebridge1 = []
        disulphidebridge2 = []
        Location_Text=[]
        text_coordinates=[]
        Domain_Text = []

        if chain_count >= 4:
            slant = True
        else:
            slant = False
        before_H = True
        Build_in = True
        Build_out = False
        righthanded = False
        Light_chain_check = False
        if dictionary == VHb_chain or dictionary == VLb_chain:
            Heavy_chain = VHb_chain
            Light_chain = VLb_chain
            righthanded = True
        elif dictionary == VHa_chain or dictionary == VLa_chain:
            Heavy_chain = VHa_chain
            Light_chain = VLa_chain
            righthanded = False
        if dictionary == VLa_chain or dictionary == VLb_chain or dictionary == fragment1 or dictionary == fragment2:
            Light_chain_check = True
        if dictionary == fragment1 or dictionary == fragment2:
            slant = False

        innie_or_outie_list = innie_or_outie(dictionary, Light_chain_check, chain_count)


        for i in range(len(dictionary)):

            keyslist = list(dictionary.keys())
            print(keyslist[i])
            V  = False
            X  = False
            direction = str(innie_or_outie_list[i])
            interaction = ""
            mod= ""
            mod_label=""
            if "V" in keyslist[i]:
                V  = True
            elif "X" in keyslist[i]:
                X  = True



            try:
                location = dictionary.get(keyslist[i])[0][0]
                disulphide_bridge = dictionary.get(keyslist[i])[1]
                location    = dictionary.get(keyslist[i])[0]
                dictionary[keyslist[i]] = location

            except:
                disulphide_bridge = ""


            if "@" in keyslist[i]:
                mod = "@"
                if Light_chain_check == True or chain_count==2:
                    direction = "outie"
                elif Light_chain_check == False:
                    if before_H == True:
                        direction = "innie"
                    elif before_H==False:
                        direction = "outie"
            elif ">" in keyslist[i]:
                mod = ">"
                if Light_chain_check == True or chain_count==2:
                    direction = "innie"
                elif Light_chain_check == False:
                    if before_H == True:
                        direction = "outie"
                    elif before_H==False:
                        direction = "innie"

            elif "+" in keyslist[i]:
                mod_label = "+"
            elif "_" in keyslist[i]:
                mod_label = "-"
            elif "*" in keyslist[i]:
                try:
                    mod_label = str(keyslist[i].split("*")[1])
                except IndexError:
                    mod = ""
            Domain_name = str(re.sub("\@|\>|\<","",str(keyslist[i])))
            print(keyslist[i],direction)



            if i == 0:
                getcoordinates = domainmaker(startx,starty, righthanded,slant,V,direction,X,mod,interaction)
                #if VHa_VLa_bond_count  < 0:
                #    VHa_VLa_bond += getcoordinates[1]
                #elif VHb_VLb_bond_count< 0:
                #    VHa_VLa_bond += getcoordinates[1]

            elif i > 0:
                previous_domain = keyslist[i-1]
                previous_chain  = chain[i-1]
                if keyslist[i-1] != "H" and "L[" not in keyslist[i-1] and "X" not in keyslist[i-1]:
                    previous_number = (dictionary.get(previous_domain)[0])+1
                elif "X" in keyslist[i-1]  and i-1 == 0:
                    previous_number = (dictionary.get(keyslist[i])[0])
                    dictionary[keyslist[i-1]] = [previous_number,0]
                    previous_domain = (keyslist[i])


                if previous_domain == "H":
                    previous_domain = keyslist[i-2]
                    previous_chain  = chain[i-2]
                    previous_number = (dictionary.get(previous_domain)[0])+2


                if keyslist[i-1] == "CL":
                    slant = False

                if keyslist[i] == "H":
                    before_H = False
                    if slant == True and righthanded == False:
                        disulphidebridge1 += [bottom_bond[0]+7,bottom_bond[1]+15]
                    elif slant == True and righthanded == True:
                        disulphidebridge1 += [bottom_bond[0]-7,bottom_bond[1]+15]
                    else:
                        disulphidebridge1 += [bottom_bond[0],bottom_bond[1]+15]
                    slant = False

                elif "X" in keyslist[i]:
                    previous_chain = chain[i-2]
                    if righthanded == False:
                        getcoordinates = domainmaker((previous_chain[0]+100),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction)
                    elif righthanded == True:
                        getcoordinates = domainmaker((previous_chain[0]-100),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction)


                elif keyslist[i-1] == "H" and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                    if righthanded == True:
                        getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction)
                    elif righthanded == False:
                        getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction)

                    top_bond = getcoordinates[2]
                    if righthanded == False:
                        disulphidebridge2 += [top_bond[0]-5,top_bond[1]-10]
                    elif  righthanded == True:
                        disulphidebridge2 += [top_bond[0]+5,top_bond[1]-10]
                    else:
                        disulphidebridge2 += [top_bond[0],top_bond[1]-10]

                elif len(dictionary.get(keyslist[i-1])) ==1:
                    if slant == True and righthanded == True:
                        getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                    elif slant == True and righthanded == False:
                        getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                    else:
                        getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                    Build_in  = False
                    Build_out = True

                elif keyslist[i] != "H" and "L[" not in keyslist[i-1] and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                    if slant == True and righthanded == True:
                        getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                    elif slant == True and righthanded == False:
                        getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                    else:
                        getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)

                    Build_in  = False
                    Build_out = True

##Linker
                elif "L[" in keyslist[i-1]:
                    previous_chain = chain[i-2]
                    previous_domain = keyslist[i-2]

                    if len(dictionary.get(keyslist[i])) == 1:
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                        else:
                            getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                        Build_in  = False
                        Build_out = True
                    elif len(dictionary.get(keyslist[i-2])) ==1:
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                        else:
                            getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                        Build_in  = False
                        Build_out = True
##Build up
                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == (dictionary.get(keyslist[0])[1]):
                        getcoordinates = domainmaker((previous_chain[0]),(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction)

##Build down
                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):

                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)
                        else:
                            getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction)


                        Build_in  = False
                        Build_out = True
##Build across
                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] == (dictionary.get(previous_domain)[1]):
                        if Build_in == True:
                            if righthanded == False:
                                getcoordinates = domainmaker((previous_chain[0]+60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction)
                            elif righthanded == True:
                                getcoordinates = domainmaker((previous_chain[0]-60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction)
                        elif Build_out == True:
                            if righthanded == False:
                                getcoordinates = domainmaker((previous_chain[0]-60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction)
                            elif righthanded == True:
                                getcoordinates = domainmaker((previous_chain[0]+60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction)

##Build diagonally
                    elif dictionary.get(keyslist[i])[0] != previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                        if righthanded  == True:
                            righthanded = False
                        elif righthanded == False:
                            righthanded = True
                        if Build_in == True:
                            if righthanded == False:
                                getcoordinates = domainmaker((previous_chain[6]-100),(previous_chain[7]+20), righthanded,slant,V,direction,X,mod,interaction)
                            elif righthanded == True:
                                getcoordinates = domainmaker((previous_chain[6]+100),(previous_chain[7]+20), righthanded,slant,V,direction,X,mod,interaction)
                        if Build_in == False:
                            if righthanded == False:
                                getcoordinates = domainmaker((previous_chain[6]+100),(previous_chain[7]+20), righthanded,slant,V,direction,X,mod,interaction)
                            elif righthanded == True:
                                getcoordinates = domainmaker((previous_chain[6]-100),(previous_chain[7]+20), righthanded,slant,V,direction,X,mod,interaction)


##append coordinates to chains
            if keyslist[i] == "H" or "L[" in keyslist[i]:
                chain.append([])
                bonds.append([])
            else:
                chain.append(getcoordinates[0])
            if i > 0:
                top_bond = getcoordinates[2]
                bonds.append(bottom_bond + top_bond)
            bottom_bond = getcoordinates[1]




###append domains to right list
            if keyslist[i] != "H" and "L[" not in keyslist[i]:
                if "a" in keyslist[i] or "a" in keyslist[i-1]:
                    if "VL" in keyslist[i] or "CL" in keyslist[i]:
                        coordinates_list_light_a.append(getcoordinates[0])
                    else:
                        coordinates_list_heavy_a.append(getcoordinates[0])
                elif "b" in keyslist[i] or "b" in keyslist[i-1]:
                    if "VL" in keyslist[i] or "CL" in keyslist[i]:
                        coordinates_list_light_b.append(getcoordinates[0])
                    else:
                        coordinates_list_heavy_b.append(getcoordinates[0])
                elif "c" in keyslist[i]:
                    if "VL" in keyslist[i] or "CL" in keyslist[i]:
                        coordinates_list_light_c.append(getcoordinates[0])
                    else:
                        coordinates_list_heavy_c.append(getcoordinates[0])
                elif "d" in keyslist[i]:
                    if "VL" in keyslist[i] or "CL" in keyslist[i]:
                        coordinates_list_light_d.append(getcoordinates[0])
                    else:
                        coordinates_list_heavy_d.append(getcoordinates[0])
                elif "CL" in keyslist[i]:
                    if dictionary == VHa_chain or dictionary == VLa_chain:
                        coordinates_list_light_a.append(getcoordinates[0])
                    elif dictionary == VHb_chain or dictionary == VLb_chain:
                        coordinates_list_light_b.append(getcoordinates[0])
                elif  dictionary == VHa_chain or dictionary == VLa_chain:
                    coordinates_list_heavy_a.append(getcoordinates[0])
                elif dictionary == VHb_chain or dictionary == VLb_chain:
                    coordinates_list_heavy_b.append(getcoordinates[0])
                elif dictionary == fragment1:
                    if "a" in keyslist[i-1] or "a" in keyslist[i-2] :
                        coordinates_list_heavy_a.append(getcoordinates[0])
                    elif "b" in keyslist[i-1]or "b" in keyslist[i-2] :
                        coordinates_list_heavy_b.append(getcoordinates[0])





##Get labels and positions
            if keyslist[i] != "H" and "L[" not in keyslist[i]:
                Label_Locations = getcoordinates[3]
                text_coordinates.append([Label_Locations])
                text = dictionary.get(keyslist[i])[0]
                Location_Text.append(str(text)+mod_label)
                All_positions_and_chains[text] = [getcoordinates[0], righthanded]
                Domain_Text.append(str(Domain_name)+mod_label)

            elif keyslist[i]=="H":
                Label_Locations = [getcoordinates[3][0],getcoordinates[3][1]+60]
                text_coordinates.append(Label_Locations)
                text = dictionary.get(keyslist[i])[0]
                Location_Text.append(str(text)+mod_label)
                All_positions_and_chains[text] = [Label_Locations, righthanded]
                Domain_Text.append(str(Domain_name)+mod_label)
##Sort extra disulphide bridges
            if disulphide_bridge != "" and keyslist[i] !="H":
                interactor = dictionary.get(keyslist[i])[1]
                extra_disulphide_bridges[interactor] = bottom_bond
                disulphidebridge_keyslist = list(extra_disulphide_bridges.keys())

                for j in range(len(extra_disulphide_bridges)):
                    if int(disulphidebridge_keyslist[j]) == int(location[0]):
                        coordinates = extra_disulphide_bridges.get(disulphidebridge_keyslist[j])+bottom_bond
                        completed_disulphidebridges.append(coordinates)


        chain = [x for x in chain if x != []]
        allbonds = bonds
        allbonds = [x for x in bonds if x != []]
        #allbonds = bonds
        return(coordinates_list_heavy_a,coordinates_list_light_a,coordinates_list_heavy_b,coordinates_list_light_b,coordinates_list_heavy_c,coordinates_list_light_c,coordinates_list_heavy_d,coordinates_list_light_d,allbonds,Location_Text,text_coordinates, disulphidebridge1, disulphidebridge2,completed_disulphidebridges,Domain_Text)

    VLa_Domains_before_CL = 0
    VLb_Domains_before_CL = 0
    for i in range(len(VLa_chain)):
        keyslist = list(VLa_chain.keys())
        if  "CL" not in keyslist[i] and "L[" not in keyslist[i]:
            VLa_Domains_before_CL += 1
        elif "CL" in keyslist[i]:
            break
    for i in range(len(VLb_chain)):
        keyslist = list(VLb_chain.keys())
        if "CL" not in keyslist[i]  and "L[" not in keyslist[i]:
            VLb_Domains_before_CL += 1
        elif "CL" in keyslist[i]:
            break

    VHa_Domains_before_H = 0
    VHb_Domains_before_H = 0
    keyslist = list(VHa_chain.keys())
    for i in range(len(VHa_chain)):
        if keyslist[i] != "H" and "L[" not in keyslist[i]:
            VHa_Domains_before_H += 1
        elif keyslist[i] == "H":
            break


    VLa_startx, VLa_starty = 75, 100
    VHa_startx, VHa_starty = 145, 100
    if chain_count == 2:
        if "L[" in keyslist[1] and VHa_chain.get(keyslist[0])[0] == (VHa_chain.get(keyslist[2])[1]):
            VHa_startx, VHa_starty = 150, 100
        else:
            VHa_startx, VHa_starty = 250, 100
    elif chain_count ==1:
        VHa_startx, VHa_starty = 250,200
    elif chain_count >= 4:
        if VHa_Domains_before_H   == 2 and VLa_Domains_before_CL == 1:
            VLa_startx, VLa_starty = 75, 100
            VHa_startx, VHa_starty = 145, 100
        elif VHa_Domains_before_H == 2 and VLa_Domains_before_CL == 2:
            VLa_startx, VLa_starty = 30,5
            VHa_startx, VHa_starty = 100,5
        elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 1:
            VLa_startx, VLa_starty = 75, 100
            VHa_startx, VHa_starty = 140,5
        elif VHa_Domains_before_H == 4 and VLa_Domains_before_CL == 1:
            VLa_startx, VLa_starty  = 75,100
            VHa_startx, VHa_starty  = 30,5
        #elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 2:
        #    VLa_startx, VLa_starty  = 75,95
        #    VHa_startx, VHa_starty  = 255,95
    if VLa_Domains_before_CL == 1:
        VLa_startx, VLa_starty = 75, 100

    keyslist = list(VHb_chain.keys())
    for i in range(len(VHb_chain)):
        if keyslist[i] != "H" and "L[" not in keyslist[i]:
            VHb_Domains_before_H += 1
        elif keyslist[i] == "H":
            break

    VLb_startx, VLb_starty = 520, 100
    VHb_startx, VHb_starty = 450, 100
    if chain_count == 2:
        if "L[" in keyslist[1] and  VHb_chain.get(keyslist[0])[0] == (VHb_chain.get(keyslist[2])[1]):
            VHb_startx, VHb_starty = 400,100
        else:
            VHb_startx, VHb_starty = 350,100

    elif chain_count >= 4:
        if VHb_Domains_before_H   == 2 and VLa_Domains_before_CL == 1:
            VLb_startx, VLb_starty = 520, 100
            VHb_startx, VHb_starty = 450, 100
        elif VHa_Domains_before_H == 2 and VLa_Domains_before_CL == 2:
            VLb_startx, VLb_starty = 560,5
            VHb_startx, VHb_starty = 480,5
        elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 1:
            VLb_startx, VLb_starty = 515, 100
            VHb_startx, VHb_starty = 440,5
        elif VHa_Domains_before_H == 4 and VLa_Domains_before_CL == 1:
            VLb_startx, VLb_starty  = 520, 100
            VHb_startx, VHb_starty  = 560,5
        #elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 2:
        #    VLb_startx, VLb_starty = 690, 95
        #    VHb_startx, VHb_starty = 490,95
    if VLa_Domains_before_CL == 1:
        VLb_startx, VLb_starty = 515, 100







    VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty)
    VLa_stats = renderchains(VLa_chain,VLa_startx, VLa_starty)
    VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty)
    VLb_stats = renderchains(VLb_chain,VLb_startx, VLb_starty)
    frag1_startx,frag1_starty,frag2_startx,frag2_starty,frag3_startx,frag3_starty,frag4_startx,frag4_starty = 0,0,0,0,0,0,0,0

    if fragment1 != {}:
        fragment1_list = list(fragment1.keys())
        fragment_inter = fragment1.get(fragment1_list[0])[1]
        frag1_start = find_the_fragment(fragment_inter,All_positions_and_chains)
        righthanded = frag1_start[1]
        if righthanded == True:
            frag1_startx=frag1_start[0][0]+50
            frag1_starty=frag1_start[0][1]
        elif righthanded==False:
            frag1_startx=frag1_start[0][0]-50
            frag1_starty=frag1_start[0][1]
    if fragment2 != {}:
        fragment2_list = list(fragment2.keys())
        fragment_inter = fragment2.get(fragment2_list[0])[1]
        frag2_start = find_the_fragment(fragment_inter,All_positions_and_chains)
        righthanded = frag2_start[1]
        if righthanded == True:
            frag2_startx=frag2_start[0][0]+50
            frag2_starty=frag2_start[0][1]
        elif righthanded==False:
            frag2_startx=frag2_start[0][0]-50
            frag2_starty=frag2_start[0][1]
    frag1_stat= renderchains(fragment1,frag1_startx,frag1_starty)
    frag2_stat= renderchains(fragment2,frag2_startx,frag2_starty)
    frag3_stat= renderchains(fragment3,frag3_startx,frag3_starty)
    frag4_stat= renderchains(fragment4,frag4_startx,frag4_starty)



    Heavy_Domains_a     = VHa_stats[0] + VLa_stats[0] + VHb_stats[0] + VLb_stats[0] + frag1_stat[0] + frag2_stat[0] + frag3_stat[0] + frag4_stat[0]
    Light_Domains_a     = VHa_stats[1] + VLa_stats[1] + VHb_stats[1] + VLb_stats[1] + frag1_stat[1] + frag2_stat[1] + frag3_stat[1] + frag4_stat[1]
    Heavy_Domains_b     = VHa_stats[2] + VLa_stats[2] + VHb_stats[2] + VLb_stats[2] + frag1_stat[2] + frag2_stat[2] + frag3_stat[2] + frag4_stat[2]
    Light_Domains_b     = VHa_stats[3] + VLa_stats[3] + VHb_stats[3] + VLb_stats[3] + frag1_stat[3] + frag2_stat[3] + frag3_stat[3] + frag4_stat[3]
    Heavy_Domains_c     = VHa_stats[4] + VLa_stats[4] + VHb_stats[4] + VLb_stats[4] + frag1_stat[4] + frag2_stat[4] + frag3_stat[4] + frag4_stat[4]
    Light_Domains_c     = VHa_stats[5] + VLa_stats[5] + VHb_stats[5] + VLb_stats[5] + frag1_stat[5] + frag2_stat[5] + frag3_stat[5] + frag4_stat[5]
    Heavy_Domains_d     = VHa_stats[6] + VLa_stats[6] + VHb_stats[6] + VLb_stats[6] + frag1_stat[6] + frag2_stat[6] + frag3_stat[6] + frag4_stat[6]
    Light_Domains_d     = VHa_stats[7] + VLa_stats[7] + VHb_stats[7] + VLb_stats[7] + frag1_stat[7] + frag2_stat[7] + frag3_stat[7] + frag4_stat[7]
    Bonds               = VLa_stats[8] + VHa_stats[8] + VLb_stats[8] + VHb_stats[8] + frag1_stat[8] + frag2_stat[8] + frag3_stat[8] + frag4_stat[8]
    Bonds               = [x for x in Bonds if x != []]
    disulphide_bridges        = [VHa_stats[11] + VHb_stats[11]] + [VHa_stats[12] + VHb_stats[12]] + completed_disulphidebridges
    disulphide_bridges        = [x for x in disulphide_bridges if x != []]
    Label_Text          = VLa_stats[9] + VHa_stats[9] + VLb_stats[9] + VHb_stats[9] + frag1_stat[9] + frag2_stat[9] + frag3_stat[9] + frag4_stat[9]
    Label_spot          = VLa_stats[10] +VHa_stats[10] +VLb_stats[10] +VHb_stats[10]+ frag1_stat[10]+ frag2_stat[10]+ frag3_stat[10]+ frag4_stat[10]
    Domain_Text         = VLa_stats[14] + VHa_stats[14] + VLb_stats[14] + VHb_stats[14] + frag1_stat[14] + frag2_stat[14] + frag3_stat[14] + frag4_stat[14]

    return(Heavy_Domains_a,Light_Domains_a,Heavy_Domains_b,Light_Domains_b,Heavy_Domains_c,Light_Domains_c,Heavy_Domains_d,Light_Domains_d,Bonds,disulphide_bridges,Label_Text,Label_spot,Domain_Text)

def render(chains_list,canvas):
    canvas.delete("all")

    Heavy_Domains_a    = chains_list[0]
    Light_Domains_a    = chains_list[1]
    Heavy_Domains_b    = chains_list[2]
    Light_Domains_b    = chains_list[3]
    Heavy_Domains_c    = chains_list[4]
    Light_Domains_c    = chains_list[5]
    Heavy_Domains_d    = chains_list[6]
    Light_Domains_d    = chains_list[7]
    Bonds              = chains_list[8]
    disulphide_bridges       = chains_list[9]
    Label_Text         = chains_list[10]
    Label_positions    = chains_list[11]
    Domain_Text        = chains_list[12]




#Bonds
    for i in range(len(Bonds)):
        canvas.create_line(Bonds[i], fill='#000000', width = 2)
#disulphide_bridge
    if disulphide_bridges != []:
        for i in range(len(disulphide_bridges)):
            canvas.create_line(disulphide_bridges[i], fill='#FF4040', width = 2)

#A domains
    for i in range(len(Heavy_Domains_a)):
        canvas.create_polygon(Heavy_Domains_a[i], outline='#000000',fill='#007ECB', width=2)
    for i in range(len(Light_Domains_a)):
        canvas.create_polygon(Light_Domains_a[i], outline='#000000',fill='#73CAFF', width=2)

#B domains
    for i in range(len(Heavy_Domains_b)):
        canvas.create_polygon(Heavy_Domains_b[i], outline='#000000',fill='#FF43EE', width=2)
    for i in range(len(Light_Domains_b)):
        canvas.create_polygon(Light_Domains_b[i], outline='#000000',fill='#F9D3F5', width=2)

#C domains
    for i in range(len(Heavy_Domains_c)):
        canvas.create_polygon(Heavy_Domains_c[i], outline='#000000',fill='#0BD05A', width=2)
    for i in range(len(Light_Domains_c)):
        canvas.create_polygon(Light_Domains_c[i], outline='#000000',fill='#B9FAD3', width=2)
#
#D domains
    for i in range(len(Heavy_Domains_d)):
        canvas.create_polygon(Heavy_Domains_d[i], outline='#000000',fill='#D9DE4A', width=2)
    for i in range(len(Light_Domains_d)):
        canvas.create_polygon(Light_Domains_d[i], outline='#000000',fill='#E2F562', width=2)

    for i in range(len(Label_positions)):
        canvas.create_text(Label_positions[i], text=Domain_Text[i])
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
frame.place(relwidth=1,relheight=1)

###Input box
textBox = tk.Text(frame, font=40)
textBox.place(relx=0.01,rely = 0.05, relwidth=0.4,relheight=0.3)

###Option box
frame2 = tk.Frame(frame, bg = '#D3D3D3')
frame2.place(relx=0.01, rely = 0.45,relheight = 0.35, relwidth = 0.4)

status_label = tk.Label(root, text='test', bd=1)
status_label.place(rely = 0.98, relheight = 0.02, relwidth = 1)

##Big button
button = tk.Button(frame, text = "Get Structure", bg = "grey", font=40, command=lambda: render_pipeline(lower_canvas))
button.place(relx=0.1,rely=0.85,relheight=0.1, relwidth=0.2)
button.bind("<Enter>", button_hover)
button.bind("<Leave>", button_hover_leave)

###Results canvas
lower_frame = tk.Frame(root, bg = '#80c1ff', bd=10)
lower_frame.place(relx=0.45, rely=0.05, relwidth=0.55,relheight=0.9)

lower_canvas = tk.Canvas(lower_frame)
lower_canvas.place(relheight=1,relwidth=1)




root.mainloop()
