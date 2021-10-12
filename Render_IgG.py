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
    elif len(splitx) <4:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a', 'X','Y']


    chain_count = len(splitx)
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


    for i in range(len(chains)):

        chain = splitx[i].split("-")
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
        #elif "VH" in chain[0] and
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

            if domain != "X" and domain !="L" and domain != "":
                for x in range(len(locationstr)):
                    location.append(int(locationstr[x]))
            elif domain == "X":
                location = re.findall("\[(.*?)\]", str(chain[j]))
            elif domain == "L":
                domain = "L"+str([j])
                location.append("")

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

    return(VHa,VLa,VHb,VLb,Salt_bridges,VHa_VLa_bonds,VHb_VLb_bonds,CH1a_CL1a_bonds,CH1b_CL1b_bonds,fragment,chain_count)

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
    CH1a_CLa_bond_count   = chains_list[7]
    CH1b_CLb_bond_count   = chains_list[8]
    fragment              = chains_list[9]
    chain_count           = chains_list[10]




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
                            if chain.get(keyslist[n+2])[0] == chain.get(keyslist[n])[0]+1 and chain.get(keyslist[n])[0] == (chain.get(keyslist[n+2])[1]):
                                if Light_chain_check == False:
                                    innie_or_outie_list.append("innie")
                                elif Light_chain_check == True:
                                    innie_or_outie_list.append("innie")
                            elif  chain.get(keyslist[n+2])[0] == chain.get(keyslist[n])[0]+1 and chain.get(keyslist[n])[0] != (chain.get(keyslist[n+2])[1]):
                                if chain_count == 2:
                                    if Light_chain_check == False:
                                        innie_or_outie_list.append("innie")
                                    elif Light_chain_check == True:
                                        innie_or_outie_list.append("innie")
                                else:
                                    innie_or_outie_list.append(default)
                        elif "L[" not in keyslist[n+1]:
                            if chain_count == 2:
                                if Light_chain_check == False:
                                    innie_or_outie_list.append("innie")
                            else:
                                innie_or_outie_list.append(default)
                    except IndexError:
                        continue



                elif n > 0 and n < len(chain):
                    if "L[" in keyslist[n-1]:
                        if chain.get(keyslist[n-2])[0]+1 == chain.get(keyslist[n])[0] and chain.get(keyslist[n-2])[0] == (chain.get(keyslist[n])[1]):

                            if innie_or_outie_list[-2] == "outie":
                                innie_or_outie_list.append("innie")
                            elif innie_or_outie_list[-2] == "innie":
                                innie_or_outie_list.append("outie")
                            else:
                                innie_or_outie_list.append(default)

                        elif  chain.get(keyslist[n-2])[0]+1 == chain.get(keyslist[n])[0] and chain.get(keyslist[n])[0] != (chain.get(keyslist[n])[1]):
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
            #print("innie")
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
            #print("outie")
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
        elif V == False  and mod == "@" and X == False:
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

        elif V == False and mod == ">" and interaction == "H_L" and X == False:
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
        elif V == False and mod == ">" and interaction != "H_L" and  X == False:
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
                sixthx     = fifthx - 20
            elif slant == False and righthanded == True:
                sixthx     = fifthx + 20
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

        top_bond    = [firstx,firsty]
        bottom_bond = [fourthx, fourthy]
        Labelbond   = [firstx, secondy/thirdy]
        return(coordinates, bottom_bond,top_bond,Labelbond)



    def renderchains(dictionary,startx,starty):
        chain                = []
        coordinates_list_heavy= []
        coordinates_list_light= []
        bonds = []
        saltbridge1 = []
        saltbridge2 = []
        H_count     = 0
        Label_Text=[]
        text_coordinates=[]

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
        if dictionary == VLa_chain or dictionary == VLb_chain:
            Light_chain_check = True

        innie_or_outie_list = innie_or_outie(dictionary, Light_chain_check, chain_count)
        #print(dictionary, innie_or_outie_list)


        for i in range(len(dictionary)):

            keyslist = list(dictionary.keys())

            V  = False
            X  = False
            direction = str(innie_or_outie_list[i])
            #print(keyslist[i],direction)
            interaction = ""
            mod= ""
            if "V" in keyslist[i]:
                V  = True
            elif "X" in keyslist[i]:
                X  = True


            if "@" in keyslist[i]:
                mod = "@"
            elif ">" in keyslist[i]:
                mod = ">"
                if "H" in str(dictionary):
                    interaction = VH_VL_check(keyslist[i],dictionary,Light_chain_check,Heavy_chain,Light_chain,">","@")
            elif "+" in keyslist[i]:
                mod = "+"
            elif "_" in keyslist[i]:
                mod = "-"
            elif "*" in keyslist[i]:
                try:
                    mod = str(keyslist[i].split("*")[1])
                except IndexError:
                    mod = ""



            if i == 0:


                getcoordinates = domainmaker(startx,starty, righthanded,slant,V,direction,X,mod,interaction)



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
                    slant = False
                    before_H == False
                    H_count +=1
                    if H_count == 1:
                        saltbridge1 += bottom_bond
                    elif H_count == 2:
                        saltbridge1 += bottom_bond
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


                    if H_count == 1:
                        saltbridge2 += top_bond
                    elif H_count == 2:
                        saltbridge2 += top_bond

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

##Build up
                    if dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == (dictionary.get(keyslist[0])[1]):
                        if "VH" in keyslist[i-2] and "VL" in keyslist[i]:
                            if righthanded == True:
                                righthanded = False
                            elif righthanded==False:
                                righthanded = True
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
                                getcoordinates = domainmaker((previous_chain[0]+50),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction)
                            elif righthanded == True:
                                getcoordinates = domainmaker((previous_chain[0]-50),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction)
                        elif Build_out == True:
                            if righthanded == False:
                                getcoordinates = domainmaker((previous_chain[0]-50),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction)
                            elif righthanded == True:
                                getcoordinates = domainmaker((previous_chain[0]+50),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction)


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
            if "VL" in keyslist[i] or "CL" in keyslist[i]:
                coordinates_list_heavy.append(getcoordinates[0])
            else:
                coordinates_list_light.append(getcoordinates[0])

        chain = [x for x in chain if x != []]
        bonds = [x for x in bonds if x != []]
        allbonds = bonds + saltbridge1 + saltbridge2
        
        return(chain,allbonds,Label_Text,text_coordinates)

    VHa_Domains_before_H = 0
    VHb_Domains_before_H = 0
    for i in range(len(VHa_chain)):
        keyslist = list(VHa_chain.keys())
        if keyslist[i] != "H" and "L[" not in keyslist[i]:
            VHa_Domains_before_H += 1
        elif keyslist[i] == "H":
            break
    for i in range(len(VHb_chain)):
        keyslist = list(VHb_chain.keys())
        if keyslist[i] != "H" and "L[" not in keyslist[i]:
            VHb_Domains_before_H += 1
        elif keyslist[i] == "H":
            break

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
    VLa_startx, VLa_starty = 75, 95
    VHa_startx, VHa_starty = 145, 95
    if VHa_Domains_before_H   == 2 and VLa_Domains_before_CL == 1:
        VLa_startx, VLa_starty = 75, 95
        VHa_startx, VHa_starty = 145, 95
    elif VHa_Domains_before_H == 2 and VLa_Domains_before_CL == 2:
        VLa_startx, VLa_starty = 30,0
        VHa_startx, VHa_starty = 100,0
    elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 1:
        VLa_startx, VLa_starty = 75, 95
        VHa_startx, VHa_starty = 140,0
    elif VHa_Domains_before_H == 4 and VLa_Domains_before_CL == 1:
        VLa_startx, VLa_starty  = 75,95
        VHa_startx, VHa_starty  = 30,0
    #elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 2:
    #    VLa_startx, VLa_starty  = 75,95
    #    VHa_startx, VHa_starty  = 255,95
    if VLa_Domains_before_CL == 1:
        VLa_startx, VLa_starty = 75, 95

    VLb_startx, VLb_starty = 520, 95
    VHb_startx, VHb_starty = 450, 95
    if VHb_Domains_before_H   == 2 and VLa_Domains_before_CL == 1:
        VLb_startx, VLb_starty = 520, 95
        VHb_startx, VHb_starty = 450, 95
    elif VHa_Domains_before_H == 2 and VLa_Domains_before_CL == 2:
        VLb_startx, VLb_starty = 560,0
        VHb_startx, VHb_starty = 480,0
    elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 1:
        VLb_startx, VLb_starty = 515, 95
        VHb_startx, VHb_starty = 440,0
    elif VHa_Domains_before_H == 4 and VLa_Domains_before_CL == 1:
        VLb_startx, VLb_starty  = 520, 95
        VHb_startx, VHb_starty  = 560,0
    #elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 2:
    #    VLb_startx, VLb_starty = 690, 95
    #    VHb_startx, VHb_starty = 490,95
    if VLa_Domains_before_CL == 1:
        VLb_startx, VLb_starty = 515, 95

    Coordinates_a = renderchains(VHa_chain,VHa_startx, VHa_starty)[0] + renderchains(VLa_chain,VLa_startx, VLa_starty)[0]
    Coordinates_b = renderchains(VLb_chain,VLb_startx, VLb_starty)[0] + renderchains(VHb_chain,VHb_startx, VHb_starty)[0]
    Bonds         = renderchains(VLa_chain,VLa_startx, VLa_starty)[1] + renderchains(VHa_chain,VHa_startx, VLa_starty)[1] + renderchains(VLb_chain,VLb_startx, VLb_starty)[1] + renderchains(VHb_chain,VHb_startx, VHb_starty)[1]
    Label_Text    = renderchains(VLa_chain,VLa_startx, VLa_starty)[2] + renderchains(VHa_chain,VHa_startx, VHa_starty)[2] + renderchains(VLb_chain,VLb_startx, VLb_starty)[2] + renderchains(VHb_chain,VHb_startx, VHb_starty)[2]
    Label_spot    = renderchains(VLa_chain,VLa_startx, VLa_starty)[3] + renderchains(VHa_chain,VHa_startx, VHa_starty)[3] + renderchains(VLb_chain,VLb_startx, VLb_starty)[3] + renderchains(VHb_chain,VHb_startx, VHb_starty)[3]

    return(Coordinates_a,Coordinates_b,Bonds,Label_Text,Label_spot)

def render(chains_list,canvas):
    canvas.delete("all")
    Coordinates_a    = chains_list[0]
    Coordinates_b    = chains_list[1]
    Bonds            = chains_list[2]
    Label_Text       = chains_list[3]
    Label_positions  = chains_list[4]





    for i in range(len(Bonds)):
        canvas.create_line(Bonds[i], fill='#000000', width = 2)
    for i in range(len(Coordinates_a)):
        canvas.create_polygon(Coordinates_a[i], outline='#000000',fill='#007ECB', width=2)
    for i in range(len(Coordinates_b)):
        canvas.create_polygon(Coordinates_b[i], outline='#000000',fill='#FF43EE', width=2)


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
